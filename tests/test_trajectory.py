"""Tests for trajectory ingestion and kinetic temperature validation."""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from ase import units
from ase.build import bulk
from ase.io import Trajectory as AseTrajectoryWriter

from md_ins.trajectory import (
    TrajectoryData,
    calculate_kinetic_temperature,
    load_trajectory,
)


def _trajectory_velocity_for_temperature(
    n_atoms: int, mass_amu: float, temperature_k: float
) -> np.ndarray:
    """Build a single-frame velocity array whose kinetic temperature equals
    ``temperature_k`` exactly (all Cartesian components equal)."""
    v0 = units.fs * np.sqrt(temperature_k * units.kB / mass_amu)
    return np.full((1, n_atoms, 3), v0)


def test_trajectory_data_initialization_validates_shapes(
    rng: np.random.Generator,
) -> None:
    """Task 2.1: TrajectoryData accepts well-shaped synthetic arrays."""
    n_steps, n_atoms = 50, 4
    velocities = rng.standard_normal((n_steps, n_atoms, 3))
    symbols = ["Ar"] * n_atoms
    masses = np.full(n_atoms, 39.948)
    cell = np.eye(3) * 5.26

    traj = TrajectoryData(
        velocities=velocities,
        symbols=symbols,
        masses=masses,
        timestep_fs=1.0,
        cell=cell,
    )
    assert traj.n_steps == n_steps
    assert traj.n_atoms == n_atoms
    assert traj.velocities.shape == (n_steps, n_atoms, 3)
    assert traj.species == ["Ar"]
    assert traj.duration_fs == pytest.approx((n_steps - 1) * 1.0)


def test_trajectory_data_rejects_mismatched_shapes(rng: np.random.Generator) -> None:
    """Task 2.1: shape mismatches raise ValueError."""
    with pytest.raises(ValueError, match="velocities must have shape"):
        TrajectoryData(
            velocities=rng.standard_normal((10, 4)),
            symbols=["Ar"] * 4,
            masses=np.full(4, 39.948),
            timestep_fs=1.0,
        )
    with pytest.raises(ValueError, match="does not match n_atoms"):
        TrajectoryData(
            velocities=rng.standard_normal((10, 4, 3)),
            symbols=["Ar"] * 3,
            masses=np.full(4, 39.948),
            timestep_fs=1.0,
        )
    with pytest.raises(ValueError, match="timestep_fs must be positive"):
        TrajectoryData(
            velocities=rng.standard_normal((10, 4, 3)),
            symbols=["Ar"] * 4,
            masses=np.full(4, 39.948),
            timestep_fs=0.0,
        )


def test_load_trajectory_with_stored_velocities(tmp_path) -> None:
    """Task 2.2: loading a .traj file with velocities converts to Ang/fs."""
    atoms = bulk("Ar", "fcc", a=5.26, cubic=True)
    atoms.calc = None
    n_frames = 20
    timestep_fs = 1.0
    rng = np.random.default_rng(0)
    path = tmp_path / "ar.traj"
    with AseTrajectoryWriter(str(path), mode="w", atoms=atoms) as writer:
        for _ in range(n_frames):
            atoms.positions = (
                atoms.positions + rng.standard_normal(atoms.positions.shape) * 0.01
            )
            atoms.set_velocities(rng.standard_normal(atoms.positions.shape) * 0.001)
            writer.write()

    traj = load_trajectory(path, timestep_fs=timestep_fs)
    assert traj.n_steps == n_frames
    assert traj.n_atoms == atoms.get_global_number_of_atoms()
    assert traj.velocities.shape == (n_frames, traj.n_atoms, 3)
    assert traj.cell is not None
    # velocities should be finite and in a sensible Ang/fs range
    assert np.all(np.isfinite(traj.velocities))


def test_load_trajectory_stride_skips_frames(tmp_path) -> None:
    """Task 2.2: the stride parameter subsamples frames."""
    atoms = bulk("Ar", "fcc", a=5.26, cubic=True)
    rng = np.random.default_rng(1)
    path = tmp_path / "ar_stride.traj"
    n_frames = 30
    with AseTrajectoryWriter(str(path), mode="w", atoms=atoms) as writer:
        for _ in range(n_frames):
            atoms.set_velocities(rng.standard_normal(atoms.positions.shape) * 0.001)
            writer.write()
    traj = load_trajectory(path, timestep_fs=1.0, stride=5)
    assert traj.n_steps == n_frames // 5


def test_load_trajectory_derives_velocities_from_positions(tmp_path) -> None:
    """Task 2.3: position-only trajectories are differentiated, trimmed, warned."""
    atoms = bulk("Ar", "fcc", a=5.26, cubic=True)
    path = tmp_path / "ar_posonly.traj"
    n_frames = 30
    rng = np.random.default_rng(2)
    with AseTrajectoryWriter(str(path), mode="w", atoms=atoms) as writer:
        for _ in range(n_frames):
            atoms.positions = (
                atoms.positions + rng.standard_normal(atoms.positions.shape) * 0.01
            )
            writer.write()  # no velocities stored
    with pytest.warns(UserWarning, match="deriving them by central finite differences"):
        traj = load_trajectory(path, timestep_fs=1.0)
    # central differences drop the first and last frames
    assert traj.n_steps == n_frames - 2
    assert traj.positions is not None
    assert traj.positions.shape[0] == n_frames - 2
    assert np.all(np.isfinite(traj.velocities))


def test_central_difference_matches_harmonic_oscillator() -> None:
    """Task 2.3: derived velocities match the analytic harmonic-oscillator
    velocity for a small timestep."""
    n_steps = 4000
    timestep_fs = 0.5
    amplitude = 0.2  # Angstrom
    omega_per_fs = 0.005  # rad / fs
    times = np.arange(n_steps) * timestep_fs
    # one atom, x-component oscillating
    positions = np.zeros((n_steps, 1, 3))
    positions[:, 0, 0] = amplitude * np.cos(omega_per_fs * times)

    from md_ins.trajectory import _derive_velocities_central

    velocities, trimmed_positions = _derive_velocities_central(
        positions, timestep_fs, cell=None
    )
    # exact velocity at interior times
    interior_times = times[1:-1]
    exact_v = -amplitude * omega_per_fs * np.sin(omega_per_fs * interior_times)
    # central-difference error is O((omega*dt)^2); assert tight agreement
    derived = velocities[:, 0, 0]
    np.testing.assert_allclose(derived, exact_v, rtol=5e-3, atol=1e-6)
    assert trimmed_positions.shape[0] == n_steps - 2


def test_temperature_matching_no_warning() -> None:
    """Task 2.4: kinetic temperature matches nominal -> no warning."""
    mass = 39.948
    velocities = _trajectory_velocity_for_temperature(4, mass, 300.0)
    traj = TrajectoryData(
        velocities=velocities,
        symbols=["Ar"] * 4,
        masses=np.full(4, mass),
        timestep_fs=1.0,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        t_traj = calculate_kinetic_temperature(traj, t_md=300.0)
    assert t_traj == pytest.approx(300.0, rel=1e-6)


def test_temperature_hot_emits_warning() -> None:
    """Task 2.4: a hot trajectory warns about the discrepancy."""
    mass = 39.948
    velocities = _trajectory_velocity_for_temperature(4, mass, 300.0) * 2.0
    traj = TrajectoryData(
        velocities=velocities,
        symbols=["Ar"] * 4,
        masses=np.full(4, mass),
        timestep_fs=1.0,
    )
    with pytest.warns(UserWarning, match="differs from nominal"):
        t_traj = calculate_kinetic_temperature(traj, t_md=300.0)
    assert t_traj == pytest.approx(1200.0, rel=1e-6)


def test_temperature_cold_emits_warning() -> None:
    """Task 2.4: a cold trajectory warns about the discrepancy."""
    mass = 39.948
    velocities = _trajectory_velocity_for_temperature(4, mass, 300.0) * 0.25
    traj = TrajectoryData(
        velocities=velocities,
        symbols=["Ar"] * 4,
        masses=np.full(4, mass),
        timestep_fs=1.0,
    )
    with pytest.warns(UserWarning, match="differs from nominal"):
        t_traj = calculate_kinetic_temperature(traj, t_md=300.0)
    # velocity scales T quadratically: 0.25**2 * 300 = 18.75 K
    assert t_traj == pytest.approx(18.75, rel=1e-6)
