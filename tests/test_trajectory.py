"""Stage A: trajectory input, pre-processing and sampling limits."""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms
from ase import units as ase_units

from tests.conftest import write_trajectory
from mdins.provenance import Provenance
from mdins.trajectory import VelocityTrajectory, read_velocities
from mdins.units import PLANCK


def synthetic(velocities, dt=0.002, masses=None, symbols=None):
    n_atoms = velocities.shape[1]
    return VelocityTrajectory(
        velocities=velocities,
        masses=np.full(n_atoms, 1.0) if masses is None else np.asarray(masses),
        symbols=["H"] * n_atoms if symbols is None else list(symbols),
        dt=dt,
        provenance=Provenance(source="synthetic"),
    )


class TestReading:
    def test_velocities_round_trip_through_ase_units(self, tmp_path):
        """The conversion that is easy to invert: ASE stores Å per ASE time unit."""
        rng = np.random.default_rng(1)
        expected = rng.normal(scale=10.0, size=(5, 5, 3))
        path = tmp_path / "traj.extxyz"
        write_trajectory(path, expected)

        traj = read_velocities(path, dt=0.01)
        # extxyz is a text format; the tolerance is its written precision, not ours.
        np.testing.assert_allclose(traj.velocities, expected, rtol=1e-5)

    def test_reads_masses_and_symbols(self, tmp_path):
        path = tmp_path / "traj.extxyz"
        write_trajectory(path, np.ones((3, 5, 3)))
        traj = read_velocities(path, dt=0.01)
        assert traj.symbols == ["C", "H", "H", "H", "H"]
        assert traj.masses[0] == pytest.approx(12.011, rel=1e-3)
        assert traj.n_atoms == 5
        assert traj.n_frames == 3

    def test_rejects_a_file_without_velocities(self, tmp_path):
        path = tmp_path / "positions.extxyz"
        write_trajectory(path, np.ones((3, 5, 3)), with_velocities=False)
        with pytest.raises(ValueError, match="no velocities"):
            read_velocities(path, dt=0.01)

    def test_rejects_all_zero_velocities(self, tmp_path):
        """A populated but empty momenta array would otherwise give a spectrum of
        zeros rather than an error."""
        path = tmp_path / "zeros.extxyz"
        write_trajectory(path, np.zeros((3, 5, 3)))
        with pytest.raises(ValueError, match="zero"):
            read_velocities(path, dt=0.01)

    def test_error_message_names_the_usable_formats(self, tmp_path):
        path = tmp_path / "positions.extxyz"
        write_trajectory(path, np.ones((3, 5, 3)), with_velocities=False)
        with pytest.raises(ValueError, match=r"GROMACS \.trr"):
            read_velocities(path, dt=0.01)

    def test_index_selects_frames(self, tmp_path):
        path = tmp_path / "traj.extxyz"
        write_trajectory(path, np.ones((10, 5, 3)))
        traj = read_velocities(path, dt=0.02, index="::2")
        assert traj.n_frames == 5

    def test_positions_are_omitted_unless_requested(self, tmp_path):
        path = tmp_path / "traj.extxyz"
        write_trajectory(path, np.ones((3, 5, 3)))
        assert read_velocities(path, dt=0.01).positions is None
        assert read_velocities(path, dt=0.01, keep_positions=True).positions is not None

    def test_provenance_records_the_source_and_settings(self, tmp_path):
        path = tmp_path / "traj.extxyz"
        write_trajectory(path, np.ones((3, 5, 3)))
        traj = read_velocities(path, dt=0.01)
        assert traj.provenance.source == str(path)
        assert any("read_velocities" in step for step in traj.provenance.steps)


class TestValidation:
    def test_rejects_a_single_frame(self):
        with pytest.raises(ValueError, match="at least two frames"):
            synthetic(np.ones((1, 3, 3)))

    def test_rejects_nonpositive_dt(self):
        with pytest.raises(ValueError, match="dt must be positive"):
            synthetic(np.ones((10, 3, 3)), dt=0.0)

    def test_rejects_mismatched_masses(self):
        with pytest.raises(ValueError, match="masses"):
            VelocityTrajectory(
                velocities=np.ones((10, 3, 3)),
                masses=np.ones(2),
                symbols=["H", "H"],
                dt=0.01,
                provenance=Provenance(source="x"),
            )

    def test_rejects_wrong_shape(self):
        with pytest.raises(ValueError, match="n_frames, n_atoms, 3"):
            synthetic(np.ones((10, 3)))


class TestSamplingLimits:
    def test_nyquist_energy(self):
        traj = synthetic(np.ones((100, 2, 3)), dt=0.002)
        assert traj.max_resolvable_energy == pytest.approx(PLANCK / 0.004)
        # A 2 fs dump interval reaches about 1034 meV.
        assert traj.max_resolvable_energy == pytest.approx(1033.9, rel=1e-3)

    def test_energy_resolution_from_duration(self):
        traj = synthetic(np.ones((10000, 2, 3)), dt=0.002)
        assert traj.duration == pytest.approx(20.0)
        assert traj.energy_resolution == pytest.approx(PLANCK / 20.0)

    def test_validate_sampling_accepts_a_reachable_request(self):
        synthetic(np.ones((10000, 2, 3)), dt=0.002).validate_sampling(500.0)

    def test_validate_sampling_rejects_aliasing(self):
        traj = synthetic(np.ones((1000, 2, 3)), dt=0.02)
        with pytest.raises(ValueError, match="Nyquist"):
            traj.validate_sampling(500.0)

    def test_aliasing_message_says_what_dump_interval_would_work(self):
        traj = synthetic(np.ones((1000, 2, 3)), dt=0.02)
        with pytest.raises(ValueError, match="fs or shorter"):
            traj.validate_sampling(500.0)

    def test_validate_sampling_rejects_too_short_a_run(self):
        traj = synthetic(np.ones((100, 2, 3)), dt=0.002)
        with pytest.raises(ValueError, match="resolution"):
            traj.validate_sampling(100.0, e_resolution=0.01)


class TestPreprocessing:
    def test_com_removal_zeroes_the_net_momentum(self):
        rng = np.random.default_rng(2)
        velocities = rng.normal(size=(50, 4, 3)) + np.array([3.0, 0.0, -1.0])
        masses = np.array([12.0, 1.0, 1.0, 1.0])
        traj = synthetic(velocities, masses=masses, symbols=["C", "H", "H", "H"])
        cleaned = traj.remove_com_velocity()
        momentum = np.einsum("a,fad->fd", masses, cleaned.velocities)
        np.testing.assert_allclose(momentum, 0.0, atol=1e-12)

    def test_com_removal_is_mass_weighted(self):
        """An unweighted mean would leave momentum behind for unequal masses."""
        velocities = np.zeros((10, 2, 3))
        velocities[:, 0, 0] = 1.0  # heavy atom moving
        traj = synthetic(velocities, masses=np.array([100.0, 1.0]), symbols=["Xe", "H"])
        cleaned = traj.remove_com_velocity()
        assert cleaned.velocities[0, 0, 0] == pytest.approx(1.0 / 101.0, rel=1e-9)

    def test_com_removal_is_recorded(self):
        traj = synthetic(np.ones((10, 3, 3))).remove_com_velocity()
        assert "remove_com_velocity" in traj.provenance.steps
        assert traj.dof_removed == 3

    def test_angular_removal_needs_positions(self):
        traj = synthetic(np.ones((10, 3, 3)))
        with pytest.raises(ValueError, match="needs positions"):
            traj.remove_angular_velocity()

    def test_angular_removal_kills_rigid_rotation(self):
        """The two removals are complementary: rotation about an origin away from the
        centre of mass carries a net translation, which is the other step's job."""
        rng = np.random.default_rng(3)
        positions = rng.normal(size=(4, 3))
        omega = np.array([0.0, 0.0, 0.7])
        velocities = np.tile(np.cross(omega, positions), (20, 1, 1))
        traj = VelocityTrajectory(
            velocities=velocities,
            masses=np.ones(4),
            symbols=["H"] * 4,
            dt=0.01,
            provenance=Provenance(source="x"),
            positions=np.tile(positions, (20, 1, 1)),
        )
        cleaned = traj.remove_com_velocity().remove_angular_velocity()
        np.testing.assert_allclose(cleaned.velocities, 0.0, atol=1e-12)

    def test_angular_removal_alone_leaves_only_translation(self):
        rng = np.random.default_rng(3)
        positions = rng.normal(size=(4, 3))
        omega = np.array([0.0, 0.0, 0.7])
        velocities = np.tile(np.cross(omega, positions), (20, 1, 1))
        traj = VelocityTrajectory(
            velocities=velocities,
            masses=np.ones(4),
            symbols=["H"] * 4,
            dt=0.01,
            provenance=Provenance(source="x"),
            positions=np.tile(positions, (20, 1, 1)),
        )
        residual = traj.remove_angular_velocity().velocities
        # Identical for every atom, i.e. a pure translation.
        np.testing.assert_allclose(
            residual, np.broadcast_to(residual[:, :1, :], residual.shape), atol=1e-12
        )

    def test_angular_removal_preserves_vibration(self):
        """Rotation goes, internal motion stays."""
        rng = np.random.default_rng(4)
        positions = rng.normal(size=(5, 3))
        vibration = rng.normal(size=(5, 3)) * 0.1
        vibration -= vibration.mean(axis=0)
        omega = np.array([0.2, -0.1, 0.4])
        velocities = np.tile(np.cross(omega, positions) + vibration, (10, 1, 1))
        traj = VelocityTrajectory(
            velocities=velocities,
            masses=np.ones(5),
            symbols=["H"] * 5,
            dt=0.01,
            provenance=Provenance(source="x"),
            positions=np.tile(positions, (10, 1, 1)),
        )
        cleaned = traj.remove_angular_velocity()
        assert np.linalg.norm(cleaned.velocities) > 0.5 * np.linalg.norm(vibration)


class TestDiagnostics:
    def test_temperature_from_maxwell_boltzmann_velocities(self):
        from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

        rng = np.random.default_rng(5)
        atoms = Atoms("Ar" * 300, positions=rng.random((300, 3)))
        frames = []
        for _ in range(40):
            MaxwellBoltzmannDistribution(atoms, temperature_K=150.0, rng=rng)
            frames.append(atoms.get_velocities() * ase_units.fs * 1e3)
        traj = synthetic(
            np.array(frames), masses=atoms.get_masses(), symbols=["Ar"] * 300
        )
        assert traj.temperature() == pytest.approx(150.0, rel=0.05)

    def test_temperature_accounts_for_removed_degrees_of_freedom(self):
        rng = np.random.default_rng(6)
        traj = synthetic(rng.normal(size=(200, 5, 3)))
        before = traj.temperature()
        after = traj.remove_com_velocity().temperature()
        # 15 dof become 12; removing drift lowers the kinetic energy but raises the
        # per-degree-of-freedom temperature back towards the original.
        assert after == pytest.approx(before, rel=0.15)

    def test_drift_is_zero_for_a_stationary_run(self):
        rng = np.random.default_rng(7)
        traj = synthetic(rng.normal(size=(3000, 4, 3)))
        assert traj.temperature_drift() < 0.15

    def test_drift_is_detected_when_the_run_heats_up(self):
        rng = np.random.default_rng(8)
        velocities = rng.normal(size=(3000, 4, 3))
        velocities *= np.linspace(1.0, 3.0, 3000)[:, None, None]
        assert synthetic(velocities).temperature_drift() > 1.0


class TestDegenerateDegreesOfFreedom:
    def test_a_single_atom_projected_onto_itself_has_no_temperature(self):
        """``3N - dof_removed`` is zero, so the kinetic temperature is ``0/0``. Left
        alone this returns NaN, which then disables the sum rule downstream instead of
        failing it."""
        traj = VelocityTrajectory(
            velocities=np.ones((16, 1, 3)),
            masses=np.array([1.008]),
            symbols=["H"],
            dt=0.002,
            provenance=Provenance(source="synthetic"),
        ).remove_com_velocity()

        assert np.abs(traj.velocities).max() == 0.0
        with pytest.raises(ValueError, match="nothing to measure a temperature from"):
            traj.temperature()

    def test_two_atoms_still_have_a_temperature(self):
        """The guard must not catch the smallest system that is actually meaningful."""
        rng = np.random.default_rng(0)
        traj = VelocityTrajectory(
            velocities=rng.normal(size=(64, 2, 3)),
            masses=np.array([1.008, 1.008]),
            symbols=["H", "H"],
            dt=0.002,
            provenance=Provenance(source="synthetic"),
        ).remove_com_velocity()
        assert traj.temperature() > 0.0
