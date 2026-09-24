# SPDX-License-Identifier: GPL-3.0-or-later
"""Precompute the grid of results behind ``docs/interactive-validation.ipynb``.

Run from the repository root::

    uv run --extra euphonic --group docs python docs/make_interactive_data.py

The notebook lets a reader move trajectory length, Welch segment length, overlap,
temperature, seed and species and watch the harmonic comparison respond. Running MD
inside a slider callback is not viable — one run is tens of seconds — so every
combination is computed once here and cached in
``docs/figures/interactive-validation-data.npz``; the notebook only ever indexes into
arrays.

Two of the axes are free, and that is what makes the grid affordable:

- **Trajectory length.** The MD is run once at the full length and the velocities are
  kept. Truncating to the first *N* frames and re-estimating is exactly the trajectory a
  shorter run of the same seed would have produced, so the whole length axis costs one
  run rather than one run per point.
- **Segment length and overlap.** These are properties of the estimator, not of the
  dynamics, so they too are re-derived from the stored velocities. The length sweep is
  therefore run at every segment length as well, at no extra MD cost.

What is *not* free is the seed, because in NVE the initial Maxwell-Boltzmann draw fixes
the normal-mode energies for the whole run (see the notebook): each seed is a separate
run, and the scatter between them is what the length axis has to be read against.

This module deliberately re-implements the MD loop of ``tests.lj_reference.run_nve``
instead of calling it. ``run_nve`` returns only the spectral density and discards the
velocities, and it is asserted on by the M2 test suite, so it is not the place to add a
return value. Everything that defines the physics — potential, timestep, dump interval,
equilibration — is imported from there, so the two cannot drift apart.

A second, independent cache — ``docs/figures/interactive-benzene-data.npz`` — holds the
same estimator sweep applied to a MACE-MP NVE trajectory of solid benzene that was
produced outside this repository. That half runs no dynamics at all, only re-estimation,
so it is a few seconds; it is kept in its own file precisely so that rebuilding it does
not drag the six minutes of Lennard-Jones MD along with it.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import time
from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from tests.lj_reference import (  # noqa: E402
    DUMP_EVERY,
    N_FRAMES,
    SEEDS,
    SEGMENT_LENGTH,
    TEMPERATURE,
    TIMESTEP_FS,
    binary_cell,
    equilibrium_lattice_constant,
    harmonic_modes,
    lennard_jones,
    moment,
    moment_uncertainty,
    reference_pdos,
    single_species_cell,
    species_weight,
)

from mdins.provenance import Provenance  # noqa: E402 - needs the path above
from mdins.spectral import velocity_spectral_density  # noqa: E402
from mdins.trajectory import VelocityTrajectory  # noqa: E402
from mdins.units import KB, PLANCK  # noqa: E402

if TYPE_CHECKING:  # pragma: no cover
    from ase import Atoms
    from numpy.typing import NDArray

    from mdins.ir import VelocitySpectralDensity

FIGURES = Path(__file__).parent / "figures"
CACHE = FIGURES / "interactive-validation-data.npz"

SUPERCELL = (3, 3, 3)
BINARY_SUPERCELL = (2, 2, 2)

#: Truncations of the stored velocities, in frames — a factor of sixteen in run length,
#: from 10 ps to 164 ps. Welch needs at least one whole segment, so the shortest usable
#: length depends on the segment length and the grid is ragged; the gaps are stored as
#: NaN.
LENGTHS = (512, 768, 1024, 1536, 2048, 3072, 4096, 5120, 6144, 7168, 8192)

#: Resolution axis. 2048 is the M2 choice; the others bracket it, from far coarser than
#: the width of a mode cluster to finer than a short run can support. The full
#: seed × length grid is evaluated at each of these, so the reader can ask whether the
#: length trend is a property of the dynamics or of the estimator.
SEGMENTS = (256, 512, 1024, 2048, 4096)

#: Overlap axis, at two segment lengths. Overlapping segments are correlated, so this
#: mostly moves the *reported* uncertainty rather than the estimate.
OVERLAPS = (0.0, 0.25, 0.5, 0.75)
OVERLAP_SEGMENTS = (1024, 2048)

#: Window axis, passed straight through to :func:`scipy.signal.get_window`. The boxcar
#: is the instructive one: it has the narrowest main lobe of the four, so it "resolves"
#: best by the reported number, and by far the worst sidelobes, so it leaks the most
#: weight out of the band — which is exactly what contaminates the second moment.
WINDOWS = ("hann", "hamming", "blackman", "boxcar")

#: Target temperatures. The M2 runs sit at 10 K (about T_melt/40); 20 K is still a cold
#: crystal but far enough up for anharmonicity to show against a harmonic reference, and
#: 5 K is the control that says the residual at 10 K is not a bug in the pipeline.
TEMPERATURES = (5.0, 10.0, 20.0)

E_MAX = 20.0
N_BINS = 400


def run_velocities(
    primitive: Atoms,
    supercell: tuple[int, int, int],
    *,
    seed: int = 0,
    temperature: float = TEMPERATURE,
    n_frames: int = N_FRAMES,
) -> VelocityTrajectory:
    """NVE molecular dynamics, keeping the velocities rather than a spectrum.

    Mirrors :func:`tests.lj_reference.run_nve` step for step — same potential, same
    timestep, same dump interval, same 500-step equilibration, same start at twice the
    target temperature because half the kinetic energy goes into potential energy as
    the system equilibrates — but returns the trajectory, so that every estimator
    setting and every truncation can be derived from one run.

    Args:
        primitive: Unit cell to repeat.
        supercell: Repetitions along each axis.
        seed: Seed of the Maxwell-Boltzmann draw. Under NVE this fixes the normal-mode
            energies for the whole run, so it is a physical parameter here and not
            merely a numerical detail.
        temperature: Target temperature in K.
        n_frames: Frames to dump.

    Returns:
        The trajectory, with centre-of-mass motion already projected out.
    """
    from ase import units as ase_units
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.verlet import VelocityVerlet

    atoms = primitive * supercell
    atoms.calc = lennard_jones()

    rng = np.random.default_rng(seed)
    MaxwellBoltzmannDistribution(atoms, temperature_K=2 * temperature, rng=rng)
    dynamics = VelocityVerlet(atoms, timestep=TIMESTEP_FS * ase_units.fs)
    dynamics.run(500)  # let the initial condition forget itself

    velocities = np.empty((n_frames, len(atoms), 3))
    for frame in range(n_frames):
        dynamics.run(DUMP_EVERY)
        velocities[frame] = atoms.get_velocities() * ase_units.fs * 1e3

    return VelocityTrajectory(
        velocities=velocities,
        masses=atoms.get_masses(),
        symbols=atoms.get_chemical_symbols(),
        dt=DUMP_EVERY * TIMESTEP_FS * 1e-3,
        provenance=Provenance(source="ase-lj-nve"),
    ).remove_com_velocity()


def estimate(
    trajectory: VelocityTrajectory,
    *,
    n_frames: int | None = None,
    segment_length: int = SEGMENT_LENGTH,
    overlap: float = 0.5,
    window: str = "hann",
) -> VelocitySpectralDensity:
    """Spectral density from the first ``n_frames`` frames of a stored trajectory.

    Truncation is applied to the velocities, not to the spectrum: the result is exactly
    what a shorter run of the same seed would have given, which is what makes the
    notebook's length slider a real experiment rather than an interpolation.
    """
    truncated = (
        trajectory
        if n_frames is None
        else replace(trajectory, velocities=trajectory.velocities[:n_frames])
    )
    return velocity_spectral_density(
        truncated,
        e_max=E_MAX,
        n_bins=N_BINS,
        segment_length=segment_length,
        overlap=overlap,
        window=window,
        temperature=truncated.temperature(),
        ensemble="NVE",
    )


#: Keys produced by :func:`summarise`, and the shape of each per estimate.
SUMMARY_SHAPES: dict[str, tuple[int, ...]] = {
    "weight": (N_BINS,),
    "error": (N_BINS,),
    "m1": (),
    "m2": (),
    "s1": (),
    "s2": (),
    "n_segments": (),
    "resolution": (),
    "temperature": (),
}


def summarise(density: VelocitySpectralDensity, species: str) -> dict[str, Any]:
    """Per-species pDOS, its Welch error, and the first two moments with uncertainties.

    The moments and their uncertainties come from the same functions the M2 tests use,
    so a number read off the notebook is the number the test suite would assert on.
    """
    weight, error = species_weight(density, species)
    energies = density.frequencies
    return {
        "weight": weight,
        "error": error,
        "m1": moment(energies, weight, 1),
        "m2": moment(energies, weight, 2),
        "s1": moment_uncertainty(energies, weight, error, 1),
        "s2": moment_uncertainty(energies, weight, error, 2),
        "n_segments": float(density.metadata.n_segments or 0),
        "resolution": float(density.metadata.energy_resolution or np.nan),
        "temperature": float(density.temperature_md),
    }


def missing() -> dict[str, Any]:
    """Placeholder for a grid point the estimator cannot deliver.

    The seed × length × segment grid is ragged — a 512-frame truncation cannot be cut
    into 2048-frame segments — and filling the holes with NaN keeps it a rectangular
    array that Plotly can index, with the gaps simply not drawn.
    """
    return {key: np.full(shape, np.nan) for key, shape in SUMMARY_SHAPES.items()}


def stack(rows: list[dict[str, Any]], prefix: str, shape: tuple[int, ...]) -> dict:
    """Fold a flat list of :func:`summarise` outputs into a keyed, rectangular grid."""
    return {
        f"{prefix}_{key}": np.array([row[key] for row in rows]).reshape(
            (*shape, *trailing)
        )
        for key, trailing in SUMMARY_SHAPES.items()
    }


def compute(tmp: Path) -> dict[str, Any]:
    """Every run and every derived estimate the notebook needs."""
    data: dict[str, Any] = {
        "energies": np.linspace(0.0, E_MAX, N_BINS),
        "seeds": np.array(SEEDS),
        "lengths": np.array(LENGTHS),
        "segments": np.array(SEGMENTS),
        "overlaps": np.array(OVERLAPS),
        "overlap_segments": np.array(OVERLAP_SEGMENTS),
        "windows": np.array(WINDOWS),
        "temperatures": np.array(TEMPERATURES),
        "dt": DUMP_EVERY * TIMESTEP_FS * 1e-3,
        "reference_segment": SEGMENT_LENGTH,
    }

    lattice_constant = equilibrium_lattice_constant()
    data["lattice_constant"] = lattice_constant
    print(f"equilibrium lattice constant {lattice_constant:.6f} Å")

    single = single_species_cell(lattice_constant)
    single_modes = harmonic_modes(single, SUPERCELL, str(tmp / "lj"))
    modes = single_modes.frequencies.to("meV").magnitude.ravel()
    data["single_modes"] = modes[modes > 1e-3]

    # -- seed x length x segment, the headline grid ----------------------------------
    trajectories = {}
    for seed in SEEDS:
        start = time.perf_counter()
        trajectories[seed] = run_velocities(single, SUPERCELL, seed=seed)
        print(f"single-species MD, seed {seed}: {time.perf_counter() - start:.1f} s")

    start = time.perf_counter()
    rows = [
        summarise(estimate(trajectories[seed], n_frames=n, segment_length=s), "Ar")
        if n >= s
        else missing()
        for seed in SEEDS
        for n in LENGTHS
        for s in SEGMENTS
    ]
    data |= stack(rows, "grid", (len(SEEDS), len(LENGTHS), len(SEGMENTS)))
    print(f"{len(rows)} Welch estimates: {time.perf_counter() - start:.1f} s")

    # -- overlap, from seed 0 at full length -----------------------------------------
    reference_trajectory = trajectories[SEEDS[0]]
    rows = [
        summarise(estimate(reference_trajectory, segment_length=s, overlap=o), "Ar")
        for s in OVERLAP_SEGMENTS
        for o in OVERLAPS
    ]
    data |= stack(rows, "overlap", (len(OVERLAP_SEGMENTS), len(OVERLAPS)))

    # -- window, the other half of the same trade-off ----------------------------------
    rows = [
        summarise(
            estimate(reference_trajectory, segment_length=s, window=w),
            "Ar",
        )
        for w in WINDOWS
        for s in SEGMENTS
    ]
    data |= stack(rows, "window", (len(WINDOWS), len(SEGMENTS)))

    # -- temperature -------------------------------------------------------------------
    rows = []
    for target in TEMPERATURES:
        if target == TEMPERATURE:
            trajectory = reference_trajectory
        else:
            start = time.perf_counter()
            trajectory = run_velocities(single, SUPERCELL, temperature=target)
            print(f"single-species MD, {target} K: {time.perf_counter() - start:.1f} s")
        rows.append(summarise(estimate(trajectory), "Ar"))
    data |= stack(rows, "temperature", (len(TEMPERATURES),))

    # -- mixed species ---------------------------------------------------------------
    binary_primitive = binary_cell(lattice_constant)
    binary_modes = harmonic_modes(
        binary_primitive, BINARY_SUPERCELL, str(tmp / "ljmix")
    )
    binary_reference = reference_pdos(binary_modes, binary_primitive)

    start = time.perf_counter()
    binary_trajectory = run_velocities(binary_primitive, BINARY_SUPERCELL)
    print(f"binary MD: {time.perf_counter() - start:.1f} s")

    for species in ("Ar", "Kr"):
        rows = [
            summarise(
                estimate(binary_trajectory, n_frames=n, segment_length=s), species
            )
            if n >= s
            else missing()
            for n in LENGTHS
            for s in SEGMENTS
        ]
        data |= stack(rows, f"binary_{species}", (len(LENGTHS), len(SEGMENTS)))
        frequencies, weights = binary_reference[species]
        data[f"binary_{species}_ref_energies"] = frequencies
        data[f"binary_{species}_ref_weights"] = weights

    return data


# =====================================================================================
# Solid benzene: the same estimator sweep on a real MD trajectory
# =====================================================================================
#
# Everything above is a Lennard-Jones crystal run inside this repository against a
# Euphonic reference. This half is the opposite situation: a trajectory somebody else
# produced, with no harmonic reference available, so nothing here can be a validation.
# It is a demonstration that the estimator behaves the same way on real data, and — as
# it turns out — a fairly sharp diagnostic of the trajectory itself.

#: Where the trajectory is expected to sit, relative to the repository root. Not
#: committed: 12.9 MB of extxyz is too much to carry in git, and unlike the caches it
#: cannot be regenerated from anything in this repository.
BENZENE_TRAJECTORY = (
    ROOT / "docs" / "data" / "benzene_solid-average-nve-T25.0-traj.extxyz"
)

BENZENE_CACHE = FIGURES / "interactive-benzene-data.npz"

#: Frame spacing in ps. **Read from the file, not assumed.** Each frame's ``info``
#: carries ``time`` and ``step``, and a ``units`` dictionary that states ``"time":
#: "fs"``; consecutive frames are ``step`` 0, 100, 200 … at ``time`` 0.0, 100.0,
#: 200.0 fs — a 1 fs integration timestep dumped every 100 steps.
#: :func:`benzene_dt` re-derives this from the file at build time and fails loudly if it
#: disagrees, because every energy in this section scales inversely with it.
BENZENE_DT = 0.1

#: Truncations of the 1001 frames, in frames: 12.8 ps to the full 100.1 ps.
BENZENE_LENGTHS = (128, 192, 256, 384, 512, 768, 1001)

#: Welch segment lengths, in frames. 0.65 meV down to 0.081 meV of nominal resolution.
BENZENE_SEGMENTS = (64, 128, 256, 512)

BENZENE_WINDOWS = ("hann", "hamming", "blackman", "boxcar")
BENZENE_OVERLAPS = (0.0, 0.25, 0.5, 0.75)
BENZENE_OVERLAP_SEGMENT = 256

#: Decimation factors for the aliasing experiment. Throwing away every other frame
#: halves the Nyquist energy, and whatever extra weight then appears below the new limit
#: is weight that was folded down — measured rather than argued.
BENZENE_DECIMATIONS = (1, 2)

BENZENE_N_BINS = 400
BENZENE_SPECIES = ("C", "H")

#: Shared grid for the decimation panel, just inside the dt = 0.2 ps Nyquist energy so
#: that both rates can be asked for the same interval.
BENZENE_LOW_E_MAX = 0.999 * PLANCK / (2.0 * 2 * BENZENE_DT)


def benzene_dt(path: Path) -> float:
    """Frame spacing in ps, recovered from the trajectory's own ``info`` fields.

    Raises:
        ValueError: If the file does not record times in fs, or if the spacing is not
            uniform, or if it differs from :data:`BENZENE_DT`. Any of those would make
            every energy in this section wrong by a constant factor with no other
            symptom, so none of them is allowed to pass silently.
    """
    from ase.io import read

    frames = read(str(path), index=":3")
    units = frames[0].info.get("units", {})
    if units.get("time") != "fs":
        raise ValueError(f"{path} does not record times in fs: units={units!r}")
    times = np.array([frame.info["time"] for frame in frames])
    spacings = np.diff(times) * 1e-3
    if not np.allclose(spacings, spacings[0], rtol=1e-6):
        raise ValueError(f"non-uniform frame spacing in {path}: {spacings} ps")
    if not np.isclose(spacings[0], BENZENE_DT, rtol=1e-6):
        raise ValueError(
            f"{path} records a {spacings[0]} ps frame spacing but BENZENE_DT is "
            f"{BENZENE_DT} ps. Every energy scales inversely with this; fix the "
            "constant."
        )
    return float(spacings[0])


def benzene_trajectory() -> VelocityTrajectory:
    """The benzene trajectory, drift removed, read through the package's own loader."""
    from mdins.trajectory import read_velocities

    if not BENZENE_TRAJECTORY.exists():
        raise FileNotFoundError(
            f"{BENZENE_TRAJECTORY} does not exist. It is a 12.9 MB MACE-MP NVE run of "
            "solid benzene and is deliberately not committed; place it at that path to "
            "rebuild the benzene half of the notebook."
        )
    dt = benzene_dt(BENZENE_TRAJECTORY)
    return read_velocities(BENZENE_TRAJECTORY, dt=dt).remove_com_velocity()


def benzene_estimate(
    trajectory: VelocityTrajectory,
    *,
    n_frames: int | None = None,
    segment_length: int = 256,
    overlap: float = 0.5,
    window: str = "hann",
    estimator: str = "welch",
    frequencies: NDArray[np.float64] | None = None,
) -> VelocitySpectralDensity:
    """One estimate from the stored benzene velocities.

    ``segment_length`` doubles as the VACF maximum lag, so that the two estimators are
    compared at the same nominal resolution rather than at their respective defaults.
    """
    truncated = (
        trajectory
        if n_frames is None
        else replace(trajectory, velocities=trajectory.velocities[:n_frames])
    )
    if frequencies is None:
        frequencies = np.linspace(
            0.0, 0.999 * truncated.max_resolvable_energy, BENZENE_N_BINS
        )
    return velocity_spectral_density(
        truncated,
        frequencies=frequencies,
        estimator=estimator,
        segment_length=segment_length if estimator == "welch" else None,
        max_lag=segment_length if estimator == "vacf" else None,
        overlap=overlap,
        window=window,
        temperature=truncated.temperature(),
        ensemble="NVE",
    )


#: Keys produced by :func:`benzene_summarise`, and the shape of each per estimate.
BENZENE_SHAPES: dict[str, tuple[int, ...]] = {
    "total": (BENZENE_N_BINS,),
    "C": (BENZENE_N_BINS,),
    "H": (BENZENE_N_BINS,),
    "C_error": (BENZENE_N_BINS,),
    "H_error": (BENZENE_N_BINS,),
    "C_integral": (),
    "H_integral": (),
    "m1": (),
    "m2": (),
    "s1": (),
    "s2": (),
    "n_segments": (),
    "resolution": (),
}


def benzene_summarise(density: VelocitySpectralDensity) -> dict[str, Any]:
    """Absolute per-species pDOS and the integrals the sum rule is stated in.

    Unlike :func:`summarise`, nothing here is normalised to unit area. The whole point
    of this section is that the absolute scale is checkable —
    ``∫ tr P_i dE = 3 k_B T / m_i`` — and normalising would throw away the only
    quantitative handle available without a harmonic reference.

    Per-atom errors within a species are averaged rather than added in quadrature, for
    the reason :func:`tests.lj_reference.species_weight` gives: atoms of one species
    sample the same modes, so their segment-to-segment fluctuations move together.
    """
    from mdins.ir import bin_widths

    energies = density.frequencies
    widths = bin_widths(energies)
    symbols = np.asarray(density.symbols)
    pdos = density.pdos()
    std = (
        density.pdos_std if density.pdos_std is not None else np.full_like(pdos, np.nan)
    )

    row: dict[str, Any] = {
        "total": density.total_pdos(),
        "n_segments": float(density.metadata.n_segments or 0),
        "resolution": float(density.metadata.energy_resolution or np.nan),
    }
    for species in BENZENE_SPECIES:
        index = symbols == species
        row[species] = pdos[index].mean(axis=0)
        row[f"{species}_error"] = std[index].mean(axis=0)
        # 3x because pdos() is tr(P)/3 and the sum rule is stated on the trace.
        row[f"{species}_integral"] = float(3.0 * (row[species] * widths).sum())

    # Moments of the hydrogen projection: the INS-relevant one, and the quantity the
    # estimator sweeps below are read against. Not a physical spectrum — see the
    # notebook on aliasing — but a perfectly good stability metric.
    weight, error = row["H"], row["H_error"]
    row["m1"] = moment(energies, weight, 1)
    row["m2"] = moment(energies, weight, 2)
    row["s1"] = moment_uncertainty(energies, weight, error, 1)
    row["s2"] = moment_uncertainty(energies, weight, error, 2)
    return row


def benzene_stack(rows: list[dict[str, Any]], prefix: str, shape: tuple[int, ...]):
    """Fold :func:`benzene_summarise` outputs into a keyed, rectangular grid."""
    return {
        f"{prefix}_{key}": np.array([row[key] for row in rows]).reshape(
            (*shape, *trailing)
        )
        for key, trailing in BENZENE_SHAPES.items()
    }


def compute_benzene() -> dict[str, Any]:
    """Every benzene estimate the notebook needs. No MD: the trajectory is on disk."""
    trajectory = benzene_trajectory()
    symbols = np.asarray(trajectory.symbols)
    energies = np.linspace(
        0.0, 0.999 * trajectory.max_resolvable_energy, BENZENE_N_BINS
    )

    data: dict[str, Any] = {
        "benzene_energies": energies,
        "benzene_lengths": np.array(BENZENE_LENGTHS),
        "benzene_segments": np.array(BENZENE_SEGMENTS),
        "benzene_windows": np.array(BENZENE_WINDOWS),
        "benzene_overlaps": np.array(BENZENE_OVERLAPS),
        "benzene_overlap_segment": BENZENE_OVERLAP_SEGMENT,
        "benzene_decimations": np.array(BENZENE_DECIMATIONS),
        "benzene_dt": trajectory.dt,
        "benzene_nyquist": trajectory.max_resolvable_energy,
        "benzene_duration": trajectory.duration,
        "benzene_n_frames": trajectory.n_frames,
        "benzene_n_atoms": trajectory.n_atoms,
        "benzene_temperature": trajectory.temperature(),
        "benzene_temperature_drift": trajectory.temperature_drift(),
        "benzene_masses": np.array(
            [trajectory.masses[symbols == s][0] for s in BENZENE_SPECIES]
        ),
        "benzene_counts": np.array(
            [int((symbols == s).sum()) for s in BENZENE_SPECIES]
        ),
    }

    # Kinetic temperature of each species separately, straight from the velocities.
    # This is the number that turns out to matter most, and it needs no estimator.
    for species in BENZENE_SPECIES:
        index = symbols == species
        square = (trajectory.velocities[:, index] ** 2).sum(axis=-1)
        mass = trajectory.masses[index][0]
        data[f"benzene_T_{species}"] = float(mass * square.mean() / (3.0 * KB))
        data[f"benzene_msv_{species}"] = float(square.mean(axis=0).mean())

    start = time.perf_counter()

    # -- length x segment --------------------------------------------------------------
    rows = [
        benzene_summarise(
            benzene_estimate(
                trajectory, n_frames=n, segment_length=s, frequencies=energies
            )
        )
        if n >= s
        else {key: np.full(shape, np.nan) for key, shape in BENZENE_SHAPES.items()}
        for n in BENZENE_LENGTHS
        for s in BENZENE_SEGMENTS
    ]
    data |= benzene_stack(
        rows, "benzene_grid", (len(BENZENE_LENGTHS), len(BENZENE_SEGMENTS))
    )

    # -- window x segment, at full length ----------------------------------------------
    rows = [
        benzene_summarise(
            benzene_estimate(
                trajectory, segment_length=s, window=w, frequencies=energies
            )
        )
        for w in BENZENE_WINDOWS
        for s in BENZENE_SEGMENTS
    ]
    data |= benzene_stack(
        rows, "benzene_window", (len(BENZENE_WINDOWS), len(BENZENE_SEGMENTS))
    )

    # -- overlap -----------------------------------------------------------------------
    rows = [
        benzene_summarise(
            benzene_estimate(
                trajectory,
                segment_length=BENZENE_OVERLAP_SEGMENT,
                overlap=o,
                frequencies=energies,
            )
        )
        for o in BENZENE_OVERLAPS
    ]
    data |= benzene_stack(rows, "benzene_overlap", (len(BENZENE_OVERLAPS),))

    # -- welch against vacf, at matched nominal resolution -----------------------------
    rows = [
        benzene_summarise(
            benzene_estimate(
                trajectory, segment_length=s, estimator=e, frequencies=energies
            )
        )
        for e in ("welch", "vacf")
        for s in BENZENE_SEGMENTS
    ]
    data |= benzene_stack(rows, "benzene_estimator", (2, len(BENZENE_SEGMENTS)))

    # -- decimation: aliasing measured rather than asserted ----------------------------
    low = np.linspace(0.0, BENZENE_LOW_E_MAX, BENZENE_N_BINS)
    data["benzene_low_energies"] = low
    rows = []
    for factor in BENZENE_DECIMATIONS:
        thinned = replace(
            trajectory,
            velocities=trajectory.velocities[::factor],
            dt=trajectory.dt * factor,
        )
        rows.append(
            benzene_summarise(
                benzene_estimate(thinned, segment_length=512 // factor, frequencies=low)
            )
        )
    data |= benzene_stack(rows, "benzene_decimated", (len(BENZENE_DECIMATIONS),))

    print(f"benzene: estimates in {time.perf_counter() - start:.1f} s")
    return data


def load() -> dict[str, NDArray[np.float64]]:
    """Read the cache, with a message pointing at this script if it is not there."""
    if not CACHE.exists():
        raise FileNotFoundError(
            f"{CACHE} does not exist. Build it with\n"
            "    uv run --extra euphonic --group docs "
            "python docs/make_interactive_data.py"
        )
    return dict(np.load(CACHE))


def load_benzene() -> dict[str, NDArray[np.float64]]:
    """Read the benzene cache. Kept separate from :func:`load` so that rebuilding it
    does not mean re-running six minutes of Lennard-Jones molecular dynamics."""
    if not BENZENE_CACHE.exists():
        raise FileNotFoundError(
            f"{BENZENE_CACHE} does not exist. Build it with\n"
            "    uv run --group docs python docs/make_interactive_data.py"
        )
    return dict(np.load(BENZENE_CACHE))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="re-run the MD instead of using the cache",
    )
    args = parser.parse_args()

    FIGURES.mkdir(exist_ok=True)

    if BENZENE_CACHE.exists() and not args.refresh:
        print(f"using cached {BENZENE_CACHE}")
    else:
        start = time.perf_counter()
        np.savez(BENZENE_CACHE, **compute_benzene())
        print(f"wrote {BENZENE_CACHE} in {time.perf_counter() - start:.0f} s")

    if CACHE.exists() and not args.refresh:
        print(f"using cached {CACHE}")
        return

    start = time.perf_counter()
    with tempfile.TemporaryDirectory() as tmp:
        data = compute(Path(tmp))
    np.savez(CACHE, **data)
    print(f"wrote {CACHE} in {time.perf_counter() - start:.0f} s")


if __name__ == "__main__":
    main()
