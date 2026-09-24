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


def load() -> dict[str, NDArray[np.float64]]:
    """Read the cache, with a message pointing at this script if it is not there."""
    if not CACHE.exists():
        raise FileNotFoundError(
            f"{CACHE} does not exist. Build it with\n"
            "    uv run --extra euphonic --group docs "
            "python docs/make_interactive_data.py"
        )
    return dict(np.load(CACHE))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="re-run the MD instead of using the cache",
    )
    args = parser.parse_args()

    FIGURES.mkdir(exist_ok=True)
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
