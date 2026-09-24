# SPDX-License-Identifier: GPL-3.0-or-later
"""Command line interface.

Deliberately thin. Every subcommand is a short call into the library, so that anything
the CLI can do is reachable from Python with the same arguments, and the CLI never
becomes the only place a piece of logic lives.

Subcommands are added as pipeline stages land. Currently:

``pdos``
    Stages A-C: trajectory in, velocity spectral density out.
"""

from __future__ import annotations

import argparse
import sys
from typing import TYPE_CHECKING

import numpy as np

from mdins import __version__
from mdins.ir import VelocitySpectralDensity
from mdins.spectral import velocity_spectral_density
from mdins.trajectory import read_velocities

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Sequence

__all__ = ["main"]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mdins",
        description="Inelastic neutron scattering spectra from MD trajectories.",
    )
    parser.add_argument("--version", action="version", version=f"mdins {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    pdos = subparsers.add_parser(
        "pdos",
        help="compute the velocity spectral density from a trajectory",
        description=(
            "Read velocities, remove centre-of-mass drift, and estimate the per-atom "
            "velocity cross-spectral density. The trace of the tensor divided by three "
            "is the atom-projected density of states."
        ),
    )
    pdos.add_argument("trajectory", help="trajectory file carrying velocities")
    pdos.add_argument(
        "--dt",
        type=float,
        required=True,
        help=(
            "interval between the frames in the file, in ps. This is the dump "
            "interval, not the MD integration timestep; a wrong value rescales every "
            "energy in the result. Never inferred."
        ),
    )
    pdos.add_argument(
        "-o", "--output", required=True, help="HDF5 file to write the result to"
    )
    pdos.add_argument("--format", default=None, help="ASE format name")
    pdos.add_argument("--index", default=":", help="ASE frame slice, e.g. '::2'")
    pdos.add_argument(
        "--e-max", type=float, default=None, help="highest output energy in meV"
    )
    pdos.add_argument(
        "--n-bins", type=int, default=1024, help="number of output energy bins"
    )
    pdos.add_argument(
        "--e-resolution",
        type=float,
        default=None,
        help=(
            "energy resolution required in meV; an error if the segment length or "
            "maximum lag cannot deliver it. Bin count is not resolution."
        ),
    )
    pdos.add_argument(
        "--estimator",
        choices=["welch", "vacf"],
        default="welch",
        help="see design.md D5",
    )
    pdos.add_argument("--segment-length", type=int, default=None, help="Welch segments")
    pdos.add_argument("--overlap", type=float, default=0.5, help="Welch overlap")
    pdos.add_argument("--max-lag", type=int, default=None, help="VACF maximum lag")
    pdos.add_argument("--window", default="hann", help="window function")
    pdos.add_argument(
        "--temperature",
        type=float,
        default=None,
        help="MD temperature in K; defaults to the trajectory's kinetic temperature",
    )
    pdos.add_argument(
        "--ensemble",
        default=None,
        help=(
            "ensemble the trajectory was generated in, e.g. NVE or NVT. Recorded in "
            "the output; never guessed, since no format stores it."
        ),
    )
    pdos.add_argument(
        "--by-species",
        action="store_true",
        help="average over atoms of each element before writing",
    )
    pdos.add_argument(
        "--remove-rotation",
        action="store_true",
        help=(
            "also project out rigid-body rotation. Only valid for a non-periodic "
            "system; for a periodic cell there is no well-defined global rotation."
        ),
    )
    return parser


def _run_pdos(args: argparse.Namespace) -> int:
    trajectory = read_velocities(
        args.trajectory,
        dt=args.dt,
        index=args.index,
        format=args.format,
        keep_positions=args.remove_rotation,
    )
    trajectory = trajectory.remove_com_velocity()
    if args.remove_rotation:
        trajectory = trajectory.remove_angular_velocity()

    drift = trajectory.temperature_drift()
    print(
        f"{trajectory.n_frames} frames x {trajectory.n_atoms} atoms, "
        f"{trajectory.duration:.4g} ps, "
        f"Nyquist {trajectory.max_resolvable_energy:.4g} meV",
        file=sys.stderr,
    )
    print(
        f"kinetic temperature {trajectory.temperature():.4g} K, drift {drift:.1%}",
        file=sys.stderr,
    )
    if drift > 0.1:
        print(
            f"warning: kinetic temperature drifts by {drift:.1%} across the run, so "
            "the spectrum averages over states the system was passing through rather "
            "than describing one.",
            file=sys.stderr,
        )

    density = velocity_spectral_density(
        trajectory,
        e_max=args.e_max,
        n_bins=args.n_bins,
        estimator=args.estimator,
        segment_length=args.segment_length,
        overlap=args.overlap,
        window=args.window,
        max_lag=args.max_lag,
        temperature=args.temperature,
        ensemble=args.ensemble,
        e_resolution=args.e_resolution,
    )
    # The estimator's resolution, not the trajectory's: for a segmented estimator the
    # segment is the binding constraint, and it is coarser by the number of segments.
    print(
        f"energy resolution {density.metadata.energy_resolution:.4g} meV "
        f"over {density.n_freq} output bins",
        file=sys.stderr,
    )
    if args.by_species:
        density = density.group_by_species()

    _report_checks(density)
    density.to_hdf5(args.output)
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


def _report_checks(density: VelocitySpectralDensity) -> None:
    """Run the invariants and report them rather than raising.

    A violated sum rule means the result is wrong, but it is more useful to write the
    file and say so than to discard an expensive calculation.
    """
    try:
        density.check_sum_rule()
    except ValueError as error:
        print(f"warning: {error}", file=sys.stderr)
    else:
        print("sum rule satisfied", file=sys.stderr)

    try:
        density.check_positive_semidefinite()
    except ValueError as error:
        print(f"warning: {error}", file=sys.stderr)

    weights = density.pdos().sum(axis=0)
    peak = float(density.frequencies[int(np.argmax(weights))])
    print(f"largest spectral weight at {peak:.4g} meV", file=sys.stderr)


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point. Returns a process exit status."""
    args = _build_parser().parse_args(argv)
    try:
        if args.command == "pdos":
            return _run_pdos(args)
    except ValueError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    raise AssertionError(f"unhandled command {args.command!r}")


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
