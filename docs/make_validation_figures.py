# SPDX-License-Identifier: GPL-3.0-or-later
"""Generate the figures in ``docs/validation.md``.

Run from the repository root::

    uv run --group docs python docs/make_validation_figures.py

Imports the same harness the M2 tests assert on (``tests/lj_reference.py``), so the
figures cannot drift away from the numbers in the test suite. Roughly two minutes, most
of it molecular dynamics; results are cached in ``docs/figures/validation-data.npz`` and
reused unless ``--refresh`` is given.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from tests.lj_reference import (  # noqa: E402 - needs the path above
    binary_cell,
    equilibrium_lattice_constant,
    harmonic_modes,
    moment,
    moment_uncertainty,
    reference_pdos,
    run_nve,
    single_species_cell,
    species_weight,
)

FIGURES = Path(__file__).parent / "figures"
CACHE = FIGURES / "validation-data.npz"

SUPERCELL = (3, 3, 3)
BINARY_SUPERCELL = (2, 2, 2)

#: The harmonic reference is a set of delta functions and the MD estimate is the true
#: spectrum seen through the Welch window, so they cannot be plotted against each other
#: raw: at any bin width the reference is a spike and the MD a peak of finite height,
#: and the eye reads the height difference as disagreement. The reference is therefore
#: convolved with the estimator's own resolution kernel before plotting, which is the
#: only way the two curves are comparable in height as well as in area. Gaussian of
#: this FWHM, taken from ``metadata.energy_resolution``.
#:
#: ``metadata.energy_resolution`` is the frequency spacing of one Welch segment,
#: ``h/(L·dt)``. The kernel the estimator actually applies is the Hann main lobe, which
#: is 1.44 of those wide at half maximum — so using the reported spacing directly as a
#: FWHM under-broadens the reference by that factor and makes the MD look worse than it
#: is. Measured MD peaks come out 1.5-2.5× the reported spacing, consistent with this
#: plus the finite spread of modes inside each cluster.
#:
#: Nothing in the test suite depends on any of this — the tests compare moments and
#: coarse bins, neither of which is affected by broadening. It is a plotting choice.
HANN_MAIN_LOBE = 1.44
FWHM_TO_SIGMA = 1.0 / 2.3548

ARGON = "#1f77b4"
KRYPTON = "#d62728"
REFERENCE = "0.45"


def compute(tmp: Path) -> dict[str, np.ndarray]:
    """Run both routes for both crystals."""
    lattice_constant = equilibrium_lattice_constant()
    print(f"equilibrium lattice constant {lattice_constant:.6f} Å")

    single = single_species_cell(lattice_constant)
    single_modes = harmonic_modes(single, SUPERCELL, str(tmp / "lj"))
    print("single-species MD ...")
    single_density = run_nve(single, SUPERCELL)

    binary_primitive = binary_cell(lattice_constant)
    binary_modes = harmonic_modes(
        binary_primitive, BINARY_SUPERCELL, str(tmp / "ljmix")
    )
    binary_reference = reference_pdos(binary_modes, binary_primitive)
    print("binary MD ...")
    binary_density = run_nve(binary_primitive, BINARY_SUPERCELL)

    single_weight, single_error = species_weight(single_density, "Ar")
    data = {
        "lattice_constant": lattice_constant,
        "single_energies": single_density.frequencies,
        "single_weight": single_weight,
        "single_error": single_error,
        "single_modes": single_modes.frequencies.to("meV").magnitude.ravel(),
        "single_temperature": single_density.temperature_md,
        "single_resolution": single_density.metadata.energy_resolution,
        "binary_energies": binary_density.frequencies,
        "binary_temperature": binary_density.temperature_md,
        "binary_resolution": binary_density.metadata.energy_resolution,
    }
    for species in ("Ar", "Kr"):
        weight, error = species_weight(binary_density, species)
        frequencies, reference = binary_reference[species]
        data[f"binary_{species}_weight"] = weight
        data[f"binary_{species}_error"] = error
        data[f"binary_{species}_ref_energies"] = frequencies
        data[f"binary_{species}_ref_weights"] = reference
    return data


def broadened(energies, weights, grid, resolution):
    """Weighted mode list as a curve, broadened by the estimator's resolution.

    Each mode becomes a Gaussian of the reported resolution width. The result is
    normalised to unit area so it shares an axis with the MD pDOS.
    """
    sigma = resolution * HANN_MAIN_LOBE * FWHM_TO_SIGMA
    delta = grid[:, None] - np.asarray(energies)[None, :]
    curve = (np.exp(-0.5 * (delta / sigma) ** 2) * np.asarray(weights)).sum(axis=1)
    return curve / np.trapezoid(curve, grid)


def density_curve(energies, weight):
    """Normalise a pDOS to unit area so it can share axes with the reference."""
    return weight / np.trapezoid(weight, energies)


def figure_single(data):
    """MD against the harmonic reference, single-species argon."""
    energies = data["single_energies"]
    weight = density_curve(energies, data["single_weight"])
    error = data["single_error"] / np.trapezoid(data["single_weight"], energies)
    modes = data["single_modes"]
    real = modes[modes > 1e-3]

    resolution = float(data["single_resolution"])
    reference = broadened(real, np.ones_like(real), energies, resolution)

    fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    ax.fill_between(
        energies,
        reference,
        color=REFERENCE,
        alpha=0.35,
        lw=0,
        label=(
            f"Euphonic, {len(real)} modes broadened to the "
            f"{HANN_MAIN_LOBE * resolution:.2f} meV instrument resolution"
        ),
    )
    ax.plot(energies, weight, color=ARGON, lw=1.4, label="mdins, MD velocity pDOS")
    ax.fill_between(
        energies, weight - error, weight + error, color=ARGON, alpha=0.3, lw=0
    )
    ax.plot(
        real,
        np.full_like(real, -0.012),
        "|",
        color=REFERENCE,
        ms=6,
        label="individual harmonic modes",
    )

    ax.set_xlim(0.0, 10.0)
    ax.set_ylim(-0.025, None)
    ax.set_xlabel("energy transfer (meV)")
    ax.set_ylabel("pDOS (1/meV, unit area)")
    ax.set_title(
        f"FCC argon, {np.prod(SUPERCELL)} atoms, {data['single_temperature']:.1f} K"
    )
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES / "single-species.png", dpi=200)
    plt.close(fig)


def figure_binary(data):
    """The mixed-species comparison: same forces, different masses."""
    energies = data["binary_energies"]
    resolution = float(data["binary_resolution"])

    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    means: list[tuple[str, float, str]] = []
    for species, colour in (("Ar", ARGON), ("Kr", KRYPTON)):
        weight = density_curve(energies, data[f"binary_{species}_weight"])
        error = data[f"binary_{species}_error"] / np.trapezoid(
            data[f"binary_{species}_weight"], energies
        )
        reference = broadened(
            data[f"binary_{species}_ref_energies"],
            data[f"binary_{species}_ref_weights"],
            energies,
            resolution,
        )
        ax.fill_between(energies, reference, color=colour, alpha=0.16, lw=0)
        ax.plot(energies, reference, color=colour, alpha=0.75, lw=1.1, ls="--")
        ax.plot(energies, weight, color=colour, lw=1.5, label=f"{species}, MD")
        ax.fill_between(
            energies, weight - error, weight + error, color=colour, alpha=0.3, lw=0
        )

        md_mean = moment(energies, data[f"binary_{species}_weight"], 1)
        ax.axvline(md_mean, color=colour, lw=0.9, ls=":")
        means.append((species, md_mean, colour))

    # Collected in a corner rather than beside each line: at these energies the labels
    # would otherwise sit on top of the peaks they describe.
    for row, (species, md_mean, colour) in enumerate(means):
        ax.annotate(
            f"{species} $\\langle E\\rangle$ = {md_mean:.2f} meV",
            (0.985, 0.95 - 0.07 * row),
            xycoords="axes fraction",
            ha="right",
            va="top",
            fontsize=8,
            color=colour,
        )

    ax.plot(
        [],
        [],
        color=REFERENCE,
        ls="--",
        lw=1.1,
        label="Euphonic, broadened to the same resolution",
    )
    ax.set_xlim(0.0, 10.0)
    ax.set_xlabel("energy transfer (meV)")
    ax.set_ylabel("pDOS (1/meV, unit area per species)")
    ax.set_title(
        f"Ordered Ar/Kr, {4 * np.prod(BINARY_SUPERCELL)} atoms, "
        f"{data['binary_temperature']:.1f} K — identical forces, mass ratio 2.10"
    )
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES / "mixed-species.png", dpi=200)
    plt.close(fig)


def figure_residuals(data):
    """Every comparison, in units of its own derived uncertainty."""
    rows = []

    energies = data["single_energies"]
    weight, error = data["single_weight"], data["single_error"]
    modes = data["single_modes"]
    real = modes[modes > 1e-3]
    for order in (1, 2):
        got = moment(energies, weight, order)
        want = float((real**order).mean())
        sigma = moment_uncertainty(energies, weight, error, order)
        rows.append((f"Ar only, moment {order}", (got - want) / sigma, ARGON))

    binary_energies = data["binary_energies"]
    top = max(data[f"binary_{s}_ref_energies"].max() for s in ("Ar", "Kr"))
    inside = binary_energies <= top + 5 * data["binary_resolution"]
    for species in ("Ar", "Kr"):
        weight = data[f"binary_{species}_weight"][inside]
        error = data[f"binary_{species}_error"][inside]
        for order in (1, 2):
            got = moment(binary_energies[inside], weight, order)
            want = moment(
                data[f"binary_{species}_ref_energies"],
                data[f"binary_{species}_ref_weights"],
                order,
            )
            sigma = moment_uncertainty(binary_energies[inside], weight, error, order)
            rows.append(
                (
                    f"Ar/Kr, {species} moment {order}",
                    (got - want) / sigma,
                    ARGON if species == "Ar" else KRYPTON,
                )
            )

    labels = [row[0] for row in rows][::-1]
    deviations = [row[1] for row in rows][::-1]
    colours = [row[2] for row in rows][::-1]

    fig, ax = plt.subplots(figsize=(7.0, 3.4), constrained_layout=True)
    ax.axvspan(-3, 3, color="0.9", label="$\\pm 3\\sigma$, the test threshold")
    ax.axvline(0, color="0.4", lw=0.8)
    ax.barh(labels, deviations, color=colours, height=0.55)
    ax.set_xlim(-4, 4)
    ax.set_xlabel("(MD $-$ Euphonic) / derived standard error")
    ax.set_title("Every M2 comparison, against its own propagated uncertainty")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.savefig(FIGURES / "residuals.png", dpi=200)
    plt.close(fig)


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
        data = dict(np.load(CACHE))
    else:
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            data = compute(Path(tmp))
        np.savez(CACHE, **data)

    figure_single(data)
    figure_binary(data)
    figure_residuals(data)
    print(f"wrote three figures to {FIGURES}")


if __name__ == "__main__":
    main()
