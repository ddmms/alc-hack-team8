# SPDX-License-Identifier: GPL-3.0-or-later
"""Generate the figures in ``docs/validation.md``.

Run from the repository root::

    uv run --group docs python docs/make_validation_figures.py

Imports the same harness the M2 tests assert on (``tests/lj_reference.py``), so the
figures cannot drift away from the numbers in the test suite. That means an ensemble of
independent NVE runs rather than a single one, because in NVE the normal-mode energies
are fixed by the initial velocity draw and the within-run Welch spread cannot see the
resulting run-to-run scatter. Roughly ten minutes, almost all of it molecular dynamics;
results are cached in ``docs/figures/validation-data.npz`` and reused unless
``--refresh`` is given.
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
    SEEDS,
    binary_cell,
    ensemble,
    ensemble_moment,
    equilibrium_lattice_constant,
    harmonic_modes,
    moment,
    reference_pdos,
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
    """Run both routes for both crystals, over the full ensemble of seeds.

    Every per-seed spectrum is kept rather than averaged here, so the figures can show
    both the ensemble mean and the scatter it was drawn from. The within-run Welch
    spread is kept alongside it for the single-species case only, because the point of
    one panel is to show how much narrower it is than the truth.
    """
    lattice_constant = equilibrium_lattice_constant()
    print(f"equilibrium lattice constant {lattice_constant:.6f} Å")

    single = single_species_cell(lattice_constant)
    single_modes = harmonic_modes(single, SUPERCELL, str(tmp / "lj"))
    print(f"single-species MD, {len(SEEDS)} seeds ...")
    single_ensemble = ensemble(single, SUPERCELL)

    binary_primitive = binary_cell(lattice_constant)
    binary_modes = harmonic_modes(
        binary_primitive, BINARY_SUPERCELL, str(tmp / "ljmix")
    )
    binary_reference = reference_pdos(binary_modes, binary_primitive)
    print(f"binary MD, {len(SEEDS)} seeds ...")
    binary_ensemble = ensemble(binary_primitive, BINARY_SUPERCELL)

    first = single_ensemble[0]
    data = {
        "lattice_constant": lattice_constant,
        "seeds": np.array(SEEDS),
        "single_energies": first.frequencies,
        # (n_seeds, n_bins): one normalised pDOS per independent run.
        "single_weights": np.array(
            [species_weight(density, "Ar")[0] for density in single_ensemble]
        ),
        # The within-run Welch spread of one run, for the comparison panel only.
        "single_welch_error": species_weight(first, "Ar")[1],
        "single_modes": single_modes.frequencies.to("meV").magnitude.ravel(),
        "single_temperatures": np.array(
            [density.temperature_md for density in single_ensemble]
        ),
        "single_resolution": first.metadata.energy_resolution,
        "binary_energies": binary_ensemble[0].frequencies,
        "binary_temperatures": np.array(
            [density.temperature_md for density in binary_ensemble]
        ),
        "binary_resolution": binary_ensemble[0].metadata.energy_resolution,
    }
    for species in ("Ar", "Kr"):
        frequencies, reference = binary_reference[species]
        data[f"binary_{species}_weights"] = np.array(
            [species_weight(density, species)[0] for density in binary_ensemble]
        )
        data[f"binary_{species}_ref_energies"] = frequencies
        data[f"binary_{species}_ref_weights"] = reference
        # Per-seed moments, restricted to the harmonic band, so the residual figure can
        # be drawn from the cache without re-running anything.
        top = max(binary_reference[s][0].max() for s in ("Ar", "Kr"))
        limit = top + 5 * data["binary_resolution"]
        for order in (1, 2):
            mean, error = ensemble_moment(
                binary_ensemble, order, species=species, limit=limit
            )
            data[f"binary_{species}_moment{order}"] = np.array([mean, error])
    for order in (1, 2):
        data[f"single_moment{order}"] = np.array(
            ensemble_moment(single_ensemble, order)
        )
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


def curves(energies, weights):
    """Per-seed pDOS as unit-area curves, and their mean."""
    normalised = np.array([density_curve(energies, weight) for weight in weights])
    return normalised, normalised.mean(axis=0)


def figure_single(data):
    """MD against the harmonic reference, single-species argon.

    Two bands are drawn deliberately. The narrow one is what a single run reports from
    its own Welch segments; the wide one is the actual spread over independent runs.
    The gap between them is the point of the figure.
    """
    energies = data["single_energies"]
    seed_curves, mean_curve = curves(energies, data["single_weights"])
    welch = data["single_welch_error"] / np.trapezoid(
        data["single_weights"][0], energies
    )
    spread = seed_curves.std(axis=0, ddof=1)
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
    ax.fill_between(
        energies,
        mean_curve - spread,
        mean_curve + spread,
        color=ARGON,
        alpha=0.18,
        lw=0,
        label=f"scatter over {len(seed_curves)} independent NVE runs",
    )
    ax.plot(
        energies,
        mean_curve,
        color=ARGON,
        lw=1.4,
        label="mdins, MD velocity pDOS (ensemble mean)",
    )
    ax.fill_between(
        energies,
        mean_curve - welch,
        mean_curve + welch,
        color=ARGON,
        alpha=0.45,
        lw=0,
        label="what one run reports from its own Welch segments",
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
    temperatures = data["single_temperatures"]
    ax.set_title(
        f"FCC argon, {np.prod(SUPERCELL)} atoms, "
        f"{temperatures.mean():.1f} ± {temperatures.std(ddof=1):.1f} K over "
        f"{len(temperatures)} runs"
    )
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES / "single-species.png", dpi=200)
    plt.close(fig)


def figure_binary(data):
    """The mixed-species comparison: same forces, different masses."""
    energies = data["binary_energies"]
    resolution = float(data["binary_resolution"])

    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    means: list[tuple[str, float, float, str]] = []
    for species, colour in (("Ar", ARGON), ("Kr", KRYPTON)):
        seed_curves, weight = curves(energies, data[f"binary_{species}_weights"])
        spread = seed_curves.std(axis=0, ddof=1)
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
            energies, weight - spread, weight + spread, color=colour, alpha=0.3, lw=0
        )

        md_mean, md_error = data[f"binary_{species}_moment{1}"]
        ax.axvline(md_mean, color=colour, lw=0.9, ls=":")
        means.append((species, float(md_mean), float(md_error), colour))

    # Collected in a corner rather than beside each line: at these energies the labels
    # would otherwise sit on top of the peaks they describe.
    for row, (species, md_mean, md_error, colour) in enumerate(means):
        ax.annotate(
            f"{species} $\\langle E\\rangle$ = {md_mean:.2f} $\\pm$ {md_error:.2f} meV",
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
    temperatures = data["binary_temperatures"]
    ax.set_title(
        f"Ordered Ar/Kr, {4 * np.prod(BINARY_SUPERCELL)} atoms, "
        f"{temperatures.mean():.1f} ± {temperatures.std(ddof=1):.1f} K over "
        f"{len(temperatures)} runs — identical forces, mass ratio 2.10"
    )
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURES / "mixed-species.png", dpi=200)
    plt.close(fig)


def figure_residuals(data):
    """Every comparison, in units of the ensemble standard error.

    The error on each bar is the scatter of that statistic over independent runs,
    divided by the root of their number. Nothing is propagated from within a run: the
    within-run Welch spread does not know about the initial velocity draw, which is the
    dominant term, and using it here would shrink every bar by roughly a factor of four.
    """
    rows = []

    modes = data["single_modes"]
    real = modes[modes > 1e-3]
    for order in (1, 2):
        got, sigma = data[f"single_moment{order}"]
        want = float((real**order).mean())
        rows.append((f"Ar only, moment {order}", (got - want) / sigma, ARGON))

    for species in ("Ar", "Kr"):
        for order in (1, 2):
            got, sigma = data[f"binary_{species}_moment{order}"]
            want = moment(
                data[f"binary_{species}_ref_energies"],
                data[f"binary_{species}_ref_weights"],
                order,
            )
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
    ax.set_xlabel("(MD $-$ Euphonic) / ensemble standard error")
    ax.set_title(
        f"Every M2 comparison, against the scatter over {len(data['seeds'])} NVE runs"
    )
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.savefig(FIGURES / "residuals.png", dpi=200)
    plt.close(fig)


def report(data) -> None:
    """Print the results table in the form ``validation.md`` carries it.

    Printed rather than written so that updating the page is a deliberate act: the
    numbers move a little with every ensemble, and a table that rewrote itself would
    make it easy to change a documented result without noticing.
    """
    modes = data["single_modes"]
    real = modes[modes > 1e-3]
    print("\n| Comparison | MD | Euphonic | Difference | 1σ | Deviation |")  # noqa: RUF001
    print("|---|---|---|---|---|---|")
    rows = [
        (
            f"Ar only, ⟨E{'²' if order == 2 else ''}⟩",
            f"single_moment{order}",
            float((real**order).mean()),
        )
        for order in (1, 2)
    ]
    for species in ("Ar", "Kr"):
        for order in (1, 2):
            rows.append(
                (
                    f"Ar/Kr, {species} ⟨E{'²' if order == 2 else ''}⟩",
                    f"binary_{species}_moment{order}",
                    moment(
                        data[f"binary_{species}_ref_energies"],
                        data[f"binary_{species}_ref_weights"],
                        order,
                    ),
                )
            )
    for label, key, want in rows:
        got, sigma = data[key]
        print(
            f"| {label} | {got:.4f} ± {sigma:.4f} | {want:.4f} | "
            f"{100 * (got - want) / want:+.2f}% | {100 * sigma / got:.2f}% | "
            f"{(got - want) / sigma:+.1f}σ |"  # noqa: RUF001
        )
    temperatures = data["single_temperatures"]
    print(
        f"\nsingle-species T {temperatures.mean():.2f} ± "
        f"{temperatures.std(ddof=1):.2f} K over {len(temperatures)} runs "
        f"(range {temperatures.min():.2f}-{temperatures.max():.2f})"
    )


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
    report(data)


if __name__ == "__main__":
    main()
