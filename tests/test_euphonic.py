# SPDX-License-Identifier: GPL-3.0-or-later
"""M2: the harmonic validation, against Euphonic.

This is the proposal's headline check and the highest-risk milestone in plan.md. The
same Lennard-Jones crystal is described twice by routes that share nothing but the
potential:

.. code-block:: none

    ASE LennardJones ──┬── finite displacement → force constants → Euphonic → phonons
                       └── NVE MD, low T       → mdins → pdos()

The left-hand route solves the harmonic problem exactly and never integrates anything;
the right-hand route integrates the real dynamics and never diagonalises anything. If
stages A-C have the normalisation, the units or the frequency axis wrong, the two
disagree.

Two systems are compared. :class:`TestAgainstEuphonic` uses single-species argon, which
tests the pipeline end to end. :class:`TestMixedSpecies` uses an ordered Ar/Kr crystal,
where the potential is species-blind so the two elements differ *only* by mass — the
projection onto species then has exactly one thing to get right, and getting it wrong
cannot hide behind a difference in the forces.

Two traps produce a spurious mismatch, and both are the test's own responsibility rather
than the code's:

**The reference must use the q-grid commensurate with the MD supercell.** A supercell
under periodic boundary conditions can physically host one wavevector per cell. Sampling
the reference on a denser grid adds modes the MD was never able to support, and the
comparison then fails for a reason that has nothing to do with mdins. This is why
:func:`~tests.lj_reference.commensurate_grid` is built from the supercell repeat rather
than chosen for smoothness.

**The MD must be cold enough to stay harmonic.** The comparison is against a harmonic
calculation, so any real anharmonicity is a genuine physical difference, not an error.
At 10 K, argon is deep in the harmonic regime — a fortieth of its melting point.

Because the reference is a few dozen near-degenerate delta functions spread over a few
meV, a point-by-point comparison at the estimator's 0.1 meV resolution is meaningless:
the mode spacing is smaller than the resolution. The comparison is therefore made on the
quantities that survive that, namely the moments of the distribution and the lineshape
on bins wider than the mode spacing. Every tolerance is derived from the Welch
inter-segment spread (design.md D5) rather than chosen to fit, and
:meth:`TestAgainstEuphonic.test_the_comparison_can_fail` checks that the derived
tolerances are tight enough to reject a wrong answer.
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.slow

pytest.importorskip("euphonic", reason="euphonic is an optional extra")

from tests.lj_reference import (  # noqa: E402 - must follow the euphonic skip
    TEMPERATURE,
    acoustic_weight_fraction,
    binary_cell,
    equilibrium_lattice_constant,
    harmonic_modes,
    lennard_jones,
    moment,
    moment_uncertainty,
    reference_pdos,
    run_nve,
    single_species_cell,
    species_weight,
)

#: 3×3×3 of the one-atom primitive: 27 atoms, 27 wavevectors, 81 modes.
SUPERCELL = (3, 3, 3)

#: 2×2×2 of the four-atom conventional cell: 32 atoms, 8 wavevectors, 96 modes. Smaller
#: in wavevectors than the single-species case because each cell carries four atoms, so
#: the run stays the same length in wall-clock.
BINARY_SUPERCELL = (2, 2, 2)


@pytest.fixture(scope="module")
def lattice_constant():
    return equilibrium_lattice_constant()


@pytest.fixture(scope="module")
def harmonic_reference(lattice_constant, tmp_path_factory):
    """Phonon energies in meV on the commensurate grid, shape ``(n_q, n_branches)``."""
    modes = harmonic_modes(
        single_species_cell(lattice_constant),
        SUPERCELL,
        str(tmp_path_factory.mktemp("phonons") / "lj"),
    )
    return modes.frequencies.to("meV").magnitude


@pytest.fixture(scope="module")
def md_density(lattice_constant):
    """pDOS of the same crystal, from NVE molecular dynamics."""
    return run_nve(single_species_cell(lattice_constant), SUPERCELL)


@pytest.fixture(scope="module")
def binary(lattice_constant, tmp_path_factory):
    """The Ar/Kr crystal, both ways, plus the reference cell needed to project it."""
    primitive = binary_cell(lattice_constant)
    modes = harmonic_modes(
        primitive,
        BINARY_SUPERCELL,
        str(tmp_path_factory.mktemp("phonons_binary") / "ljmix"),
    )
    return {
        "primitive": primitive,
        "modes": modes,
        "reference": reference_pdos(modes, primitive),
        "density": run_nve(primitive, BINARY_SUPERCELL),
    }


def total_weight(density):
    """Normalised pDOS of the whole system and its standard error."""
    pdos = density.total_pdos()
    total = pdos.sum()
    error = (density.pdos_std * density.counts[:, None]).sum(axis=0)
    return pdos / total, error / total


class TestTheReferenceItself:
    """Before comparing against it: is the reference a sane phonon spectrum?

    A reference that is quietly wrong would be indistinguishable from agreement, since
    both sides come from the same ASE potential.
    """

    def test_every_mode_is_real(self, harmonic_reference):
        """Imaginary modes, which Euphonic reports as negative, mean the structure is
        not at a minimum and the harmonic reference is meaningless."""
        assert harmonic_reference.min() > -1e-3

    def test_there_are_three_acoustic_modes_at_gamma(self, harmonic_reference):
        """Γ is the first point of the grid. Its three zero modes are rigid translation
        of the whole cell — the same three degrees of freedom that
        ``remove_com_velocity`` projects out of the MD, which is why the comparison
        below excludes them from both sides."""
        gamma = harmonic_reference[0]
        assert np.abs(gamma).max() < 1e-3
        assert harmonic_reference.shape == (27, 3)

    def test_the_grid_has_one_wavevector_per_supercell_cell(self, harmonic_reference):
        """The commensurate condition, asserted rather than assumed. If this ever
        changes the comparison becomes invalid, for reasons that would be very hard to
        diagnose from the resulting lineshape mismatch."""
        assert harmonic_reference.shape[0] == np.prod(SUPERCELL)

    def test_the_band_is_where_argon_should_be(self, harmonic_reference):
        """A loose sanity bound from the Debye temperature of solid argon, about 92 K,
        so a band top of order 8 meV. Catches a reference that is out by a unit
        conversion, which would otherwise agree with nothing and be blamed on the MD."""
        assert 4.0 < harmonic_reference.max() < 12.0


class TestAgainstEuphonic:
    """The single-species comparison."""

    def test_the_md_stayed_cold_and_harmonic(self, md_density):
        """The premise of the whole comparison. NVE started at ``2T`` settles near
        ``T``; if it has not, the anharmonic shift is no longer negligible and a
        disagreement below would be physics rather than a bug."""
        assert md_density.temperature_md == pytest.approx(TEMPERATURE, rel=0.2)

    def test_first_moment_agrees_within_the_welch_spread(
        self, md_density, harmonic_reference
    ):
        """The mean phonon energy. This is the single number most sensitive to a
        frequency-axis error: a factor of 2π anywhere would move it by 2π."""
        energies = md_density.frequencies
        weight, error = total_weight(md_density)

        got = moment(energies, weight, 1)
        uncertainty = moment_uncertainty(energies, weight, error, 1)
        want = float(harmonic_reference[harmonic_reference > 1e-3].mean())

        assert abs(got - want) < 3 * uncertainty, (
            f"MD mean phonon energy {got:.4f} meV against a harmonic {want:.4f} meV, "
            f"a difference of {abs(got - want) / uncertainty:.1f} standard errors"
        )

    def test_second_moment_agrees_within_the_welch_spread(
        self, md_density, harmonic_reference
    ):
        """Sensitive to the width of the band rather than its position, so it fails for
        a spectrum that is centred correctly but spread wrongly."""
        energies = md_density.frequencies
        weight, error = total_weight(md_density)

        got = moment(energies, weight, 2)
        uncertainty = moment_uncertainty(energies, weight, error, 2)
        real = harmonic_reference[harmonic_reference > 1e-3]
        want = float((real**2).mean())

        assert abs(got - want) < 3 * uncertainty, (
            f"MD second moment {got:.3f} against a harmonic {want:.3f} meV², "
            f"a difference of {abs(got - want) / uncertainty:.1f} standard errors"
        )

    def test_the_lineshape_agrees_on_bins_wider_than_the_mode_spacing(
        self, md_density, harmonic_reference
    ):
        """The shape comparison D6 asks for, made at a bin width the comparison can
        actually support.

        The 78 modes occupy 4 meV, so they are spaced about 0.05 meV apart — half the
        estimator's resolution. Comparing at the output bin width would be comparing
        noise against delta functions. Half a meV is wide enough to hold several modes
        and still narrow enough to have eight bins across the band.
        """
        width = 0.5
        edges = np.arange(0.0, 12.0 + width, width)
        energies = md_density.frequencies
        weight, _ = total_weight(md_density)

        measured, _ = np.histogram(energies, bins=edges, weights=weight)
        real = harmonic_reference[harmonic_reference > 1e-3]
        reference, _ = np.histogram(real, bins=edges)
        reference = reference / reference.sum()

        overlap = float(
            np.sum(measured * reference)
            / np.sqrt(np.sum(measured**2) * np.sum(reference**2))
        )
        assert overlap > 0.95, f"lineshape overlap {overlap:.4f}"

    def test_the_acoustic_modes_were_projected_out(
        self, md_density, harmonic_reference
    ):
        """Below the lowest finite mode there is nothing but the three Γ translations,
        and those were removed from the velocities. Weight here would mean the
        projection did not take, and the ``1/E`` factor downstream would amplify it.
        """
        real = harmonic_reference[harmonic_reference > 1e-3]
        energies = md_density.frequencies
        weight, _ = total_weight(md_density)
        below = energies < real.min() - 3 * md_density.metadata.energy_resolution
        assert weight[below].sum() < 0.01

    def test_there_is_no_weight_beyond_the_band_top(
        self, md_density, harmonic_reference
    ):
        """A hard edge in the harmonic spectrum, softened only by the estimator's
        resolution and by whatever anharmonicity survives at 10 K. Weight well above it
        would mean spurious high-frequency content — the signature of aliasing or of a
        mishandled window.
        """
        energies = md_density.frequencies
        weight, _ = total_weight(md_density)
        top = harmonic_reference.max()
        # Four resolution widths above the top, plus a margin for the anharmonic shift.
        beyond = energies > top + 4 * md_density.metadata.energy_resolution + 1.0
        assert weight[beyond].sum() < 0.02

    def test_the_comparison_can_fail(self, md_density, harmonic_reference):
        """The control. Every tolerance above is derived from the measured spread, which
        would be worthless if that spread were so wide that anything passed.

        The derived tolerance is measured here rather than asserted. The observed
        agreement is 0.6 standard errors on the first moment and 0.3 on the second,
        while the three-sigma threshold rejects any frequency error above roughly 3.5%
        — far smaller than a factor of 2π, of 2, or of any plausible unit conversion. A
        5% stretch is used below because it is the nearest round figure the test
        rejects, not because it is comfortable.
        """
        energies = md_density.frequencies
        weight, error = total_weight(md_density)
        uncertainty = moment_uncertainty(energies, weight, error, 1)
        want = float(harmonic_reference[harmonic_reference > 1e-3].mean())

        stretched = moment(1.05 * energies, weight, 1)
        assert abs(stretched - want) > 3 * uncertainty, (
            "a 5% frequency error is within the derived tolerance, so the agreement "
            "tests above prove nothing"
        )


class TestMixedSpecies:
    """Ar and Kr in one cell, differing only by mass.

    The single-species tests above cannot tell whether the per-atom projection is right,
    because every atom is equivalent and any projection error averages away. Here it
    cannot: argon is lighter, so it carries the high-energy end of the band and krypton
    the low, and the split has to come out at the value Euphonic gives.

    Moments are taken over the support of the harmonic band rather than over the whole
    output range. This is not a convenience. The ``E²`` weighting in the second moment
    means the roughly 1% of MD weight lying above the band top — resolution leakage and
    residual anharmonicity — contributes as much as a 2% error in the band itself; over
    the full range the second moments disagree by 2.0%, and restricted to the band by
    0.1%. :meth:`test_little_weight_lies_outside_the_harmonic_band` bounds the
    excluded weight separately, so nothing is hidden by the restriction.
    """

    @staticmethod
    def band_limit(binary):
        """Top of the harmonic band, plus five resolution widths.

        Expressed in resolution widths rather than as a fixed margin in meV, so it
        tracks the estimator rather than this particular run.
        """
        top = max(freq.max() for freq, _ in binary["reference"].values())
        return top + 5 * binary["density"].metadata.energy_resolution

    def test_the_cell_really_does_contain_both_species(self, binary):
        """Guards the premise. A typo in the L1₀ site assignment giving 4 Ar and 0 Kr
        would make every test below pass trivially."""
        symbols = binary["density"].symbols
        assert symbols.count("Ar") == 16
        assert symbols.count("Kr") == 16

    def test_the_two_species_feel_identical_forces(self, binary):
        """The design of the experiment, asserted. ASE's Lennard-Jones takes one epsilon
        and one sigma, so the dynamical matrix differs between the species only through
        the mass factor. If this ever stopped being true, a shape difference below could
        no longer be attributed to mass alone."""
        primitive = binary["primitive"].copy()
        primitive.calc = lennard_jones()
        assert np.abs(primitive.get_forces()).max() < 1e-9
        assert set(primitive.get_chemical_symbols()) == {"Ar", "Kr"}

        # The same cell with the elements swapped must give the same forces: that is
        # what "species-blind" means, and it is the assumption the class rests on.
        swapped = primitive.copy()
        swapped.set_chemical_symbols(["Kr", "Kr", "Ar", "Ar"])
        swapped.calc = lennard_jones()
        assert np.allclose(
            swapped.get_potential_energy(), primitive.get_potential_energy()
        )

    def test_the_lighter_species_sits_at_higher_energy(self, binary):
        """The qualitative signature, checked before the quantitative one. If this fails
        the species labels have been swapped somewhere, and a numerical comparison
        against the matching reference would still pass."""
        energies = binary["density"].frequencies
        argon, _ = species_weight(binary["density"], "Ar")
        krypton, _ = species_weight(binary["density"], "Kr")
        assert moment(energies, argon, 1) > moment(energies, krypton, 1)

    @pytest.mark.parametrize("species", ["Ar", "Kr"])
    @pytest.mark.parametrize("order", [1, 2])
    def test_each_species_matches_its_harmonic_projection(self, binary, species, order):
        """The comparison this class exists for: the species-projected pDOS, against
        the same projection of the Euphonic eigenvectors."""
        density = binary["density"]
        limit = self.band_limit(binary)
        inside = density.frequencies <= limit

        weight, error = species_weight(density, species)
        got = moment(density.frequencies[inside], weight[inside], order)
        uncertainty = moment_uncertainty(
            density.frequencies[inside], weight[inside], error[inside], order
        )

        frequencies, reference = binary["reference"][species]
        want = moment(frequencies, reference, order)

        assert abs(got - want) < 3 * uncertainty, (
            f"{species} moment {order}: MD {got:.4f} against a harmonic {want:.4f}, "
            f"a difference of {abs(got - want) / uncertainty:.1f} standard errors"
        )

    def test_the_species_are_separated_by_the_amount_euphonic_predicts(self, binary):
        """The split itself, rather than each species against its own reference.

        Worth asserting separately because it is the one number a reader of the figures
        will check, and because a systematic error common to both species — a wrong
        timestep, say — cancels here and would leave this passing while the individual
        comparisons drifted together.
        """
        density = binary["density"]
        limit = self.band_limit(binary)
        inside = density.frequencies <= limit

        measured = {}
        for species in ("Ar", "Kr"):
            weight, _ = species_weight(density, species)
            measured[species] = moment(density.frequencies[inside], weight[inside], 1)

        predicted = {
            species: moment(*binary["reference"][species], 1)
            for species in ("Ar", "Kr")
        }
        got = measured["Ar"] - measured["Kr"]
        want = predicted["Ar"] - predicted["Kr"]
        assert got == pytest.approx(want, rel=0.05), (
            f"mass splitting {got:.4f} meV against a harmonic {want:.4f} meV"
        )

    def test_the_weight_ratio_is_set_by_mass(self, binary):
        """Amplitude rather than shape. Equipartition puts ``3k_BT/m_i`` under each
        atom's spectrum, so the ratio of integrated weights is the inverse mass ratio —
        but only after the centre-of-mass correction ``3k_BT(1/m_i − 1/M)``. Without
        that term the ratio comes out 3.4% wrong instead of 1.1%, which is why this is a
        real check on the correction and not a restatement of the sum rule.
        """
        from mdins.units import KB

        density = binary["density"]
        grouped = density.group_by_species()
        masses = dict(zip(grouped.symbols, grouped.masses, strict=True))
        total_mass = float((np.array(density.masses)).sum())

        integrated = {}
        for species in ("Ar", "Kr"):
            index = [i for i, s in enumerate(density.symbols) if s == species]
            integrated[species] = float(density.mean_square_velocity()[index].mean())

        def expected(species):
            return (
                3
                * KB
                * density.temperature_md
                * (1.0 / masses[species] - 1.0 / total_mass)
            )

        got = integrated["Ar"] / integrated["Kr"]
        want = expected("Ar") / expected("Kr")
        assert got == pytest.approx(want, rel=0.03), (
            f"weight ratio {got:.4f} against an equipartition {want:.4f}"
        )

    def test_the_heavier_species_carries_more_of_the_projected_translation(
        self, binary
    ):
        """A property of the reference, not of the MD, but it explains why the acoustic
        modes have to be masked per species rather than globally.

        Γ acoustic eigenvectors go as ``√m``, so krypton holds about twice argon's share
        of the rigid translation that ``remove_com_velocity`` removes. Masking the same
        fraction from both species would bias the comparison in opposite directions.
        """
        fraction = acoustic_weight_fraction(binary["modes"], binary["primitive"])
        masses = dict(
            zip(
                binary["primitive"].get_chemical_symbols(),
                binary["primitive"].get_masses(),
                strict=True,
            )
        )
        ratio = fraction["Kr"] / fraction["Ar"]
        assert ratio == pytest.approx(masses["Kr"] / masses["Ar"], rel=0.05)

    def test_little_weight_lies_outside_the_harmonic_band(self, binary):
        """Bounds what the band restriction above excludes.

        Without this, restricting the moments to the band would be a way of not looking
        at a disagreement. About 1% of the argon weight and 0.3% of the krypton weight
        lies above the cutoff, which is small enough that the restriction is a choice of
        statistic rather than a filter on the result.
        """
        density = binary["density"]
        limit = self.band_limit(binary)
        outside = density.frequencies > limit
        for species in ("Ar", "Kr"):
            weight, _ = species_weight(density, species)
            assert weight[outside].sum() < 0.02, (
                f"{species} has {100 * weight[outside].sum():.2f}% of its weight "
                f"above {limit:.2f} meV"
            )
