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
on bins wider than the mode spacing.

**Where the tolerances come from.** Every comparison of a *moment* is made against the
run-to-run standard error over :data:`~tests.lj_reference.SEEDS` independent
trajectories, at three sigma. That is a measured sampling distribution rather than a
propagated one, which matters here more than it usually does: these runs are NVE on a
harmonic crystal, so the normal-mode energies are constants of motion fixed by the
initial velocity draw, and the inter-segment spread that the estimator reports is
structurally blind to them. It comes out about four times too small — small enough that
the comparisons passed on one seed and failed on most others. The long comment in
:mod:`tests.lj_reference` sets this out in full.

Not every tolerance is derived, and the ones that are not should not be read as if they
were. The lineshape overlap, the out-of-band weight bounds and the band limit are
judgement calls, chosen to be loose enough not to fail on noise; they bound the size of
an effect rather than measure it.
:meth:`TestAgainstEuphonic.test_the_comparison_can_fail` records what the derived
tolerances can actually discriminate, which is a frequency-axis error of about 3.5% —
not the 1% the within-run spread would have implied.
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.slow

pytest.importorskip("euphonic", reason="euphonic is an optional extra")

from tests.lj_reference import (  # noqa: E402 - must follow the euphonic skip
    SEEDS,
    SEGMENT_LENGTH,
    TEMPERATURE,
    acoustic_weight_fraction,
    binary_cell,
    ensemble,
    ensemble_moment,
    equilibrium_lattice_constant,
    estimate,
    harmonic_modes,
    lennard_jones,
    md_trajectory,
    moment,
    reference_pdos,
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
def md_ensemble(lattice_constant):
    """pDOS of the same crystal, from one NVE run per seed in :data:`SEEDS`.

    An ensemble rather than a single run because a single NVE run cannot measure its own
    uncertainty here — see the long comment in :mod:`tests.lj_reference`. This is the
    expensive fixture in the module; everything that only needs *a* trajectory takes
    element zero rather than paying for a second one.
    """
    return ensemble(single_species_cell(lattice_constant), SUPERCELL)


@pytest.fixture(scope="module")
def md_density(md_ensemble):
    """One representative run, for tests about shape rather than about scatter."""
    return md_ensemble[0]


@pytest.fixture(scope="module")
def binary(lattice_constant, tmp_path_factory):
    """The Ar/Kr crystal, both ways, plus the reference cell needed to project it."""
    primitive = binary_cell(lattice_constant)
    modes = harmonic_modes(
        primitive,
        BINARY_SUPERCELL,
        str(tmp_path_factory.mktemp("phonons_binary") / "ljmix"),
    )
    runs = ensemble(primitive, BINARY_SUPERCELL)
    return {
        "primitive": primitive,
        "modes": modes,
        "reference": reference_pdos(modes, primitive),
        "ensemble": runs,
        "density": runs[0],
    }


@pytest.fixture(scope="module")
def md_velocities(lattice_constant):
    """One trajectory, kept so that both estimators can be run over it."""
    return md_trajectory(single_species_cell(lattice_constant), SUPERCELL)


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

    def test_the_md_stayed_cold_and_harmonic(self, md_ensemble):
        """The premise of the whole comparison. NVE started at ``2T`` settles near
        ``T``; if it has not, the anharmonic shift is no longer negligible and a
        disagreement below would be physics rather than a bug.

        Asserted on the ensemble mean, not on a single run. An individual run lands
        anywhere in roughly 7-13 K, because the temperature is set by the initial
        Maxwell-Boltzmann draw over 27 atoms and then conserved; that scatter is the
        expected behaviour of a small NVE system, not a failure to equilibrate. What
        would be a real failure is the mean sitting away from ``T``.
        """
        temperatures = np.array([density.temperature_md for density in md_ensemble])
        error = temperatures.std(ddof=1) / np.sqrt(len(temperatures))
        assert abs(temperatures.mean() - TEMPERATURE) < 3 * error, (
            f"ensemble mean temperature {temperatures.mean():.2f} K against a target "
            f"{TEMPERATURE} K, with a standard error of {error:.2f} K"
        )

    @pytest.mark.parametrize("order", [1, 2])
    def test_the_moments_agree_with_the_harmonic_reference(
        self, md_ensemble, harmonic_reference, order
    ):
        """The headline comparison.

        The first moment is the number most sensitive to a frequency-axis error: a
        factor of 2π anywhere moves it by 2π. The second is sensitive to the width of
        the band rather than its position, so it catches a spectrum that is centred
        correctly but spread wrongly.

        Both are compared as ensemble means against the run-to-run standard error. The
        tolerance is therefore the measured sampling distribution of the statistic, with
        no assumption about how bins or atoms are correlated.
        """
        got, uncertainty = ensemble_moment(md_ensemble, order)
        real = harmonic_reference[harmonic_reference > 1e-3]
        want = float((real**order).mean())

        assert abs(got - want) < 3 * uncertainty, (
            f"MD moment {order} of {got:.4f} ± {uncertainty:.4f} against a "
            f"harmonic {want:.4f}, a difference of "
            f"{abs(got - want) / uncertainty:.1f} "
            f"standard errors over {len(SEEDS)} seeds"
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

    def test_the_comparison_can_fail(self, md_ensemble, harmonic_reference):
        """The control. Every tolerance above is derived from the measured spread, which
        would be worthless if that spread were so wide that anything passed.

        The discriminating power is a property of the system and the ensemble size, not
        something to be asserted into existence: with ten seeds the standard error on
        the first moment is about 1.2%, so the three-sigma threshold rejects a
        frequency-axis error above roughly 3.5%. That is the honest figure. It is far
        smaller than a factor of 2π, of 2, or of any plausible unit conversion, and it
        is what makes the agreement above meaningful — but it is also much weaker than
        the 1% that the *within-run* Welch spread would have suggested, which is exactly
        the overclaim this ensemble replaced.

        A 5% stretch is used below: the nearest round figure comfortably outside the
        measured power, chosen so that the control does not sit on its own threshold.
        """
        got, uncertainty = ensemble_moment(md_ensemble, 1)
        want = float(harmonic_reference[harmonic_reference > 1e-3].mean())

        stretched = 1.05 * got
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
        no longer be attributed to mass alone.

        Getting this test right is harder than it looks, and two natural forms of it are
        worthless. Asserting that the forces vanish proves nothing: every L1₀ site is an
        inversion centre of the underlying FCC lattice, so the forces are zero by
        symmetry for *any* pair potential at *any* lattice constant. Nor does swapping
        the two labels: Ar and Kr sit on sublattices related by the translation
        ``(0, ½, ½)``, so the swapped cell is the original one shifted, and its energy
        is unchanged for a species-*aware* potential too.

        What does discriminate is relabelling every atom as a single element, which is
        not a symmetry of anything, evaluated away from the equilibrium geometry so that
        the forces are not fixed by the lattice. A potential with per-species parameters
        fails both assertions below.
        """
        primitive = binary["primitive"].copy()
        assert set(primitive.get_chemical_symbols()) == {"Ar", "Kr"}

        rng = np.random.default_rng(0)
        displaced = primitive.copy()
        displaced.positions += rng.normal(scale=0.05, size=displaced.positions.shape)

        mixed = displaced.copy()
        mixed.calc = lennard_jones()
        relabelled = displaced.copy()
        relabelled.set_chemical_symbols(["Ar"] * len(relabelled))
        relabelled.calc = lennard_jones()

        assert mixed.get_potential_energy() == pytest.approx(
            relabelled.get_potential_energy(), abs=1e-12
        )
        assert np.allclose(mixed.get_forces(), relabelled.get_forces(), atol=1e-12)

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
        limit = self.band_limit(binary)
        got, uncertainty = ensemble_moment(
            binary["ensemble"], order, species=species, limit=limit
        )
        want = moment(*binary["reference"][species], order)

        assert abs(got - want) < 3 * uncertainty, (
            f"{species} moment {order}: MD {got:.4f} ± {uncertainty:.4f} "
            f"against a harmonic {want:.4f}, a difference of "
            f"{abs(got - want) / uncertainty:.1f} standard errors"
        )

    def test_the_species_are_separated_by_the_amount_euphonic_predicts(self, binary):
        """The split itself, rather than each species against its own reference.

        Worth asserting separately because it is the one number a reader of the figures
        will check, and because a systematic error common to both species — a wrong
        timestep, say — cancels here and would leave this passing while the individual
        comparisons drifted together.

        The splitting is a difference of two numbers that are each uncertain and, it
        turns out, anti-correlated: a draw that puts extra energy into the high modes
        raises the argon mean and lowers the krypton one. That makes it the noisiest
        statistic in the module, and the one where a fixed percentage tolerance was
        least defensible. It is formed per run and compared through the spread of those,
        so the anti-correlation is measured rather than assumed away.
        """
        limit = self.band_limit(binary)
        splitting = []
        for density in binary["ensemble"]:
            inside = density.frequencies <= limit
            means = {
                species: moment(
                    density.frequencies[inside],
                    species_weight(density, species)[0][inside],
                    1,
                )
                for species in ("Ar", "Kr")
            }
            splitting.append(means["Ar"] - means["Kr"])

        sample = np.array(splitting)
        got = float(sample.mean())
        uncertainty = float(sample.std(ddof=1) / np.sqrt(len(sample)))
        want = moment(*binary["reference"]["Ar"], 1) - moment(
            *binary["reference"]["Kr"], 1
        )

        assert abs(got - want) < 3 * uncertainty, (
            f"mass splitting {got:.4f} ± {uncertainty:.4f} meV against a harmonic "
            f"{want:.4f} meV, a difference of {abs(got - want) / uncertainty:.1f} "
            f"standard errors"
        )

    def test_the_weight_ratio_is_set_by_mass(self, binary):
        """Amplitude rather than shape. Equipartition puts ``3k_BT/m_i`` under each
        atom's spectrum, so the ratio of integrated weights is the inverse mass ratio —
        but only after the centre-of-mass correction ``3k_BT(1/m_i − 1/M)``.

        The correction is worth about 2.3% in this ratio, and a single run scatters by
        rather more than that, so on one trajectory this test cannot honestly claim to
        be checking the correction at all: an earlier version asserted a fixed 3%
        tolerance, which happened to pass on one seed and failed on most. Over the
        ensemble the mean is precise enough for the corrected and uncorrected
        expectations to be distinguished, and both are checked below.
        """
        ratios = []
        for density in binary["ensemble"]:
            velocity = density.mean_square_velocity()
            mean = {
                species: float(
                    velocity[
                        [i for i, s in enumerate(density.symbols) if s == species]
                    ].mean()
                )
                for species in ("Ar", "Kr")
            }
            ratios.append(mean["Ar"] / mean["Kr"])

        sample = np.array(ratios)
        got = float(sample.mean())
        uncertainty = float(sample.std(ddof=1) / np.sqrt(len(sample)))

        masses = dict(
            zip(binary["density"].symbols, binary["density"].masses, strict=True)
        )
        total_mass = float(np.sum(binary["density"].masses))
        corrected = (1.0 / masses["Ar"] - 1.0 / total_mass) / (
            1.0 / masses["Kr"] - 1.0 / total_mass
        )
        uncorrected = masses["Kr"] / masses["Ar"]

        assert abs(got - corrected) < 3 * uncertainty, (
            f"weight ratio {got:.4f} ± {uncertainty:.4f} against an equipartition "
            f"{corrected:.4f}, a difference of "
            f"{abs(got - corrected) / uncertainty:.1f} standard errors"
        )
        # And the correction is what makes that work: the uncorrected ratio is further
        # out than the corrected one, which is the whole claim of this test. Asserted as
        # a comparison rather than against a fixed percentage, because the 2.3% size of
        # the correction is smaller than a single run's scatter and only the ensemble
        # can see it at all.
        assert abs(got - corrected) < abs(got - uncorrected), (
            f"the centre-of-mass correction moves the expectation from "
            f"{uncorrected:.4f} to {corrected:.4f}, and the measured {got:.4f} is not "
            f"closer to the "
            "corrected value, so this test does not show what it claims to"
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
        limit = self.band_limit(binary)
        for species in ("Ar", "Kr"):
            for density in binary["ensemble"]:
                outside = density.frequencies > limit
                weight, _ = species_weight(density, species)
                assert weight[outside].sum() < 0.02, (
                    f"{species} has {100 * weight[outside].sum():.2f}% of its weight "
                    f"above {limit:.2f} meV"
                )


class TestTheTwoEstimatorsAgree:
    """Welch against VACF, on molecular dynamics rather than on a synthetic signal.

    design.md D5 chooses Welch as the default and keeps the windowed velocity
    autocorrelation as a second route, and plan.md §5 lists their agreement as an
    invariant. :class:`tests.test_spectral.TestEstimatorAgreement` already checks it on
    pure tones and on white noise, where the answer is known analytically.

    This is the harder version of the same check, and the reason it belongs in M2 rather
    than beside the others. A real phonon band is neither a tone nor noise: it is a few
    dozen closely spaced lines with a hard edge at the band top, which is exactly the
    structure that makes the two routes diverge. Welch's leakage is set by the segment
    window; the VACF route's is set by truncating the correlation at a maximum lag and
    tapering it, a completely different approximation. Agreement here means the answer
    does not depend on which of those was chosen.

    Both estimators see the *same* velocities, so the run-to-run scatter that dominates
    everything else in this module cancels exactly and the comparison is sharp. That is
    also why this class does not use the ensemble.
    """

    @pytest.fixture(scope="class")
    def welch(self, md_velocities):
        return estimate(md_velocities)

    @pytest.fixture(scope="class")
    def vacf(self, md_velocities):
        return estimate(md_velocities, estimator="vacf", max_lag=SEGMENT_LENGTH)

    @pytest.mark.parametrize("order", [1, 2])
    def test_the_moments_agree_far_inside_the_ensemble_scatter(
        self, welch, vacf, md_ensemble, order
    ):
        """The tolerance is the same one the harmonic comparison uses, which makes the
        two numbers directly comparable.

        The point is not merely that the estimators agree, but that they agree by a wide
        margin on the scale that matters. If the choice of estimator moved a moment by
        anything approaching the run-to-run spread, every tolerance in this module would
        have to carry a second term for it.
        """
        _, scatter = ensemble_moment(md_ensemble, order)

        got = moment(welch.frequencies, welch.total_pdos(), order)
        other = moment(vacf.frequencies, vacf.total_pdos(), order)

        assert abs(got - other) < scatter, (
            f"Welch moment {order} of {got:.4f} against a VACF {other:.4f}: the "
            f"estimators differ by {abs(got - other) / scatter:.2f} of the run-to-run "
            "standard error, so the choice of estimator is not negligible"
        )

    def test_both_recover_the_harmonic_reference(
        self, welch, vacf, harmonic_reference, md_ensemble
    ):
        """Agreement between two estimators is not correctness — they could be wrong
        together, and on a shared trajectory a normalisation error would be shared too.
        Anchoring both against Euphonic is what rules that out."""
        real = harmonic_reference[harmonic_reference > 1e-3]
        want = float(real.mean())
        _, scatter = ensemble_moment(md_ensemble, 1)

        for name, density in (("welch", welch), ("vacf", vacf)):
            got = moment(density.frequencies, density.total_pdos(), 1)
            assert abs(got - want) < 3 * scatter, (
                f"{name} mean phonon energy {got:.4f} meV against a harmonic "
                f"{want:.4f} meV"
            )

    def test_the_vacf_route_reports_no_uncertainty(self, welch, vacf):
        """Not an agreement check but the reason the default is what it is: the VACF
        route has no segments to disagree with each other, so it cannot produce the
        spread this module's tolerances are built on. An estimator that cannot say how
        wrong it might be is a poor default even when it is right."""
        assert vacf.pdos_std is None
        assert welch.pdos_std is not None
        assert np.all(np.isfinite(welch.pdos_std))
