"""The intermediate representation: validation, invariants and persistence."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from mdins import units
from mdins.ir import (
    NORMALISATION,
    SpectralMetadata,
    VelocitySpectralDensity,
    bin_widths,
)
from mdins.provenance import Provenance


def make_density(
    symbols=("H", "H", "C"),
    temperature=300.0,
    centre=40.0,
    width=5.0,
    n_freq=512,
    e_max=200.0,
    anisotropy=None,
):
    """A synthetic IR obeying the equipartition sum rule exactly.

    A Gaussian band per atom, normalised so that ``∫ tr P dE = 3 k_B T / m``.
    """
    from ase.data import atomic_masses, atomic_numbers

    masses = np.array([atomic_masses[atomic_numbers[s]] for s in symbols])
    energies = np.linspace(0.0, e_max, n_freq)

    shape = np.exp(-0.5 * ((energies - centre) / width) ** 2)
    shape = shape / np.trapezoid(shape, x=energies)

    weights = np.array(anisotropy if anisotropy is not None else [1 / 3, 1 / 3, 1 / 3])
    weights = weights / weights.sum()

    target = units.thermal_velocity_squared(masses, temperature)  # ∫ tr P dE
    density = (
        target[:, None, None, None]
        * shape[None, :, None, None]
        * np.diag(weights)[None, None, :, :]
    )

    metadata = SpectralMetadata(
        estimator="welch",
        dt=0.005,
        n_frames=20000,
        segment_length=4096,
        n_segments=9,
        overlap=0.5,
        provenance=Provenance(source="synthetic"),
    )
    return VelocitySpectralDensity(
        frequencies=energies,
        density=density,
        symbols=list(symbols),
        masses=masses,
        temperature_md=temperature,
        metadata=metadata,
    )


class TestValidation:
    """Malformed input is rejected at construction, not discovered downstream."""

    def test_accepts_well_formed_input(self):
        sd = make_density()
        assert sd.n_entity == 3
        assert sd.n_freq == 512
        assert sd.metadata.normalisation == NORMALISATION

    def test_rejects_negative_frequencies(self):
        sd = make_density()
        freqs = sd.frequencies.copy()
        freqs[0] = -1.0
        with pytest.raises(ValueError, match="one-sided"):
            VelocitySpectralDensity(
                freqs, sd.density, sd.symbols, sd.masses, 300.0, sd.metadata
            )

    def test_rejects_unsorted_frequencies(self):
        sd = make_density()
        freqs = sd.frequencies.copy()
        freqs[3], freqs[4] = freqs[4], freqs[3]
        with pytest.raises(ValueError, match="ascending"):
            VelocitySpectralDensity(
                freqs, sd.density, sd.symbols, sd.masses, 300.0, sd.metadata
            )

    def test_rejects_wrong_density_shape(self):
        sd = make_density()
        with pytest.raises(ValueError, match="n_entity, n_freq, 3, 3"):
            VelocitySpectralDensity(
                sd.frequencies,
                sd.density[..., 0],
                sd.symbols,
                sd.masses,
                300.0,
                sd.metadata,
            )

    def test_rejects_asymmetric_tensor(self):
        sd = make_density()
        density = sd.density.copy()
        density[0, 10, 0, 1] = 1.0  # break symmetry without its transpose
        with pytest.raises(ValueError, match="symmetric"):
            VelocitySpectralDensity(
                sd.frequencies, density, sd.symbols, sd.masses, 300.0, sd.metadata
            )

    def test_rejects_symbol_count_mismatch(self):
        sd = make_density()
        with pytest.raises(ValueError, match="symbols"):
            VelocitySpectralDensity(
                sd.frequencies,
                sd.density,
                ["H", "H"],
                sd.masses,
                300.0,
                sd.metadata,
            )

    def test_rejects_nonpositive_temperature(self):
        sd = make_density()
        with pytest.raises(ValueError, match="temperature"):
            VelocitySpectralDensity(
                sd.frequencies, sd.density, sd.symbols, sd.masses, 0.0, sd.metadata
            )


class TestDerivedQuantities:
    def test_pdos_is_one_third_of_the_trace(self):
        sd = make_density()
        np.testing.assert_allclose(
            sd.pdos(), np.trace(sd.density, axis1=-2, axis2=-1) / 3.0
        )

    def test_pdos_is_invariant_under_rotation(self):
        """The isotropic method's input must not depend on the lab frame."""
        sd = make_density(anisotropy=[0.6, 0.3, 0.1])
        rotation = np.linalg.qr(np.random.default_rng(3).normal(size=(3, 3)))[0]
        rotated = np.einsum("ai,nfij,bj->nfab", rotation, sd.density, rotation)
        other = VelocitySpectralDensity(
            sd.frequencies, rotated, sd.symbols, sd.masses, 300.0, sd.metadata
        )
        np.testing.assert_allclose(sd.pdos(), other.pdos(), rtol=1e-12)

    def test_total_pdos_sums_over_entities(self):
        sd = make_density()
        np.testing.assert_allclose(sd.total_pdos(), sd.pdos().sum(axis=0))

    def test_heavier_atoms_carry_less_spectral_weight(self):
        sd = make_density(symbols=("H", "C"))
        weight = sd.mean_square_velocity()
        assert weight[0] / weight[1] == pytest.approx(
            sd.masses[1] / sd.masses[0], rel=1e-6
        )


class TestInvariants:
    def test_sum_rule_holds_for_synthetic_input(self):
        make_density().check_sum_rule(rtol=1e-6)

    def test_sum_rule_catches_a_scaling_error(self):
        """The failure this exists to catch: right shape, wrong overall factor."""
        sd = make_density()
        wrong = VelocitySpectralDensity(
            sd.frequencies,
            sd.density * 2.0 * np.pi,
            sd.symbols,
            sd.masses,
            sd.temperature_md,
            sd.metadata,
        )
        with pytest.raises(ValueError, match="sum rule violated"):
            wrong.check_sum_rule()

    def test_sum_rule_catches_a_wrong_temperature(self):
        sd = make_density(temperature=300.0)
        wrong = VelocitySpectralDensity(
            sd.frequencies, sd.density, sd.symbols, sd.masses, 10.0, sd.metadata
        )
        with pytest.raises(ValueError, match="sum rule violated"):
            wrong.check_sum_rule()

    def test_positive_semidefinite_accepts_physical_input(self):
        make_density(anisotropy=[0.7, 0.2, 0.1]).check_positive_semidefinite()

    def test_positive_semidefinite_rejects_a_negative_eigenvalue(self):
        """Blackman-Tukey can produce these; Welch cannot (design.md D5)."""
        sd = make_density()
        density = sd.density.copy()
        density[0, 100, 0, 0] = -density[0, :, 0, 0].max()
        bad = VelocitySpectralDensity(
            sd.frequencies, density, sd.symbols, sd.masses, 300.0, sd.metadata
        )
        with pytest.raises(ValueError, match="positive semi-definite"):
            bad.check_positive_semidefinite()


class TestCentreOfMassCorrection:
    """The sum rule after the centre-of-mass velocity has been projected out.

    Removing ``V = Σ m_j v_j / M`` is not a small correction on a small cell, and it is
    exactly computable: at equilibrium the momentum distribution factorises across atoms
    whatever the potential, so ``⟨v_i·V⟩ = ⟨|V|²⟩ = 3 k_B T / M`` and therefore
    ``⟨|v_i − V|²⟩ = 3 k_B T (1/m_i − 1/M)``. Nothing here is asymptotic in N.
    """

    def test_uncorrected_expectation_is_plain_equipartition(self):
        sd = make_density()
        np.testing.assert_allclose(
            sd.expected_mean_square_velocity(),
            units.thermal_velocity_squared(sd.masses, sd.temperature_md),
        )

    def test_equal_masses_lose_exactly_one_nth(self):
        """For N identical atoms the factor collapses to ``1 − 1/N``, which is the one
        case that can be checked without trusting the mass algebra."""
        n_atoms = 8
        sd = make_density(symbols=("H",) * n_atoms)
        projected = replace(sd, com_projected_mass=float(sd.masses.sum()))
        ratio = (
            projected.expected_mean_square_velocity()
            / sd.expected_mean_square_velocity()
        )
        np.testing.assert_allclose(ratio, 1.0 - 1.0 / n_atoms)

    def test_methane_carbon_loses_a_quarter_of_its_kinetic_energy(self):
        """Not a rounding effect: the carbon in an isolated methane keeps only 75% of
        ``3 k_B T / m_C`` once the molecular centre of mass is frozen. A sum rule that
        ignored this would be wrong by a factor of four on the hydrogen side too."""
        sd = make_density(symbols=("C", "H", "H", "H", "H"))
        total = float(sd.masses.sum())
        projected = replace(sd, com_projected_mass=total)
        expected = projected.expected_mean_square_velocity()
        plain = sd.expected_mean_square_velocity()
        assert expected[0] / plain[0] == pytest.approx(1.0 - sd.masses[0] / total)
        assert expected[0] / plain[0] == pytest.approx(0.251, abs=0.002)
        assert expected[1] / plain[1] == pytest.approx(1.0 - sd.masses[1] / total)

    def test_sum_rule_uses_the_correction(self):
        """A density built to satisfy *uncorrected* equipartition must fail the sum rule
        once it is labelled as centre-of-mass projected, and vice versa."""
        sd = make_density(symbols=("C", "H", "H", "H", "H"))
        sd.check_sum_rule()

        mislabelled = replace(sd, com_projected_mass=float(sd.masses.sum()))
        with pytest.raises(ValueError, match="sum rule"):
            mislabelled.check_sum_rule()

    def test_a_correctly_projected_density_passes(self):
        """The positive half: scale a synthetic density by the exact factor and the
        check accepts it. Without this the test above would also pass if
        ``check_sum_rule`` simply always failed when the mass was set."""
        sd = make_density(symbols=("C", "H", "H", "H", "H"))
        total = float(sd.masses.sum())
        factor = 1.0 - sd.masses / total
        projected = replace(
            sd,
            density=sd.density * factor[:, None, None, None],
            com_projected_mass=total,
        )
        projected.check_sum_rule()

    def test_rejects_a_nonpositive_projected_mass(self):
        with pytest.raises(ValueError, match="com_projected_mass"):
            replace(make_density(), com_projected_mass=0.0)


class TestEntityCounts:
    """One entity can stand for many atoms after grouping, and the total has to stay
    extensive when it does."""

    def test_counts_default_to_one_per_entity(self):
        sd = make_density()
        np.testing.assert_array_equal(sd.counts, [1, 1, 1])

    def test_total_pdos_weights_by_count(self):
        sd = make_density(symbols=("H", "H", "H", "H"))
        grouped = sd.group_by_species()
        assert grouped.n_entity == 1
        np.testing.assert_array_equal(grouped.counts, [4])
        np.testing.assert_allclose(grouped.total_pdos(), sd.total_pdos(), rtol=1e-12)

    def test_grouping_without_counts_would_lose_atoms(self):
        """Guards the specific regression: if ``entity_counts`` were dropped, the total
        DOS of a grouped density would be the *mean* over each species rather than the
        sum, silently shrinking the spectrum by the multiplicity."""
        sd = make_density(symbols=("C", "H", "H", "H", "H"))
        grouped = sd.group_by_species()
        unweighted = grouped.pdos().sum(axis=0)
        assert not np.allclose(unweighted, sd.total_pdos())
        np.testing.assert_allclose(grouped.total_pdos(), sd.total_pdos(), rtol=1e-12)

    def test_rejects_the_wrong_number_of_counts(self):
        with pytest.raises(ValueError, match="entity_counts"):
            replace(make_density(), entity_counts=np.array([1, 2]))

    def test_rejects_a_zero_count(self):
        with pytest.raises(ValueError, match="positive"):
            replace(make_density(), entity_counts=np.array([1, 0, 1]))


class TestGrouping:
    def test_grouping_by_species_preserves_the_sum_rule(self):
        sd = make_density(symbols=("H", "H", "H", "C"))
        grouped = sd.group_by_species()
        assert grouped.symbols == ["C", "H"]
        grouped.check_sum_rule(rtol=1e-6)

    def test_grouping_is_recorded_in_provenance(self):
        grouped = make_density().group_by_species()
        assert "group_by_species" in grouped.metadata.provenance.steps

    def test_rejects_an_isotope_mixture(self):
        """Two atoms labelled H with different masses must not be merged.

        The sum rule is linear in ``1/m``, not in ``m``, so the mean density of a
        protium/deuterium pair does not satisfy it for either mass. Refusing is the only
        answer that cannot be silently wrong.
        """
        sd = make_density(symbols=("H", "H", "C"))
        mixed = replace(sd, masses=np.array([1.008, 2.014, 12.011]))
        with pytest.raises(ValueError, match="differing masses"):
            mixed.group_by_species()

    def test_the_error_names_the_offending_species(self):
        sd = make_density(symbols=("H", "H", "C"))
        mixed = replace(sd, masses=np.array([1.008, 2.014, 12.011]))
        with pytest.raises(ValueError, match="'H'"):
            mixed.group_by_species()


class TestIntegrationRule:
    """Which quadrature rule the IR promises, and why it is not the trapezoid rule."""

    def test_widths_sum_to_the_full_range(self):
        """Bin widths tile the grid: the trapezoid rule over bin centres instead loses
        half a bin at each end, a ``1/(2 (n_freq - 1))`` bias that masquerades as a
        normalisation error."""
        centres = np.linspace(0.0, 200.0, 65)
        spacing = centres[1] - centres[0]
        widths = bin_widths(centres)
        assert widths.sum() == pytest.approx(65 * spacing)
        np.testing.assert_allclose(widths, spacing)

    def test_widths_follow_a_nonuniform_grid(self):
        centres = np.array([0.0, 1.0, 3.0, 7.0])
        np.testing.assert_allclose(bin_widths(centres), [1.0, 1.5, 3.0, 4.0])

    def test_sum_rule_does_not_drift_with_bin_count(self):
        """The trapezoid rule would make the measured integral depend on ``n_freq``;
        this is the test that pins the choice down."""
        coarse = make_density(n_freq=64)
        fine = make_density(n_freq=4096)
        np.testing.assert_allclose(
            coarse.mean_square_velocity(), fine.mean_square_velocity(), rtol=2e-3
        )

    def test_a_single_bin_grid_is_rejected(self):
        with pytest.raises(ValueError, match="at least two"):
            bin_widths(np.array([1.0]))


class TestPersistence:
    def test_round_trip_preserves_everything(self, tmp_path):
        sd = make_density(anisotropy=[0.5, 0.3, 0.2])
        path = tmp_path / "sd.h5"
        sd.to_hdf5(path)
        back = VelocitySpectralDensity.from_hdf5(path)

        np.testing.assert_allclose(back.frequencies, sd.frequencies)
        np.testing.assert_allclose(back.density, sd.density)
        np.testing.assert_allclose(back.masses, sd.masses)
        assert back.symbols == sd.symbols
        assert back.temperature_md == sd.temperature_md
        assert back.metadata.estimator == sd.metadata.estimator
        assert back.metadata.n_segments == sd.metadata.n_segments
        assert back.metadata.normalisation == NORMALISATION

    def test_round_trip_preserves_provenance(self, tmp_path):
        sd = make_density()
        provenance = (
            Provenance(source="traj.extxyz")
            .with_step("remove_com_velocity")
            .with_note("velocities finite-differenced from positions")
        )
        sd = VelocitySpectralDensity(
            sd.frequencies,
            sd.density,
            sd.symbols,
            sd.masses,
            sd.temperature_md,
            SpectralMetadata(
                estimator="vacf",
                dt=0.01,
                n_frames=5000,
                max_lag=1024,
                provenance=provenance,
            ),
        )
        path = tmp_path / "sd.h5"
        sd.to_hdf5(path)
        back = VelocitySpectralDensity.from_hdf5(path)

        assert back.metadata.provenance == provenance
        assert back.metadata.max_lag == 1024
        assert back.metadata.segment_length is None

    def test_off_diagonal_survives_voigt_packing(self, tmp_path):
        """Catches a transposed or mis-ordered Voigt convention."""
        sd = make_density()
        density = sd.density.copy()
        for n, (a, b) in enumerate(((1, 2), (0, 2), (0, 1))):
            density[0, :, a, b] = density[0, :, b, a] = (n + 1) * 1e-3
        sd = VelocitySpectralDensity(
            sd.frequencies, density, sd.symbols, sd.masses, 300.0, sd.metadata
        )
        path = tmp_path / "sd.h5"
        sd.to_hdf5(path)
        np.testing.assert_allclose(
            VelocitySpectralDensity.from_hdf5(path).density, density
        )

    def test_rejects_a_foreign_file(self, tmp_path):
        import h5py

        path = tmp_path / "other.h5"
        with h5py.File(path, "w") as handle:
            handle.create_dataset("something", data=[1, 2, 3])
        with pytest.raises(ValueError, match="not an mdins"):
            VelocitySpectralDensity.from_hdf5(path)

    def test_round_trip_preserves_the_new_optional_fields(self, tmp_path):
        """``pdos_std``, ``entity_counts`` and ``com_projected_mass`` all change what
        downstream code computes, so losing one on the way through HDF5 would change a
        result rather than only a label."""
        sd = make_density(symbols=("C", "H", "H", "H", "H"))
        rng = np.random.default_rng(0)
        sd = replace(
            sd,
            pdos_std=rng.random((sd.n_entity, sd.n_freq)),
            com_projected_mass=float(sd.masses.sum()),
        ).group_by_species()

        path = tmp_path / "ir.h5"
        sd.to_hdf5(path)
        loaded = VelocitySpectralDensity.from_hdf5(path)

        np.testing.assert_array_equal(loaded.entity_counts, sd.entity_counts)
        np.testing.assert_allclose(loaded.pdos_std, sd.pdos_std)
        assert loaded.com_projected_mass == pytest.approx(sd.com_projected_mass)
        np.testing.assert_allclose(loaded.total_pdos(), sd.total_pdos())
        np.testing.assert_allclose(
            loaded.expected_mean_square_velocity(), sd.expected_mean_square_velocity()
        )

    def test_round_trip_preserves_the_achieved_resolution(self, tmp_path):
        sd = make_density()
        sd = replace(sd, metadata=replace(sd.metadata, energy_resolution=0.25))
        path = tmp_path / "ir.h5"
        sd.to_hdf5(path)
        loaded = VelocitySpectralDensity.from_hdf5(path)
        assert loaded.metadata.energy_resolution == 0.25

    def test_absent_optional_fields_stay_absent(self, tmp_path):
        """A file written without them must load as ``None`` rather than as zeros,
        which would quietly disable the centre-of-mass correction."""
        path = tmp_path / "ir.h5"
        make_density().to_hdf5(path)
        loaded = VelocitySpectralDensity.from_hdf5(path)
        assert loaded.pdos_std is None
        assert loaded.entity_counts is None
        assert loaded.com_projected_mass is None


class TestDegenerateInput:
    """Cases where the sum rule has nothing to say, and must say so rather than pass.

    Every guard here exists because the natural spelling of the check silently succeeds:
    ``nan > rtol`` is False, ``nan <= 0`` is False, and ``0/0`` is ``nan``. A validation
    that cannot fail is worse than none, because it is reported as having passed.
    """

    def test_rejects_a_nan_temperature(self):
        with pytest.raises(ValueError, match="finite"):
            replace(make_density(), temperature_md=float("nan"))

    def test_rejects_an_infinite_temperature(self):
        with pytest.raises(ValueError, match="finite"):
            replace(make_density(), temperature_md=float("inf"))

    def test_rejects_a_projected_mass_below_the_heaviest_entity(self):
        """The correction is ``1 - m_i/M``; an ``M`` smaller than some ``m_i`` makes the
        expected mean square velocity negative, and every comparison against it
        meaningless."""
        sd = make_density(symbols=("C", "H"))
        with pytest.raises(ValueError, match="heaviest entity"):
            replace(sd, com_projected_mass=1.0)

    def test_a_nan_residual_is_a_failure(self):
        """Constructed directly, since the validation above now makes it hard to reach
        by accident. If a NaN ever does get into the density, the sum rule must not
        report success."""
        sd = make_density()
        density = sd.density.copy()
        density[0, 10] = np.nan
        broken = replace(sd, density=density)
        with pytest.raises(ValueError, match="sum rule violated"):
            broken.check_sum_rule()

    def test_a_single_atom_projected_onto_itself_must_carry_no_velocity(self):
        """``M = m`` makes the expectation exactly zero. The check then has to be
        absolute: a relative one would divide zero by zero and pass."""
        sd = make_density(symbols=("H",))
        projected = replace(sd, com_projected_mass=float(sd.masses[0]))
        assert projected.expected_mean_square_velocity()[0] == 0.0
        with pytest.raises(ValueError, match="no velocity at all"):
            projected.check_sum_rule()

    def test_a_genuinely_empty_spectrum_passes_that_check(self):
        """The other side of it: zero expected and zero measured is consistent."""
        sd = make_density(symbols=("H",))
        projected = replace(
            sd,
            density=np.zeros_like(sd.density),
            com_projected_mass=float(sd.masses[0]),
        )
        projected.check_sum_rule()
