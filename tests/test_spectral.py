"""Layer 1: spectral estimation against analytically known answers.

Velocities are generated directly rather than by running MD, so the expected spectrum is
exact and the tests run in milliseconds. This layer is what catches normalisation, unit
and factor-of-2π errors; by the time a real trajectory is involved, those are expensive
to distinguish from physics.
"""

from __future__ import annotations

import re

import numpy as np
import pytest

from mdins import units
from mdins.provenance import Provenance
from mdins.spectral import (
    CHUNK_BYTES,
    _atom_chunks,
    velocity_spectral_density,
)
from mdins.trajectory import VelocityTrajectory
from mdins.units import HBAR, PLANCK

ESTIMATORS = ["welch", "vacf"]


def make_trajectory(velocities, dt=0.002, masses=None, symbols=None):
    """Wrap a synthetic velocity array.

    Deliberately does *not* claim ``remove_com_velocity`` in provenance. These arrays
    are constructed drift-free, so the real check passes on its own, and asserting a
    projection that was never applied would make the sum-rule correction wrong — the
    correction is keyed off that provenance step.
    """
    n_atoms = velocities.shape[1]
    return VelocityTrajectory(
        velocities=np.asarray(velocities, dtype=np.float64),
        masses=np.full(n_atoms, 1.0) if masses is None else np.asarray(masses),
        symbols=["H"] * n_atoms if symbols is None else list(symbols),
        dt=dt,
        provenance=Provenance(source="synthetic"),
    )


def oscillating(energy_mev, amplitude=1.0, axis=0, n_frames=8192, dt=0.002, phase=0.0):
    """A single atom oscillating along one axis at a known energy."""
    nu = energy_mev / PLANCK  # 1/ps
    t = np.arange(n_frames) * dt
    velocities = np.zeros((n_frames, 1, 3))
    velocities[:, 0, axis] = amplitude * np.cos(2 * np.pi * nu * t + phase)
    return velocities


def integrate(sd):
    """The same rule the estimator conserves; see mdins.ir.bin_widths."""
    return sd.mean_square_velocity()


def _quadrature_fraction(sd):
    """The A2 diagnostic, recovered from the provenance note that carries it."""
    notes = [n for n in sd.metadata.provenance.notes if "ambiguity A2" in n]
    assert len(notes) == 1, f"expected exactly one A2 note, got {notes}"
    return float(re.search(r"([\d.]+)%", notes[0]).group(1)) / 100.0


class TestPeakPosition:
    """Where the peak lands. A 2π error moves it by a factor of 6.28."""

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    @pytest.mark.parametrize("energy", [5.0, 40.0, 150.0])
    def test_peak_is_at_the_input_energy(self, estimator, energy):
        traj = make_trajectory(oscillating(energy))
        sd = velocity_spectral_density(
            traj, e_max=400.0, n_bins=2048, estimator=estimator
        )
        peak = sd.frequencies[np.argmax(sd.pdos()[0])]
        assert peak == pytest.approx(energy, abs=2.0)

    def test_peak_is_not_at_the_angular_frequency(self):
        """Guards the specific slip of storing rad/ps where meV is expected."""
        traj = make_trajectory(oscillating(40.0))
        sd = velocity_spectral_density(traj, e_max=400.0, n_bins=2048)
        peak = sd.frequencies[np.argmax(sd.pdos()[0])]
        assert abs(peak - 2 * np.pi * 40.0) > 10.0


class TestNormalisation:
    """That ``∫ tr P dE`` is the mean square velocity, in absolute units (D6)."""

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_single_oscillator_integrates_to_half_the_amplitude_squared(
        self, estimator
    ):
        amplitude = 3.0
        traj = make_trajectory(oscillating(40.0, amplitude=amplitude))
        sd = velocity_spectral_density(
            traj, e_max=400.0, n_bins=2048, estimator=estimator
        )
        assert integrate(sd)[0] == pytest.approx(amplitude**2 / 2, rel=0.02)

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_white_noise_recovers_its_variance(self, estimator):
        rng = np.random.default_rng(0)
        sigma = 2.5
        velocities = rng.normal(scale=sigma, size=(16384, 4, 3))
        traj = make_trajectory(velocities)
        sd = velocity_spectral_density(traj, n_bins=256, estimator=estimator)
        np.testing.assert_allclose(integrate(sd), 3 * sigma**2, rtol=0.05)

    def test_white_noise_spectrum_is_flat(self):
        rng = np.random.default_rng(1)
        traj = make_trajectory(rng.normal(size=(32768, 2, 3)))
        sd = velocity_spectral_density(traj, n_bins=64)
        pdos = sd.pdos()[0][2:-2]
        assert pdos.std() / pdos.mean() < 0.2

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_normalisation_is_independent_of_the_output_grid(self, estimator):
        """Rebinning conserves the integral, so the sum rule cannot be an artefact of
        bin count."""
        traj = make_trajectory(oscillating(40.0, amplitude=2.0))
        coarse = velocity_spectral_density(
            traj, e_max=400.0, n_bins=128, estimator=estimator
        )
        fine = velocity_spectral_density(
            traj, e_max=400.0, n_bins=4096, estimator=estimator
        )
        assert integrate(coarse)[0] == pytest.approx(integrate(fine)[0], rel=0.02)

    def test_normalisation_is_independent_of_segment_length(self):
        traj = make_trajectory(oscillating(40.0, amplitude=2.0, n_frames=16384))
        short = velocity_spectral_density(traj, e_max=400.0, segment_length=512)
        long = velocity_spectral_density(traj, e_max=400.0, segment_length=4096)
        assert integrate(short)[0] == pytest.approx(integrate(long)[0], rel=0.02)

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_thermal_velocities_satisfy_the_sum_rule(self, estimator):
        """The end-to-end invariant: Maxwell-Boltzmann velocities in, 3kT/m out."""
        from ase import Atoms
        from ase import units as ase_units
        from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

        rng = np.random.default_rng(4)
        atoms = Atoms("CH4", positions=rng.random((5, 3)))
        frames = []
        for _ in range(4096):
            MaxwellBoltzmannDistribution(atoms, temperature_K=300.0, rng=rng)
            frames.append(atoms.get_velocities() * ase_units.fs * 1e3)

        traj = make_trajectory(
            np.array(frames), masses=atoms.get_masses(), symbols=["C", *["H"] * 4]
        )
        sd = velocity_spectral_density(traj, n_bins=128, estimator=estimator)
        assert sd.temperature_md == pytest.approx(300.0, rel=0.1)
        sd.check_sum_rule(rtol=0.1)


class TestTensorStructure:
    """The part the isotropic method throws away and the anisotropic method needs."""

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_axis_aligned_motion_populates_one_diagonal_element(self, estimator):
        traj = make_trajectory(oscillating(40.0, axis=1))
        sd = velocity_spectral_density(
            traj, e_max=400.0, n_bins=512, estimator=estimator
        )
        total = np.trapezoid(sd.density[0], x=sd.frequencies, axis=0)
        assert total[1, 1] == pytest.approx(0.5, rel=0.02)
        assert abs(total[0, 0]) < 1e-6
        assert abs(total[2, 2]) < 1e-6

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_anisotropic_amplitudes_give_the_right_diagonal_ratio(self, estimator):
        amplitudes = np.array([3.0, 2.0, 1.0])
        velocities = sum(
            oscillating(40.0, amplitude=a, axis=d) for d, a in enumerate(amplitudes)
        )
        sd = velocity_spectral_density(
            make_trajectory(velocities), e_max=400.0, n_bins=512, estimator=estimator
        )
        total = np.trapezoid(sd.density[0], x=sd.frequencies, axis=0)
        np.testing.assert_allclose(np.diag(total), amplitudes**2 / 2, rtol=0.02)

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_in_phase_motion_gives_a_positive_off_diagonal(self, estimator):
        velocities = oscillating(40.0, axis=0) + oscillating(40.0, axis=1)
        sd = velocity_spectral_density(
            make_trajectory(velocities), e_max=400.0, n_bins=512, estimator=estimator
        )
        total = np.trapezoid(sd.density[0], x=sd.frequencies, axis=0)
        assert total[0, 1] == pytest.approx(0.5, rel=0.02)

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_antiphase_motion_flips_the_sign(self, estimator):
        """Catches a dropped sign or an absolute value in the cross term."""
        velocities = oscillating(40.0, axis=0) + oscillating(40.0, axis=1, phase=np.pi)
        sd = velocity_spectral_density(
            make_trajectory(velocities), e_max=400.0, n_bins=512, estimator=estimator
        )
        total = np.trapezoid(sd.density[0], x=sd.frequencies, axis=0)
        assert total[0, 1] == pytest.approx(-0.5, rel=0.02)

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_quadrature_motion_has_no_cospectrum(self, estimator):
        """Circular motion: all the cross-correlation is in the imaginary part, which
        the co-spectrum convention discards (ambiguity A2)."""
        velocities = oscillating(40.0, axis=0) + oscillating(
            40.0, axis=1, phase=np.pi / 2
        )
        sd = velocity_spectral_density(
            make_trajectory(velocities), e_max=400.0, n_bins=512, estimator=estimator
        )
        total = np.trapezoid(sd.density[0], x=sd.frequencies, axis=0)
        assert abs(total[0, 1]) < 0.02 * total[0, 0]

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_quadrature_discard_is_recorded(self, estimator):
        velocities = oscillating(40.0, axis=0) + oscillating(
            40.0, axis=1, phase=np.pi / 2
        )
        sd = velocity_spectral_density(
            make_trajectory(velocities), e_max=400.0, estimator=estimator
        )
        assert any("ambiguity A2" in note for note in sd.metadata.provenance.notes)

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_quadrature_diagnostic_is_large_for_circular_motion(self, estimator):
        """Asserted by value, not by the presence of a note.

        Circular motion puts all of its cross-correlation in the discarded component,
        so the reported fraction has to be of order one. A diagnostic that always read
        zero would satisfy a string match and tell nobody anything (design.md A2).
        """
        velocities = oscillating(40.0, axis=0) + oscillating(
            40.0, axis=1, phase=np.pi / 2
        )
        sd = velocity_spectral_density(
            make_trajectory(velocities), e_max=400.0, estimator=estimator
        )
        assert _quadrature_fraction(sd) > 0.3

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_quadrature_diagnostic_is_small_for_in_phase_motion(self, estimator):
        """The other half: linear motion has nothing in quadrature, so a diagnostic
        stuck at one would be just as useless as one stuck at zero."""
        velocities = oscillating(40.0, axis=0) + oscillating(40.0, axis=1)
        sd = velocity_spectral_density(
            make_trajectory(velocities), e_max=400.0, estimator=estimator
        )
        assert _quadrature_fraction(sd) < 0.05

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_tensor_transforms_correctly_under_rotation(self, estimator):
        """P is a tensor: rotating the velocities must rotate the result."""
        rotation = np.linalg.qr(np.random.default_rng(7).normal(size=(3, 3)))[0]
        velocities = sum(
            oscillating(40.0, amplitude=a, axis=d)
            for d, a in enumerate([3.0, 2.0, 1.0])
        )
        plain = velocity_spectral_density(
            make_trajectory(velocities), e_max=400.0, n_bins=256, estimator=estimator
        )
        rotated = velocity_spectral_density(
            make_trajectory(velocities @ rotation.T),
            e_max=400.0,
            n_bins=256,
            estimator=estimator,
        )
        expected = np.einsum("ai,nfij,bj->nfab", rotation, plain.density, rotation)
        np.testing.assert_allclose(rotated.density, expected, atol=1e-8)


class TestPositiveSemidefiniteness:
    """The reason Welch is the default (design.md D5)."""

    @pytest.mark.parametrize("n_bins", [256, 1024, 4096])
    def test_vacf_fails_where_welch_succeeds(self, n_bins):
        """The whole content of D5, as an executable statement.

        Welch averages outer products ``V*Vᵀ``, so each tensor is PSD by construction.
        Blackman-Tukey has no such guarantee, and the Hann lag window's sidelobes drive
        eigenvalues negative on a signal as ordinary as one sinusoid. If this ever stops
        failing for VACF, the estimator changed and the justification for the default
        needs rechecking — it does not mean the problem went away.
        """
        traj = make_trajectory(oscillating(40.0, amplitude=2.0))
        welch = velocity_spectral_density(traj, e_max=400.0, n_bins=n_bins)
        vacf = velocity_spectral_density(
            traj, e_max=400.0, n_bins=n_bins, estimator="vacf"
        )

        welch.check_positive_semidefinite()
        assert welch.eigenvalues().min() >= -1e-12 * welch.eigenvalues().max()

        with pytest.raises(ValueError, match="not positive semi-definite"):
            vacf.check_positive_semidefinite()

    def test_welch_stays_positive_semidefinite_on_a_degenerate_tensor(self):
        """Purely one-dimensional motion: two eigenvalues are exactly zero, so any
        normalisation slip that overshoots shows up as a negative one."""
        traj = make_trajectory(oscillating(40.0, axis=2, amplitude=5.0))
        velocity_spectral_density(
            traj, e_max=400.0, n_bins=512
        ).check_positive_semidefinite()

    def test_welch_is_positive_semidefinite_on_noise(self):
        rng = np.random.default_rng(2)
        traj = make_trajectory(rng.normal(size=(8192, 3, 3)))
        velocity_spectral_density(traj, n_bins=512).check_positive_semidefinite()

    def test_welch_is_positive_semidefinite_on_correlated_motion(self):
        rng = np.random.default_rng(5)
        base = rng.normal(size=(8192, 2, 1))
        velocities = np.concatenate([base, 0.7 * base, -0.3 * base], axis=2)
        velocities += 0.1 * rng.normal(size=velocities.shape)
        traj = make_trajectory(velocities)
        velocity_spectral_density(traj, n_bins=256).check_positive_semidefinite()


class TestEstimatorAgreement:
    """Welch and VACF are independent implementations of the same quantity."""

    def test_agreement_on_a_multimode_signal(self):
        velocities = (
            oscillating(20.0, amplitude=1.0, axis=0, n_frames=16384)
            + oscillating(60.0, amplitude=2.0, axis=1, n_frames=16384)
            + oscillating(120.0, amplitude=0.5, axis=2, n_frames=16384)
        )
        traj = make_trajectory(velocities)
        welch = velocity_spectral_density(traj, e_max=300.0, n_bins=128)
        vacf = velocity_spectral_density(
            traj, e_max=300.0, n_bins=128, estimator="vacf"
        )
        assert integrate(welch)[0] == pytest.approx(integrate(vacf)[0], rel=0.05)
        # Compare shapes, since the two smooth differently.
        overlap = np.sum(welch.pdos()[0] * vacf.pdos()[0]) / np.sqrt(
            np.sum(welch.pdos()[0] ** 2) * np.sum(vacf.pdos()[0] ** 2)
        )
        assert overlap > 0.95

    def test_agreement_is_within_the_welch_inter_segment_spread(self):
        """The tolerance comes from the data, not from me.

        design.md D5 and plan.md §5 both require the two estimators to agree "within the
        Welch inter-segment spread". That spread is computed and calibrated
        (:meth:`TestUncertainty.test_std_error_matches_the_scatter_across_independent_runs`),
        so the comparison can use it instead of a hand-picked percentage. This is the
        pattern M2's Euphonic comparison has to inherit.
        """
        rng = np.random.default_rng(21)
        traj = make_trajectory(rng.normal(scale=1.3, size=(32768, 2, 3)))
        welch = velocity_spectral_density(
            traj, n_bins=32, segment_length=2048, overlap=0.0
        )
        vacf = velocity_spectral_density(
            traj, n_bins=32, max_lag=2048, estimator="vacf"
        )

        # Interior bins only: the end bins are half-width and edge-affected, and the
        # VACF route leaks differently there.
        deviation = np.abs(welch.pdos()[:, 1:-1] - vacf.pdos()[:, 1:-1])
        spread = welch.pdos_std[:, 1:-1]
        worst = float((deviation / spread).max())
        assert worst < 3.0, (
            f"estimators differ by up to {worst:.1f} standard errors, which is more "
            "than sampling noise explains"
        )

    def test_agreement_on_thermal_noise(self):
        rng = np.random.default_rng(11)
        traj = make_trajectory(rng.normal(scale=1.7, size=(16384, 3, 3)))
        welch = velocity_spectral_density(traj, n_bins=64)
        vacf = velocity_spectral_density(traj, n_bins=64, estimator="vacf")
        np.testing.assert_allclose(integrate(welch), integrate(vacf), rtol=0.05)


class TestUncertainty:
    """Welch's inter-segment spread, used to set validation tolerances."""

    def test_std_error_is_reported(self):
        rng = np.random.default_rng(6)
        traj = make_trajectory(rng.normal(size=(8192, 2, 3)))
        sd = velocity_spectral_density(traj, n_bins=64, segment_length=512)
        assert sd.pdos_std is not None
        assert sd.pdos_std.shape == sd.pdos().shape
        assert np.all(sd.pdos_std > 0.0)

    def test_std_error_falls_with_more_data(self):
        rng = np.random.default_rng(8)
        velocities = rng.normal(size=(65536, 1, 3))
        short = velocity_spectral_density(
            make_trajectory(velocities[:16384]),
            n_bins=32,
            segment_length=1024,
            overlap=0.0,
        )
        long = velocity_spectral_density(
            make_trajectory(velocities), n_bins=32, segment_length=1024, overlap=0.0
        )
        # Four times the data, so half the standard error.
        assert long.pdos_std.mean() == pytest.approx(
            0.5 * short.pdos_std.mean(), rel=0.25
        )

    def test_std_error_is_set_by_duration_not_by_segmentation(self):
        """Non-obvious, and worth pinning: at fixed total length, shortening the
        segments makes each one noisier but yields proportionally more of them, and the
        two effects cancel. Segment length trades resolution against variance, not
        against total information."""
        rng = np.random.default_rng(8)
        traj = make_trajectory(rng.normal(size=(32768, 1, 3)))
        few = velocity_spectral_density(
            traj, n_bins=32, segment_length=8192, overlap=0.0
        )
        many = velocity_spectral_density(
            traj, n_bins=32, segment_length=1024, overlap=0.0
        )
        assert many.pdos_std.mean() == pytest.approx(few.pdos_std.mean(), rel=0.3)

    def test_std_error_matches_the_scatter_across_independent_runs(self):
        """The error bar has to be calibrated, not merely present: validation
        tolerances are derived from it (plan.md §5)."""
        rng = np.random.default_rng(12)
        estimates = []
        reported = []
        for _ in range(40):
            traj = make_trajectory(rng.normal(size=(4096, 1, 3)))
            sd = velocity_spectral_density(
                traj, n_bins=16, segment_length=512, overlap=0.0
            )
            estimates.append(sd.pdos()[0])
            reported.append(sd.pdos_std[0])

        empirical = np.std(estimates, axis=0, ddof=1)
        predicted = np.mean(reported, axis=0)
        # Compare interior bins; the first and last are half-width and edge-affected.
        ratio = (predicted[1:-1] / empirical[1:-1]).mean()
        assert ratio == pytest.approx(1.0, abs=0.25)

    def test_vacf_reports_no_std_error(self):
        traj = make_trajectory(oscillating(40.0))
        sd = velocity_spectral_density(traj, e_max=200.0, estimator="vacf")
        assert sd.pdos_std is None


class TestGuards:
    def test_rejects_an_unknown_estimator(self):
        traj = make_trajectory(oscillating(40.0))
        with pytest.raises(ValueError, match="unknown estimator"):
            velocity_spectral_density(traj, estimator="periodogram")

    def test_rejects_a_grid_beyond_nyquist(self):
        traj = make_trajectory(oscillating(40.0), dt=0.01)
        with pytest.raises(ValueError, match="Nyquist"):
            velocity_spectral_density(traj, e_max=1000.0)

    def test_rejects_a_drifting_trajectory(self):
        velocities = oscillating(40.0) + 5.0
        traj = make_trajectory(velocities)
        with pytest.raises(ValueError, match="centre-of-mass"):
            velocity_spectral_density(traj, e_max=200.0)

    def test_accepts_a_trajectory_whose_drift_was_removed(self):
        rng = np.random.default_rng(9)
        velocities = rng.normal(size=(4096, 6, 3)) + 5.0
        traj = make_trajectory(velocities).remove_com_velocity()
        velocity_spectral_density(traj, n_bins=32)

    def test_oscillating_single_atom_is_not_mistaken_for_drift(self):
        """A lone oscillator has a large instantaneous centre-of-mass speed and no
        drift; the check must distinguish them."""
        traj = make_trajectory(oscillating(40.0))
        velocity_spectral_density(traj, e_max=200.0)

    def test_rejects_too_long_a_segment(self):
        traj = make_trajectory(oscillating(40.0, n_frames=512))
        with pytest.raises(ValueError, match="exceeds"):
            velocity_spectral_density(traj, e_max=200.0, segment_length=1024)

    def test_metadata_records_the_estimator_settings(self):
        traj = make_trajectory(oscillating(40.0, n_frames=8192))
        sd = velocity_spectral_density(
            traj, e_max=200.0, segment_length=1024, overlap=0.5
        )
        assert sd.metadata.estimator == "welch"
        assert sd.metadata.segment_length == 1024
        assert sd.metadata.n_segments == 15
        assert any("welch" in step for step in sd.metadata.provenance.steps)


class TestResolution:
    """Resolution is set by the segment, not by the output grid or the total duration.

    This is the most frequently misunderstood knob in the pipeline: a long run cut into
    short segments resolves poorly no matter how many output bins are asked for, and
    interpolating onto a finer grid makes a spectrum that *looks* sharper without being
    so. The requested resolution is therefore validated against what the estimator will
    actually deliver.
    """

    def test_achieved_resolution_is_recorded(self):
        traj = make_trajectory(oscillating(40.0, n_frames=8192))
        sd = velocity_spectral_density(traj, e_max=200.0, segment_length=1024)
        assert sd.metadata.energy_resolution == pytest.approx(PLANCK / (1024 * 0.002))

    def test_achieved_resolution_tracks_max_lag_for_the_vacf_route(self):
        traj = make_trajectory(oscillating(40.0, n_frames=8192))
        sd = velocity_spectral_density(
            traj, e_max=200.0, max_lag=2048, estimator="vacf"
        )
        assert sd.metadata.energy_resolution == pytest.approx(PLANCK / (2048 * 0.002))

    def test_achieved_resolution_does_not_depend_on_the_output_grid(self):
        traj = make_trajectory(oscillating(40.0, n_frames=8192))
        coarse = velocity_spectral_density(
            traj, e_max=200.0, segment_length=1024, n_bins=32
        )
        fine = velocity_spectral_density(
            traj, e_max=200.0, segment_length=1024, n_bins=8192
        )
        assert coarse.metadata.energy_resolution == fine.metadata.energy_resolution

    @pytest.mark.parametrize("estimator", ESTIMATORS)
    def test_rejects_a_resolution_the_segment_cannot_deliver(self, estimator):
        """1 meV is comfortably within what 16 ps of trajectory supports (0.25 meV) and
        far outside what a 256-frame segment supports (8 meV). Checking the duration
        alone, which is what "energy resolution" usually means, would accept this."""
        traj = make_trajectory(oscillating(40.0, n_frames=8192))
        assert traj.energy_resolution < 1.0
        with pytest.raises(ValueError, match="Resolution is set by the segment"):
            velocity_spectral_density(
                traj,
                e_max=200.0,
                segment_length=256,
                max_lag=256,
                estimator=estimator,
                e_resolution=1.0,
            )

    def test_rejects_a_resolution_the_whole_trajectory_cannot_deliver(self):
        """Caught earlier, by the trajectory rather than the estimator, and with the
        advice that matters in that case: run for longer."""
        traj = make_trajectory(oscillating(40.0, n_frames=1024))
        with pytest.raises(ValueError, match="Run for at least"):
            velocity_spectral_density(traj, e_max=200.0, e_resolution=0.01)

    def test_accepts_a_resolution_the_segment_can_deliver(self):
        traj = make_trajectory(oscillating(40.0, n_frames=8192))
        sd = velocity_spectral_density(
            traj, e_max=200.0, segment_length=2048, e_resolution=2.0
        )
        assert sd.metadata.energy_resolution <= 2.0

    def test_a_finer_output_grid_does_not_buy_resolution(self):
        """The peak of a pure tone is broadened by the segment window. Quadrupling the
        bin count must not narrow it, or the resolution figure would be a lie."""
        traj = make_trajectory(oscillating(40.0, n_frames=8192))
        widths = []
        for n_bins in (512, 2048):
            sd = velocity_spectral_density(
                traj, e_max=200.0, segment_length=512, n_bins=n_bins
            )
            pdos = sd.pdos()[0]
            above = pdos > 0.5 * pdos.max()
            energies = sd.frequencies[above]
            widths.append(energies.max() - energies.min())
        assert widths[1] == pytest.approx(widths[0], abs=2.0)


class TestMemoryBounding:
    """The estimator accumulates onto the output grid, so its footprint is bounded by
    the chunk size rather than by the number of atoms (design.md §2[B])."""

    def test_chunks_cover_every_atom_exactly_once(self):
        chunks = _atom_chunks(1000, 4096)
        covered = [a for lo, hi in chunks for a in range(lo, hi)]
        assert covered == list(range(1000))

    def test_chunk_size_shrinks_as_the_native_grid_grows(self):
        small = _atom_chunks(100_000, 1024)[0]
        large = _atom_chunks(100_000, 65_536)[0]
        assert small[1] > large[1]

    def test_a_single_atom_is_never_split(self):
        """Even when one atom's cross-spectrum exceeds the budget: splitting it is not
        possible, and returning an empty chunk list would drop it silently."""
        assert _atom_chunks(3, 10**9) == [(0, 1), (1, 2), (2, 3)]

    def test_the_budget_is_respected(self):
        n_native = 8192
        lo, hi = _atom_chunks(10_000, n_native)[0]
        assert (hi - lo) * n_native * 9 * 16 * 3 <= CHUNK_BYTES

    def test_chunking_does_not_change_the_answer(self, monkeypatch):
        """The whole point: the result is identical whether the atoms are processed in
        one pass or in many."""
        rng = np.random.default_rng(11)
        traj = make_trajectory(rng.normal(size=(4096, 9, 3)))
        whole = velocity_spectral_density(traj, n_bins=64)

        monkeypatch.setattr("mdins.spectral.CHUNK_BYTES", 1)
        chunked = velocity_spectral_density(traj, n_bins=64)
        assert len(_atom_chunks(9, 2048)) == 9

        np.testing.assert_allclose(chunked.density, whole.density, rtol=1e-12)


class TestCentreOfMassProvenance:
    """The sum-rule correction is keyed off provenance, so the wiring between stage A
    and the IR has to be exact."""

    def test_projection_is_propagated_to_the_ir(self):
        rng = np.random.default_rng(12)
        traj = make_trajectory(rng.normal(size=(4096, 6, 3))).remove_com_velocity()
        sd = velocity_spectral_density(traj, n_bins=32)
        assert sd.com_projected_mass == pytest.approx(traj.total_mass)

    def test_an_unprojected_trajectory_leaves_it_unset(self):
        traj = make_trajectory(oscillating(40.0))
        assert velocity_spectral_density(traj, e_max=200.0).com_projected_mass is None

    def test_projected_thermal_velocities_satisfy_the_corrected_sum_rule(self):
        """The end-to-end statement, on a cell small enough that the correction is 17%
        and an uncorrected check would fail outright.

        Independent Maxwell-Boltzmann draws per frame, so the only thing between the
        input and the sum rule is the projection and the estimator.
        """
        rng = np.random.default_rng(13)
        n_atoms = 6
        mass = 39.948
        sigma = np.sqrt(units.KB * 300.0 / mass)
        velocities = rng.normal(scale=sigma, size=(8192, n_atoms, 3))
        traj = make_trajectory(
            velocities, masses=np.full(n_atoms, mass), symbols=["Ar"] * n_atoms
        ).remove_com_velocity()

        sd = velocity_spectral_density(traj, n_bins=64)
        ratio = sd.mean_square_velocity() / units.thermal_velocity_squared(
            sd.masses, sd.temperature_md
        )
        np.testing.assert_allclose(ratio, 1.0 - 1.0 / n_atoms, rtol=0.05)
        sd.check_sum_rule(rtol=0.05)


class TestLineshape:
    """Not just where the peak is, but how wide it is.

    Every other test here uses an undamped tone, whose measured width is entirely the
    window function — so nothing yet checks that a *real* linewidth survives the
    estimator. The velocity spectral density of a damped oscillator
    ``v(t) = e^{-t/τ} cos(ω₀ t)`` is a Lorentzian of half-width ``ħ/τ`` in energy. This
    is what makes the instrument stage [F] meaningful later: convolving with a
    resolution function is only honest if the intrinsic width underneath it is right.

    The width is obtained by fitting, not by thresholding at half the peak bin. The peak
    bin of a noisy spectrum is biased high, which biases a thresholded width low — by
    34% for the narrowest case here, enough to hide a real error.
    """

    #: Energy of the oscillator, in meV. Chosen well away from zero so the wings are not
    #: cut off by the lower edge of the grid.
    CENTRE = 40.0

    @staticmethod
    def damped(tau_ps, n_frames=65536, n_atoms=2, dt=0.002, seed=0):
        """An ensemble of randomly-timed, randomly-phased ringdowns.

        A single decaying transient is not a stationary process and has no spectral
        density; Welch would report its segment-to-segment decay as structure. Filtered
        Poisson noise is stationary and its spectral density is exactly the squared
        response of one ringdown, which is the Lorentzian being tested for. The kick
        rate scales with ``1/τ`` so that the estimate is equally well converged at
        every damping.
        """
        rng = np.random.default_rng(seed)
        nu = TestLineshape.CENTRE / PLANCK
        t = np.arange(n_frames) * dt
        velocities = np.zeros((n_frames, n_atoms, 3))
        for atom in range(n_atoms):
            signal = np.zeros(n_frames)
            for start in rng.integers(
                0, n_frames, size=int(4 * n_frames * dt / tau_ps)
            ):
                lag = t[: n_frames - start]
                signal[start:] += np.exp(-lag / tau_ps) * np.cos(
                    2 * np.pi * nu * lag + rng.uniform(0, 2 * np.pi)
                )
            velocities[:, atom, 0] = signal
        return velocities

    @staticmethod
    def fit_width(sd):
        """Half-width at half maximum of the peak, in meV, by least squares.

        The constant baseline absorbs the multi-kick background; without it the fit
        trades width against offset.
        """
        from scipy.optimize import curve_fit

        def lorentzian(energy, amplitude, centre, gamma, baseline):
            return amplitude / (1.0 + ((energy - centre) / gamma) ** 2) + baseline

        pdos = sd.pdos().mean(axis=0)
        near = np.abs(sd.frequencies - TestLineshape.CENTRE) < 20.0
        (_, centre, gamma, _), _ = curve_fit(
            lorentzian,
            sd.frequencies[near],
            pdos[near],
            p0=[pdos[near].max(), TestLineshape.CENTRE, 1.0, 0.0],
        )
        assert centre == pytest.approx(TestLineshape.CENTRE, abs=1.0)
        return abs(gamma)

    def estimate(self, tau_ps, seed=0):
        traj = make_trajectory(self.damped(tau_ps, seed=seed))
        return velocity_spectral_density(
            traj, e_max=120.0, n_bins=512, segment_length=8192
        )

    @pytest.mark.parametrize("tau_ps", [0.3, 0.6])
    def test_width_matches_the_damping_rate(self, tau_ps):
        """HWHM is ``ħ/τ``. A missing 2π would be out by 6.28 and confusing the
        amplitude decay with the energy decay by 2, so the band is wide enough to be
        insensitive to estimator bias and still far narrower than either error.

        The measured width runs a little high because the Hann window on an 8192-frame
        segment adds about 0.25 meV of instrumental broadening on top of the intrinsic
        line; both damping times here are chosen several times wider than that.
        """
        gamma = self.fit_width(self.estimate(tau_ps))
        expected = HBAR / tau_ps
        assert expected > 4 * PLANCK / (8192 * 0.002)
        assert 0.85 < gamma / expected < 1.35, (
            f"fitted HWHM {gamma:.3f} meV against an expected {expected:.3f} meV"
        )

    def test_a_faster_decay_is_a_broader_line(self):
        """The monotonic statement, free of any window bias since it affects both."""
        assert self.fit_width(self.estimate(0.2)) > 2 * self.fit_width(
            self.estimate(1.0)
        )

    def test_the_line_is_lorentzian_rather_than_gaussian(self):
        """Measured in the wings, where the two shapes differ by orders of magnitude.

        Four half-widths out, a Lorentzian sits at 1/17 of its peak and a Gaussian at
        3e-4 of it. For INS this is the part that matters most: the wings are where the
        multiphonon background lives.
        """
        tau_ps = 0.3
        sd = self.estimate(tau_ps)
        pdos = sd.pdos().mean(axis=0)
        gamma = HBAR / tau_ps
        offset = 4 * gamma
        peak = np.interp(self.CENTRE, sd.frequencies, pdos)

        wing = np.interp(self.CENTRE + offset, sd.frequencies, pdos)
        lorentzian = peak / (1 + (offset / gamma) ** 2)
        gaussian = peak * np.exp(-0.5 * (offset / gamma) ** 2)

        assert wing == pytest.approx(lorentzian, rel=0.6)
        assert wing > 100 * gaussian


class TestEnsembleLabel:
    """The one piece of metadata nothing can infer."""

    def test_is_recorded_when_given(self):
        traj = make_trajectory(oscillating(40.0))
        sd = velocity_spectral_density(traj, e_max=200.0, ensemble="NVT")
        assert sd.metadata.ensemble == "NVT"

    def test_is_not_guessed(self):
        traj = make_trajectory(oscillating(40.0))
        assert velocity_spectral_density(traj, e_max=200.0).metadata.ensemble is None

    def test_survives_a_round_trip(self, tmp_path):
        traj = make_trajectory(oscillating(40.0))
        sd = velocity_spectral_density(traj, e_max=200.0, ensemble="NVE")
        path = tmp_path / "ir.h5"
        sd.to_hdf5(path)
        from mdins.ir import VelocitySpectralDensity

        assert VelocitySpectralDensity.from_hdf5(path).metadata.ensemble == "NVE"
