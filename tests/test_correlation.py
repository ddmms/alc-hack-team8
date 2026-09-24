"""Tests for the VACF and atom-projected pDOS correlation engine."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from md_ins.correlation import (
    H_MEV_FS,
    _apply_window,
    calculate_pdos,
    calculate_vacf,
    invcm_to_mev,
    mev_to_invcm,
    mev_to_thz,
    thz_to_mev,
)


def test_unit_conversions_round_trip() -> None:
    """Energy conversion helpers round-trip correctly."""
    energies = np.array([0.0, 5.0, 8.0])
    np.testing.assert_allclose(mev_to_invcm(energies) / 8.0655440045, energies)
    np.testing.assert_allclose(mev_to_thz(energies) / 0.2417990502, energies)
    np.testing.assert_allclose(invcm_to_mev(mev_to_invcm(energies)), energies)
    np.testing.assert_allclose(thz_to_mev(mev_to_thz(energies)), energies)


def _oscillator_velocities(
    omega_per_fs: float, n_steps: int, timestep_fs: float, n_atoms: int = 1
) -> np.ndarray:
    """Deterministic single-frequency oscillator velocities (all atoms equal)."""
    times = np.arange(n_steps) * timestep_fs
    component = np.cos(omega_per_fs * times)
    velocities = np.zeros((n_steps, n_atoms, 3))
    velocities[:, :, :] = component[:, None, None]
    return velocities


def test_vacf_zero_lag_is_one() -> None:
    """Task 3.1: the normalised VACF equals 1.0 at zero lag."""
    velocities = _oscillator_velocities(
        omega_per_fs=0.01, n_steps=2048, timestep_fs=1.0
    )
    vacf = calculate_vacf(velocities)
    assert vacf.shape == (2048, 1)
    np.testing.assert_allclose(vacf[0, 0], 1.0, atol=1e-10)


def test_vacf_smooth_oscillatory_decay() -> None:
    """Task 3.1: the VACF of a single-frequency oscillator is a smooth cosine."""
    omega = 0.02  # rad / fs
    n_steps = 4096
    timestep_fs = 1.0
    velocities = _oscillator_velocities(omega, n_steps, timestep_fs)
    vacf = calculate_vacf(velocities)[:, 0]
    # Normalised VACF ~ cos(omega * tau) for an oscillator
    # check the first several extrema are smooth and bounded
    assert np.all(np.abs(vacf) <= 1.0 + 1e-9)
    # zero crossing near quarter period
    quarter = int(round(np.pi / (2 * omega * timestep_fs)))
    assert abs(vacf[quarter]) < 0.05
    # near half period it should be negative
    half = int(round(np.pi / (omega * timestep_fs)))
    assert vacf[half] < -0.8


def test_window_endpoints_and_low_frequency_preservation() -> None:
    """Task 3.2: Hann/Blackman taper to zero at the end, rectangular does not,
    and low-frequency content is preserved by the Hann window."""
    n = 1024
    vacf = np.ones((n, 2))
    hann = _apply_window(vacf, "hann")
    blackman = _apply_window(vacf, "blackman")
    rectangular = _apply_window(vacf, "rectangular")
    assert hann[0, 0] == pytest.approx(1.0)
    assert hann[-1, 0] == pytest.approx(0.0, abs=1e-12)
    assert blackman[0, 0] == pytest.approx(1.0)
    assert blackman[-1, 0] == pytest.approx(0.0, abs=1e-12)
    assert np.allclose(rectangular, vacf)


def test_hann_resolves_low_frequency_modes() -> None:
    """Task 3.2: the default Hann window resolves a low-frequency (~20 cm^-1)
    mode rather than suppressing it."""
    energy_target_invcm = 20.0
    energy_target_mev = invcm_to_mev(energy_target_invcm)
    frequency_per_fs = energy_target_mev / H_MEV_FS
    omega = 2.0 * np.pi * frequency_per_fs
    n_steps = 8192
    timestep_fs = 1.0
    velocities = _oscillator_velocities(omega, n_steps, timestep_fs)
    result = calculate_pdos(velocities, timestep_fs, window="hann")
    peak_energy = result.energies_mev[np.argmax(result.total_dos)]
    assert peak_energy == pytest.approx(energy_target_mev, rel=0.05)


def test_pdos_normalises_to_one() -> None:
    """Task 3.3: the per-atom pDOS integrates to 1.0."""
    omega = 0.01
    n_steps = 4096
    timestep_fs = 1.0
    velocities = _oscillator_velocities(omega, n_steps, timestep_fs, n_atoms=3)
    result = calculate_pdos(
        velocities, timestep_fs, symbols=["Ar", "Ar", "Ar"], window="hann"
    )
    for atom_idx in range(3):
        integral = np.trapezoid(result.atom_dos[:, atom_idx], result.energies_mev)
        assert integral == pytest.approx(1.0, abs=0.01)
    total_integral = np.trapezoid(result.total_dos, result.energies_mev)
    assert total_integral == pytest.approx(1.0, abs=0.01)


def test_pdos_multi_frequency_peaks() -> None:
    """Task 3.3: a multi-frequency trajectory produces peaks at the right
    energies."""
    e1_mev = 4.0
    e2_mev = 7.5
    f1 = e1_mev / H_MEV_FS
    f2 = e2_mev / H_MEV_FS
    n_steps = 16384
    timestep_fs = 1.0
    times = np.arange(n_steps) * timestep_fs
    signal = np.cos(2.0 * np.pi * f1 * times) + 0.7 * np.cos(2.0 * np.pi * f2 * times)
    velocities = np.zeros((n_steps, 1, 3))
    velocities[:, 0, :] = signal[:, None]
    result = calculate_pdos(velocities, timestep_fs, window="hann")
    # locate the two strongest peaks
    spectrum = result.total_dos
    # split into low (<5.75) and high (>5.75) regions and find peaks
    energies = result.energies_mev
    low_mask = energies < 5.75
    high_mask = energies > 5.75
    low_peak = energies[low_mask][np.argmax(spectrum[low_mask])]
    high_peak = energies[high_mask][np.argmax(spectrum[high_mask])]
    resolution = H_MEV_FS / ((n_steps - 1) * timestep_fs)
    assert low_peak == pytest.approx(e1_mev, abs=max(0.2, 3 * resolution))
    assert high_peak == pytest.approx(e2_mev, abs=max(0.2, 3 * resolution))


def test_frequency_resolution_warning() -> None:
    """Task 3.4: requesting finer resolution than h/T_tot warns the user."""
    n_steps = 256
    timestep_fs = 1.0
    total_duration = (n_steps - 1) * timestep_fs
    physical_resolution = H_MEV_FS / total_duration
    velocities = _oscillator_velocities(0.01, n_steps, timestep_fs)
    with pytest.warns(UserWarning, match="finer than the physical limit"):
        calculate_pdos(
            velocities, timestep_fs, energy_resolution_mev=physical_resolution * 0.1
        )
    # requesting coarser resolution does not warn
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        calculate_pdos(
            velocities, timestep_fs, energy_resolution_mev=physical_resolution * 5
        )
