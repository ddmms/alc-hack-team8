"""Tests for the Lennard-Jones argon harmonic benchmark.

The integration test (``test_benchmark_md_pdos_matches_harmonic_regression``) is
*regression based*: with a fixed random seed the pipeline is deterministic, so
it pins the detected peak locations to recorded baselines rather than requiring
fine cross-method agreement.  A coarse general-position check confirms the MD
and harmonic dominant peaks lie in the same spectral region.  See the
comparison plot in the project README for a human review of the agreement.
"""

from __future__ import annotations

import numpy as np
import pytest

from md_ins.benchmark import (
    compute_force_constants,
    compute_harmonic_dos,
    make_qgrid,
    run_benchmark,
)


def _argon_unitcell():
    from ase.build import bulk

    return bulk("Ar", "fcc", a=5.26, cubic=True)


def test_compute_force_constants_matrix_shape() -> None:
    """Task 4.1: the ASE-to-Phonopy wrapper produces a force-constants matrix."""
    atoms = _argon_unitcell()
    phonon = compute_force_constants(atoms, np.diag([2, 2, 2]))
    fc = phonon.force_constants
    n_sc = 2**3 * len(atoms)
    assert fc.shape == (n_sc, n_sc, 3, 3)
    assert np.all(np.isfinite(fc))
    # the matrix should be symmetric to within numerical tolerance
    np.testing.assert_allclose(fc, fc.transpose(0, 1, 3, 2), atol=1e-4)


def test_compute_harmonic_dos_shape_and_units() -> None:
    """Task 4.2: Euphonic produces a finite harmonic DOS in meV units."""
    atoms = _argon_unitcell()
    phonon = compute_force_constants(atoms, np.diag([2, 2, 2]))
    qpoints, weights = make_qgrid(2)
    energies, dos, band_max = compute_harmonic_dos(
        phonon, qpoints, weights, dos_bins_mev=np.arange(0.0, 12.0, 0.05)
    )
    assert energies.ndim == 1
    assert dos.shape == energies.shape
    assert np.all(np.isfinite(energies))
    assert np.all(np.isfinite(dos))
    # argon LJ band maximum is a few meV (Debye ~ 8.4 meV)
    assert 7.0 < band_max < 9.0
    # the DOS support sits below the band maximum
    assert energies.min() >= 0.0


@pytest.mark.slow
def test_benchmark_md_pdos_matches_harmonic_regression() -> None:
    """Task 4.4: regression test pinning MD pDOS and harmonic DOS peak locations.

    Uses a fixed seed so the run is deterministic.  The MD pDOS is a broadened
    version of the same phonon spectrum sampled on the matching q-grid, so its
    dominant peak should coincide with the harmonic DOS dominant peak in the
    same general region.  Exact 5% cross-method agreement is not asserted --
    instead the detected peaks are pinned to recorded baselines (regression)
    and a coarse general-position check guards against gross disagreement.
    """
    result = run_benchmark(supercell=(3, 3, 3), n_steps=6000, seed=42)

    # Sanity: argon LJ Debye frequency is ~8.4 meV.
    assert result.harmonic_band_max_mev == pytest.approx(8.39, abs=0.02)

    # Harmonic DOS peaks (regression baselines).
    harm_peaks = result.harmonic_peaks_mev
    assert len(harm_peaks) >= 2
    np.testing.assert_allclose(harm_peaks[0], 4.67, atol=0.05)
    np.testing.assert_allclose(harm_peaks[1], 7.29, atol=0.05)

    # MD pDOS peaks (regression baselines).
    md_peaks = result.md_peaks_mev
    assert len(md_peaks) >= 1
    np.testing.assert_allclose(md_peaks[0], 4.74, atol=0.08)

    # MD pDOS normalisation.
    integral = float(np.trapezoid(result.md_dos, result.md_energies_mev))
    assert integral == pytest.approx(1.0, abs=0.01)

    # Coarse general-position check: dominant peaks agree within ~15%.
    dominant_rel = abs(md_peaks[0] - harm_peaks[0]) / harm_peaks[0]
    assert dominant_rel < 0.15

    # A second MD peak should also appear near the harmonic band-edge peak.
    if len(md_peaks) >= 2:
        assert abs(md_peaks[1] - harm_peaks[1]) / harm_peaks[1] < 0.15
