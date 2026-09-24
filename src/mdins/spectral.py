# SPDX-License-Identifier: GPL-3.0-or-later
"""Stage B: velocity spectral density estimators.

Two estimators for the same quantity, both producing a
:class:`~mdins.ir.VelocitySpectralDensity`:

``welch``
    Segment-averaged cross-periodogram. The default, because it averages outer products
    and so yields a positive semi-definite tensor by construction. A negative eigenvalue
    is a negative mean-square displacement along some direction, which is not a physical
    tensor, and the anisotropic method depends on the tensor being one.

``vacf``
    Windowed velocity autocorrelation, transformed (Blackman-Tukey). What the source
    papers do. Positive semi-definiteness is not guaranteed.

They agree wherever both are valid, and requiring that agreement is an independent check
on this stage. See design.md D5.

The output energy grid is an *input* to estimation, not a post-hoc downsample: each
segment is rebinned onto it as it is computed, so the full-resolution per-atom array is
never accumulated. Without this a 3000-atom system at 25000 native bins needs 5.4 GB
before anything useful has happened.

That bounds what is *kept*, not what is touched: forming one segment's cross-spectrum is
itself an ``(n_native, n_atoms, 3, 3)`` intermediate. Atoms are therefore processed in
chunks sized by :data:`CHUNK_BYTES`, so peak memory is bounded by the chunk rather than
by the system. Resolution is a separate axis from bin count and is recorded as
``metadata.energy_resolution``; a finer output grid than that is interpolation.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import fft as sp_fft
from scipy.integrate import cumulative_trapezoid
from scipy.signal import get_window

from mdins.ir import SpectralMetadata, VelocitySpectralDensity, bin_edges
from mdins.units import PLANCK

if TYPE_CHECKING:  # pragma: no cover
    from mdins.trajectory import VelocityTrajectory

__all__ = ["velocity_spectral_density"]


def _output_grid(
    frequencies: ArrayLike | None,
    e_max: float | None,
    n_bins: int,
    nyquist: float,
) -> NDArray[np.float64]:
    if frequencies is not None:
        grid = np.ascontiguousarray(frequencies, dtype=np.float64)
        if grid.ndim != 1 or grid.size < 2:
            raise ValueError("frequencies must be a 1-D grid with at least two points")
        return grid
    top = nyquist if e_max is None else float(e_max)
    return np.linspace(0.0, top, int(n_bins))


def _rebin(
    values: NDArray[np.float64],
    native: NDArray[np.float64],
    target: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Rebin a spectral density from ``native`` onto ``target`` energies.

    Integral-conserving: interpolates the cumulative integral at the target bin edges
    and differences it, so ``∫ P dE`` survives, which is what the sum rule tests.

    Args:
        values: Shape ``(..., n_native)``, density per unit energy.
        native: Ascending native energy grid.
        target: Ascending target energy grid, treated as bin centres.

    Returns:
        Shape ``(..., n_target)``.
    """
    cumulative = np.concatenate(
        [
            np.zeros((*values.shape[:-1], 1)),
            cumulative_trapezoid(values, x=native, axis=-1),
        ],
        axis=-1,
    )
    edges = bin_edges(target)
    lower, weight = _interpolation_weights(native, edges)
    sampled = (
        cumulative[..., lower] * (1.0 - weight) + cumulative[..., lower + 1] * weight
    )
    return np.diff(sampled, axis=-1) / np.diff(edges)


def _interpolation_weights(
    native: NDArray[np.float64], edges: NDArray[np.float64]
) -> tuple[NDArray[np.intp], NDArray[np.float64]]:
    """Linear interpolation indices and weights, computed once and reused.

    Equivalent to :func:`numpy.interp`, including its clamping outside ``native``, but
    vectorised over leading axes — which matters because rebinning happens once per
    Welch segment.
    """
    lower = np.clip(np.searchsorted(native, edges) - 1, 0, native.size - 2)
    span = native[lower + 1] - native[lower]
    weight = np.clip((edges - native[lower]) / span, 0.0, 1.0)
    return lower, weight


#: Target size of the per-segment working set, in bytes. Peak memory in
#: :func:`_welch` scales with this rather than with the number of atoms.
CHUNK_BYTES = 64 * 1024 * 1024


def _atom_chunks(n_atoms: int, n_native: int) -> list[tuple[int, int]]:
    """Half-open atom ranges whose cross-spectra fit in :data:`CHUNK_BYTES`.

    The intermediate is ``(n_native, n_chunk, 3, 3)`` complex, at 16 bytes per element,
    and roughly three such arrays are live at once (the cross-spectrum, its real view
    and the cumulative integral inside the rebinning).
    """
    per_atom = n_native * 9 * 16 * 3
    size = max(int(CHUNK_BYTES // max(per_atom, 1)), 1)
    return [(lo, min(lo + size, n_atoms)) for lo in range(0, n_atoms, size)]


def _check_resolution(
    achieved: float, requested: float | None, knob: str, value: int
) -> None:
    """Refuse a resolution the estimator cannot deliver.

    design.md §2[A].3 requires validating the duration against the requested resolution
    as well as the dump interval against the requested energy range. For a segmented
    estimator the binding constraint is the *segment*, not the whole trajectory: a 1 ns
    run cut into 1 ps segments resolves 4 meV, not 0.004 meV. Checking the total
    duration would validate a quantity nothing computes.
    """
    if requested is None or achieved <= requested:
        return
    raise ValueError(
        f"requested {requested:.4g} meV resolution but {knob}={value} gives "
        f"{achieved:.4g} meV. Resolution is set by the segment, not by the total "
        f"duration: raise {knob} (at the cost of fewer averages) or run for longer and "
        f"raise it. Increasing the output bin count interpolates and does not help."
    )


def _check_com_removed(trajectory: VelocityTrajectory) -> None:
    """Refuse to estimate from a drifting trajectory.

    Centre-of-mass drift is a spurious zero-energy feature, and the ``1/E`` factor in
    the displacement model amplifies exactly that region. This checks rather than
    silently fixing, so the pre-processing applied stays visible in provenance.
    """
    if any("remove_com_velocity" in step for step in trajectory.provenance.steps):
        return
    masses = trajectory.masses
    com = np.einsum("a,fad->fd", masses, trajectory.velocities) / masses.sum()
    # Net translation, i.e. the time average — not the instantaneous centre-of-mass
    # speed, which fluctuates harmlessly at 1/sqrt(N) of the thermal scale even in a
    # perfectly drift-free run.
    drift = float(np.linalg.norm(com.mean(axis=0)))
    thermal = float(np.sqrt((trajectory.velocities**2).sum(axis=2).mean()))
    if thermal > 0.0 and drift > 0.05 * thermal:
        raise ValueError(
            f"net centre-of-mass motion is {drift:.4g} Å/ps, {drift / thermal:.1%} of "
            "the per-atom thermal speed. Call trajectory.remove_com_velocity() first, "
            "or the drift will appear as a spurious low-energy feature."
        )


def _welch(
    velocities: NDArray[np.float64],
    dt: float,
    grid: NDArray[np.float64],
    segment_length: int,
    overlap: float,
    window: str,
) -> tuple[NDArray[np.float64], NDArray[np.float64], int, float]:
    """Segment-averaged cross-periodogram.

    Atoms are processed in chunks: the per-segment cross-spectrum is an
    ``(n_native, n_chunk, 3, 3)`` complex intermediate, which for a large system is far
    bigger than either the trajectory or the result. Chunking bounds it without changing
    the answer, since atoms are independent here.

    Returns:
        ``(density, trace_std_error, n_segments, quadrature)``. ``density`` is shaped
        ``(n_atoms, n_target, 3, 3)`` in Å²·ps⁻²·meV⁻¹. ``quadrature`` is the largest
        imaginary part of the cross-spectrum relative to the largest real part — the
        component the co-spectrum convention discards (ambiguity A2).
    """
    n_frames, n_atoms, _ = velocities.shape
    step = max(round(segment_length * (1.0 - overlap)), 1)
    starts = range(0, n_frames - segment_length + 1, step)
    if not starts:
        raise ValueError(
            f"segment_length={segment_length} exceeds the {n_frames} frames available"
        )

    taper = get_window(window, segment_length, fftbins=True)
    native = PLANCK * sp_fft.rfftfreq(segment_length, d=dt)

    # One-sided PSD per unit frequency (1/ps), then per unit energy (meV).
    scale = 2.0 * dt / (taper**2).sum() / PLANCK
    one_sided = np.full(native.size, scale)
    one_sided[0] /= 2.0
    if segment_length % 2 == 0:
        one_sided[-1] /= 2.0

    total = np.zeros((n_atoms, grid.size, 3, 3))
    trace_sum = np.zeros((n_atoms, grid.size))
    trace_sq = np.zeros((n_atoms, grid.size))
    n_segments = len(starts)
    real_scale = 0.0
    imaginary_scale = 0.0

    for lo, hi in _atom_chunks(n_atoms, native.size):
        for start in starts:
            block = velocities[start : start + segment_length, lo:hi]
            spectrum = sp_fft.rfft(block * taper[:, None, None], axis=0)
            cross = np.einsum("fai,faj->faij", spectrum.conj(), spectrum)
            cross = np.moveaxis(cross, 0, 1) * one_sided[None, :, None, None]

            # Diagnostic before the discard, so A2 is a decision and not an accident.
            real_scale = max(real_scale, float(np.abs(cross.real).max(initial=0.0)))
            imaginary_scale = max(
                imaginary_scale, float(np.abs(cross.imag).max(initial=0.0))
            )

            rebinned = _rebin(np.moveaxis(cross.real, 1, -1), native, grid)
            rebinned = np.moveaxis(rebinned, -1, 1)

            total[lo:hi] += rebinned
            trace = np.trace(rebinned, axis1=-2, axis2=-1)
            trace_sum[lo:hi] += trace
            trace_sq[lo:hi] += trace**2

    mean = total / n_segments
    if n_segments > 1:
        variance = np.maximum(
            trace_sq / n_segments - (trace_sum / n_segments) ** 2, 0.0
        )
        std_error = np.sqrt(variance / (n_segments - 1))
    else:
        std_error = np.full_like(trace_sum, np.nan)
    quadrature = imaginary_scale / real_scale if real_scale > 0.0 else 0.0
    return mean, std_error / 3.0, n_segments, quadrature


def _vacf(
    velocities: NDArray[np.float64],
    dt: float,
    grid: NDArray[np.float64],
    max_lag: int,
    window: str,
) -> tuple[NDArray[np.float64], float]:
    """Windowed velocity autocorrelation, transformed.

    Returns:
        ``(density, quadrature_fraction)``. The second value is the size of the
        antisymmetric part of the correlation tensor relative to the symmetric part —
        the component discarded by taking the co-spectrum (ambiguity A2). It is a
        diagnostic, not a correction.
    """
    n_frames, n_atoms, _ = velocities.shape
    if max_lag >= n_frames:
        raise ValueError(f"max_lag={max_lag} needs fewer than {n_frames} frames")

    n_fft = sp_fft.next_fast_len(2 * n_frames)
    transform = sp_fft.rfft(velocities, n=n_fft, axis=0)
    # Correlation for lags 0..max_lag, biased (divide by n_frames rather than by the
    # number of overlapping samples). The biased estimator is the one whose transform is
    # non-negative; the unbiased one is noisier at long lag and less well behaved here.
    correlation = (
        sp_fft.irfft(
            np.einsum("fai,faj->faij", transform.conj(), transform), n=n_fft, axis=0
        )[: max_lag + 1]
        / n_frames
    )

    symmetric = 0.5 * (correlation + np.swapaxes(correlation, -1, -2))
    antisymmetric = 0.5 * (correlation - np.swapaxes(correlation, -1, -2))
    scale = float(np.abs(symmetric).max(initial=0.0))
    quadrature = (
        float(np.abs(antisymmetric).max(initial=0.0) / scale) if scale > 0.0 else 0.0
    )

    taper = get_window(window, 2 * max_lag + 1, fftbins=False)[max_lag:]
    tapered = symmetric * taper[:, None, None, None]

    n_out = sp_fft.next_fast_len(2 * max_lag + 1)
    padded = np.zeros((n_out, n_atoms, 3, 3))
    padded[: max_lag + 1] = tapered
    padded[n_out - max_lag :] = tapered[:0:-1]

    native = PLANCK * sp_fft.rfftfreq(n_out, d=dt)
    spectrum = sp_fft.rfft(padded, axis=0).real * (2.0 * dt / PLANCK)
    spectrum[0] /= 2.0
    if n_out % 2 == 0:
        spectrum[-1] /= 2.0

    spectrum = np.moveaxis(spectrum, 0, 1)  # (n_atoms, n_native, 3, 3)
    rebinned = _rebin(np.moveaxis(spectrum, 1, -1), native, grid)
    return np.moveaxis(rebinned, -1, 1), quadrature


def velocity_spectral_density(
    trajectory: VelocityTrajectory,
    *,
    frequencies: ArrayLike | None = None,
    e_max: float | None = None,
    n_bins: int = 1024,
    e_resolution: float | None = None,
    estimator: str = "welch",
    segment_length: int | None = None,
    overlap: float = 0.5,
    window: str = "hann",
    max_lag: int | None = None,
    temperature: float | None = None,
    ensemble: str | None = None,
) -> VelocitySpectralDensity:
    """Estimate the per-atom velocity cross-spectral density.

    Args:
        trajectory: Source velocities. Centre-of-mass drift must already have been
            removed; this is checked, not done silently.
        frequencies: Output energy grid in meV. Defaults to ``n_bins`` points from zero
            to ``e_max``.
        e_max: Upper energy in meV when ``frequencies`` is not given. Defaults to the
            Nyquist energy of the dump interval.
        n_bins: Number of output bins when ``frequencies`` is not given. This sets the
            sampling of the output, not the resolution of the estimate: see
            ``e_resolution``.
        e_resolution: Energy resolution required, in meV. Checked against what the
            estimator can actually deliver — the Welch segment length or the VACF
            maximum lag — and an error if it cannot (design.md §2[A].3). Leave as
            ``None`` to accept whatever the settings give; the achieved value is always
            recorded in ``metadata.energy_resolution``.
        estimator: ``"welch"`` or ``"vacf"``.
        segment_length: Welch segment length in frames. Defaults to the largest power of
            two not exceeding an eighth of the trajectory, which with the default 50%
            overlap gives at least fifteen segments.
        overlap: Welch fractional segment overlap. Note that overlapping segments are
            correlated, so the reported ``pdos_std`` understates the true uncertainty;
            use ``overlap=0.0`` where the error bar matters.
        window: Window function, applied to the data for Welch and to the lag for VACF.
        max_lag: VACF maximum lag in frames. Defaults to an eighth of the trajectory.
        temperature: MD temperature in K. Defaults to the trajectory's kinetic
            temperature.
        ensemble: Ensemble the trajectory was generated in, recorded verbatim in the
            metadata. Never inferred: no trajectory format records it, and it changes
            how the result should be read — a strongly coupled thermostat perturbs the
            dynamics it is sampling, which shows up in the low-energy pDOS, and a
            velocity-rescaling one can break the equipartition the sum rule assumes.

    Returns:
        The intermediate representation, in absolute units.

    Raises:
        ValueError: If the requested grid exceeds the Nyquist energy, if the requested
            resolution is unachievable, if drift has not been removed, or if the
            estimator name is unknown.
    """
    _check_com_removed(trajectory)

    grid = _output_grid(frequencies, e_max, n_bins, trajectory.max_resolvable_energy)
    trajectory.validate_sampling(float(grid[-1]), e_resolution)

    velocities = trajectory.velocities
    provenance = trajectory.provenance
    pdos_std = None
    metadata_extra: dict[str, Any] = {}

    if estimator == "welch":
        if segment_length is None:
            exponent = int(np.floor(np.log2(max(trajectory.n_frames // 8, 2))))
            segment_length = max(2**exponent, 2)
        achieved = PLANCK / (segment_length * trajectory.dt)
        _check_resolution(achieved, e_resolution, "segment_length", segment_length)
        density, pdos_std, n_segments, quadrature = _welch(
            velocities, trajectory.dt, grid, segment_length, overlap, window
        )
        provenance = provenance.with_step(
            f"welch(segment_length={segment_length}, overlap={overlap}, "
            f"window={window!r}, n_segments={n_segments})"
        ).with_note(
            f"co-spectrum retained; discarded imaginary part of the cross-spectrum was "
            f"{quadrature:.2%} of the real part (ambiguity A2)"
        )
        metadata_extra = {
            "segment_length": segment_length,
            "n_segments": n_segments,
            "overlap": overlap,
            "energy_resolution": achieved,
        }
    elif estimator == "vacf":
        if max_lag is None:
            max_lag = max(trajectory.n_frames // 8, 2)
        achieved = PLANCK / (max_lag * trajectory.dt)
        _check_resolution(achieved, e_resolution, "max_lag", max_lag)
        density, quadrature = _vacf(velocities, trajectory.dt, grid, max_lag, window)
        provenance = provenance.with_step(
            f"vacf(max_lag={max_lag}, window={window!r})"
        ).with_note(
            "co-spectrum retained; discarded antisymmetric part was "
            f"{quadrature:.2%} of the symmetric part (ambiguity A2)"
        )
        metadata_extra = {"max_lag": max_lag, "energy_resolution": achieved}
    else:
        raise ValueError(f"unknown estimator {estimator!r}, expected 'welch' or 'vacf'")

    # Symmetric by construction, but floating point leaves a residue that the IR's
    # validation would reject.
    density = 0.5 * (density + np.swapaxes(density, -1, -2))

    return VelocitySpectralDensity(
        frequencies=grid,
        density=density,
        symbols=trajectory.symbols,
        masses=trajectory.masses,
        temperature_md=(
            trajectory.temperature() if temperature is None else float(temperature)
        ),
        metadata=SpectralMetadata(
            estimator=estimator,  # type: ignore[arg-type]
            dt=trajectory.dt,
            n_frames=trajectory.n_frames,
            window=window,
            provenance=provenance,
            ensemble=ensemble,
            **metadata_extra,
        ),
        pdos_std=pdos_std,
        # Recorded so the sum rule can correct for the projection exactly rather than
        # tolerate it approximately.
        com_projected_mass=(trajectory.total_mass if trajectory.com_removed else None),
    )
