"""Velocity autocorrelation and atom-projected phonon density of states.

The correlation engine computes the velocity autocorrelation function (VACF)
via the Wiener-Khinchin FFT theorem and the atom-/species-projected phonon
density of states (pDOS) :math:`g_d(\\omega)` by Fourier transformation of the
windowed VACF.  Frequencies are handled internally in meV (standard for INS)
with conversion helpers to cm^-1 and THz.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
from scipy import fft

# Planck constant in meV * fs.  One inverse-femtosecond equals 1000 THz, and
# 1 THz = 4.135667697 meV, so 1/fs = 4135.667697 meV.  E[meV] = h * f[1/fs].
H_MEV_FS: float = 4135.667697

# Energy unit conversions (energy <-> frequency).
MEV_TO_INVCM: float = 8.0655440045  # 1 meV  -> cm^-1
MEV_TO_THZ: float = 0.2417990502  # 1 meV  -> THz
INVCM_TO_MEV: float = 1.0 / MEV_TO_INVCM  # 1 cm^-1 -> meV
THZ_TO_MEV: float = 1.0 / MEV_TO_THZ  # 1 THz   -> meV


def mev_to_invcm(energies_mev: np.ndarray | float) -> np.ndarray | float:
    """Convert energies from meV to cm^-1."""
    return np.asarray(energies_mev) * MEV_TO_INVCM


def mev_to_thz(energies_mev: np.ndarray | float) -> np.ndarray | float:
    """Convert energies from meV to THz."""
    return np.asarray(energies_mev) * MEV_TO_THZ


def invcm_to_mev(energies_invcm: np.ndarray | float) -> np.ndarray | float:
    """Convert energies from cm^-1 to meV."""
    return np.asarray(energies_invcm) * INVCM_TO_MEV


def thz_to_mev(energies_thz: np.ndarray | float) -> np.ndarray | float:
    """Convert energies from THz to meV."""
    return np.asarray(energies_thz) * THZ_TO_MEV


def _next_fast_len(n: int) -> int:
    """Return a fast FFT length ``>= 4 * n`` for a smooth frequency grid."""
    target = max(4 * n, 16)
    return fft.next_fast_len(target)


def _correlation_lag0(velocities: np.ndarray) -> np.ndarray:
    """Return ``<|v(t)|^2>`` per atom, shape ``(n_atoms,)``.

    This is the zero-lag (normalisation) value of the per-atom VACF.
    """
    speed_sq = np.sum(velocities**2, axis=2)  # (n_steps, n_atoms)
    return np.mean(speed_sq, axis=0)


def calculate_vacf(velocities: np.ndarray) -> np.ndarray:
    """Compute the per-atom normalised velocity autocorrelation function.

    Uses the Wiener-Khinchin theorem: the autocorrelation is the inverse FFT
    of the squared magnitude of the velocity FFT, with zero-padding to avoid
    circular wrap-around.

    Parameters
    ----------
    velocities : np.ndarray
        Shape ``(n_steps, n_atoms, 3)`` in Angstrom / fs (any consistent unit
        works; the VACF is normalised).

    Returns
    -------
    np.ndarray
        Normalised per-atom VACF of shape ``(n_steps, n_atoms)`` with
        ``vacf[0] == 1.0`` for every atom.
    """
    velocities = np.asarray(velocities, dtype=float)
    if velocities.ndim != 3 or velocities.shape[2] != 3:
        raise ValueError(
            f"velocities must have shape (n_steps, n_atoms, 3), got {velocities.shape}."
        )
    n_steps = velocities.shape[0]
    n_fft = 2 * n_steps

    velocity_fft = fft.fft(velocities, n=n_fft, axis=0)
    power = velocity_fft * np.conj(velocity_fft)
    acf = fft.ifft(power, axis=0).real
    acf = acf[:n_steps]  # positive lags only

    vacf = acf.sum(axis=2)  # dot product over Cartesian, (n_steps, n_atoms)

    normalisation = vacf[0].copy()
    # Guard against stationary atoms producing a zero denominator.
    normalisation[normalisation == 0.0] = 1.0
    return vacf / normalisation


def _apply_window(vacf: np.ndarray, window: str) -> np.ndarray:
    """Apply a window along the time (lag) axis to the per-atom VACF.

    Windows are *one-sided*: equal to 1 at zero lag and tapering to 0 at the
    maximum lag, so the physically important zero-lag value is preserved and
    spectral leakage at long correlation times is suppressed.  This keeps
    low-frequency vibrational modes (down to tens of cm^-1) resolvable.
    """
    n = vacf.shape[0]
    if n == 1:
        return vacf
    index = np.arange(n)
    if window == "rectangular":
        weights = np.ones(n)
    elif window == "hann":
        weights = 0.5 * (1.0 + np.cos(np.pi * index / (n - 1)))
    elif window == "blackman":
        arg = np.pi * index / (n - 1)
        weights = 0.42 + 0.5 * np.cos(arg) + 0.08 * np.cos(2.0 * arg)
    else:
        raise ValueError(
            f"Unknown window {window!r}; use 'hann', 'blackman', or 'rectangular'."
        )
    return vacf * weights[:, None]


def _energy_grid(n_fft: int, timestep_fs: float) -> np.ndarray:
    """Return the positive-frequency energy grid in meV for an FFT of length
    ``n_fft`` over a signal sampled every ``timestep_fs`` femtoseconds.
    """
    frequencies_per_fs = np.arange(n_fft // 2 + 1) / (n_fft * timestep_fs)
    return H_MEV_FS * frequencies_per_fs


@dataclass
class PDOSResult:
    """Result of a pDOS calculation.

    Attributes
    ----------
    energies_mev : np.ndarray
        Energy grid in meV (positive frequencies, zero included).
    atom_dos : np.ndarray
        Per-atom pDOS, shape ``(len(energies_mev), n_atoms)``, each column
        normalised so the integral over energy equals 1.
    species_dos : dict[str, np.ndarray]
        Species-averaged pDOS keyed by chemical symbol, each integrating to 1.
    total_dos : np.ndarray
        pDOS averaged over all atoms, integrating to 1.
    window : str
        Window function used.
    timestep_fs : float
        Sampling timestep.
    n_steps : int
        Number of frames in the input VACF.
    """

    energies_mev: np.ndarray
    atom_dos: np.ndarray
    species_dos: dict[str, np.ndarray]
    total_dos: np.ndarray
    window: str
    timestep_fs: float
    n_steps: int


def calculate_pdos(
    velocities: np.ndarray,
    timestep_fs: float,
    symbols: list[str] | None = None,
    *,
    window: str = "hann",
    n_fft: int | None = None,
    energy_resolution_mev: float | None = None,
) -> PDOSResult:
    """Compute the atom- and species-projected phonon density of states.

    The pDOS is the Fourier transform of the windowed, normalised VACF, taken
    as the (real) cosine transform.  Each per-atom spectrum is normalised so
    that the integral over positive frequencies equals 1.0.

    Parameters
    ----------
    velocities : np.ndarray
        Shape ``(n_steps, n_atoms, 3)``.
    timestep_fs : float
        Sampling timestep in femtoseconds.
    symbols : list[str] or None
        Chemical symbol per atom for species projection.  Required for
        ``species_dos`` to be populated.
    window : str, default 'hann'
        Window applied to the VACF: ``'hann'``, ``'blackman'``, or
        ``'rectangular'``.
    n_fft : int or None
        FFT length (zero-padding) for the frequency grid.  Defaults to a fast
        length ``>= 4 * n_steps``.
    energy_resolution_mev : float or None
        Requested energy resolution in meV.  If finer than the physical limit
        ``h / T_tot`` a :class:`UserWarning` is issued.

    Returns
    -------
    PDOSResult
    """
    velocities = np.asarray(velocities, dtype=float)
    n_steps = velocities.shape[0]
    total_duration_fs = (n_steps - 1) * timestep_fs

    if energy_resolution_mev is not None:
        physical_resolution = (
            H_MEV_FS / total_duration_fs if total_duration_fs > 0 else np.inf
        )
        if energy_resolution_mev < physical_resolution:
            warnings.warn(
                f"Requested energy resolution {energy_resolution_mev:.4g} meV "
                f"is finer than the physical limit h/T_tot = "
                f"{physical_resolution:.4g} meV set by the trajectory duration "
                f"({total_duration_fs:.2f} fs). The spectrum will be "
                "interpolated, not truly resolved.",
                UserWarning,
                stacklevel=2,
            )

    vacf = calculate_vacf(velocities)  # (n_steps, n_atoms)
    vacf_windowed = _apply_window(vacf, window)

    length = n_fft if n_fft is not None else _next_fast_len(n_steps)
    spectrum = fft.rfft(vacf_windowed, n=length, axis=0).real  # cosine transform
    energies = _energy_grid(length, timestep_fs)

    # Normalize each per-atom spectrum so the integral over energy equals 1.
    integrals = np.trapezoid(spectrum, energies, axis=0)
    integrals = np.where(integrals == 0.0, 1.0, integrals)
    atom_dos = spectrum / integrals

    total_dos = np.mean(atom_dos, axis=1)

    species_dos: dict[str, np.ndarray] = {}
    if symbols is not None:
        unique_symbols: list[str] = []
        for sym in symbols:
            if sym not in unique_symbols:
                unique_symbols.append(sym)
        for sym in unique_symbols:
            mask = np.array([s == sym for s in symbols], dtype=bool)
            species_dos[sym] = np.mean(atom_dos[:, mask], axis=1)

    return PDOSResult(
        energies_mev=energies,
        atom_dos=atom_dos,
        species_dos=species_dos,
        total_dos=total_dos,
        window=window,
        timestep_fs=timestep_fs,
        n_steps=n_steps,
    )
