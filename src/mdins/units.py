# SPDX-License-Identifier: GPL-3.0-or-later
"""Canonical internal units and the conversions between them.

The internal unit system is **Å, ps, amu, meV**.

Two rules, which exist because this is where calculations of this kind most often go
quietly wrong:

1. Every frequency that crosses a function boundary, enters a dataclass, or is written
   to a file is an **energy in meV**. Angular frequency in rad/ps may be used inside a
   function body; it must never appear in a public signature. This makes most
   factor-of-2π errors unrepresentable rather than merely testable.
2. Conversions live here and nowhere else.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import constants

__all__ = [
    "AMU",
    "ANGSTROM",
    "HBAR",
    "KB",
    "KB_MEV",
    "MEV_PER_AMU_ANGSTROM2_PS2",
    "PICOSECOND",
    "PLANCK",
    "energy_to_angular_frequency",
    "energy_to_frequency",
    "energy_to_wavenumber",
    "frequency_to_energy",
    "thermal_velocity_squared",
    "wavenumber_to_energy",
]

#: Length of the canonical length unit, in metres.
ANGSTROM = 1e-10
#: Length of the canonical time unit, in seconds.
PICOSECOND = 1e-12
#: Mass of the canonical mass unit, in kilograms.
AMU = constants.value("atomic mass constant")

_MEV_IN_JOULES = constants.eV * 1e-3
_VELOCITY_SI = ANGSTROM / PICOSECOND  # 1 Å/ps in m/s

#: Energy of 1 amu·Å²/ps², in meV. Equivalently, the factor converting an energy
#: expressed in the canonical mass/length/time units into meV.
MEV_PER_AMU_ANGSTROM2_PS2 = AMU * _VELOCITY_SI**2 / _MEV_IN_JOULES

#: Reduced Planck constant, in meV·ps.
HBAR = constants.hbar / _MEV_IN_JOULES / PICOSECOND
#: Planck constant, in meV·ps. Energy in meV divided by this is a frequency in 1/ps.
PLANCK = constants.h / _MEV_IN_JOULES / PICOSECOND

#: Boltzmann constant, in meV/K.
KB_MEV = constants.value("Boltzmann constant in eV/K") * 1e3
#: Boltzmann constant, in amu·Å²·ps⁻²·K⁻¹ — the form needed alongside MD velocities.
KB = KB_MEV / MEV_PER_AMU_ANGSTROM2_PS2

_MEV_PER_WAVENUMBER = (
    constants.h * constants.c * 100.0 / _MEV_IN_JOULES
)  # 1 cm⁻¹ in meV


def energy_to_angular_frequency(energy: ArrayLike) -> NDArray[np.float64]:
    """Convert energy in meV to angular frequency in rad/ps.

    For internal use only: angular frequencies must not cross a public boundary.
    """
    return np.asarray(energy, dtype=np.float64) / HBAR


def angular_frequency_to_energy(omega: ArrayLike) -> NDArray[np.float64]:
    """Convert angular frequency in rad/ps to energy in meV."""
    return np.asarray(omega, dtype=np.float64) * HBAR


def energy_to_frequency(energy: ArrayLike) -> NDArray[np.float64]:
    """Convert energy in meV to frequency in 1/ps (THz)."""
    return np.asarray(energy, dtype=np.float64) / PLANCK


def frequency_to_energy(frequency: ArrayLike) -> NDArray[np.float64]:
    """Convert frequency in 1/ps (THz) to energy in meV."""
    return np.asarray(frequency, dtype=np.float64) * PLANCK


def energy_to_wavenumber(energy: ArrayLike) -> NDArray[np.float64]:
    """Convert energy in meV to wavenumber in cm⁻¹."""
    return np.asarray(energy, dtype=np.float64) / _MEV_PER_WAVENUMBER


def wavenumber_to_energy(wavenumber: ArrayLike) -> NDArray[np.float64]:
    """Convert wavenumber in cm⁻¹ to energy in meV."""
    return np.asarray(wavenumber, dtype=np.float64) * _MEV_PER_WAVENUMBER


def thermal_velocity_squared(
    masses: ArrayLike, temperature: float
) -> NDArray[np.float64]:
    """Classical equipartition value of ``<|v|²>`` in (Å/ps)².

    This is ``3 k_B T / m``, the invariant that anchors the absolute normalisation of
    the spectral density (design.md D6).

    Args:
        masses: Atomic masses in amu.
        temperature: Temperature in K.

    Returns:
        Mean square speed for each mass, in (Å/ps)².
    """
    return 3.0 * KB * float(temperature) / np.asarray(masses, dtype=np.float64)
