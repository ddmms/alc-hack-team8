# SPDX-License-Identifier: GPL-3.0-or-later
"""Inelastic neutron scattering spectra from molecular dynamics trajectories.

Implements two published methods that route from an MD trajectory to a simulated INS
spectrum via a derived phonon density of states. See ``design.md`` for the architecture
and ``docs/method-review.md`` for the source-paper analysis.
"""

from __future__ import annotations

from mdins.ir import NORMALISATION, SpectralMetadata, VelocitySpectralDensity
from mdins.provenance import Provenance
from mdins.scattering import anisotropic_spectrum, isotropic_spectrum
from mdins.spectral import velocity_spectral_density
from mdins.trajectory import VelocityTrajectory, read_velocities

__all__ = [
    "NORMALISATION",
    "Provenance",
    "SpectralMetadata",
    "VelocitySpectralDensity",
    "VelocityTrajectory",
    "__version__",
    "anisotropic_spectrum",
    "isotropic_spectrum",
    "read_velocities",
    "velocity_spectral_density",
]

__version__ = "0.1.0.dev0"
