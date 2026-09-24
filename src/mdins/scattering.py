# SPDX-License-Identifier: GPL-3.0-or-later
"""Stages D-E: displacements and scattering intensity.

Not implemented. These entry points exist so that the instrument stage [F] and the test
scaffolding can be written against a stable signature while stages D-E are in progress
(plan.md §4, M0). Both consume the single intermediate representation produced by stage
B and never touch a trajectory.

The two methods differ only in what they do with the tensor:

``isotropic_spectrum``
    Method 1 (Cheng et al. 2020). Uses ``pdos = tr P / 3`` and discards the off-diagonal
    structure, a fully isotropic approximation.

``anisotropic_spectrum``
    Method 2 (Harrelson et al. 2021). Retains the cross-correlation terms, contracting
    the tensor against **Q**. This is why the IR stores the full ``(3, 3)`` tensor
    rather than only its trace.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    import numpy as np
    from numpy.typing import NDArray

    from mdins.ir import VelocitySpectralDensity

__all__ = ["anisotropic_spectrum", "isotropic_spectrum"]

_NOT_YET = (
    "not implemented yet: stages D-E land in milestones M3 (isotropic) and M4 "
    "(anisotropic). See plan.md §7. Stages A-C are usable now via "
    "mdins.spectral.velocity_spectral_density."
)


def isotropic_spectrum(
    density: VelocitySpectralDensity,
    *,
    temperature: float | None = None,
    max_order: int = 10,
) -> NDArray[np.float64]:
    """Method 1: INS intensity in the isotropic approximation.

    Args:
        density: The intermediate representation from stage B.
        temperature: Temperature to evaluate the spectrum at, in K. Distinct from
            ``density.temperature_md``, which is a property of the trajectory; defaults
            to it (design.md D8).
        max_order: Highest phonon order to include in the iterative convolution.

    Raises:
        NotImplementedError: Always, for now.
    """
    raise NotImplementedError(_NOT_YET)


def anisotropic_spectrum(
    density: VelocitySpectralDensity,
    *,
    temperature: float | None = None,
    max_order: int = 10,
) -> NDArray[np.float64]:
    """Method 2: INS intensity retaining the cross-correlation terms.

    Args:
        density: The intermediate representation from stage B. The full ``(3, 3)``
            tensor is used, not only its trace.
        temperature: Temperature to evaluate the spectrum at, in K.
        max_order: Highest phonon order to include in the iterative convolution.

    Raises:
        NotImplementedError: Always, for now.
    """
    raise NotImplementedError(_NOT_YET)
