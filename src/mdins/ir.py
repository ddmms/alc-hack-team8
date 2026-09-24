# SPDX-License-Identifier: GPL-3.0-or-later
"""The intermediate representation: velocity spectral density.

This module is the contract between the trajectory side of the pipeline and the
scattering side. Stages A–B produce a :class:`VelocitySpectralDensity`; stages D–G
consume one and never touch a trajectory.

The stored quantity is the one-sided velocity cross-spectral density tensor

.. math::

    P_{i,\\alpha\\beta}(E) \\;=\\;
    \\langle v_{i,\\alpha}^{*}(E)\\, v_{i,\\beta}(E)\\rangle

in **absolute** units of Å²·ps⁻²·meV⁻¹, per atom, with energy transfer in meV. Absolute
rather than a normalised lineshape (design.md D6), so the equipartition sum rule

.. math::

    \\int \\mathrm{tr}\\, P_i(E)\\, \\mathrm{d}E \\;=\\; 3 k_B T / m_i

applies directly to the stored array and is checkable without external context.

The scalar atom-projected DOS used by the isotropic method is the trace of this tensor,
so there is one intermediate representation and two consumers, not two pipelines.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
from numpy.typing import NDArray

from mdins.provenance import Provenance
from mdins.units import thermal_velocity_squared

if TYPE_CHECKING:  # pragma: no cover
    from os import PathLike

__all__ = [
    "NORMALISATION",
    "SpectralMetadata",
    "VelocitySpectralDensity",
    "bin_edges",
    "bin_widths",
]

#: The only normalisation convention this package writes. Recorded in every artefact so
#: that the convention travels with the data.
NORMALISATION = "absolute-spectral-density"

#: Voigt ordering used when packing the symmetric tensor for storage.
_VOIGT = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))

Estimator = Literal["welch", "vacf"]


def bin_edges(centres: NDArray[np.float64]) -> NDArray[np.float64]:
    """Edges of the bins whose centres are ``centres``.

    The end bins are given the same width as their neighbour, so the covered range
    extends half a bin beyond the first and last centre.
    """
    if centres.size < 2:
        raise ValueError("need at least two bin centres to infer edges")
    mid = 0.5 * (centres[:-1] + centres[1:])
    return np.concatenate(
        [
            [centres[0] - (mid[0] - centres[0])],
            mid,
            [centres[-1] + (centres[-1] - mid[-1])],
        ]
    )


def bin_widths(centres: NDArray[np.float64]) -> NDArray[np.float64]:
    """Width of each bin, shape ``(n_freq,)``.

    Every integral over the energy axis in this package uses these weights, because
    that is the rule :func:`mdins.spectral._rebin` conserves. Integrating the same array
    with the trapezoid rule instead loses half a bin at each end of the grid, a
    ``1/(2 (n_freq - 1))`` bias that shows up as a grid-dependent sum rule.
    """
    return np.diff(bin_edges(centres))


@dataclass(frozen=True)
class SpectralMetadata:
    """How a spectral density was estimated.

    Attributes:
        estimator: ``"welch"`` or ``"vacf"``.
        dt: Interval between trajectory frames used, in ps.
        n_frames: Number of frames consumed.
        window: Window function name.
        segment_length: Welch segment length in frames, or ``None`` for the VACF route.
        n_segments: Number of Welch segments averaged.
        overlap: Fractional segment overlap.
        max_lag: VACF maximum lag in frames, or ``None`` for the Welch route.
        energy_resolution: The resolution the estimator actually achieved, in meV. Set
            by the Welch segment length or the VACF maximum lag, **not** by the output
            bin width and **not** by the trajectory duration. A finer output grid than
            this is interpolation, so any downstream comparison of lineshapes has to be
            made at this scale or coarser.
        ensemble: Ensemble the trajectory was generated in, if the caller said so
            (``"NVE"``, ``"NVT"``, ...). Never inferred — no trajectory format records
            it — and carried because it qualifies the result: a thermostat coupled
            tightly enough to control the temperature also perturbs the dynamics being
            measured.
        normalisation: Convention identifier; see :data:`NORMALISATION`.
        provenance: What produced this artefact.
    """

    estimator: Estimator
    dt: float
    n_frames: int
    window: str = "hann"
    segment_length: int | None = None
    n_segments: int | None = None
    overlap: float | None = None
    max_lag: int | None = None
    energy_resolution: float | None = None
    ensemble: str | None = None
    normalisation: str = NORMALISATION
    provenance: Provenance = field(default_factory=lambda: Provenance(source="unknown"))


@dataclass(frozen=True, eq=False)
class VelocitySpectralDensity:
    """Per-atom velocity cross-spectral density — the intermediate representation.

    Attributes:
        frequencies: Energy transfer grid in meV, ascending and non-negative, shape
            ``(n_freq,)``.
        density: Spectral density tensor in Å²·ps⁻²·meV⁻¹, shape
            ``(n_entity, n_freq, 3, 3)``, real and symmetric in the trailing axes.
        symbols: Chemical symbol of each entity, length ``n_entity``.
        masses: Mass of each entity in amu, shape ``(n_entity,)``.
        temperature_md: Temperature the trajectory sampled the potential energy surface
            at, in K. This is a property of the trajectory, distinct from the
            temperature a spectrum is later evaluated at (design.md D8).
        metadata: How the estimate was made.
        pdos_std: Optional standard error of ``pdos()``, shape ``(n_entity, n_freq)``,
            from the spread across Welch segments. Validation tolerances are derived
            from this rather than chosen by hand.
        entity_counts: Number of atoms each entity stands for, shape ``(n_entity,)``.
            All ones for a per-atom density; set by :meth:`group_by_species`, where the
            density is a per-atom *average* over the group. Carrying the counts is what
            keeps :meth:`total_pdos` extensive after grouping.
        com_projected_mass: Total mass of the system whose centre-of-mass velocity was
            projected out, in amu, or ``None`` if it was not. The projection changes the
            equipartition sum rule by a known amount; see :meth:`check_sum_rule`.

    The trailing ``(3, 3)`` axes are stored in full rather than as six symmetric
    components. That costs 1.5x the memory and buys the ability to write
    ``np.einsum("...ij,i,j->...", density, q, q)`` directly, with no packing convention
    to get wrong. Storage on disk packs to six components.
    """

    frequencies: NDArray[np.float64]
    density: NDArray[np.float64]
    symbols: list[str]
    masses: NDArray[np.float64]
    temperature_md: float
    metadata: SpectralMetadata
    pdos_std: NDArray[np.float64] | None = None
    entity_counts: NDArray[np.int64] | None = None
    com_projected_mass: float | None = None

    def __post_init__(self) -> None:
        freqs = np.ascontiguousarray(self.frequencies, dtype=np.float64)
        density = np.ascontiguousarray(self.density, dtype=np.float64)
        masses = np.ascontiguousarray(self.masses, dtype=np.float64)

        if freqs.ndim != 1:
            raise ValueError(f"frequencies must be 1-D, got shape {freqs.shape}")
        if np.any(freqs < 0.0):
            raise ValueError(
                "frequencies must be non-negative: the spectral density is one-sided"
            )
        if np.any(np.diff(freqs) <= 0.0):
            raise ValueError("frequencies must be strictly ascending")

        expected = (len(masses), freqs.size, 3, 3)
        if density.shape != expected:
            raise ValueError(
                f"density has shape {density.shape}, expected {expected} "
                "(n_entity, n_freq, 3, 3)"
            )
        if len(self.symbols) != len(masses):
            raise ValueError(
                f"got {len(self.symbols)} symbols for {len(masses)} masses"
            )
        if np.any(masses <= 0.0):
            raise ValueError("masses must be positive")
        # ``not > 0`` rather than ``<= 0``: NaN fails every comparison, so the obvious
        # spelling would let a NaN temperature through, and a NaN then disables the sum
        # rule silently instead of failing it.
        if not self.temperature_md > 0.0:
            raise ValueError(
                f"temperature_md must be positive and finite, got {self.temperature_md}"
            )
        if not np.isfinite(self.temperature_md):
            raise ValueError("temperature_md must be finite")

        asymmetry = np.abs(density - np.swapaxes(density, -1, -2)).max(initial=0.0)
        scale = np.abs(density).max(initial=0.0)
        if scale > 0.0 and asymmetry > 1e-10 * scale:
            raise ValueError(
                "density must be symmetric in its trailing axes; "
                f"maximum asymmetry {asymmetry:.3e} against scale {scale:.3e}"
            )

        if self.pdos_std is not None:
            std = np.ascontiguousarray(self.pdos_std, dtype=np.float64)
            if std.shape != (len(masses), freqs.size):
                raise ValueError(
                    f"pdos_std has shape {std.shape}, "
                    f"expected {(len(masses), freqs.size)}"
                )
            object.__setattr__(self, "pdos_std", std)

        if self.entity_counts is not None:
            counts = np.ascontiguousarray(self.entity_counts, dtype=np.int64)
            if counts.shape != (len(masses),):
                raise ValueError(
                    f"entity_counts has shape {counts.shape}, expected {(len(masses),)}"
                )
            if np.any(counts < 1):
                raise ValueError("entity_counts must be positive")
            object.__setattr__(self, "entity_counts", counts)

        if self.com_projected_mass is not None:
            if not self.com_projected_mass > 0.0:
                raise ValueError(
                    f"com_projected_mass must be positive and finite, got "
                    f"{self.com_projected_mass}"
                )
            if self.com_projected_mass < masses.max():
                raise ValueError(
                    f"com_projected_mass {self.com_projected_mass} amu is less than "
                    f"the heaviest entity ({masses.max()} amu), which would make the "
                    "centre-of-mass correction to the sum rule negative. The projected "
                    "mass is the mass of the whole system and cannot be smaller than "
                    "one of its parts."
                )

        object.__setattr__(self, "frequencies", freqs)
        object.__setattr__(self, "density", density)
        object.__setattr__(self, "masses", masses)
        object.__setattr__(self, "symbols", list(self.symbols))

    # -- shape -------------------------------------------------------------------

    @property
    def n_entity(self) -> int:
        """Number of atoms (or groups) the density is resolved over."""
        return self.density.shape[0]

    @property
    def n_freq(self) -> int:
        """Number of energy bins."""
        return self.frequencies.size

    # -- derived quantities ------------------------------------------------------

    def pdos(self) -> NDArray[np.float64]:
        """Atom-projected density of states, shape ``(n_entity, n_freq)``.

        This is ``tr P_i(E) / 3`` — the scalar quantity the isotropic method uses in
        place of the displacement tensors. Units Å²·ps⁻²·meV⁻¹.
        """
        return np.trace(self.density, axis1=-2, axis2=-1) / 3.0

    @property
    def counts(self) -> NDArray[np.int64]:
        """Atoms represented by each entity, defaulting to one each."""
        if self.entity_counts is None:
            return np.ones(self.n_entity, dtype=np.int64)
        return self.entity_counts

    def total_pdos(self) -> NDArray[np.float64]:
        """Density of states of the whole system, shape ``(n_freq,)``.

        Weighted by :attr:`counts`, so this stays extensive across
        :meth:`group_by_species`: grouping stores a per-atom average and the count
        restores the total. Summing :meth:`pdos` directly would not.
        """
        return (self.pdos() * self.counts[:, None]).sum(axis=0)

    def mean_square_velocity(self) -> NDArray[np.float64]:
        """Integrated ``<|v_i|²>`` in (Å/ps)², shape ``(n_entity,)``.

        The left-hand side of the equipartition sum rule. Integrated with
        :func:`bin_widths`, matching the rule the rebinning conserves.
        """
        trace = np.trace(self.density, axis1=-2, axis2=-1)
        return (trace * bin_widths(self.frequencies)).sum(axis=-1)

    def eigenvalues(self) -> NDArray[np.float64]:
        """Ascending eigenvalues of each tensor, shape ``(n_entity, n_freq, 3)``."""
        return np.linalg.eigvalsh(self.density)

    # -- invariants --------------------------------------------------------------

    def expected_mean_square_velocity(self) -> NDArray[np.float64]:
        """Right-hand side of the equipartition sum rule, shape ``(n_entity,)``.

        Classical equipartition gives ``<|v_i|²> = 3 k_B T / m_i``. If the
        centre-of-mass velocity was projected out, the sampled velocities are
        ``v_i - Σ_j m_j v_j / M`` and the expectation shifts by a known amount. At
        equilibrium the momentum distribution factorises across atoms whatever the
        potential, so

        .. math::

            \\langle |v_i - V|^2 \\rangle
            = 3 k_B T \\left( \\frac{1}{m_i} - \\frac{1}{M} \\right)

        exactly, using ``<v_i·V> = <|V|²> = 3 k_B T / M``. This is not a small
        correction: for equal masses it is a factor ``1 - 1/N``, and for the carbon in
        methane ``1 - m_C/M = 0.251``. Applying it is what lets the sum rule work as a
        normalisation check on a small cell instead of only asymptotically.
        """
        want = thermal_velocity_squared(self.masses, self.temperature_md)
        if self.com_projected_mass is None:
            return want
        return want * (1.0 - self.masses / self.com_projected_mass)

    def check_sum_rule(self, rtol: float = 5e-2) -> None:
        """Assert classical equipartition against
        :meth:`expected_mean_square_velocity`.

        This is the anchor for the absolute normalisation (design.md D6, A1). A failure
        here means the estimator normalisation is wrong, and every intensity downstream
        will be wrong by the same factor.

        Two things this does not correct for. Removing global *rotation* is not a linear
        projection with a mass-independent expectation, so no closed form is applied and
        a rotation-removed trajectory will fail this at the ``3/(3N)`` level. And the
        integral is truncated at the grid's upper edge, so any spectral weight above it
        is missing here as well as everywhere downstream.

        Args:
            rtol: Relative tolerance. The default allows for finite sampling error in a
                real trajectory; tighten it for synthetic tests.

        Raises:
            ValueError: If any entity violates the sum rule beyond ``rtol``.
        """
        got = self.mean_square_velocity()
        want = self.expected_mean_square_velocity()

        # A single atom whose own centre of mass was projected out has no velocity left,
        # so the expectation is exactly zero and a relative tolerance is meaningless.
        # Comparing absolutely is the only thing that can fail here; dividing would give
        # 0/0 and pass whatever the data said.
        degenerate = want <= 0.0
        if degenerate.any() and np.abs(got[degenerate]).max() > 0.0:
            raise ValueError(
                f"{int(degenerate.sum())} entities should carry no velocity at all "
                "once the centre-of-mass projection is accounted for, but their "
                f"spectra integrate to up to {np.abs(got[degenerate]).max():.6g} "
                "(Å/ps)². This means the projection recorded in the metadata is not "
                "the one that was applied."
            )

        with np.errstate(invalid="ignore", divide="ignore"):
            rel = np.where(degenerate, 0.0, np.abs(got - want) / want)
        # A non-finite residual is a failure, not a pass: ``nan > rtol`` is False.
        bad = np.flatnonzero(~(rel <= rtol))
        if bad.size:
            worst = int(bad[np.argmax(rel[bad])])
            raise ValueError(
                f"equipartition sum rule violated for {bad.size} of {self.n_entity} "
                f"entities; worst is entity {worst} ({self.symbols[worst]}): "
                f"got {got[worst]:.6g}, expected {want[worst]:.6g} (Å/ps)², "
                f"relative error {rel[worst]:.3%}"
            )

    def check_positive_semidefinite(self, atol_fraction: float = 1e-8) -> None:
        """Assert every tensor has non-negative eigenvalues.

        A negative eigenvalue is a negative mean-square displacement along some
        direction, which is not a physical tensor. The Welch estimator guarantees this
        by construction and the VACF estimator does not, which is why Welch is the
        default (design.md D5).

        Args:
            atol_fraction: Tolerance as a fraction of the largest eigenvalue present.

        Raises:
            ValueError: If any eigenvalue is negative beyond tolerance.
        """
        eigs = self.eigenvalues()
        scale = float(np.abs(eigs).max(initial=0.0))
        tol = -atol_fraction * scale
        if np.any(eigs < tol):
            worst = float(eigs.min())
            count = int(np.count_nonzero(eigs < tol))
            raise ValueError(
                f"spectral density is not positive semi-definite: {count} negative "
                f"eigenvalues, most negative {worst:.3e} against scale {scale:.3e}"
            )

    # -- reshaping ---------------------------------------------------------------

    def group_by_species(self) -> VelocitySpectralDensity:
        """Average over entities sharing a chemical symbol.

        The stored density is the per-atom *mean* over the group, so the sum rule still
        reads ``3 k_B T / m`` per entity. :attr:`entity_counts` records how many atoms
        each entry stands for, which is what keeps :meth:`total_pdos` extensive.

        A memory escape hatch, not a default: the Debye-Waller factor is genuinely
        per-atom in a disordered system, which is the point of these methods. Grouping
        discards that.

        Raises:
            ValueError: If a group's masses are not all equal. Averaging them would be
                wrong — the sum rule is linear in ``1/m``, not in ``m`` — so an isotope
                mixture has to be separated by the caller rather than silently merged.
        """
        order = sorted(set(self.symbols))
        symbols = np.asarray(self.symbols)
        index = [np.flatnonzero(symbols == s) for s in order]

        for symbol, i in zip(order, index, strict=True):
            if not np.allclose(self.masses[i], self.masses[i][0]):
                raise ValueError(
                    f"entities labelled {symbol!r} have differing masses "
                    f"({np.unique(self.masses[i])}); grouping them would average a "
                    "quantity the sum rule is not linear in. Label the isotopes "
                    "distinctly and group those."
                )

        density = np.stack([self.density[i].mean(axis=0) for i in index])
        masses = np.array([self.masses[i][0] for i in index])
        counts = np.array([self.counts[i].sum() for i in index], dtype=np.int64)
        std = None
        if self.pdos_std is not None:
            std = np.stack(
                [np.sqrt((self.pdos_std[i] ** 2).sum(axis=0)) / len(i) for i in index]
            )
        return replace(
            self,
            density=density,
            symbols=order,
            masses=masses,
            pdos_std=std,
            entity_counts=counts,
            metadata=replace(
                self.metadata,
                provenance=self.metadata.provenance.with_step("group_by_species"),
            ),
        )

    # -- persistence -------------------------------------------------------------

    def to_hdf5(self, path: str | PathLike[str]) -> None:
        """Write to HDF5, packing the symmetric tensor to six Voigt components."""
        import h5py

        packed = np.stack(
            [self.density[..., a, b] for a, b in _VOIGT], axis=-1
        )  # (n_entity, n_freq, 6)

        with h5py.File(Path(path), "w") as handle:
            handle.attrs["format"] = "mdins-velocity-spectral-density"
            handle.attrs["format_version"] = 1
            handle.attrs["voigt_order"] = "xx,yy,zz,yz,xz,xy"

            handle.create_dataset("frequencies", data=self.frequencies)
            handle.create_dataset("density_voigt", data=packed, compression="gzip")
            handle.create_dataset("masses", data=self.masses)
            handle.create_dataset(
                "symbols", data=np.array(self.symbols, dtype=h5py.string_dtype())
            )
            if self.pdos_std is not None:
                handle.create_dataset(
                    "pdos_std", data=self.pdos_std, compression="gzip"
                )
            if self.entity_counts is not None:
                handle.create_dataset("entity_counts", data=self.entity_counts)

            meta = handle.create_group("metadata")
            meta.attrs["temperature_md"] = self.temperature_md
            if self.com_projected_mass is not None:
                meta.attrs["com_projected_mass"] = self.com_projected_mass
            for key, value in vars(self.metadata).items():
                if key == "provenance":
                    meta.attrs["provenance"] = value.to_json()
                elif value is not None:
                    meta.attrs[key] = value
                else:
                    meta.attrs[key] = json.dumps(None)

    @classmethod
    def from_hdf5(cls, path: str | PathLike[str]) -> VelocitySpectralDensity:
        """Read a file written by :meth:`to_hdf5`."""
        import h5py

        with h5py.File(Path(path), "r") as handle:
            fmt = handle.attrs.get("format")
            if fmt != "mdins-velocity-spectral-density":
                raise ValueError(f"not an mdins spectral density file (format={fmt!r})")

            packed = np.asarray(handle["density_voigt"])
            density = np.empty((*packed.shape[:-1], 3, 3), dtype=np.float64)
            for n, (a, b) in enumerate(_VOIGT):
                density[..., a, b] = packed[..., n]
                density[..., b, a] = packed[..., n]

            meta_attrs: dict[str, Any] = dict(handle["metadata"].attrs)
            temperature = float(meta_attrs.pop("temperature_md"))
            provenance = Provenance.from_json(meta_attrs.pop("provenance"))
            projected = meta_attrs.pop("com_projected_mass", None)
            fields: dict[str, Any] = {
                key: (None if value == "null" else value)
                for key, value in meta_attrs.items()
            }
            for key in ("n_frames", "segment_length", "n_segments", "max_lag"):
                value = fields.get(key)
                if value is not None:
                    fields[key] = int(value)

            pdos_std = np.asarray(handle["pdos_std"]) if "pdos_std" in handle else None
            counts = (
                np.asarray(handle["entity_counts"])
                if "entity_counts" in handle
                else None
            )

            return cls(
                frequencies=np.asarray(handle["frequencies"]),
                density=density,
                symbols=[s.decode() for s in handle["symbols"]],
                masses=np.asarray(handle["masses"]),
                temperature_md=temperature,
                metadata=SpectralMetadata(provenance=provenance, **fields),
                pdos_std=pdos_std,
                entity_counts=counts,
                com_projected_mass=None if projected is None else float(projected),
            )
