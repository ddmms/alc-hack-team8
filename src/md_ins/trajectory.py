"""Trajectory data container and ASE-backed ingestion routines.

All quantities are stored in *standard* units:

* velocities  -- Angstrom / femtosecond
* positions   -- Angstrom
* masses      -- atomic mass units (amu)
* timestep    -- femtoseconds
* cell        -- Angstrom (unit-cell row vectors)
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from ase import units
from ase.atoms import Atoms
from ase.io import read

if TYPE_CHECKING:
    from collections.abc import Sequence

# ASE stores velocities in units of Angstrom per ASE-time-unit, where the
# kinetic energy satisfies ``KE = 0.5 * sum(masses * v**2)`` in eV.  One
# femtosecond equals ``ase.units.fs`` ASE-time-units, so multiplying an ASE
# velocity by this factor yields Angstrom / fs.
_ASE_VEL_TO_ANG_PER_FS: float = units.fs

# Conversion from (amu * Angstrom**2 / fs**2) to eV.  Because ASE velocities in
# Ang/ASE-time-unit satisfy ``KE[eV] = 0.5 * m[amu] * v_ASE**2`` and
# ``v[Ang/fs] = v_ASE * units.fs``, kinetic energy in eV is
# ``0.5 * m * v[Ang/fs]**2 / units.fs**2``.
_KE_TO_EV: float = 1.0 / units.fs**2


def _minimum_image(dr: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Apply the minimum-image convention to displacement vectors ``dr``.

    ``dr`` has shape ``(..., 3)`` and ``cell`` is a ``(3, 3)`` matrix of row
    vectors.  Displacements are wrapped into the Wigner-Seitz cell so that
    central differences across periodic boundaries are correct.
    """
    inv_cell = np.linalg.inv(cell)
    frac = dr @ inv_cell
    frac -= np.round(frac)
    return frac @ cell


def _derive_velocities_central(
    positions: np.ndarray, timestep_fs: float, cell: np.ndarray | None
) -> tuple[np.ndarray, np.ndarray]:
    """Derive velocities from positions via central finite differences.

    Uses ``v(t) = (r(t + dt) - r(t - dt)) / (2 dt)`` and drops the first and
    last frames (which lack neighbours on both sides).  For periodic systems
    the minimum-image convention is applied to the displacement.

    Returns ``(velocities, trimmed_positions)`` both of shape
    ``(n_steps - 2, n_atoms, 3)``.
    """
    if positions.shape[0] < 3:
        raise ValueError(
            "Need at least 3 frames to derive velocities by central "
            f"differences, got {positions.shape[0]}."
        )
    dr = positions[2:] - positions[:-2]
    if cell is not None:
        dr = _minimum_image(dr, cell)
    velocities = dr / (2.0 * timestep_fs)
    trimmed_positions = positions[1:-1]
    return velocities, trimmed_positions


@dataclass
class TrajectoryData:
    """In-memory container for a molecular-dynamics trajectory.

    Attributes
    ----------
    velocities : np.ndarray
        Velocities with shape ``(n_steps, n_atoms, 3)`` in Angstrom / fs.
    symbols : list[str]
        Chemical symbol per atom (length ``n_atoms``).
    masses : np.ndarray
        Atomic masses in amu, shape ``(n_atoms,)``.
    timestep_fs : float
        Simulation timestep in femtoseconds.
    cell : np.ndarray or None
        Unit-cell row vectors ``(3, 3)`` in Angstrom, or ``None`` for
        non-periodic systems.
    positions : np.ndarray or None
        Positions with shape ``(n_steps, n_atoms, 3)`` in Angstrom, or
        ``None`` when only velocities are available.
    """

    velocities: np.ndarray
    symbols: list[str]
    masses: np.ndarray
    timestep_fs: float
    cell: np.ndarray | None = None
    positions: np.ndarray | None = None

    def __post_init__(self) -> None:
        self._validate()

    def _validate(self) -> None:
        velocities = np.asarray(self.velocities, dtype=float)
        if velocities.ndim != 3 or velocities.shape[2] != 3:
            raise ValueError(
                "velocities must have shape (n_steps, n_atoms, 3), got "
                f"{velocities.shape}."
            )
        n_atoms = velocities.shape[1]

        if len(self.symbols) != n_atoms:
            raise ValueError(
                f"len(symbols)={len(self.symbols)} does not match n_atoms={n_atoms}."
            )

        masses = np.asarray(self.masses, dtype=float)
        if masses.shape != (n_atoms,):
            raise ValueError(
                f"masses must have shape ({n_atoms},), got {masses.shape}."
            )

        if self.timestep_fs <= 0:
            raise ValueError(f"timestep_fs must be positive, got {self.timestep_fs}.")

        cell: np.ndarray | None
        if self.cell is None:
            cell = None
        else:
            cell = np.asarray(self.cell, dtype=float)
            if cell.shape != (3, 3):
                raise ValueError(f"cell must have shape (3, 3), got {cell.shape}.")
        self.cell = cell

        positions: np.ndarray | None
        if self.positions is None:
            positions = None
        else:
            positions = np.asarray(self.positions, dtype=float)
            if positions.shape != velocities.shape:
                raise ValueError(
                    "positions must match velocities shape "
                    f"{velocities.shape}, got {positions.shape}."
                )
        self.positions = positions
        self.velocities = velocities
        self.masses = masses

    @property
    def n_steps(self) -> int:
        """Number of stored frames."""
        return self.velocities.shape[0]

    @property
    def n_atoms(self) -> int:
        """Number of atoms."""
        return self.velocities.shape[1]

    @property
    def duration_fs(self) -> float:
        """Total trajectory duration in femtoseconds ``(n_steps - 1) * dt``."""
        return (self.n_steps - 1) * self.timestep_fs

    @property
    def species(self) -> list[str]:
        """Unique chemical species in order of first appearance."""
        seen: list[str] = []
        for sym in self.symbols:
            if sym not in seen:
                seen.append(sym)
        return seen

    def atom_indices(self, symbol: str) -> np.ndarray:
        """Indices of atoms matching ``symbol``."""
        return np.array(
            [idx for idx, sym in enumerate(self.symbols) if sym == symbol],
            dtype=int,
        )


def _frame_has_velocities(atoms: Atoms) -> bool:
    """Whether an ASE ``Atoms`` frame carries velocity/momentum data."""
    return "momenta" in atoms.arrays or "velocities" in atoms.arrays


def load_trajectory(
    source: str | Path,
    timestep_fs: float,
    *,
    stride: int = 1,
    index: slice | str = ":",
    derive_velocities: bool | None = None,
) -> TrajectoryData:
    """Load a trajectory from an ASE-readable file into a :class:`TrajectoryData`.

    Parameters
    ----------
    source : str or Path
        Path to a trajectory file readable by :func:`ase.io.read`
        (``.traj``, ExtXYZ, LAMMPS, VASP, ...).
    timestep_fs : float
        Simulation timestep in femtoseconds.  ASE trajectory files do not
        store the timestep, so it must be supplied by the caller.
    stride : int, default 1
        Take every ``stride``-th frame (``1`` keeps all frames).
    index : slice or str, default ":"
        Frame selector forwarded to :func:`ase.io.read`.
    derive_velocities : bool or None, default None
        Control velocity handling.  ``None`` (auto) uses stored velocities
        when present in every frame and otherwise derives them from
        positions.  ``True`` always derives from positions.  ``False``
        requires stored velocities and raises if they are absent.

    Returns
    -------
    TrajectoryData
        Velocities are converted to Angstrom / fs, masses to amu, and the
        cell is kept in Angstrom.
    """
    path = Path(source)
    if stride < 1:
        raise ValueError(f"stride must be >= 1, got {stride}.")

    frames: Sequence[Atoms] = read(str(path), index=index)
    if isinstance(frames, Atoms):
        frames = [frames]
    frames = list(frames)[::stride]
    if len(frames) == 0:
        raise ValueError(f"No frames read from {path!s}.")

    first = frames[0]
    symbols = list(first.get_chemical_symbols())
    masses = np.asarray(first.get_masses(), dtype=float)
    cell_array = np.asarray(first.get_cell(), dtype=float)
    pbc = np.asarray(first.pbc)
    cell = cell_array if pbc.any() else None

    positions_all = np.asarray([f.get_positions() for f in frames], dtype=float)
    has_vel = all(_frame_has_velocities(f) for f in frames)

    use_stored: bool
    if derive_velocities is True:
        use_stored = False
    elif derive_velocities is False:
        if not has_vel:
            raise ValueError(
                "derive_velocities=False but no velocities are stored in the "
                f"trajectory {path!s}."
            )
        use_stored = True
    else:
        use_stored = has_vel

    if use_stored:
        velocities_ase = np.asarray([f.get_velocities() for f in frames], dtype=float)
        velocities = velocities_ase * _ASE_VEL_TO_ANG_PER_FS
        positions = positions_all
    else:
        if not np.all(positions_all == positions_all):  # pragma: no cover
            pass
        velocities, positions = _derive_velocities_central(
            positions_all, timestep_fs, cell
        )
        warnings.warn(
            "Velocities were not found in the trajectory; deriving them by "
            "central finite differences from positions. The first and last "
            "frames have been trimmed.",
            UserWarning,
            stacklevel=2,
        )

    return TrajectoryData(
        velocities=velocities,
        symbols=symbols,
        masses=masses,
        timestep_fs=timestep_fs,
        cell=cell,
        positions=positions,
    )


def calculate_kinetic_temperature(
    trajectory: TrajectoryData,
    t_md: float | None = None,
    *,
    degrees_of_freedom: int | None = None,
) -> float:
    """Return the average kinetic temperature of ``trajectory`` in kelvin.

    Computes the per-frame kinetic energy

    .. math:: E_k(t) = \\tfrac12 \\sum_i m_i \\lVert\\mathbf{v}_i(t)\\rVert^2

    from velocities in Angstrom / fs and masses in amu, then the kinetic
    temperature :math:`T_{traj} = 2\\langle E_k\\rangle / (3 N k_B)`.

    When ``t_md`` is provided and ``|T_traj - T_md| / T_md > 0.20`` a
    :class:`UserWarning` is emitted (execution proceeds).

    Parameters
    ----------
    trajectory : TrajectoryData
    t_md : float or None
        Nominal simulation temperature in kelvin for validation.
    degrees_of_freedom : int or None
        Number of degrees of freedom used in the denominator.  Defaults to
        ``3 * n_atoms``.
    """
    velocities = trajectory.velocities
    masses = trajectory.masses
    speed_sq = np.sum(velocities**2, axis=2)  # (n_steps, n_atoms)
    kinetic_per_frame = 0.5 * np.sum(masses * speed_sq, axis=1) * _KE_TO_EV
    kinetic_mean = float(np.mean(kinetic_per_frame))

    dof = 3 * trajectory.n_atoms if degrees_of_freedom is None else degrees_of_freedom
    t_traj = 2.0 * kinetic_mean / (dof * units.kB)

    if t_md is not None:
        relative = abs(t_traj - t_md) / abs(t_md)
        if relative > 0.20:
            warnings.warn(
                f"Kinetic temperature {t_traj:.2f} K differs from nominal "
                f"T_MD={t_md:.2f} K by {relative * 100:.1f}% (tolerance 20%).",
                UserWarning,
                stacklevel=2,
            )

    return t_traj
