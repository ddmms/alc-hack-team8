# SPDX-License-Identifier: GPL-3.0-or-later
"""Stage A: trajectory input.

Reads **velocities** through ASE. Positions are not sufficient, which is the tightest
constraint in the pipeline and is enforced rather than discovered: a format that does
not carry velocities produces an error naming the problem, not a spectrum of zeros.

ASE holds velocities in Å per ASE time unit; everything here is Å/ps.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from mdins.provenance import Provenance
from mdins.units import KB, PLANCK

if TYPE_CHECKING:  # pragma: no cover
    from os import PathLike

__all__ = ["VelocityTrajectory", "read_velocities"]

_FORMAT_ADVICE = (
    "Formats that carry velocities and that ASE can read include ASE .traj, extxyz "
    "with a velocities/momenta array, and LAMMPS dumps written with vx vy vz. "
    "GROMACS .xtc and DCD do not store velocities at all; GROMACS .trr stores them but "
    "ASE cannot read it. Convert to extxyz or .traj first."
)


@dataclass(frozen=True, eq=False)
class VelocityTrajectory:
    """Atomic velocities on a uniform time grid.

    Attributes:
        velocities: Shape ``(n_frames, n_atoms, 3)``, in Å/ps.
        masses: Shape ``(n_atoms,)``, in amu.
        symbols: Chemical symbol of each atom.
        dt: Interval between the frames actually held, in ps. This is the trajectory's
            dump interval, not the MD integration timestep, and confusing the two is the
            easiest way to get every energy wrong by a constant factor.
        provenance: What produced this trajectory.
        positions: Optional ``(n_frames, n_atoms, 3)`` in Å, only needed for angular
            velocity removal.
        dof_removed: Degrees of freedom already projected out, used when inferring
            temperature.
    """

    velocities: NDArray[np.float64]
    masses: NDArray[np.float64]
    symbols: list[str]
    dt: float
    provenance: Provenance
    positions: NDArray[np.float64] | None = None
    dof_removed: int = 0

    def __post_init__(self) -> None:
        velocities = np.ascontiguousarray(self.velocities, dtype=np.float64)
        masses = np.ascontiguousarray(self.masses, dtype=np.float64)

        if velocities.ndim != 3 or velocities.shape[2] != 3:
            raise ValueError(
                f"velocities must have shape (n_frames, n_atoms, 3), "
                f"got {velocities.shape}"
            )
        if velocities.shape[1] != masses.size:
            raise ValueError(
                f"{velocities.shape[1]} atoms in velocities but {masses.size} masses"
            )
        if len(self.symbols) != masses.size:
            raise ValueError(f"{len(self.symbols)} symbols for {masses.size} masses")
        if self.dt <= 0.0:
            raise ValueError(f"dt must be positive, got {self.dt}")
        if velocities.shape[0] < 2:
            raise ValueError("need at least two frames to form a spectrum")

        object.__setattr__(self, "velocities", velocities)
        object.__setattr__(self, "masses", masses)
        object.__setattr__(self, "symbols", list(self.symbols))

    @property
    def n_frames(self) -> int:
        """Number of frames."""
        return self.velocities.shape[0]

    @property
    def n_atoms(self) -> int:
        """Number of atoms."""
        return self.velocities.shape[1]

    @property
    def duration(self) -> float:
        """Total sampled time, in ps."""
        return self.n_frames * self.dt

    @property
    def total_mass(self) -> float:
        """Mass of the whole system, in amu."""
        return float(self.masses.sum())

    @property
    def com_removed(self) -> bool:
        """Whether the centre-of-mass velocity has been projected out.

        Read from provenance rather than from a flag, so it stays true of trajectories
        that arrived already cleaned and said so.
        """
        return any("remove_com_velocity" in step for step in self.provenance.steps)

    # -- sampling limits ---------------------------------------------------------

    @property
    def max_resolvable_energy(self) -> float:
        """Nyquist energy in meV.

        Motion above this aliases back into the spectrum rather than being filtered out.
        """
        return PLANCK / (2.0 * self.dt)

    @property
    def energy_resolution(self) -> float:
        """Finest resolvable energy spacing in meV, set by the total duration."""
        return PLANCK / self.duration

    def validate_sampling(
        self, e_max: float, e_resolution: float | None = None
    ) -> None:
        """Check the trajectory can support the requested energy range.

        Raises rather than warns. Aliasing is silent and unrecoverable after the fact:
        an aliased spectrum looks like a real one, so a warning would be read as
        advisory when it is fatal.

        Args:
            e_max: Highest energy transfer required, in meV.
            e_resolution: Required energy resolution in meV, if any.

        Raises:
            ValueError: If the dump interval or the duration is insufficient.
        """
        if e_max > self.max_resolvable_energy:
            needed = PLANCK / (2.0 * e_max)
            raise ValueError(
                f"requested energies up to {e_max:.4g} meV but the Nyquist limit for a "
                f"{self.dt * 1e3:.4g} fs dump interval is "
                f"{self.max_resolvable_energy:.4g} meV. Motion above the limit aliases "
                f"back into the spectrum and cannot be removed afterwards. Dump at "
                f"{needed * 1e3:.4g} fs or shorter, or lower e_max."
            )
        if e_resolution is not None and e_resolution < self.energy_resolution:
            needed = PLANCK / e_resolution
            raise ValueError(
                f"requested {e_resolution:.4g} meV resolution but "
                f"{self.duration:.4g} ps of trajectory supports only "
                f"{self.energy_resolution:.4g} meV. "
                f"Run for at least {needed:.4g} ps."
            )

    # -- pre-processing ----------------------------------------------------------

    def remove_com_velocity(self) -> VelocityTrajectory:
        """Subtract the centre-of-mass velocity from every frame.

        Drift puts a spurious feature at zero energy, which the ``1/E`` factor in the
        displacement model then amplifies.
        """
        total_mass = self.masses.sum()
        com = (
            np.einsum("a,fad->fd", self.masses, self.velocities) / total_mass
        )  # (n_frames, 3)
        return replace(
            self,
            velocities=self.velocities - com[:, None, :],
            dof_removed=self.dof_removed + 3,
            provenance=self.provenance.with_step("remove_com_velocity"),
        )

    def remove_angular_velocity(self) -> VelocityTrajectory:
        """Subtract rigid-body rotation from every frame.

        Only meaningful for a non-periodic system: for a periodic cell there is no
        well-defined global rotation, and applying this would be wrong rather than
        merely unnecessary. Requires positions.

        Raises:
            ValueError: If the trajectory was read without positions.
        """
        if self.positions is None:
            raise ValueError(
                "angular velocity removal needs positions; "
                "read the trajectory with keep_positions=True"
            )

        masses = self.masses
        total_mass = masses.sum()
        centre = np.einsum("a,fad->fd", masses, self.positions) / total_mass
        rel = self.positions - centre[:, None, :]

        angular_momentum = np.einsum(
            "a,fad->fd", masses, np.cross(rel, self.velocities)
        )
        # Inertia tensor per frame: I = Σ m (r·r δ - r ⊗ r)
        r2 = np.einsum("fad,fad->fa", rel, rel)
        inertia = np.einsum("a,fa,de->fde", masses, r2, np.eye(3)) - np.einsum(
            "a,fad,fae->fde", masses, rel, rel
        )
        # Trailing axis kept explicit: numpy reads a 2-D right-hand side as a single
        # matrix rather than a batch of vectors.
        omega = np.linalg.solve(inertia, angular_momentum[..., None])[..., 0]
        return replace(
            self,
            velocities=self.velocities - np.cross(omega[:, None, :], rel),
            dof_removed=self.dof_removed + 3,
            provenance=self.provenance.with_step("remove_angular_velocity"),
        )

    # -- diagnostics -------------------------------------------------------------

    def temperature(self) -> float:
        """Kinetic temperature in K, averaged over frames.

        Uses ``3 N - dof_removed`` degrees of freedom, so it stays correct after
        centre-of-mass removal.
        """
        kinetic = np.einsum(
            "a,fad,fad->", self.masses, self.velocities, self.velocities
        )
        kinetic /= self.n_frames
        dof = 3 * self.n_atoms - self.dof_removed
        if dof <= 0:
            raise ValueError(
                f"{self.n_atoms} atoms with {self.dof_removed} degrees of freedom "
                "projected out leaves nothing to measure a temperature from. A single "
                "atom with its own centre-of-mass motion removed has no velocity left, "
                "and no spectrum can be computed from it."
            )
        return float(kinetic / (dof * KB))

    def temperature_drift(self) -> float:
        """Fractional difference in kinetic temperature between the first and last
        thirds of the trajectory.

        A large value means the run was not equilibrated, and the spectral density is an
        average over states the system was passing through rather than a property of
        one.
        """
        third = max(self.n_frames // 3, 1)
        squared = np.einsum(
            "a,fad,fad->f", self.masses, self.velocities, self.velocities
        )
        early = squared[:third].mean()
        late = squared[-third:].mean()
        return float(abs(late - early) / early)


def read_velocities(
    path: str | PathLike[str],
    *,
    dt: float,
    index: str = ":",
    format: str | None = None,  # noqa: A002 - matches ase.io.read
    keep_positions: bool = False,
) -> VelocityTrajectory:
    """Read a trajectory's velocities through ASE.

    Args:
        path: Trajectory file.
        dt: Interval between the frames being read, in ps. Required, and not inferred:
            most trajectory formats do not record it, and a wrong value rescales every
            energy in the final spectrum without any other symptom. If ``index`` skips
            frames, this is the interval after skipping.
        index: ASE slice string selecting frames.
        format: ASE format name; inferred from the extension when omitted.
        keep_positions: Also retain positions, needed only for
            :meth:`VelocityTrajectory.remove_angular_velocity`.

    Returns:
        The trajectory, with velocities in Å/ps.

    Raises:
        ValueError: If the file carries no velocities.
    """
    from ase.io import read
    from ase.units import fs as ase_fs

    frames = read(str(path), index=index, format=format)
    if not isinstance(frames, list):
        frames = [frames]
    if not frames:
        raise ValueError(f"no frames read from {path}")

    first = frames[0]
    if not first.has("momenta"):
        raise ValueError(
            f"{path} contains no velocities, so no spectrum can be computed from it. "
            + _FORMAT_ADVICE
        )

    # ASE velocities are Å per ASE time unit; ase.units.fs is a femtosecond in those
    # units, so multiplying gives Å/fs.
    to_angstrom_per_ps = ase_fs * 1e3
    velocities = (
        np.array([frame.get_velocities() for frame in frames], dtype=np.float64)
        * to_angstrom_per_ps
    )

    if not np.any(velocities):
        raise ValueError(
            f"all velocities in {path} are zero. This usually means the file has a "
            "momenta array that was never populated. " + _FORMAT_ADVICE
        )

    positions = (
        np.array([frame.get_positions() for frame in frames], dtype=np.float64)
        if keep_positions
        else None
    )

    provenance = Provenance(source=str(path)).with_step(
        f"read_velocities(index={index!r}, dt={dt} ps, n_frames={len(frames)})"
    )
    return VelocityTrajectory(
        velocities=velocities,
        masses=np.asarray(first.get_masses(), dtype=np.float64),
        symbols=list(first.get_chemical_symbols()),
        dt=dt,
        provenance=provenance,
        positions=positions,
    )
