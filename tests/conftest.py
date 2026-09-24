"""Shared test helpers."""

from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms
from ase import units as ase_units
from ase.io import write


def write_trajectory(
    path, velocities_angstrom_per_ps, symbols="CH4", with_velocities=True
):
    """Write an extxyz trajectory carrying the given velocities in Å/ps."""
    frames = []
    rng = np.random.default_rng(0)
    for velocities in velocities_angstrom_per_ps:
        atoms = Atoms(symbols, positions=rng.random((len(velocities), 3)) * 5.0)
        if with_velocities:
            # ase.units.fs converts Å/fs into ASE's internal velocity unit.
            atoms.set_velocities(velocities / 1e3 / ase_units.fs)
        frames.append(atoms)
    write(str(path), frames, format="extxyz")
    return frames


@pytest.fixture
def extxyz_writer():
    """The ``write_trajectory`` helper, as a fixture."""
    return write_trajectory
