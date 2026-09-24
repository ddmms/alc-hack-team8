"""Deeper check on atom projection: mass-disordered argon.

Two inequivalent sublattices are created by assigning alternating atom masses
(light = 40 amu, heavy = 160 amu) in a 2x2x2 FCC argon supercell.  Heavy atoms
should contribute more to low-frequency (acoustic) modes and light atoms to
high-frequency modes.  The harmonic atom-projected DOS (Euphonic, via phonon
eigenvectors) and the MD atom-projected pDOS (this pipeline) are computed for
the same system and compared per sublattice.

Run:  uv run python docs/make_mass_projection_plot.py
"""

from __future__ import annotations

import tempfile
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms, units
from ase.build import bulk
from ase.md.velocitydistribution import (
    Stationary,
    ZeroRotation,
    force_temperature,
    thermalize_momenta,
)
from ase.md.verlet import VelocityVerlet
from euphonic import ForceConstants, ureg
from phonopy import Phonopy
from phonopy.file_IO import write_FORCE_CONSTANTS
from phonopy.structure.atoms import PhonopyAtoms

from md_ins.benchmark import make_lj_calculator, make_qgrid
from md_ins.correlation import calculate_pdos

warnings.filterwarnings("ignore")

LATTICE = 5.26
SUPERCELL = (2, 2, 2)
MASS_LIGHT = 40.0
MASS_HEAVY = 160.0
TEMPERATURE_K = 10.0
TIMESTEP_FS = 2.0
N_STEPS = 12000
SEED = 42

THIS_DIR = Path(__file__).resolve().parent


def _unitcell_masses() -> np.ndarray:
    """Mass pattern for the 4-atom cubic FCC cell: [L, H, L, H]."""
    return np.array([MASS_LIGHT, MASS_HEAVY, MASS_LIGHT, MASS_HEAVY])


def _supercell_masses(n_atoms: int) -> np.ndarray:
    return _unitcell_masses()[np.arange(n_atoms) % 4]


def _light_mask(masses: np.ndarray) -> np.ndarray:
    return np.isclose(masses, MASS_LIGHT)


# ---------------------------------------------------------------------------
# Harmonic atom-projected DOS
# ---------------------------------------------------------------------------
unit_atoms = bulk("Ar", "fcc", a=LATTICE, cubic=True)
cell_masses = _unitcell_masses()
unitcell = PhonopyAtoms(
    symbols=unit_atoms.get_chemical_symbols(),
    scaled_positions=unit_atoms.get_scaled_positions(),
    cell=unit_atoms.get_cell().array,
    masses=cell_masses,
)
phonon = Phonopy(unitcell, supercell_matrix=np.diag(SUPERCELL), primitive_matrix="P")
phonon.generate_displacements(distance=0.01)

forces_set = []
for displaced in phonon.supercells_with_displacements:
    calc = make_lj_calculator()
    a = Atoms(
        symbols=displaced.symbols,
        scaled_positions=displaced.scaled_positions,
        cell=displaced.cell,
        pbc=True,
    )
    a.calc = calc
    forces_set.append(a.get_forces())

phonon.forces = forces_set
phonon.produce_force_constants(show_drift=False)
phonon.symmetrize_force_constants()

tmp = Path(tempfile.mkdtemp(prefix="md_ins_mass_"))
write_FORCE_CONSTANTS(phonon.force_constants, filename=str(tmp / "FORCE_CONSTANTS"))
phonon.save(filename=str(tmp / "phonopy.yaml"), settings={"force_constants": True})
fc = ForceConstants.from_phonopy(path=str(tmp))

nq = SUPERCELL[0]
qpoints, weights = make_qgrid(nq)
qpm = fc.calculate_qpoint_phonon_modes(
    qpoints, weights=weights, asr="reciprocal", use_c=True
)

dos_bins = np.arange(0.0, 14.0, 0.04) * ureg("meV")
pdos_collection = qpm.calculate_pdos(dos_bins)
harm_energies = pdos_collection[0].get_bin_centres().to("meV").magnitude

euph_masses = np.array([m.to("amu").magnitude for m in qpm.crystal.atom_mass])
light_mask = _light_mask(euph_masses)
heavy_mask = ~light_mask

harm_light = np.zeros_like(harm_energies)
harm_heavy = np.zeros_like(harm_energies)
for spectrum in pdos_collection:
    idx = int(spectrum.metadata["index"])
    y = spectrum.y_data.to("1/meV").magnitude
    if light_mask[idx]:
        harm_light += y
    else:
        harm_heavy += y
harm_light /= light_mask.sum()
harm_heavy /= heavy_mask.sum()
band_max = float(qpm.frequencies.to("meV").magnitude.max())
print("harmonic band max:", round(band_max, 3), "meV")

# ---------------------------------------------------------------------------
# MD atom-projected pDOS
# ---------------------------------------------------------------------------
atoms = bulk("Ar", "fcc", a=LATTICE, cubic=True).repeat(SUPERCELL)
atoms.set_masses(_supercell_masses(len(atoms)))
atoms.calc = make_lj_calculator()

rng = np.random.default_rng(SEED)
thermalize_momenta(atoms, TEMPERATURE_K, exact_temperature=True, rng=rng)
Stationary(atoms)
ZeroRotation(atoms)
force_temperature(atoms, TEMPERATURE_K)

dynamics = VelocityVerlet(atoms, TIMESTEP_FS * units.fs, logfile=None)
velocities_ase = np.empty((N_STEPS, len(atoms), 3))
for step in range(N_STEPS):
    dynamics.run(1)
    velocities_ase[step] = atoms.get_velocities()
velocities = velocities_ase * units.fs  # Ang/fs

md_masses = np.asarray(atoms.get_masses())
md_light = _light_mask(md_masses)
md_heavy = ~md_light

# use symbols only for species grouping; we project per atom manually here
symbols = ["Ar"] * len(atoms)
result = calculate_pdos(velocities, TIMESTEP_FS, symbols=symbols, window="hann")
md_energies = result.energies_mev
md_light_dos = result.atom_dos[:, md_light].mean(axis=1)
md_heavy_dos = result.atom_dos[:, md_heavy].mean(axis=1)
print("MD light peak:", round(md_energies[np.argmax(md_light_dos)], 3), "meV")
print("MD heavy peak:", round(md_energies[np.argmax(md_heavy_dos)], 3), "meV")
print("harm light peak:", round(harm_energies[np.argmax(harm_light)], 3), "meV")
print("harm heavy peak:", round(harm_energies[np.argmax(harm_heavy)], 3), "meV")


# ---------------------------------------------------------------------------
# Plot: per-sublattice comparison (MD solid, harmonic dashed)
# ---------------------------------------------------------------------------
def _norm(e, y):
    return y / np.trapezoid(y, e)


fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True, sharex=True)

for ax, title, md_d, h_d in [
    (axes[0], f"Light sublattice ({int(MASS_LIGHT)} amu)", md_light_dos, harm_light),
    (axes[1], f"Heavy sublattice ({int(MASS_HEAVY)} amu)", md_heavy_dos, harm_heavy),
]:
    ax.plot(
        md_energies,
        _norm(md_energies, md_d),
        color="tab:blue",
        lw=1.6,
        ls="-",
        label="MD pDOS (this pipeline)",
    )
    ax.plot(
        harm_energies,
        _norm(harm_energies, h_d),
        color="tab:red",
        lw=1.4,
        ls="--",
        label="Harmonic pDOS (Euphonic)",
    )
    ax.set_title(title)
    ax.set_xlabel("Energy (meV)")
    ax.set_xlim(0, max(band_max * 1.1, 10))

axes[0].set_ylabel("Projected DOS (normalised)")
axes[0].legend(loc="upper right")
fig.suptitle(
    "Atom-projected pDOS of mass-disordered FCC argon — MD vs harmonic",
    y=1.02,
)
fig.tight_layout()
fig.savefig(THIS_DIR / "mass_projection_subplots.png", dpi=130, bbox_inches="tight")
plt.close(fig)

# Overlay: both sublattices together, MD solid / harmonic dashed
fig, ax = plt.subplots(figsize=(7.5, 4.8))
ax.plot(
    md_energies,
    _norm(md_energies, md_light_dos),
    color="tab:blue",
    lw=1.6,
    ls="-",
    label="MD light",
)
ax.plot(
    md_energies,
    _norm(md_energies, md_heavy_dos),
    color="tab:orange",
    lw=1.6,
    ls="-",
    label="MD heavy",
)
ax.plot(
    harm_energies,
    _norm(harm_energies, harm_light),
    color="tab:blue",
    lw=1.3,
    ls="--",
    label="Harmonic light",
)
ax.plot(
    harm_energies,
    _norm(harm_energies, harm_heavy),
    color="tab:orange",
    lw=1.3,
    ls="--",
    label="Harmonic heavy",
)
ax.set_xlabel("Energy (meV)")
ax.set_ylabel("Projected DOS (normalised)")
ax.set_title(
    "Mass-disordered argon atom-projected pDOS\n(solid = MD, dashed = harmonic)"
)
ax.set_xlim(0, max(band_max * 1.1, 10))
ax.legend(loc="upper right")
fig.tight_layout()
fig.savefig(THIS_DIR / "mass_projection_overlay.png", dpi=130)
plt.close(fig)
print("plots saved")
