"""Deeper check on atom projection: mass-disordered argon.

Two inequivalent sublattices are created by assigning alternating atom masses
(light = 40 amu, heavy = 160 amu) in a 2x2x2 FCC argon supercell. Heavy atoms
should contribute more to low-frequency (acoustic) modes and light atoms to
high-frequency modes. The harmonic atom-projected DOS (Euphonic, via phonon
eigenvectors) and the MD atom-projected pDOS (this pipeline) are computed for
the same system and compared per sublattice.

Gaussian broadening is applied to the harmonic pDOS via Euphonic's
``Spectrum1DCollection.broaden()`` method to match finite-temperature MD
linewidths and facilitate comparison.

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
from euphonic import ForceConstants, Quantity, Spectrum1D, Spectrum1DCollection, ureg
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
BROADENING_FWHM = 0.35 * ureg("meV")

THIS_DIR = Path(__file__).resolve().parent


def get_unitcell_masses() -> np.ndarray:
    """Return alternating mass pattern for 4-atom cubic FCC cell: [L, H, L, H]."""
    return np.array([MASS_LIGHT, MASS_HEAVY, MASS_LIGHT, MASS_HEAVY])


def get_supercell_masses(n_atoms: int) -> np.ndarray:
    """Tile the alternating mass pattern across all atoms in the supercell."""
    return get_unitcell_masses()[np.arange(n_atoms) % 4]


def is_light_atom(masses: np.ndarray) -> np.ndarray:
    """Return a boolean mask that is True for light sublattice atoms."""
    return np.isclose(masses, MASS_LIGHT)


def build_harmonic_force_constants(
    lattice: float = LATTICE,
    supercell: tuple[int, int, int] = SUPERCELL,
    displacement: float = 0.01,
) -> Phonopy:
    """Build mass-disordered FCC argon and compute force constants with ASE LJ."""
    unit_atoms = bulk("Ar", "fcc", a=lattice, cubic=True)
    cell_masses = get_unitcell_masses()
    unitcell = PhonopyAtoms(
        symbols=unit_atoms.get_chemical_symbols(),
        scaled_positions=unit_atoms.get_scaled_positions(),
        cell=unit_atoms.get_cell().array,
        masses=cell_masses,
    )
    phonon = Phonopy(
        unitcell, supercell_matrix=np.diag(supercell), primitive_matrix="P"
    )
    phonon.generate_displacements(distance=displacement)

    forces_set = []
    for displaced in phonon.supercells_with_displacements:
        calc = make_lj_calculator()
        atoms = Atoms(
            symbols=displaced.symbols,
            scaled_positions=displaced.scaled_positions,
            cell=displaced.cell,
            pbc=True,
        )
        atoms.calc = calc
        forces_set.append(atoms.get_forces())

    phonon.forces = forces_set
    phonon.produce_force_constants(show_drift=False)
    phonon.symmetrize_force_constants()
    return phonon


def compute_harmonic_pdos(
    phonon: Phonopy,
    qgrid_n: int = SUPERCELL[0],
    dos_bins: Quantity | np.ndarray | None = None,
    broadening_fwhm: Quantity | float | None = BROADENING_FWHM,
) -> tuple[Spectrum1DCollection, float]:
    """Compute harmonic atom-projected DOS via Euphonic, optionally broadened.

    Returns a ``Spectrum1DCollection`` containing the light and heavy sublattice
    harmonic pDOS spectra, along with the maximum phonon frequency (meV).
    """
    if dos_bins is None:
        dos_bins = np.arange(0.0, 14.0, 0.04) * ureg("meV")
    elif not isinstance(dos_bins, Quantity):
        dos_bins = dos_bins * ureg("meV")

    with tempfile.TemporaryDirectory(prefix="md_ins_mass_") as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        write_FORCE_CONSTANTS(
            phonon.force_constants, filename=str(tmp_dir / "FORCE_CONSTANTS")
        )
        phonon.save(
            filename=str(tmp_dir / "phonopy.yaml"), settings={"force_constants": True}
        )
        fc = ForceConstants.from_phonopy(path=str(tmp_dir))

    qpoints, weights = make_qgrid(qgrid_n)
    qpm = fc.calculate_qpoint_phonon_modes(
        qpoints, weights=weights, asr="reciprocal", use_c=True
    )
    raw_pdos: Spectrum1DCollection = qpm.calculate_pdos(dos_bins)
    band_max = float(qpm.frequencies.to("meV").magnitude.max())

    euph_masses = np.array([m.to("amu").magnitude for m in qpm.crystal.atom_mass])
    light_mask = is_light_atom(euph_masses)
    heavy_mask = ~light_mask

    # Average per-atom spectra for each sublattice
    light_y = raw_pdos.y_data[light_mask].mean(axis=0)
    heavy_y = raw_pdos.y_data[heavy_mask].mean(axis=0)

    light_spectrum = Spectrum1D(
        x_data=raw_pdos.x_data,
        y_data=light_y,
        metadata={"sublattice": "light", "label": "Harmonic light"},
    )
    heavy_spectrum = Spectrum1D(
        x_data=raw_pdos.x_data,
        y_data=heavy_y,
        metadata={"sublattice": "heavy", "label": "Harmonic heavy"},
    )
    sublattice_collection = Spectrum1DCollection.from_spectra(
        [light_spectrum, heavy_spectrum]
    )

    if broadening_fwhm is not None:
        if not isinstance(broadening_fwhm, Quantity):
            broadening_fwhm = broadening_fwhm * ureg("meV")
        sublattice_collection = sublattice_collection.broaden(
            broadening_fwhm, shape="gauss"
        )

    return sublattice_collection, band_max


def run_md_simulation(
    lattice: float = LATTICE,
    supercell: tuple[int, int, int] = SUPERCELL,
    temperature_k: float = TEMPERATURE_K,
    timestep_fs: float = TIMESTEP_FS,
    n_steps: int = N_STEPS,
    seed: int = SEED,
) -> tuple[np.ndarray, np.ndarray]:
    """Run an NVE MD trajectory for mass-disordered FCC argon.

    Returns ``(velocities, masses)`` where velocities has shape
    ``(n_steps, n_atoms, 3)`` in units of Angstrom/fs.
    """
    atoms = bulk("Ar", "fcc", a=lattice, cubic=True).repeat(supercell)
    masses = get_supercell_masses(len(atoms))
    atoms.set_masses(masses)
    atoms.calc = make_lj_calculator()

    rng = np.random.default_rng(seed)
    thermalize_momenta(atoms, temperature_k, exact_temperature=True, rng=rng)
    Stationary(atoms)
    ZeroRotation(atoms)
    force_temperature(atoms, temperature_k)

    dynamics = VelocityVerlet(atoms, timestep_fs * units.fs, logfile=None)
    velocities_ase = np.empty((n_steps, len(atoms), 3))
    for step in range(n_steps):
        dynamics.run(1)
        velocities_ase[step] = atoms.get_velocities()

    velocities = velocities_ase * units.fs  # convert ASE velocities to Angstrom/fs
    return velocities, masses


def compute_md_pdos(
    velocities: np.ndarray,
    masses: np.ndarray,
    timestep_fs: float = TIMESTEP_FS,
) -> Spectrum1DCollection:
    """Compute atom-projected MD pDOS, aggregated into sublattice spectra.

    Returns a ``Spectrum1DCollection`` containing the light and heavy
    sublattice MD pDOS spectra.
    """
    n_atoms = velocities.shape[1]
    symbols = ["Ar"] * n_atoms
    result = calculate_pdos(velocities, timestep_fs, symbols=symbols, window="hann")

    light_mask = is_light_atom(masses)
    heavy_mask = ~light_mask

    md_energies = result.energies_mev * ureg("meV")
    md_light_y = result.atom_dos[:, light_mask].mean(axis=1) * ureg("1/meV")
    md_heavy_y = result.atom_dos[:, heavy_mask].mean(axis=1) * ureg("1/meV")

    light_spectrum = Spectrum1D(
        x_data=md_energies,
        y_data=md_light_y,
        metadata={"sublattice": "light", "label": "MD light"},
    )
    heavy_spectrum = Spectrum1D(
        x_data=md_energies,
        y_data=md_heavy_y,
        metadata={"sublattice": "heavy", "label": "MD heavy"},
    )
    return Spectrum1DCollection.from_spectra([light_spectrum, heavy_spectrum])


def get_peak_energy(spectrum: Spectrum1D) -> float:
    """Return the energy (in meV) corresponding to the maximum intensity peak."""
    energies = spectrum.get_bin_centres().to("meV").magnitude
    intensities = np.asarray(spectrum.y_data.magnitude)
    return float(energies[np.argmax(intensities)])


def normalize_spectrum(spectrum: Spectrum1D) -> tuple[np.ndarray, np.ndarray]:
    """Return bin centres (meV) and area-normalised intensities for plotting."""
    energies = spectrum.get_bin_centres().to("meV").magnitude
    intensities = np.asarray(spectrum.y_data.magnitude)
    area = float(np.trapezoid(intensities, energies))
    normalized = intensities / area if area > 0 else intensities
    return energies, normalized


def plot_subplots(
    md_spectra: Spectrum1DCollection,
    harm_spectra: Spectrum1DCollection,
    band_max: float,
    output_path: Path,
) -> None:
    """Plot per-sublattice comparison as side-by-side subplots."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True, sharex=True)

    sublattices = [
        ("light", f"Light sublattice ({int(MASS_LIGHT)} amu)", axes[0]),
        ("heavy", f"Heavy sublattice ({int(MASS_HEAVY)} amu)", axes[1]),
    ]

    for key, title, ax in sublattices:
        md_spec = md_spectra.select(sublattice=key)[0]
        harm_spec = harm_spectra.select(sublattice=key)[0]

        md_e, md_y = normalize_spectrum(md_spec)
        harm_e, harm_y = normalize_spectrum(harm_spec)

        ax.plot(
            md_e,
            md_y,
            color="tab:blue",
            lw=1.6,
            ls="-",
            label="MD pDOS (this pipeline)",
        )
        ax.plot(
            harm_e,
            harm_y,
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
    fig.savefig(output_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_overlay(
    md_spectra: Spectrum1DCollection,
    harm_spectra: Spectrum1DCollection,
    band_max: float,
    output_path: Path,
) -> None:
    """Plot all four spectra overlaid on a single axis."""
    md_light = md_spectra.select(sublattice="light")[0]
    md_heavy = md_spectra.select(sublattice="heavy")[0]
    harm_light = harm_spectra.select(sublattice="light")[0]
    harm_heavy = harm_spectra.select(sublattice="heavy")[0]

    md_light_e, md_light_y = normalize_spectrum(md_light)
    md_heavy_e, md_heavy_y = normalize_spectrum(md_heavy)
    harm_light_e, harm_light_y = normalize_spectrum(harm_light)
    harm_heavy_e, harm_heavy_y = normalize_spectrum(harm_heavy)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(
        md_light_e,
        md_light_y,
        color="tab:blue",
        lw=1.6,
        ls="-",
        label="MD light",
    )
    ax.plot(
        md_heavy_e,
        md_heavy_y,
        color="tab:orange",
        lw=1.6,
        ls="-",
        label="MD heavy",
    )
    ax.plot(
        harm_light_e,
        harm_light_y,
        color="tab:blue",
        lw=1.3,
        ls="--",
        label="Harmonic light",
    )
    ax.plot(
        harm_heavy_e,
        harm_heavy_y,
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
    fig.savefig(output_path, dpi=130)
    plt.close(fig)


def main() -> None:
    # 1. Harmonic pDOS via Phonopy + Euphonic (broadened for MD comparison)
    phonon = build_harmonic_force_constants()
    harm_spectra, band_max = compute_harmonic_pdos(
        phonon,
        qgrid_n=SUPERCELL[0],
        broadening_fwhm=BROADENING_FWHM,
    )
    print("harmonic band max:", round(band_max, 3), "meV")

    # 2. MD simulation and atom-projected pDOS
    velocities, masses = run_md_simulation()
    md_spectra = compute_md_pdos(velocities, masses)

    # 3. Peak analysis
    md_light = md_spectra.select(sublattice="light")[0]
    md_heavy = md_spectra.select(sublattice="heavy")[0]
    harm_light = harm_spectra.select(sublattice="light")[0]
    harm_heavy = harm_spectra.select(sublattice="heavy")[0]

    print("MD light peak:", round(get_peak_energy(md_light), 3), "meV")
    print("MD heavy peak:", round(get_peak_energy(md_heavy), 3), "meV")
    print("harm light peak:", round(get_peak_energy(harm_light), 3), "meV")
    print("harm heavy peak:", round(get_peak_energy(harm_heavy), 3), "meV")

    # 4. Save comparison figures
    plot_subplots(
        md_spectra,
        harm_spectra,
        band_max,
        THIS_DIR / "mass_projection_subplots.png",
    )
    plot_overlay(
        md_spectra,
        harm_spectra,
        band_max,
        THIS_DIR / "mass_projection_overlay.png",
    )
    print("plots saved")


if __name__ == "__main__":
    main()
