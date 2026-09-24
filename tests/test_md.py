# SPDX-License-Identifier: GPL-3.0-or-later
"""Layer 3: the whole A-C chain driven by real molecular dynamics.

Everything in `test_spectral.py` feeds the estimator velocity arrays built to have a
known spectrum. That isolates the estimator, but never exercises the parts most likely
to be wrong in practice: ASE's integrator, ASE's unit system, trajectory file I/O, and
how centre-of-mass removal interacts with a real thermostat-free ensemble.

The oracle here is the Einstein crystal — every atom bound independently to its own site
by the same spring, so every mode sits at exactly ``ħ √(k/m)`` and the pDOS is a single
delta. It is the one many-atom system whose spectrum can be written down without solving
anything, which makes it the right first contact with real MD: a disagreement is a bug
in the code, not a question about the physics.

Marked ``slow``: these run MD and take seconds rather than milliseconds.
"""

from __future__ import annotations

import numpy as np
import pytest

from mdins.spectral import velocity_spectral_density
from mdins.trajectory import read_velocities
from mdins.units import HBAR, KB, PLANCK

pytestmark = pytest.mark.slow

#: MD integration timestep in fs. Ten steps per period of the stiffest mode here, which
#: keeps the Verlet energy drift well below the tolerances used below.
TIMESTEP_FS = 1.0

#: Frames are dumped every this many steps. The dump interval, not the timestep, sets
#: the Nyquist energy, and confusing the two is the mistake this whole layer exists to
#: catch.
DUMP_EVERY = 4

N_FRAMES = 4096
TEMPERATURE = 50.0

#: Welch segment length in frames, giving 0.5 meV resolution at the 4 fs dump interval.
#: Not left to the default (an eighth of the trajectory, so 2 meV here): the argon and
#: krypton modes below are only 6.2 meV apart, and at the default they would not be
#: cleanly separated. Resolution is set by this, not by the run length or the bin count.
SEGMENT_LENGTH = 2048


def einstein_energy(spring_constant: float, mass: float) -> float:
    """Exact phonon energy of an Einstein oscillator, in meV.

    ``ω = √(k/m)``, with ``k`` in eV/Å² and ``m`` in amu. The conversion
    ``√(eV / amu) = 98.22 Å/ps`` is the same physical constant as ASE's velocity unit,
    written out here independently so that this oracle does not inherit an error from
    the conversion in :mod:`mdins.trajectory` that it is meant to check.
    """
    angstrom_per_ps = np.sqrt(1.602176634e-19 / 1.66053906660e-27) * 1e-2
    return float(HBAR * angstrom_per_ps * np.sqrt(spring_constant / mass))


def run_einstein_md(path, symbols, spring_constant, *, seed=0, n_frames=N_FRAMES):
    """Integrate an Einstein crystal under NVE and write the trajectory.

    No thermostat: one would add its own coupling to the spectrum at low energy, which
    is precisely the region a spurious feature would hide in. NVE with a fixed seed is
    also reproducible to the bit.
    """
    from ase import Atoms
    from ase import units as ase_units
    from ase.calculators.harmonic import HarmonicCalculator, HarmonicForceField
    from ase.io.trajectory import Trajectory
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.verlet import VelocityVerlet

    rng = np.random.default_rng(seed)
    n_atoms = len(symbols)
    # Sites on a loose cubic grid. The spacing is irrelevant to the physics — the atoms
    # do not interact — but keeping them apart makes the trajectory readable.
    side = int(np.ceil(n_atoms ** (1 / 3)))
    sites = (
        np.array(
            [(i, j, k) for i in range(side) for j in range(side) for k in range(side)],
            dtype=float,
        )[:n_atoms]
        * 3.0
    )

    atoms = Atoms(symbols=list(symbols), positions=sites, pbc=False)
    # An Einstein crystal is a diagonal Hessian: every Cartesian coordinate is bound to
    # its own site by the same spring, with no coupling between atoms or directions.
    hessian = spring_constant * np.eye(3 * n_atoms)
    atoms.calc = HarmonicCalculator(
        HarmonicForceField(ref_atoms=atoms, hessian_x=hessian)
    )

    # Initialise at twice the target: equipartition hands half of it back to the
    # potential as the system settles.
    MaxwellBoltzmannDistribution(atoms, temperature_K=2 * TEMPERATURE, rng=rng)

    dynamics = VelocityVerlet(atoms, timestep=TIMESTEP_FS * ase_units.fs)
    dynamics.run(200)  # discard the initial transient

    with Trajectory(str(path), "w", atoms) as trajectory:
        for _ in range(n_frames):
            dynamics.run(DUMP_EVERY)
            trajectory.write(atoms)

    return DUMP_EVERY * TIMESTEP_FS * 1e-3  # dump interval in ps


@pytest.fixture(scope="module")
def einstein_argon(tmp_path_factory):
    """A 27-atom argon Einstein crystal with its mode at about 20 meV."""
    path = tmp_path_factory.mktemp("md") / "einstein.traj"
    spring_constant = 3.824  # eV/Å², chosen to put the mode near 20 meV
    dt = run_einstein_md(path, ["Ar"] * 27, spring_constant)
    return path, dt, einstein_energy(spring_constant, 39.948)


class TestEinsteinCrystal:
    """The exact case: one spring constant, one mass, one line."""

    def test_the_trajectory_survives_a_round_trip_through_ase(self, einstein_argon):
        """Before trusting any spectrum from it: the file really does carry velocities,
        and they arrive in Å/ps at a physically sane magnitude."""
        path, dt, _ = einstein_argon
        trajectory = read_velocities(path, dt=dt)

        assert trajectory.n_frames == N_FRAMES
        assert trajectory.n_atoms == 27
        assert trajectory.dt == pytest.approx(0.004)

        # Argon near 50 K: ⟨|v|²⟩ = 3kT/m is about 3 (Å/ps)². A factor-of-two band,
        # which is loose as a physics statement but tight enough to catch the failure
        # that actually threatens here — ASE's velocity unit being mistaken for Å/fs or
        # Å/ase-time, which would be out by 10³ or 10².
        mean_square = (trajectory.velocities**2).sum(axis=-1).mean()
        expected = 3 * KB * TEMPERATURE / 39.948
        assert 0.5 * expected < mean_square < 2.0 * expected

    def test_peak_is_at_the_einstein_energy(self, einstein_argon):
        """The headline. ``ħ√(k/m)`` from first principles against the measured peak,
        with no fitted or borrowed constant in between."""
        path, dt, expected = einstein_argon
        trajectory = read_velocities(path, dt=dt).remove_com_velocity()
        density = velocity_spectral_density(
            trajectory, e_max=60.0, n_bins=512, segment_length=SEGMENT_LENGTH
        )

        pdos = density.total_pdos()
        peak = density.frequencies[np.argmax(pdos)]
        resolution = density.metadata.energy_resolution
        assert peak == pytest.approx(expected, abs=2 * resolution), (
            f"Einstein mode expected at {expected:.3f} meV, found at {peak:.3f} meV"
        )

    def test_the_line_is_resolution_limited(self, einstein_argon):
        """An Einstein crystal is exactly harmonic, so the only width the line can have
        is the estimator's own. A broader line would mean the integrator is pumping
        energy between modes that are not supposed to be coupled.
        """
        path, dt, _ = einstein_argon
        trajectory = read_velocities(path, dt=dt).remove_com_velocity()
        density = velocity_spectral_density(
            trajectory, e_max=60.0, n_bins=1024, segment_length=SEGMENT_LENGTH
        )

        pdos = density.total_pdos()
        half = density.frequencies[pdos > 0.5 * pdos.max()]
        fwhm = half.max() - half.min()
        assert fwhm < 4 * density.metadata.energy_resolution

    def test_the_spectrum_is_empty_away_from_the_mode(self, einstein_argon):
        """There is one mode and nothing else. In particular nothing at zero energy,
        which is where drift and a badly handled window both show up, and where the
        ``1/E`` factor in the displacement model will amplify anything left behind.
        """
        path, dt, expected = einstein_argon
        trajectory = read_velocities(path, dt=dt).remove_com_velocity()
        density = velocity_spectral_density(
            trajectory, e_max=60.0, n_bins=512, segment_length=SEGMENT_LENGTH
        )

        pdos = density.total_pdos()
        far = np.abs(density.frequencies - expected) > 5.0
        assert pdos[far].max() < 0.02 * pdos.max()

    def test_the_sum_rule_holds_on_a_real_trajectory(self, einstein_argon):
        """Equipartition, end to end, including the centre-of-mass correction on a cell
        small enough (27 atoms) for it to matter at the 3.7% level.

        Averaged over the species, not per atom, and for a physical reason rather than a
        statistical one. The atoms of an Einstein crystal are uncoupled by construction,
        so under NVE each one conserves its own energy exactly and keeps whatever the
        initial Maxwell-Boltzmann draw gave it forever. The system is not ergodic per
        atom: a time average over one atom converges to that atom's energy, not to
        ``3 k_B T / m``. Equipartition is recovered by the average over atoms, which is
        what ``group_by_species`` forms. A coupled system would not need this.
        """
        path, dt, _ = einstein_argon
        trajectory = read_velocities(path, dt=dt).remove_com_velocity()
        density = velocity_spectral_density(
            trajectory, e_max=60.0, n_bins=512, segment_length=SEGMENT_LENGTH
        )

        assert density.com_projected_mass == pytest.approx(27 * 39.948)
        density.check_positive_semidefinite()
        density.group_by_species().check_sum_rule(rtol=0.05)

    def test_individual_atoms_do_not_satisfy_it(self, einstein_argon):
        """The other half of the statement above, asserted rather than assumed.

        If the per-atom spread were small, the grouping in the previous test would be
        pointless caution and should be removed. It is not: the uncoupled atoms retain
        their initial energies, which are χ²-distributed with three degrees of freedom,
        a spread of order 80%.
        """
        path, dt, _ = einstein_argon
        trajectory = read_velocities(path, dt=dt).remove_com_velocity()
        density = velocity_spectral_density(
            trajectory, e_max=60.0, n_bins=512, segment_length=SEGMENT_LENGTH
        )

        ratio = density.mean_square_velocity() / density.expected_mean_square_velocity()
        assert ratio.std() > 0.3
        assert ratio.mean() == pytest.approx(1.0, rel=0.05)
        with pytest.raises(ValueError, match="sum rule violated"):
            density.check_sum_rule(rtol=0.05)

    def test_the_measured_temperature_is_the_md_temperature(self, einstein_argon):
        """Half the initial kinetic energy goes into the springs, so a run started at
        ``2T`` settles at ``T``. Getting this wrong by the factor of two would rescale
        every intensity downstream."""
        path, dt, _ = einstein_argon
        trajectory = read_velocities(path, dt=dt).remove_com_velocity()
        assert trajectory.temperature() == pytest.approx(TEMPERATURE, rel=0.25)

    def test_the_result_does_not_depend_on_the_estimator(self, einstein_argon):
        """Welch and the VACF route are independent code paths; on an exactly harmonic
        system they have no excuse to disagree about where the line is."""
        path, dt, expected = einstein_argon
        trajectory = read_velocities(path, dt=dt).remove_com_velocity()
        peaks = []
        for estimator in ("welch", "vacf"):
            density = velocity_spectral_density(
                trajectory,
                e_max=60.0,
                n_bins=512,
                estimator=estimator,
                segment_length=SEGMENT_LENGTH,
                max_lag=SEGMENT_LENGTH,
            )
            pdos = density.total_pdos()
            peaks.append(density.frequencies[np.argmax(pdos)])
        assert peaks[0] == pytest.approx(peaks[1], abs=1.0)
        assert peaks[0] == pytest.approx(expected, abs=1.5)


class TestMassDependence:
    """Two species in one cell, bound by the same spring. Their modes must separate as
    ``1/√m``, which no normalisation error can fake and no fitted constant can absorb.
    """

    @pytest.fixture(scope="class")
    def mixed_crystal(self, tmp_path_factory):
        path = tmp_path_factory.mktemp("md") / "mixed.traj"
        spring_constant = 3.824
        symbols = ["Ar"] * 13 + ["Kr"] * 14
        dt = run_einstein_md(path, symbols, spring_constant, seed=1)
        return path, dt, spring_constant

    def test_each_species_peaks_at_its_own_energy(self, mixed_crystal):
        from ase.data import atomic_masses, atomic_numbers

        path, dt, spring_constant = mixed_crystal
        trajectory = read_velocities(path, dt=dt).remove_com_velocity()
        density = velocity_spectral_density(
            trajectory, e_max=60.0, n_bins=1024, segment_length=SEGMENT_LENGTH
        ).group_by_species()

        resolution = density.metadata.energy_resolution
        for index, symbol in enumerate(density.symbols):
            mass = atomic_masses[atomic_numbers[symbol]]
            expected = einstein_energy(spring_constant, mass)
            peak = density.frequencies[np.argmax(density.pdos()[index])]
            assert peak == pytest.approx(expected, abs=2 * resolution), (
                f"{symbol} expected at {expected:.3f} meV, found at {peak:.3f} meV"
            )

    def test_the_two_peaks_are_resolved(self, mixed_crystal):
        """Guards against the test above passing on a single blended peak that happens
        to sit between the two expected energies."""
        from ase.data import atomic_masses, atomic_numbers

        path, dt, spring_constant = mixed_crystal
        energies = [
            einstein_energy(spring_constant, atomic_masses[atomic_numbers[s]])
            for s in ("Ar", "Kr")
        ]
        trajectory = read_velocities(path, dt=dt).remove_com_velocity()
        density = velocity_spectral_density(
            trajectory, e_max=60.0, n_bins=1024, segment_length=SEGMENT_LENGTH
        )
        assert abs(energies[0] - energies[1]) > 4 * density.metadata.energy_resolution

        between = np.abs(density.frequencies - np.mean(energies)) < 0.5
        total = density.total_pdos()
        assert total[between].max() < 0.5 * total.max()


class TestSamplingLimits:
    """The guards in stage A, exercised against a trajectory where the aliasing they
    prevent is real rather than hypothetical."""

    def test_an_undersampled_dump_interval_is_refused(self, einstein_argon):
        """Reading the same trajectory at every thirty-second frame pushes the Nyquist
        energy below the Einstein mode. The mode does not vanish — it folds back to a
        lower energy and looks entirely plausible — so this has to be an error."""
        path, dt, expected = einstein_argon
        coarse = read_velocities(path, dt=32 * dt, index="::32").remove_com_velocity()
        assert coarse.max_resolvable_energy < expected

        with pytest.raises(ValueError, match="Nyquist"):
            velocity_spectral_density(coarse, e_max=60.0)

    def test_aliasing_is_what_the_guard_prevents(self, einstein_argon):
        """Demonstrates the harm, by asking only for energies below the folded-back
        Nyquist limit so the guard stays quiet. The peak moves; nothing about the
        spectrum says it has.
        """
        path, dt, expected = einstein_argon
        coarse = read_velocities(path, dt=32 * dt, index="::32").remove_com_velocity()
        density = velocity_spectral_density(
            coarse, e_max=coarse.max_resolvable_energy, n_bins=256
        )
        peak = density.frequencies[np.argmax(density.total_pdos())]
        assert abs(peak - expected) > 2.0

    def test_the_full_rate_trajectory_is_accepted(self, einstein_argon):
        path, dt, expected = einstein_argon
        trajectory = read_velocities(path, dt=dt)
        assert trajectory.max_resolvable_energy > expected
        assert trajectory.max_resolvable_energy == pytest.approx(PLANCK / (2 * dt))
