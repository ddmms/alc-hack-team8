"""Unit conversions and physical constants.

Cheap tests, but they guard the factor-of-2π and factor-of-mass errors that are the
characteristic failure mode of this kind of calculation: a spectrum that looks entirely
plausible and is wrong by a constant.
"""

from __future__ import annotations

import numpy as np
import pytest

from mdins import units


class TestConstants:
    """Constants against independently known values (CODATA/NIST)."""

    def test_hbar_in_mev_ps(self):
        assert units.HBAR == pytest.approx(0.6582119569, rel=1e-9)

    def test_planck_in_mev_ps(self):
        assert units.PLANCK == pytest.approx(4.135667696, rel=1e-9)

    def test_boltzmann_in_mev_per_kelvin(self):
        assert units.KB_MEV == pytest.approx(0.08617333262, rel=1e-9)

    def test_boltzmann_in_md_units(self):
        # Equals the gas constant in kJ/mol/K, 0.0083144626, scaled to amu Å²/ps²/K.
        assert units.KB == pytest.approx(0.831446262, rel=1e-8)

    def test_energy_of_one_amu_angstrom2_per_ps2(self):
        assert units.MEV_PER_AMU_ANGSTROM2_PS2 == pytest.approx(0.103642696, rel=1e-8)

    def test_planck_is_two_pi_hbar(self):
        assert units.PLANCK == pytest.approx(2.0 * np.pi * units.HBAR, rel=1e-12)


class TestConversions:
    """Known conversion factors, and round-trips."""

    def test_one_terahertz_in_mev(self):
        assert units.frequency_to_energy(1.0) == pytest.approx(4.135667696, rel=1e-9)

    def test_one_mev_in_wavenumbers(self):
        assert units.energy_to_wavenumber(1.0) == pytest.approx(8.065543937, rel=1e-8)

    def test_angular_frequency_carries_the_two_pi(self):
        energy = 10.0
        omega = units.energy_to_angular_frequency(energy)
        nu = units.energy_to_frequency(energy)
        assert omega == pytest.approx(2.0 * np.pi * nu, rel=1e-12)

    @pytest.mark.parametrize(
        ("forward", "backward"),
        [
            (units.energy_to_frequency, units.frequency_to_energy),
            (units.energy_to_wavenumber, units.wavenumber_to_energy),
            (units.energy_to_angular_frequency, units.angular_frequency_to_energy),
        ],
    )
    def test_round_trip(self, forward, backward):
        energies = np.array([0.0, 1.0, 12.5, 500.0])
        np.testing.assert_allclose(backward(forward(energies)), energies, rtol=1e-12)


class TestEquipartition:
    """The invariant that anchors the absolute normalisation (design.md D6)."""

    def test_hydrogen_at_room_temperature(self):
        # <|v|²> = 3kT/m; for 1 amu at 300 K this is ~748 (Å/ps)², i.e. ~27 Å/ps rms.
        got = units.thermal_velocity_squared(1.008, 300.0)
        assert got == pytest.approx(3 * units.KB * 300.0 / 1.008, rel=1e-12)
        assert np.sqrt(got) == pytest.approx(27.2, rel=1e-2)

    def test_scales_inversely_with_mass(self):
        light, heavy = units.thermal_velocity_squared([1.0, 100.0], 300.0)
        assert light / heavy == pytest.approx(100.0, rel=1e-12)

    def test_scales_linearly_with_temperature(self):
        cold = units.thermal_velocity_squared([12.0], 10.0)
        hot = units.thermal_velocity_squared([12.0], 300.0)
        assert (hot / cold)[0] == pytest.approx(30.0, rel=1e-12)

    def test_matches_ase_maxwell_boltzmann(self):
        """Cross-check against ASE, which uses its own unit system throughout.

        ASE holds velocities in Å per ASE time unit. ``ase.units.fs`` is the size of a
        femtosecond in those units, so multiplying converts to Å/fs.
        """
        from ase import Atoms
        from ase import units as ase_units
        from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

        rng = np.random.default_rng(0)
        atoms = Atoms("Ar" * 2000, positions=rng.random((2000, 3)))
        MaxwellBoltzmannDistribution(atoms, temperature_K=300.0, rng=rng)

        velocities = atoms.get_velocities() * ase_units.fs * 1e3  # Å/ps
        got = (velocities**2).sum(axis=1).mean()
        want = units.thermal_velocity_squared([atoms[0].mass], 300.0)[0]
        assert got == pytest.approx(want, rel=0.05)
