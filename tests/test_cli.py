"""The command line interface.

Thin by design, so most of these check wiring and failure reporting rather than physics:
that the entry point exists, that arguments reach the library, and that a bad input
produces a non-zero status with a message rather than a traceback. One end-to-end case
does run the invariants, because a smoke test that never checks the answer is not one.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.conftest import write_trajectory
from mdins.cli import main
from mdins.ir import VelocitySpectralDensity
from mdins.units import KB

TEMPERATURE = 300.0
ARGON_MASS = 39.948


def write_thermal_argon(path, n_atoms=64, n_frames=1024, seed=0):
    """Uncorrelated thermal velocities: a flat spectrum with the right total weight.

    Many atoms rather than a handful, so that removing the three centre-of-mass degrees
    of freedom is a small projection and the per-atom equipartition sum rule survives
    it. See :meth:`VelocitySpectralDensity.check_sum_rule` for why that matters.
    """
    rng = np.random.default_rng(seed)
    sigma = np.sqrt(KB * TEMPERATURE / ARGON_MASS)
    velocities = rng.normal(size=(n_frames, n_atoms, 3)) * sigma
    write_trajectory(path, velocities, symbols=f"Ar{n_atoms}")
    return path


@pytest.fixture
def trajectory_file(tmp_path):
    return write_thermal_argon(tmp_path / "traj.extxyz"), 0.005


class TestPdos:
    def test_writes_a_readable_file(self, trajectory_file, tmp_path, capsys):
        path, dt = trajectory_file
        out = tmp_path / "pdos.h5"
        assert main(["pdos", str(path), "--dt", str(dt), "-o", str(out)]) == 0

        density = VelocitySpectralDensity.from_hdf5(out)
        assert density.n_entity == 64
        assert density.metadata.estimator == "welch"
        assert density.metadata.dt == pytest.approx(dt)
        assert "wrote" in capsys.readouterr().err

    def test_end_to_end_result_satisfies_the_invariants(
        self, trajectory_file, tmp_path, capsys
    ):
        path, dt = trajectory_file
        # Per species, not per atom: a single atom's mean square velocity over 1024
        # frames carries a 2.6% standard error, so a few of 64 would cross the 5%
        # tolerance by chance alone. Averaging over the species is the meaningful check.
        main(
            [
                "pdos", str(path), "--dt", str(dt), "-o", str(tmp_path / "p.h5"),
                "--temperature", str(TEMPERATURE), "--by-species",
            ]
        )  # fmt: skip
        err = capsys.readouterr().err
        assert "sum rule satisfied" in err
        assert "warning" not in err

        density = VelocitySpectralDensity.from_hdf5(tmp_path / "p.h5")
        density.check_sum_rule()
        density.check_positive_semidefinite()

    def test_by_species_collapses_atoms(self, trajectory_file, tmp_path):
        path, dt = trajectory_file
        out = tmp_path / "pdos.h5"
        main(["pdos", str(path), "--dt", str(dt), "-o", str(out), "--by-species"])
        density = VelocitySpectralDensity.from_hdf5(out)
        assert density.symbols == ["Ar"]
        assert density.n_entity == 1

    def test_grid_options_reach_the_estimator(self, trajectory_file, tmp_path):
        path, dt = trajectory_file
        out = tmp_path / "pdos.h5"
        main(
            [
                "pdos", str(path), "--dt", str(dt), "-o", str(out),
                "--e-max", "100", "--n-bins", "64",
            ]
        )  # fmt: skip
        density = VelocitySpectralDensity.from_hdf5(out)
        assert density.n_freq == 64
        assert density.frequencies[-1] == pytest.approx(100.0)

    def test_vacf_estimator_is_selectable(self, trajectory_file, tmp_path):
        path, dt = trajectory_file
        out = tmp_path / "pdos.h5"
        main(
            ["pdos", str(path), "--dt", str(dt), "-o", str(out), "--estimator", "vacf"]
        )
        assert VelocitySpectralDensity.from_hdf5(out).metadata.estimator == "vacf"

    def test_reports_the_sampling_limits(self, trajectory_file, tmp_path, capsys):
        path, dt = trajectory_file
        main(["pdos", str(path), "--dt", str(dt), "-o", str(tmp_path / "p.h5")])
        err = capsys.readouterr().err
        assert "Nyquist" in err
        assert "413.6 meV" in err  # h / (2 * 5 fs)

    def test_com_drift_is_removed_before_estimating(self, tmp_path):
        """Without the removal step the estimator refuses the trajectory outright."""
        rng = np.random.default_rng(1)
        velocities = rng.normal(size=(256, 8, 3)) + np.array([20.0, 0.0, 0.0])
        path = tmp_path / "drift.extxyz"
        write_trajectory(path, velocities, symbols="Ar8")
        out = tmp_path / "p.h5"
        assert main(["pdos", str(path), "--dt", "0.005", "-o", str(out)]) == 0

    def test_heating_is_reported(self, tmp_path, capsys):
        rng = np.random.default_rng(2)
        velocities = rng.normal(size=(600, 8, 3))
        velocities *= np.linspace(1.0, 3.0, 600)[:, None, None]
        path = tmp_path / "heating.extxyz"
        write_trajectory(path, velocities, symbols="Ar8")
        main(["pdos", str(path), "--dt", "0.005", "-o", str(tmp_path / "p.h5")])
        assert "drifts by" in capsys.readouterr().err


class TestErrors:
    def test_missing_velocities_reports_rather_than_raises(self, tmp_path, capsys):
        path = tmp_path / "positions.extxyz"
        write_trajectory(path, np.ones((4, 5, 3)), with_velocities=False)
        status = main(["pdos", str(path), "--dt", "0.01", "-o", str(tmp_path / "p.h5")])
        assert status == 1
        assert "no velocities" in capsys.readouterr().err

    def test_aliasing_reports_rather_than_raises(
        self, trajectory_file, tmp_path, capsys
    ):
        path, dt = trajectory_file
        status = main(
            [
                "pdos", str(path), "--dt", str(dt), "-o", str(tmp_path / "p.h5"),
                "--e-max", "5000",
            ]
        )  # fmt: skip
        assert status == 1
        assert "Nyquist" in capsys.readouterr().err

    def test_dt_is_required(self, trajectory_file, tmp_path):
        path, _ = trajectory_file
        with pytest.raises(SystemExit):
            main(["pdos", str(path), "-o", str(tmp_path / "p.h5")])

    def test_version(self, capsys):
        with pytest.raises(SystemExit) as exit_info:
            main(["--version"])
        assert exit_info.value.code == 0
        assert "mdins" in capsys.readouterr().out


class TestEnsembleLabel:
    """``--ensemble`` is a passthrough, and has to reach the file to be worth having."""

    def test_reaches_the_output_file(self, trajectory_file, tmp_path):
        path, dt = trajectory_file
        main(
            [
                "pdos", str(path), "--dt", str(dt), "-o", str(tmp_path / "p.h5"),
                "--ensemble", "NVT",
            ]
        )  # fmt: skip
        density = VelocitySpectralDensity.from_hdf5(tmp_path / "p.h5")
        assert density.metadata.ensemble == "NVT"

    def test_is_absent_unless_asked_for(self, trajectory_file, tmp_path):
        """No format records the ensemble, so the CLI must not invent one."""
        path, dt = trajectory_file
        main(["pdos", str(path), "--dt", str(dt), "-o", str(tmp_path / "p.h5")])
        density = VelocitySpectralDensity.from_hdf5(tmp_path / "p.h5")
        assert density.metadata.ensemble is None
