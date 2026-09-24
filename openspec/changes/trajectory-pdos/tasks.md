# Tasks

## 1. Project Setup & Package Scaffolding

- [ ] 1.1 Create `pyproject.toml` in `alc-hack-team8` defining the `md_ins` package and runtime/dev dependencies (`numpy`, `scipy`, `ase`, `euphonic`, `phonopy`, `pytest`), and verify dependency resolution via `uv sync`
- [ ] 1.2 Create package directory layout `src/md_ins/` and `tests/` with `__init__.py` files, and verify `uv run python -c "import md_ins"` succeeds

## 2. Trajectory Ingestion & Temperature Validation

- [ ] 2.1 Implement `TrajectoryData` dataclass in `src/md_ins/trajectory.py` representing velocities, positions, symbols, masses, timestep, and unit cell, and verify initialization with unit tests on synthetic arrays in `tests/test_trajectory.py`
- [ ] 2.2 Implement ASE trajectory loader function `load_trajectory` supporting `.traj`, ExtXYZ, and other ASE-readable formats with striding and unit conversion to standard units ($\text{\AA}$, $\text{fs}$, $\text{amu}$), verified by loading a sample ASE trajectory
- [ ] 2.3 Implement central-difference velocity derivation fallback when velocities are absent in trajectory coordinates, trimming the boundary frames ($t=0$ and $t=t_{end}$) and raising a visible `UserWarning`, verified by unit tests comparing derived velocities against harmonic oscillator positions
- [ ] 2.4 Implement `calculate_kinetic_temperature` and validation against nominal $T_{MD}$, raising a `UserWarning` if $|T_{traj} - T_{MD}| / T_{MD} > 0.20$, verified by unit tests on hot, cold, and matching trajectory fixtures

## 3. Correlation Engine & Atom-Projected pDOS

- [ ] 3.1 Implement vectorised velocity autocorrelation function (VACF) calculation in `src/md_ins/correlation.py` using Wiener-Khinchin FFT, and verify that $\text{VACF}(0) = 1.0$ and decay is smooth via `tests/test_correlation.py`
- [ ] 3.2 Implement default Hann windowing configured for modes down to tens of $\text{cm}^{-1}$ (with optional Blackman and rectangular filters), and verify window attenuation and low-frequency preservation with unit tests
- [ ] 3.3 Implement Fourier transform to atom-projected and species-projected pDOS $g_d(\omega)$ with internal meV units, frequency unit conversion helpers ($\text{cm}^{-1}$, $\text{THz}$), and area normalization ($\int_0^{\infty} g_d(\omega) d\omega = 1.0$), verified by unit tests on synthetic multi-frequency trajectories
- [ ] 3.4 Implement frequency resolution validation warning when trajectory duration $T_{tot}$ is too short for user-requested frequency resolution ($\Delta\omega < 2\pi / T_{tot}$), verified by unit tests

## 4. Harmonic Benchmark Validation with Phonopy & Euphonic

- [ ] 4.1 Implement in-house ASE-to-Phonopy wrapper in `src/md_ins/benchmark.py` that computes force constants for an ASE Atoms structure and calculator (Lennard-Jones) without external Calorine dependency, and verify force constants matrix generation
- [ ] 4.2 Ingest Phonopy force constants into Euphonic via `ForceConstants.from_phonopy()` and calculate reference harmonic phonon DOS over reciprocal space grid, verified by checking output Euphonic spectrum shape and units
- [ ] 4.3 Implement automated FCC Argon Lennard-Jones simulation runner in `src/md_ins/benchmark.py` generating a low-temperature ($10\text{ K}$) NVE trajectory and computing MD pDOS
- [ ] 4.4 Implement benchmark integration test in `tests/test_benchmark.py` verifying that MD pDOS peak frequencies match Euphonic harmonic phonon DOS within 5% relative error, and verify entire test suite passes with `uv run pytest`
