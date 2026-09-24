# Design

## Context

This design covers **Stage 1** of the MD-to-INS simulation workflow. The broader project implements simulation of inelastic neutron scattering (INS) from MD trajectories across three sequential changes:
1. **`trajectory-pdos` (Stage 1 - Current)**: Trajectory ingestion, in-house velocity autocorrelation, atom-projected pDOS generation, and harmonic validation against Euphonic.
2. **`ins-isotropic-sim` (Stage 2)**: Isotropic INS cross-section calculation (Cheng et al., 2020), TOSCA kinematics, multiphonon convolutions, and benchmarking against `abinslib`.
3. **`ins-tensor-anisotropic` (Stage 3)**: Outer-product velocity cross-correlation tensors (Harrelson et al., 2021), almost-isotropic powder average, tensor overtone convolutions, and benchmarking against `abinslib`.

*Note on Benchmarking Tools*: Euphonic provides the reference harmonic phonon density of states (via dynamical matrix / force constants) for validating the pDOS in this stage. For subsequent INS spectrum stages, `abinslib` is the appropriate reference comparison code (rather than Euphonic) because `abinslib` implements the almost-isotropic approximation from harmonic input data.

## Goals / Non-Goals

**Goals:**
- **Trajectory Ingestion**: Provide an `ase.io`-backed loader supporting standard formats (`.traj`, ExtXYZ, LAMMPS, VASP) along with an in-memory `TrajectoryData` container (`velocities`, `symbols`, `masses`, `timestep`, `cell`).
- **In-House Correlation Engine**: Vectorised NumPy/SciPy engine to compute:
  - Velocity autocorrelation function (VACF) per atom and averaged per chemical species.
  - Atom-projected phonon density of states $g_d(\omega)$ via Fourier transform with configurable windowing (Hann, Blackman).
  - Proper normalization: $\int_0^{\infty} g_d(\omega) d\omega = 1$.
- **Kinetic Energy & Temperature Validation**: Compute trajectory temperature $T_{traj} = \frac{2 \langle E_k \rangle}{3 N k_B}$ and log a warning if $|T_{traj} - T_{MD}| / T_{MD} > 20\%$.
- **Harmonic Validation Suite**: Automated benchmark comparing MD-derived pDOS for an FCC Lennard-Jones Argon crystal simulated with ASE at low temperature ($T \approx 10\text{ K}$) against Euphonic's harmonic phonon DOS derived from the identical LJ potential.

**Non-Goals:**
- INS neutron cross sections, Debye-Waller factors, or instrument kinematics (deferred to Stage 2).
- Cartesian velocity cross-correlation tensors (deferred to Stage 3).
- Running production MD simulations (the library provides ingestion and analysis tools, plus a small LJ fixture for testing).

## Decisions

### 1. Trajectory Ingestion via ASE with In-Memory Fallback
- **Decision**: Read trajectories using `ase.io`. Represent the in-memory state with a dataclass:
  ```python
  @dataclass
  class TrajectoryData:
      velocities: np.ndarray  # Shape: (n_steps, n_atoms, 3) in Angstrom/fs or m/s
      symbols: list[str]      # Length: n_atoms
      masses: np.ndarray       # Shape: (n_atoms,) in amu
      timestep_fs: float       # Timestep in femtoseconds
      cell: np.ndarray | None  # Shape: (3, 3) unit cell vectors if periodic
  ```
- **Rationale**: `ase.io` provides transparent support for nearly all atomistic file formats. The decoupled dataclass allows fast unit testing without touching disk.

### 2. In-House Vectorised Correlation Engine
- **Decision**: Implement VACF via Wiener-Khinchin FFT theorem using `scipy.fft`:
  1. Zero-pad velocities along the time axis to length $2 N_{steps}$.
  2. Compute FFT of velocities, take squared magnitude, and inverse FFT to get the autocorrelation.
  3. Apply time-domain window (Hann or Blackman) to damp spectral leakage before Fourier transforming to frequency space.
- **Rationale**: Keeps dependencies minimal, achieves $O(N \log N)$ performance, and directly prepares the codebase for computing cross-correlation components in Stage 3.

### 3. Energy/Frequency Units & Nyquist Grid Management
- **Decision**: Internally manage frequencies in meV (standard for INS) with conversion utilities to $\text{cm}^{-1}$ and THz:
  - Maximum frequency bounded by Nyquist frequency: $f_{max} = \frac{1}{2 \Delta t}$.
  - Frequency resolution bounded by total trajectory duration: $\Delta f = \frac{1}{T_{tot}}$.
  - The module validates whether the trajectory has sufficient duration and time resolution for user-requested frequency ranges.

### 4. Temperature Validation Rationale
- **Decision**: Accept explicit $T_{MD}$ parameter from the user. Compute:
  $$E_k(t) = \frac{1}{2} \sum_{i=1}^N m_i |\mathbf{v}_i(t)|^2, \quad T_{traj} = \frac{2 \langle E_k \rangle}{3 N k_B}$$
  If $|T_{traj} - T_{MD}| / T_{MD} > 0.2$, emit a `UserWarning`. Do not throw an error, as effective temperatures may deliberately differ from kinetic temperatures in certain sampling schemes.

### 5. Harmonic Validation Strategy (Lennard-Jones + Euphonic)
- **Decision**: Implement a self-contained test fixture:
  1. Set up an FCC Argon supercell in ASE ($a \approx 5.26\text{ \AA}$).
  2. Evaluate harmonic force constants using finite differences with ASE's Lennard-Jones calculator ($\epsilon = 0.01042\text{ eV}$, $\sigma = 3.4\text{ \AA}$).
  3. Ingest force constants into Euphonic's `ForceConstants` model and calculate the harmonic phonon DOS over a fine $q$-point grid.
  4. Run a brief ASE Velocity-Verlet NVE trajectory at $T = 10\text{ K}$.
  5. Compute MD pDOS and assert:
     - Integral normalization: $\int g(\omega) d\omega \approx 1$.
     - Acoustic peak location matches Euphonic harmonic pDOS within 5%.

## Component Architecture

```
alc-hack-team8/
├── pyproject.toml
└── src/
    └── md_ins/
        ├── __init__.py
        ├── trajectory.py       # TrajectoryData container & ASE loader
        ├── correlation.py      # VACF & atom-projected pDOS calculation
        └── benchmark.py        # Euphonic + ASE Lennard-Jones benchmark
```

## Risks / Trade-offs

- **[Trajectory memory overhead]** → Support a `stride` parameter during loading to skip frames if trajectory sampling rate is higher than needed for the vibrational frequency range.
- **[Thermal broadening in MD vs delta peaks in harmonic LD]** → At $10\text{ K}$, thermal broadening in MD is minimal; Euphonic's harmonic pDOS can be slightly broadened with a narrow Gaussian ($\sigma \approx 0.2\text{ meV}$) to enable direct Pearson correlation or peak alignment matching.
