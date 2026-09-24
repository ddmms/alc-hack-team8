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
      velocities: np.ndarray  # Shape: (n_steps, n_atoms, 3) in Angstrom/fs
      symbols: list[str]  # Length: n_atoms
      masses: np.ndarray  # Shape: (n_atoms,) in amu
      timestep_fs: float  # Timestep in femtoseconds
      cell: np.ndarray | None  # Shape: (3, 3) unit cell vectors if periodic
  ```
  If precomputed velocities exist in the trajectory, use them directly. If velocities are missing, calculate them via central finite differences from positions $\mathbf{v}(t) = \frac{\mathbf{r}(t+\Delta t) - \mathbf{r}(t-\Delta t)}{2 \Delta t}$ (applying minimum image convention for periodic cells), drop the first and last steps from the trajectory, and emit a visible `UserWarning`.
- **Rationale**: `ase.io` provides transparent support for nearly all atomistic file formats. Providing a velocity derivation fallback ensures that position-only trajectory files (common in lightweight simulations) remain usable, while visibly notifying the user of the numerical differentiation and frame trimming.

### 2. In-House Vectorised Correlation Engine
- **Decision**: Implement VACF via Wiener-Khinchin FFT theorem using `scipy.fft`:
  1. Zero-pad velocities along the time axis to length $2 N_{steps}$.
  2. Compute FFT of velocities, take squared magnitude, and inverse FFT to get the autocorrelation.
  3. Apply a robust default Hann window to damp spectral leakage before Fourier transforming to frequency space. The default window length and shape are configured to cleanly resolve low-frequency vibrational and librational modes down to tens of $\text{cm}^{-1}$ (approximately 2–5 meV) without requiring user tuning, while preserving optional customization (Blackman, rectangular).
- **Rationale**: Keeps dependencies minimal, achieves $O(N \log N)$ performance, protects end-users from windowing subtleties while ensuring high spectral fidelity in the low-frequency acoustic regime, and directly prepares the codebase for computing cross-correlation components in Stage 3.

### 3. Energy/Frequency Units & Nyquist Grid Management
- **Decision**: Internally manage frequencies in meV (standard for INS) with conversion utilities to $\text{cm}^{-1}$ and THz:
  - Maximum frequency bounded by Nyquist frequency: $f_{max} = \frac{1}{2 \Delta t}$.
  - Frequency resolution bounded by total trajectory duration: $\Delta f = \frac{1}{T_{tot}}$.
  - The module validates whether the trajectory has sufficient duration and time resolution for user-requested frequency ranges.

### 4. Temperature Validation Rationale
- **Decision**: Accept explicit $T_{MD}$ parameter from the user. Compute:
  $$E_k(t) = \frac{1}{2} \sum_{i=1}^N m_i |\mathbf{v}_i(t)|^2, \quad T_{traj} = \frac{2 \langle E_k \rangle}{3 N k_B}$$
  If $|T_{traj} - T_{MD}| / T_{MD} > 0.2$, emit a `UserWarning`. Do not throw an error, as effective temperatures may deliberately differ from kinetic temperatures in certain sampling schemes.

### 5. Harmonic Validation Strategy (Lennard-Jones + Phonopy + Euphonic)
- **Decision**: Implement a self-contained test benchmark:
  1. Set up an FCC Argon supercell in ASE ($a \approx 5.26\text{ \AA}$).
  2. Because Euphonic supports CASTEP and Phonopy formats for force constants but not ASE's native `ase.phonons`, implement a lightweight helper function connecting ASE calculators to `phonopy` (following the pattern in Calorine, without introducing Calorine as a dependency):
     - Convert the ASE structure to `PhonopyAtoms`.
     - Generate finite displacements using `Phonopy(structure, supercell_matrix)`.
     - Evaluate supercell forces using ASE's Lennard-Jones calculator ($\epsilon = 0.01042\text{ eV}$, $\sigma = 3.4\text{ \AA}$).
     - Compute force constants in Phonopy and export to Euphonic via `ForceConstants.from_phonopy()`.
  3. Calculate the reference harmonic phonon DOS over a reciprocal-space $q$-point grid using Euphonic's `calculate_qpoint_phonon_modes().calculate_dos()`.
  4. Run a brief ASE Velocity-Verlet NVE trajectory of the same system at $T = 10\text{ K}$.
  5. Compute MD pDOS with default Hann windowing and assert:
     - Integral normalization: $\int g(\omega) d\omega \approx 1$.
     - Acoustic peak location matches Euphonic harmonic pDOS within 5%.
- **Rationale**: Keeps the test suite fully autonomous and exact while cleanly bridging ASE calculators, Phonopy force constants, and Euphonic spectral calculation without adding unnecessary heavy dependencies.

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
