# Spec Delta

## Purpose

Extracts atom-projected phonon density of states (pDOS) and velocity autocorrelation functions from molecular dynamics trajectories, with kinetic temperature validation and harmonic lattice dynamics benchmarks against Euphonic.

## ADDED Requirements

### Requirement: Trajectory Data Ingestion
The system SHALL ingest molecular dynamics trajectory data providing positions, velocities, atomic species, masses, unit cell dimensions, and simulation timestep. The system SHALL support loading from standard atomistic file formats via ASE as well as from in-memory array data. When velocities are present in the input, the system SHALL prefer and use them directly. When velocities are missing from the input, the system SHALL compute velocities by central differences from atomic positions, remove the first and last steps from the trajectory, and emit a visible warning.

#### Scenario: Load trajectory from supported file
- **WHEN** a valid trajectory file (such as `.traj` or ExtXYZ) with velocity information is loaded
- **THEN** the system extracts velocity arrays, atomic symbols, masses, timestep, and unit cell geometry without loss of precision

#### Scenario: Ingest trajectory from in-memory arrays
- **WHEN** raw velocity arrays of shape `(n_steps, n_atoms, 3)` along with masses and timestep are provided directly in memory
- **THEN** the system successfully initializes the trajectory representation without requiring a filesystem read

#### Scenario: Trajectory missing velocity data
- **WHEN** a trajectory containing atomic positions but no velocity information is loaded
- **THEN** the system computes velocities using central finite differences, drops the first and last frames from the trajectory, and emits a visible warning informing the user that velocities were derived from positions

### Requirement: Trajectory Kinetic Temperature Validation
The system SHALL calculate the average instantaneous kinetic temperature of the trajectory from atomic velocities and masses. The system SHALL allow users to specify a nominal simulation temperature and SHALL emit a warning if the calculated kinetic temperature diverges significantly from the nominal value.

#### Scenario: Kinetic temperature matches nominal temperature
- **WHEN** an equilibrated trajectory at nominal temperature $T_{MD}$ is validated and the calculated kinetic temperature is within 20% of $T_{MD}$
- **THEN** the system records the validated temperature without emitting warnings

#### Scenario: Kinetic temperature diverges from nominal temperature
- **WHEN** the trajectory kinetic temperature differs from the user-specified $T_{MD}$ by more than 20%
- **THEN** the system emits a warning detailing the discrepancy while allowing execution to proceed

### Requirement: Velocity Autocorrelation Function Calculation
The system SHALL calculate the normalized velocity autocorrelation function (VACF) for individual atoms and averaged across chemical species over specified lag times using vectorised Fourier transform methods.

#### Scenario: VACF calculation on equilibrium trajectory
- **WHEN** an equilibrium velocity trajectory is provided
- **THEN** the calculated VACF at zero lag equals 1.0 and decays smoothly with increasing correlation time

### Requirement: Atom-Projected Phonon Density of States
The system SHALL compute the atom-projected and species-projected vibrational density of states $g_d(\omega)$ by Fourier transformation of the velocity autocorrelation function. The system SHALL provide a robust default Hann window that accurately resolves low-frequency vibrational modes down to tens of $\text{cm}^{-1}$ (approximately 2 to 5 meV) without requiring user window configuration, while optionally supporting alternative windowing filters (Blackman, rectangular). The system SHALL normalize the pDOS such that the integral over positive frequencies equals 1.0 for each atomic degree of freedom.

#### Scenario: Default Hann windowing on low-frequency modes
- **WHEN** a trajectory containing low-frequency modes down to tens of $\text{cm}^{-1}$ is processed with default parameters
- **THEN** the system applies the default Hann window, suppresses spectral leakage, and resolves the low-frequency features without user window configuration

#### Scenario: pDOS calculation and normalization
- **WHEN** a velocity trajectory is processed with a specified energy grid and windowing filter
- **THEN** the output pDOS is evaluated on the frequency grid and its numerical integral over positive frequencies equals 1.0 within numerical precision (±1%)

#### Scenario: Trajectory duration insufficient for requested resolution
- **WHEN** the total trajectory duration provides a frequency resolution coarser than requested by the user
- **THEN** the system issues a warning informing the user that the requested energy resolution exceeds the physical resolution limit of the trajectory

### Requirement: Harmonic Benchmark Validation against Euphonic
The system SHALL provide an automated benchmark suite that computes the pDOS of an FCC Lennard-Jones crystal from an ASE molecular dynamics trajectory and verifies that the vibrational frequency peaks match the harmonic phonon density of states calculated by Euphonic from the identical Lennard-Jones potential using force constants generated via a lightweight Phonopy wrapper.

#### Scenario: Closed-loop Lennard-Jones benchmark execution
- **WHEN** the benchmark simulation of an FCC Argon crystal is executed at low temperature (approximately 10 K)
- **THEN** the peak frequencies of the MD-derived pDOS agree with Euphonic's harmonic phonon DOS within 5% relative error
