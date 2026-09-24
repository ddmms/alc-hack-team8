# Proposal

## Why

Simulating inelastic neutron scattering (INS) spectra directly from molecular dynamics (MD) trajectories requires an accurate intermediate representation of system dynamics: the atom-projected phonon density of states (pDOS). Generating pDOS from classical trajectories via the velocity autocorrelation function (VACF) bridges atomic simulation engines and vibrational spectroscopy. 

Before implementing downstream INS cross-section calculations, a robust, validated trajectory-to-pDOS pipeline is needed that integrates with modern Python simulation tools (ASE) and is validated against harmonic lattice dynamics (Euphonic).

## What Changes

- Implement trajectory ingestion routines supporting standard simulation formats via `ase.io` (`.traj`, ExtXYZ, LAMMPS, etc.) and in-memory NumPy arrays.
- Implement an in-house vectorised NumPy/SciPy Fourier transform engine to calculate the velocity autocorrelation function (VACF) and atom-projected vibrational density of states (pDOS) $g_d(\omega)$ with windowing (Hann/Blackman).
- Implement kinetic energy and temperature validation: calculate $T_{traj} = \frac{2 \langle E_k \rangle}{3 N k_B}$, compare against nominal simulation temperature $T_{MD}$, and emit a warning if discrepancy exceeds tolerance.
- Construct an automated validation test suite simulating an FCC Lennard-Jones crystal (Argon) in ASE and benchmarking the resulting MD pDOS against the harmonic phonon DOS computed with Euphonic from the identical LJ potential.

This change serves as **Stage 1** of the 3-stage MD-INS pipeline, providing the intermediate pDOS representation used by subsequent INS spectrum calculation changes (`ins-isotropic-sim` and `ins-tensor-anisotropic`).

## Capabilities

### New Capabilities
- `trajectory-pdos`: Extraction of atom-projected phonon density of states (pDOS) from MD trajectories via velocity autocorrelation, with harmonic validation benchmarks using Euphonic and ASE.

### Modified Capabilities
<!-- None: Greenfield capability. -->

## Impact

- Introduces core modules: `trajectory.py` (ingestion/data container), `correlation.py` (VACF & pDOS engine), and `benchmark.py` (Lennard-Jones validation).
- Adds project dependencies on `numpy`, `scipy`, `ase`, and `euphonic`.
- Provides automated test suites validating pDOS normalization, frequency grids, and peak alignment against harmonic Euphonic calculations.
