# Proposal

## Why

Simulating inelastic neutron scattering (INS) spectra for large, disordered, or amorphous molecular systems via ab initio lattice dynamics is computationally prohibitive due to $O(N^3)$ scaling and the breakdown of crystalline symmetry. While molecular dynamics (MD) can model complex systems at scale, existing workflows often either compare raw classical velocity power spectra directly without quantum Debye-Waller or overtone corrections, or require cumbersome specialized pipelines. 

By implementing published methods that transform MD trajectories into intermediate dynamic representations—specifically atom-projected phonon density of states (pDOS) and velocity Cartesian cross-correlation tensors—we can accurately calculate 1-phonon and multiphonon INS spectra incorporating quantum Bose-Einstein populations, Debye-Waller attenuation, and instrument kinematics.

## What Changes

- Implement trajectory analysis routines to extract atom-projected phonon density of states (pDOS) and velocity autocorrelation functions (VACF) from MD trajectories using robust, modern Python libraries.
- Establish a validation pipeline verifying MD-derived pDOS against harmonic lattice dynamics calculated with Euphonic and ASE force calculators (e.g., Lennard-Jones).
- Implement isotropic INS spectrum calculation using atom-projected pDOS in place of explicit displacement tensors (Cheng et al., *J. Chem. Theory Comput.* 2020), including quantum thermal populations, isotropic Debye-Waller factors, instrument momentum-energy relationships ($Q(\omega)$), and iterative scalar overtone/combination convolutions.
- Implement anisotropic/cross-correlation INS spectrum calculation using the outer-product velocity correlation tensor and the "almost isotropic" powder-averaging approximation (Harrelson et al., *Sci. Rep.* 2021), including tensor-contracted overtone convolutions.
- Provide benchmarks against reference calculations and existing INS simulation codes (e.g., OCLIMAX and Euphonic).

## Capabilities

### New Capabilities
- `trajectory-pdos`: Extraction of atom-projected phonon density of states (pDOS) from MD trajectories with harmonic validation benchmarks using Euphonic and ASE.
- `ins-isotropic-sim`: Calculation of INS spectra from pDOS using the fully isotropic approximation, quantum Debye-Waller factors, and multiphonon convolutions.
- `ins-tensor-cross-correlation`: Calculation of INS spectra using velocity cross-correlation tensors, "almost isotropic" powder averaging, and tensor overtone convolutions.

### Modified Capabilities
<!-- None: This is a greenfield implementation with no existing spec-level capabilities modified. -->

## Impact

- Adds new core modules for trajectory velocity processing, vibrational correlation tensors, and INS cross-section calculations.
- Integrates key scientific Python ecosystem dependencies: `numpy`, `scipy`, `ase`, and `euphonic`.
- Adds test suites and benchmark datasets validating pDOS against harmonic lattice dynamics and comparing simulated spectra against published reference results.
