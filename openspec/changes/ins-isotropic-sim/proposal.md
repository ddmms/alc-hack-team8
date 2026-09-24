# Proposal

## Why

Following the generation of atom-projected phonon density of states (pDOS) from trajectories in Stage 1 (`trajectory-pdos`), the next capability is simulating inelastic neutron scattering (INS) spectra. The method established by Cheng et al. (*J. Chem. Theory Comput.* 2020, 16, 7702) uses the atom-projected pDOS in place of full lattice dynamics displacement tensors in a fully-isotropic, incoherent powder approximation. This decouples classical MD sampling from quantum thermal populations, enabling accurate 1-phonon and multiphonon spectral simulation for complex materials.

## What Changes

- Implement neutron scattering cross-section resolution supporting `'incoherent'`, `'coherent'`, and `'total'` options using Euphonic's `IsotopeData` (Sears reference).
- Implement spectrometer kinematics starting with a prototype for TOSCA-like indirect geometry (135° backscattering bank, $E_f \approx 3.32\text{ meV}$) as well as fixed-$Q$ modes, following `abinslib` conventions.
- Implement isotropic 1-phonon INS spectrum calculation from atom-projected pDOS $g_d(\omega)$:
  $$S_1(Q, \omega) = \sum_d \frac{3 \sigma_d}{2 M_d} \exp(-2 W_d) \frac{Q^2 g_d(\omega)}{\omega} (n(\omega, T_{INS}) + 1)$$
  with quantum Debye-Waller factor $W_d = \frac{1}{3} Q^2 \int \frac{\hbar g_d(\omega)}{2 M_d \omega} (2 n(\omega, T_{INS}) + 1) d\omega$.
- Implement multiphonon combination and overtone evaluation via iterative 1D scalar convolutions: $S_n(\omega) = S_1 * S_{n-1}$ with overtone prefactors.
- Implement decoupled spectral resolution broadening (configurable fixed-width Gaussian and energy-dependent instrument $\Delta E / E$).
- Implement a validation and benchmark suite comparing simulated spectra against `abinslib`'s implementation from harmonic input data.

## Capabilities

### New Capabilities
- `ins-isotropic-sim`: Calculation of 1-phonon and multiphonon inelastic neutron scattering (INS) spectra from atom-projected pDOS using the fully isotropic approximation, quantum Debye-Waller factors, and scalar convolutions, benchmarked against `abinslib`.

### Modified Capabilities
<!-- None: Greenfield capability. -->

## Impact

- Builds directly upon the `trajectory-pdos` intermediate representation.
- Introduces modules: `cross_sections.py`, `kinematics.py`, `isotropic.py`, and `broadening.py`.
- Integrates `abinslib` into the benchmark suite for reference harmonic INS calculations.
