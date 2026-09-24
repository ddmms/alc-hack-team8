# Design

## Context

This design covers **Stage 2** of the 3-stage MD-to-INS simulation workflow:
1. `trajectory-pdos` (Stage 1): Trajectory ingestion, VACF, atom-projected pDOS, and Euphonic harmonic pDOS benchmark.
2. **`ins-isotropic-sim` (Stage 2 - Current)**: Isotropic INS cross section (Cheng et al., 2020), TOSCA kinematics, 1D multiphonon convolutions, and `abinslib` benchmark.
3. `ins-tensor-anisotropic` (Stage 3): Outer-product velocity correlation tensors (Harrelson et al., 2021), almost-isotropic powder average, and tensor overtone convolutions.

See `proposal.md` for motivation. In this stage, the intermediate atom-projected pDOS $g_d(\omega)$ generated in Stage 1 is converted into 1-phonon and multiphonon INS spectra under the incoherent isotropic approximation. `abinslib` serves as the reference benchmark code because it implements the isotropic and almost-isotropic powder INS approximations from harmonic input data.

## Goals / Non-Goals

**Goals:**
- **Cross-Section Resolution**: Retrieve `'incoherent'`, `'coherent'`, or `'total'` scattering cross-sections and atomic masses using `euphonic.data.isotopes.IsotopeData` (Sears reference).
- **Spectrometer Kinematics**:
  - Implement TOSCA-like indirect geometry (135° backscattering bank, $E_f \approx 3.32\text{ meV}$) kinematics:
    $$Q^2(\omega) = \frac{m_n}{\hbar^2} \left(E_i + E_f - 2\sqrt{E_i E_f}\cos(135^\circ)\right), \quad E_i = E_f + \hbar\omega$$
  - Support fixed-$Q$ evaluation mode.
- **Isotropic 1-Phonon INS**:
  - Calculate 1-phonon intensity per atom/species:
    $$S_1(Q, \omega) = \sum_d \frac{3 \sigma_d}{2 M_d} \exp(-2 W_d) \frac{Q^2 g_d(\omega)}{\omega} (n(\omega, T_{INS}) + 1)$$
  - Calculate quantum Debye-Waller factor:
    $$W_d = \frac{1}{3} Q^2 \int_0^{\infty} \frac{\hbar g_d(\omega)}{2 M_d \omega} \coth\left(\frac{\hbar\omega}{2 k_B T_{INS}}\right) d\omega$$
- **Multiphonon Iterative Convolutions**:
  - Compute higher-order overtones ($N=2 \dots N_{max}$, default 4) via iterative 1D scalar convolutions:
    $$S_n(\omega) = \int S_1(\omega - \omega') S_{n-1}(\omega') d\omega'$$
  - Weight by appropriate prefactors.
- **Decoupled Broadening**: Provide standalone post-processing functions for fixed-width Gaussian and energy-dependent instrument resolution ($\Delta E / E$), referencing `abinslib` and `resins`.
- **Validation**: Benchmark simulated spectra against `abinslib`'s powder INS calculation using consistent harmonic input data.

**Non-Goals:**
- $3 \times 3$ outer-product Cartesian velocity correlation tensors and almost-isotropic tensor contractions (deferred to Stage 3).
- Single-crystal coherent dispersion curves $S(\mathbf{Q}, \omega)$.

## Decisions

### 1. Cross-Section Lookup via Euphonic IsotopeData
- **Decision**: Query `euphonic.data.isotopes.IsotopeData` to obtain bound coherent cross section ($\sigma_{coh}$), incoherent cross section ($\sigma_{inc}$), total cross section ($\sigma_{tot} = \sigma_{coh} + \sigma_{inc}$), and atomic masses.
- **Rationale**: Reuses a well-tested, standard-compliant database based on Sears (1992).

### 2. Spectrometer Kinematics Architecture
- **Decision**: Create an abstract `Kinematics` base class with implementations:
  - `TOSCAKinematics`: $E_f = 3.32\text{ meV}$, $\theta = 135^\circ$, calculating $Q(\omega)$ along the indirect trajectory.
  - `FixedQKinematics`: Constant $Q$ value across all energy transfers.
- **Rationale**: Directly aligns with the primary experimental geometry of interest (TOSCA at ISIS) while allowing straightforward extension.

### 3. Low-Frequency Regularization & DWF Integration
- **Decision**: The factor $1/\omega$ diverges as $\omega \to 0$. Implement a low-frequency cutoff (default $0.5\text{ meV}$ or user-specified):
  - For $\omega < \omega_{cutoff}$, extrapolate linearly from zero to avoid $1/\omega$ singularity.
  - Compute $W_d$ by numerical Simpson/trapezoidal quadrature of the quantum-weighted pDOS.

### 4. Multiphonon 1D Scalar Convolutions
- **Decision**: Resample $S_1(\omega)$ onto a fine, uniform frequency grid with uniform $\Delta\omega$. Compute successive convolutions using `scipy.signal.fftconvolve`:
  $$S_n(\omega) = \Delta\omega \cdot (S_1 * S_{n-1})(\omega)$$
  Truncate or zero-pad back to the analysis window. Sum $S_{total}(\omega) = \sum_{n=1}^{N_{max}} S_n(\omega)$.

### 5. Decoupled Instrument Broadening
- **Decision**: INS simulation outputs unbroadened spectra. A separate `broadening` module provides:
  - `apply_gaussian_broadening(frequencies, intensity, fwhm)`
  - `apply_instrument_broadening(frequencies, intensity, resolution_func)` where resolution can be a constant fraction $\Delta E / E$ (e.g. $1.5\%$ for TOSCA/VISION).

### 6. Benchmark against Abinslib
- **Decision**: Compare simulated 1-phonon and multiphonon spectra against `abinslib` calculations using harmonic phonon modes from the Lennard-Jones Argon reference.
- **Rationale**: `abinslib` provides the established Python implementation of the powder INS approximation from harmonic input data, ensuring our intermediate pDOS approach matches standard harmonic INS baselines.

## Component Architecture

```
alc-hack-team8/
└── src/
    └── md_ins/
        ├── cross_sections.py   # Euphonic IsotopeData query interface
        ├── kinematics.py       # TOSCA indirect kinematics & fixed-Q
        ├── isotropic.py        # Cheng et al. 1-phonon & scalar overtone engine
        └── broadening.py       # Gaussian and instrument ΔE/E broadening
```

## Risks / Trade-offs

- **[Convolution grid boundary effects]** → Zero-pad the frequency grid by at least $N_{max} \cdot \omega_{max}$ before running FFT convolutions, then slice back to the requested display range.
- **[Discrepancy at low energy due to MD acoustic phonons]** → Trajectory finite size limits acoustic mode sampling near $\Gamma$ ($\omega \to 0$). The low-frequency cutoff and comparison with `abinslib` will document this expected boundary behavior.
