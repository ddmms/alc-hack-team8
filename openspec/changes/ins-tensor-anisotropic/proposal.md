# Proposal

## Why

In molecular solids, polymers, and low-symmetry materials, atomic vibrations are strongly anisotropic. The isotropic approximation (Cheng et al., Stage 2) assumes spherically symmetric displacements, which can miss directional orientation effects and alter overtone intensities. The method established by Harrelson et al. (*Sci. Rep.* 2021, 11, 86771) captures vibrational anisotropy directly from MD trajectories by computing outer-product Cartesian velocity cross-correlation tensors, evaluating the Sears/Tomkinson "almost isotropic" powder average for the 1-phonon spectrum, and contracting matrix-multiplied convolutions for higher overtones.

## What Changes

- Implement outer-product Cartesian velocity correlation tensor calculation in the correlation engine:
  $$\mathbf{v}\mathbf{v}^T(\omega) = \mathcal{F}\{\langle \mathbf{v}_i^*(t) \mathbf{v}_i^T(t+\tau) \rangle\}$$
  yielding a $3 \times 3$ Hermitian matrix per atom across all frequency bins.
- Implement mapping via equipartition to frequency-dependent atomic displacement tensors $\overline{\overline{B}}_i(\omega)$ and total displacement tensors $\overline{\overline{A}}_i = \int \overline{\overline{B}}_i(\omega) d\omega$.
- Implement the "almost isotropic" powder-averaged 1-phonon INS spectrum:
  $$S_{0 \to 1}(Q, \omega) \approx \sum_i \frac{Q^2}{3} \text{Tr}(\overline{\overline{B}}_i(\omega)) \exp(-Q^2 \alpha_i(\omega))$$
  where $\alpha_i(\omega) = \frac{1}{5} \left[\text{Tr}(\overline{\overline{A}}_i) + 2 \frac{\overline{\overline{B}}_i(\omega) : \overline{\overline{A}}_i}{\text{Tr}(\overline{\overline{B}}_i(\omega))}\right]$.
- Implement higher overtone calculation via matrix-multiplication frequency convolutions:
  $$\overline{\overline{B}}_n(\omega) = \overline{\overline{B}} * \overline{\overline{B}}_{n-1}$$
  combined with isotropic Debye-Waller factors and prefactors $\frac{3^{n-2}}{n! 5^{n-1}}$.
- Implement benchmark and consistency suites:
  - Verify exact numerical convergence to the Stage 2 isotropic result on isotropic cubic systems (FCC Argon).
  - Compare almost-isotropic spectrum against `abinslib` calculations on anisotropic molecular crystals.

## Capabilities

### New Capabilities
- `ins-tensor-cross-correlation`: Calculation of INS spectra incorporating Cartesian velocity cross-correlation tensors, "almost isotropic" powder averaging, and tensor-contracted overtone convolutions.

### Modified Capabilities
<!-- None: Greenfield capability. -->

## Impact

- Expands `correlation.py` to produce $3 \times 3$ Cartesian tensor representations $\mathbf{v}\mathbf{v}^T(\omega)$.
- Introduces `anisotropic.py` implementing tensor powder averaging and matrix-multiplied convolutions.
- Reuses kinematics, cross-section lookup, and resolution broadening modules developed in Stage 2.
- Verifies consistency against both Stage 2 isotropic calculations and `abinslib`'s almost-isotropic implementation.
