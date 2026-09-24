# Design

## Context

This design covers **Stage 3** of the 3-stage MD-to-INS simulation workflow:
1. `trajectory-pdos` (Stage 1): Trajectory ingestion, VACF, atom-projected pDOS, and Euphonic harmonic pDOS benchmark.
2. `ins-isotropic-sim` (Stage 2): Isotropic INS cross section (Cheng et al., 2020), TOSCA kinematics, 1D multiphonon convolutions, and `abinslib` benchmark.
3. **`ins-tensor-anisotropic` (Stage 3 - Current)**: Outer-product velocity correlation tensors (Harrelson et al., 2021), almost-isotropic powder average, tensor overtone convolutions, and `abinslib` benchmark.

See `proposal.md` for motivation. While Stage 2 treats vibrations as spherically isotropic, Stage 3 captures directional vibrational motions and cross-correlations between Cartesian components directly from MD trajectories.

## Goals / Non-Goals

**Goals:**
- **Outer-Product Velocity Correlation Tensor**:
  - Compute the $3 \times 3$ matrix $\langle \mathbf{v}_i^*(t) \mathbf{v}_i^T(t+\tau) \rangle$ for each atom $i$.
  - Transform via FFT to obtain $\mathbf{v}\mathbf{v}^T(\omega)$ as a Hermitian matrix per frequency bin.
- **Quantum Displacement Tensors**:
  - Map to mode displacement tensor:
    $$\overline{\overline{B}}_i(\omega) = \frac{\hbar}{\omega k_B T_{MD}} \mathbf{v}\mathbf{v}^T(\omega)$$
  - Compute total atomic displacement tensor by frequency integration:
    $$\overline{\overline{A}}_i = \int_0^{\infty} \overline{\overline{B}}_i(\omega) d\omega$$
- **Almost-Isotropic 1-Phonon Powder Average**:
  - Implement Sears/Tomkinson formula for the fundamental excitation:
    $$S_{0 \to 1}(Q, \omega) = \sum_i \frac{Q^2}{3} \text{Tr}(\overline{\overline{B}}_i(\omega)) \exp(-Q^2 \alpha_i(\omega))$$
    $$\alpha_i(\omega) = \frac{1}{5} \left[\text{Tr}(\overline{\overline{A}}_i) + 2 \frac{\sum_{jk} B_{i,jk}(\omega) A_{i,jk}}{\text{Tr}(\overline{\overline{B}}_i(\omega))}\right]$$
- **Multiphonon Matrix Convolutions**:
  - Evaluate higher overtones ($n=2 \dots N_{max}$) by frequency-domain matrix-multiplication convolutions:
    $$(\overline{\overline{B}} * \overline{\overline{B}})_{ik}(\omega) = \sum_j \int B_{ij}(\omega') B_{jk}(\omega - \omega') d\omega'$$
  - Combine with isotropic Debye-Waller factor $\exp(-Q^2 \text{Tr}(\overline{\overline{A}}_i)/3)$ and prefactor $\frac{3^{n-2}}{n! 5^{n-1}}$.
- **Benchmarking & Validation**:
  - Assert numerical equivalence with Stage 2 (`ins-isotropic-sim`) on cubic isotropic systems (LJ Argon) where $\mathbf{v}\mathbf{v}^T(\omega) \propto \mathbf{I}_3$.
  - Compare almost-isotropic intensities against `abinslib` on an anisotropic molecular crystal.

**Non-Goals:**
- Single-crystal angle-resolved coherent scattering (focus remains on powder averaging).
- Fitting anharmonic damping rates $k_j$ from complex MD lineshapes (retains the harmonic ground-state displacement mapping as established in Harrelson et al.).

## Decisions

### 1. Vectorised Outer-Product Correlation Engine
- **Decision**: Extend `correlation.py` to evaluate the 6 independent elements of the symmetric/Hermitian $3 \times 3$ velocity outer product simultaneously using FFT Wiener-Khinchin correlation across all atoms.
- **Rationale**: Computing the full tensor takes approximately $2\text{--}3\times$ the operations of scalar VACF, but vectorisation in NumPy keeps execution fast and avoids Python-level atom loops.

### 2. Double-Dot Tensor Contraction & Regularization
- **Decision**: In computing $\alpha_i(\omega)$, evaluate the tensor contraction term $\overline{\overline{B}}_i(\omega) : \overline{\overline{A}}_i = \sum_{j,k} B_{i,jk}(\omega) A_{i,jk}$. If $\text{Tr}(\overline{\overline{B}}_i(\omega))$ falls below a numerical threshold (e.g., $10^{-12}$), regularize $\alpha_i(\omega) \to \frac{1}{3} \text{Tr}(\overline{\overline{A}}_i)$ (reverting smoothly to the isotropic Debye-Waller factor).
- **Rationale**: Prevents division-by-zero instability at frequency bins where vibrational intensity is negligible.

### 3. Matrix Convolution Algorithm
- **Decision**: Perform the matrix convolution by convolving the component functions using `scipy.signal.fftconvolve`:
  $$C_{ik}(\omega) = \Delta\omega \sum_j (B_{ij} * B_{jk})(\omega)$$
  Ensures the resulting tensor $\overline{\overline{B}}_n(\omega)$ remains Hermitian at every overtone order.

### 4. Convergence Test Suite (Cubic Argon)
- **Decision**: Include a regression test asserting that for an FCC Argon crystal:
  1. The off-diagonal terms of $\overline{\overline{B}}_i(\omega)$ are zero within simulation noise ($< 1\%$ of trace).
  2. $\text{Tr}(\overline{\overline{B}}_i(\omega)) / 3$ matches the scalar $g_i(\omega)$ from Stage 1.
  3. The resulting $S(Q, \omega)$ matches Stage 2's isotropic $S(Q, \omega)$ within 1%.

## Component Architecture

```
alc-hack-team8/
└── src/
    └── md_ins/
        ├── correlation.py      # Expanded with outer_product_velocity_tensor()
        └── anisotropic.py      # Harrelson et al. almost-isotropic powder & matrix overtones
```

## Risks / Trade-offs

- **[Memory footprint of $3 \times 3$ tensors across frequency grid]** → For $N_{atoms}$ and $N_{\omega}$ frequency bins, a full tensor is $(N_{atoms}, N_{\omega}, 3, 3)$ floats. For large systems ($>10^4$ atoms), accumulate species-averaged tensors or stream atom batches rather than allocating all atoms at once.
- **[Computational cost of matrix convolutions]** → Matrix convolution requires $3 \times 3 \times 3 = 27$ 1D FFT convolutions per overtone iteration per species. For a typical species count ($1\text{--}4$ species), this completes in seconds on standard CPUs.
