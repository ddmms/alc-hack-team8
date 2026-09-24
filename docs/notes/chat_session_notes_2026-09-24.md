# Comprehensive Project & Session Notes: MD-to-INS Simulation

**Date:** 2026-09-24  
**Project:** Simulating Inelastic Neutron Scattering (INS) from Molecular Dynamics (MD) via Intermediate Dynamical Representations  
**Context:** ALC Hack Team 8  
**Reference Papers:**
1. **Cheng et al.** (*J. Chem. Theory Comput.* 2020, 16, 7702–7708): Isotropic partial phonon density of states (pDOS) bridging MD to quantum 1-phonon cross sections and scalar multiphonon series.
2. **Harrelson et al.** (*Sci. Rep.* 2021, 11, 7938): Anisotropic Cartesian velocity cross-correlation tensors, "almost-isotropic" powder averaging, and overtone convolutions.

---

## 1. Project Motivation & Scope Clarification

### Motivation
- **The Challenge:** Simulating INS spectra for large, disordered, amorphous, or biologically complex systems via ab initio lattice dynamics (LD) scales as $\mathcal{O}(N_a^3)$ due to dynamical matrix diagonalization and fails in the absence of crystalline periodicity.
- **The Solution:** Molecular dynamics scales favorably ($\mathcal{O}(N_a)$ to $\mathcal{O}(N_a \log N_a)$) and naturally samples structural disorder and anharmonicity.
- **The Core Strategy:** Instead of direct 4-point coordinate space-time correlations $S(\mathbf{Q}, \omega) = \mathcal{F}_{\mathbf{r}, t}\{G(\mathbf{r}, t)\}$, use **intermediate dynamical representations** (atom-projected vibrational density of states and $3 \times 3$ Cartesian velocity cross-correlation tensors) to decouple the classical simulation from the quantum scattering statistics.

### Explicit User Constraints & Scoping
- **Drop MLIP Parts:** Drop machine-learned interatomic potential (MLIP) training/fine-tuning, on-the-fly dynamics, and inverse potential fitting.
- **Trajectory Provided:** Assume MD trajectories (ExtXYZ, NetCDF, HDF5, ASE Atoms) are already computed and provided as input (e.g. `trajectories/small-sto-T50-nve.extxyz`).
- **Focus:** Ground strictly in the mathematical formulations of the two papers, derive all equations from first principles, audit the open-source reference code (`MolDyINS`), and design a clean, modular Python architecture.

---

## 2. Critical Audit of Published Papers and Open-Source Code (`MolDyINS`)

### Theoretical Inconsistencies Identified
1. **Missing Zero-Point Factor in Cheng (2020) Eq. (2):**
   Cheng defines $B_{ds} = e_{ds} e_{ds}^T$ and $A_d = \sum_s B_{ds} (2n_s + 1)$ with dimensionless polarization vectors $e_{ds}$, causing $A_d$ to be dimensionless and the Debye–Waller exponent $W_d = \frac{1}{3} Q^2 \operatorname{Tr}[A_d]$ to have unphysical dimensions of $\text{Length}^{-2}$. The physical displacement tensor requires the quantum amplitude factor $\frac{\hbar}{2 M_d \omega_s}$.
2. **Prefactor Ambiguity in Cheng (2020) Eq. (8):**
   An unphysical factor of 3 appears in the numerator ($\frac{3\sigma_d}{2M_d}$ vs $\frac{\sigma_d}{2M_d}$) depending on whether $g_d(\omega)$ is normalized to 1 or 3 (3 degrees of freedom per atom).
3. **Multiphonon Debye–Waller Decoupling:**
   Directly convolving the full single-phonon cross section $S_1(Q, \omega)$ iteratively produces an unphysical decay factor $\exp(-2n W_d)$ instead of the physical global prefactor $\exp(-2W_d)$. The convolution must operate on the normalized spectral distribution.

### Audit of Harrelson Reference Code (`tfharrelson/MolDyINS`)
Our inspection of `MolDyINS.py`, `src/atom.py`, and `src/utils.py` uncovered 6 critical bugs and omissions:

| Flaw / Bug in `MolDyINS` | Code Location | Physical Consequence | Rectification in Our Codebase |
| :--- | :--- | :--- | :--- |
| **Missing Cross Section ($\sigma_i$)** | `MolDyINS.py:182, 185` | Atomic cross section is never multiplied into `currAtom.Slaw`. In multi-element systems, H, C, and O are weighted identically! | Multiply each species by $\frac{\sigma_i}{4\pi}$ (from NIST neutron scattering length tables). |
| **Missing Atomic Mass ($M_i$)** | `src/atom.py:82` | Displacement tensor `f1` denominator omits atom mass $M_i$. | Divide displacement by atomic mass: $\mathbf{B}_i(\omega) \propto \frac{\hbar}{2 M_i \omega} \mathbf{\Gamma}_i(\omega)$. |
| **Hardcoded Proton Mass** | `MolDyINS.py:156` | Equipartition verification test hardcodes proton mass `m_p`. | Use element-specific mass $M_d$ loaded from trajectory metadata. |
| **Non-Hermitian Truncation Suppression** | `src/atom.py:63` | Uses `np.absolute` on trace to suppress complex noise from finite FFT truncation. | Symmetrize $\frac{1}{2}(\mathbf{\Gamma} + \mathbf{\Gamma}^\dagger)$ upfront; eigenvalues and trace are guaranteed real. |
| **Hardcoded Spectrometer Kinematics** | `src/atom.py:97–100` | Hardcodes $E_f = 32\text{ meV}$ and $\theta = 135^\circ$, violating both TOSCA ($3.32\text{ meV}$) and VISION ($3.9\text{ meV}$). | Parameterized instrument classes for VISION, TOSCA, and direct geometry. |
| **Omission of Finite-$T$ Multi-State Convolutions** | `MolDyINS.py:181–186` | Paper Eqs. 18–21 for finite-temperature detailed balance were never implemented in code (strictly 0 K ground-state). | Implement detailed balance factor $\exp(-\beta\hbar\omega/2)$ and Boltzmann multi-state convolution. |

---

## 3. First-Principles Derivation of Key Equations

### A. Intermediate Displacement Tensors from MD Trajectories
From classical trajectories $v_{d,\alpha}(t)$, windowed discrete Fourier transforms give:
$$\tilde{v}_{d,\alpha}(\omega_k) = \sum_{n=0}^{N_t-1} w(t_n) v_{d,\alpha}(t_n) e^{-i\omega_k t_n} \Delta t$$
The cross-spectral density matrix:
$$\mathbf{\Gamma}_{d,\alpha\beta}(\omega) = \frac{1}{2\pi T_f} \tilde{v}_{d,\alpha}^*(\omega) \tilde{v}_{d,\beta}(\omega) \in \mathbb{C}^{3 \times 3}, \quad \mathbf{\Gamma}_d^{\text{sym}}(\omega) = \frac{1}{2}(\mathbf{\Gamma}_d(\omega) + \mathbf{\Gamma}_d(\omega)^\dagger)$$
Equipartition theorem (Smith, *Elements of Molecular Dynamics*):
$$\frac{1}{2\pi} \int_0^\infty \operatorname{Tr}[\mathbf{\Gamma}_d(\omega)] d\omega = \frac{3 k_B T_{\text{MD}}}{M_d}$$
Mapping to quantum harmonic displacements (Balucani & Zoppi, *Dynamics of the Liquid State*):
$$\mathbf{B}_d(\omega) = \frac{\hbar}{2 M_d \omega} \coth\left(\frac{\hbar \omega}{2 k_B T_{\text{exp}}}\right) \frac{\operatorname{Re}[\mathbf{\Gamma}_d^{\text{sym}}(\omega)]}{k_B T_{\text{MD}} / M_d}$$
Total atomic displacement parameter (ADP) tensor (Frenkel & Smit, *Understanding Molecular Simulation*):
$$\mathbf{A}_d = \int_0^\infty \mathbf{B}_d(\omega) d\omega = \langle \mathbf{u}_d \mathbf{u}_d^T \rangle$$

---

### B. Derivation of Harrelson Equation (11) and (12)
The powder-averaged single-phonon cross section over the unit sphere $S^2$ is:
$$\langle S_{0\to 1}(q, \omega_j) \rangle = \frac{1}{4\pi} \int_{S^2} S_{0\to 1}(q\hat{\mathbf{n}}, \omega_j) d\Omega = \sum_i \sigma_i q^2 \langle (\hat{\mathbf{n}}^T \mathbf{B}_{ij} \hat{\mathbf{n}}) \exp(-q^2 \hat{\mathbf{n}}^T \mathbf{A}_i \hat{\mathbf{n}}) \rangle_{S^2}$$

#### Three Independent Proof Methods:
1. **Method 1: Cartesian Moment Contraction & Cumulant Resummation**
   - The $\mathbf{B}$-weighted expectation value is:
     $$\langle e^{-q^2 \hat{\mathbf{n}}^T \mathbf{A} \hat{\mathbf{n}}} \rangle_{\mathbf{B}} \equiv \frac{\langle (\hat{\mathbf{n}}^T \mathbf{B} \hat{\mathbf{n}}) e^{-q^2 \hat{\mathbf{n}}^T \mathbf{A} \hat{\mathbf{n}}} \rangle_{S^2}}{\frac{1}{3}\operatorname{Tr}(\mathbf{B})}$$
   - The unweighted 2nd moment is $\langle n_\alpha n_\beta \rangle_{S^2} = \frac{1}{3}\delta_{\alpha\beta}$.
   - The 4th moment by $SO(3)$ isotropic invariance is:
     $$\langle n_\alpha n_\beta n_\gamma n_\delta \rangle_{S^2} = \frac{1}{15} (\delta_{\alpha\beta}\delta_{\gamma\delta} + \delta_{\alpha\gamma}\delta_{\beta\delta} + \delta_{\alpha\delta}\delta_{\beta\gamma})$$
   - Contracting with symmetric tensors $B_{\alpha\beta} A_{\gamma\delta}$:
     $$\langle (\hat{\mathbf{n}}^T \mathbf{B} \hat{\mathbf{n}})(\hat{\mathbf{n}}^T \mathbf{A} \hat{\mathbf{n}}) \rangle_{S^2} = \frac{1}{15} [\operatorname{Tr}(\mathbf{B})\operatorname{Tr}(\mathbf{A}) + 2 (\mathbf{B}:\mathbf{A})]$$
   - Dividing by $\frac{1}{3}\operatorname{Tr}(\mathbf{B})$ yields the first cumulant $\alpha$:
     $$\alpha = \frac{1}{5} \left[ \operatorname{Tr}(\mathbf{A}) + 2 \left( \frac{\mathbf{B}:\mathbf{A}}{\operatorname{Tr}(\mathbf{B})} \right) \right] \quad \text{\textbf{(Harrelson Eq. 12)}}$$
   - Re-exponentiating gives:
     $$\langle S_{0\to 1}(q, \omega_j) \rangle \approx \sum_i \frac{q^2}{3} \operatorname{Tr}(\mathbf{B}_{ij}) \exp(-q^2 \alpha_{ij}) \quad \text{\textbf{(Harrelson Eq. 11)}}$$

2. **Method 2: Irreducible Spherical Tensor Decomposition (Wigner–Eckart)**
   - Decompose symmetric tensors into $L=0$ monopoles and $L=2$ quadrupoles:
     $$\hat{\mathbf{n}}^T \mathbf{T} \hat{\mathbf{n}} = \frac{1}{3}\operatorname{Tr}(\mathbf{T}) + \sqrt{\frac{8\pi}{15}} \sum_{m=-2}^2 T^{(2)}_m Y_{2m}^*(\hat{\mathbf{n}})$$
   - Orthogonality $\int Y_{lm} Y_{l'm'}^* d\Omega = \delta_{ll'}\delta_{mm'}$ forces all cross-terms between $l=0$ and $l=2$ to vanish, proving why only the scalar trace and the quadrupole contraction $\frac{2}{15}(\widetilde{\mathbf{B}}:\widetilde{\mathbf{A}})$ appear.

3. **Method 3: Generating Functional & Source-Derivative Method**
   - Differentiate the sphere partition function $\mathcal{Z}(\mathbf{J}) = \frac{1}{4\pi}\int_{S^2} e^{-\hat{\mathbf{n}}^T \mathbf{J} \hat{\mathbf{n}}} d\Omega$ with respect to source tensor $\mathbf{J}$ at $\mathbf{J} = q^2 \mathbf{A}$.

---

### C. Exact Error Analysis & Physical Breakdown (The Second Cumulant $\kappa_2$)
- Calculating the 6th-order spherical moment tensor:
  $$\langle n_\alpha n_\beta n_\gamma n_\delta n_\mu n_\nu \rangle_{S^2} = \frac{1}{105} \sum_{15 \text{ pairs}} \delta_{\cdot\cdot}\delta_{\cdot\cdot}\delta_{\cdot\cdot}$$
  yields the exact variance / second cumulant:
  $$\kappa_2 = \frac{1}{35} \left[ \operatorname{Tr}(\mathbf{A})^2 + 2 \operatorname{Tr}(\mathbf{A}^2) + 4 \operatorname{Tr}(\mathbf{A}) \frac{\mathbf{B}:\mathbf{A}}{\operatorname{Tr}(\mathbf{B})} + 8 \frac{\operatorname{Tr}(\mathbf{B}\mathbf{A}^2)}{\operatorname{Tr}(\mathbf{B})} \right] - \alpha^2$$
- **Isotropic Limit:** $\kappa_2 = 0$ identically (the isotropic powder average has zero variance).
- **Anisotropic Materials (P3HT):** When $\Delta A \sim 0.1\text{ \AA}^2$, the approximation is highly accurate up to $q \approx 4\text{ \AA}^{-1}$. At high momentum transfers ($q > 6\text{ \AA}^{-1}$ on VISION/TOSCA), higher-order cumulants produce non-negligible damping.

---

### D. Multiphonon Series & Spherical Contraction Factor $(3/5)^{n-1}$
The $n$-th order overtone excitation is given by:
$$S_{n,d}(Q, \omega) = \frac{\sigma_d}{4\pi} \frac{Q^{2n}}{3 \cdot n!} \left(\frac{3}{5}\right)^{n-1} \operatorname{Tr}[\mathbf{b}_d^{*n}(\omega)] \exp\left( -\frac{1}{3} Q^2 \operatorname{Tr}(\mathbf{A}_d) \right)$$
- The prefactor $(3/5)^{n-1} / (3 \cdot n!)$ is derived from the $n$-th order spherical contraction of $(\hat{\mathbf{n}}^T \mathbf{B} \hat{\mathbf{n}})^n$.
- Evaluated efficiently via fast Fourier transform:
  $$\mathbf{b}_d^{*n}(\omega) = \mathcal{F}^{-1}\left\{ (\mathcal{F}\{\mathbf{B}_d(\omega)\})^n \right\}$$

---

## 4. Software Architecture Blueprint (4-Phase Implementation)

```
md_ins/
├── io/
│   ├── trajectory.py        # Stream ExtXYZ, NetCDF, HDF5, ASE Atoms
│   └── cross_sections.py    # NIST nuclear scattering lengths (sigma_inc, sigma_coh, M)
├── intermediate/
│   ├── vacf.py              # Windowed FFT VACF & isotropic pDOS g_d(w) (Cheng 2020)
│   ├── tensors.py           # Strictly Hermitian Gamma_d(w) & displacement B_d(w) (Harrelson 2021)
│   └── debye_waller.py      # Isotropic <u^2> and anisotropic tensor A_d
├── scattering/
│   ├── detailed_balance.py  # Quantum Bose & detailed-balance factors
│   ├── powder_average.py    # Isotropic (Cheng) & Almost-Isotropic (Harrelson) contractions
│   ├── cumulants.py         # 2nd cumulant kappa_2 calculation & error diagnostics
│   └── multiphonon.py       # FFT polynomial overtone convolutions (n=1..10)
├── instruments/
│   ├── kinematics.py        # Q(w) relations for TOSCA, VISION, and direct-geometry
│   └── resolution.py        # Energy-dependent instrumental broadening kernels (Delta E / E)
└── pipeline.py              # End-to-end runner: Trajectory -> Intermediate -> S(Q,w) -> Spectrum
```

---

## 5. Artifacts and Generated Reference Files

All mathematical derivations, reports, and compiled PDFs generated:

| File | Type | Description |
| :--- | :--- | :--- |
| [`derivation_eq_11_take3.pdf`](file:///home/drFaustroll/playground/ml/alc-hack-team8/derivation_eq_11_take3.pdf) | PDF (13 pages) | **Take 3 (Definitive)**: Unified theoretical and numerical treatise with Lebedev spherical quadrature benchmarking, $\text{SrTiO}_3$ vs. P3HT case studies. |
| [`derivation_eq_11_take3.tex`](file:///home/drFaustroll/playground/ml/alc-hack-team8/derivation_eq_11_take3.tex) | LaTeX Source | Complete LaTeX source for the 13-page Take 3 treatise. |
| [`docs/prds/derivation-equation-11-2026-09-24.md`](file:///home/drFaustroll/playground/ml/alc-hack-team8/docs/prds/derivation-equation-11-2026-09-24.md) | PRD Markdown | Comprehensive Product Requirements Document generated via `eg-prd`. |
| [`derivation_eq_11_take2.pdf`](file:///home/drFaustroll/playground/ml/alc-hack-team8/derivation_eq_11_take2.pdf) | PDF (12 pages) | Comprehensive multi-perspective derivation of Eq. 11 & 12 (Take 2). |
| [`derivation_eq_11_take2.tex`](file:///home/drFaustroll/playground/ml/alc-hack-team8/derivation_eq_11_take2.tex) | LaTeX Source | LaTeX source for Take 2. |
| [`docs/derivation_equation_11.pdf`](file:///home/drFaustroll/playground/ml/alc-hack-team8/docs/derivation_equation_11.pdf) | PDF (6 pages) | First iteration derivation of Harrelson Eq. (11). |
| [`docs/derivations_and_concrete_steps.md`](file:///home/drFaustroll/playground/ml/alc-hack-team8/docs/derivations_and_concrete_steps.md) | Markdown | 4-part technical blueprint and code audit of `MolDyINS`. |
| [`docs/proposals/original_proposal_full.md`](file:///home/drFaustroll/playground/ml/alc-hack-team8/docs/proposals/original_proposal_full.md) | Markdown | Preserved original proposal before simplification. |

