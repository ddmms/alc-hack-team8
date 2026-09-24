# Proposal: Simulation of Inelastic Neutron Scattering (INS) from Molecular Dynamics Trajectories via Intermediate Dynamical Representations

## 1. Executive Summary & Project Context

Inelastic neutron scattering (INS) provides an exceptionally direct probe of atomic and molecular vibrations because neutrons interact via nuclear potentials rather than electronic charge distributions, yielding spectra free from electromagnetic selection rules. However, simulating INS spectra for large, disordered, amorphous, or biologically complex systems using conventional ab initio lattice dynamics (LD) is computationally prohibitive due to $\mathcal{O}(N_a^3)$ scaling of dynamical matrix diagonalization and the loss of periodic crystalline symmetry.

Classical and ab initio molecular dynamics (MD) scale favorably ($\mathcal{O}(N_a)$ to $\mathcal{O}(N_a \log N_a)$) and naturally capture structural disorder, large unit cells, and anharmonic sampling of the potential energy surface. However, generating INS spectra directly from MD trajectories presents two major hurdles:
1. **Classical vs. Quantum Discrepancy**: Classical trajectories follow the equipartition theorem, underestimating zero-point motion and producing mode amplitudes proportional to the classical kinetic temperature $T_{\text{MD}}$ rather than following quantum Bose-Einstein statistics ($n(\omega) + 1$).
2. **Kinematic & Multiphonon Modeling**: Standard velocity autocorrelation power spectra fail to model neutron momentum transfer ($Q$-dependence), the Debye-Waller attenuation factor ($\exp(-2W)$), instrument trajectories ($Q(\omega)$ on inverted-geometry spectrometers like VISION or TOSCA), and higher-order multiphonon transitions (overtones and combination bands).

This project implements an open-source, modular Python framework to simulate 1-phonon and multiphonon INS spectra from MD trajectories by using **intermediate dynamical representations**—specifically, atom-projected phonon densities of states (pDOS) and Cartesian velocity correlation tensors—as formulated in two foundational papers:
- **Cheng et al. (JCTC 2020)**: Fully isotropic approximation bridging MD-derived atom-projected pDOS to quantum-corrected 1-phonon cross-sections and scalar multiphonon convolutions.
- **Harrelson et al. (Sci. Rep. 2021)**: Anisotropic framework employing a $3 \times 3$ Cartesian outer-product velocity correlation tensor per atom, the "almost isotropic" powder averaging approximation, and matrix-valued overtone convolutions.

---

## 2. Critical Scientific Assessment of the Two Papers

### 2.1 Cheng et al. (*J. Chem. Theory Comput.* 2020, 16, 7702–7708)

#### Core Philosophy & Merits
Cheng et al. integrate MD trajectories into the theoretical framework of the **OCLIMAX** code. The core idea is to extract the atom-projected density of states (pDOS) $g_d(\omega)$ from the trace of the velocity autocorrelation function (VACF) and map it into the standard lattice dynamics incoherent scattering formula.
- **Decoupling MD Temperature from Experimental Temperature**: A major scientific strength of Cheng et al. is recognizing that classical MD at a "fictive" temperature $T_{\text{MD}}$ is only used to survey the potential energy landscape and obtain the vibrational spectrum $g_d(\omega)$. The quantum thermal occupation $n(\omega) = [\exp(\hbar\omega / k_B T_{\text{exp}}) - 1]^{-1}$ and the quantum Debye-Waller factor $W_d(Q)$ are then evaluated analytically at the experimental temperature $T_{\text{exp}}$, correctly recovering zero-point fluctuations.
- **Demonstrated Applicability**: Tested successfully against experimental data from the VISION (indirect geometry) and SEQUOIA (direct geometry) spectrometers at ORNL on hydrogen-disordered ice Ih (384–1296 atoms), MOF ZIF-8 (276 atoms), and amorphous silica glass (3000 atoms).

#### Theoretical & Mathematical Issues / Ambiguities
1. **Dimensional Inconsistency & Missing $\hbar / (2 M_d \omega)$ Factor in Displacement Tensors**:
   - In Equation (2), Cheng et al. define:
     $$B_{ds} = e_{ds} e_{ds}^T, \quad A_d = \sum_s B_{ds} (2 n_s + 1)$$
     and in Equation (4):
     $$W_d = \frac{1}{3} Q^2 \text{Tr}[A_d]$$
   - If the polarization vectors $e_{ds}$ are dimensionless unit vectors, $B_{ds}$ and $A_d$ are strictly dimensionless. Consequently, $W_d \propto Q^2$ has dimensions of $\text{Length}^{-2}$ (e.g., $\text{Å}^{-2}$), rendering the Debye-Waller exponent $\exp(-2W_d)$ dimensionally non-sensical!
   - In standard lattice dynamics (Squires, Lovesey, and Cheng's own 2019 OCLIMAX paper), the displacement tensor includes the quantum zero-point amplitude factor $\frac{\hbar}{2 M_d \omega_s}$:
     $$A_d = \sum_s \frac{\hbar}{2 M_d \omega_s} e_{ds} e_{ds}^T (2 n_s + 1)$$
     This factor was inadvertently omitted in Eq. (2) of Cheng (2020).
2. **Prefactor Ambiguity in Eq. (8) vs Eq. (3) & (7)**:
   - In Eq. (3), the mode sum contains $\sum_s \frac{Q^2 \text{Tr}(B_{ds})}{3 \omega_s}$.
   - In Eq. (7), the partial DOS is defined via $g_d(\omega) d\omega = \sum_{\omega \le \omega_s < \omega + d\omega} \frac{\text{Tr}(B_{ds})}{3}$, which normalizes $\int g_d(\omega) d\omega = 1$.
   - However, in Eq. (8), Cheng writes:
     $$S_{\text{inc}\pm 1}(Q, \omega) = \sum_d \frac{3 \sigma_d}{2 M_d} \exp(-2 W_d) \frac{Q^2 g_d(\omega)}{\omega} \left(n + \frac{1}{2} \pm \frac{1}{2}\right)$$
     An extra factor of 3 appears in the numerator ($\frac{3\sigma_d}{2 M_d}$ vs $\frac{\sigma_d}{2 M_d}$). Depending on whether $g_d(\omega)$ is normalized to 1 or 3 (3 degrees of freedom per atom), this represents an internal scaling inconsistency that must be standardized in implementation.
3. **Unspecified Multiphonon Formulation**:
   - Eq. (9) presents an iterative convolution:
     $$S_n(\omega) = \int_{-\infty}^\infty S_1(\omega - \omega') S_{n-1}(\omega') d\omega'$$
   - Convolving the already-evaluated $S_1(Q, \omega)$ without decoupling the Debye-Waller factor would result in an erroneous factor of $\exp(-2n W_d)$ instead of the physical overall factor $\exp(-2W_d)$. The actual multiphonon series expansion requires convolving the normalized spectral distribution $[\frac{g_d(\omega)}{\omega (1 - e^{-\beta\hbar\omega})}]$ with prefactors $\frac{Q^{2n}}{n!}$, which is handled inside OCLIMAX Fortran routines but glossed over in the text.
4. **Isotropic Limitation**:
   - The method assumes an isotropic displacement tensor ($\text{Tr}(B_{ds}) / 3 \cdot \mathbf{I}$). For highly anisotropic systems (such as aligned conjugated polymers or layered 2D materials), orientation-dependent excitation probabilities cannot be captured.

---

### 2.2 Harrelson et al. (*Sci. Rep.* 2021, 11, 7938)

#### Core Philosophy & Merits
Harrelson et al. aim to eliminate the isotropic approximation by retaining the directional Cartesian information from MD trajectories without performing full normal mode analysis:
- **Cartesian Velocity Correlation Tensor**: Rather than reducing the trajectory to a scalar VACF, they compute the $3 \times 3$ outer-product velocity correlation tensor per atom:
  $$\mathbf{C}_i(t) = \langle \vec{v}_i(t) \vec{v}_i^T(t + \tau) \rangle \xrightarrow{\mathcal{F}} \vec{v}\vec{v}^T(\omega)$$
- **"Almost Isotropic" Powder Averaging**: Adopts the analytical powder-averaging expansion of Ramirez-Cuesta (aCLIMAX 2004), which accounts for the projection of the mode displacement tensor $B_{ij}$ along the total anisotropic thermal displacement ellipsoid $A_i$.
- **Equipartition-to-Quantum Mapping**: Connects classical velocity spectral power to quantum ground-state displacement via:
  $$u_i^2(\omega) = \frac{\hbar v_i^2(f)}{\omega k_B T_{\text{MD}}}$$
- **Demonstrated on P3HT**: Applied to semicrystalline and amorphous poly(3-hexylthiophene) (P3HT), demonstrating that MD captures low-frequency disorder and conformational breadth better than single-conformation DFT.

#### Analysis of Supporting Information (SI)
- **Cheng et al. (2020)**: Published with **no Supporting Information** (the article contains no Associated Content; all calculations relied directly on OCLIMAX Fortran internals).
- **Harrelson et al. (2021)**: Published an official 8-page Supporting Information document (`harrelson2021_SI.pdf`) comprising six detailed sections:
  1. *Section S1 (Eqs. S1–S16)*: Normal mode decomposition of MD velocities, proving random phase cancellation for degenerate modes.
  2. *Section S2 (Eqs. S17–S24)*: Operator derivation connecting the convolution series $\hat{\mathcal{C}}^n \mathbf{B}(\omega)$ to the quantum intermediate scattering function $I_i(\vec{q}, t) = \sigma_i \langle e^{-i\vec{q}\cdot\vec{r}_i(0)} e^{i\vec{q}\cdot\vec{r}_i(t)} \rangle$.
  3. *Section S3 (Eqs. S25–S35)*: Second-order perturbation theory relating classical mode relaxation rates $\langle R_i \rangle$ to quantum anharmonic energy level shifts $\Delta E_i^{(2)} = \frac{\pi\hbar}{4}\langle R_i \rangle$.
  4. *Section S4 (Figs. S2–S3)*: Assessment of thermostat and barostat time constants (recommends $\tau = 10$ ps to avoid unphysical spectral spikes at low frequency).
  5. *Section S5 & S6 (Figs. S4–S5)*: Experimental P3HT comparisons and verification of the overtone background in the 1600–3000 cm$^{-1}$ gap.

#### Reconciliation with Open-Source Code (`MolDyINS`)
Direct inspection of the authors' codebase (`tfharrelson/MolDyINS`: `MolDyINS.py`, `src/atom.py`, `src/utils.py`) yielded critical reconciliations:
1. **Critical Code Bug: Nuclear Scattering Cross Section ($\sigma_i$) Omitted**:
   In `MolDyINS.py` (line 88), `currAtom.assignXS(...)` reads $\sigma_i$, but in the scattering summation loop (`MolDyINS.py:181–186`), `currAtom.Slaw` is updated **without ever multiplying by `self.XS`**! For pure P3HT ($>99\%$ scattering from $^1$H), this functioned as an arbitrary global scale factor, but for multi-element materials (ice, hydrates, MOFs), this omission produces invalid relative element intensities. Our implementation strictly multiplies by $\frac{\sigma_i}{2M_i}$.
2. **Resolution of Typo in Eq. (12)**:
   In `src/atom.py` (line 63), the code computes `alpha = 0.2*(np.trace(A) + 2*np.absolute(np.trace(np.dot(A,B))) / np.absolute(np.trace(B)))`. This confirms the true intention was $\mathbf{B}_i(\omega) : \mathbf{A}_i = \text{Tr}(\mathbf{A}_i \mathbf{B}_i(\omega))$ normalized by $\text{Tr}(\mathbf{B}_i(\omega))$, resolving the typographical ambiguities in the published paper.
3. **Omission of Finite-Temperature Boltzmann Convolution**:
   There is **no code** implementing Eqs. (18–21) in `MolDyINS`. The code executes strictly the 0 K ground-state formulation.
4. **Omission of SI Anharmonic Energy Shifts**:
   The second-order perturbation energy corrections described in SI Section S3 were never integrated into `MolDyINS.py`; overtones remain harmonic convolutions.
5. **Non-Hermitian Convolutions & Absolute Value Workaround**:
   Matrix multiplication of frequency functions in `src/atom.py:93` generates non-Hermitian complex components, which the author discarded in `MolDyINS.py:182` by taking `np.absolute(...)` of the trace. Our framework rectifies this by enforcing matrix symmetrization $\frac{1}{2}(\mathbf{B}_1\mathbf{B}_2 + \mathbf{B}_2\mathbf{B}_1)$.

---

## 3. Comparison and Synthesis of Approaches

| Feature / Property | Cheng et al. (JCTC 2020) | Harrelson et al. (Sci. Rep. 2021) | Proposed Implementation Standard |
| :--- | :--- | :--- | :--- |
| **Intermediate Representation** | Scalar atom-projected pDOS $g_d(\omega)$ (1D vector per atom) | $3 \times 3$ Cartesian velocity correlation tensor $\mathbf{B}_i(\omega)$ | Modular support for both scalar pDOS and $3 \times 3$ tensor |
| **Vibrational Anisotropy** | Fully isotropic ($\text{Tr}(\mathbf{B})/3$) | Directional / "Almost Isotropic" | Evaluates both isotropic and almost-isotropic powder averages |
| **Temperature Decoupling** | Rigorous: classical MD for potential survey; quantum Bose statistics & DWF at $T_{\text{exp}}$ | Approximate: ground-state zero-point mapping ($T \to 0$ K); ad-hoc finite-$T$ convolution | Rigorous quantum thermal population $(n(\omega)+1)$ and analytic quantum DWF |
| **Debye-Waller Factor** | Analytic quantum expectation $W_d = \frac{1}{3} Q^2 \text{Tr}(\mathbf{A}_d)$ | Powder-averaged directional exponent $\alpha_i(\omega)$ for $n=1$, isotropic for $n \ge 2$ | Corrected analytic quantum $\mathbf{A}_d(T_{\text{exp}})$ with corrected $\alpha_i(\omega)$ formula |
| **Multiphonon Handling** | Scalar 1D iterative convolutions | $3 \times 3$ matrix convolutions with inner-index contraction | Scalar 1D FFT convolutions + optional symmetrized tensor convolutions |
| **Instrument Kinematics** | General direct/indirect ($Q(\omega)$ parameterization) | Specific inverted geometry (VISION backscattering: $E_f \approx 32\text{ cm}^{-1}$) | Explicit inverted-geometry kinematics ($E_f = 3.97$ meV) & energy-dependent resolution |

---

## 4. Unified Mathematical Framework for Implementation

### Stage 1: Trajectory Processing & Intermediate Representations
From an MD trajectory with velocity time series $\vec{v}_d(t)$ for atom $d$:
1. **Cartesian Velocity Autocorrelation Tensor**:
   $$\mathbf{C}_d(\tau) = \langle \vec{v}_d(t) \vec{v}_d^T(t + \tau) \rangle_t$$
   Fourier transform (Wiener-Khinchin):
   $$\tilde{\mathbf{C}}_d(\omega) = \int_{-\infty}^\infty \mathbf{C}_d(\tau) e^{-i\omega \tau} d\tau$$
   By time-reversal symmetry of equilibrium classical MD, $\tilde{\mathbf{C}}_d(\omega)$ is real and symmetric. To suppress finite-time noise, windowing (Hann or Welch) and symmetrization $\frac{1}{2}(\tilde{\mathbf{C}} + \tilde{\mathbf{C}}^T)$ are applied.
2. **Scalar Projected pDOS**:
   $$g_d(\omega) = \frac{\text{Tr}[\tilde{\mathbf{C}}_d(\omega)]}{\int_0^\infty \text{Tr}[\tilde{\mathbf{C}}_d(\omega')] d\omega'}$$
   normalized such that $\int_0^\infty g_d(\omega) d\omega = 1$.

### Stage 2: Quantum Corrections & Thermal Displacement Tensors
The classical-to-quantum bridge maps MD power to quantum displacements at target experimental temperature $T_{\text{exp}}$:
1. **Quantum Thermal Amplitude**:
   For mode energy $\hbar\omega$, the mean-square displacement in quantum statistical mechanics is:
   $$\langle u^2(\omega) \rangle = \frac{\hbar}{2 M_d \omega} \coth\left( \frac{\hbar\omega}{2 k_B T_{\text{exp}}} \right) = \frac{\hbar}{M_d \omega} \left( n(\omega, T_{\text{exp}}) + \frac{1}{2} \right)$$
2. **Frequency-Resolved Displacement Tensor $\mathbf{B}_d(\omega)$**:
   $$\mathbf{B}_d(\omega) = \frac{\hbar}{2 M_d \omega} \frac{\tilde{\mathbf{C}}_d(\omega)}{\frac{1}{3} \int_0^\infty \text{Tr}[\tilde{\mathbf{C}}_d(\omega')] d\omega'}$$
3. **Total Displacement Tensor $\mathbf{A}_d$ & Debye-Waller Factor**:
   $$\mathbf{A}_d = \int_0^\infty \mathbf{B}_d(\omega) (2 n(\omega, T_{\text{exp}}) + 1) d\omega$$
   Isotropic mean-square displacement:
   $$\langle u_d^2 \rangle = \text{Tr}(\mathbf{A}_d)$$
   Isotropic Debye-Waller factor:
   $$W_{d,\text{iso}}(Q) = \frac{1}{6} Q^2 \langle u_d^2 \rangle = \frac{1}{6} Q^2 \text{Tr}(\mathbf{A}_d)$$

### Stage 3: Fundamental (1-Phonon) INS Cross Section
1. **Method A: Fully Isotropic (Cheng et al.)**:
   $$S_{1,\text{iso}}(Q, \omega) = \sum_d \frac{\sigma_d}{2 M_d} \frac{Q^2 g_d(\omega)}{\omega} (n(\omega, T_{\text{exp}}) + 1) \exp\left( - \frac{1}{3} Q^2 \text{Tr}(\mathbf{A}_d) \right)$$
2. **Method B: Anisotropic / "Almost Isotropic" (Harrelson et al., Corrected)**:
   $$S_{1,\text{aniso}}(Q, \omega) = \sum_d \frac{\sigma_d}{2 M_d} \frac{Q^2 \text{Tr}(\mathbf{B}_d(\omega))}{\omega} (n(\omega, T_{\text{exp}}) + 1) \exp\left( - Q^2 \alpha_d(\omega) \right)$$
   where:
   $$\alpha_d(\omega) = \frac{1}{5} \left[ \text{Tr}(\mathbf{A}_d) + 2 \frac{\mathbf{B}_d(\omega) : \mathbf{A}_d}{\text{Tr}(\mathbf{B}_d(\omega))} \right]$$

### Stage 4: Multiphonon Excitations (Overtones & Combinations)
Incoherent multiphonon transitions are evaluated via convolution of the normalized 1-phonon spectral profile:
$$\tilde{S}_{1,d}(\omega) = \frac{g_d(\omega)}{\omega} (n(\omega, T_{\text{exp}}) + 1)$$
Higher-order transitions:
$$\tilde{S}_{n,d}(\omega) = \int_0^\omega \tilde{S}_{1,d}(\omega') \tilde{S}_{n-1,d}(\omega - \omega') d\omega'$$
Total scattering function:
$$S(Q, \omega) = \sum_d \frac{\sigma_d}{2 M_d} \exp(-2 W_d(Q)) \sum_{n=1}^{N_{\text{max}}} \frac{Q^{2n}}{n!} \tilde{S}_{n,d}(\omega)$$

### Stage 5: Inverted Geometry Spectrometer Kinematics & Resolution
For instruments like VISION (SNS) or TOSCA (ISIS):
- Fixed analyzer energy: $E_f = 3.97$ meV ($32\text{ cm}^{-1}$)
- Incident energy: $E_i = E_f + \hbar\omega$
- Wavevectors: $k_i = \sqrt{2 m_n E_i} / \hbar$, $k_f = \sqrt{2 m_n E_f} / \hbar$
- Momentum transfer trajectory (average scattering angle $2\theta \approx 90^\circ$ or backscattering $135^\circ$):
  $$Q^2(\omega) = k_i^2 + k_f^2 - 2 k_i k_f \cos(2\theta)$$
- Instrument resolution convolution:
  $$I_{\text{conv}}(\omega) = \int S(Q(\omega'), \omega') \mathcal{G}(\omega - \omega'; \sigma(\omega')) d\omega'$$
  with $\sigma(E) = a + b E$ (for VISION, $\sigma = 1.21\text{ cm}^{-1} + 0.01 E$).

---

## 5. Software Architecture & Implementation Plan

### 5.1 Proposed Package Layout (`md_ins/`)
```
md_ins/
├── __init__.py
├── trajectory/
│   ├── __init__.py
│   ├── reader.py          # ASE-compatible trajectory reader (ExtXYZ, NetCDF, Gromacs, LAMMPS)
│   └── vacf.py            # Fast FFT-based VACF and 3x3 Cartesian velocity correlation tensor
├── intermediate/
│   ├── __init__.py
│   ├── pdos.py            # Scalar projected pDOS computation & normalization
│   └── tensor.py          # Cartesian B(omega) and total thermal displacement A tensors
├── physics/
│   ├── __init__.py
│   ├── cross_section.py   # 1-phonon cross section (isotropic and almost-isotropic)
│   ├── multiphonon.py     # Multiphonon expansion via FFT convolutions
│   └── bose.py            # Bose-Einstein populations & quantum Debye-Waller factors
├── instruments/
│   ├── __init__.py
│   ├── vision.py          # VISION spectrometer Q(omega) trajectory and resolution function
│   └── tosca.py           # TOSCA spectrometer parameters
└── validation/
    ├── __init__.py
    └── euphonic_bench.py  # Benchmark harness comparing MD pDOS against Euphonic lattice dynamics
```

### 5.2 Phased Milestones

#### Phase 1: Trajectory Processing & Harmonic Lattice Dynamics Validation
- Implement `md_ins.trajectory.vacf` and `md_ins.intermediate.pdos`.
- Generate harmonic MD trajectories of crystalline Argon (Lennard-Jones FCC) using ASE.
- Compute the exact harmonic pDOS and phonon dispersion using **Euphonic** with the matching LJ potential.
- Benchmark MD-derived pDOS against Euphonic's analytical lattice dynamics result to validate convergence with trajectory length and windowing.

#### Phase 2: Fully Isotropic INS Engine (Cheng 2020 Formulation)
- Implement `md_ins.physics.bose` and `md_ins.physics.cross_section` (Method A).
- Implement scalar multiphonon FFT convolution (`md_ins.physics.multiphonon`).
- Implement VISION inverted-geometry kinematics ($Q(\omega)$) and energy-dependent Gaussian broadening.
- Benchmark simulated 1-phonon and multiphonon spectra against Euphonic's native INS module.

#### Phase 3: Anisotropic & "Almost Isotropic" Tensor Engine (Harrelson 2021 Formulation)
- Implement `md_ins.intermediate.tensor` for full $3 \times 3$ Cartesian outer-product correlation tensors.
- Implement the corrected "almost isotropic" powder average formula.
- Implement matrix-valued overtone convolutions with Hermiticity preservation.
- Compare isotropic vs. almost-isotropic predictions on an anisotropic test system (e.g., planar molecule or polymer chain).

#### Phase 4: Benchmarking & Integration
- Benchmark computational scaling and memory usage with system size $N_a$.
- Test on published benchmark datasets (e.g., Ice Ih or P3HT).
- Package as an installable Python library with CLI and OpenSpec compliance.

