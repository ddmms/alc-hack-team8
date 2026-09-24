# Comprehensive Derivation & Concrete Implementation Blueprint
## Simulation of Inelastic Neutron Scattering from Molecular Dynamics via Intermediate Dynamic Representations

**References:**
- Cheng, Kolesnikov, & Ramirez-Cuesta, *J. Chem. Theory Comput.* 16, 3190–3197 (2020) (`Cheng2020.pdf`)
- Harrelson et al., *Sci. Rep.* 11, 7938 (2021) (`harrelson2021.pdf` and `harrelson2021_SI.pdf`)
- Open-Source Reference Code: `MolDyINS` (`tfharrelson/MolDyINS`)

---

## 1. Executive Mathematical Overview & Definitions

The central premise of both papers is to replace the direct $O(N^2)$ evaluation of the dynamic structure factor:
$$S(\mathbf{Q}, \omega) = \frac{1}{2\pi \hbar} \int_{-\infty}^{\infty} dt \, e^{-i\omega t} \sum_{j, l} b_j b_l \langle e^{-i \mathbf{Q} \cdot \hat{\mathbf{r}}_j(0)} e^{i \mathbf{Q} \cdot \hat{\mathbf{r}}_l(t)} \rangle$$
with **single-site intermediate representations**:
1. **Isotropic pDOS:** $g_d(\omega)$ (Cheng et al. 2020).
2. **Anisotropic Cartesian Cross-Spectral Tensors:** $\mathbf{\Gamma}_d(\omega)$ and displacement tensors $\mathbf{B}_d(\omega)$ (Harrelson et al. 2021).

These intermediate representations compress $10^5$-frame trajectories into compact 1D spectral arrays, enabling exact analytical powder averaging, multiphonon overtone convolutions, and instrument resolution broadening in seconds.

---

## 2. Rigorous Formula Derivations

### 2.1 Velocity Autocorrelation & Cross-Spectral Density Tensor $\mathbf{\Gamma}_d(\omega)$
Let $\mathbf{v}_d(t) \in \mathbb{R}^3$ be the velocity of atom $d$ at time $t$ along Cartesian axes $\alpha, \beta \in \{x, y, z\}$.

The time-domain Cartesian velocity cross-correlation tensor is:
$$\mathbf{C}_{d, \alpha\beta}(t) = \langle v_{d,\alpha}(0) v_{d,\beta}(t) \rangle = \lim_{T \to \infty} \frac{1}{T} \int_0^T v_{d,\alpha}(\tau) v_{d,\beta}(\tau + t) d\tau$$

By the Wiener–Khinchin theorem, its Fourier transform yields the velocity cross-spectral density tensor:
$$\mathbf{\Gamma}_{d, \alpha\beta}(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \mathbf{C}_{d, \alpha\beta}(t) e^{-i\omega t} dt$$

#### Fast Fourier Transform Formulation (Welch Windowing):
For a discrete trajectory of $N_t$ frames with timestep $\Delta t$ ($T_f = N_t \Delta t$), the discrete Fourier transform of the windowed velocity is:
$$\tilde{v}_{d, \alpha}(\omega_k) = \sum_{n=0}^{N_t-1} w(t_n) v_{d, \alpha}(t_n) e^{-i \omega_k t_n} \Delta t$$
where $w(t_n)$ is a normalized window (e.g. Hann window with $\sum_n w^2(t_n) / N_t = 1$).

The spectral tensor is computed via the outer product:
$$\mathbf{\Gamma}_{d, \alpha\beta}(\omega_k) = \frac{1}{2\pi T_f} \tilde{v}_{d, \alpha}^*(\omega_k) \tilde{v}_{d, \beta}(\omega_k)$$

#### Mathematical Properties:
1. **Hermiticity:** $\mathbf{\Gamma}_d(\omega) = \mathbf{\Gamma}_d(\omega)^\dagger$, meaning $\mathbf{\Gamma}_{d, \beta\alpha}(\omega) = \mathbf{\Gamma}_{d, \alpha\beta}^*(\omega)$.
2. **Positive Semi-Definiteness:** $\mathbf{z}^\dagger \mathbf{\Gamma}_d(\omega) \mathbf{z} \ge 0$ for any vector $\mathbf{z} \in \mathbb{C}^3$.
3. **Classical Equipartition Sum Rule:**
   $$\frac{1}{2\pi} \int_0^\infty \operatorname{Tr}[\mathbf{\Gamma}_d(\omega)] d\omega = \langle |\mathbf{v}_d|^2 \rangle = \frac{3 k_B T}{M_d}$$

---

### 2.2 Quantum Scaling & The Displacement Tensor $\mathbf{B}_d(\omega)$
Classical molecular dynamics satisfies classical equipartition where mode energy is $k_B T$. In a quantum harmonic oscillator of frequency $\omega$, the average kinetic energy is:
$$\langle E_{\text{kin}}(\omega) \rangle = \frac{1}{2} \hbar \omega \left[ n(\omega, T) + \frac{1}{2} \right] = \frac{1}{4} \hbar \omega \coth\left(\frac{\hbar \omega}{2 k_B T}\right)$$

To correct for quantum zero-point motion, the classical spectral tensor $\mathbf{\Gamma}_d(\omega)$ is scaled by the harmonic quantum correction factor:
$$Q(\omega, T) = \frac{\hbar \omega}{2 k_B T} \coth\left(\frac{\hbar \omega}{2 k_B T}\right)$$

Relating velocity $\mathbf{v}$ to displacement $\mathbf{u}$ via $\mathbf{v}(\omega) = -i\omega \mathbf{u}(\omega)$, the frequency-dependent atomic displacement tensor is:
$$\mathbf{B}_d(\omega) = \frac{\operatorname{Re}[\mathbf{\Gamma}_d(\omega)]}{\omega^2} \cdot Q(\omega, T) = \frac{\hbar}{2 M_d \omega} \coth\left(\frac{\hbar \omega}{2 k_B T}\right) \operatorname{Re}[\mathbf{\Gamma}_d(\omega)]$$
*Units:* $\text{\AA}^2 \cdot \text{s}$ (or $\text{\AA}^2 / \text{meV}$).

The total mean-squared displacement (Debye-Waller) tensor is:
$$\mathbf{A}_d = \int_0^\infty \mathbf{B}_d(\omega) d\omega = \langle \mathbf{u}_d \mathbf{u}_d^T \rangle \in \mathbb{R}^{3 \times 3}$$

---

### 2.3 Exact Derivation of Harrelson's "Almost-Isotropic" Powder Averaging
In powder neutron scattering, the sample is macroscopically unoriented. The scattering vector $\mathbf{Q} = Q \hat{\mathbf{n}}$ samples all orientations of the unit vector $\hat{\mathbf{n}} \in S^2$ with equal probability measure $d\Omega / 4\pi$.

The powder-averaged 1-phonon cross section requires evaluating:
$$I_1(Q, \omega) = \langle (\mathbf{Q}^T \mathbf{B}_d(\omega) \mathbf{Q}) \exp(-\mathbf{Q}^T \mathbf{A}_d \mathbf{Q}) \rangle_{S^2} = Q^2 \langle (\hat{\mathbf{n}}^T \mathbf{B}_d \hat{\mathbf{n}}) \exp(-Q^2 \hat{\mathbf{n}}^T \mathbf{A}_d \hat{\mathbf{n}}) \rangle_{S^2}$$

#### Step-by-Step Spherical Integration:
Recall the spherical moments of a random unit vector $\hat{\mathbf{n}}$ on $S^2$:
$$\langle \hat{n}_\alpha \hat{n}_\beta \rangle = \frac{1}{3} \delta_{\alpha\beta}$$
$$\langle \hat{n}_\alpha \hat{n}_\beta \hat{n}_\gamma \hat{n}_\delta \rangle = \frac{1}{15} \left( \delta_{\alpha\beta}\delta_{\gamma\delta} + \delta_{\alpha\gamma}\delta_{\beta\delta} + \delta_{\alpha\delta}\delta_{\beta\gamma} \right)$$

Taylor-expanding the Debye-Waller exponential to first order in $Q^2$:
$$\exp(-Q^2 \hat{\mathbf{n}}^T \mathbf{A}_d \hat{\mathbf{n}}) \approx 1 - Q^2 (\hat{\mathbf{n}}^T \mathbf{A}_d \hat{\mathbf{n}}) + \mathcal{O}(Q^4)$$

Substituting this into the expectation value:
$$\langle (\hat{\mathbf{n}}^T \mathbf{B}_d \hat{\mathbf{n}}) \left( 1 - Q^2 \hat{\mathbf{n}}^T \mathbf{A}_d \hat{\mathbf{n}} \right) \rangle = \langle \hat{\mathbf{n}}^T \mathbf{B}_d \hat{\mathbf{n}} \rangle - Q^2 \langle (\hat{\mathbf{n}}^T \mathbf{B}_d \hat{\mathbf{n}})(\hat{\mathbf{n}}^T \mathbf{A}_d \hat{\mathbf{n}}) \rangle$$

1. **First Term:**
   $$\langle \hat{\mathbf{n}}^T \mathbf{B}_d \hat{\mathbf{n}} \rangle = \sum_{\alpha\beta} B_{d,\alpha\beta} \langle \hat{n}_\alpha \hat{n}_\beta \rangle = \frac{1}{3} \operatorname{Tr}(\mathbf{B}_d)$$

2. **Second Term:**
   $$\langle (\hat{\mathbf{n}}^T \mathbf{B}_d \hat{\mathbf{n}})(\hat{\mathbf{n}}^T \mathbf{A}_d \hat{\mathbf{n}}) \rangle = \sum_{\alpha\beta\gamma\delta} B_{d,\alpha\beta} A_{d,\gamma\delta} \langle \hat{n}_\alpha \hat{n}_\beta \hat{n}_\gamma \hat{n}_\delta \rangle$$
   $$= \frac{1}{15} \sum_{\alpha\beta\gamma\delta} B_{d,\alpha\beta} A_{d,\gamma\delta} \left[ \delta_{\alpha\beta}\delta_{\gamma\delta} + \delta_{\alpha\gamma}\delta_{\beta\delta} + \delta_{\alpha\delta}\delta_{\beta\gamma} \right]$$
   Since $\mathbf{A}_d$ and $\mathbf{B}_d$ are symmetric:
   $$= \frac{1}{15} \left[ \operatorname{Tr}(\mathbf{B}_d)\operatorname{Tr}(\mathbf{A}_d) + 2 \operatorname{Tr}(\mathbf{B}_d \mathbf{A}_d) \right] = \frac{1}{15} \left[ \operatorname{Tr}(\mathbf{B}_d)\operatorname{Tr}(\mathbf{A}_d) + 2 \mathbf{B}_d : \mathbf{A}_d \right]$$

3. **Re-factorization:**
   Factoring out the leading $\frac{1}{3} \operatorname{Tr}(\mathbf{B}_d)$:
   $$I_1(Q, \omega) \approx \frac{1}{3} Q^2 \operatorname{Tr}(\mathbf{B}_d) \left[ 1 - Q^2 \cdot \frac{1}{5} \left( \operatorname{Tr}(\mathbf{A}_d) + 2 \frac{\mathbf{B}_d : \mathbf{A}_d}{\operatorname{Tr}(\mathbf{B}_d)} \right) \right]$$

4. **Re-Exponentiation (Cumulant Resummation):**
   Using the cumulant identity $1 - x \approx e^{-x}$:
   $$I_1(Q, \omega) = \frac{1}{3} Q^2 \operatorname{Tr}(\mathbf{B}_d(\omega)) \exp\left( -Q^2 \alpha_d(\omega) \right)$$
   where:
   $$\alpha_d(\omega) \equiv \frac{1}{5} \left[ \operatorname{Tr}(\mathbf{A}_d) + 2 \frac{\mathbf{B}_d(\omega) : \mathbf{A}_d}{\operatorname{Tr}(\mathbf{B}_d(\omega))} \right]$$
   **Q.E.D.** This proves Eq. (12) of Harrelson et al. (2021) and validates `src/atom.py:63` of `MolDyINS`.

---

### 2.4 Multiphonon Overtone Convolution Series
For higher-order vibrational excitations ($n \ge 2$), multiphonon processes arise from combinations of $n$ vibrational quanta.

In Harrelson's formalism:
$$S_{n, d}(Q, \omega) = \frac{\sigma_d}{4\pi} \frac{Q^{2n}}{3 \cdot n!} \left(\frac{3}{5}\right)^{n-1} \operatorname{Tr}[\mathbf{b}_d^{*n}(\omega)] \exp\left( -\frac{1}{3} Q^2 \operatorname{Tr}(\mathbf{A}_d) \right)$$
where $\mathbf{b}_d^{*n}(\omega)$ is the $n$-fold recursive Cartesian convolution:
$$\mathbf{b}_d^{*n}(\omega) = \int_{-\infty}^{\infty} \mathbf{b}_d^{*(n-1)}(\omega') \mathbf{B}_d(\omega - \omega') d\omega'$$

Computed via 1D Fast Fourier Transform (FFT) along the frequency axis:
$$\mathcal{F}\{\mathbf{b}_d^{*n}\} = \left( \mathcal{F}\{\mathbf{B}_d\} \right)^n$$
$$\mathbf{b}_d^{*n}(\omega) = \mathcal{F}^{-1}\left\{ \left( \mathcal{F}\{\mathbf{B}_d\} \right)^n \right\}$$

---

### 2.5 Detailed Balance & Finite-Temperature Boltzmann Population
In Harrelson et al. (2021), Eqs. (18)–(21) state that transitions from an initial vibrational state $\nu$ to a final state $\nu'$ are thermally weighted by the canonical Boltzmann distribution:
$$P_\nu = \frac{e^{-\beta E_\nu}}{\mathcal{Z}}, \quad \mathcal{Z} = \sum_\nu e^{-\beta E_\nu}$$
The transition rate involves the matrix element $|\langle \nu' | e^{i\mathbf{Q}\cdot\hat{\mathbf{u}}} | \nu \rangle|^2$.

For down-scattering (neutron energy loss, creating phonons, $\omega > 0$):
$$S_{\text{quantum}}(Q, \omega) = S_{\text{sym}}(Q, \omega) \cdot e^{\frac{\hbar \omega}{2 k_B T}}$$
Or in standard asymmetric form:
$$S(Q, -\omega) = e^{-\frac{\hbar \omega}{k_B T}} S(Q, \omega)$$
*Implementation Note:* In `MolDyINS`, the author omitted finite-temperature multi-state convolutions, calculating solely the $0 \to n$ transition. In our clean pipeline, detailed balance must be rigorously enforced via $e^{-\beta \hbar \omega / 2}$ symmetry.

---

### 2.6 Spectrometer Kinematics & Resolution Broadening
1. **Inverted-Geometry Kinematics (TOSCA & VISION):**
   Neutrons scatter off the sample and pass through an analyzer crystal (pyrolytic graphite) that fixes the final energy $E_f$:
   * **TOSCA (ISIS):** $E_f \approx 3.32\text{ meV}$ ($\lambda_f \approx 4.96\text{ \AA}$)
   * **VISION (SNS):** $E_f \approx 3.90\text{ meV}$
   
   The incident energy is $E_i = E_f + \hbar \omega$. The neutron wavevectors are:
   $$k_f = \sqrt{\frac{2 m_n E_f}{\hbar^2}}, \quad k_i = \sqrt{\frac{2 m_n (E_f + \hbar \omega)}{\hbar^2}}$$
   For a scattering angle $2\theta$:
   $$Q(\omega) = \sqrt{k_i^2 + k_f^2 - 2 k_i k_f \cos(2\theta)}$$
   * Forward bank: $2\theta \approx 45^\circ$
   * Backward bank: $2\theta \approx 135^\circ$

2. **Instrument Resolution Kernel:**
   Spectrometers have finite energy resolution that broadens delta-function peaks:
   $$I_{\text{conv}}(\omega) = \int S(Q(\omega'), \omega') R(\omega - \omega', \omega') d\omega'$$
   For TOSCA:
   $$\frac{\Delta E}{E} \approx 0.0125 \implies \sigma(\omega) = \frac{0.0125 \cdot \omega}{2\sqrt{2\ln 2}} + \sigma_0$$
   For VISION:
   $$\frac{\Delta E}{E} \approx 0.0100 \implies \sigma(\omega) = 0.010 \cdot \omega + 1.21\text{ cm}^{-1}$$

---

## 3. Dissected Code Flaws in `MolDyINS` & Their Rectifications

| # | File & Line in `MolDyINS` | Code Flaw / Bug | Physical Impact | Rectification in Clean Pipeline |
|---|---|---|---|---|
| 1 | `MolDyINS.py:182, 185` | `currAtom.Slaw = currAtom.Slaw + np.absolute(pre*q_sq*np.trace(currAtom.fn)*DWF)` | Missing bound neutron cross-section $\sigma_i$. All atoms weighted identically! | Multiply by $\frac{\sigma_i}{4\pi}$ (or $\sigma_i^{\text{inc}} + \sigma_i^{\text{coh}}$) for each atom species. |
| 2 | `src/atom.py:82` | `f1 = 2*10**12*const.hbar/(2*self.freqs*const.Boltzmann*self.T)*self.powspec` | Mass $M_i$ missing in denominator of $f_1$. Dimensions of displacement are wrong. | Divide by $M_i$: $\mathbf{B}_i(\omega) \propto \frac{\hbar}{2 M_i \omega} \operatorname{Re}[\mathbf{\Gamma}_i(\omega)]$. |
| 3 | `MolDyINS.py:156` | `kT_testx = ... * const.m_p ...` | Hardcoded proton mass `m_p` instead of element mass $M_i$. | Use element mass array $M_d$ loaded from trajectory / periodic table. |
| 4 | `src/atom.py:63` | `alpha = 0.2*(np.trace(A)+2*np.true_divide(np.absolute(np.trace(np.dot(A,B))),np.absolute(np.trace(B))))` | Takes `np.absolute` of tensor trace to suppress imaginary parts from FFT noise. | Enforce Hermiticity $\frac{1}{2}(\mathbf{\Gamma} + \mathbf{\Gamma}^\dagger)$ on spectral tensor upfront; traces are strictly real. |
| 5 | `src/atom.py:97-100` | $E_f = 32\text{ meV}$, $\theta = 135^\circ$ hardcoded | Inverted geometry kinematics do not match TOSCA ($3.32\text{ meV}$) or VISION ($3.9\text{ meV}$). | Parameterized instrument class (`InstrumentModel.TOSCA()`, `InstrumentModel.VISION()`). |
| 6 | `MolDyINS.py:181-186` | Only $0 \to n$ ground-state harmonic transitions evaluated | Paper Eqs. (18)–(21) for finite-$T$ Boltzmann convolutions were never coded. | Implement exact detailed balance factor $e^{-\beta \hbar \omega / 2}$ and multi-state convolution. |

---

## 4. Concrete Implementation Steps (Phased Roadmap)

### Phase 1: Core Physics & Intermediate Data Structures
1. **Element & Cross-Section Database:** NIST bound coherent/incoherent neutron scattering lengths and cross-sections ($\sigma_{\text{inc}}$, $\sigma_{\text{coh}}$, $M_d$).
2. **Trajectory Ingestor:** Clean reader for ASE `Atoms` objects, ExtXYZ, NetCDF (`.nc`, `.mdt`), and HDF5, outputting positions $\mathbf{R}(t)$ and velocities $\mathbf{V}(t)$.
3. **Spectral Tensor Calculator:**
   * Welch-windowed FFT yielding strictly Hermitian $\mathbf{\Gamma}_d(\omega) \in \mathbb{C}^{N_{\text{atoms}} \times N_\omega \times 3 \times 3}$.
   * Equipartition sum-rule unit validation: $\frac{1}{2\pi}\int \operatorname{Tr}(\mathbf{\Gamma}_d) d\omega = \frac{3 k_B T}{M_d}$.
   * Displacement tensor $\mathbf{B}_d(\omega)$ with exact $\hbar / (2 M_d \omega)$ scaling.

### Phase 2: Powder Averaging & Multiphonon Engine
1. **Debye-Waller Exponents:**
   * Isotropic: $W_d(Q) = \frac{1}{6} Q^2 \langle u_d^2 \rangle$.
   * Almost-Isotropic: $\alpha_d(\omega) = \frac{1}{5}[\operatorname{Tr}(\mathbf{A}_d) + 2 \frac{\mathbf{B}_d(\omega) : \mathbf{A}_d}{\operatorname{Tr}(\mathbf{B}_d(\omega))}]$.
2. **Multiphonon Series ($n = 1 \dots 10$):**
   * Dual-space FFT polynomial convolution $\mathbf{b}_d^{*n}(\omega) = \mathcal{F}^{-1}\{(\mathcal{F}\{\mathbf{B}_d\})^n\}$.
   * Prefactor weighting $\frac{1}{3 \cdot n!} (\frac{3}{5})^{n-1} Q^{2n}$.

### Phase 3: Instrument Kinematics & Resolution Broadening
1. **Spectrometer Trajectories $Q(\omega)$:**
   * TOSCA ($E_f = 3.32\text{ meV}$, forward bank $45^\circ$, backward bank $135^\circ$).
   * VISION ($E_f = 3.90\text{ meV}$).
2. **Broadening Convolution:**
   * FFT convolution with energy-dependent Gaussian/Lorentzian FWHM $\Delta E(E)$.

### Phase 4: Verification Gate & Unit Test Suite
1. **Harmonic Solid Benchmark:** Compare against Euphonic on an analytical force-constant model.
2. **SrTiO3 ExtXYZ Trajectory Test:** Run against the existing `small-sto-T50-nve.extxyz` in this repository and benchmark against MDANSE DOS.
3. **Cross-Section Weighting Test:** Verify that hydrogen dominates over carbon in P3HT or polymer models according to NIST $\sigma_{\text{inc}}$ tables.

