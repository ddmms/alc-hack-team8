# Brainstorm Brief: Trajectory-to-INS via Intermediate Dynamical Representations

**Date:** 2026-09-23  
**Grounding:** Cheng et al. (*JCTC* 2020), Harrelson et al. (*Sci. Rep.* 2021 & SI).  
**Premise:** Trajectory data is assumed to be provided (e.g. ExtXYZ / NetCDF / HDF5). MLIP training/inversion is excluded; focus is 100% on rigorous physics and mathematical implementation of intermediate representations bridging trajectories to INS spectra.

---

## 1. Core Physics & Equation Synthesis from the Two Papers

### A. Trajectory Ingestion & Intermediate Dynamic Representation
Given atomic coordinates and velocities $\mathbf{v}_d(t)$ for each atom $d$:

1. **Cheng et al. (2020) — Isotropic Scalar Representation:**
   * Atom-projected velocity autocorrelation function (VACF):
     $$\gamma_d(t) = \frac{\langle \mathbf{v}_d(0) \cdot \mathbf{v}_d(t) \rangle}{\langle |\mathbf{v}_d(0)|^2 \rangle}$$
   * Partial vibrational density of states (pDOS):
     $$g_d(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \gamma_d(t) e^{-i\omega t} dt$$
   * Mean squared displacement and scalar Debye-Waller exponent:
     $$\langle u_d^2 \rangle = \frac{3 k_B T}{M_d} \int_0^\infty \frac{g_d(\omega)}{\omega^2} d\omega, \quad W_d(Q) = \frac{1}{6} Q^2 \langle u_d^2 \rangle$$

2. **Harrelson et al. (2021) — Anisotropic Cartesian Tensor Representation:**
   * Velocity cross-spectral density tensor $\mathbf{\Gamma}_d(\omega) \in \mathbb{C}^{3 \times 3}$:
     $$\mathbf{\Gamma}_{d, \alpha\beta}(\omega) = \frac{1}{2\pi} \int_{-\infty}^{\infty} \langle v_{d,\alpha}(0) v_{d,\beta}(t) \rangle e^{-i\omega t} dt$$
     Guaranteed Hermitian and positive semi-definite: $\mathbf{\Gamma}_d(\omega) = \mathbf{\Gamma}_d(\omega)^\dagger$, $\text{Tr}(\mathbf{\Gamma}_d(\omega)) \ge 0$.
   * Frequency-dependent displacement tensor $\mathbf{B}_d(\omega) \in \mathbb{R}^{3 \times 3}$:
     $$\mathbf{B}_d(\omega) = \frac{\hbar}{2 M_d \omega} \coth\left(\frac{\hbar \omega}{2 k_B T}\right) \text{Re}[\mathbf{\Gamma}_d(\omega)]$$
   * Total atomic displacement tensor $\mathbf{A}_d$:
     $$\mathbf{A}_d = \int_0^\infty \mathbf{B}_d(\omega) d\omega = \langle \mathbf{u}_d \mathbf{u}_d^T \rangle$$

---

## 2. Quantum Corrections & Powder Averaging

1. **Detailed Balance & Zero-Point Factor:**
   * Classical MD trajectories obey classical equipartition $\langle v^2 \rangle = \frac{3 k_B T}{M}$. To match quantum neutron scattering, the spectrum is scaled by the harmonic quantum factor:
     $$Q(\omega, T) = \frac{\hbar \omega}{2 k_B T} \coth\left(\frac{\hbar \omega}{2 k_B T}\right)$$
   * Detailed balance is enforced by multiplying down-scattering (neutron energy loss) by $\exp\left(-\frac{\hbar \omega}{2 k_B T}\right)$ or the standard asymmetric factor $1 / [1 - \exp(-\beta \hbar \omega)]$.

2. **Directional Powder Averaging on the Unit Sphere:**
   * **Cheng et al. (2020):** Assumes complete isotropic orientation:
     $$\langle (\mathbf{Q} \cdot \mathbf{e})^2 \exp(-(\mathbf{Q} \cdot \mathbf{u})^2) \rangle \approx \frac{1}{3} Q^2 \exp\left(-\frac{1}{3} Q^2 \langle u^2 \rangle\right)$$
   * **Harrelson et al. (2021) — "Almost-Isotropic" Averaging:**
     For an anisotropic displacement tensor $\mathbf{A}_d$, the spherical average of $(\mathbf{Q}^T \mathbf{B}_d \mathbf{Q}) \exp(-\mathbf{Q}^T \mathbf{A}_d \mathbf{Q})$ contracts to:
     $$\alpha_d(\omega) = \frac{1}{5} \left[ \text{Tr}(\mathbf{A}_d) + 2 \frac{\mathbf{B}_d(\omega) : \mathbf{A}_d}{\text{Tr}(\mathbf{B}_d(\omega))} \right]$$
     $$S_{1, d}(Q, \omega) = \frac{\sigma_d}{4\pi} \frac{Q^2}{3} \text{Tr}(\mathbf{B}_d(\omega)) \exp\left( -Q^2 \alpha_d(\omega) \right)$$
     *(Correcting the omission in MolDyINS where cross-section $\sigma_d$ was never multiplied).*

---

## 3. Multiphonon Expansion via Fast Fourier Convolution

For higher-order overtones ($n \ge 2$):
1. **Cheng et al. (2020) Recursive Convolution:**
   $$S_{n, d}(Q, \omega) = \frac{1}{n} \left[ S_{1, d}(Q, \omega) \circledast_\omega S_{n-1, d}(Q, \omega) \right]$$
   Computed via 1D Fast Fourier Transform (FFT) in frequency space:
   $$\mathcal{F}\{S_{n, d}\} = \frac{1}{n!} \left( \mathcal{F}\{S_{1, d}\} \right)^n$$
   Summed up to convergence (typically $n = 1 \dots 10$).

---

## 4. Instrument Kinematics & Resolution Broadening

1. **Spectrometer Kinematics Trajectories $Q(\omega)$:**
   * **Inverted-Geometry (TOSCA at ISIS, VISION at SNS):**
     Fixed final neutron energy $E_f$ ($E_f \approx 3.32\text{ meV}$ for TOSCA, $\sim 3.9\text{ meV}$ for VISION):
     $$Q(\omega) = \sqrt{\frac{2 m_n}{\hbar^2}} \sqrt{2 E_f + \hbar \omega - 2 \sqrt{E_f (E_f + \hbar \omega)} \cos(2\theta)}$$
     where $2\theta$ is the scattering angle (forward bank $\approx 45^\circ$, backward bank $\approx 135^\circ$).
   * **Direct-Geometry (MAPS, MARI, SEQUOIA):**
     Fixed incident neutron energy $E_i$, measuring $S(Q, \omega)$ across a 2D kinematic envelope $Q_{\text{min}}(\omega) \le Q \le Q_{\text{max}}(\omega)$.

2. **Resolution Convolution:**
   * Instrumental Gaussian/Lorentzian broadening with energy-dependent FWHM:
     $$\frac{\Delta E}{E} \approx 1.25\% \text{ (TOSCA)}$$
     $$I_{\text{inst}}(\omega) = \int S(Q(\omega'), \omega') R(\omega - \omega', \omega') d\omega'$$

---

## 5. Architectural Blueprint for Implementation

```
md_ins/
├── io/
│   ├── trajectory.py        # Reads ExtXYZ, NetCDF, HDF5, ASE Atoms list
│   └── cross_sections.py    # NIST neutron scattering cross-sections (sigma_inc, sigma_coh, M)
├── intermediate/
│   ├── vacf.py              # Windowed FFT VACF and isotropic pDOS g_d(w) (Cheng 2020)
│   ├── tensors.py           # 3x3 Cartesian velocity spectral tensors Gamma_d(w) & B_d(w) (Harrelson 2021)
│   └── debye_waller.py      # Isotropic <u^2> and anisotropic displacement tensors A_d
├── scattering/
│   ├── detailed_balance.py  # Quantum Bose/detailed-balance factors
│   ├── powder_average.py    # Isotropic (Cheng) & Almost-Isotropic (Harrelson) contractions
│   └── multiphonon.py       # FFT-accelerated multiphonon overtone convolutions (n=1..10)
├── instruments/
│   ├── kinematics.py        # Q(w) relations for TOSCA, VISION, and direct-geometry spectrometers
│   └── resolution.py        # Instrument-specific resolution broadening kernels
└── pipeline.py              # End-to-end driver: trajectory -> intermediate -> S(Q,w) -> spectrum
```

