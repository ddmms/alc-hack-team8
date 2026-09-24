# Comparison of `opus` and `openspec`: Implementation and Physics Review

**Date:** September 2026  
**Subject:** Technical, mathematical, and physical comparison of the `opus` and `openspec` branches in `alc-hack-team8`.

---

## 1. Executive Summary & Verdict

Both branches set out to solve the same problem: converting molecular dynamics (MD) trajectories into atom-projected phonon densities of states (pDOS) and subsequent inelastic neutron scattering (INS) spectra (following Cheng et al. 2020 and Harrelson et al. 2021).

* **The Reality of Implemented Code:** Neither branch has implemented the actual scattering spectra (Stages 2 & 3: isotropic and anisotropic INS intensities). In `opus`, `src/mdins/scattering.py` contains typed `NotImplementedError` stubs; in `openspec`, Stages 2 and 3 exist only as text proposal documents. **Both branches are currently Stage 1 (Trajectory $\to$ pDOS) codebases.**
* **Recommendation:** **`opus` is strongly recommended.** Its Stage 1 implementation is vastly superior in physical rigor, statistical mechanics foundation, streaming memory management, tensor mathematics, persistence, and verification.

---

## 2. Actual Code Implementation: Reality vs. Promises

| Capability / Stage | `openspec` (`src/md_ins`) | `opus` (`src/mdins`) |
| :--- | :--- | :--- |
| **Stage 1: Trajectory Ingestion** | Implemented. Ingests via ASE. Derives velocities from positions if missing. | Implemented. Ingests via ASE. Strictly enforces explicit velocities. Removes COM translation and angular drift. |
| **Stage 1: Spectral Estimation** | Implemented (~250 lines). Scalar VACF via Wiener–Khinchin. Unit-normalized. | Implemented (~1,500 lines). Full $3 \times 3$ tensor. Dual Welch / VACF estimators. Guaranteed PSD. |
| **Stage 1: Memory & Streaming** | In-memory only. Allocates full FFT array across all atoms at once. | Actively chunks atoms into 64 MB blocks. Direct integral-conserving rebinning. |
| **Stage 1: Persistence & CLI** | None (ephemeral `PDOSResult`). No CLI. | HDF5 serialization with 6-component Voigt packing and provenance. Working `mdins pdos` CLI. |
| **Stages 2 & 3: Scattering Spectra** | **0 lines of Python** (proposals in `openspec/changes/`). | **0 lines of Python** (typed `NotImplementedError` stubs in `scattering.py`). |

---

## 3. Physical Analysis: The Underlying Physics Compared

### 3.1 Normalisation and the Equipartition Theorem
* **`openspec`**: Normalises every atom's pDOS to unit integral:
  $$\int_0^\infty g_i(E)\, dE = 1.0$$
  This discards the absolute velocity amplitude $\langle |v_i|^2 \rangle$. Consequently, downstream scattering intensities cannot recover absolute cross sections without ad-hoc scaling, and all species contribute equal weight to total DOS regardless of mass.
* **`opus`**: Stores the absolute velocity cross-spectral density tensor $P_{i,\alpha\beta}(E)$ in units of $\text{Å}^2\cdot\text{ps}^{-2}\cdot\text{meV}^{-1}$, adhering directly to classical equipartition:
  $$\int_0^\infty \frac{1}{3}\mathrm{tr}\, P_i(E)\, dE = \frac{k_B T}{m_i}$$
  The physical scale required for displacement tensors $\mathbf{B}_i$ and Debye–Waller factors is preserved, and lighter atoms naturally carry larger velocity spectral densities ($\propto 1/m_i$).

---

### 3.2 Statistical Mechanics of NVE Trajectories: The "Frozen Draw"
Both pipelines validate MD pDOS against harmonic phonon references (Phonopy/Euphonic) in low-temperature NVE argon.
* **`openspec`**: Assumes that running an NVE trajectory longer will converge the pDOS onto the true harmonic answer.
* **`opus`** (`docs/interactive-validation.ipynb` §3): Recognizes that in an NVE simulation of a near-harmonic crystal, **the energy of every normal mode is a constant of motion**. The initial velocity draw from the Maxwell–Boltzmann distribution freezes the energy partition across normal modes:
  * Running the MD longer measures the *same frozen microstate* with higher precision; it does **not** converge to the ensemble harmonic DOS.
  * In a 27-atom cell (78 modes), mode energy fluctuations are $\sim 1/\sqrt{78} \approx 11\%$, and seed-to-seed scatter is $\sim 4.5\times$ larger than the single-run Welch uncertainty.
  * `opus` models this by testing across multi-seed ensembles rather than mistaking single-run precision for physical accuracy.

---

### 3.3 Deriving Velocities from Positions: The $\operatorname{sinc}$ Transfer Function
When input trajectories lack explicit velocities:
* **`openspec`**: Automatically computes central differences:
  $$v(t) \approx \frac{r(t + \Delta t) - r(t - \Delta t)}{2\Delta t}$$
  **Physics Failure**: Central differencing applies an uncorrected frequency-dependent transfer function:
  $$\hat{v}(\omega) = i\omega \operatorname{sinc}(\omega \Delta t)\, \hat{r}(\omega)$$
  Near the Nyquist boundary, $\operatorname{sinc}(\omega \Delta t)$ attenuates high-frequency modes by up to 36%. Furthermore, periodic boundary condition (PBC) box jumps create unphysical delta-function velocity spikes.
* **`opus`**: Refuses to derive velocities from positions, requiring true velocities from the integrator.

---

### 3.4 Finite-Size Center-of-Mass Corrections
Removing center-of-mass velocity drift changes the single-atom velocity expectation in a finite box of mass $M$:
* **`openspec`**: Ignores COM correction.
* **`opus`**: Derives and tests the exact analytical finite-cell identity:
  $$\langle |\mathbf{v}_i - \mathbf{V}_{\text{COM}}|^2 \rangle = 3 k_B T \left(\frac{1}{m_i} - \frac{1}{M}\right)$$
  For equal-mass atoms, this introduces an exact $(1 - 1/N)$ correction; for light atoms in small cells (e.g., carbon in methane), it shifts the expectation by 25%. `opus` validates its normalisation against this exact finite-size bound.

---

### 3.5 Real MD Trajectories, Dump Intervals, and Nyquist Aliasing
In `docs/interactive-validation.ipynb` §9, a real MACE-MP machine-learning MD trajectory of solid benzene was evaluated:
* Standard MD trajectories often dump frames every $100\text{ fs}$ ($\Delta t = 0.1\text{ ps}$), giving a Nyquist energy limit of $E_{\text{Nyquist}} = h/(2\Delta t) \approx 20.7\text{ meV}$ ($167\text{ cm}^{-1}$).
* Benzene's $\text{C-H}$ stretches sit near $3000\text{ cm}^{-1}$ ($370\text{ meV}$), which is $18\times$ higher than the Nyquist cutoff.
* In any Fourier transform, mode energy above $E_{\text{Nyquist}}$ folds back (aliases) into the low-frequency region.
* **`opus`** uses the equipartition sum rule as a physical detector: total kinetic power is conserved under aliasing, exposing that hydrogen and carbon were at different effective temperatures ($23.7\text{ K}$ vs $15.8\text{ K}$) and that 83% of the degrees of freedom were aliased. `openspec`'s unit-normalized pDOS disguises this aliasing completely.

---

### 3.6 Downstream INS Scattering: Series Truncation vs. Gaussian Approximation
Looking ahead to Stages 2 and 3 (Cheng et al. and Harrelson et al.):
* **`openspec`**: Proposes evaluating iterative self-convolutions with combinatorial prefactors ($1/n!$ and $3^{n-2}/(n! 5^{n-1})$) truncated at $N_{\max} = 4$.
* **`opus`** (`docs/method-review.md` §2): Demonstrates that the infinite convolution series resums to the Gaussian-approximation incoherent intermediate scattering function:
  $$I_i(\mathbf{q}, t) = \sigma_i \exp\left(-\frac{1}{2}\langle [\mathbf{q}\cdot \mathbf{u}_i(0)]^2\rangle\right) \exp\left(\langle [\mathbf{q}\cdot \mathbf{u}_i(0)][\mathbf{q}\cdot \mathbf{u}_i(t)]\rangle\right)$$
  Direct exponentiation in the time domain yields the Debye–Waller factor, all multiphonon orders, and exact detailed balance without truncation or ad-hoc prefactors.

---

## 4. Software Engineering and Usability

| Dimension | `openspec` | `opus` |
| :--- | :--- | :--- |
| **Tensor Storage** | Scalar only (dot product). | Full $3 \times 3$ cross-spectral tensor; Voigt 6-component disk packing. |
| **Spectral Estimator** | Blackman–Tukey VACF only (can produce non-PSD tensors). | Welch (guaranteed positive semi-definite) + Blackman–Tukey VACF. |
| **Memory Control** | Unbounded (full FFT arrays). | Chunks atoms in 64 MB blocks; integral-conserving rebinning. |
| **File Format** | None. | Structured HDF5 with complete metadata and provenance. |
| **Command Line** | None. | `mdins pdos <traj> --dt <dt> -o <out.h5>`. |
| **Test Suite** | ~415 lines; heuristic tolerances. | ~3,500 lines; analytical models (harmonic oscillator, Langevin Brownian motion). |
| **Code Hygiene** | No type checking, no CI. | Strict `mypy` typing, `ruff` linting, GitHub Actions CI. |

---

## 5. Recommended Action Plan

1. **Adopt `opus` as the codebase**: It has the correct physical foundation, rigorous math, robust data contracts, and complete Stage 1 execution.
2. **Borrow two ergonomic features from `openspec`**:
   * Add a nominal kinetic temperature sanity check ($|T_{\text{traj}} - T_{\text{MD}}| / T_{\text{MD}} > 0.2$) with a descriptive warning.
   * Add an in-memory trajectory ingestion helper (`from_arrays(velocities, masses, symbols, dt)`) to streamline Jupyter notebook prototyping.
3. **Proceed with Stages 2 & 3 in `opus`**: Implement `isotropic_spectrum` and `anisotropic_spectrum` in `src/mdins/scattering.py` using the Gaussian-approximation time-domain route proven in `docs/method-review.md`.

