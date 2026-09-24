# Simulating INS spectra from MD trajectories — method review

Review of the two target papers, the literature they depend on, and the information
that is missing if you want to reimplement the methods generally (MD-engine agnostic,
instrument agnostic, open source).

## 1. The two papers

### Paper A — Cheng, Kolesnikov & Ramirez-Cuesta (2020)

*Simulation of Inelastic Neutron Scattering Spectra Directly from Molecular Dynamics
Trajectories*, J. Chem. Theory Comput. **16**, 7702–7708.
[doi:10.1021/acs.jctc.0c00937](https://doi.org/10.1021/acs.jctc.0c00937).
Closed source (OCLIMAX, distributed as a Docker image).

Free accepted manuscript: <https://www.osti.gov/servlets/purl/1845806> (OSTI 1845806).

Method, in full:

1. Normalised velocity autocorrelation per atom,
   `C_d(t) = <v_d(t+t0)·v_d(t0)> / <v_d(t0)·v_d(t0)>`   (Eq. 5)
2. Partial (atom-projected) phonon DOS, `g_d(ω) = (1/2π) ∫ exp(-2πiωt) C_d(t) dt`  (Eq. 6)
3. Bridge to lattice dynamics: `∫_ω^{ω+dω} g_d dω = Σ_s tr(B_ds)/3`  (Eq. 7),
   where `B_ds = e_ds e_ds^T` is the displacement tensor and `A_d = Σ_s B_ds (2n_s+1)`.
4. Incoherent one-phonon scattering
   `S±1(Q,ω) = Σ_d (3σ_d / 2M_d) Q² exp(-2W_d) [g_d(ω)/ω] (n + ½ ± ½)`  (Eq. 8),
   with `W_d = Q² tr(A_d) / 6`  (Eq. 4).
5. Higher orders by repeated convolution, `S_i(ω) = ∫ S_1(ω-ω') S_{i-1}(ω') dω'`  (Eq. 9),
   performed **at atomic level, before** applying cross-section and Debye–Waller factor.
6. Bose–Einstein statistics throughout, so the MD temperature is decoupled from the
   experimental temperature — it becomes a *fictive* parameter chosen only to survey the
   right region of the PES (they run ice and ZIF-8 at 180 K to mimic proton zero-point
   motion classically), while population and DWF come from quantum statistics at the real T.

Everything else — powder averaging, multiphonon prefactors, instrument Q–ω trajectories,
resolution functions, coherent scattering — is inherited from the earlier OCLIMAX paper.
**Paper A on its own is about six equations; the real specification is Paper A′.**

### Paper A′ — Cheng, Daemen, Kolesnikov & Ramirez-Cuesta (2019)

*Simulation of Inelastic Neutron Scattering Spectra Using OCLIMAX*, JCTC **15**, 1974–1982.
[doi:10.1021/acs.jctc.8b01250](https://doi.org/10.1021/acs.jctc.8b01250).
Free accepted manuscript: <https://www.osti.gov/servlets/purl/1845813> (OSTI 1845813).

Contains what Paper A omits: coherent `S(Q,ω)` (Eq. 1), incoherent form (Eq. 2), low-T
overtone expansion (Eq. 3), the analytic "almost isotropic" powder average (Eq. 4), the
exact two-quantum expression (Eq. 5), phonon wings (Eqs. 9–10), and resolution functions
(Eqs. 11–12).

### Paper B — Harrelson, Dettmann, Scherer, Andrienko, Moulé & Faller (2021)

*Computing inelastic neutron scattering spectra from molecular dynamics trajectories*,
Sci. Rep. **11**, 7938. [doi:10.1038/s41598-021-86771-5](https://doi.org/10.1038/s41598-021-86771-5).
Open access; full text and SI on PMC (PMC8041884).

**The code is public**: <https://github.com/tfharrelson/MolDyINS> (GPL-3.0). See §5.

Same skeleton, three differences:

- Keeps the **full 3×3 tensor** `<v_α*(ω) v_β(ω)>` (Eq. 16), not just its trace, so
  vibrational anisotropy survives into the powder average.
- Eliminates the mode reduced mass via equipartition:
  `u_i²(ω) = ħ v_i²(ω) / (ω k_B T)`  (Eq. 8).
- Overtones by repeated self-convolution of `B` (Eq. 17) with an assumed `1/n!` prefactor;
  0 K base spectrum, with an optional convolution against a Boltzmann-weighted DOS to go
  to finite T (Eqs. 18–20).

## 2. The most useful thing in either paper

Supplementary §2 of Paper B (Eqs. S17–S24) proves that the infinite series of convolutions
sums to an exponential:

```
I_i(q,t) = σ_i · exp(-∫F_i(ω)dω) · exp(f_i(t))
F_i(ω)   = ħ |q·v_i(ω)|² / (2 ω k_B T)
f_i(t)   = <[q·r_i(0)][q·r_i(t)]>
```

which is just the standard **Gaussian-approximation incoherent intermediate scattering
function**, `I_i(q,t) = σ_i <exp(-iq·r_i(0)) exp(iq·r_i(t))>`.

Implication for implementation: do not implement an N-term convolution loop with
hand-guessed prefactors. Exponentiate once. You get the Debye–Waller factor, every
multiphonon order, and the correct `1/n!` combinatorics for free and self-consistently,
and the resulting object is the same one MDANSE and nMOLDYN already compute — which
gives you an independent cross-check. The convolution formulation is then an optimisation
/ order-decomposition view, not the definition.

## 3. Missing information — Paper A (OCLIMAX)

1. **Absolute normalisation is undetermined.** Eq. 5 defines a *normalised* VACF (divided
   by `v(0)·v(0)`), so `∫g_d dω` is fixed by that normalisation and is dimensionless.
   Eq. 7 equates it to `Σ_s tr(B_ds)/3`, a squared displacement. These cannot both hold as
   written. The intensity scale has to be recovered from a sum rule (e.g.
   `Σ_s tr(B_ds)(2n_s+1) = 3<u_d²>`) — the paper never states which convention is used.
2. **Anisotropy is discarded and the consequence is unstated.** Eq. 7 delivers only
   `tr(B_ds)`, but the almost-isotropic powder average (A′ Eq. 4) needs `B_s:A` and
   `tr(B_s)` separately. Paper A says the reduced information "is sufficient", implying a
   silent fall back to the fully isotropic DWF for MD input. Not confirmed anywhere.
3. **The `tclimax` trajectory format is undocumented**, and the format converter ships only
   inside the Docker image. This is the concrete blocker: you must define your own
   intermediate representation.
4. **No VACF numerics at all.** Nothing on window/apodisation function, correlation-window
   length vs trajectory length, zero padding, direct FT of velocities vs Wiener–Khinchin,
   block averaging, centre-of-mass drift removal, or how the `g(ω)/ω` divergence as ω→0
   (and quasi-elastic/diffusive contamination) is regularised.
5. **Detailed balance after convolution.** Eq. 2 carries `∓`, but how both energy-gain and
   energy-loss sides are generated and how detailed balance is preserved through the
   multiphonon convolution is not described.
6. **Quantum correction factor is implicit.** Classical MD amplitudes already carry
   classical occupancy; Eq. 8 then multiplies by `(n+½±½)`. This is only consistent if `g`
   is a pure normalised lineshape — which loops back to (1). Neither paper names the QCF
   choice it is making; see Ramírez et al. in §6.
7. **Resolution constants disagree with the open implementation.** A′ Eq. 12 gives
   `σ = 0.25 + 0.005·ΔE + 1e-7·ΔE²` (wavenumbers) for TOSCA/VISION. Mantid AbINS uses the
   same `a` and `b` but `c = 2.5`, a factor of ten larger. One is a typo; check with an
   instrument scientist before trusting either. Also, VISION ≠ TOSCA and VISION's
   parameters are not published anywhere we found.
8. **Coherent scattering from MD does not exist in either paper.** OCLIMAX does coherent
   `S(Q,ω)`, but only from lattice dynamics. If you want coherent-from-MD you are past the
   literature and into `F(Q,t)` territory (see dynasor, §5).

## 4. Missing information — Paper B (MolDyINS)

1. **Eq. 8 is dimensionally inconsistent** and the paper says so ("the correct units are
   recovered with an integral over ω"). The actual bookkeeping is only in SI Eq. S13, in
   terms of discrete Fourier coefficients. Follow S13, not Eq. 8.
2. **The `1/n!` overtone prefactor is asserted, not derived** — "it is a general relation
   that closely resembles the known prefactors for the first four overtones". The
   exponential form in §2 above removes the need for it entirely.
3. **The rank-4 → rank-2 contraction in Eq. 17 is an uncontrolled approximation**
   ("approximated as a simple matrix multiplication"), with no error estimate.
4. **Eq. 19 looks wrong as printed.** `ρ(ω;T) = Σ_a exp(+ħω_a/k_B T) δ(ω+ω_a)` has a
   positive exponent and no partition function, so it grows without bound with mode energy.
   Almost certainly should be `exp(-ħω_a/k_BT)/Z`. Verify before implementing.
5. **DWF Eq. 21 is a heuristic**: `∫ max[(q·u_i(ω))², (q·ω v_i(ω))²] dω`, taking whichever
   of the ground-state or thermal displacement is larger at each frequency. No derivation,
   no continuity guarantee at the crossover.
6. **The anharmonic correction (Eqs. 22–23, SI §3) was never used.** It is derived, then
   dropped — and the authors note the classical decay rates correspond to the thermalised
   initial state, not the vibrational ground state the theory needs.
7. **Thermostat constraints are load-bearing but buried in SI §4.** Production runs need
   NVE or a weak decorrelating thermostat (Andersen-massive, τ ≈ 10 ps). τ = 1 ps corrupts
   the spectrum; barostats damp the dynamics. Any general implementation should validate
   and warn about this.
8. **Instrument handling is hardcoded to VISION.** `E_f = 32 cm⁻¹` is assumed; the `Q(ω)`
   relation is said to be "determined by the instrument setup" but never written down; the
   resolution `σ = 0.01·E` is admitted to be chosen by eye to match the C–H stretch width.

### Sampling requirements (from Paper B, generalise these)

2 fs velocity dump over 100 ps → Nyquist ≈ 8333 cm⁻¹, frequency spacing 0.33 cm⁻¹. In
general: `Δt_dump < 1/(2 c ν_max)` and `T_total > 1/(c Δν)`. Aliasing is not filtered —
modes above Nyquist fold back and contaminate the spectrum, so the dump interval is a
correctness requirement, not a performance knob.

## 5. State of the released code

`tfharrelson/MolDyINS`, GPL-3.0, last commit June 2021, ~2 stars. Six Python files
(~40 kB), plus a 13 MB `.trr` test trajectory at `data/25K_NVE/`.

It is GROMACS-only by construction: it shells out to `gmx`/`gmx_mpi` and reads `.trr`,
`.tpr` and `.ndx`. It also uses `np.float`, removed in NumPy 1.24, so it will not run on a
current stack without edits, and it relies on module-level globals for configuration.

Treat it as a reference and a source of regression fixtures — the bundled 25 K trajectory
is a ready-made test case — rather than a base to build on.

OCLIMAX itself is Docker-only and binary, so it cannot be a base either, but it is usable
as a numerical oracle for validation.

## 6. Literature you need

### Theory

| Work | Why |
|---|---|
| Squires, *Introduction to the Theory of Thermal Neutron Scattering* (Dover, 1996) | Definitions of `S(Q,ω)`; cited as the source equation by both papers |
| Mitchell, Parker, Ramirez-Cuesta & Tomkinson, *Vibrational Spectroscopy with Neutrons* (World Scientific, 2005) | **The** INS-simulation reference. Source of the almost-isotropic powder average and phonon wings. Paper B's ref 11, cited 8×, and never reproduced |
| Ramirez-Cuesta, *aCLIMAX 4.0.1*, Comput. Phys. Commun. **157**, 226 (2004), [doi:10.1016/S0010-4655(03)00520-4](https://doi.org/10.1016/S0010-4655(03)00520-4) | The powder-average and multiphonon derivations both papers cite but omit |
| Cheng et al., *OCLIMAX*, JCTC **15**, 1974 (2019) | Actual spec for Paper A (see §1) |
| Cheng & Ramirez-Cuesta, JCTC **16**, 5212 (2020), doi:10.1021/acs.jctc.0c00569 | Thermal neutron cross-sections; free at OSTI 1845808 |
| Maradudin & Fein, Phys. Rev. **128**, 2589 (1962) | Anharmonic crystal scattering; the detailed-balance result Paper B Eq. 20 reduces to |
| Fultz et al., *Experimental INS with a Chopper Spectrometer and Virtual Neutron Scattering with a Computer* (2018) | Origin of the multiphonon convolution both codes use |
| Ramírez, López-Ciudad, Kumar & Marx, J. Chem. Phys. **121**, 3973 (2004), [doi:10.1063/1.1774986](https://doi.org/10.1063/1.1774986) | Quantum correction factors for classical time-correlation functions. This is the choice both papers make implicitly and differently — read before picking one |

### Open-source implementations

| Tool | Licence | What it gives you |
|---|---|---|
| **Mantid AbINS / Abins2** (`scripts/abins/`) | GPL-3.0 | Open implementation of exactly the aCLIMAX/OCLIMAX maths: `powdercalculator.py`, `frequencypowdergenerator.py`, `spowdersemiempiricalcalculator.py`, plus `instruments/` with calibrated TOSCA, Lagrange, PANTHER and MAPS/MARI/MERLIN (via PyChop). Documented numerical parameters (`autoconvolution`: `max_order=10`, `min_order=3`, `fine_bin_factor=10`; `s_relative_threshold=0.01`). **Its loaders are all ab initio phonon codes — there is no MD trajectory input.** That gap is precisely the contribution being proposed |
| **MDANSE** ([ISISNeutronMuon/MDANSE](https://github.com/ISISNeutronMuon/MDANSE), active) | GPL-3.0 | Multi-format trajectory reader, VACF/DOS, incoherent and coherent dynamic structure factors under the Gaussian approximation. Closest existing general engine; lacks the quantum multiphonon and instrument layer. Goret, Aoun & Pellegrini, JCIM **57**, 1 (2017) |
| **nMOLDYN 3** | | Predecessor; the original Gaussian-approximation `I(q,t)` implementation. Hinsen et al., JCC **33**, 2043 (2012) |
| **Euphonic** ([pace-neutrons/Euphonic](https://github.com/pace-neutrons/Euphonic), active) | GPL-3.0 | Coherent/incoherent `S(Q,ω)` from force constants; already wired in as an AbINS loader |
| **dynasor** ([gitlab materials-modeling](https://gitlab.com/materials-modeling/dynasor), v2.5) | | Coherent and incoherent dynamic structure factors straight from MD. The natural starting point if you want the coherent branch neither paper covers |
| **LiquidLib** | | Walter et al., CPC **228**, 209 (2018) |
| **OCLIMAX** | binary/Docker | Not a base, but a validation oracle |

### Validation systems (all with published experimental spectra)

- Ice Ih at VISION and SEQUOIA, ZIF-8, silica glass — Paper A. Built with GenIce, MB-pol/MBX + i-PI, CHIK/LAMMPS.
- P3HT crystalline and amorphous — Paper B, plus the 25 K trajectory in the MolDyINS repo.
- Toluene, MgH₂ (strong multiphonon), LiAlH₄ (temperature series), graphite (coherent) — Paper A′.

## 7. Suggested shape of a general implementation

```
trajectory (any engine, via ASE / MDAnalysis)
   → atom-resolved velocity spectrum  v_i,α(ω)   [full 3×3 tensor, Paper B]
   → I_i(Q,t) via the Gaussian-approximation exponential   [SI §2 — not a convolution loop]
   → S_inc(Q,ω), all orders, DWF consistent by construction
   → powder average (almost-isotropic; isotropic as a fallback/option)
   → instrument Q(ω) trajectory + resolution broadening   [reuse AbINS]
```

Two decisions worth settling early, because they change the architecture:

- **Coherent scattering?** Absent from both papers. Including it means carrying
  `F(Q,t) = Σ_{i,j}` rather than per-atom terms — a different cost class, and where dynasor
  would come in.
- **Fictive-temperature decoupling (Paper A) or not (Paper B)?** Paper A's separation of MD
  temperature from the quantum statistics is genuinely useful for nuclear quantum effects,
  but it only works if the normalisation in §3(1) is pinned down first.
