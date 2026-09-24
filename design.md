# Design: INS spectra from MD trajectories

Technical blueprint and architecture decisions. Companion to [`proposal.md`](proposal.md)
(why/what) and [`docs/method-review.md`](docs/method-review.md) (source-paper analysis).

**Status: settled.** All eight decisions in §6 are resolved. What remains open are the
physics and numerics ambiguities in §5, which have working resolutions and are settled
during implementation rather than up front.

---

## 1. Pipeline

Seven stages with explicit data contracts between them. The stage boundaries are the
architecture; everything else is an implementation detail we can change later.

```
[A] Trajectory source     velocities v_i(t), masses, species, cell, dt, T_MD
        ↓
[B] Spectral estimator    per-atom velocity cross-spectral density  P_i,αβ(ω)
        ↓
[C] Intermediate repr.    numpy-backed dataclass, HDF5 on request       ← the IR
        ↓
[D] Displacement model    quantum amplitudes u²_i,αβ(ω), Debye–Waller factor A_i
        ↓
[E] Scattering model      S_inc(Q,ω), resolved by excitation order
        ↓
[F] Instrument model      Q(ω) trajectory, resolution broadening, binning
        ↓
[G] Spectrum              output container
```

### The fork point is B→C, and there is only one IR

Method 1 (Cheng) needs only the scalar atom-projected DOS; Method 2 (Harrelson) needs the
full cross-correlation tensor. These are not two pipelines — the tensor IR contains the
scalar as its trace:

```
g_i(ω) = tr P_i(ω) / 3
```

So we compute and store the tensor once, and Method 1 is a consumer that takes the trace
at stage D. This is the single most important structural decision in the design: it means
the two methods differ only in stages D and E, share A–C and F–G, and can be compared on
byte-identical input.

### Stage independence

Stages A–C are useful without the rest — atom-projected DOS from MD has value beyond INS,
and is the first deliverable. Stages D–G consume the IR and never touch a trajectory.
Anyone with phonon data from another source can enter at C.

---

## 2. Stage detail

### [A] Trajectory source — ASE only (D2)

Input is **velocities**, not positions. Combined with the ASE-only decision this is the
tightest constraint in the design, and it should be stated plainly in user-facing docs
rather than discovered:

| Source | Velocities via ASE |
|---|---|
| ASE `.traj` | yes |
| extxyz (with a velocity/momenta array) | yes |
| LAMMPS dump, if velocities were requested | yes |
| GROMACS `.trr` | **no native reader** |
| GROMACS `.xtc`, DCD | not stored in the format at all |

Users of other engines convert to extxyz or `.traj` first. That is a real usability cost
and should be acknowledged, not hidden — see the consequence noted under D2 in §6.

Deriving velocities by finite-differencing positions is possible but introduces a
frequency-dependent transfer function (a `sinc` factor) that must be divided out, and
requires periodic-image unwrapping. It is a fallback, not a default, and if implemented
must be loudly flagged in provenance.

Pre-processing, in order:
1. Remove centre-of-mass velocity per frame (drift produces a spurious ω→0 feature).
2. Optionally remove global angular velocity — only meaningful for non-periodic systems.
3. Validate sampling: assert `dt_dump < 1/(2 c ν_max)` for the requested `ν_max`, and
   `t_total > 1/(c Δν)` for the requested resolution. **Aliasing is silent and
   unrecoverable**, so this is an error, not a warning.

### [B] Spectral estimator — Welch default, VACF as second estimator (D5)

Both routes to `P_i,αβ(ω) = ⟨v*_i,α(ω) v_i,β(ω)⟩` are implemented; they agree wherever both
are valid, and requiring that agreement is an independent check on the stage most likely to
be quietly wrong.

- **Welch-averaged cross-periodogram** (default). Segment the trajectory, window, FFT each
  segment, average the cross-products. Well-defined normalisation
  (`scipy.signal.csd(..., scaling='density')`), variance `∝ 1/K` controlled explicitly by
  segment count, streams naturally, and gives real and imaginary parts separately.
- **VACF → lag window → FFT** (Blackman–Tukey). What Paper A does, so it is the
  paper-faithful route, and the VACF is a useful diagnostic in its own right — whether the
  correlation has actually decayed within `τ_max` is visible there and not in a spectrum.

Welch is the default for one decisive reason: it averages outer products `V*V^T`, so every
`P_i(ω)` is **positive semi-definite by construction**. Blackman–Tukey does not guarantee
this. For Method 1 a negative excursion is a cosmetic wiggle; for Method 2 a non-PSD
`P_i,αβ` means a negative mean-square displacement along some direction, which is not a
physical tensor. PSD-ness is asserted on the IR, not assumed.

The off-diagonal terms motivate the same choice. `C_αβ(t)` is not even in `t` for α≠β, so a
symmetric-lag VACF implementation silently discards the quadrature part — turning ambiguity
A2 into an accident rather than a decision. Welch keeps both parts explicit.

Segment length sets frequency resolution (`1/(L·dt)`); overlap and window (Hann by default)
set variance and spectral leakage. For the VACF route the corresponding knob is `τ_max`.
Both are documented, testable parameters recorded in the IR, rather than the undocumented
choices they are in both source codes.

Welch also yields the spread across segments, which gives an uncertainty estimate on the
spectrum and a stationarity check for free. §4's agreement tolerances should be stated
relative to it rather than picked by hand.

**The Hermitian question.** `P_αβ` is Hermitian in general, but the displacement tensor
`B` it stands in for is real symmetric. Paper B prints the full Hermitian matrix and never
says what happens to the imaginary (quadrature) part. We take the real part — the
co-spectrum — and assert the imaginary part is negligible in the harmonic limit as a
diagnostic. See §5.

**Sum rule.** Classical equipartition gives a hard, testable invariant:

```
∫ tr P_i(ω) dω = ⟨|v_i|²⟩ = 3 k_B T_MD / m_i
```

This is the normalisation anchor the papers lack (see method-review §3.1), and it becomes
a unit test that runs on every estimator change.

### [C] Intermediate representation — numpy dataclass, HDF5 on request (D3)

The IR is a lightweight dataclass holding numpy arrays, with `to_hdf5` / `from_hdf5` for
persistence. The dataclass *is* the schema: no separate schema document to drift out of
sync, and the in-memory object is the thing tests exercise.

Fields:

- `frequencies` — 1-D grid, with units
- `spectral_density` — per-entity real symmetric tensor, **absolute** units (D6)
- `species`, `masses`, `entity_counts` — entity → physical meaning
- `temperature_md`, `dt`, `n_frames`, estimator and segment/window parameters, ensemble
- `normalisation` — convention identifier, `"absolute-spectral-density"` by default
- `provenance` — source identity, code version, pre-processing applied

Stored as absolute spectral density rather than a normalised lineshape plus a separate
amplitude (D6): the sum rule in [B] then applies directly to the stored array, the units
are dimensionally checkable, and the IR cannot reproduce the failure mode that makes
Paper A's absolute scale irrecoverable. Conversion to a normalised lineshape is a
one-line helper, and is needed anyway for §4 — Euphonic and OCLIMAX each emit pDOS under
their own conventions, so comparison happens in normalised space even though storage does
not.

Granularity is **per-atom**, with optional grouping applied at write time. Per-species
alone is insufficient: the Debye–Waller factor is genuinely per-atom in a disordered
system, which is the entire point of the method.

The memory ceiling is real and unchanged by the container choice — 3000 atoms × 25000
bins × 6 components in float64 is ~3.6 GB, and both source papers hit exactly this wall.
Mitigations (frequency rebinning for the stored IR, entity grouping, float32 storage) are
options on the dump path, not changes to the in-memory contract.

### [D] Displacement model

Convert classical spectral density to quantum harmonic displacement amplitudes:

```
u²_i,αβ(ω) = ħ P_i,αβ(ω) / (ω k_B T_MD)
```

(Paper B Eq. 8; follow the discrete-transform bookkeeping in its SI Eq. S13, not Eq. 8 as
printed — see method-review §4.1.)

Method 1 takes `tr(u²)/3` here and proceeds isotropically. Method 2 keeps the tensor.

Debye–Waller factor `A_i = ∫ u²_i(ω) dω`, with Bose–Einstein occupation at the
experimental temperature.

**Two temperatures, defaulting to equal (D8).** `T_MD` is a property of the IR — it is the
temperature at which the PES was sampled, and it is what appears in the classical-to-quantum
conversion above. `T_experiment` is an argument to stages D and E, setting Bose–Einstein
occupation and the Debye–Waller factor. They default to equal, and the API should make
supplying only one the easy path. Decoupling them is one of Paper A's real contributions —
it lets a "fictive" MD temperature be chosen to sample the right region of the PES while the
spectrum is evaluated at the temperature the measurement was actually made at — but it is
an approximation with no error estimate attached, so using it must be a deliberate act,
recorded in provenance and surfaced in the output metadata.

The `1/ω` factor diverges as ω→0. Needs an explicit low-frequency cutoff, and diffusive or
quasi-elastic contributions must be excluded rather than silently integrated. This is a
correctness issue neither paper addresses.

### [E] Scattering model — iterative convolution (D7)

Incoherent `S(Q,ω)`, resolved by excitation order, built as `S_n = S_1 * S_{n-1}`
exactly as both papers do. Convolution happens **at atomic level, before** cross-section
and Debye–Waller factor are applied (Paper A is explicit about this ordering).

Two consequences of choosing convolution as the only route, both of which need handling
rather than just noting:

- **The order prefactor becomes a first-class choice, not an implementation detail.**
  Paper B asserts `1/n!` on the grounds that it "closely resembles the known prefactors";
  the exact aCLIMAX prefactors differ. With no exponential resummation to derive them
  from, the prefactor must be selected explicitly, documented, and made configurable so
  the two conventions can be compared.
- **We lose the internal cross-check.** The exponential form would have provided an
  independent, exact result to validate the convolution against, including quantifying
  ambiguity A8. Verification now rests entirely on the Tier 1 analytic tests and external
  oracles, which raises the importance of both.

Maximum order is a parameter; AbINS uses 10 as its ceiling, which is a reasonable default.

Powder averaging is a strategy chosen per method:
- **isotropic** (Method 1)
- **almost-isotropic** for the fundamental, isotropic above (Method 2, following aCLIMAX)

Both are small, well-defined functions over `tr(A)`, `tr(B)` and `B:A`.

### [F] Instrument model — Euphonic in core, AbINS optional (D4)

Two pieces, both small:

- **Kinematics.** For indirect geometry, `Q² = (2m_n/ħ²)(E_i + E_f − 2√(E_i E_f) cos θ)`
  with `E_i = E_f + ħω`. Direct geometry is the same relation with `E_i` fixed. Roughly a
  hundred lines from Squires.
- **Resolution.** Gaussian with an energy-dependent width, in practice a polynomial in
  energy transfer.

Broadening and Sears cross-section data come from Euphonic (`broadening`, `isotopes`),
which is a core dependency. Instrument definitions belong in data files (TOML/YAML), not
code, so adding a spectrometer needs no release; parameter values may be taken from AbINS
with attribution, now that the licences are compatible. AbINS itself is an optional extra
used to cross-check our broadening, not to supply it.

The published resolution constants for TOSCA disagree between OCLIMAX and AbINS by a
factor of ten (method-review §3.7). That has to be resolved explicitly when the data file
is written — not settled by whichever source we happened to copy from.

### [G] Output

A spectrum container with units and metadata, per-order decomposition retained (both
papers show order-resolved plots, and it is diagnostically essential). Export to a form
the neutron ecosystem already reads.

---

## 3. Cross-cutting rules

**Units.** One canonical internal system — Å, ps, amu, meV — converted at the boundaries.
Unit confusion is the most likely source of a plausible-looking wrong answer in this
domain, so it is defended against by documenting the unit of every stored array and by
testing absolute normalisation against equipartition (D6), not by a unit library.

`pint` was considered for the public API and rejected during implementation. It does not
survive contact with the arrays that matter: the IR's `(n_entity, n_freq, 3, 3)` density
is handed straight to `numpy.linalg.eigvalsh` and written to HDF5, neither of which takes
a quantity, so units would be stripped at exactly the boundary where they were supposed
to help. It would also put a hard dependency in every downstream import for a guarantee
the sum-rule check already provides end to end. Units are instead recorded in the IR's
metadata and asserted by tests.

**Optional dependencies.** `euphonic` is used only to validate M2 against a phonon
calculation and is an extra rather than a runtime requirement; nothing on the main path
imports it.

**Provenance is not optional.** Every artefact records what produced it. The IR is meant
to be shared and re-analysed; an unlabelled pDOS with an unknown normalisation convention
is worse than useless.

**Per-atom independence.** Stages B, D and E are embarrassingly parallel over atoms. The
design should not gratuitously destroy that. Prefer chunked shared-memory parallelism over
the MPI approach MolDyINS took — this is a library, not a batch job.

**Streaming.** Welch segmentation means stage B can process a trajectory in chunks and
never hold it in memory. Keep that property.

---

## 4. Validation

### Tier 1 — analytic

A single 1D harmonic oscillator: `P(ω)` is a delta at `ω₀` with amplitude fixed by
equipartition. A diatomic adds a second known mode. These catch normalisation, units and
factor-of-2π errors immediately and run in milliseconds.

Given D7 removed the exponential cross-check, this tier now carries more of the
verification load than originally planned. It should be extended to cover the multiphonon
path specifically: for a single oscillator the exact *n*-quantum intensities are known in
closed form, so the convolution and its prefactor can be checked against analytic values
rather than only against other codes.

### Tier 2 — harmonic reference via Euphonic

The core validation, as proposed:

```
ASE LennardJones calculator
   ├── force constants (finite displacement) → Euphonic ForceConstants
   │        → calculate_pdos → reference atom-projected DOS
   └── MD (same calculator, same cell, low T) → our pipeline → pDOS
```

Euphonic's `calculate_pdos` returns a `Spectrum1DCollection` with optional neutron
cross-section weighting, which makes it a direct comparison target. The MD side stays
within ASE throughout, so this tier is unaffected by D2.

Two subtleties that will otherwise produce a spurious mismatch:

- **Brillouin-zone sampling must match.** An MD supercell samples only the q-points
  commensurate with it. The Euphonic reference must use that same commensurate q-grid, not
  a dense mesh, or the reference will be smoother than anything MD can produce.
- **Stay in the harmonic limit.** Low temperature, small displacements. LJ is anharmonic;
  at high T the two will legitimately disagree and the test tells us nothing.

Agreement criterion: integrated intensity per species, peak positions, and a normalised
spectral overlap. Tolerances are derived from the inter-segment spread of the Welch
estimate (D5) rather than chosen by hand — a mismatch only counts as a failure if it
exceeds the sampling noise of the MD side. Comparison is done on normalised lineshapes,
since Euphonic's convention differs from the IR's absolute one (D6).

### Tier 3 — reference implementations and experiment

- **MolDyINS**: its bundled 25 K trajectory is a ready-made regression fixture, but it is
  a GROMACS `.trr`, which ASE cannot read. It needs one-off conversion to extxyz or
  `.traj` — a developer-side step using MDAnalysis or `gmx`, checked in as a fixture
  rather than added as a runtime dependency.
- **OCLIMAX**: numerical oracle for the isotropic method.
- **Published spectra**: ice Ih, MgH₂ (strong multiphonon), P3HT.

---

## 5. Known ambiguities

Physics and numerics questions the papers leave open. These need resolving during
implementation but are **not** the human decisions in §6 — they have defensible answers we
can adopt and document.

| # | Ambiguity | Working resolution |
|---|---|---|
| A1 | Absolute normalisation of the pDOS is not determined by Paper A (its Eq. 5 normalises the VACF, its Eq. 7 equates the integral to a squared displacement) | Settled by D6: store absolute spectral density, anchor with the equipartition sum rule, make it a unit test |
| A2 | Paper B's `P_αβ` is Hermitian but the displacement tensor is real symmetric; fate of the imaginary part unstated | Take the real part (co-spectrum); assert the imaginary part is negligible as a diagnostic. D5's Welch default keeps this an explicit discard rather than an implicit one |
| A3 | Paper B Eq. 19's Boltzmann factor has a positive exponent and no partition function; diverges as printed | Treat as a typo, use `exp(−ħω/kT)/Z`, flag in docs |
| A4 | Paper B Eq. 21's Debye–Waller `max[...]` heuristic has no derivation | Implement as published, behind a flag; prefer the consistent BE-weighted form as default |
| A5 | Paper A never states what OCLIMAX does with anisotropy given only `tr(B)` from MD | Assume isotropic fall-back; confirm against the OCLIMAX oracle |
| A6 | `g(ω)/ω` divergence and quasi-elastic contamination at low ω | Explicit cutoff, surfaced as a required parameter, not a hidden constant |
| A7 | Multiphonon order prefactor: Paper B asserts `1/n!`, aCLIMAX's exact prefactors differ | Promoted to a configurable, documented parameter (see [E]); validate against analytic single-oscillator intensities |
| A8 | Rank-4 → rank-2 contraction in Paper B Eq. 17 "approximated as matrix multiplication" | Implement as published. With the exponential route ruled out by D7, the error can only be bounded against an external oracle, not internally |

---

## 6. Decisions

### Resolved

**D1 — Licence: relicense to GPL-3.0.**
The core may depend on Euphonic, AbINS and MDANSE directly, and may reuse their code.

Done with the copyright holder's approval: `LICENSE` now holds the verbatim GPL-3.0 text.

*Outstanding housekeeping:* the GPL `LICENSE` is the licence text only and carries no
project copyright line, so the `Copyright (c) 2026, Alin Marin Elena` statement that the
BSD file used to hold must be restated elsewhere — the standard GPL header on each source
file, and a licence section in the README. Packaging metadata needs
`License-Expression: GPL-3.0-or-later` (or `GPL-3.0-only`; the "or later" choice is worth a
moment's thought, as it is not recoverable once distributed).

*Consequence:* unblocked D4 below — Euphonic became a permissible core dependency.

**D2 — Trajectory reader: ASE only.**
One dependency, tight integration with the calculator ecosystem used for validation, no
abstraction layer.

*Consequence:* GROMACS users — the audience of Paper B — cannot bring `.trr` files
directly, and neither can the MolDyINS regression fixture. Conversion to extxyz or `.traj`
becomes a documented prerequisite, and the fixture is converted once offline and checked
in. Worth revisiting if external users hit it; a reader protocol can be retrofitted, just
less cheaply than building it now.

**D3 — IR: lightweight numpy-backed dataclass, with HDF5 load/dump.**
The dataclass is the schema. HDF5 is I/O, not the definition. Granularity stays per-atom
with optional grouping on write.

**D7 — Multiphonon: iterative convolution only.**
Matches both papers and the reference codes, so any discrepancy against OCLIMAX or
MolDyINS is unambiguously ours rather than a method difference.

*Consequence:* the order prefactor becomes an explicit configurable choice, and Tier 1
analytic tests must be extended to cover multiphonon intensities — see [E] and §4.

**D5 — Spectral estimator: Welch default, VACF as a second estimator.**
Chosen on positive semi-definiteness of the tensor rather than on accuracy — the two agree
where both are valid.

*Consequence:* PSD-ness becomes an assertion on the IR, not an assumption. The two
estimators run against each other as a cross-check on stage B. Welch's inter-segment spread
supplies the uncertainty estimate that §4's tolerances should be stated against, replacing
hand-picked thresholds.

**D6 — IR normalisation: absolute spectral density.**
Equipartition sum rule as the invariant. `normalisation` field retained so the convention
travels with the data.

*Consequence:* a normalised-lineshape conversion helper is needed regardless, because
Euphonic and OCLIMAX each normalise differently and Tier 2/3 comparison happens in that
space.

**D8 — Temperature: `T_MD` and `T_experiment` supported separately, defaulting to equal.**
`T_MD` belongs to the IR; `T_experiment` is an argument to stages D and E.

*Consequence:* decoupling is an approximation without an error estimate, so it must be
opt-in, recorded in provenance, and visible in output metadata — never a silent default.

**D4 — Instrument layer: Euphonic in core, own kinematics, AbINS optional.**
Euphonic supplies broadening and Sears cross-section data; the instrument kinematics are
ours; AbINS sits behind an optional extra. The licence objection disappeared with D1, but
Mantid remains a very heavy dependency for the amount of it we would use, and it is
awkward to install outside a conda environment.

*Consequence:* AbINS's calibrated instrument parameters are the most valuable thing we are
declining to depend on, and our instrument data files will partly duplicate them. Since
AbINS is now licence-compatible, those values can be copied into our TOML/YAML definitions
with attribution rather than re-derived — which is also the cleanest way to confront the
factor-of-ten TOSCA resolution disagreement (method-review §3.7) rather than silently
inheriting one side of it. The optional AbINS extra then exists to cross-check our
broadening against theirs, not to supply it.

### Open

None. All eight decisions are settled; the remaining open items are the ambiguities in §5,
which are resolved during implementation rather than by decision.

---

## 7. What is deliberately not here

Module layout, class names, function signatures and package structure. Those follow from
§1 and are cheap to change; the stage contracts, the single-IR decision, and §6 are not.
