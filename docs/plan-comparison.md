# Plan comparison: `main` vs `origin/openspec`

**Produced:** 2026-09-24, by a read-only review of both branches. Nothing was checked
out, merged or modified; all branch content was read with `git show`. The numbers and
quotations below were accurate at the tip commits recorded in the appendix, and will go
stale as either branch moves.

## Verdict — the four things that actually matter

1. **These are not two versions of one plan; they are two independent plans with no
   common ancestor text.** The merge base (`8b6848d`) contains only `LICENSE`. Every
   planning document on each branch was written from scratch against the same seed
   prompt (preserved verbatim on openspec as `proposal-prompt.md`).
   `git diff main..origin/openspec -- '*.md'` therefore shows pure add/delete:
   `proposal.md`, `design.md`, `plan.md`, `docs/method-review.md`, `docs/validation.md`
   deleted; the `openspec/` tree added. No shared file except `README.md`, which is
   wholly rewritten.

2. **The single most consequential technical divergence is the intermediate
   representation.** `main` design.md §1 ("The fork point is B→C, and there is only one
   IR") computes the full 3×3 velocity cross-spectral tensor **once**, stores it in
   absolute units (D6, equipartition sum rule `∫ tr P dω = 3k_BT/m`), and makes the
   isotropic method a consumer that takes the trace. openspec splits the tensor out into
   Stage 3 (`ins-tensor-anisotropic`, Decision 1: "Extend `correlation.py` to evaluate
   the 6 independent elements…"), and its Stage 1 spec **mandates the opposite
   normalisation**: "SHALL normalize the pDOS such that the integral over positive
   frequencies equals 1.0" (`openspec/specs/trajectory-pdos/spec.md`, *Atom-Projected
   Phonon Density of States*). `main`'s method-review §3.1 identifies exactly that
   normalised-VACF convention as the defect in Paper A that makes its absolute scale
   irrecoverable. openspec has adopted the paper's bug as a requirement.

3. **`main` carries a large body of analysis that openspec has no counterpart for, and
   would simply lose.** `docs/method-review.md` (223 lines: the two papers plus Paper A′,
   eight enumerated gaps in Paper A, eight in Paper B, the literature table, the state of
   MolDyINS/OCLIMAX) and design.md §5's ambiguity register A1–A8 have no analogue
   anywhere in `openspec/`. Most notably method-review §2 ("The most useful thing in
   either paper") — the SI Eqs. S17–S24 proof that the convolution series resums to an
   exponential, i.e. the Gaussian-approximation `I(q,t)` — is a genuine finding that
   openspec never mentions, and whose absence means its Stage 2/3 multiphonon prefactors
   (`1/n!`, `3^{n-2}/(n!5^{n-1})`) are taken on faith from the papers rather than treated
   as an open question.

4. **openspec's form is machine-checkable where `main`'s is prose.** openspec has a
   living capability spec in requirement/scenario (SHALL + WHEN/THEN) form, a per-change
   `tasks.md` with ticked checkboxes traceable to test files, and an archive convention
   (`openspec/changes/archive/2026-09-24-trajectory-pdos/`) that snapshots the delta that
   produced the current spec. `main` has richer reasoning but its only task decomposition
   (`docs/issues/01..17`) is **untracked in git** — present in the working tree, not on
   the branch. Structurally, openspec is ahead on traceability; substantively, `main` is
   far ahead on physics.

Both plans are for the same three capabilities in the same order, so the disagreement is
not about *what* to build. It is about normalisation, about the IR, and about what counts
as proof.

---

## Scope and goals

Identical at the top level, and both trace to the same prompt: pDOS from trajectories →
Cheng isotropic → Harrelson anisotropic.

| | `main` | `openspec` |
|---|---|---|
| Stage 1 | proposal.md "What" §1, *Atom-projected density of states from trajectory data* | `openspec/changes/archive/2026-09-24-trajectory-pdos/` (implemented, archived) |
| Stage 2 | §2 *INS intensity calculation — isotropic method (Cheng et al.)* | `openspec/changes/ins-isotropic-sim/` (proposal + design only, no tasks) |
| Stage 3 | §3 *…anisotropic method (Harrelson et al.)* | `openspec/changes/ins-tensor-anisotropic/` (proposal + design only, no tasks) |

Explicit non-goals agree where they overlap: both exclude coherent single-crystal
scattering and running production MD (main proposal.md "Not in scope"; openspec Stage 1
design "Non-Goals", Stage 3 "Non-Goals").

**Only on `main`:**

- A stated non-goal of *"A direct trajectory-to-spectrum route"* and of *"New physics"* —
  including the deliberate decision to leave Paper B's derived-but-unused anharmonic
  correction alone. openspec Stage 3 makes the same call in passing ("Fitting anharmonic
  damping rates… Non-Goals") but does not name it as a scope boundary for the project.
- Deliverables that are not physics at all: a CLI (`mdins pdos`, `mdins spectrum`,
  plan.md §3), HDF5 persistence of the IR, provenance records ("Provenance is not
  optional", design.md §3), instrument definitions as data files (TOML/YAML, design.md
  [F]), and CI (plan.md §1: "Python 3.11–3.13, ruff + mypy + fast tests on every PR, slow
  and reference tests nightly"). None of these appear anywhere in the openspec plan;
  openspec's `pyproject.toml` has no mypy config and the branch has no `.github/`.
- *"What success looks like"* includes a documentation deliverable — "The ambiguities we
  had to resolve are written down" — restated at the foot of plan.md §7 as "a
  deliverable, not a byproduct".

**Only on openspec:**

- A first-class *velocity-derivation-from-positions* feature. Stage 1 spec, *Trajectory
  Data Ingestion*: "When velocities are missing… the system SHALL compute velocities by
  central differences… remove the first and last steps… and emit a visible warning."
  `main` design.md [A] considers and demotes this: "Deriving velocities by
  finite-differencing positions is possible but introduces a frequency-dependent transfer
  function (a `sinc` factor) that must be divided out, and requires periodic-image
  unwrapping. It is a fallback, not a default." openspec's design and spec mention neither
  the sinc correction nor image unwrapping. This is a scope difference *and* a correctness
  difference.
- Kinetic-temperature validation as a named requirement (`T_traj = 2⟨E_k⟩/3Nk_B`, warn at
  20%). `main` has no equivalent requirement; its IR carries `temperature_md` but the plan
  never proposes checking it against the trajectory. This is a small, real thing `main`
  lacks.

---

## Architecture and decisions

### Do D1–D8 survive? Mostly not, and three are resolved differently

`main`'s design.md §6 has eight numbered decisions, all marked resolved ("Open: None").
openspec has three separate unnumbered `## Decisions` lists, one per change (5 + 6 + 4 =
15 decisions), scoped to the stage rather than the project.

| `main` | openspec counterpart | Same answer? |
|---|---|---|
| **D1** Relicense to GPL-3.0 | **none** | **No decision at all.** openspec's `LICENSE` is still "BSD 3-Clause License / Copyright (c) 2026, Alin Marin Elena" and its `pyproject.toml` has no licence field. `main` relicensed with the copyright holder's approval and tracked the outstanding housekeeping. Since openspec plans to depend on `abinslib` (GPL-3.0) for Stage 2/3 benchmarking, this is a latent licence conflict its plan never confronts. |
| **D2** ASE only, with the velocity-format table and the "GROMACS users cannot bring `.trr`" consequence | Stage 1 Decision 1, "Trajectory Ingestion via ASE with In-Memory Fallback" | **Partly.** Same reader choice; openspec adds an in-memory array path (reasonable, and `main` has no equivalent) but omits the format-capability table and answers the velocity-availability problem with finite differences instead of documenting a conversion prerequisite. |
| **D3** IR = frozen numpy dataclass + HDF5, per-atom granularity | Stage 1 Decision 1's `TrajectoryData` dataclass | **No.** openspec's dataclass is the *trajectory* container, not a persisted spectral IR. There is no HDF5, no schema, no persisted intermediate; Stage 2 consumes `g_d(ω)` in memory. `main`'s claim that the IR "is meant to be shared and re-analysed" has no openspec analogue. |
| **D4** Euphonic in core, own kinematics, AbINS optional extra | Stage 2 Decisions 1 and 6: Euphonic `IsotopeData` for cross sections; `abinslib` as the benchmark | **Different weighting.** Both take Sears data from Euphonic and write their own kinematics (main design.md [F]; openspec Stage 2 Decision 2, abstract `Kinematics` base class with `TOSCAKinematics`/`FixedQKinematics`). But openspec makes `abinslib` the *primary reference oracle* for both intensity stages. `main` demotes AbINS to "cross-check our broadening, not to supply it" and plan.md §1 records "`abins` is conda-only and cannot be an extra" — which, if true, undercuts openspec's whole Stage 2/3 validation plan. |
| **D5** Welch default, VACF second, chosen on positive semi-definiteness | Stage 1 Decision 2, "In-House Vectorised Correlation Engine": Wiener–Khinchin VACF, zero-padded, Hann window | **Resolved oppositely.** openspec implements only the Blackman–Tukey-style route that `main` demotes. `main`'s reason is explicit and physical: "Welch averages outer products `V*V^T`, so every `P_i(ω)` is positive semi-definite by construction. Blackman–Tukey does not guarantee this… a non-PSD `P_i,αβ` means a negative mean-square displacement along some direction". openspec Stage 3 will need exactly that tensor to be PSD and its plan contains no PSD assertion — only a Hermiticity claim. `main` additionally uses the Welch inter-segment spread as its uncertainty estimator; openspec has no uncertainty estimate anywhere, which is why its tolerances are hand-picked. |
| **D6** Absolute spectral density, equipartition sum rule | Stage 1 spec requirement: `∫₀^∞ g_d(ω) dω = 1.0` (±1%) | **Resolved oppositely**, and this is the deepest disagreement. See verdict item 2. |
| **D7** Multiphonon by iterative convolution only, prefactor a configurable documented parameter | Stage 2 Decision 4 (`scipy.signal.fftconvolve`), Stage 3 prefactor `3^{n-2}/(n!5^{n-1})` | **Same mechanism, different status.** Both convolve. `main` flags the prefactor as ambiguity A7 — "Paper B asserts `1/n!`… the exact aCLIMAX prefactors differ" — makes it configurable, and adds an analytic oracle to decide it. openspec states prefactors as settled fact. openspec Stage 2 also defaults `N_max = 4`; `main` design.md [E] uses AbINS's ceiling of 10. |
| **D8** `T_MD` and `T_experiment` separate, defaulting to equal, opt-in and recorded in provenance | implicit only | **Not decided.** openspec Stage 2 uses `T_INS` in the Bose factor and DWF; Stage 3 uses `T_MD` in `B̄̄ᵢ(ω)`. The two temperatures therefore exist, but the decoupling is never named, never justified, never flagged as an approximation without an error estimate, and Stage 1's spec offers no field to carry `T_MD` into Stage 2. |

### Decisions on openspec with no counterpart on `main`

Three are worth keeping regardless of which branch wins:

- **Stage 1 Decision 5, the ASE→Phonopy→Euphonic bridge.** "Because Euphonic supports
  CASTEP and Phonopy formats for force constants but not ASE's native `ase.phonons`,
  implement a lightweight helper… following the pattern in Calorine, without introducing
  Calorine as a dependency", exporting via `ForceConstants.from_phonopy()`. `main`
  plan.md §5 Layer 3 just draws an arrow "finite displacement → Euphonic ForceConstants"
  and does not say how the bridge is built. `main` evidently solved this
  (`tests/lj_reference.py`) but never recorded the route as a decision.
- **Stage 3 Decision 2, regularisation of the double-dot contraction.** "If `Tr(B̄̄ᵢ(ω))`
  falls below a numerical threshold (e.g. 10⁻¹²), regularize `αᵢ(ω) → ⅓Tr(Āᵢ)`." `main`
  handles the `1/ω` divergence at low frequency (A6, an explicit cutoff) but has no
  equivalent for the `B:A/Tr(B)` denominator in the almost-isotropic average — a distinct
  division-by-zero in a different place.
- **Stage 2 Decision 5, broadening decoupled from the calculation.** "INS simulation
  outputs unbroadened spectra. A separate `broadening` module…". `main` puts broadening
  inside stage [F]; openspec's separation is arguably cleaner and makes the resolution
  question testable in isolation.

Conversely, openspec has nothing corresponding to `main`'s cross-cutting rules in
design.md §3: the canonical Å/ps/amu/meV unit system and the explicit, reasoned rejection
of `pint` ("It does not survive contact with the arrays that matter… units would be
stripped at exactly the boundary where they were supposed to help"); per-atom parallelism;
streaming. openspec's unit story is one line in its README ("velocities in Å/fs, positions
in Å, masses in amu, timesteps in fs, energies in meV") — note **fs, not ps**, so the two
branches do not even share a canonical time unit.

### Pipeline stages A–G

`main` design.md §1's seven lettered stages with named data contracts have no openspec
equivalent. openspec's decomposition is by *module within change*: Stage 1 →
`trajectory.py`/`correlation.py`/`benchmark.py`; Stage 2 →
`cross_sections.py`/`kinematics.py`/`isotropic.py`/`broadening.py`; Stage 3 →
`correlation.py` (extended) + `anisotropic.py`. This is a finer and more concrete
decomposition than `main`'s §1, but it is a *file* decomposition, not a contract
decomposition — and note that `main` design.md §7 explicitly refuses to put module layout
in the design document, deferring it to plan.md §1. So the two branches split the same
material along a different line: openspec's per-change design.md ≈ `main`'s design.md
§§1–2 *plus* plan.md §1.

The consequence of openspec having no stage contracts is visible in Stage 3, which must
*reopen* `correlation.py` ("Expands `correlation.py` to produce 3×3 Cartesian tensor
representations"). `main`'s single-IR decision exists precisely to avoid that: "the two
methods differ only in stages D and E, share A–C and F–G, and can be compared on
byte-identical input."

---

## Milestones, sequencing and dependencies

Both plans run Stage 1 → isotropic → anisotropic. The differences are granularity,
parallelism and what "done" means.

**`main`** — plan.md §4, six milestones M0–M5 (Scaffolding, pDOS from MD, Harmonic
validation, Method 1, Method 2, Benchmarks) with relative sizes and an explicit commitment
line: **"M0–M2 is the honest first cut… Treating M5 as in-scope for a short effort is how
this ends up with a spectrum nobody has verified."** §4 also names three parallel tracks
after M0 (A–C critical path; instrument stage F, "completely self-contained"; validation
harness) and states the dependency rule that makes them work: "D–E depend on C existing
but not on C being correct, so they can start once the dataclass is frozen and fed
synthetic data." That is why plan.md §3 is titled "The IR is the contract — write it
first", and why `src/mdins/scattering.py` is 78 lines of documented `NotImplementedError`
stubs.

**openspec** — three sequential changes, no parallelism, no sizing, no risk ranking, no
committed cut. The dependency statement is a single sentence per change. Within Stage 1
there is finer sequencing than `main` offers anywhere on-branch: `tasks.md` has 14
numbered tasks in 4 groups, each carrying its own verification clause (e.g. "2.3 Implement
central-difference velocity derivation fallback… verified by unit tests comparing derived
velocities against harmonic oscillator positions"), all ticked `[x]`. Stages 2 and 3 have
**no `tasks.md` at all** — they are proposals awaiting decomposition.

**Definition of done.** `main` plan.md §7 gives one per milestone, in testable terms ("M0 —
CI green, IR round-trips through HDF5, entry points importable and raising
`NotImplementedError`"). openspec's equivalent is the task checkbox plus the spec
scenarios. `main`'s M2 entry has been rewritten in place to record the achieved result —
i.e. `main` uses the plan document itself as the progress record, where openspec uses
checkbox state plus an archive directory.

`main`'s §6 risk table (7 risks with mitigations, including "Scope creep into M5 → M0–M2
is the committed cut") has no project-level counterpart; openspec's risks are per-change
"Risks / Trade-offs" pairs and are all technical (memory footprint, convolution grid edge
effects, FFT cost), never programmatic.

---

## Testing and validation strategy — the sharpest substantive gap

**`main`** defines three validation tiers in design.md §4 and five test layers in plan.md
§5, and is unambiguous about which comparison is decisive.

- Layer 1 (analytic, synthetic velocities, always run): seven named cases including
  "Correlated x/y motion → known off-diagonal term, and known *sign*" and "Circular motion
  → nonzero quadrature part, which is where ambiguity A2 becomes visible", plus four
  asserted invariants (equipartition, non-negative eigenvalues of every `P_i(ω)`,
  Welch-vs-VACF agreement within the inter-segment spread, exact HDF5 round-trip).
- Layer 2: **"the single most valuable test in the suite and should be written before the
  convolution code"** — the closed-form isotropic 3D oscillator
  `S(Q,ω) = e^{-2W} Σ_n [(Q²u²)ⁿ/n!] δ(ω − nω₀)`, used as an analytic oracle that
  *decides* the multiphonon prefactor (A7) rather than arguing about it.
- Layer 3: Einstein crystal first (exact answer, isolates code bugs from physics), then the
  LJ/Euphonic Tier 2 comparison.
- **Tolerances are derived, not chosen.** design.md §4: "Tolerances are derived from the
  inter-segment spread of the Welch estimate (D5) rather than chosen by hand — a mismatch
  only counts as a failure if it exceeds the sampling noise of the MD side."
  `docs/validation.md` carries this through with a table of MD vs Euphonic moments, 1σ
  columns, a 3σ threshold and a negative control asserting that a 5% frequency stretch is
  rejected. It also documents its own weak spots.
- Layer 5 oracles: OCLIMAX (Method 1), MolDyINS (Method 2, "Needs repair first — it fails
  on modern NumPy (`np.float`)"), AbINS broadening, and published spectra (ice Ih, MgH₂,
  P3HT) — explicitly *not* CI, "results written up rather than asserted".

**openspec** has one validation idea, executed well, and hand-set tolerances everywhere
else.

- Stage 1's decisive comparison is the same LJ-argon/Euphonic closed loop, on the same
  commensurate q-grid, and it is genuinely done — plus a mass-disordered check (light 40
  amu / heavy 160 amu on one FCC lattice) that is the direct analogue of `main`'s Ar/Kr
  L1₀ case and rests on the same insight (a species-blind LJ potential makes mass the only
  asymmetry).
- But the acceptance criterion moved during implementation and the plan records the retreat
  honestly. The spec scenario *Closed-loop Lennard-Jones benchmark execution* still says
  "within 5% relative error", while `tasks.md` 4.4 notes: "cross-method spectral comparison
  is hard to assert tightly; the test is **regression-based** (fixed seed → deterministic
  baselines) with a coarse general-position peak check (~15%)… the assertion is
  regression-pinned rather than a strict 5% cross-method gate." So the living spec and the
  implementation now disagree — the exact drift OpenSpec's structure is supposed to prevent.
- The comparison statistic is **peak positions** (Gaussian smoothing + prominence
  thresholding); `main` compares **moments** (⟨E⟩, ⟨E²⟩) per species with propagated
  uncertainties, plus coarse bins, plus a bound on out-of-band weight. `main`'s statistic is
  the stronger one.
- No analytic layer. There is no synthetic-velocity oracle, no equipartition sum rule
  (foreclosed by the ∫g = 1 normalisation), no PSD assertion, no 2π/unit trap test.
  openspec's Stage 1 correctness rests on `VACF(0) = 1.0`, `∫g dω = 1`, and peak positions
  — all three of which are invariant under a global frequency-axis or mass error of the
  kind `main`'s Layer 1 exists to catch. (VACF(0)=1 holds by construction of the
  normalisation, so it proves nothing about scale.)
- Stage 2/3 validation is entirely "compare against `abinslib`", plus one good internal
  consistency test (Stage 3 Decision 4: off-diagonals < 1% of trace on cubic argon,
  `Tr(B)/3` matches Stage 1's scalar `g_i(ω)`, and `S(Q,ω)` matches Stage 2 within 1%).
  `main`'s M4 done-criterion is the same idea in one line.
- No experimental comparison anywhere. `main` names ice Ih, MgH₂ and P3HT; openspec's plan
  never proposes checking against a measured spectrum.

---

## Form: what the OpenSpec structure buys and costs

openspec uses the OpenSpec spec-driven convention:

```
openspec/config.yaml                                   (schema: spec-driven; all guidance commented out)
openspec/specs/<capability>/spec.md                    living spec: Purpose + Requirements + Scenarios
openspec/changes/<change-id>/{proposal,design,tasks}.md + specs/<cap>/spec.md   (a "Spec Delta")
openspec/changes/archive/<date>-<change-id>/...        change snapshot after it lands
```

`main` uses four hand-written long-form documents with explicit cross-references
(`README.md` → `proposal.md` / `design.md` / `plan.md` / `docs/method-review.md`;
design.md header: "Companion to proposal.md (why/what) and docs/method-review.md
(source-paper analysis)"; plan.md opens "Concrete work breakdown for the architecture in
design.md. That document deliberately stopped short of module layout and signatures; this
one supplies them"). Plus `docs/validation.md`, which is a written-up result rather than a
plan, and the 17 untracked `docs/issues/*.md`.

**What openspec's structure buys.**

- *A living spec distinct from the change that produced it.*
  `openspec/specs/trajectory-pdos/spec.md` is the current contract; the delta that created
  it is preserved under `changes/archive/…` with `## ADDED Requirements`. `main` has no
  object that answers "what does the system currently guarantee?" — it has a design
  (intent) and a validation write-up (evidence), with plan.md §7 patched in place to bridge
  them.
- *Near-machine-readable acceptance criteria.* SHALL/WHEN/THEN scenarios map close to
  one-test-each.
- *Explicit change lifecycle.* Proposal → design → tasks → implement → archive. Review can
  happen at proposal stage, before design effort is spent; stages 2 and 3 are sitting in
  exactly that state.
- *Task-level completion state in the repo*, on-branch and diffable.

**What it costs, as actually used here.**

- *Triplication.* Each change repeats its stage-numbered context header nearly verbatim.
  The archived spec delta is byte-identical to the living spec, so the "delta" carries no
  information on a greenfield change.
- *Nowhere for cross-cutting or non-functional material.* Units, licence, provenance,
  parallelism, memory strategy, CI, and the source-paper analysis have no home in a
  per-capability spec, and consequently do not exist on openspec. `config.yaml` has a
  `context:` field designed for exactly this and it is left as commented-out boilerplate.
- *The spec can drift anyway.* The 5% vs 15% mismatch above is already live.
- *Requirement prose is weaker than the reasoning it replaces.* "SHALL provide a robust
  default Hann window" is a requirement; "Welch is the default for one decisive reason: it
  averages outer products `V*V^T`, so every `P_i(ω)` is positive semi-definite by
  construction" is an argument you can audit. The scenario form tends to record the *what*
  and drop the *why*.

`main`'s cost is the mirror image: 23 kB of design.md and 13 kB of plan.md with no
mechanical link to tests, decisions numbered by hand, and an issue set that exists only in
someone's working tree. Its `docs/issues/*.md` are in fact the closest thing to openspec's
`tasks.md` and are *better* — each cites its authority — but they are untracked and would
vanish with the working tree.

---

## Divergence risk

**Who is ahead.** From the common ancestor `8b6848d`: `main` +2 commits, openspec +6
(`git rev-list --left-right --count main...origin/openspec` → `2 6`). Commit count is
misleading; by content:

- *Planning*: `main` is substantially ahead — ~48 kB of planning prose plus a 223-line
  method review and a 158-line validation write-up, against ~33 kB of openspec artefacts of
  which a third is duplicated between the live and archived specs.
- *Implementation*: comparable in volume, different in reach. `main`: ~1,850 lines of
  `src/mdins` (units, provenance, trajectory, spectral, ir with HDF5, cli) + ~3,180 lines
  of tests, milestones M0–M2 complete, D–G stubbed. openspec: ~975 lines of `src/md_ins`
  (trajectory, correlation, benchmark) + ~480 lines of tests, Stage 1 complete, Stages 2–3
  not started.

**The branches also disagree on the package name** (`mdins` vs `md_ins`, distributions
`mdins` vs `md-ins`), on the Python floor (≥3.11 vs ≥3.10), on the time unit (ps vs fs), on
the dependency set (`h5py` vs `phonopy`+`matplotlib`), and on the licence. There is no
cheap merge here; adopting one means porting from the other.

**What `main` would lose if it were abandoned** (i.e. things worth carrying to openspec):

1. `docs/method-review.md` in its entirety — in particular §2's exponential resummation,
   §3's eight Paper-A gaps, §3.7's factor-of-ten TOSCA resolution disagreement between
   OCLIMAX and AbINS, §4.4's sign error in Paper B Eq. 19, §4.7's thermostat constraint
   (NVE or Andersen-massive with τ ≈ 10 ps; "τ = 1 ps corrupts the spectrum"). None of this
   exists on openspec, and §4.7 in particular is an input-validation requirement openspec's
   spec should have.
2. design.md §5's ambiguity register A1–A8 — a project-level object OpenSpec has no
   artefact type for, but which could live in `config.yaml`'s `context:` or as a
   capability-level spec section.
3. D6 + the equipartition sum rule, and D5's PSD argument. Both are load-bearing for
   openspec's own Stage 3 and both are currently decided the other way there.
4. The derived-tolerance methodology and the negative control in `docs/validation.md`,
   which is the difference between "the peaks are in about the right place" and "the
   discrepancy is a fraction of a measured uncertainty and a 5% error would be caught".
5. The unit-boundary rule (plan.md §2, meV at every function boundary) and the memory rule
   ("The output frequency grid is an input to stage B, not a post-hoc downsample… Both
   source papers hit this wall; the fix is an API decision, not an optimisation").

**What openspec would lose if it were abandoned** (i.e. things worth carrying to `main`):

1. The document *structure* itself — specifically a living requirements spec, which `main`
   has no equivalent of and which would give its issue files an authority to cite that
   outlives them.
2. Stage 1 Decision 5's ASE→Phonopy→Euphonic force-constant bridge, written down as a
   decision with a rationale.
3. The kinetic-temperature validation requirement (`T_traj` vs nominal `T_MD`, 20% warning)
   — a real check `main`'s plan does not contain.
4. Stage 3 Decision 2's regularisation of `α_i(ω)` when `Tr(B)` → 0.
5. Stage 2's argument for `abinslib` as the right oracle for the isotropic/almost-isotropic
   maths — worth weighing against `main`'s D4, which declines AbINS on dependency-weight
   grounds and then relies on a closed binary (OCLIMAX) and an unmaintained script
   (MolDyINS) instead.
6. Stage 2 Decision 5's decoupling of broadening from the spectrum calculation.

**Immediate risk items, independent of which branch wins.**

- openspec's ∫g dω = 1 normalisation is baked into a shipped, archived spec *and* the
  implementation. If Stage 2's absolute intensity scale is ever to mean anything, that
  requirement has to be amended, and OpenSpec's change process makes that a visible,
  reviewable act — which is the structure working as intended.
- openspec's spec says 5%, its test asserts 15% and regression baselines. Someone should
  reconcile them.
- openspec's LICENSE is still BSD-3 while its plan names GPL-3.0 dependencies (`abinslib`)
  as core to validation. `main` resolved this deliberately as D1.
- `main`'s `docs/issues/*.md` are untracked. That is the single highest-value thing on
  `main` that is not actually on `main`.

---

## Appendix — what was compared

**Branches**

| Branch | Tip | Date | Subject |
|---|---|---|---|
| `main` (local, checked out, dirty) | `f1a12f9` | 2026-09-24 | Add second stage of implementation |
| `origin/openspec` | `7fce47d` | 2026-09-24 | improved pdos demo plot; broader |

Merge base `8b6848d` (2026-09-23, "Initial commit") contains `LICENSE` only.
`git rev-list --left-right --count main...origin/openspec` → **main +2, openspec +6**.
Neither is an ancestor of the other; there is no shared planning content to merge.

Other branches present but not compared: `origin/main` (`e66673f`, a different lineage
carrying MDANSE scripts, papers and SrTiO₃ data), `origin/proposal` (`0f09b98`, openspec's
ancestor), `origin/elephant_goldfish` (`e8e8a88`).

**Read on `main`** (all tracked and clean except where noted): `proposal.md`; `design.md`;
`plan.md`; `docs/method-review.md`; `docs/validation.md`; `README.md`. Also inspected,
**untracked**: `docs/issues/01-project-skeleton.md` … `17-euphonic-comparison.md` (17
files) — plus `src/mdins/*.py` and `tests/` line counts and `src/mdins/scattering.py` for
stub status.

**Read on `origin/openspec`** (via `git show`): `openspec/config.yaml`;
`openspec/specs/trajectory-pdos/spec.md`;
`openspec/changes/archive/2026-09-24-trajectory-pdos/{proposal,design,tasks}.md` and its
`specs/trajectory-pdos/spec.md`; `openspec/changes/ins-isotropic-sim/{proposal,design}.md`;
`openspec/changes/ins-tensor-anisotropic/{proposal,design}.md`; the three `.openspec.yaml`
files; `proposal-prompt.md`; `README.md`; `pyproject.toml`; `LICENSE` (first lines).
