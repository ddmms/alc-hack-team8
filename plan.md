# Implementation and test plan

Concrete work breakdown for the architecture in [`design.md`](design.md). That document
deliberately stopped short of module layout and signatures; this one supplies them, plus
sequencing and the test strategy.

---

## 1. Package and tooling

Proposed distribution name **`mdins`**, import name `mdins`. Check PyPI before committing
to it; the name is cheap to change now and expensive later.

```
pyproject.toml          hatchling; requires-python >=3.11
src/mdins/
  units.py              canonical units; plain floats, no pint (see design.md §3)
  provenance.py         provenance record construction
  trajectory.py         [A] ASE loading, COM removal, sampling validation
  spectral.py           [B] Welch and VACF estimators
  ir.py                 [C] VelocitySpectralDensity dataclass, HDF5 I/O
  displacement.py       [D] u²(ω), Debye–Waller factor
  scattering/
    convolution.py      multiphonon orders, prefactors
    powder.py           isotropic and almost-isotropic averaging
    isotropic.py        Method 1 (Cheng)
    anisotropic.py      Method 2 (Harrelson)
  instrument/
    kinematics.py       Q(ω) for indirect and direct geometry
    resolution.py       energy-dependent broadening
    data/*.toml         TOSCA, VISION, MARI, …
  spectrum.py           [G] output container
  cli.py                thin argparse/typer wrapper
tests/                  mirrors the above, plus tests/reference/ and tests/data/
```

Runtime dependencies: `numpy`, `scipy`, `ase`, `h5py`.
Optional extras: `euphonic` (M2 validation only, so not a runtime requirement),
`mdanalysis` (dev-only, fixture conversion). `abins` is conda-only and cannot be an
extra; see the note in `pyproject.toml`. `pint` is not used — see design.md §3.

Dev tooling: `uv` for environments and locking, `ruff` for lint and format, `mypy`
(strict on `ir.py`, `units.py`, `spectral.py` at minimum), `pytest`.

CI on GitHub Actions: Python 3.11–3.13, ruff + mypy + fast tests on every PR, slow and
reference tests nightly.

---

## 2. Two conventions to fix before any code is written

These cause the failure mode this domain is notorious for — a plausible-looking spectrum
that is wrong by 2π or by a factor of mass.

**Frequencies are energies, in meV, everywhere a value crosses a function boundary.**
Angular frequency `ω` in rad/ps may exist inside a function body and must never appear in
a signature, a dataclass field, or a file. This single rule eliminates most of the 2π
bugs by construction. `units.py` provides the conversions and nothing else may.

**The output frequency grid is an input to stage B, not a post-hoc downsample.**
Passing the target grid into the estimator lets it rebin per segment, so the
full-resolution per-atom array never materialises. Without this, a 3000-atom system at
25000 bins is 5.4 GB before anything useful has happened. Both source papers hit this
wall; the fix is an API decision, not an optimisation.

---

## 3. The IR is the contract — write it first

Everything else is either upstream or downstream of this object, so freezing it on day one
is what lets work proceed in parallel.

```python
@dataclass(frozen=True)
class VelocitySpectralDensity:
    frequencies: NDArray[np.float64]  # (n_freq,)              meV
    density: NDArray[np.float64]  # (n_entity, n_freq, 3, 3)
    symbols: list[str]  # (n_entity,)
    masses: NDArray[np.float64]  # (n_entity,)            amu
    temperature_md: float  # K
    metadata: SpectralMetadata  # estimator, window, segments, normalisation,
    # ensemble, provenance

    def pdos(self) -> NDArray[np.float64]: ...  # trace/3 → (n_entity, n_freq)
    def total_pdos(self) -> NDArray[np.float64]: ...
    def check_sum_rule(self, rtol: float = 1e-3) -> None: ...
    def to_hdf5(self, path) -> None: ...
    @classmethod
    def from_hdf5(cls, path) -> "VelocitySpectralDensity": ...
```

Trailing `(3, 3)` axes are stored in full rather than as six symmetric components: 1.5×
the memory, but `np.einsum('...ij,i,j->...', density, q, q)` works directly and the
symmetric packing is an easy place to introduce an off-by-a-factor-of-two. Pack to six
components on HDF5 write only.

Public entry points, so callers can be written against them before internals exist:

```python
sd = velocity_spectral_density(traj, frequencies=grid, estimator="welch")
spec = isotropic_spectrum(sd, instrument="TOSCA", temperature=10.0, max_order=10)
spec = anisotropic_spectrum(sd, instrument="TOSCA", temperature=10.0, max_order=10)
```

CLI mirrors these: `mdins pdos traj.extxyz -o sd.h5` and
`mdins spectrum sd.h5 --instrument tosca --temperature 10`.

---

## 4. Milestones

Each produces something checkable. Sizes are relative, not calendar estimates — I don't
know the timebox.

| | Milestone | Deliverable | Size |
|---|---|---|---|
| **M0** | Scaffolding | repo, CI, `units.py`, `ir.py`, stub entry points | S |
| **M1** | pDOS from MD | stages A–C, analytic tests passing | M |
| **M2** | Harmonic validation | Euphonic Tier 2 comparison agrees | M, highest risk |
| **M3** | Method 1 | stages D, E-isotropic, F, G — end-to-end spectrum | L |
| **M4** | Method 2 | anisotropic path | M |
| **M5** | Benchmarks | OCLIMAX, MolDyINS, published spectra | M, unbounded |

**M0–M2 is the honest first cut.** It delivers the capability the proposal lists first,
it is the part with an unambiguous right answer to check against, and it is the foundation
under everything else. M3 is a reasonable stretch. Treating M5 as in-scope for a short
effort is how this ends up with a spectrum nobody has verified.

### Parallelisation after M0

Three tracks that barely touch:

- **A–C** (trajectory → IR). The critical path.
- **F** (instrument). Completely self-contained — takes an energy grid, returns `Q(ω)` and
  a broadening kernel. Testable against published TOSCA curves with no MD involved.
- **Validation harness** (Euphonic reference, ASE LJ setup, fixture generation). Can be
  built against the frozen IR before the IR has a real implementation behind it.

D–E depend on C existing but not on C being correct, so they can start once the dataclass
is frozen and fed synthetic data.

---

## 5. Test strategy

Five layers, decreasing in speed and increasing in what they actually prove.

### Layer 1 — analytic, on synthetic velocities (fast, always run)

Generate velocity arrays directly. No MD, no files, exact expected answers, millisecond
runtimes. This layer catches essentially every normalisation, unit and 2π error.

- Single sinusoid → delta at known energy, amplitude fixed by the sum rule
- Sum of sinusoids at different amplitudes → correct relative weights
- Damped oscillator → Lorentzian of known width
- White noise → flat spectrum; Parseval holds
- Anisotropic signal (different amplitude per Cartesian axis) → correct diagonal tensor
- Correlated x/y motion → known off-diagonal term, and known *sign*
- Circular motion → nonzero quadrature part, which is where ambiguity A2 becomes visible

Invariants asserted as properties across randomised inputs:

- `∫ tr P dω == 3k_BT/m` (equipartition, D6's anchor)
- every `P_i(ω)` has non-negative eigenvalues (D5's justification — this test is the
  reason Welch is the default, so it must exist)
- Welch and VACF agree within the Welch inter-segment spread (D5's cross-check)
- HDF5 round-trip is exact, including metadata

### Layer 2 — analytic physics (fast)

- **Multiphonon against closed form.** For an isotropic 3D harmonic oscillator the exact
  incoherent result is `S(Q,ω) = e^{-2W} Σ_n [(Q²u²)ⁿ/n!] δ(ω − nω₀)`. This is a direct
  analytic oracle for the convolution chain, the Debye–Waller factor, *and* the order
  prefactor — which is exactly the question D7 left open as ambiguity A7. It is the single
  most valuable test in the suite and should be written before the convolution code.
- **Instrument kinematics** against hand-computed `Q(ω)` at a few energies, and against the
  published TOSCA trajectory curve.
- **Powder averaging**: an isotropic input must give identical results through the
  isotropic and almost-isotropic paths.

### Layer 3 — end-to-end through real MD (slow, marked)

- **Einstein crystal.** ASE atoms on harmonic springs, NVE, fixed seed. The pDOS is a
  delta at a known frequency. Proves the whole A–C chain against an exact answer while
  exercising real trajectory I/O.
- **LJ crystal, Tier 2.** The proposal's headline validation:

  ```
  ASE LennardJones ──┬── finite displacement → Euphonic ForceConstants → calculate_pdos
                     └── NVE MD, low T      → mdins → pdos()
  ```

  Two traps that will produce a spurious mismatch, both of which belong in the test's
  docstring: the Euphonic reference must use the q-grid commensurate with the MD
  supercell, and the MD must be cold enough to stay harmonic. Comparison is on normalised
  lineshapes (D6), with tolerance derived from the Welch inter-segment spread (D5).

Determinism: fixed RNG seeds, NVE rather than Langevin where possible. Where a thermostat
is unavoidable, assert statistically with a stated confidence rather than pinning numbers.

### Layer 4 — regression (fast)

Golden HDF5 outputs for small checked-in fixtures, compared with tolerance. A
`scripts/regenerate_fixtures.py` regenerates them, and regeneration must be a reviewed,
deliberate commit — a golden file silently updated to match a bug is worse than no test.

Fixtures stay under ~1 MB. The MolDyINS 25 K GROMACS trajectory is converted once offline
(MDAnalysis, dev-only) to extxyz and checked in, since ASE cannot read `.trr` (D2).

### Layer 5 — cross-code benchmarks (not CI)

Scripted, run on demand, results written up rather than asserted:

- **OCLIMAX** as a numerical oracle for Method 1. Docker-only binary, so the harness needs
  Docker and cannot run in CI.
- **MolDyINS** for Method 2. Needs repair first — it fails on modern NumPy (`np.float`).
  Effort here is genuinely unknown and should not block anything.
- **AbINS** broadening against ours, via the optional extra (D4).
- **Published spectra** — ice Ih, MgH₂ (strong multiphonon, so it exercises the
  convolution chain hardest), P3HT. Figures in docs, not assertions.

---

## 6. Risks, and what to do about them

| Risk | Mitigation |
|---|---|
| Silent factor errors (2π, mass, ħ) | Layer 1 is written *first*, and the meV-only boundary rule in §2 |
| Tier 2 mismatch is ambiguous — is it us, anharmonicity, or q-sampling? | Einstein crystal first: it has an exact answer, so it isolates our bugs from physics before LJ introduces both |
| Memory on realistic systems | Target grid as estimator input (§2); entity grouping; chunking over atoms |
| Multiphonon prefactor genuinely undetermined | Layer 2's closed-form oscillator decides it, rather than a convention argument |
| TOSCA resolution constants disagree ×10 | Resolve explicitly when writing the data file; document the choice and the disagreement |
| MolDyINS repair is open-ended | Layer 5 only, blocks nothing |
| Scope creep into M5 | M0–M2 is the committed cut |

---

## 7. Definition of done, per milestone

- **M0** — CI green, IR round-trips through HDF5, entry points importable and raising
  `NotImplementedError`.
- **M1** — all of Layer 1 passes; pDOS from a real ASE trajectory; sum rule holds.
- **M2** — Einstein crystal exact; LJ crystal within stated tolerance of Euphonic, with
  the tolerance justified rather than tuned.
- **M3** — TOSCA spectrum for a published system, order-resolved, Layer 2's analytic
  multiphonon test passing.
- **M4** — anisotropic path reproduces the isotropic one on isotropic input, and differs
  in the expected direction on anisotropic input.
- **M5** — written-up comparison against at least one external code, including where it
  disagrees.

Throughout: every ambiguity in design.md §5 that gets resolved is recorded with the
evidence that resolved it. Per the proposal, that record is a deliverable, not a byproduct.
