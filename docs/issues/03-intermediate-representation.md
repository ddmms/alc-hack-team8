# Freeze the VelocitySpectralDensity intermediate representation

## Context

`plan.md` §3 says to write this first, and the reason is scheduling as much as design:
everything in the pipeline is either upstream or downstream of this object, so freezing
it on day one is what allows the three tracks in `plan.md` §4 to run in parallel. Stages
D–G can be written and tested against synthetic IR instances long before stage B produces
a correct one.

The structural claim behind it is `design.md` §1: Method 1 and Method 2 are not two
pipelines. The scalar atom-projected DOS is the trace of the tensor, `g_i = tr P_i / 3`,
so the tensor is computed once and Method 1 is a consumer that takes the trace at stage
D. That is what lets the two methods be compared on byte-identical input, and it is the
single most important thing this dataclass encodes.

Storage is **absolute spectral density**, not a normalised lineshape (D6). The
equipartition sum rule then applies directly to the stored array, which is the
normalisation anchor both source papers lack (method-review §3.1).

## Scope

- `src/mdins/provenance.py`: a `Provenance` record capturing source identity, package
  version, and the ordered list of pre-processing steps applied. Serialisable to and from
  plain JSON-compatible types.
- `src/mdins/ir.py`: `SpectralMetadata` (estimator, window, segment count and length,
  normalisation identifier, ensemble, energy resolution, provenance) and the frozen
  `VelocitySpectralDensity` dataclass with the fields in `plan.md` §3.
- The trailing `(3, 3)` axes stored in full, not packed to six symmetric components:
  1.5× the memory, but `np.einsum('...ij,i,j->...', density, q, q)` works directly and
  symmetric packing is an easy place to introduce a factor of two.
- Methods: `pdos()` (trace/3), `total_pdos()`, `mean_square_velocity()`,
  `check_sum_rule(rtol)`, and a `normalised_lineshape()` helper.
- A `NORMALISATION` constant recording the one convention this package writes, stored in
  every artefact so the convention travels with the data.
- Bin-edge and bin-width helpers, so the integration rule the estimator conserves and the
  rule the sum-rule check applies are literally the same function.

## Out of scope

- HDF5 I/O — `04-hdf5-serialisation.md`.
- Producing an instance from a trajectory — `07-welch-estimator.md`.
- Asserting the invariants across estimators — `09-analytic-test-suite.md`. This issue
  provides `check_sum_rule`; it does not prove any estimator satisfies it.

## Acceptance criteria

- [ ] The dataclass is frozen, and mutating a field raises. Arrays are the shapes
      documented in `plan.md` §3; a constructor validation test feeds mismatched
      `symbols`/`masses`/`density` lengths and asserts a clear error naming the offending
      field, rather than a `ValueError` from broadcasting three stages later.
- [ ] `pdos()` equals `trace(density)/3` on a randomised tensor, and `total_pdos()` is the
      sum over entities.
- [ ] `check_sum_rule` passes on a synthetic density constructed to satisfy
      `∫ tr P dE = 3 k_B T / m` and fails with a message naming the worst-offending atom
      and the observed ratio when one atom's density is scaled by 1.05. A check that only
      says "sum rule violated" is not enough to debug a 3000-atom system.
- [ ] `normalised_lineshape()` integrates to one per entity, and the round trip through
      the stored absolute scale recovers the original array.
- [ ] `Provenance` round-trips through JSON, including an empty step list.
- [ ] Every array field has its unit stated in the docstring. A reader should not have to
      infer that `density` is Å²·ps⁻²·meV⁻¹ from the sum rule.
- [ ] mypy strict passes on both modules.

## Depends on

- `02-units-and-mev-boundary.md`

## Open questions / risks

- Grouping (per-species rather than per-atom) is deliberately a write-time option, not a
  field on the in-memory object (D3). If it later needs to be a first-class field, that
  is a breaking change to the contract everything else was written against — so resist it
  here and revisit only with evidence from a real system.

**Labels:** `milestone:M0`, `stage:C`, `architecture`
