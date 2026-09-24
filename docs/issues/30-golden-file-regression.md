# Add the golden-file regression layer

## Context

`plan.md` §5 Layer 4. The analytic layers prove correctness against known answers; this
layer catches the different failure of a refactor that quietly changes a number nobody was
asserting. Golden HDF5 outputs for small checked-in fixtures, compared with tolerance.

The governance matters more than the mechanism. `plan.md` is explicit that regeneration
must be a reviewed, deliberate commit, because a golden file silently updated to match a
bug is worse than no test at all — it converts a regression into a documented expectation.

## Scope

- Small checked-in input fixtures, under ~1 MB total, covering at least: a synthetic
  multi-species trajectory, an IR produced from it, and a spectrum from each method.
- Golden artefacts stored as HDF5 alongside a manifest recording the code version, the
  parameters, and the date each was generated.
- A comparison helper with per-field tolerances, reporting the worst-disagreeing field and
  its magnitude rather than a bare assertion failure.
- `scripts/regenerate_fixtures.py`, which regenerates everything and prints a diff summary
  against the existing goldens before writing.
- A CONTRIBUTING note, or a docstring at the top of the script, stating the review
  requirement.

## Out of scope

- External fixtures such as the MolDyINS trajectory — `32-moldyins-benchmark.md`.
- Anything requiring Docker or an optional dependency; this layer is fast and runs on
  every PR.

## Acceptance criteria

- [ ] The regression tests run in the fast suite and complete in seconds.
- [ ] A deliberate 1% perturbation to any stage's output fails at least one golden
      comparison, and the failure message names the artefact and the field. Verify by
      perturbing and observing, not by inspection.
- [ ] Tolerances are per-field and justified: floating-point reproducibility across
      platforms and NumPy versions is not exact, so the tolerance has to be wide enough to
      survive that and narrow enough to catch a real change. State which of the two set
      each number.
- [ ] `regenerate_fixtures.py` is idempotent — running it twice produces byte-identical
      files — and prints what changed before overwriting.
- [ ] Regenerating requires an explicit flag; running the script with no arguments reports
      the diff and exits without writing.
- [ ] Total fixture size stays under 1 MB, asserted by a test so it cannot creep.

## Depends on

- `24-method-1-assembly.md`

**Labels:** `milestone:M5`, `testing`, `infrastructure`
