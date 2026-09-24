# Maintain the ambiguity resolution log

## Context

`plan.md` §7 closes with a requirement that runs through every other issue: "every
ambiguity in design.md §5 that gets resolved is recorded with the evidence that resolved
it. Per the proposal, that record is a deliverable, not a byproduct."

That is not documentation hygiene. The eight ambiguities in `design.md` §5 are places where
two published papers are silent, contradictory, or wrong, and where this project therefore
had to choose. Anyone comparing our output against OCLIMAX, MolDyINS or an experiment needs
to know which choices we made and on what basis, or a disagreement is uninterpretable. The
review in `docs/method-review.md` did this for the papers; this does it for us.

## Scope

- A document listing A1–A8 from `design.md` §5, each with: the ambiguity as stated, the
  working resolution adopted, the **evidence** that supports it, the issue and commit where
  it was resolved, and its current status (resolved, provisional, or open).
- The same treatment for the non-§5 discrepancies that also needed a decision: the TOSCA
  resolution constants (method-review §3.7), the implicit quantum correction factor
  (method-review §3.6), and detailed balance through the multiphonon convolution
  (method-review §3.5).
- A convention that an issue resolving an ambiguity updates this file in the same PR, noted
  in CONTRIBUTING.
- Entries that remain open stay in the list marked open, with what evidence would close
  them. An ambiguity that was never settled is as important to record as one that was.

## Out of scope

- Re-reviewing the papers. `docs/method-review.md` already does that and is the input to
  this, not a duplicate of it.
- General user documentation — `35-user-documentation.md`, which should link here for the
  parameters that are choices.

## Acceptance criteria

- [ ] All eight of A1–A8 appear, including those whose status at the end is "open" or
      "provisional". A log that only records successes would misrepresent the project's
      confidence in its own output.
- [ ] Each resolved entry cites something checkable — a test name, an equation, an external
      comparison — not a preference. "We chose `1/n!`" is not a resolution; "the
      single-oscillator closed form in `tests/test_multiphonon.py::test_order_ratios`
      selects `1/n!` and rejects the aCLIMAX prefactors at the 3σ level" is.
- [ ] A test or CI check asserts that every ambiguity identifier referenced in a source
      docstring (`A1`…`A8`) has an entry here, so a code comment can never point at a
      log entry that does not exist.
- [ ] Entries whose resolution is provisional say what would change it.
- [ ] The log is linked from the README and from the user documentation, since it is one of
      the project's stated deliverables rather than an internal note.

## Depends on

- `03-intermediate-representation.md` (A1 is settled by D6's sum rule, which is the first
  ambiguity to close)

## Open questions / risks

- The entries that are likeliest to stay open are A7 (multiphonon prefactor, if the
  analytic oracle cannot distinguish the conventions) and A8 (the rank-4 contraction,
  whose error has no internal bound at all under D7). Plan for the log to end with genuine
  open items and resist the pressure to write them up as resolved.

**Labels:** `docs`, `physics`, `open-question`
