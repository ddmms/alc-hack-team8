# Resolve the TOSCA and VISION resolution constants and write their data files

## Context

This issue exists because of a genuine, unresolved contradiction in the literature, and it
is deliberately separated from the mechanics of the data format so that it cannot be
settled by accident.

method-review §3.7: Paper A′ Eq. 12 gives the TOSCA/VISION resolution as
`σ = 0.25 + 0.005·ΔE + 1e-7·ΔE²` in wavenumbers. Mantid AbINS uses the same `a` and `b`
but `c = 2.5` — a factor of ten larger. One of them is a typo. `design.md` §2[F] is
explicit that this "has to be resolved explicitly when the data file is written — not
settled by whichever source we happened to copy from", and `plan.md` §6 lists it as a
standing risk.

There is a second, quieter problem in the same review item: **VISION is not TOSCA**, and
VISION's own parameters are not published anywhere we found. Paper B's VISION handling
does not help — its `σ = 0.01·E` is admitted to have been chosen by eye to match the C–H
stretch width (method-review §4.8).

## Scope

- Evaluate both candidate `c` values against something external: the published TOSCA
  resolution curve, AbINS's own output via the optional extra, and if possible an
  instrument scientist's confirmation.
- Write `tosca.toml` with the chosen constants and a provenance block that records **both**
  published values, which was chosen, and on what evidence.
- Write `vision.toml`, or document explicitly why it is not shipped. Shipping TOSCA's
  constants relabelled as VISION would be the worst outcome: it looks like support and is
  not.
- Add the resolution entry to the ambiguity resolution log
  (`36-ambiguity-resolution-log.md`).
- A user-facing note in the docs that the TOSCA resolution carries a known literature
  disagreement, with the magnitude of its effect.

## Out of scope

- The broadening machinery that consumes these constants —
  `14-resolution-broadening.md`.
- Other instruments — `12-instrument-data-format.md`.

## Acceptance criteria

- [ ] A test evaluates the shipped TOSCA resolution at several energies and compares
      against the published curve; it must be able to *distinguish* the two candidate `c`
      values. If the quadratic term contributes less than the comparison's own uncertainty
      over the instrument's range, say so in the write-up rather than claiming the test
      resolved anything — that outcome is informative and must not be dressed up.
- [ ] The chosen value and the rejected one both appear in the data file's provenance, so
      a future reader can re-litigate the choice without rediscovering the disagreement.
- [ ] If AbINS is installed, a cross-check test compares our broadened lineshape against
      theirs and reports the discrepancy; it is marked `reference` and skipped otherwise.
- [ ] `docs/` states the disagreement in one paragraph, with the effect on a
      representative spectrum quantified rather than described.

## Depends on

- `12-instrument-data-format.md`
- `14-resolution-broadening.md`

## Open questions / risks

- **This may not be resolvable from the literature alone.** Both published values are in
  print and neither source shows its derivation. The acceptable outcomes are: resolved
  with evidence; or resolved provisionally with the alternative shipped as a commented
  option and the uncertainty documented. The unacceptable outcome is a data file that
  states one value as fact with no trace of the other.
- The quadratic term is small over TOSCA's usual range, so the factor of ten may be
  practically irrelevant at low energy transfer and matter only above a few hundred meV.
  Quantify where it starts to matter; that bound is more useful to a user than a verdict.

**Labels:** `milestone:M3`, `stage:F`, `physics`, `open-question`
