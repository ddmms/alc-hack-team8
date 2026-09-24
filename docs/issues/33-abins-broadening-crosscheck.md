# Cross-check broadening against AbINS

## Context

`design.md` D4 keeps AbINS behind an optional extra rather than in the core: Mantid is a
very heavy dependency for the amount of it we would use, and it is awkward to install
outside conda. But D4 is equally clear that AbINS's calibrated instrument parameters are
the most valuable thing we are declining to depend on, and that the optional extra exists
specifically **to cross-check our broadening against theirs, not to supply it**.

This is also the second half of the TOSCA resolution question. `13-tosca-resolution-constants.md`
chooses between two published constants; this issue checks the resulting lineshape against
the implementation that uses the other one.

## Scope

- A `scripts/` harness comparing our broadened lineshape against AbINS's, for a delta input
  and for a representative spectrum, on TOSCA and on at least one other shared instrument.
- Comparison of the resolution *widths* AbINS computes against ours, energy by energy,
  separately from the comparison of the broadened result — a width discrepancy and a
  convolution discrepancy have different causes and should not be conflated.
- A `reference`-marked test that runs the width comparison when AbINS is importable and
  skips otherwise.
- A documented conda-based setup, since `abins` cannot be a pip extra
  (`plan.md` §1).
- The result fed back into `13-tosca-resolution-constants.md`'s provenance block.

## Out of scope

- Depending on AbINS at runtime. D4 settled that.
- Comparing full `S(Q,ω)` — AbINS has no MD input path (method-review §6), which is the
  gap this project exists to fill, so there is nothing to compare end to end.

## Acceptance criteria

- [ ] Widths are compared at a spread of energies covering the range where the quadratic
      term matters, so the comparison can actually see the factor-of-ten disagreement
      rather than sampling only where it is negligible.
- [ ] If the widths differ by the known factor of ten in the quadratic coefficient, the
      write-up says so and connects it to method-review §3.7 rather than reporting an
      unexplained mismatch.
- [ ] With the same width model forced on both sides, the broadened lineshapes agree to
      within stated tolerance. This separates "we disagree about the constants" from "we
      disagree about the convolution", which is the whole value of the exercise.
- [ ] The setup instructions are followed successfully by someone other than the author,
      or the harness is marked as unverified.
- [ ] The test skips cleanly without AbINS and never appears in the default PR run.

## Depends on

- `14-resolution-broadening.md`
- `13-tosca-resolution-constants.md`

**Labels:** `milestone:M5`, `stage:F`, `validation`
