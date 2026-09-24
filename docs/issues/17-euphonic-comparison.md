# Compare LJ pDOS against Euphonic and write up the harmonic validation

## Context

The headline validation from the proposal, and the harder half of M2's definition of done
(`plan.md` §7). The reference side exists by now; this issue supplies the MD side, the
comparison, the tolerance, and the write-up.

The tolerance is the part that needs care. `design.md` §4 is explicit that it is **derived
from the inter-segment spread of the Welch estimate** (D5) rather than chosen by hand: a
mismatch only counts as a failure if it exceeds the sampling noise of the MD side.
Anything else is a number tuned until the test went green, which proves nothing. The same
section names the second trap — stay in the harmonic limit. LJ is anharmonic, and at high
temperature the two will legitimately disagree and the test tells us nothing.

A single-species crystal cannot check the per-species projection, which is the feature the
entire method rests on. An ordered Ar/Kr crystal fixes that: the potential is
species-blind, so mass is the only asymmetry, and any per-species difference we report has
a known origin.

## Scope

- An LJ argon crystal fixture (a supercell such as 3×3×3), NVE MD at low temperature with
  a fixed seed, run through stages A–C.
- Comparison against the Euphonic reference on the commensurate q-grid, in normalised
  lineshape space.
- Tolerances propagated from the Welch inter-segment spread through each compared
  quantity, expressed in σ rather than in percent.
- Agreement criteria per `design.md` §4: integrated intensity per species, peak positions
  (first and second spectral moments are a more robust proxy than peak picking), and a
  normalised spectral overlap.
- A negative control asserting the same threshold *rejects* a deliberately wrong spectrum.
- A second fixture: ordered Ar/Kr, checking the per-species projection.
- A written-up validation document with figures — both spectra overlaid, residuals, and the
  mixed-species per-species panels — plus the script that regenerates them.
- Both marked `slow` and `reference`.

## Out of scope

- The reference computation — `16-euphonic-reference-harness.md`.
- Anything past stage C. This compares pDOS, not spectra.
- Cross-code benchmarks — M5.

## Acceptance criteria

- [ ] Every stated tolerance is computed from the inter-segment spread at runtime, not
      written as a literal. A reviewer should be able to see where each number came from.
- [ ] The observed discrepancy in each compared moment is reported in units of that σ, so
      the test result is a measurement rather than a pass/fail.
- [ ] A negative control injects a 5% frequency error and asserts the same threshold
      rejects it. Without this, agreement could be an artefact of a tolerance wide enough
      to admit anything, and the whole milestone would be unsupported.
- [ ] The MD temperature is low enough that anharmonicity is below the tolerance, and this
      is *demonstrated* — run at two temperatures and show the discrepancy grows with the
      higher one — not asserted in a comment.
- [ ] The Euphonic reference uses the commensurate q-grid; a test that deliberately uses a
      dense mesh instead shows a visibly worse match, which documents why the constraint
      exists.
- [ ] Per-species moments for the Ar/Kr crystal agree within the propagated tolerance,
      species by species. The single-species case cannot fail in a way that catches a
      mass mix-up in the projection; this one can.
- [ ] The two traps from `design.md` §4 appear in the test docstrings, so the next person
      to see a red run reads the explanation before debugging the wrong thing.
- [ ] The write-up states the tolerance derivation, the observed discrepancies, and the
      negative control result, and is regenerable from a checked-in script.

## Depends on

- `15-einstein-crystal.md`
- `16-euphonic-reference-harness.md`

## Open questions / risks

- Highest-risk item in the plan (`plan.md` §4, M2 is marked "M, highest risk"). If the
  comparison disagrees beyond tolerance after the Einstein crystal passes, the remaining
  candidates are q-sampling, anharmonicity, and the finite-displacement reference itself —
  in that order of likelihood. Budget for that investigation rather than assuming it is
  our estimator.

**Labels:** `milestone:M2`, `testing`, `physics`, `docs`
