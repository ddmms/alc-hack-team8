# Build the OCLIMAX oracle harness for Method 1

## Context

`plan.md` §5 Layer 5 and `design.md` §4 Tier 3. OCLIMAX is closed source and distributed
as a Docker image, so it cannot be a base — but it is the reference implementation of
Paper A's method and therefore the numerical oracle for Method 1.

It carries disproportionate weight because of D7. Ruling out the exponential resummation
removed the internal cross-check on the multiphonon chain, and `design.md` §2[E] says
verification now rests entirely on the Tier 1 analytic tests and external oracles. This is
the external oracle for the isotropic method.

It is also the only way to settle ambiguity A5: Paper A never states what OCLIMAX does with
anisotropy given only `tr(B)` from MD, and `design.md` §5 says to assume an isotropic
fall-back and **confirm against the OCLIMAX oracle**. That confirmation happens here or
nowhere.

## Scope

- A scripted harness under `scripts/` (not `tests/`) that runs OCLIMAX in Docker on a
  chosen system and captures its output.
- Conversion of our IR into whatever OCLIMAX needs, or generation of the same input for
  both codes. Note method-review §3.3: the `tclimax` trajectory format is undocumented and
  its converter ships only inside the Docker image, so plan for the input path to be the
  hard part.
- Comparison of our Method 1 spectrum against OCLIMAX's, order-resolved where OCLIMAX
  exposes orders, on a common grid and in a common normalisation.
- A written-up result: where we agree, where we disagree, and by how much. Results are
  written up, not asserted — this cannot run in CI.
- Explicit checks on A5 (isotropic fall-back) and, if the orders are separable, on A7 (the
  multiphonon prefactor).

## Out of scope

- Method 2 — `32-moldyins-benchmark.md`.
- Any CI integration. `plan.md` §5 is explicit that the Docker requirement puts this
  outside CI.
- Reverse-engineering the OCLIMAX binary.

## Acceptance criteria

- [ ] The harness runs end to end on at least one system and produces both spectra on a
      common grid, reproducibly, from a single documented command.
- [ ] The comparison is quantitative: integrated intensity, peak positions, and per-order
      intensities where available, with the disagreement reported as a number.
- [ ] A5 is addressed explicitly. Either the comparison is consistent with an isotropic
      fall-back — state the evidence and the sensitivity — or it is not, in which case the
      finding is recorded and `22-isotropic-powder-average.md`'s assumption is revisited.
- [ ] Every disagreement is attributed, or explicitly recorded as unattributed. "Close
      enough" without a candidate explanation is not an acceptable write-up for a
      comparison whose entire purpose is to find our bugs.
- [ ] The write-up states what the comparison *cannot* distinguish, given that OCLIMAX's
      internals are not inspectable.

## Depends on

- `24-method-1-assembly.md`

## Open questions / risks

- **The input path may be the blocker.** The undocumented `tclimax` format with a converter
  only inside the image means the practical route may be to drive OCLIMAX's own converter
  from within the container. Timebox this; if the input cannot be matched, the comparison
  is on a system both codes can construct independently, which weakens it but does not
  make it worthless.
- Unbounded in effort by `plan.md` §4's own assessment of M5. Nothing may depend on this.

**Labels:** `milestone:M5`, `validation`, `open-question`
