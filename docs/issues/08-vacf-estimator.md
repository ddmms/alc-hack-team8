# Implement the VACF (Blackman–Tukey) estimator

## Context

The second estimator in `design.md` D5. It exists for two reasons that are both worth
more than the redundancy suggests.

It is the **paper-faithful** route: Paper A's Eqs. 5–6 are a normalised VACF transformed
to a DOS, so any discrepancy against OCLIMAX can be attributed to our method choice or
ruled out. And the VACF is a genuinely useful diagnostic in its own right — whether the
correlation has actually decayed within `τ_max` is visible in the VACF and invisible in
the spectrum it produces.

Most importantly, running the two estimators against each other is an independent check on
the stage of the pipeline most likely to be quietly wrong. `design.md` D5 makes their
agreement, within the Welch inter-segment spread, an explicit cross-check rather than a
nice-to-have.

## Scope

- A `estimator="vacf"` path through `velocity_spectral_density`, returning the same
  `VelocitySpectralDensity` type with `metadata.estimator` set accordingly.
- VACF computed by Wiener–Khinchin (FFT) rather than a direct double loop, with `τ_max`
  as the documented resolution knob corresponding to Welch's segment length.
- A lag window (Blackman–Harris or similar, configurable) applied before transforming,
  with the choice and its parameters recorded in metadata.
- **Full-lag, not symmetric-lag, correlation for the off-diagonal components.** `C_αβ(t)`
  is not even in `t` for α≠β, so a symmetric-lag implementation silently discards the
  quadrature part — turning ambiguity A2 into an accident rather than a decision
  (`design.md` §2[B]).
- The same absolute normalisation as the Welch path, so the two are directly comparable
  without a convention conversion.
- Exposure of the VACF itself as a diagnostic return, not just its transform.

## Out of scope

- Making this the default. It is not, and the reason (positive semi-definiteness) is D5.
- The cross-estimator agreement test — `09-analytic-test-suite.md`.

## Acceptance criteria

- [ ] On synthetic velocities the two estimators agree within the Welch inter-segment
      spread, component by component, including off-diagonals. A test that compares only
      the trace would pass while the off-diagonal sign convention was inverted.
- [ ] The VACF of a single undamped sinusoid is a cosine at the same frequency with
      `C(0) = ⟨|v|²⟩`, checked directly against the closed form before any transform.
- [ ] For a correlated x/y signal, `C_xy(t) != C_xy(-t)`, and a test asserts the
      implementation keeps both halves. Force the symmetric-lag shortcut in a test double
      and confirm it loses the quadrature part — this documents why the full-lag version
      is required.
- [ ] Positive semi-definiteness is checked and allowed to fail: the test asserts that a
      short `τ_max` with an aggressive window can produce a small negative eigenvalue, and
      that the IR's PSD assertion reports it rather than the estimator hiding it. This is
      the concrete evidence for D5's choice of default.
- [ ] `τ_max` longer than the trajectory raises; `τ_max` short enough that the correlation
      has not decayed produces a warning naming the residual correlation at `τ_max`.
- [ ] The sum rule holds on the VACF path to the same tolerance as on the Welch path.

## Depends on

- `07-welch-estimator.md`

**Labels:** `milestone:M1`, `stage:B`, `physics`
