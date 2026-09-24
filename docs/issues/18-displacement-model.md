# Implement the displacement model and the low-frequency cutoff

## Context

Stage D (`design.md` §2[D]): the conversion from a classical spectral density to quantum
harmonic displacement amplitudes,

```
u²_i,αβ(E) = ħ P_i,αβ(E) / (ω k_B T_MD)
```

This is Paper B Eq. 8, but `design.md` says to follow the discrete-transform bookkeeping in
its SI Eq. S13 rather than Eq. 8 as printed — method-review §4.1 records that Eq. 8 is
dimensionally inconsistent and the paper admits as much.

The `1/ω` factor diverges as ω→0, and this is a correctness issue neither paper addresses
(ambiguity A6). Diffusive and quasi-elastic contributions sit exactly where the divergence
is, so integrating through them does not produce a slightly wrong Debye–Waller factor, it
produces an arbitrary one. `design.md` §5 requires an explicit cutoff **surfaced as a
required parameter, not a hidden constant**, which is a deliberate contrast with both
source implementations.

## Scope

- `src/mdins/displacement.py` with `mean_square_displacement(sd, ...)` returning the
  `u²(E)` tensor per atom in the IR's layout.
- The low-frequency cutoff as an explicit argument with no silent default, applied by
  zeroing or excluding bins below it, and recorded in the output metadata.
- A diagnostic reporting how much spectral weight lies below the cutoff, so a user can see
  whether they have excluded a rounding error or a third of the spectrum.
- Method 1's entry into this stage: `tr(u²)/3`, taken here rather than earlier, so that
  both methods consume the same tensor.
- Documentation of which discrete-transform convention (SI Eq. S13) is implemented, in the
  module docstring, with the factor written out.

## Out of scope

- The Debye–Waller factor and the two-temperature question — `19-debye-waller-factor.md`.
- Multiphonon — `21-multiphonon-convolution.md`.
- Detecting diffusion. We exclude the low-frequency region; we do not attempt to model it.

## Acceptance criteria

- [ ] For a single harmonic oscillator, `∫ u²(E) dE` equals the analytic harmonic
      mean-square displacement `k_B T / (m ω₀²)` in the classical limit. This one test
      pins the entire ħ/k_B/ω bookkeeping, which is where Paper B's printed equation goes
      wrong.
- [ ] Calling without a cutoff raises, with a message explaining why there is no safe
      default. A default of zero would produce an infinite DWF on any real trajectory; a
      default of "something small" would hide the problem.
- [ ] The excluded-weight diagnostic is exact: on synthetic input with a known fraction of
      weight below the cutoff, the reported fraction matches.
- [ ] Changing the cutoff across a region with no spectral weight changes `∫ u² dE` by less
      than the quadrature error; changing it across a real peak changes it substantially.
      Asserting both directions shows the parameter is doing what it claims.
- [ ] `tr(u²)/3` from the tensor equals the scalar path computed from `sd.pdos()`, to
      machine precision. If these ever differ, the single-IR claim in `design.md` §1 is
      broken.
- [ ] Units: `u²` comes out in Å², checked dimensionally by a test that varies mass and
      temperature and asserts the expected scaling.

## Depends on

- `03-intermediate-representation.md`
- `09-analytic-test-suite.md`

## Open questions / risks

- **A6 has a shape but not a number.** `design.md` settles that the cutoff is explicit; it
  does not settle what a user should pick. The docs need guidance — tie it to the
  estimator's frequency resolution and to the lowest physical mode — and the guidance
  should be revisited once a real system with diffusive motion has been run.
- **The quantum correction factor is implicit in both papers and they choose differently**
  (method-review §3.6, and Ramírez et al. 2004). The conversion above is one QCF choice.
  Name it in the docstring rather than leaving the reader to infer it, and record it in
  the ambiguity log.

**Labels:** `milestone:M3`, `stage:D`, `physics`, `open-question`
