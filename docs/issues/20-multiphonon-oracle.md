# Write the closed-form multiphonon oracle test

## Context

`plan.md` §5 Layer 2 calls this "the single most valuable test in the suite" and says it
should be written **before** the convolution code. That ordering is not a style preference,
it is the consequence of `design.md` D7.

D7 chose iterative convolution as the only route to multiphonon orders, matching both
papers and both reference codes, so that any discrepancy against OCLIMAX or MolDyINS is
unambiguously ours rather than a method difference. The cost, stated plainly in D7's
consequence, is that we lose the internal cross-check the exponential resummation would
have provided. Verification now rests entirely on analytic tests and external oracles.

For an isotropic 3D harmonic oscillator the exact incoherent result is known:

```
S(Q,ω) = e^{-2W} Σ_n [(Q²u²)ⁿ / n!] δ(ω − nω₀)
```

That is a direct analytic oracle for the convolution chain, the Debye–Waller factor, *and*
the order prefactor — which is exactly the question D7 left open as ambiguity A7. Writing
it first means the prefactor is decided by a closed-form result rather than by a convention
argument between two papers.

## Scope

- A test module constructing the exact single-oscillator `S(Q,ω)` expansion on a discrete
  energy grid, at several `Q` and several temperatures.
- A synthetic `VelocitySpectralDensity` for that oscillator, so the oracle can be fed
  through stages D and E as they are built.
- Assertions on: the `n = 0` elastic weight (`e^{-2W}`), the ratio of successive orders
  (`Q²u²/n`), the position of each order at `nω₀`, and the total summed intensity.
- A parametrisation over candidate order prefactors, so the test *selects* one rather than
  confirming a choice already made in the code.
- Marked fast and unmarked, so it runs on every PR.

## Out of scope

- The convolution implementation itself — `21-multiphonon-convolution.md`. This issue
  delivers the oracle and, initially, a set of failing or skipped tests.
- Powder averaging — `22-isotropic-powder-average.md`.
- Anything with a real trajectory.

## Acceptance criteria

- [ ] The oracle is written from the closed form and checked against itself in a limit
      where it is trivially known: as `Q → 0` all weight is elastic; as `Q²u²` grows the
      distribution over `n` approaches a Poisson with mean `Q²u²`. Both are asserted, so
      the oracle is not merely the code's own output copied into a fixture.
- [ ] The order-ratio assertion distinguishes the candidate prefactors. If two candidate
      conventions give indistinguishable ratios at every accessible `Q`, that finding is
      recorded in the ambiguity log — it would mean A7 cannot be settled this way, which
      is important and must not be quietly ignored.
- [ ] Discretisation error is bounded and stated: the delta functions are represented on a
      grid, so the test's own tolerance is derived from the bin width rather than tuned.
- [ ] Temperature dependence is exercised, not just `T = 0`. The DWF and the occupation
      enter differently, and a `T = 0`-only test cannot separate them.
- [ ] The test fails, clearly and with a message naming the order that disagreed, when the
      prefactor is changed away from the selected one.

## Depends on

- `19-debye-waller-factor.md`

## Open questions / risks

- **A7 is genuinely undetermined in the literature.** Paper B asserts `1/n!` on the grounds
  that it "closely resembles the known prefactors for the first four overtones"; aCLIMAX's
  exact prefactors differ. This test is the plan's designated way of settling it
  (`plan.md` §6). If it does settle it, the evidence goes in the ambiguity log; if it
  settles it only for the isotropic oscillator, say so, because the anisotropic case in
  `27-anisotropic-path.md` then inherits an open question.

**Labels:** `milestone:M3`, `stage:E`, `testing`, `physics`, `open-question`
