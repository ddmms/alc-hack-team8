# Implement iterative multiphonon convolution and the order prefactor

## Context

Stage E's core, per `design.md` §2[E] and D7. Orders are built as `S_n = S_1 * S_{n−1}`,
exactly as both papers do, and the convolution happens **at atomic level, before**
cross-section and Debye–Waller factor are applied — Paper A is explicit about that ordering
and getting it wrong changes the result, not just the arithmetic path.

Two consequences of D7 need handling here rather than noting. The order prefactor becomes
a first-class configurable choice rather than an implementation detail, because there is no
exponential resummation to derive it from and the two published conventions differ
(ambiguity A7). And per-order results must be retained, not summed away: `design.md` §2[G]
notes both papers show order-resolved plots and that the decomposition is diagnostically
essential.

## Scope

- `src/mdins/scattering/convolution.py`: iterative self-convolution to a configurable
  `max_order`, defaulting to 10 (AbINS's ceiling).
- The order prefactor as an explicit, documented, configurable parameter, with the
  conventions named — `1/n!` (Paper B) and the exact aCLIMAX prefactors — rather than a
  bare float.
- Per-order output retained through to stage G.
- An intensity threshold below which higher orders are dropped, with the cutoff recorded;
  AbINS's `s_relative_threshold = 0.01` and `min_order = 3` are a reasonable starting
  point and should be cited as such.
- Grid handling: convolution on a discrete grid needs a fine enough bin to represent the
  fundamental before it is convolved. AbINS's `fine_bin_factor = 10` exists for this
  reason; the internal grid refinement should be explicit and recorded.
- Detailed balance across the energy-gain and energy-loss sides, and its preservation
  through the convolution chain.

## Out of scope

- The analytic oracle — `20-multiphonon-oracle.md`, which is written first.
- Powder averaging — `22-isotropic-powder-average.md` and `26-almost-isotropic-powder-average.md`.
- The anisotropic rank-4 contraction — `27-anisotropic-path.md`.

## Acceptance criteria

- [ ] The oracle test from `20-multiphonon-oracle.md` passes, order by order, at multiple
      `Q` and temperatures.
- [ ] Convolution is associative and order-independent to numerical tolerance: computing
      `S_4` as `S_1*S_3` and as `S_2*S_2` agrees. A discrepancy here means the prefactor
      is being applied inside the recursion rather than once per order.
- [ ] Total intensity summed over orders is conserved against the analytic total, so the
      threshold and the `max_order` truncation are shown to be discarding what they claim.
- [ ] Changing `max_order` from 10 to 20 changes the result by less than the recorded
      threshold on a representative spectrum, which is the evidence that 10 is adequate
      rather than merely inherited from AbINS.
- [ ] Reducing the internal grid refinement factor visibly degrades the higher orders,
      demonstrating that the refinement is load-bearing and not cargo-culted.
- [ ] Detailed balance holds between the gain and loss sides at each order:
      `S(Q,−ω)/S(Q,+ω) = exp(−ħω/k_BT)` to tolerance. method-review §3.5 records that
      Paper A never describes how this is preserved through the convolution, so it has to
      be asserted rather than assumed.
- [ ] Per-order arrays reach the output container; a test asserts the summed spectrum
      equals the sum of the retained orders.

## Depends on

- `20-multiphonon-oracle.md`

## Open questions / risks

- **A7 (the prefactor) is settled by the oracle, not here.** If `20-multiphonon-oracle.md`
  could not distinguish the conventions, this issue ships both behind the parameter with
  the default documented as provisional, and the comparison against OCLIMAX in
  `31-oclimax-oracle.md` becomes the deciding evidence.
- **Detailed balance through the convolution is not described in either paper.** If
  asserting it turns out to be inconsistent with the published construction, that is a
  finding about the method, not a bug to paper over — record it.

**Labels:** `milestone:M3`, `stage:E`, `physics`, `open-question`
