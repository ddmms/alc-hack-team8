# Build the Layer 1 analytic test suite and assert the IR invariants

## Context

`plan.md` §6 lists silent factor errors — 2π, mass, ħ — as the first risk, and the
mitigation is that this layer is written *first*. Velocities are generated directly rather
than by running MD, so the expected spectrum is exact and the tests run in milliseconds.
By the time a real trajectory is involved, a normalisation error is expensive to
distinguish from physics; here it is unambiguous.

`plan.md` §5 Layer 1 enumerates the cases and the invariants. The invariants matter as
much as the cases: `design.md` D5 and D6 are both justified by assertions that have to
actually exist, or the decisions rest on nothing.

## Scope

- `tests/test_spectral.py` (and helpers) covering, for both estimators:
  - single sinusoid → delta at a known energy, amplitude fixed by the sum rule;
  - sum of sinusoids at different amplitudes → correct relative weights;
  - damped oscillator → Lorentzian of known width;
  - white noise → flat spectrum, Parseval holds;
  - anisotropic signal, different amplitude per Cartesian axis → correct diagonal tensor,
    zero off-diagonals;
  - correlated x/y motion → known off-diagonal magnitude **and sign**;
  - circular motion → nonzero quadrature part, which is where ambiguity A2 is visible.
- Property-style invariants over randomised inputs:
  - `∫ tr P dE == 3 k_B T / m` (equipartition, D6's anchor);
  - every `P_i(E)` has non-negative eigenvalues (D5's justification);
  - Welch and VACF agree within the Welch inter-segment spread (D5's cross-check);
  - HDF5 round-trip is exact, including metadata.
- A shared synthetic-trajectory helper, so a test that constructs drift-free velocities
  does not have to claim a COM removal that never ran — provenance is load-bearing for
  the sum-rule correction.
- A quantitative threshold for "the imaginary part is negligible", derived from the
  circular-motion case where it is large by construction, rather than asserted.

## Out of scope

- Anything involving a real trajectory — `15-einstein-crystal.md` onwards.
- Multiphonon analytic tests — `20-multiphonon-oracle.md`, which is Layer 2 and belongs
  with the scattering work.

## Acceptance criteria

- [ ] Every test states its exact expected value in the docstring, with the derivation in
      one line. A test whose expected value is "whatever the code produced" is not part of
      this layer.
- [ ] The single-sinusoid amplitude test fails if the estimator's window power correction
      is removed, and fails if the one-sided factor of two is removed. Verify by
      deliberately breaking each and confirming a red test, then reverting — a
      normalisation test that passes under both conventions is not testing normalisation.
- [ ] The off-diagonal sign test distinguishes `+C_xy` from `-C_xy`; a test on `|C_xy|`
      does not count.
- [ ] The PSD invariant runs on randomised multi-atom, multi-frequency input and reports
      the most negative eigenvalue and the atom it belongs to on failure.
- [ ] The Welch/VACF agreement test states its tolerance as a multiple of the inter-segment
      spread, not as a literal. A hand-picked `rtol` here would defeat the entire point of
      D5's uncertainty estimate.
- [ ] The whole layer runs in under a few seconds and is unmarked, so it runs on every PR.

## Depends on

- `04-hdf5-serialisation.md`
- `07-welch-estimator.md`
- `08-vacf-estimator.md`

**Labels:** `milestone:M1`, `stage:B`, `testing`
