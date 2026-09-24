# Implement the almost-isotropic powder average

## Context

Method 2's powder average for the fundamental, following aCLIMAX (`design.md` §2[E]).
Above the fundamental it falls back to isotropic, which is what Paper A′ does and what
makes the anisotropic method tractable at all.

The source situation is awkward and worth knowing before starting. The analytic
almost-isotropic average is Paper A′ Eq. 4, derived in aCLIMAX (Ramirez-Cuesta 2004) and
in Mitchell et al.'s *Vibrational Spectroscopy with Neutrons* — which method-review §6
notes is Paper B's ref 11, cited eight times and never reproduced. So the expression has
to be taken from A′ and checked, not copied from the paper that uses it most.

The average needs `B:A` and `tr(B)` separately (method-review §3.2), which is exactly the
information Paper A's trace-only pipeline throws away. Retaining it is the whole point of
the tensor IR.

## Scope

- An `almost_isotropic` strategy in `src/mdins/scattering/powder.py`, sharing the interface
  established by `22-isotropic-powder-average.md`.
- Implementation of A′ Eq. 4 in terms of `tr(A)`, `tr(B)` and `B:A`, with the expression
  written out in the docstring and its source cited by equation number.
- Dispatch: almost-isotropic for the fundamental, isotropic for orders above it, with the
  switch explicit and configurable.

## Out of scope

- The anisotropic displacement tensor path that feeds it — `27-anisotropic-path.md`.
- The isotropic average — `22-isotropic-powder-average.md`.

## Acceptance criteria

- [ ] **On isotropic input, the almost-isotropic path returns exactly what the isotropic
      path returns.** `plan.md` §5 Layer 2 names this test specifically, and it is the
      strongest cheap check available: the two expressions are structurally different and
      agreeing to machine precision on the isotropic special case means the contraction
      indices are right.
- [ ] Rotational invariance of the result under a random rotation applied to both `A` and
      `B` together, to machine precision.
- [ ] On a deliberately anisotropic tensor, the result differs from the isotropic one in
      the direction and by roughly the magnitude the closed form predicts for a uniaxial
      case worked out by hand in the test docstring. "Differs" alone is not an acceptance
      criterion.
- [ ] The result is bounded: it must not go negative for any physically admissible
      (positive semi-definite) `A` and `B`. Randomised PSD inputs assert this.
- [ ] Order dispatch is tested: the second order goes through the isotropic path even when
      the strategy is almost-isotropic.

## Depends on

- `22-isotropic-powder-average.md`

**Labels:** `milestone:M4`, `stage:E`, `physics`
