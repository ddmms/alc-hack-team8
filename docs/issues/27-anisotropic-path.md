# Implement the anisotropic tensor path and its rank-4 contraction

## Context

Method 2's distinguishing feature, and the reason the IR stores a tensor at all. Where
Method 1 takes `tr(u²)/3` at stage D and proceeds with scalars, Method 2 keeps the full
3×3 tensor through the displacement model and into the multiphonon chain, so vibrational
anisotropy survives into the powder average.

The difficulty is Paper B Eq. 17. Self-convolving `B` to get overtones produces a rank-4
object, which the paper reduces to rank 2 "approximated as a simple matrix multiplication",
with no derivation and no error estimate (method-review §4.3, ambiguity A8).
`design.md` §5 resolves this as "implement as published" — and then states the cost
plainly: with the exponential resummation ruled out by D7, the error **can only be bounded
against an external oracle, not internally.** That is a known, accepted limitation of this
issue, not something its tests can fix.

## Scope

- `src/mdins/scattering/anisotropic.py`: the tensor path through stages D and E.
- Tensor self-convolution with the rank-4 → rank-2 contraction implemented as published,
  isolated in one clearly named function with the approximation documented at the point it
  is made.
- The contraction kept swappable, so an alternative can be tried later without restructuring
  the caller — the approximation is uncontrolled and may need revisiting.
- Method 2's Debye–Waller factor using the full tensor rather than its trace.
- A diagnostic reporting the anisotropy of each atom's tensor (e.g. the ratio of largest to
  smallest eigenvalue), so a user can see whether the anisotropic path is buying them
  anything on their system.

## Out of scope

- The powder average that consumes the result — `26-almost-isotropic-powder-average.md`.
- Finite-temperature weighting and the DWF heuristic —
  `28-finite-temperature-and-dwf-variants.md`.
- Bounding the contraction error against MolDyINS or OCLIMAX — `32-moldyins-benchmark.md`.

## Acceptance criteria

- [ ] **On isotropic input, the anisotropic path reproduces the isotropic path exactly**,
      order by order. This is half of M4's definition of done (`plan.md` §7) and the only
      strong internal check available. A discrepancy at a single order localises the bug
      to the contraction.
- [ ] On anisotropic input the two differ in the expected direction — the other half of
      `plan.md` §7's M4 criterion. "Expected direction" must be made concrete: for a
      uniaxial tensor, state which way the fundamental's intensity should move and why,
      in the test docstring.
- [ ] Rotational covariance: rotating the input tensor rotates the output tensor, to
      machine precision, at every order. This is what catches a transposed index inside
      the contraction, which isotropic input cannot expose.
- [ ] The contracted tensor stays symmetric and positive semi-definite at every order for
      PSD input. If the published contraction does not preserve this, that is a finding
      about the approximation and must be reported and recorded, not clipped away.
- [ ] Memory and time scale acceptably: the rank-4 intermediate is never materialised for
      all atoms at once.
- [ ] The docstring at the contraction names Paper B Eq. 17, calls the approximation
      uncontrolled, and points at the ambiguity log entry.

## Depends on

- `21-multiphonon-convolution.md`
- `26-almost-isotropic-powder-average.md`

## Open questions / risks

- **A8 is unbounded by design.** No test written here can tell us how wrong the
  approximation is; the best available evidence is a comparison against MolDyINS (which
  makes the same approximation, so agreement proves only that we reproduced it) and
  against OCLIMAX (which does not use it, so a discrepancy confounds the approximation
  with the method difference). State this limitation in the user documentation; a user
  choosing the anisotropic method deserves to know its error is uncharacterised.

**Labels:** `milestone:M4`, `stage:D`, `stage:E`, `physics`, `open-question`
