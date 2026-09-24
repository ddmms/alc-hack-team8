# Implement isotropic powder averaging

## Context

`design.md` §2[E] makes powder averaging a strategy chosen per method: fully isotropic for
Method 1, almost-isotropic for the fundamental in Method 2. Both are small, well-defined
functions over `tr(A)`, `tr(B)` and `B:A`, which is why they belong in one module with a
common interface rather than being inlined into each method.

Method 1 is isotropic by construction. Paper A's Eq. 7 delivers only `tr(B_ds)`, and the
paper says the reduced information "is sufficient" — implying a silent fall back to the
fully isotropic Debye–Waller factor for MD input, which method-review §3.2 records as
never confirmed anywhere (ambiguity A5). Implementing the isotropic average explicitly,
with its assumptions named, is what makes that assumption testable later against the
OCLIMAX oracle instead of buried.

## Scope

- `src/mdins/scattering/powder.py` with a common signature over `tr(A)`, `tr(B)` and
  `B:A`, so the isotropic and almost-isotropic strategies are interchangeable at the call
  site.
- The fully isotropic average for the fundamental and for all higher orders.
- The isotropic Debye–Waller form `W = Q² tr(A)/6` (Paper A Eq. 4).
- A docstring stating the assumption Paper A leaves implicit, with the A5 reference, so
  the reader knows what has been assumed on their behalf.

## Out of scope

- The almost-isotropic average — `26-almost-isotropic-powder-average.md`. Keeping them
  separate keeps Method 1 shippable without Method 2.
- Coherent scattering. Absent from both papers and out of scope for the whole project
  (method-review §3.8).

## Acceptance criteria

- [ ] For an isotropic input tensor, the average reproduces the closed-form
      `⟨(Q·u)²⟩ = Q² tr(u²)/3` exactly, not approximately. The isotropic case has an exact
      answer and should be asserted as one.
- [ ] Rotational invariance: applying a random rotation to the input tensor leaves the
      result unchanged to machine precision. This is the property a powder average must
      have, and it catches an index error in the contraction that a scalar-only test
      cannot.
- [ ] The DWF form matches Paper A Eq. 4 at a hand-computed `Q` and tensor.
- [ ] The function is pure and array-shaped over atoms and energies, so stage E's
      per-atom parallelism (`design.md` §3) is not destroyed by a Python loop here.

## Depends on

- `19-debye-waller-factor.md`

**Labels:** `milestone:M3`, `stage:E`, `physics`
