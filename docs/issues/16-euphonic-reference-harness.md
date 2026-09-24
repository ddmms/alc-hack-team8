# Build the Euphonic harmonic reference harness

## Context

Tier 2 is the core validation (`design.md` §4), and its reference side is a phonon
calculation on the same potential and the same cell that the MD uses:

```
ASE LennardJones calculator
   ├── force constants (finite displacement) → Euphonic ForceConstants
   │        → calculate_pdos → reference atom-projected DOS
   └── MD (same calculator, same cell, low T) → our pipeline → pDOS
```

This issue is the top branch only. `plan.md` §4 lists the validation harness as one of the
three tracks that can be built in parallel immediately after M0, because it can be written
against the frozen IR before there is a real implementation behind it — the reference does
not depend on our pipeline at all.

`design.md` §4 flags the one thing that will otherwise produce a spurious mismatch:
**Brillouin-zone sampling must match.** An MD supercell samples only the q-points
commensurate with it. A reference computed on a dense mesh will be smoother than anything
MD can produce, and the comparison will fail for a reason that has nothing to do with our
code.

## Scope

- A helper that builds force constants by finite displacement from an ASE calculator, for
  a given unit cell and supercell.
- Conversion to a Euphonic `ForceConstants` object.
- Computation of the commensurate q-grid for a given supercell, as an explicit, separately
  testable function — not a line buried inside the comparison.
- `calculate_pdos` on that grid, returning a `Spectrum1DCollection` with per-species
  projections.
- Conversion of both sides into a common comparison space: normalised lineshapes on a
  shared energy grid, since Euphonic's convention differs from the IR's absolute one (D6).
- Everything behind the `euphonic` extra and marked `reference`, so `import mdins` stays
  free of it.

## Out of scope

- The comparison and its tolerance — `17-euphonic-comparison.md`.
- Running the MD side.
- Neutron cross-section weighting; the comparison is on pDOS.

## Acceptance criteria

- [ ] The commensurate q-grid function returns exactly the `N1·N2·N3` points commensurate
      with an `N1×N2×N3` supercell, verified against a hand-enumerated 2×2×2 case. Getting
      this wrong is the failure mode `design.md` §4 warns about, and it fails silently —
      the reference just looks smoother.
- [ ] Finite-displacement force constants reproduce the analytic dynamical matrix for a
      1D chain of harmonic springs, so the displacement magnitude and the symmetrisation
      are checked before LJ introduces anharmonicity.
- [ ] The acoustic sum rule holds on the computed force constants, and the three acoustic
      modes at Γ come out at zero within tolerance. A nonzero Γ acoustic frequency
      indicates a broken sum rule and would shift the whole low-energy comparison.
- [ ] The reference pDOS integrates to the same value per species under the normalisation
      helper, so the comparison space is well defined before anything is compared in it.
- [ ] Importing the harness without `euphonic` installed produces a skip, not an error.
- [ ] The harness is usable from a script as well as from pytest, since the write-up in
      `17-euphonic-comparison.md` needs to generate figures from it.

## Depends on

- `03-intermediate-representation.md`
- `05-continuous-integration.md`

## Open questions / risks

- The finite-displacement step size is a free parameter that trades truncation error
  against numerical noise. Pick it by convergence on the Γ frequencies rather than by
  copying a default, and record the choice — an under-converged reference would be
  attributed to our MD side.

**Labels:** `milestone:M2`, `testing`, `parallel-track`
