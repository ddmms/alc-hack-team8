# Add the Einstein-crystal end-to-end test

## Context

`plan.md` §6 identifies the central difficulty of M2: when the LJ comparison against
Euphonic disagrees, is it our bug, anharmonicity, or a q-sampling mismatch? Three
candidate explanations and no way to separate them is not a validation, it is an argument.

The Einstein crystal removes two of the three. ASE atoms on independent harmonic springs,
NVE, fixed seed: the pDOS is a delta at a known frequency, there is no dispersion so
q-sampling cannot matter, and the potential is exactly harmonic so anharmonicity cannot
matter. Anything that disagrees is ours. `plan.md` §5 Layer 3 puts it first for exactly
this reason, and `plan.md` §7 makes it half of M2's definition of done.

It is also the first test that exercises real trajectory I/O rather than synthetic arrays,
so it is where stage A's unit conversions get proved against a physical answer.

## Scope

- A test fixture building an ASE `Atoms` object with a harmonic spring calculator tethering
  each atom to its lattice site, with a spring constant chosen to give a convenient
  frequency in the tens of meV.
- NVE dynamics with a fixed RNG seed for the initial velocity distribution, written to a
  real trajectory file and read back through `read_velocities`, so stage A is in the loop.
- The full A–C chain to a `VelocitySpectralDensity`.
- Assertions on peak position, peak width relative to the estimator's own resolution, the
  equipartition sum rule against the actual MD temperature, and isotropy of the tensor.
- Marked `slow`.

## Out of scope

- Euphonic. This test's whole value is that it needs no external reference.
- LJ and dispersion — `17-euphonic-comparison.md`.
- Anything downstream of the IR.

## Acceptance criteria

- [ ] The recovered peak sits at the analytic Einstein frequency `√(k/m)` expressed in meV,
      within the estimator's frequency resolution `1/(L·dt)`. The test states that
      resolution explicitly rather than using a bare `rtol`, so a failure distinguishes "we
      are wrong" from "the grid is coarse".
- [ ] The sum rule holds against the *measured* MD temperature, not the requested one. In
      NVE with a thermalised start the two differ by roughly the expected equipartition
      split, and using the requested temperature would build a small systematic error into
      the benchmark everything else is compared against.
- [ ] The tensor is isotropic to within the inter-segment spread: the three diagonal
      components agree and the off-diagonals are consistent with zero. An anisotropic
      result from an isotropic potential means a bug in the tensor construction or in the
      Å/ps conversion per axis.
- [ ] A negative control: a spring constant 5% different fails the same assertion. Without
      it, the tolerance may simply be wide enough to admit anything.
- [ ] Both estimators pass, with the same tolerances.
- [ ] The test is deterministic across runs and across platforms to the stated tolerance;
      the seed and the integrator settings are in the fixture, not in a comment.

## Depends on

- `09-analytic-test-suite.md`

**Labels:** `milestone:M2`, `stage:A`, `stage:B`, `stage:C`, `testing`, `physics`
