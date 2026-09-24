# Implement the Welch cross-spectral density estimator

## Context

Stage B's default route to `P_{i,αβ}(E) = ⟨v*_{i,α} v_{i,β}⟩`, per `design.md` D5. The
choice of Welch over the VACF route the source papers use is made on positive
semi-definiteness, not accuracy: Welch averages outer products `V*V^T`, so every `P_i(E)`
is PSD by construction. Blackman–Tukey gives no such guarantee, and a non-PSD tensor
means a negative mean-square displacement along some direction — cosmetic for Method 1,
physically meaningless for Method 2.

The second reason is `plan.md` §2's other fixed convention: **the output energy grid is an
input to this stage, not a post-hoc downsample.** Passing the target grid into the
estimator lets it rebin each segment as it is formed, so the full-resolution per-atom
array never materialises. Without it, 3000 atoms at 25000 native bins is 5.4 GB before
anything useful has happened, which is the wall both source papers hit.

## Scope

- `src/mdins/spectral.py` and the `velocity_spectral_density(traj, frequencies=grid,
  estimator="welch", ...)` entry point, returning a `VelocitySpectralDensity`.
- Segmentation with configurable segment length, overlap and window (Hann by default),
  with every parameter recorded in `SpectralMetadata`. These are undocumented choices in
  both source codes; here they are documented, testable arguments.
- Per-segment rebinning onto the target grid, conserving the integral.
- Chunking over atoms so peak memory is bounded by the chunk, not the system — forming one
  segment's cross-spectrum is itself an `(n_native, n_atoms, 3, 3)` intermediate even when
  the kept array is small.
- Normalisation such that `∫ tr P dE = 3 k_B T_MD / m_i` holds directly (`design.md` D6),
  including the window's power correction and the one-sided-spectrum factor of two.
- The real part (co-spectrum) taken for the stored tensor, with the discarded imaginary
  (quadrature) part returned or recorded as a diagnostic magnitude rather than thrown
  away silently — ambiguity A2.
- Inter-segment spread returned alongside the mean, since `design.md` D5 makes it the
  basis for every tolerance in `plan.md` §5 and a stationarity check for free.
- Recording `energy_resolution = 1/(L·dt)` separately from the bin count, so a user who
  asks for a finer grid than the segment length supports is told they are interpolating.

## Out of scope

- The VACF estimator — `08-vacf-estimator.md`.
- The analytic tests that prove the normalisation — `09-analytic-test-suite.md`. They
  should be written first, but they are their own deliverable.
- Trajectory reading and pre-processing — `06-trajectory-reader.md`.

## Acceptance criteria

- [ ] For synthetic velocities the estimator reproduces `scipy.signal.csd(...,
      scaling='density')` on a single pair of components, so the normalisation is anchored
      to an external implementation and not only to our own sum rule.
- [ ] Every returned `P_i(E)` has non-negative eigenvalues to within numerical tolerance,
      asserted over randomised multi-atom input. This is the assertion the whole choice of
      Welch exists to make true, so its absence would make D5 unjustified.
- [ ] `∫ tr P dE` equals `3 k_B T_MD / m_i` per atom for synthetic thermal velocities, with
      the integral computed by the same bin-width helper the IR uses.
- [ ] Changing the target grid from 25000 to 500 bins changes peak memory by roughly the
      expected factor and leaves the integral unchanged. Measure allocation, not wall
      time; the point of the design is memory.
- [ ] Changing the window from Hann to a rectangular one changes the linewidth but not the
      integral, to the stated tolerance. A missing window power correction fails this and
      passes everything else.
- [ ] Doubling the segment count halves the variance of the estimate across repeats, to
      within sampling error — the `∝ 1/K` behaviour D5 claims.
- [ ] Atom chunking is exercised by a test that sets the chunk size small enough to force
      several chunks and asserts the result is bit-identical to the unchunked path.
- [ ] Requesting an energy grid finer than `1/(L·dt)` either raises or records the
      interpolation explicitly; it does not quietly return a smooth-looking spectrum with
      more bins than information.

## Depends on

- `03-intermediate-representation.md`
- `06-trajectory-reader.md`

## Open questions / risks

- **A2, the fate of the quadrature part.** `P_αβ` is Hermitian in general but the
  displacement tensor it stands in for is real symmetric, and Paper B prints the full
  Hermitian matrix without saying what happens to the imaginary part. The working
  resolution (`design.md` §5) is to take the co-spectrum and assert the imaginary part is
  negligible in the harmonic limit. What "negligible" means quantitatively is not settled
  and should be fixed by the circular-motion test in `09-analytic-test-suite.md`, where
  the quadrature part is large by construction, rather than by picking a threshold here.

**Labels:** `milestone:M1`, `stage:B`, `physics`
