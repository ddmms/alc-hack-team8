# Assemble `anisotropic_spectrum`, the Method 2 end-to-end path

## Context

The counterpart to `24-method-1-assembly.md`, and the payoff of `design.md` §1's single-IR
decision: the two methods share stages A–C and F–G, differ only in D and E, and can
therefore be compared on byte-identical input. That comparison is the strongest statement
this project can make about either method, and it is only possible because the two
assembly functions consume the same object.

## Scope

- `anisotropic_spectrum(sd, instrument, temperature, max_order, cutoff, ...) -> Spectrum`,
  matching `plan.md` §3.
- The pipeline: IR → full tensor `u²` → tensor DWF → tensor multiphonon convolution with
  the rank-4 contraction → almost-isotropic average for the fundamental, isotropic above →
  cross-section → `Q(ω)` → broadening → `Spectrum`.
- Shared code with Method 1 wherever the stages are genuinely the same: F and G must not
  be reimplemented here.
- A comparison helper that runs both methods on one IR and returns both spectra plus their
  difference, since that is the operation the whole architecture was arranged to make
  cheap.
- The anisotropic diagnostic from `27-anisotropic-path.md` surfaced in the output metadata.

## Out of scope

- The contraction and the powder average themselves — `27-anisotropic-path.md`,
  `26-almost-isotropic-powder-average.md`.
- External benchmarking — `32-moldyins-benchmark.md`.

## Acceptance criteria

- [ ] On an isotropic system the two methods produce the same spectrum to within
      numerical tolerance, end to end — not just at the powder-average step. This closes
      `plan.md` §7's M4 criterion and is a much stronger test than the unit-level version,
      because it also proves D, F and G are shared rather than duplicated.
- [ ] A test asserts that the two assembly functions call the *same* stage F and G code
      objects. Duplicated broadening would drift, and the divergence would be attributed
      to the physics.
- [ ] On an anisotropic system the difference is nonzero, localised to the fundamental
      (since higher orders fall back to isotropic), and in the direction the uniaxial
      analysis in `26-almost-isotropic-powder-average.md` predicts.
- [ ] Both methods run on the same IR object without either mutating it — the IR is
      frozen, so this should be structurally guaranteed, but assert it.
- [ ] Output metadata distinguishes the method used, so two spectra in a directory are
      never ambiguous.
- [ ] The comparison helper is used to generate a figure for the documentation, so the
      method difference is shown rather than described.

## Depends on

- `24-method-1-assembly.md`
- `27-anisotropic-path.md`
- `28-finite-temperature-and-dwf-variants.md`

**Labels:** `milestone:M4`, `stage:D`, `stage:E`, `physics`
