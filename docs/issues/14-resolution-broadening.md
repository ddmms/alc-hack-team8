# Implement energy-dependent resolution broadening

## Context

Stage F's second half (`design.md` §2[F]): a Gaussian whose width varies with energy
transfer, in practice a polynomial in it. The implementation itself is small; the reason
it is its own issue is the dependency decision behind it.

`design.md` D4 makes Euphonic a core dependency and takes `broadening` and the Sears
cross-section data from it rather than reimplementing either. That became permissible only
after the relicence in D1. Variable-width convolution is easy to get subtly wrong —
normalisation drift across the energy axis, or a kernel truncated too aggressively at the
wings — and Euphonic's implementation is already used by the neutron community.

## Scope

- `src/mdins/instrument/resolution.py`: evaluate the polynomial width model from an
  instrument definition, and apply variable-width Gaussian broadening to a spectrum on a
  given energy grid, via Euphonic's `broadening`.
- Per-bank broadening where the instrument defines several banks, combined with the bank
  weights.
- A fixed-width mode for testing and for instruments defined that way.
- Explicit handling of the grid: the broadening kernel must be resolved by the grid it is
  applied on, and a warning or error when the finest width is comparable to the bin
  spacing — an under-resolved Gaussian silently becomes a delta.
- Documentation of which Euphonic broadening routine is used and why, since the choice
  affects the wings.

## Out of scope

- The instrument constants — `13-tosca-resolution-constants.md`.
- Cross-checking against AbINS — `33-abins-broadening-crosscheck.md`. That is Layer 5 and
  needs a spectrum to compare.
- Kinematics — `11-instrument-kinematics.md`.

## Acceptance criteria

- [ ] Broadening a unit delta with a constant width returns a Gaussian of that width and
      unit integral, to the grid's quadrature accuracy.
- [ ] Broadening conserves total integrated intensity under a *variable* width, not only a
      constant one. This is the test that catches the most likely bug: a kernel normalised
      once at the mean width rather than per bin, which transfers intensity from the
      sharp end of the spectrum to the broad end and looks entirely plausible.
- [ ] Two deltas separated by less than the local width merge into one peak; separated by
      more, they resolve. Assert peak count, so the test says something about the physics
      rather than about array values.
- [ ] Broadening at TOSCA's width reproduces a hand-computed FWHM at a stated energy.
- [ ] A grid too coarse to resolve the narrowest width raises or warns, naming both the
      width and the bin spacing.
- [ ] Multi-bank broadening with equal weights and identical banks equals single-bank
      broadening exactly.

## Depends on

- `12-instrument-data-format.md`

**Labels:** `milestone:M3`, `stage:F`, `physics`, `parallel-track`
