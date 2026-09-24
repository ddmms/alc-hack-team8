# Implement finite-temperature Boltzmann weighting and the Debye–Waller variants

## Context

Paper B builds a 0 K base spectrum and reaches finite temperature by convolving against a
Boltzmann-weighted DOS (Eqs. 18–20). Two of the equations involved are, as printed,
wrong or unjustified, and this issue exists to implement them with that stated rather than
silently corrected.

**Eq. 19 diverges as printed.** `ρ(ω;T) = Σ_a exp(+ħω_a/k_BT) δ(ω+ω_a)` has a positive
exponent and no partition function, so it grows without bound with mode energy
(method-review §4.4). `design.md` §5 A3 treats this as a typo and uses `exp(−ħω/k_BT)/Z`,
flagged in the docs. That is a judgement about someone else's paper and must be visible in
the code, not buried.

**Eq. 21's DWF is a heuristic.** `∫ max[(q·u_i(ω))², (q·ω v_i(ω))²] dω` takes whichever of
the ground-state or thermal displacement is larger at each frequency, with no derivation
and no continuity guarantee at the crossover (method-review §4.5, ambiguity A4).
`design.md` §5 says implement it as published, behind a flag, and prefer the consistent
BE-weighted form as the default.

## Scope

- Boltzmann-weighted DOS construction with the corrected exponent and an explicit
  partition function, plus the convolution taking the 0 K spectrum to finite temperature.
- The `max[...]` DWF from Eq. 21 as a selectable, non-default option alongside the
  BE-weighted default from `19-debye-waller-factor.md`.
- A warning when the published variant is selected, naming it as a heuristic without
  derivation.
- Module docstrings recording both the printed equation and what was implemented instead,
  with the method-review section numbers.
- Entries in the ambiguity log for A3 and A4 with the evidence.

## Out of scope

- The anisotropic tensor machinery — `27-anisotropic-path.md`.
- Paper B's anharmonic correction (Eqs. 22–23). method-review §4.6 records that it was
  derived, then dropped by its own authors, who note the classical decay rates correspond
  to the thermalised initial state rather than the vibrational ground state the theory
  needs. Out of scope for this project entirely.

## Acceptance criteria

- [ ] The Boltzmann-weighted DOS is normalised: it integrates to one, and `Z` is computed
      rather than assumed. A test evaluates it at a high mode energy and asserts it decays;
      the literal Eq. 19 grows, so this test is precisely the one that distinguishes the
      correction from the published form.
- [ ] At low temperature the finite-temperature convolution returns the 0 K spectrum to
      within tolerance, so the correction reduces correctly.
- [ ] Detailed balance is checked on the result. method-review §6 notes Paper B Eq. 20
      reduces to the Maradudin–Fein detailed-balance result, which gives an independent
      target rather than only self-consistency.
- [ ] The `max[...]` DWF is exercised by a test that places the crossover inside the
      integration range and asserts the discontinuity in the integrand is real and
      measurable. The point is to characterise the heuristic's behaviour, not to hide it.
- [ ] The two DWF variants differ on a representative system, and the difference is
      quantified in the docs. If they agree everywhere we tried, say that — it would mean
      A4 does not matter in practice, which is a useful result.
- [ ] Selecting the published variant emits a warning and records the choice in
      provenance and output metadata.

## Depends on

- `27-anisotropic-path.md`

## Open questions / risks

- **A3 is a judgement call on a published equation.** We are confident Eq. 19 is a typo
  because it diverges, but we have not confirmed it with the authors. The code implements
  the corrected form and the docs must say so plainly; if anyone does get confirmation,
  it goes in the ambiguity log.
- **A4 has no derivation to check against.** The acceptance criteria above characterise
  the heuristic rather than validating it, because there is nothing to validate it
  against. That is the honest position.

**Labels:** `milestone:M4`, `stage:D`, `stage:E`, `physics`, `open-question`
