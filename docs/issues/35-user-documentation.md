# Document the velocity-input constraint and the user workflow

## Context

`design.md` §2[A] asks for one thing explicitly of the documentation: that the
velocities-only constraint be "stated plainly in user-facing docs rather than discovered".
That is the sharpest usability consequence of D2, and it lands hardest on exactly the
audience Paper B was written for — GROMACS users, who cannot bring a `.trr` at all.

The rest of the documentation exists because this package has several parameters with no
safe default (the low-frequency cutoff), several with defaults that are choices rather than
facts (the order prefactor, the DWF variant), and one approximation that must never be used
accidentally (temperature decoupling, D8). A user who does not know these exist will get a
plausible-looking number.

## Scope

- A "what you need before you start" page: velocities in the trajectory, the supported
  formats table from `design.md` §2[A], the conversion route for everything else, and the
  sampling requirements `dt < 1/(2 c ν_max)` and `t_total > 1/(c Δν)` expressed in
  practical terms.
- A thermostat note: NVE or a weakly decorrelating thermostat for production runs; τ ≈ 1 ps
  corrupts the spectrum; barostats damp the dynamics (method-review §4.7).
- A worked end-to-end example, CLI and Python, from trajectory to spectrum.
- A page on the parameters that are choices, not defaults: the low-frequency cutoff, the
  multiphonon prefactor, the DWF variant, `max_order`, and the two temperatures — each
  with what it does, what happens if it is wrong, and where to read more.
- A short statement of scope limits: incoherent scattering only, harmonic quantum
  treatment, no coherent branch.
- A licence section restating the copyright holder (the outstanding housekeeping under
  D1).
- API reference generated from docstrings.

## Out of scope

- The validation write-ups — those belong to `17-euphonic-comparison.md` and
  `34-published-spectra.md`.
- The ambiguity log — `36-ambiguity-resolution-log.md`.

## Acceptance criteria

- [ ] The formats table appears before the installation instructions, not after the
      tutorial. A user with a DCD file should find out in the first screenful.
- [ ] The worked example is executed in CI (as a doctest or a notebook run), so it cannot
      rot into referring to an API that no longer exists.
- [ ] Every parameter with no safe default is documented with the failure mode of getting
      it wrong, in concrete terms — "the Debye–Waller factor diverges and the elastic line
      swallows the spectrum", not "choose carefully".
- [ ] Temperature decoupling is documented as an approximation with no error estimate, per
      D8, and the documentation says when it is appropriate rather than only how to invoke
      it.
- [ ] A reader who has phonon data from another source can find out how to enter the
      pipeline at stage C, since `design.md` §1 advertises that as supported.
- [ ] Docs build in CI and a broken internal link fails the build.

## Depends on

- `25-spectrum-cli-subcommand.md`

**Labels:** `milestone:M3`, `docs`
