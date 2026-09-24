# Implement the indirect- and direct-geometry kinematics

## Context

Stage F, first half (`design.md` §2[F]). This is the most self-contained piece of work in
the project: it takes an energy grid and returns `Q(ω)`, with no MD, no IR and no
scattering model involved. `plan.md` §4 names the instrument layer as one of the three
tracks that can run in parallel immediately after M0, and this is its entry point.

It is also a gap in the source material. Paper B hardcodes VISION with `E_f = 32 cm⁻¹`
and never writes the `Q(ω)` relation down at all (method-review §4.8), so there is nothing
to copy — but the physics is standard and roughly a hundred lines from Squires.

## Scope

- `src/mdins/instrument/kinematics.py`.
- Indirect geometry: `Q² = (2m_n/ħ²)(E_i + E_f − 2√(E_i E_f) cos θ)` with `E_i = E_f + ħω`,
  parameterised by final energy and scattering angle.
- Direct geometry: the same relation with `E_i` fixed and `E_f = E_i − ħω`.
- Multi-bank support: a detector bank is an angle (or an angle range), and an instrument
  may have several. The return should make the per-bank trajectory available, not only an
  average.
- The neutron-mass grouping expressed through `units.py`, not inlined.
- Validity handling: in direct geometry `E_f` goes negative above `E_i`, and the
  accessible energy range must be reported rather than producing NaNs downstream.

## Out of scope

- Resolution broadening — `14-resolution-broadening.md`.
- Where the instrument parameters come from — `12-instrument-data-format.md`.
- Applying `Q(ω)` to a scattering function — `24-method-1-assembly.md`.

## Acceptance criteria

- [ ] `Q(ω)` matches hand-computed values at a few energies for both geometries, with the
      hand calculation written out in the test docstring so a reviewer can check it
      without re-deriving Squires.
- [ ] The computed TOSCA trajectory reproduces the published `Q` vs energy-transfer curve
      for the forward and backscattering banks to within reading error of the figure. This
      is the test that catches a wrong `E_f` or a degrees/radians slip, neither of which
      the hand-computed points at low energy would expose.
- [ ] At zero energy transfer, indirect geometry gives `Q² = (2m_n/ħ²)·2E_f(1 − cos θ)`,
      checked analytically.
- [ ] Energy transfers outside the accessible range raise or return a masked result with
      the accessible range named; a test asserts direct geometry at `ω > E_i` does not
      silently return `nan`.
- [ ] `Q` is monotonic in energy transfer for indirect geometry over the instrument's
      range — a cheap invariant that catches a sign error in `E_i = E_f + ħω`.
- [ ] No MD, IR or scattering import appears in this module, so the track stays
      independent.

## Depends on

- `02-units-and-mev-boundary.md`

**Labels:** `milestone:M3`, `stage:F`, `physics`, `parallel-track`
