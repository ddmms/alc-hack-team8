# Implement canonical units and enforce the meV boundary rule

## Context

`plan.md` §2 names this the first of two conventions to fix before any code is written,
and it is not a stylistic preference. The failure mode this domain is notorious for is a
spectrum that looks entirely plausible and is wrong by 2π, or by a mass, or by ħ — and
which no amount of downstream testing localises, because every stage is consistent with
the one before it. `design.md` §3 explains why the defence is a single conversion module
plus tests rather than a unit library: `pint` was tried and rejected because the arrays
that matter go straight to `numpy.linalg.eigvalsh` and to HDF5, neither of which carries
a quantity, so units would be stripped at exactly the boundary they were meant to guard.

The rule that actually buys safety is structural: **energies in meV cross every function
boundary**. Angular frequency in rad/ps may exist inside a function body and must never
appear in a signature, a dataclass field, or a file. That makes most 2π errors
unrepresentable instead of merely testable.

## Scope

- `src/mdins/units.py` holding the canonical system (Å, ps, amu, meV) and nothing else.
- Physical constants derived from `scipy.constants` rather than typed in: `HBAR`, `KB`,
  `PLANCK`, and the neutron mass grouping needed by the kinematics.
- Conversions: energy ↔ frequency (THz), energy ↔ angular frequency, energy ↔
  wavenumber, and the amu·Å²·ps⁻² ↔ meV energy conversion.
- `thermal_velocity_squared(temperature, masses) -> 3 k_B T / m` in Å²·ps⁻², the
  right-hand side of the equipartition sum rule, so that every caller of the sum rule
  gets it from one place.
- `__all__` listing the full public surface.

## Out of scope

- Neutron cross-sections and scattering lengths — those come from Euphonic
  (`design.md` D4), not from here.
- Instrument-specific constants — `12-instrument-data-format.md`.
- Any conversion that only one caller needs; those belong with the caller.

## Acceptance criteria

- [ ] Round-trip tests: `wavenumber_to_energy(energy_to_wavenumber(E)) == E` to machine
      precision, and likewise for frequency and angular frequency, over a randomised
      range spanning 0.1–1000 meV.
- [ ] Anchor values are checked against published numbers, not against the module's own
      arithmetic: 1 meV = 8.06554 cm⁻¹ = 0.2417989 THz, and k_B = 0.0861733 meV/K, each
      to the precision printed. This is the test that fails if someone "simplifies" a
      constant.
- [ ] `energy_to_angular_frequency` differs from `energy_to_frequency` by exactly 2π. A
      test asserts the ratio, because the whole point of the module is that this factor
      lives in one identifiable place.
- [ ] `thermal_velocity_squared` reproduces a hand-computed value for hydrogen at 300 K.
- [ ] A test greps the public signatures of every module in `mdins` for parameter and
      field names matching `omega`, `angular`, or `rad_per_ps`, and fails if any is
      found. Without this the boundary rule is a convention that decays; with it, the
      first violation is a red CI run.
- [ ] `units.py` imports nothing from `mdins`, so it can never participate in a cycle.

## Depends on

- `01-project-skeleton.md`

**Labels:** `milestone:M0`, `physics`, `infrastructure`
