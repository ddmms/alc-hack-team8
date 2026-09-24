# Implement the ASE velocity trajectory reader, pre-processing and sampling validation

## Context

Stage A, per `design.md` §2[A]. Two things make this more than a file reader.

The first is that the input is **velocities, not positions**, and combined with the
ASE-only decision (D2) that is the tightest constraint in the whole design. ASE `.traj`,
extxyz with a velocity or momenta array, and LAMMPS dumps that requested velocities all
work; GROMACS `.trr` has no native ASE reader, and `.xtc` and DCD do not store velocities
at all. Users of other engines convert first. `design.md` is explicit that this cost
should be stated plainly rather than discovered, so the error message when a file has no
velocities is part of the deliverable, not an afterthought.

The second is sampling validation. Aliasing is silent and unrecoverable — modes above
Nyquist fold back and contaminate the spectrum with no signature that anything went wrong
(method-review §4, "Sampling requirements"). `design.md` §2[A] therefore makes this an
error, not a warning.

## Scope

- `src/mdins/trajectory.py` with a `VelocityTrajectory` container holding velocities,
  masses, symbols, `dt`, cell, and a `Provenance` record.
- `read_velocities(path, ...)` dispatching through `ase.io`, with an explicit check that
  velocities are present and an error naming the format and the likely fix when they are
  not.
- Chunked/streaming access, so stage B can consume a trajectory it cannot hold in memory
  (`design.md` §3, "Streaming").
- Pre-processing in the order `design.md` §2[A] gives: remove centre-of-mass velocity per
  frame; optionally remove global angular velocity, meaningful only for non-periodic
  systems. Each applied step is appended to provenance.
- `validate_sampling(dt, n_frames, max_energy, energy_resolution)` asserting
  `dt < 1/(2 c ν_max)` and `t_total > 1/(c Δν)`, raising with the actual and required
  values.
- A thermostat advisory: warn if provenance or user input indicates a strongly coupled
  thermostat or a barostat (method-review §4.7 — τ = 1 ps corrupts the spectrum, barostats
  damp the dynamics).

## Out of scope

- Deriving velocities by finite-differencing positions. `design.md` §2[A] calls this a
  fallback, not a default: it needs a `sinc` transfer function divided out and periodic-
  image unwrapping. If it is ever implemented it gets its own issue and must be loudly
  flagged in provenance.
- Format conversion of external fixtures — `32-moldyins-benchmark.md`.
- Spectral estimation — `07-welch-estimator.md`.

## Acceptance criteria

- [ ] A round-trip test writes a small extxyz with a velocity array via ASE, reads it
      back, and recovers the velocities to machine precision including the unit
      conversion from ASE's Å/(Å√(amu/eV)) convention to Å/ps. Getting that conversion
      wrong rescales every spectrum by a constant and is invisible without this test.
- [ ] Reading a positions-only file raises with a message that names the format and says
      what to do; a test asserts on the message text, because this is the error most new
      users will hit.
- [ ] COM removal drives `Σ m_i v_i` to zero to machine precision on a trajectory with a
      deliberately injected drift, and the injected drift shows up as a spurious ω→0
      feature when removal is skipped — assert both directions, so the test proves the
      step does something.
- [ ] Angular-velocity removal conserves total kinetic energy minus the rotational part
      on a rigid rotating cluster, and refuses (or warns) on a periodic cell.
- [ ] `validate_sampling` raises when `dt` is too coarse for the requested `max_energy`,
      and the error states both the supplied `dt` and the largest admissible one. A test
      constructs a trajectory sampling a 200 meV mode at 10 fs and asserts the failure.
- [ ] Streaming a trajectory in chunks yields the same concatenated velocities as reading
      it whole.
- [ ] Provenance after reading and pre-processing lists the steps applied, in order.

## Depends on

- `03-intermediate-representation.md`

## Open questions / risks

- The thermostat advisory can only act on what the file tells us, and extxyz usually tells
  us nothing. It may end up as a documentation item plus an optional user-supplied
  declaration rather than a real check; decide during implementation and do not pretend to
  a detection that is not happening.

**Labels:** `milestone:M1`, `stage:A`, `physics`
