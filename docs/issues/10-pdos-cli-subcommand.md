# Add the `mdins pdos` subcommand

## Context

`plan.md` §3 gives the CLI as `mdins pdos traj.extxyz -o sd.h5`, mirroring the library
entry point. `design.md` §1 notes that stages A–C are useful without the rest — an
atom-projected DOS from MD has value beyond INS and is the first deliverable — so this is
the point at which the project becomes usable by someone who is not importing it.

The CLI is deliberately thin. Every subcommand is a short call into the library, so
anything the CLI can do is reachable from Python with the same arguments, and the CLI
never becomes the only place a piece of logic lives.

## Scope

- A `pdos` subcommand taking a trajectory path, an output HDF5 path, and the estimator
  parameters: energy grid (max energy and bin count, or an explicit grid), estimator
  choice, segment length, overlap, window, and `T_MD`.
- Pre-processing flags mirroring stage A: COM removal on by default, angular-velocity
  removal opt-in.
- Write-path options from `04-hdf5-serialisation.md`: entity grouping, storage dtype.
- Provenance recording the command line that produced the file.
- A short textual summary to stderr on completion: number of atoms, grid, resolution,
  sum-rule residual per species.

## Out of scope

- The `spectrum` subcommand — `25-spectrum-cli-subcommand.md`.
- Plotting. The IR is the output; visualisation is a downstream concern.

## Acceptance criteria

- [ ] `mdins pdos` on a checked-in small extxyz fixture produces an HDF5 file that
      `VelocitySpectralDensity.from_hdf5` reads and that passes `check_sum_rule`.
- [ ] Every CLI option maps onto a named library argument, verified by a test that
      introspects both signatures. If the CLI gains an option the library cannot express,
      that test fails — which is the property that keeps the CLI thin.
- [ ] A sampling violation (dt too coarse for the requested max energy) exits non-zero
      with the stage A error message, rather than producing an aliased file.
- [ ] `--help` output names the units of every physical argument. A `--temperature 10`
      that could be °C is exactly the kind of ambiguity `design.md` §3 is guarding
      against.
- [ ] The recorded provenance reproduces the run: a test round-trips the stored command
      line through the parser and gets the same arguments back.

## Depends on

- `04-hdf5-serialisation.md`
- `07-welch-estimator.md`

**Labels:** `milestone:M1`, `stage:A`, `stage:C`, `cli`
