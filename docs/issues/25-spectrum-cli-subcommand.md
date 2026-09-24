# Add the `mdins spectrum` subcommand

## Context

`plan.md` §3 gives the second half of the CLI: `mdins spectrum sd.h5 --instrument tosca
--temperature 10`. It consumes an IR file rather than a trajectory, which is the visible
consequence of `design.md` §1's stage independence — stages D–G never touch a trajectory,
and anyone with phonon data from another source can enter the pipeline at C.

## Scope

- A `spectrum` subcommand reading a `VelocitySpectralDensity` from HDF5 and writing a
  spectrum.
- Options for method (`isotropic` / `anisotropic`), instrument (by name or path),
  `T_experiment`, `max_order`, the low-frequency cutoff, and the order prefactor
  convention.
- Order-resolved output: write the per-order decomposition, not only the total.
- Provenance chaining: the spectrum's provenance must include the IR's, so a spectrum
  file traces back to the trajectory that produced it through the intermediate file.
- A warning on stderr when `T_experiment` differs from the IR's `T_MD` (D8).

## Out of scope

- The anisotropic method's own physics — `29-method-2-assembly.md`. This issue exposes the
  flag; the flag may raise `NotImplementedError` until that lands.
- Plotting.

## Acceptance criteria

- [ ] `mdins pdos` followed by `mdins spectrum` on a checked-in fixture produces a
      spectrum file end to end, and the two-step result is identical to the equivalent
      in-process library call.
- [ ] The spectrum's provenance contains the trajectory path and the estimator settings
      from the IR, not just the second command. A test asserts the chain: if provenance
      only records the most recent step, the IR's role as a shareable artefact is lost.
- [ ] Every CLI option maps onto a named library argument, verified by signature
      introspection, as for `mdins pdos`.
- [ ] Omitting `--temperature` uses the IR's `T_MD` and says so in the summary output;
      supplying a different one warns.
- [ ] Reading an IR with an unrecognised `normalisation` field exits non-zero rather than
      assuming ours.
- [ ] `--help` names the units of every physical argument.

## Depends on

- `24-method-1-assembly.md`
- `10-pdos-cli-subcommand.md`

**Labels:** `milestone:M3`, `cli`
