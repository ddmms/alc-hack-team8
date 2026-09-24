# Assemble `isotropic_spectrum`, the Method 1 end-to-end path

## Context

The point at which stages D, E, F and G are wired into one call and the project delivers
what the proposal promised: an INS spectrum from an MD trajectory, by Paper A's method.
Everything it needs exists by now, so this issue is assembly and ordering, not new physics
— but the ordering is where the mistakes are.

`design.md` §2[E] and method-review §1 are explicit that convolution happens at atomic
level, **before** cross-section and Debye–Waller factor are applied. Applying the
cross-section first and convolving afterwards is a natural-looking refactor that produces
a subtly wrong multiphonon series, and nothing in the individual unit tests would catch
it.

## Scope

- `isotropic_spectrum(sd, instrument, temperature, max_order, cutoff, ...) -> Spectrum`,
  matching the signature advertised in `plan.md` §3.
- The pipeline: IR → `tr(u²)/3` → DWF → per-atom multiphonon convolution → isotropic
  powder average → Sears incoherent cross-section from Euphonic → `Q(ω)` from the
  instrument → resolution broadening → binning → `Spectrum`.
- Sears cross-section and isotope data pulled from Euphonic (D4), not tabulated here.
- Chunking over atoms, preserving the per-atom independence `design.md` §3 asks not to
  destroy gratuitously.
- Full provenance and metadata propagation from the IR through to the spectrum.

## Out of scope

- The anisotropic path — `29-method-2-assembly.md`.
- The CLI — `25-spectrum-cli-subcommand.md`.
- External benchmarks — M5.

## Acceptance criteria

- [ ] End to end on the Einstein-crystal fixture, the spectrum's fundamental sits at the
      known Einstein energy after broadening, and the order-`n` features sit at `n` times
      it. This is the first test that exercises the whole chain against an exact answer.
- [ ] A test asserts the operation order: cross-section and DWF are applied after
      convolution. Implement it by scaling one species' cross-section by a factor and
      checking that every order scales by exactly that factor and not by its power — the
      wrong ordering gives `factor^n` at order `n`, which is an unmistakable signature.
- [ ] Doubling an atom's scattering cross-section doubles its contribution to every order
      and leaves the peak positions unchanged.
- [ ] Running with two identical atoms gives twice the intensity of one, per order. Cheap,
      and it catches a normalisation applied per system instead of per atom.
- [ ] Chunked and unchunked runs agree to machine precision.
- [ ] The returned spectrum's metadata names the instrument, both temperatures, the
      cutoff, the prefactor convention and the maximum order.
- [ ] Peak memory on a several-hundred-atom system stays bounded by the chunk size, not by
      `n_atoms × n_energies × max_order`. Measure it; this is the wall `design.md` §2[C]
      says both source papers hit.

## Depends on

- `14-resolution-broadening.md`
- `21-multiphonon-convolution.md`
- `22-isotropic-powder-average.md`
- `23-spectrum-container.md`

**Labels:** `milestone:M3`, `stage:D`, `stage:E`, `stage:F`, `stage:G`, `physics`
