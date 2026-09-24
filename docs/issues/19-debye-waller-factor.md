# Implement the Debye–Waller factor and the two-temperature split

## Context

Also stage D. The Debye–Waller factor `A_i = ∫ u²_i(E) dE` with Bose–Einstein occupation,
and with it the decision `design.md` D8 records: `T_MD` and `T_experiment` are separate.

`T_MD` is a property of the IR — the temperature at which the PES was sampled, and what
appears in the classical-to-quantum conversion. `T_experiment` is an argument to stages D
and E, setting the occupation and the DWF. They default to equal.

Decoupling them is one of Paper A's real contributions: a fictive MD temperature can be
chosen to sample the right region of the PES — they run ice and ZIF-8 at 180 K to mimic
proton zero-point motion classically — while the spectrum is evaluated at the temperature
the measurement was actually made at. But it is an approximation with no error estimate
attached, so `design.md` D8 requires that using it be a deliberate act: opt-in, recorded in
provenance, visible in output metadata, never a silent default.

Note also method-review §3.1: Paper A's fictive-temperature trick only works if the
absolute normalisation is pinned down first, which is what D6 and the sum rule do.

## Scope

- `debye_waller(u_squared, temperature, ...)` returning the per-atom tensor `A_i`, with
  the isotropic `W_i = Q² tr(A_i)/6` form available as a derived quantity (Paper A Eq. 4).
- Bose–Einstein occupation `n(E, T)` in `units.py`-consistent terms, with the ω→0 limit
  handled rather than divided.
- The `T_MD` / `T_experiment` API: `temperature` defaults to `sd.temperature_md`, and
  supplying only one is the easy path.
- A provenance step and an output-metadata field recorded whenever the two differ, plus a
  warning naming both values.
- The BE-weighted DWF as the default; see `28-finite-temperature-and-dwf-variants.md` for
  Paper B's alternative.

## Out of scope

- Paper B's `max[...]` DWF heuristic (ambiguity A4) — that is Method 2's variant and lives
  in `28-finite-temperature-and-dwf-variants.md`.
- Applying the DWF in a scattering calculation — `24-method-1-assembly.md`.

## Acceptance criteria

- [ ] For a single oscillator, `tr(A)/3` reproduces the analytic
      `(ħ/2mω₀)·coth(ħω₀/2k_BT)`, checked at both the classical (`k_BT ≫ ħω₀`) and
      zero-point (`k_BT ≪ ħω₀`) limits. Getting only the classical limit right is
      consistent with having dropped the `½` and would be invisible at room temperature.
- [ ] `n(E, T)` matches the closed form and tends to `k_BT/E − ½` as `E → 0` without
      overflowing; a test evaluates it at `E = 0` and at 1e-12 meV.
- [ ] Passing no `temperature` gives exactly the same result as passing
      `sd.temperature_md`, to machine precision.
- [ ] Passing a different `temperature` emits a warning naming both values and writes the
      decoupling into provenance and into the output metadata. A test asserts all three;
      D8's whole point is that this cannot happen quietly.
- [ ] Raising `T_experiment` increases `tr(A)` monotonically. Cheap, and it catches a sign
      error in the occupation that the limit tests might not.
- [ ] The DWF is finite for any input that passed stage D's cutoff check, and a test feeds
      a spectrum with weight at the cutoff edge to confirm no divergence leaks through.

## Depends on

- `18-displacement-model.md`

**Labels:** `milestone:M3`, `stage:D`, `physics`
