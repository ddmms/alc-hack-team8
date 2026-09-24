# Reproduce published spectra and write up the comparison

## Context

`design.md` §4 Tier 3 and `plan.md` §5 Layer 5. Comparison against measured spectra is the
only test that asks whether the whole chain — MD, quantum conversion, multiphonon,
instrument — produces something a neutron scientist would recognise. Everything before this
checks internal consistency or agreement with another simulation.

`plan.md` §7's M5 definition of done asks for a written-up comparison against at least one
external source "including where it disagrees". These are figures in docs, not assertions,
because a disagreement with experiment is information, not a regression.

## Scope

- At least one system from each of the two papers' validation sets, chosen for what it
  exercises:
  - **Ice Ih** (Paper A, VISION and SEQUOIA) — the fictive-temperature trick's headline
    case, so it exercises D8's decoupling on the system it was designed for.
  - **MgH₂** (Paper A′) — strong multiphonon, so it exercises the convolution chain hardest
    and is the best available experimental check on the order prefactor.
  - **P3HT** (Paper B) — anisotropic, and the system MolDyINS was written for.
- MD setup scripts for each, with the potential, thermostat and sampling parameters
  recorded. Note method-review §4.7: production runs need NVE or a weak decorrelating
  thermostat; τ = 1 ps corrupts the spectrum and barostats damp the dynamics.
- Digitised or sourced experimental reference data where licensing permits, or a clear
  pointer to it where it does not.
- A documentation page per system with the overlay, the residual, and a paragraph on what
  agrees and what does not.

## Out of scope

- Turning any of this into a CI assertion.
- Improving the force fields. If the MD is wrong, the spectrum is wrong, and saying so is
  a legitimate outcome.
- Coherent scattering. Absent from both papers and out of scope for the project
  (method-review §3.8); graphite and other coherent cases are not candidates here.

## Acceptance criteria

- [ ] Each comparison states which parts of the disagreement are attributable to the MD
      (force field, sampling, temperature) and which to our implementation. An overlay
      with no attribution is a picture, not a validation.
- [ ] The MgH₂ case reports the relative weight of each multiphonon order against the
      measured spectrum, since that is the experimentally visible consequence of the A7
      prefactor choice and the only real-world evidence available on it.
- [ ] The ice Ih case runs with `T_MD ≠ T_experiment` and reports what the decoupling
      bought, so D8's approximation is exercised on its intended use case rather than only
      supported in the API.
- [ ] Every MD setup is reproducible from a checked-in script, with the trajectory itself
      not checked in.
- [ ] The write-up is honest about failures. `plan.md` §4 warns that treating M5 as
      in-scope for a short effort is how this ends up with a spectrum nobody has verified;
      a comparison that quietly omits the system that did not work reproduces exactly that
      failure.

## Depends on

- `29-method-2-assembly.md`
- `13-tosca-resolution-constants.md`

## Open questions / risks

- Genuinely unbounded (`plan.md` §4 marks M5 "M, unbounded"). Pick one system, finish it,
  write it up, and only then consider a second.
- Experimental data availability and licensing varies by source; some spectra may only be
  referenceable by figure, which limits the comparison to qualitative.

**Labels:** `milestone:M5`, `validation`, `docs`
