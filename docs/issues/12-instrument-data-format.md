# Define the instrument data file format and loader

## Context

`design.md` D4 puts instrument definitions in data files rather than in code, so that
adding a spectrometer needs no release. That is a small decision with a large consequence
for the project's usefulness: the set of instruments is open-ended and mostly not our
expertise, and a contributor who knows MARI should be able to add MARI by writing a TOML
file and a test.

The licence change in D1 also changed what can go in those files. AbINS's calibrated
instrument parameters are the most valuable thing the project is declining to depend on,
and since AbINS is now licence-compatible those values can be copied with attribution
rather than re-derived.

## Scope

- A TOML schema under `src/mdins/instrument/data/` covering: instrument name, geometry
  (direct or indirect), fixed energy, detector banks with angles and weights, accessible
  energy range, resolution model and its coefficients, and a provenance block naming the
  source of every number and its licence attribution.
- A loader returning a validated `Instrument` object, with a clear error for an unknown
  key or a missing required field.
- Instrument lookup by case-insensitive name, and by path so a user can supply their own
  file without installing it.
- A minimal example file and a documented procedure for adding an instrument.
- At least one instrument beyond the TOSCA/VISION pair so the schema is exercised by more
  than one geometry — a direct-geometry chopper instrument such as MARI.

## Out of scope

- The TOSCA and VISION numbers themselves, and the disagreement in them —
  `13-tosca-resolution-constants.md`.
- Evaluating the resolution model — `14-resolution-broadening.md`.

## Acceptance criteria

- [ ] Every shipped data file loads, and a parametrised test instantiates each and
      computes `Q(ω)` across its stated range without error. A new file added without a
      test still gets covered by the parametrisation, which is the point.
- [ ] A file with an unrecognised key fails loudly rather than being ignored. Silently
      dropping a misspelled `resolution_coefficients` would give a default-broadened
      spectrum that looks fine.
- [ ] A file missing its provenance block fails validation. This is a hard requirement,
      not a lint: `design.md` §3 makes provenance non-optional, and instrument constants
      are exactly the numbers whose origin needs to be recoverable years later.
- [ ] A user-supplied file outside the package loads by path.
- [ ] Round-trip: an `Instrument` written back to TOML and reloaded compares equal.

## Depends on

- `11-instrument-kinematics.md`

**Labels:** `milestone:M3`, `stage:F`, `infrastructure`, `parallel-track`
