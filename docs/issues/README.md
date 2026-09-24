# Issue set

The text of the GitHub issues that take this project from an empty repository to a
complete implementation of both methods, in dependency order. Each file is one issue: a
plausible unit of work, with its scope bounded against its neighbours so the set does not
overlap. They are planning documentation, not status tracking — nothing here records what
has or has not been built. The numbering is the dependency order, so a directory listing
reads as a build sequence, and each issue states its dependencies by filename rather than
relying on that ordering alone. The architecture they implement is in
[`design.md`](../../design.md); the work breakdown and test strategy they expand on is in
[`plan.md`](../../plan.md); the source-paper discrepancies several of them have to confront
are in [`method-review.md`](../method-review.md).

## The issues

| # | Title | Milestone | Depends on |
|---|---|---|---|
| 01 | [Set up the package skeleton, packaging metadata and stub entry points](01-project-skeleton.md) | M0 | — |
| 02 | [Implement canonical units and enforce the meV boundary rule](02-units-and-mev-boundary.md) | M0 | 01 |
| 03 | [Freeze the VelocitySpectralDensity intermediate representation](03-intermediate-representation.md) | M0 | 02 |
| 04 | [Implement HDF5 serialisation for the intermediate representation](04-hdf5-serialisation.md) | M0 | 03 |
| 05 | [Set up continuous integration and the fast/slow test split](05-continuous-integration.md) | M0 | 01 |
| 06 | [Implement the ASE velocity trajectory reader, pre-processing and sampling validation](06-trajectory-reader.md) | M1 | 03 |
| 07 | [Implement the Welch cross-spectral density estimator](07-welch-estimator.md) | M1 | 03, 06 |
| 08 | [Implement the VACF (Blackman–Tukey) estimator](08-vacf-estimator.md) | M1 | 07 |
| 09 | [Build the Layer 1 analytic test suite and assert the IR invariants](09-analytic-test-suite.md) | M1 | 04, 07, 08 |
| 10 | [Add the `mdins pdos` subcommand](10-pdos-cli-subcommand.md) | M1 | 04, 07 |
| 11 | [Implement the indirect- and direct-geometry kinematics](11-instrument-kinematics.md) | M3 | 02 |
| 12 | [Define the instrument data file format and loader](12-instrument-data-format.md) | M3 | 11 |
| 13 | [Resolve the TOSCA and VISION resolution constants](13-tosca-resolution-constants.md) | M3 | 12, 14 |
| 14 | [Implement energy-dependent resolution broadening](14-resolution-broadening.md) | M3 | 12 |
| 15 | [Add the Einstein-crystal end-to-end test](15-einstein-crystal.md) | M2 | 09 |
| 16 | [Build the Euphonic harmonic reference harness](16-euphonic-reference-harness.md) | M2 | 03, 05 |
| 17 | [Compare LJ pDOS against Euphonic and write up the harmonic validation](17-euphonic-comparison.md) | M2 | 15, 16 |
| 18 | [Implement the displacement model and the low-frequency cutoff](18-displacement-model.md) | M3 | 03, 09 |
| 19 | [Implement the Debye–Waller factor and the two-temperature split](19-debye-waller-factor.md) | M3 | 18 |
| 20 | [Write the closed-form multiphonon oracle test](20-multiphonon-oracle.md) | M3 | 19 |
| 21 | [Implement iterative multiphonon convolution and the order prefactor](21-multiphonon-convolution.md) | M3 | 20 |
| 22 | [Implement isotropic powder averaging](22-isotropic-powder-average.md) | M3 | 19 |
| 23 | [Implement the spectrum output container](23-spectrum-container.md) | M3 | 03 |
| 24 | [Assemble `isotropic_spectrum`, the Method 1 end-to-end path](24-method-1-assembly.md) | M3 | 14, 21, 22, 23 |
| 25 | [Add the `mdins spectrum` subcommand](25-spectrum-cli-subcommand.md) | M3 | 10, 24 |
| 26 | [Implement the almost-isotropic powder average](26-almost-isotropic-powder-average.md) | M4 | 22 |
| 27 | [Implement the anisotropic tensor path and its rank-4 contraction](27-anisotropic-path.md) | M4 | 21, 26 |
| 28 | [Implement finite-temperature Boltzmann weighting and the Debye–Waller variants](28-finite-temperature-and-dwf-variants.md) | M4 | 27 |
| 29 | [Assemble `anisotropic_spectrum`, the Method 2 end-to-end path](29-method-2-assembly.md) | M4 | 24, 27, 28 |
| 30 | [Add the golden-file regression layer](30-golden-file-regression.md) | M5 | 24 |
| 31 | [Build the OCLIMAX oracle harness for Method 1](31-oclimax-oracle.md) | M5 | 24 |
| 32 | [Convert the MolDyINS fixture and benchmark Method 2 against it](32-moldyins-benchmark.md) | M5 | 29 |
| 33 | [Cross-check broadening against AbINS](33-abins-broadening-crosscheck.md) | M5 | 13, 14 |
| 34 | [Reproduce published spectra and write up the comparison](34-published-spectra.md) | M5 | 13, 29 |
| 35 | [Document the velocity-input constraint and the user workflow](35-user-documentation.md) | M3 | 25 |
| 36 | [Maintain the ambiguity resolution log](36-ambiguity-resolution-log.md) | — | 03 |

Milestones are those in `plan.md` §4. Issues 11–14 carry `milestone:M3` because that is
when their output is first needed, but `plan.md` §4 lists the instrument layer as
independent work that can be done any time after M0 — see below.

## The critical path

```
01 → 02 → 03 → 04 → 07 → 09 → 15 → 17        (M0–M2, the committed cut)
                              ↘ 18 → 19 → 20 → 21 → 24   (M3)
                                              ↘ 27 → 29  (M4)
```

Stated in words: the units module, then the IR, then HDF5 round-tripping, then the Welch
estimator, then the Layer 1 analytic suite. Those five are strictly sequential and
everything else in the project sits behind them. From there the path runs through the
Einstein crystal to the Euphonic comparison — **issue 17 is the highest-risk item in the
plan** (`plan.md` §4 marks M2 "highest risk") and the point at which the pipeline is
either right or not. The M3 branch then runs displacement → Debye–Waller → multiphonon
oracle → convolution → Method 1 assembly, and M4 hangs off that.

Two things off the critical path are worth starting early, because `plan.md` §4 says they
barely touch it:

- **The instrument layer, issues 11–14.** Completely self-contained: it takes an energy
  grid and returns `Q(ω)` and a broadening kernel, with no MD and no IR involved. It can
  be built and tested against published TOSCA curves the day after issue 02 lands, and it
  only rejoins the path at issue 24.
- **The validation harness, issue 16.** The Euphonic reference side depends on the frozen
  IR but not on there being a real implementation behind it, so it can be built in
  parallel with issues 06–09.

Likewise, issues 18–23 depend on the IR *existing*, not on it being *correct* — they can be
developed against synthetic IR instances as soon as issue 03 lands, which is the practical
payoff of `plan.md` §3's instruction to write the IR first.

Three ordering constraints are deliberate and should not be reversed:

- **09 before anything with real MD.** `plan.md` §6 makes the analytic layer the mitigation
  for silent factor errors, and it only works if it is written first.
- **15 before 17.** The Einstein crystal has an exact answer and no dispersion, so it
  isolates our bugs from anharmonicity and q-sampling before LJ introduces both. Debugging
  17 without 15 passing means three candidate explanations and no way to separate them.
- **20 before 21.** `design.md` D7 removed the internal cross-check on the multiphonon
  chain, so the closed-form oracle is what *selects* the order prefactor rather than
  confirming a choice already made in the code.

M5 (issues 30–34) blocks nothing and is explicitly not part of the committed cut.
`plan.md` §4: "M0–M2 is the honest first cut… Treating M5 as in-scope for a short effort is
how this ends up with a spectrum nobody has verified."

## Issues with unresolved questions

Several issues depend on questions the literature leaves open. Each carries an
`open-question` label and an **Open questions / risks** section; issue 36 is where their
eventual resolutions are recorded.

| Issue | Question |
|---|---|
| 13 | The TOSCA resolution constants disagree between OCLIMAX and AbINS by a factor of ten (method-review §3.7), and VISION's own parameters are unpublished |
| 18 | The low-frequency cutoff (A6) has a required shape but no recommended value, and the implicit quantum correction factor (method-review §3.6) is a choice both papers make differently |
| 20, 21 | The multiphonon order prefactor (A7): Paper B asserts `1/n!`, aCLIMAX's exact prefactors differ. Also detailed balance through the convolution, which method-review §3.5 records as undescribed |
| 27 | The rank-4 → rank-2 contraction (A8) is an uncontrolled approximation whose error cannot be bounded internally once D7 ruled out the exponential route |
| 28 | Paper B Eq. 19 diverges as printed and is treated as a typo (A3); the Eq. 21 Debye–Waller heuristic has no derivation (A4) |
| 31 | Whether OCLIMAX falls back to an isotropic Debye–Waller factor given only `tr(B)` (A5) is unconfirmed anywhere |
| 01 | `GPL-3.0-only` versus `GPL-3.0-or-later` is not recoverable once distributed, and `design.md` D1 leaves it open |
