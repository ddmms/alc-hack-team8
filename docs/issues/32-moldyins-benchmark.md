# Convert the MolDyINS fixture and benchmark Method 2 against it

## Context

MolDyINS is Paper B's released implementation (GPL-3.0), and its bundled 25 K trajectory is
a ready-made regression fixture (`design.md` §4 Tier 3). It is the only external check
available on the anisotropic method, which matters because `27-anisotropic-path.md`'s
rank-4 contraction (A8) has no internal bound at all.

Two obstacles, both known. The trajectory is a GROMACS `.trr`, which ASE cannot read (D2),
so it needs one-off conversion to extxyz or `.traj` — a developer-side step using
MDAnalysis or `gmx`, checked in as a fixture rather than added as a runtime dependency.
And MolDyINS itself uses `np.float`, removed in NumPy 1.24, so it will not run on a current
stack without edits (method-review §5). `plan.md` §6 is explicit that the repair effort is
genuinely unknown and must block nothing.

## Scope

- One-off conversion of the MolDyINS 25 K trajectory to a checked-in fixture, with the
  conversion script and its provenance recorded. Trim or subsample to stay near the ~1 MB
  fixture ceiling if the full trajectory is too large, and record what was trimmed.
- Our Method 2 spectrum computed from that fixture.
- A minimal, documented patch set to make MolDyINS run on a current NumPy, kept in
  `scripts/` and clearly marked as ours, not upstream's.
- Comparison of the two spectra, written up.
- An explicit statement of what agreement would and would not prove.

## Out of scope

- Upstreaming the MolDyINS repair. Worth offering, not worth blocking on.
- Method 1 — `31-oclimax-oracle.md`.
- Adding MDAnalysis as anything but a dev extra.

## Acceptance criteria

- [ ] The converted fixture reads through `read_velocities` and produces a physically
      sensible pDOS: the P3HT C–H stretch region is populated, the sum rule holds against
      the stated 25 K, and the spectrum is not aliased.
- [ ] The conversion is reproducible from the checked-in script given the original `.trr`,
      and the script records the original file's checksum so the fixture's ancestry is
      verifiable.
- [ ] The comparison is quantitative, per order where MolDyINS exposes orders.
- [ ] The write-up states plainly that agreement with MolDyINS on the rank-4 contraction
      proves only that we reproduced the same approximation, not that the approximation is
      good. A8 remains open either way, and the write-up must not be read as closing it.
- [ ] Any MolDyINS behaviour that contradicts the paper is recorded — the paper and the
      code are two sources and they may disagree.
- [ ] Nothing in the main test suite depends on this issue.

## Depends on

- `29-method-2-assembly.md`

## Open questions / risks

- **Effort is open-ended.** `plan.md` §6 already flags this. If the repair proves deep,
  the fallback is to use the converted trajectory as our own regression fixture — which is
  useful on its own — and abandon the head-to-head comparison.
- The trajectory is 13 MB as distributed; keeping the fixture small without destroying the
  frequency resolution needs thought. Subsampling in time raises the Nyquist problem
  directly, so trim in *duration* or in atom count, not in `dt`.

**Labels:** `milestone:M5`, `validation`, `open-question`
