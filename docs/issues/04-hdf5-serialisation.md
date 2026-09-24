# Implement HDF5 serialisation for the intermediate representation

## Context

`design.md` D3 splits the IR into two things deliberately: the dataclass is the schema,
and HDF5 is I/O. There is no separate schema document to drift out of sync with the code.
But the IR is meant to be shared and re-analysed by people who did not produce it, and
`design.md` §3 is blunt about the consequence — an unlabelled pDOS with an unknown
normalisation convention is worse than useless. So the file has to carry the metadata and
the provenance, not just the numbers.

This is also where the memory mitigations in `design.md` §2[C] live. 3000 atoms × 25000
bins × 6 components in float64 is ~3.6 GB, and both source papers hit exactly that wall.
Frequency rebinning, entity grouping and float32 storage are options on the dump path;
they are not changes to the in-memory contract.

## Scope

- `to_hdf5(path)` / `from_hdf5(path)` on `VelocitySpectralDensity`.
- Symmetric packing to six Voigt components **on write only**, unpacked on read.
- Metadata and provenance written as attributes, with a format-version attribute so a
  future reader can reject or migrate an old file rather than misread it.
- Write-path options: `dtype` (float32 for storage), entity grouping (sum over atoms
  sharing a label, with masses and counts handled consistently), and frequency rebinning
  onto a coarser grid.
- `mdins` CLI plumbing is not here, but the options must be expressible as keyword
  arguments so the CLI can expose them later.

## Out of scope

- Reading anyone else's format. Importing Euphonic or OCLIMAX output is not part of the
  IR contract.
- Golden-file regression comparison — `30-golden-file-regression.md`.

## Acceptance criteria

- [ ] Round-trip is **exact** for float64: `np.array_equal` on every array, and equality
      on every metadata and provenance field including the estimator parameters. A test
      that only checks `density` would pass while silently dropping the window name, which
      is precisely the information that makes a shared IR reproducible.
- [ ] Voigt packing round-trips a randomised symmetric tensor exactly, and a separate test
      asserts the off-diagonal components come back in the right slots — feed a tensor
      whose six components are 1..6 and check each position by hand. An off-by-one in the
      Voigt order produces a plausible tensor with `xy` and `xz` swapped, which nothing
      downstream would notice.
- [ ] Writing a non-symmetric density raises rather than silently discarding the
      antisymmetric part.
- [ ] float32 round-trip agrees to float32 precision and the file is close to half the
      size; the read-back object reports its storage dtype in metadata so a later
      tolerance failure can be attributed.
- [ ] Entity grouping conserves `∫ tr P dE` summed over all entities, and the grouped
      object still passes `check_sum_rule` against the grouped masses.
- [ ] Frequency rebinning conserves the integral to within the quadrature error of the
      coarser grid, and records the rebinning in provenance.
- [ ] Reading a file whose format version is newer than the reader raises a clear error
      naming both versions.

## Depends on

- `03-intermediate-representation.md`

**Labels:** `milestone:M0`, `stage:C`, `infrastructure`
