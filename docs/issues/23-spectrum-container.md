# Implement the spectrum output container

## Context

Stage G (`design.md` §2[G]). A small deliverable with two non-negotiable properties: the
per-order decomposition is retained, because both papers show order-resolved plots and
`design.md` calls the decomposition diagnostically essential; and the output carries units
and metadata, because `design.md` §3 makes provenance non-optional and an unlabelled
spectrum is as useless as an unlabelled pDOS.

It also needs to export to something the neutron ecosystem already reads, or the result of
the whole pipeline is a numpy array that nobody can open next to their data.

## Scope

- `src/mdins/spectrum.py` with a frozen container holding the energy grid, the total
  intensity, the per-order and per-species decompositions, the `Q(ω)` trajectory used, the
  instrument identity, both temperatures, and a `Provenance` record.
- Derived views: sum over orders, sum over species, a single order, a single species.
- Export to at least one ecosystem format — a Euphonic `Spectrum1D`/`Spectrum1DCollection`
  is the natural target given D4 already makes Euphonic a core dependency — plus plain
  text with a metadata header for the case where someone just wants the numbers.
- Metadata fields surfacing the things `design.md` requires to be visible: the
  normalisation, the low-frequency cutoff, the order prefactor convention, and whether
  `T_MD` and `T_experiment` were decoupled (D8).

## Out of scope

- Producing one — `24-method-1-assembly.md`.
- Plotting.
- Reading other codes' spectra; comparison scripts convert on their own terms.

## Acceptance criteria

- [ ] Summing the per-order arrays reproduces the total exactly, and a test asserts it
      rather than trusting the producer. If these drift apart, an order-resolved figure
      and the total in the same paper disagree.
- [ ] Likewise for the per-species decomposition.
- [ ] Export round-trips the energy grid and intensities without unit change, verified by
      re-importing and comparing.
- [ ] The metadata contains the four items listed above; a test asserts each key is
      present and non-null. These are the parameters a reader needs to reproduce or
      distrust the result, and they are the ones easiest to forget to propagate.
- [ ] Constructing a spectrum whose per-order arrays do not match the declared orders
      raises at construction time, not at plot time.
- [ ] The container is frozen and its arrays are not aliases of the producer's internals.

## Depends on

- `03-intermediate-representation.md`

**Labels:** `milestone:M3`, `stage:G`, `architecture`
