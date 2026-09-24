# Proposal: open tooling for simulating INS spectra from molecular dynamics

## Summary

We propose an open, well-tested Python implementation of two published methods for
simulating inelastic neutron scattering (INS) spectra from molecular dynamics (MD)
trajectories. Both methods work not from raw trajectories but from a derived
**phonon density-of-states (pDOS)**, used as an intermediate representation of the
system's dynamics.

Neither method is currently usable in practice by the wider community: one is implemented
only in a closed, binary-distributed code, the other only in an unmaintained script bound
to a single MD engine. We intend to make both available as engine-agnostic, validated,
reusable components.

## Why

### INS is uniquely well-suited to validating simulations, and underused for it

Neutrons scatter from nuclei. For non-magnetic excitations, atomic structure and dynamics
are sufficient to predict both peak positions *and* intensities — no electronic-structure
model is required, and there are no optical selection rules, so every mode is visible.
That makes INS one of the most stringent available tests of an interatomic potential, and
an unusually direct point of contact between experiment and simulation.

The conventional route to a simulated spectrum — density functional theory plus lattice
dynamics — is restricted to small, ordered, periodic cells. Most materials of real
interest are not: they are disordered, amorphous, defective, semicrystalline, or simply
too large. For those systems MD is the standard tool, and the bridge from an MD trajectory
to a comparable INS spectrum is the missing piece.

### The methods exist but the tooling does not

Two papers describe that bridge:

- **Cheng, Kolesnikov & Ramirez-Cuesta (2020)**, *JCTC* **16**, 7702.
  Uses the atom-projected pDOS in place of the atomic displacement tensors, within a fully
  isotropic approximation. Implemented in OCLIMAX — **closed source**, distributed as a
  binary Docker image.
- **Harrelson, Dettmann, Scherer, Andrienko, Moulé & Faller (2021)**, *Sci. Rep.* **11**, 7938.
  Goes further, retaining the cross-correlation terms in the energy-dependent input, so
  vibrational anisotropy survives into the spectrum. Implemented in
  [MolDyINS](https://github.com/tfharrelson/MolDyINS) — **open, but GROMACS-only,
  unmaintained since 2021, and incompatible with current NumPy**.

So a researcher today must either accept a black box they cannot inspect, extend or cite
reproducibly, or adopt a script tied to one MD engine that no longer runs.

### The surrounding ecosystem is ready, and has a hole in exactly this shape

Mature, actively maintained open tools already cover the neighbouring problems: Euphonic
and Mantid's AbINS handle phonon data and instrument modelling; ASE, MDAnalysis and
phonopy handle structures, trajectories and force constants. What is missing is the step
that turns MD-derived dynamics into the phonon-like quantities those tools consume.

Notably, **every AbINS input loader is an ab initio phonon code** — there is no route in
from MD. Filling that gap makes an existing, well-supported instrument-modelling stack
available to the large class of systems that only MD can reach.

### Why the pDOS intermediate matters

Routing through a derived pDOS rather than working directly from trajectories is a
deliberate choice, and is what both papers do. It gives a compact, physically meaningful,
inspectable intermediate that:

- decouples the MD engine from the spectrum calculation entirely;
- is independently checkable against harmonic reference calculations;
- is the same class of object that Euphonic and AbINS already work with, so it composes
  with existing tooling rather than duplicating it;
- is orders of magnitude smaller than the trajectory it came from, making results
  shareable and re-analysable without redistributing terabytes of velocities.

## What

Three capabilities, in dependency order. Each is independently useful and independently
testable.

### 1. Atom-projected density of states from trajectory data

Compute the atom-resolved pDOS from MD velocity data, in a form suitable as input to the
INS calculations that follow. Input goes through ASE, so any engine whose output ASE reads
— or which can be converted to extxyz — is supported, with no coupling between the MD code
and the spectrum calculation.

We will build on well-maintained external libraries where they fit a modern Python
ecosystem, rather than reimplementing trajectory I/O or correlation machinery that already
exists in good form.

This capability stands alone — atom-projected DOS from MD is broadly useful beyond INS.

**Validation.** Compute the pDOS for a system via this route, and independently via a
harmonic calculation in Euphonic, using a *consistent force calculator* on both sides —
the ASE Lennard-Jones implementation is a natural choice, being simple, exactly
specifiable, and free of force-field ambiguity. In the harmonic limit the two must agree.
This is a genuine end-to-end check of the correlation, normalisation and unit handling,
which is where this class of calculation most often goes quietly wrong.

### 2. INS intensity calculation — isotropic method (Cheng et al.)

Convert the atom-projected pDOS into a simulated INS spectrum under the fully isotropic
approximation: quantum displacement amplitudes, Debye–Waller factor, higher-order
(multiphonon) excitations, powder averaging, and the instrument's energy/momentum-transfer
trajectory and resolution.

### 3. INS intensity calculation — anisotropic method (Harrelson et al.)

Extend the energy-dependent input to retain cross-correlation terms, so that anisotropy of
atomic motion is carried through to the powder average rather than averaged away at the
outset. This is the same pipeline with a richer intermediate, which is why it follows
naturally from (2) rather than duplicating it.

### Benchmarking

Beyond the harmonic validation of (1), we will benchmark the intensity calculations
against existing codes — OCLIMAX and MolDyINS for the methods themselves, AbINS and
Euphonic for the shared phonon and instrument machinery — and against published
experimental spectra for systems used in the source papers.

Establishing agreement with a closed reference implementation is itself a contribution:
it converts an unverifiable black box into a result the community can reproduce.

## Scope

**In scope**

- The two published methods, implemented faithfully and documented honestly, including
  where the papers are ambiguous or appear to contain errors.
- Trajectory input via ASE, and a clearly specified pDOS intermediate that is independent
  of which MD engine produced the data.
- Validation against harmonic reference calculations, existing codes, and published
  experimental data.

**Not in scope, at least initially**

- **Coherent scattering from MD.** Absent from both source papers. It is a different and
  substantially more expensive calculation, and would change the shape of the work.
- **New physics.** We are implementing published methods, not extending them. The
  anharmonic-correction scheme derived but never used in the second paper is noted and
  deliberately left alone.
- **Running the MD.** We consume trajectories; we do not produce them.
- **A direct trajectory-to-spectrum route.** The intermediate scattering function can be
  computed without going via a pDOS, and other tools do so. Using the pDOS as the
  intermediate is the defining feature of these two methods and of this work.

## What success looks like

- pDOS from MD agrees with the Euphonic harmonic reference, under matched force
  calculators, to within a stated tolerance.
- Simulated spectra reproduce those from the reference implementations for the systems
  published in the source papers.
- Someone with an MD trajectory from an engine none of us has tested can produce a
  credible INS spectrum without reading either paper.
- Velocities, not just positions, are recognised as a hard input requirement up front
  rather than discovered after a run has finished.
- The ambiguities we had to resolve are written down, so the next person does not have to
  rediscover them.

## Known risks and open questions

- **The papers under-specify their own methods.** Several quantities needed for an exact
  reimplementation — most importantly the absolute intensity normalisation — are not
  unambiguously determined by what is published. Some equations appear to contain
  typographical errors. Resolving these against reference implementations and published
  spectra is real project work, not a preliminary.
- **Instrument parameters are partly unpublished**, and published resolution constants
  disagree between sources. Contact with instrument scientists may be needed.
- **Validation depth is bounded by reference availability.** OCLIMAX is usable as a
  numerical oracle but cannot be inspected; MolDyINS requires repair before it will run.
- **Trajectory sampling is a correctness constraint, not a performance knob.** The
  velocity dump interval sets a hard ceiling on the resolvable energy, and modes above it
  alias back into the spectrum rather than being filtered out. Input validation matters
  more here than usual.

A detailed review of both papers, the literature they depend on, and the specific gaps
identified so far is in [`docs/method-review.md`](docs/method-review.md).
