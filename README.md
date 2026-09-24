# mdins

Simulate inelastic neutron scattering (INS) spectra from molecular dynamics trajectories.

Two published methods route from an MD trajectory to a spectrum via a derived phonon
density of states, rather than working directly from the trajectory:

- **Isotropic** — Cheng, Kolesnikov & Ramirez-Cuesta, *J. Chem. Theory Comput.* **16**,
  7702 (2020). Uses the atom-projected DOS in place of the atomic displacement tensors.
- **Anisotropic** — Harrelson *et al.*, *Sci. Rep.* **11**, 7938 (2021). Retains the
  velocity cross-correlation terms, so vibrational anisotropy survives into the spectrum.

Both are implemented here as engine-neutral, tested components. The atom-projected DOS
from MD is useful in its own right and is available without the scattering machinery.

**Status: early development.** The intermediate representation and unit conventions are in
place; see [`plan.md`](plan.md) for what is implemented and what is not.

## Install

```sh
uv sync --extra euphonic
```

## Documentation

| | |
|---|---|
| [`proposal.md`](proposal.md) | why this exists and what it will do |
| [`design.md`](design.md) | architecture, and the decisions behind it |
| [`plan.md`](plan.md) | implementation and test plan |
| [`docs/method-review.md`](docs/method-review.md) | analysis of the source papers, and their gaps |

## Input requirements

Trajectories are read through ASE, and the pipeline needs **velocities**, not positions.
That rules out formats which do not store them (GROMACS `.xtc`, DCD) and formats ASE
cannot read (GROMACS `.trr`); convert to extxyz or ASE `.traj` first. The velocity dump
interval sets a hard ceiling on the resolvable energy, and modes above it alias back into
the spectrum rather than being filtered out, so this is validated as an error rather than
a warning.

## Licence

Copyright (C) 2026 Alin Marin Elena

This program is free software: you can redistribute it and/or modify it under the terms of
the GNU General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the [GNU General Public License](LICENSE) for more details.
