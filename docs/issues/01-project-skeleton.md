# Set up the package skeleton, packaging metadata and stub entry points

## Context

Nothing else can be reviewed until there is a package to put it in, and two of the
choices made here are expensive to reverse later: the distribution name, and the exact
GPL variant. `plan.md` §1 fixes the module layout; `design.md` D1 records the relicence
to GPL-3.0 and the housekeeping it left outstanding — the GPL `LICENSE` file is the
licence text only and carries no copyright line, so the `Copyright (c) 2026, Alin Marin
Elena` statement the old BSD file held has to be restated in per-file headers and in the
README.

Stub entry points are part of this issue rather than a later one because `plan.md` §3
depends on callers being writable against the public API before the internals exist.
That is what makes the parallel tracks in §4 possible at all.

## Scope

- `pyproject.toml` with hatchling, `requires-python >= 3.11`, `src/` layout, and
  `License-Expression` metadata.
- Runtime dependencies `numpy`, `scipy`, `ase`, `h5py`; optional extras `euphonic`
  (validation) and `mdanalysis` (dev-only fixture conversion). A comment recording why
  `abins` is not an extra (conda-only).
- `uv` lockfile and a documented dev workflow (`uv sync`, `uv run pytest`).
- `ruff` and `mypy` configuration, with mypy strict on `ir.py`, `units.py` and
  `spectral.py` at minimum.
- The empty module tree of `plan.md` §1, each file carrying an SPDX header and a module
  docstring saying which pipeline stage it implements.
- `src/mdins/py.typed`.
- Public entry points `velocity_spectral_density`, `isotropic_spectrum` and
  `anisotropic_spectrum`, exported from `mdins`, with the signatures in `plan.md` §3 and
  bodies raising `NotImplementedError`.
- A `mdins` console script whose `--version` works and whose subcommands are declared but
  not implemented.
- A README licence section restating the copyright holder.

## Out of scope

- Any physics. The stubs raise; they do not approximate.
- CI configuration — `05-continuous-integration.md`.
- The IR dataclass itself — `03-intermediate-representation.md`.
- User-facing documentation of the velocity constraint — `35-user-documentation.md`.

## Acceptance criteria

- [ ] `pip install -e .` succeeds on Python 3.11, 3.12 and 3.13 in a clean environment.
- [ ] `import mdins` pulls in no optional dependency. A test that stubs `euphonic` and
      `abins` out of `sys.modules` and then imports every module under `mdins` passes;
      it fails if anyone adds a top-level `import euphonic` to the main path, which
      `design.md` §3 forbids.
- [ ] `from mdins import velocity_spectral_density, isotropic_spectrum,
      anisotropic_spectrum` works, and calling each raises `NotImplementedError` rather
      than `TypeError` — i.e. the signatures accept the arguments `plan.md` §3 advertises.
- [ ] `mdins --version` prints the version from package metadata, not a hardcoded string.
- [ ] Every file under `src/mdins/` begins with
      `# SPDX-License-Identifier: GPL-3.0-or-later` (or `-only`, once decided). A test
      walks the package and asserts this, so a new module without a header fails CI
      rather than being noticed at release.
- [ ] `ruff check` and `ruff format --check` pass on an empty tree.

## Depends on

Nothing. This is the root of the graph.

## Open questions / risks

- **`GPL-3.0-only` or `GPL-3.0-or-later`.** `design.md` D1 flags this as not recoverable
  once distributed and does not settle it. Needs the copyright holder's answer before the
  first release, not before the first commit — but the metadata should not be guessed and
  then quietly changed.
- **The distribution name `mdins` is unverified on PyPI** (`plan.md` §1). Check before the
  first upload; renaming is cheap now and expensive after anyone has depended on it.

**Labels:** `milestone:M0`, `infrastructure`, `licensing`
