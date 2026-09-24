# Set up continuous integration and the fast/slow test split

## Context

`plan.md` §1 asks for GitHub Actions on Python 3.11–3.13, with ruff, mypy and the fast
tests on every PR and the slow and reference tests nightly. The split matters more than
usual here because of the shape of the test strategy in `plan.md` §5: Layer 1 and Layer 2
run in milliseconds and catch essentially every unit and normalisation error, while Layer
3 runs real MD and Layer 5 needs Docker. If they all live in one job, the fast layers stop
being run on every change, which defeats the reason they exist.

## Scope

- A PR workflow: matrix over Python 3.11–3.13, `ruff check`, `ruff format --check`,
  `mypy`, and `pytest -m "not slow"`.
- A nightly workflow running the full suite including `slow` and the Euphonic-dependent
  reference tests, on one Python version.
- `pytest` markers registered in `pyproject.toml`: `slow` (real MD, seconds to minutes)
  and `reference` (needs an optional dependency such as `euphonic`).
- A `tests/conftest.py` providing the shared fixtures and a skip mechanism for tests whose
  optional dependency is absent, so a contributor without `euphonic` gets skips rather
  than errors.
- Dependency installation via `uv` with the lockfile, so CI and local runs resolve to the
  same versions.
- Coverage reporting, reported but not gated — a coverage gate on a package this
  physics-heavy rewards testing the easy modules.

## Out of scope

- The tests themselves. This issue provides the harness; every later issue supplies the
  content.
- The Docker-based OCLIMAX harness, which `plan.md` §5 explicitly says cannot run in CI —
  `31-oclimax-oracle.md`.

## Acceptance criteria

- [ ] A PR run completes with the fast suite in under a minute on an empty-ish tree, and
      an unmarked test that takes minutes is caught by a `--durations` report in the job
      log rather than by someone noticing CI got slow.
- [ ] `pytest -m "not slow"` and `pytest -m slow` between them run every collected test.
      A test asserting the two sets partition the suite catches a marker typo, which
      otherwise silently removes a test from both jobs.
- [ ] An environment without `euphonic` installed produces skips with a reason naming the
      missing extra, not import errors.
- [ ] mypy runs in strict mode on `ir.py`, `units.py` and `spectral.py`, and non-strict
      elsewhere, with the configuration in `pyproject.toml` rather than command-line flags
      so local runs match CI.
- [ ] The nightly workflow is dispatchable manually, so a reference failure can be
      reproduced without waiting a day.

## Depends on

- `01-project-skeleton.md`

**Labels:** `milestone:M0`, `infrastructure`
