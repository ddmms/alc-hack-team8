# SPDX-License-Identifier: GPL-3.0-or-later
"""Provenance records.

Every artefact records what produced it. The intermediate representation is meant to be
shared and re-analysed, and an unlabelled spectral density with an unknown normalisation
convention is worse than useless.
"""

from __future__ import annotations

import json
import platform
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version


def _mdins_version() -> str:
    try:
        return version("mdins")
    except PackageNotFoundError:  # pragma: no cover - editable/source checkouts
        return "unknown"


@dataclass(frozen=True)
class Provenance:
    """What produced an artefact, and under what conditions.

    Attributes:
        source: Identity of the input, typically a trajectory path.
        mdins_version: Version of this package.
        created: ISO 8601 timestamp, UTC.
        python_version: Interpreter version.
        steps: Ordered record of processing applied, e.g. centre-of-mass removal.
        notes: Free-form entries. Anything approximate or unusual is recorded here, so
            that approximations cannot be applied silently.
    """

    source: str
    mdins_version: str = field(default_factory=_mdins_version)
    created: str = field(
        default_factory=lambda: datetime.now(UTC).isoformat(timespec="seconds")
    )
    python_version: str = field(default_factory=platform.python_version)
    steps: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    def with_step(self, step: str) -> Provenance:
        """Return a copy with ``step`` appended to the processing record."""
        return Provenance(
            source=self.source,
            mdins_version=self.mdins_version,
            created=self.created,
            python_version=self.python_version,
            steps=(*self.steps, step),
            notes=self.notes,
        )

    def with_note(self, note: str) -> Provenance:
        """Return a copy with ``note`` appended."""
        return Provenance(
            source=self.source,
            mdins_version=self.mdins_version,
            created=self.created,
            python_version=self.python_version,
            steps=self.steps,
            notes=(*self.notes, note),
        )

    def to_json(self) -> str:
        """Serialise to a JSON string, as stored in HDF5 attributes."""
        return json.dumps(asdict(self), sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> Provenance:
        """Inverse of :meth:`to_json`."""
        raw = json.loads(text)
        raw["steps"] = tuple(raw.get("steps", ()))
        raw["notes"] = tuple(raw.get("notes", ()))
        return cls(**raw)
