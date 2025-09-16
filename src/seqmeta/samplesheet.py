"""Utilities for parsing simple sequencing sample sheets."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional


@dataclass
class SampleEntry:
    """Single row in the sample sheet."""

    sample_id: str
    r1_path: Path
    r2_path: Optional[Path] = None

    def validate(self) -> None:
        """Ensure referenced FASTQ files exist."""

        if not self.r1_path.exists():
            raise FileNotFoundError(f"R1 FASTQ not found: {self.r1_path}")
        if self.r2_path is not None and not self.r2_path.exists():
            raise FileNotFoundError(f"R2 FASTQ not found: {self.r2_path}")


class SampleSheet:
    """Collection of :class:`SampleEntry` objects parsed from disk."""

    def __init__(self, entries: Iterable[SampleEntry]):
        self.entries: List[SampleEntry] = list(entries)

    @classmethod
    def from_file(cls, path: Path, validate_paths: bool = True) -> "SampleSheet":
        """Load a sample sheet from ``path``.

        Parameters
        ----------
        path:
            Location of the tab separated sample sheet. The format is expected to
            be ``sample_id`` followed by one or two FASTQ file paths. Comments
            starting with ``#`` and blank lines are ignored.
        validate_paths:
            Whether to ensure that the referenced FASTQ files exist on disk.
        """

        root = path.parent
        entries: List[SampleEntry] = []
        with path.open("rt", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 2:
                    raise ValueError(
                        f"Expected at least two tab separated columns on line {line_number}"
                    )
                sample_id = fields[0].strip()
                r1_path = (root / fields[1].strip()).resolve()
                r2_path = (root / fields[2].strip()).resolve() if len(fields) >= 3 else None
                entry = SampleEntry(sample_id=sample_id, r1_path=r1_path, r2_path=r2_path)
                if validate_paths:
                    entry.validate()
                entries.append(entry)
        if not entries:
            raise ValueError("Sample sheet does not contain any entries")
        return cls(entries)

    def __iter__(self):
        return iter(self.entries)

    def __len__(self) -> int:
        return len(self.entries)


__all__ = ["SampleEntry", "SampleSheet"]
