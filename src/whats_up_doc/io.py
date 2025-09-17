"""I/O helpers for working with FASTQ files."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Generator, Iterable, Iterator, Optional, Tuple


def _is_gzip_magic(magic: bytes) -> bool:
    return len(magic) >= 2 and magic[0] == 0x1F and magic[1] == 0x8B


def open_maybe_gzip(path: Path):
    """Open ``path`` transparently handling gzip-compressed files."""

    with path.open("rb") as raw:
        magic = raw.read(2)
    if _is_gzip_magic(magic):
        return gzip.open(path, mode="rt", encoding="utf-8", newline="\n")
    return path.open("rt", encoding="utf-8", newline="\n")


def read_fastq(path: Path, max_reads: Optional[int] = None) -> Iterator[Tuple[str, str, str]]:
    """Yield ``(header, sequence, quality)`` tuples from ``path``.

    Parameters
    ----------
    path:
        Location of the FASTQ file.
    max_reads:
        Optional cap on the number of reads to stream.
    """

    count = 0
    with open_maybe_gzip(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not quality:
                break
            if not header.startswith("@"):
                raise ValueError(f"Malformed FASTQ header: {header!r}")
            yield header.rstrip("\n"), sequence.rstrip("\n"), quality.rstrip("\n")
            count += 1
            if max_reads is not None and count >= max_reads:
                break


__all__ = ["open_maybe_gzip", "read_fastq"]
