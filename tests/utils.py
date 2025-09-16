from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable, Sequence, Tuple


def write_fastq(path: Path, records: Iterable[Tuple[str, str]]) -> Path:
    """Write FASTQ records to ``path`` with constant high quality."""

    opener = gzip.open if path.suffix.endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8") as handle:
        for header, sequence in records:
            quality = "I" * len(sequence)
            handle.write(f"{header}\n{sequence}\n+\n{quality}\n")
    return path


def example_paired_reads() -> Tuple[Sequence[Tuple[str, str]], Sequence[Tuple[str, str]]]:
    """Return synthetic paired-end reads with Illumina-style headers."""

    r1_records = [
        (
            "@NS500123:45:HGTT2BGX2:1:11101:10000:100000 1:N:0:ACGTACGT+TGCATGCA",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC",
        ),
        (
            "@NS500123:45:HGTT2BGX2:1:11101:10000:100001 1:N:0:ACGTACGT+TGCATGCA",
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGATCGGAAGAGCACACGTCT",
        ),
        (
            "@NS500123:45:HGTT2BGX2:1:11101:10000:100002 1:N:0:ACGTACGT+TGCATGCA",
            "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC",
        ),
    ]
    r2_records = [
        (
            "@NS500123:45:HGTT2BGX2:1:11101:10000:100000 2:N:0:ACGTACGT+TGCATGCA",
            "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
        ),
        (
            "@NS500123:45:HGTT2BGX2:1:11101:10000:100001 2:N:0:ACGTACGT+TGCATGCA",
            "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT",
        ),
        (
            "@NS500123:45:HGTT2BGX2:1:11101:10000:100002 2:N:0:ACGTACGT+TGCATGCA",
            "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT",
        ),
    ]
    return r1_records, r2_records

