"""Command line interface for the :mod:`whats_up_doc` toolkit."""

from __future__ import annotations

import argparse
import json
import sys
from os import PathLike
from pathlib import Path
from typing import Any, Dict, List, Union

from .analysis import analyze_sample
from .samplesheet import SampleEntry, SampleSheet


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect FASTQ files referenced by a simple tab-separated sample sheet "
            "and extract sequencing metadata."
        )
    )
    parser.add_argument("samplesheet", type=Path, help="Path to the sample sheet (TSV)")
    parser.add_argument(
        "--max-reads",
        type=int,
        default=50000,
        help=(
            "Maximum number of reads to sample from each FASTQ file."
            " Larger values increase accuracy at the cost of runtime."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional location to write the JSON report (defaults to stdout)",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON with indentation for easier inspection",
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Skip verification that FASTQ files referenced in the sample sheet exist",
    )
    return parser


def analyze_samplesheet(
    samplesheet: Union[str, PathLike[str], SampleSheet],
    *,
    max_reads: int = 50000,
    validate_paths: bool = True,
) -> Dict[str, Any]:
    """Generate a metadata report for the provided sample sheet.

    Parameters
    ----------
    samplesheet:
        Either a filesystem path (``str`` or ``os.PathLike``) pointing to a tab
        separated sample sheet on disk or an in-memory
        :class:`~whats_up_doc.samplesheet.SampleSheet` instance.
    max_reads:
        Maximum number of reads to sample from each FASTQ file during analysis.
    validate_paths:
        When ``True`` (the default) FASTQ paths referenced by the sample sheet are
        checked for existence. When a :class:`SampleSheet` instance is supplied the
        contained entries are validated directly.
    """

    if isinstance(samplesheet, (str, PathLike)):
        sheet = SampleSheet.from_file(Path(samplesheet), validate_paths=validate_paths)
    else:
        sheet = samplesheet
        if validate_paths:
            for entry in sheet:
                if isinstance(entry, SampleEntry):
                    entry.validate()

    results: List[Dict[str, Any]] = []
    for entry in sheet:
        metadata = analyze_sample(entry, max_reads=max_reads)
        results.append(metadata)

    return {"samples": results}


def main(argv: List[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    report = analyze_samplesheet(
        args.samplesheet,
        max_reads=args.max_reads,
        validate_paths=not args.no_validate,
    )
    indent = 2 if args.pretty or args.output else None
    output_text = json.dumps(report, indent=indent)

    if args.output:
        args.output.write_text(output_text + "\n", encoding="utf-8")
    else:
        sys.stdout.write(output_text)
        if not output_text.endswith("\n"):
            sys.stdout.write("\n")


if __name__ == "__main__":  # pragma: no cover
    main()
