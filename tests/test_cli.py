from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from whats_up_doc import analyze_samplesheet

from .utils import example_paired_reads, write_fastq


def test_cli_generates_report(tmp_path):
    r1_records, r2_records = example_paired_reads()
    r1_path = write_fastq(tmp_path / "sample_cli_R1.fastq", r1_records)
    r2_path = write_fastq(tmp_path / "sample_cli_R2.fastq.gz", r2_records)

    sheet_path = tmp_path / "sheet.tsv"
    sheet_path.write_text(
        f"Sample_cli\t{Path(r1_path).name}\t{Path(r2_path).name}\n",
        encoding="utf-8",
    )

    cmd = [
        sys.executable,
        "-m",
        "whats_up_doc.cli",
        str(sheet_path),
        "--max-reads",
        "50",
        "--pretty",
    ]
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    payload = json.loads(completed.stdout)

    assert "samples" in payload
    assert payload["samples"][0]["sample_id"] == "Sample_cli"
    assert payload["samples"][0]["sequencing_run"]["platform"] == "Illumina NextSeq"


def test_library_entry_point_matches_cli(tmp_path):
    r1_records, r2_records = example_paired_reads()
    r1_path = write_fastq(tmp_path / "sample_lib_R1.fastq", r1_records)
    r2_path = write_fastq(tmp_path / "sample_lib_R2.fastq.gz", r2_records)

    sheet_path = tmp_path / "sheet.tsv"
    sheet_path.write_text(
        f"Sample_lib\t{Path(r1_path).name}\t{Path(r2_path).name}\n",
        encoding="utf-8",
    )

    report = analyze_samplesheet(sheet_path, max_reads=50)

    assert "samples" in report
    assert report["samples"][0]["sample_id"] == "Sample_lib"
    assert report["samples"][0]["sequencing_run"]["platform"] == "Illumina NextSeq"
