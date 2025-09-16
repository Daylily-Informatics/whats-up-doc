from __future__ import annotations

from pathlib import Path

from seqmeta.samplesheet import SampleSheet

from .utils import example_paired_reads, write_fastq


def test_samplesheet_parsing(tmp_path):
    r1_records, r2_records = example_paired_reads()
    r1_path = write_fastq(tmp_path / "sample_R1.fastq", r1_records)
    r2_path = write_fastq(tmp_path / "sample_R2.fastq.gz", r2_records)

    sheet_path = tmp_path / "sheet.tsv"
    sheet_path.write_text(
        f"SampleA\t{Path(r1_path).name}\t{Path(r2_path).name}\n",
        encoding="utf-8",
    )

    sheet = SampleSheet.from_file(sheet_path)
    assert len(sheet) == 1
    entry = sheet.entries[0]
    assert entry.sample_id == "SampleA"
    assert entry.r1_path == r1_path.resolve()
    assert entry.r2_path == r2_path.resolve()
