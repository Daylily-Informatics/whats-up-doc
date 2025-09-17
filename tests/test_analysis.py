from __future__ import annotations

from whats_up_doc.analysis import analyze_sample
from whats_up_doc.samplesheet import SampleEntry

from .utils import example_paired_reads, write_fastq


def test_analyze_sample_reports_metadata(tmp_path):
    r1_records, r2_records = example_paired_reads()
    r1_path = write_fastq(tmp_path / "sample_human_rna_R1.fastq", r1_records)
    r2_path = write_fastq(tmp_path / "sample_human_rna_R2.fastq.gz", r2_records)

    entry = SampleEntry("Sample_Human_RNA", r1_path, r2_path)
    metadata = analyze_sample(entry, max_reads=100)

    assert metadata["paired_end"] is True
    assert metadata["sequencing_run"]["platform"] == "Illumina NextSeq"
    assert metadata["sequencing_run"]["index_type"] == "dual"
    assert "Homo sapiens" in metadata["organism"]["candidates"]
    assert "RNA-seq" in metadata["library"]["library_type_guesses"]
    assert metadata["read_metrics"]["r1"]["read_count"] == len(r1_records)
    assert metadata["library"]["adapter_hits"]
    assert metadata["library"]["adapter_hit_rate"] > 0  # Should have detected adapter sequence
