"""Core routines for extracting metadata from FASTQ files."""

from __future__ import annotations

import math
import statistics
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from .inference import (
    classify_index_type,
    detect_umi_from_headers,
    infer_enrichments,
    infer_library_types,
    infer_organisms,
    infer_tissues,
    instrument_to_platform,
    summarise_adapter_hits,
)
from .constants import ADAPTER_SEQUENCES
from .io import read_fastq
from .samplesheet import SampleEntry


@dataclass
class ReadStats:
    """Container accumulating summary statistics for a FASTQ file."""

    read_count: int = 0
    total_bases: int = 0
    gc_bases: int = 0
    lengths: List[int] = field(default_factory=list)
    read_mean_qualities: List[float] = field(default_factory=list)
    min_length: Optional[int] = None
    max_length: Optional[int] = None
    min_read_mean_quality: Optional[float] = None
    max_read_mean_quality: Optional[float] = None
    sum_quality_scores: int = 0
    base_counts: Counter = field(default_factory=Counter)
    cycle_quality_sums: List[float] = field(default_factory=list)
    cycle_quality_counts: List[int] = field(default_factory=list)
    headers_sampled: List[str] = field(default_factory=list)
    index_sequences: Counter = field(default_factory=Counter)
    adapter_hits: Counter = field(default_factory=Counter)
    first_header_info: Optional[Dict[str, str]] = None

    def add_read(self, header: str, sequence: str, quality: str) -> None:
        sequence = sequence.strip()
        quality = quality.strip()
        if len(sequence) == 0:
            return
        if len(sequence) != len(quality):
            # Skip malformed records.
            return

        length = len(sequence)
        self.read_count += 1
        self.total_bases += length
        self.lengths.append(length)
        self.min_length = length if self.min_length is None else min(self.min_length, length)
        self.max_length = length if self.max_length is None else max(self.max_length, length)

        seq_upper = sequence.upper()
        self.base_counts.update(seq_upper)
        self.gc_bases += sum(1 for base in seq_upper if base in {"G", "C"})

        scores = [ord(char) - 33 for char in quality]
        mean_quality = sum(scores) / length if length else 0.0
        self.read_mean_qualities.append(mean_quality)
        self.sum_quality_scores += sum(scores)
        self.min_read_mean_quality = (
            mean_quality
            if self.min_read_mean_quality is None
            else min(self.min_read_mean_quality, mean_quality)
        )
        self.max_read_mean_quality = (
            mean_quality
            if self.max_read_mean_quality is None
            else max(self.max_read_mean_quality, mean_quality)
        )

        if len(self.cycle_quality_sums) < length:
            extend_by = length - len(self.cycle_quality_sums)
            self.cycle_quality_sums.extend([0.0] * extend_by)
            self.cycle_quality_counts.extend([0] * extend_by)
        for idx, score in enumerate(scores):
            self.cycle_quality_sums[idx] += score
            self.cycle_quality_counts[idx] += 1

        if len(self.headers_sampled) < 100:
            self.headers_sampled.append(header)

        header_info = parse_read_header(header)
        if self.first_header_info is None:
            self.first_header_info = header_info

        index_sequence = header_info.get("index_sequence")
        if index_sequence:
            self.index_sequences[index_sequence] += 1

        for adapter_name, adapter_seq in ADAPTER_SEQUENCES.items():
            if adapter_seq and adapter_seq in seq_upper:
                self.adapter_hits[adapter_name] += 1

    def gc_fraction(self) -> Optional[float]:
        if self.total_bases == 0:
            return None
        return self.gc_bases / self.total_bases

    def mean_quality(self) -> Optional[float]:
        if self.total_bases == 0:
            return None
        return self.sum_quality_scores / self.total_bases

    def summary(self) -> Dict[str, object]:
        if self.read_count == 0:
            return {
                "read_count": 0,
                "average_length": 0,
                "median_length": 0,
                "gc_fraction": None,
            }
        average_length = self.total_bases / self.read_count if self.read_count else 0
        gc_fraction = self.gc_fraction()
        n_bases = self.base_counts.get("N", 0)
        n_fraction = n_bases / self.total_bases if self.total_bases else 0.0
        per_cycle_quality = []
        for total, count in zip(self.cycle_quality_sums, self.cycle_quality_counts):
            per_cycle_quality.append(total / count if count else math.nan)
        return {
            "read_count": self.read_count,
            "average_length": average_length,
            "median_length": statistics.median(self.lengths) if self.lengths else 0,
            "length_range": [self.min_length, self.max_length],
            "gc_fraction": gc_fraction,
            "n_fraction": n_fraction,
            "base_composition": dict(self.base_counts),
            "quality": {
                "mean_phred": self.mean_quality(),
                "median_read_mean_phred": statistics.median(self.read_mean_qualities)
                if self.read_mean_qualities
                else None,
                "min_read_mean_phred": self.min_read_mean_quality,
                "max_read_mean_phred": self.max_read_mean_quality,
                "per_cycle_mean": per_cycle_quality,
            },
            "index_sequences": self.index_sequences.most_common(10),
            "adapter_hits": dict(self.adapter_hits),
            "header_examples": self.headers_sampled[:5],
            "first_header": self.first_header_info,
        }


def parse_read_header(header: str) -> Dict[str, str]:
    header = header.strip()
    if header.startswith("@"):
        header = header[1:]
    parts = header.split(" ", 1)
    first_part = parts[0]
    second_part = parts[1] if len(parts) > 1 else ""
    data: Dict[str, str] = {}
    tokens = first_part.split(":")
    if len(tokens) >= 4:
        data["instrument"] = tokens[0]
        data["run_number"] = tokens[1]
        data["flowcell_id"] = tokens[2]
        data["lane"] = tokens[3]
    if len(tokens) >= 5:
        data["tile"] = tokens[4]
    if len(tokens) >= 6:
        data["x_pos"] = tokens[5]
    if len(tokens) >= 7:
        data["y_pos"] = tokens[6]

    if second_part:
        second_tokens = second_part.split(":")
        if len(second_tokens) >= 1:
            data["read"] = second_tokens[0]
        if len(second_tokens) >= 2:
            data["is_filtered"] = second_tokens[1]
        if len(second_tokens) >= 3:
            data["control_number"] = second_tokens[2]
        if len(second_tokens) >= 4:
            data["index_sequence"] = second_tokens[3]
    return data


def collect_read_statistics(path: Path, max_reads: int) -> ReadStats:
    stats = ReadStats()
    for header, sequence, quality in read_fastq(path, max_reads=max_reads):
        stats.add_read(header, sequence, quality)
    return stats


def _combine_counters(counters: Iterable[Counter]) -> Counter:
    combined = Counter()
    for counter in counters:
        combined.update(counter)
    return combined


def _aggregate_headers(stats_list: Iterable[ReadStats]) -> List[str]:
    headers: List[str] = []
    for stats in stats_list:
        headers.extend(stats.headers_sampled)
    return headers


def analyze_sample(entry: SampleEntry, max_reads: int = 50000) -> Dict[str, object]:
    """Inspect FASTQ files for ``entry`` and return metadata."""

    r1_stats = collect_read_statistics(entry.r1_path, max_reads=max_reads)
    r2_stats = collect_read_statistics(entry.r2_path, max_reads=max_reads) if entry.r2_path else None

    stats_list = [r1_stats] + ([r2_stats] if r2_stats else [])
    total_sampled_reads = sum(stats.read_count for stats in stats_list)

    combined_index_counter = _combine_counters(stats.index_sequences for stats in stats_list)
    combined_adapter_counter = _combine_counters(stats.adapter_hits for stats in stats_list)

    gc_fraction_values = [stats.gc_fraction() for stats in stats_list if stats.gc_fraction() is not None]
    gc_fraction = sum(gc_fraction_values) / len(gc_fraction_values) if gc_fraction_values else None

    instrument = r1_stats.first_header_info.get("instrument") if r1_stats.first_header_info else None
    flowcell_id = r1_stats.first_header_info.get("flowcell_id") if r1_stats.first_header_info else None
    run_number = r1_stats.first_header_info.get("run_number") if r1_stats.first_header_info else None
    lane = r1_stats.first_header_info.get("lane") if r1_stats.first_header_info else None

    platform = instrument_to_platform(instrument)

    adapter_summary = summarise_adapter_hits(combined_adapter_counter)
    adapter_names = [item["adapter"] for item in adapter_summary]
    top_adapter_fraction = 0.0
    if combined_adapter_counter and total_sampled_reads:
        top_adapter_fraction = (
            combined_adapter_counter.most_common(1)[0][1] / total_sampled_reads
        )

    organism_info = infer_organisms(entry.sample_id, gc_fraction)
    tissues = infer_tissues(entry.sample_id)
    library_types = infer_library_types(entry.sample_id, adapter_names)
    enrichments = infer_enrichments(entry.sample_id)
    index_type = classify_index_type(combined_index_counter.keys())
    umi_present = detect_umi_from_headers(_aggregate_headers(stats_list))

    r1_summary = r1_stats.summary()

    read_metrics = {
        "r1": r1_summary,
        "r2": r2_stats.summary() if r2_stats else None,
    }

    sequencing_run = {
        "instrument": instrument,
        "platform": platform,
        "flowcell_id": flowcell_id,
        "run_number": run_number,
        "lane": lane,
        "index_type": index_type,
        "index_sequences": combined_index_counter.most_common(10),
        "umi_in_headers": umi_present,
    }

    library = {
        "adapter_hits": adapter_summary,
        "library_type_guesses": library_types,
        "enrichment_guesses": enrichments,
        "adapter_hit_rate": top_adapter_fraction,
    }

    warnings: List[str] = []
    primary_quality = r1_stats.mean_quality()
    if primary_quality is not None and primary_quality < 25:
        warnings.append("Mean R1 quality score is below Q25")
    if top_adapter_fraction > 0.2:
        warnings.append(
            "Adapter sequence present in more than 20% of sampled reads (consider trimming)"
        )
    if r1_summary["n_fraction"] > 0.05:
        warnings.append("More than 5% ambiguous bases in R1 reads")

    metadata = {
        "sample_id": entry.sample_id,
        "files": {
            "r1": str(entry.r1_path),
            "r2": str(entry.r2_path) if entry.r2_path else None,
        },
        "paired_end": entry.r2_path is not None,
        "sequencing_run": sequencing_run,
        "library": library,
        "organism": organism_info,
        "tissue_guesses": tissues,
        "read_metrics": read_metrics,
        "warnings": warnings,
    }
    return metadata


__all__ = ["analyze_sample", "collect_read_statistics", "parse_read_header", "ReadStats"]
