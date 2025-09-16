"""Heuristic inference helpers."""

from __future__ import annotations

import re
from collections import Counter
from typing import Dict, Iterable, List, Optional

from .constants import (
    ADAPTER_SEQUENCES,
    ENRICHMENT_KEYWORDS,
    GC_CONTENT_HEURISTICS,
    INSTRUMENT_PREFIXES,
    LIBRARY_KEYWORDS,
    ORGANISM_KEYWORDS,
    TISSUE_KEYWORDS,
)


def normalize_text(value: str) -> str:
    return re.sub(r"\s+", " ", value.strip().lower())


def find_keyword_matches(value: str, keyword_map) -> List[str]:
    value_norm = normalize_text(value)
    matches: List[str] = []
    for label, keywords in keyword_map.items():
        for keyword in keywords:
            if keyword and keyword.lower() in value_norm:
                matches.append(label)
                break
    return matches


def infer_organisms(sample_id: str, gc_fraction: Optional[float] = None) -> Dict[str, List[str]]:
    """Infer candidate organisms from text and GC content."""

    candidates = find_keyword_matches(sample_id, ORGANISM_KEYWORDS)
    annotations: List[str] = []
    if gc_fraction is not None:
        for threshold, description in GC_CONTENT_HEURISTICS:
            if gc_fraction >= threshold:
                annotations.append(description)
                break
        else:
            annotations.append("Very AT rich content")
    return {"candidates": candidates, "gc_annotation": annotations}


def infer_tissues(sample_id: str) -> List[str]:
    return find_keyword_matches(sample_id, TISSUE_KEYWORDS)


def infer_library_types(sample_id: str, adapter_hits: Iterable[str]) -> List[str]:
    annotations = set(find_keyword_matches(sample_id, LIBRARY_KEYWORDS))
    for adapter_name in adapter_hits:
        if "Nextera" in adapter_name:
            annotations.add("Tagmentation library")
        if "TruSeq" in adapter_name:
            annotations.add("TruSeq-based library")
        if "Ion Torrent" in adapter_name:
            annotations.add("Ion Torrent library")
    return sorted(annotations)


def infer_enrichments(sample_id: str) -> List[str]:
    return find_keyword_matches(sample_id, ENRICHMENT_KEYWORDS)


def instrument_to_platform(instrument: Optional[str]) -> Optional[str]:
    if not instrument:
        return None
    instrument_upper = instrument.upper()
    for prefixes, platform in INSTRUMENT_PREFIXES.items():
        for prefix in prefixes:
            if instrument_upper.startswith(prefix):
                return platform
    return None


def classify_index_type(index_sequences: Iterable[str]) -> str:
    unique_indexes = {index for index in index_sequences if index}
    if not unique_indexes:
        return "none"
    contains_dual = any("+" in index for index in unique_indexes)
    if contains_dual:
        return "dual"
    return "single"


def detect_umi_from_headers(headers: Iterable[str]) -> bool:
    pattern = re.compile(r"umi[_:-]?[acgt]+", re.IGNORECASE)
    for header in headers:
        if pattern.search(header):
            return True
    return False


def summarise_adapter_hits(hit_counter: Counter) -> List[Dict[str, float | int]]:
    total_hits = sum(hit_counter.values())
    summaries: List[Dict[str, float | int]] = []
    if total_hits == 0:
        return summaries
    for adapter_name, count in hit_counter.most_common():
        summaries.append(
            {
                "adapter": adapter_name,
                "count": count,
                "fraction_of_hits": count / total_hits,
            }
        )
    return summaries


def adapter_catalog() -> Dict[str, str]:
    return dict(ADAPTER_SEQUENCES)
