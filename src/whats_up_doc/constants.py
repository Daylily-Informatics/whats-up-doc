"""Reference data used by the heuristic inference engine."""

from __future__ import annotations

from collections import OrderedDict

# Illumina adapter sequences compiled from vendor documentation and common kits.
ADAPTER_SEQUENCES = OrderedDict(
    [
        ("TruSeq Universal Adapter", "AGATCGGAAGAG"),
        ("TruSeq LT Adapter", "AGATCGGAAGAGCACACGTCT"),
        ("Nextera Transposase Adapter", "CTGTCTCTTATACACATCT"),
        ("Nextera Read 2", "CTGTCTCTTATACACATCTGACGCTGCCGACGA"),
        ("Nextera Read 1", "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"),
        ("NEBNext Universal Adapter", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"),
        ("Ion Torrent Adapter A", "CCATCTCATCCCTGCGTGTCTCCGACTCAG"),
        ("Ion Torrent Adapter P1", "CCTCTCTATGGGCAGTCGGTGAT"),
        ("Swift Accel-NGS", "GTTCAGAGTTCTACAGTCCGACGATC"),
    ]
)

# Keywords describing common organisms. Values are ordered to allow deterministic output.
ORGANISM_KEYWORDS = OrderedDict(
    [
        ("Homo sapiens", ["human", "homo", "hg19", "hg38", "grch", "hs" ]),
        ("Mus musculus", ["mouse", "mm", "mus", "grcm", "mm10"]),
        ("Rattus norvegicus", ["rat", "rn" ]),
        ("Danio rerio", ["zebrafish", "dr", "danio"]),
        ("Drosophila melanogaster", ["dros", "dm", "fruitfly"]),
        ("Arabidopsis thaliana", ["arabidopsis", "at"]),
        ("Saccharomyces cerevisiae", ["yeast", "sc", "s288c"]),
        ("Escherichia coli", ["ecoli", "e.coli", "mg1655", "ec" ]),
        ("Mycobacterium tuberculosis", ["mtb", "h37rv", "mycobacter" ]),
        ("Plasmodium falciparum", ["plasmodium", "pfal", "pf3d7"]),
    ]
)

# Keywords describing tissues or sample origins.
TISSUE_KEYWORDS = OrderedDict(
    [
        ("Blood", ["blood", "plasma", "serum", "pbmc"]),
        ("Brain", ["brain", "cortex", "hippocampus", "neuron"]),
        ("Liver", ["liver", "hepato"]),
        ("Heart", ["heart", "cardiac"]),
        ("Lung", ["lung", "pulmo"]),
        ("Kidney", ["kidney", "renal"]),
        ("Tumor", ["tumor", "tumour", "cancer", "carcinoma"]),
        ("Plant", ["leaf", "root", "seed", "flower"]),
        ("Microbial culture", ["culture", "isolate", "strain"]),
    ]
)

LIBRARY_KEYWORDS = OrderedDict(
    [
        ("RNA-seq", ["rna", "mrna", "transcript", "transcriptome"]),
        ("Single cell RNA-seq", ["scrna", "10x", "cellranger"]),
        ("Total RNA-seq", ["totalrna", "ribo", "ribodeplete"]),
        ("ChIP-seq", ["chip", "h3k", "h3"]),
        ("ATAC-seq", ["atac", "tn5"]),
        ("Whole genome sequencing", ["wgs", "genome", "illumina dna prep"]),
        ("Whole exome sequencing", ["wes", "exome", "target"]),
        ("Metagenomic", ["metagenome", "metagenomic", "microbiome"]),
        ("Amplicon", ["amplicon", "16s", "its"]),
        ("Methylation", ["methyl", "bsseq", "bisulfite"]),
    ]
)

ENRICHMENT_KEYWORDS = OrderedDict(
    [
        ("Poly-A selection", ["polya", "poly-a", "mrna"]),
        ("Ribosomal depletion", ["ribodeplete", "rrna", "ribominus"]),
        ("Exome capture", ["exome", "sureselect", "twist", "xgen"]),
        ("PCR amplicon", ["amplicon", "pcr"]),
        ("Immune repertoire", ["vdj", "igh", "tcr"]),
    ]
)

# Mapping instrument prefixes to sequencing platforms.
INSTRUMENT_PREFIXES = OrderedDict(
    [
        (("NS", "NB"), "Illumina NextSeq"),
        (("MN",), "Illumina MiniSeq"),
        (("M0", "M1", "MI", "MS"), "Illumina MiSeq"),
        (("A", "D", "E"), "Illumina NovaSeq"),
        (("J",), "Illumina HiSeq X"),
        (("K", "L"), "Illumina NovaSeq X"),
        (("HWI", "HWUSI", "HS", "HI", "HC"), "Illumina HiSeq"),
        (("ST", "SL", "SN"), "Illumina HiSeq 3000/4000"),
        (("VH", "V3", "V4"), "Illumina NextSeq 2000"),
    ]
)

# GC content heuristics for coarse organism class prediction.
GC_CONTENT_HEURISTICS = [
    (0.65, "High GC content (possible bacteria or GC-rich organism)"),
    (0.50, "Moderate GC content (many microbial or plant genomes)"),
    (0.40, "Mammalian-like GC content"),
    (0.30, "AT-rich content (possible parasites or amplicons)"),
]

__all__ = [
    "ADAPTER_SEQUENCES",
    "ORGANISM_KEYWORDS",
    "TISSUE_KEYWORDS",
    "LIBRARY_KEYWORDS",
    "ENRICHMENT_KEYWORDS",
    "INSTRUMENT_PREFIXES",
    "GC_CONTENT_HEURISTICS",
]
