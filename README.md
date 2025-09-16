# whats-up-doc

`seqmeta` is a small toolkit for extracting informative metadata from FASTQ files
referenced by a simple tab-separated sample sheet. It inspects the sequencing
headers and read content to produce a JSON report that summarises attributes
useful for downstream analysis and quality control such as the likely organism,
sequencing platform, indexing scheme, library preparation hints and read quality
metrics.

## Features

* Parses sample sheets with ``sample_id`` and R1/R2 FASTQ file paths.
* Supports plain or gzip-compressed FASTQ files (detected automatically by
  reading the file magic bytes).
* Samples the first N reads of each file (configurable) and reports:
  * Read length distribution, GC content, N content and mean/median quality
    scores.
  * Flowcell, lane and instrument information parsed from Illumina-style
    headers with heuristic mapping to the sequencing platform.
  * Index sequences and whether single or dual indexing is present.
  * Adapter contamination estimates for a curated catalogue of common library
    kits.
  * Keyword-based guesses for organism, tissue/source, library type and
    enrichment strategy using both the sample identifier and observed adapters.
  * Simple quality warnings (e.g., low average quality or excessive adapter
    signal).

The heuristics are intentionally conservative and designed to provide guidance
rather than definitive answers. They can be extended easily by editing
``seqmeta/constants.py``.

## Installation

The project is packaged as a standard Python module. Install it into a virtual
environment or use it directly via ``python -m``:

```bash
pip install .
```

## Usage

Create a tab-separated sample sheet with the following format:

```text
sample_id\t/path/to/sample_R1.fastq.gz\t/path/to/sample_R2.fastq.gz
```

The second column is mandatory; provide the third column when a paired-end file
is available. Paths may be relative to the location of the sample sheet.

Run the CLI on the sample sheet to produce a JSON report:

```bash
seqmeta path/to/samples.tsv --max-reads 75000 --pretty --output metadata.json
```

* ``--max-reads`` controls how many reads are sampled from each FASTQ file (the
  default is 50,000).
* ``--pretty`` toggles human-readable JSON indentation.
* ``--output`` writes the report to a file instead of printing to stdout.

## Development

This repository includes a pytest suite that exercises the sample sheet parser,
core analyser and CLI. After making changes run:

```bash
pip install -e .[test]
pytest
```

(The tests dynamically construct their own FASTQ fixtures, so there is no large data bundle tracked in version control.)

## License

MIT
