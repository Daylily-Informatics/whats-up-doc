"""Top level package for the :mod:`whats_up_doc` toolkit."""

from .analysis import analyze_sample, collect_read_statistics
from .cli import analyze_samplesheet
from .samplesheet import SampleEntry, SampleSheet

__all__ = [
    "__version__",
    "analyze_sample",
    "analyze_samplesheet",
    "collect_read_statistics",
    "SampleEntry",
    "SampleSheet",
]

__version__ = "0.1.0"
