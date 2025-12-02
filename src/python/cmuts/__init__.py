from importlib.metadata import version as _get_version

# Core data types
from .internal import (
    DataGroups,
    Opts,
    NormScheme,
    ProbingData,
)

# Core functions
from .internal import compute_reactivity, stats, title

# Visualization submodule
from . import visualize

__version__ = _get_version("cmuts")

__all__ = [
    "__version__",
    "DataGroups",
    "Opts",
    "NormScheme",
    "ProbingData",
    "compute_reactivity",
    "stats",
    "title",
    "visualize",
]
