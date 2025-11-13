from importlib.metadata import version
from .internal import (
    title,
    DataGroups,
    Opts,
    NormScheme,
    ProbingData,
    normalize,
    stats,
)
from . import visualize

__version__ = version("cmuts")

normalize = internal.normalize
stats = internal.stats
