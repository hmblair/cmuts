from importlib.metadata import version
from .vis import visualize, visualize_atom
from .plotting import generate
from .internal import (
    title,
    DataGroups,
    Opts,
    NormScheme,
    ProbingData,
    normalize,
    stats,
)

__version__ = version("cmuts")

normalize = internal.normalize
stats = internal.stats
