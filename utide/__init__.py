import pkg_resources

from ._solve import solve
from ._reconstruct import reconstruct
from ._ut_constants import (
    ut_constants,
    constit_index_dict,
    hours_per_cycle,
    cycles_per_hour,
)

try:
    __version__ = pkg_resources.get_distribution("utide").version
except Exception:
    __version__ = "unknown"

__all__ = [
    "solve",
    "reconstruct",
    "ut_constants",
    "constit_index_dict",
    "hours_per_cycle",
    "cycles_per_hour",
]
