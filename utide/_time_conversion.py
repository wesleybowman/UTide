"""
Utility for allowing flexible time input.
"""

from __future__ import (absolute_import, division, print_function)

import warnings
from datetime import date, datetime

try:
    from datetime import timezone
    have_tz = True
except ImportError:
    have_tz = False


def _normalize_time(t, epoch):
    if epoch == 'python':
        return t
    if epoch == 'matlab':
        return t - 366
    try:
        epoch = datetime.strptime(epoch, "%Y-%m-%d")
    except (TypeError, ValueError):
        pass
    if isinstance(epoch, date):
        if not hasattr(epoch, 'time'):
            return t + epoch.toordinal()
        # It must be a datetime, which is also an instance of date.
        if epoch.tzinfo is not None:
            if have_tz:
                epoch = epoch.astimezone(timezone.utc)
            else:
                warnings.warn("Timezone info in epoch is being ignored;"
                              " UTC is assumed.")
        ofs = (epoch.toordinal() + epoch.hour / 24 +
               epoch.minute / 1440 + epoch.second / 86400)
        return t + ofs
    raise ValueError("Cannot parse epoch as string or date or datetime")
