from __future__ import (absolute_import, division, print_function)

import datetime

import numpy as np
from utide._time_conversion import _normalize_time


def test_formats():
    forms = [(np.array([693595.1]), 'python'),
             (np.array([693961.1]), 'matlab'),
             (np.array([2.1]), datetime.date(1899, 12, 29)),
             (np.array([3.1]), datetime.datetime(1899, 12, 28)),
             (np.array([2.6]), datetime.datetime(1899, 12, 28, 12)),
             (np.array([2.1]), '1899-12-29')]

    expected = _normalize_time(*forms[0])

    for form in forms[1:]:
        np.testing.assert_almost_equal(_normalize_time(*form), expected)
