import datetime

import numpy as np
import pandas as pd

from utide._time_conversion import _date2num, _normalize_time


def test_formats():
    forms = [
        (np.array([693595.1]), "python"),
        (np.array([693961.1]), "matlab"),
        (np.array([2.1]), datetime.date(1899, 12, 29)),
        (np.array([3.1]), datetime.datetime(1899, 12, 28)),
        (np.array([2.6]), datetime.datetime(1899, 12, 28, 12)),
        (np.array([2.1]), "1899-12-29"),
    ]

    expected = _normalize_time(*forms[0])

    for form in forms[1:]:
        np.testing.assert_almost_equal(_normalize_time(*form), expected)


def test_datenum():
    assert _date2num("1970-01-01") == 0
    assert _date2num("1970-01-01", epoch="0001-01-01") == 719162
    assert _date2num("0001-01-01") == -719162
    assert _date2num("0001-01-01", epoch="0000-12-31") == 1
    assert _date2num(datetime.datetime(1, 1, 1)) == -719162
    assert _date2num(pd.to_datetime("1970-01-01")) == 0


def test_normalize_time_numpy():
    tin = np.arange("0001-01-01", "0001-01-10", dtype="datetime64[D]")  # till 09
    texp = np.arange(1, 10)
    assert np.array_equal(_normalize_time(tin, epoch=None), texp)

    tin = np.arange("1970-01-01", "1970-01-10", dtype="datetime64[D]")  # till 09
    texp = np.arange(719163, 719172)
    assert np.array_equal(_normalize_time(tin, epoch=None), texp)


def test_normalize_time_pandas():
    tin = pd.date_range("1970-01-01", "1970-01-10")  # till 10
    texp = np.arange(719163, 719173)
    assert np.array_equal(_normalize_time(tin, epoch=None), texp)
