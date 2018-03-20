"""
Full example test.

"""

from matplotlib.dates import date2num

import numpy as np

from pandas import date_range

import pytest

from utide import reconstruct, solve
from utide.utilities import Bunch


@pytest.fixture
def fake_tide(t, M2amp, M2phase):
    """
    Generate a minimally realistic-looking fake semidiurnal tide.

    t is time in hours
    phases are in radians

    Modified from:
    http://currents.soest.hawaii.edu/ocn760_4/_static/plotting.html

    """
    return M2amp * np.sin(2 * np.pi * t / 12.42 - M2phase)


@pytest.fixture
def make_data():
    N = 500
    np.random.seed(1234)
    t = date_range(start='2016-03-29', periods=N, freq='H')
    # Signal + some noise.
    u = fake_tide(np.arange(N), M2amp=2, M2phase=0) + np.random.randn(N)
    v = fake_tide(np.arange(N), M2amp=1, M2phase=np.pi) + np.random.randn(N)
    time = date2num(t.to_pydatetime())
    return time, u, v


def test_solve():
    time, u, v = make_data()
    coef = solve(time, u, v,
                 lat=-42.5,
                 nodal=False,
                 trend=False,
                 method='ols',
                 conf_int='linear',
                 Rayleigh_min=0.95,)
    assert isinstance(coef, Bunch)

    tide = reconstruct(time, coef)
    assert isinstance(tide, Bunch)
