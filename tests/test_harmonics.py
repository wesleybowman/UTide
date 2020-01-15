"""
Test for FUV in harmonics.py.

The test data are generated using octave with a slight
modification of ut_FUV extracted from ut_solv.m.  The
data-generating script is make_FUV_data.py.
"""

import os

from numpy.testing import assert_array_almost_equal

from utide._ut_constants import _base_dir
from utide.harmonics import FUV
from utide.utilities import loadbunch

fname = os.path.join(_base_dir, 'FUV0.npz')


def test_FUV():
    x = loadbunch(fname, masked=False)
    # Switch epoch from Matlab to Python
    x.t -= 366
    x.t0 -= 366

    for i, flag in enumerate(x.flags):
        F, U, V = FUV(x.t, x.t0, x.lind-1, x.lat, flag)
        print('i: {} ngflags: {}'.format(i, flag))

        # We use broadcasting instead of replication, so
        # we need to test only against the first row of
        # the octave output in such cases.
        if F.shape[0] == 1:
            sub = (i, slice(0, 1))
        else:
            sub = (i,)
        assert_array_almost_equal(F, x.Fo[sub])
        assert_array_almost_equal(U, x.Uo[sub])

        if V.shape[0] == 1:
            sub = (i, slice(0, 1))
        else:
            sub = (i,)
        assert_array_almost_equal(V, x.Vo[sub])
