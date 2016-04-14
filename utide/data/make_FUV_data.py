"""
Script to generate data for testing FUV.

It must be run in a directory containing ut_FUV.m and
ut_constants.mat. ut_FUV.m needs to be modified to load
ut_constants.mat at the beginning, in addition to all the
intermediate points where it is already being loaded.
"""

import numpy as np
from scipy.io.matlab import savemat

from oct2py import octave, Oct2PyError
octave.convert_to_float = False


t0 = octave.datenum(1950., 1., 2.)    # floats are required
print('t0 = ', t0)

t = np.linspace(t0, t0+300, 5)[:, None]
# Probably an octave oddity: the following *must* be a float
lat = 30.0
# Otherwise, in the expression "pi * lat / 180" or any variation
# of it, all operations seem to involve conversion of arguments
# to integers.  This is not an oct2py problem; when lat is saved
# in a matfile, and then loaded in an interactive octave session,
# the problem persists.


linds = [1 + np.arange(146, dtype=int)[:, None],
         [7, 8],
         [12, 13, 14],
         ]

for ilind, lind in enumerate(linds):
    shape = (7, len(t), len(lind))
    Fo = np.zeros(shape, dtype=float)
    Uo = Fo.copy()
    Vo = Fo.copy()

    flags = [[0, 0, 0, 0],
             [0, 1, 0, 0],
             [0, 0, 0, 1],
             [0, 1, 0, 1],
             [1, 0, 0, 0],
             [0, 0, 1, 0],
             [1, 0, 1, 0]]

    for i, flag in enumerate(flags):
        print(flag)
        try:
            F, U, V = octave.ut_FUV(t, t0, lind, lat, flag)
            Fo[i] = F
            Uo[i] = U
            Vo[i] = V
        except Oct2PyError:
            print('failed')

    save_args = dict(t=t, t0=t0, lat=lat, lind=lind, flags=flags,
                     Fo=Fo, Uo=Uo, Vo=Vo)

    np.savez('FUV%d.npz' % ilind, **save_args)

    savemat('FUV%d.mat' % ilind, save_args)
