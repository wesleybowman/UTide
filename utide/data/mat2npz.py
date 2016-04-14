from __future__ import (absolute_import, division, print_function)

import os
import numpy as np
from utide.utilities import loadbunch


datadir = os.path.join(os.path.dirname(__file__))

fname = 'ut_constants'
data = loadbunch(
    os.path.join(datadir, '{}.mat'.format(fname)), masked=False)
np.savez('{}.npz'.format(fname), **data)

fname = 'FUV0'
data = loadbunch(
    os.path.join(datadir, '{}.mat'.format(fname)), masked=False)
np.savez('{}.npz'.format(fname), **data)
