from __future__ import absolute_import

import os

from utide.utilities import loadmatbunch

_base_dir = os.path.join(os.path.dirname(__file__), 'data')
_ut_constants_fname = os.path.join(_base_dir, 'ut_constants.mat')
# At least for now, use NaNs rather than masked arrays.
ut_constants = loadmatbunch(_ut_constants_fname, masked=False)

from ._solve import solve
from ._reconstruct import reconstruct


__version__ = '0.1b0.dev0'
