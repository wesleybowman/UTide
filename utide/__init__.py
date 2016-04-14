from __future__ import (absolute_import, division, print_function)

import os

from utide.utilities import loadbunch

__version__ = '0.1.0'


_base_dir = os.path.join(os.path.dirname(__file__), 'data')
_ut_constants_fname = os.path.join(_base_dir, 'ut_constants.npz')

# At least for now, use NaNs rather than masked arrays.
ut_constants = loadbunch(_ut_constants_fname, masked=False)

# A list of strings is much easier to work with.
constit_names = [n.strip() for n in ut_constants.const.name]

# Make a dictionary for index lookups.
constit_index_dict = dict([(name, i) for (i, name) in
                           enumerate(constit_names)])

from ._solve import solve  # noqa
from ._reconstruct import reconstruct  # noqa

__all__ = ['solve',
           'reconstruct']
