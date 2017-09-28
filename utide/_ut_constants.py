from __future__ import (absolute_import, division, print_function)

import os

from .utilities import convert_unicode_arrays, loadbunch


_base_dir = os.path.join(os.path.dirname(__file__), 'data')
_ut_constants_fname = os.path.join(_base_dir, 'ut_constants.npz')

# At least for now, use NaNs rather than masked arrays.
ut_constants = loadbunch(_ut_constants_fname, masked=False)
ut_constants = convert_unicode_arrays(ut_constants)

constit_names = list(ut_constants.const.name)

# Make a dictionary for index lookups.
constit_index_dict = dict([(name, i) for (i, name) in
                           enumerate(constit_names)])
