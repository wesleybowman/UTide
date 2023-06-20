import os

import numpy as np

from utide.utilities import loadbunch

datadir = os.path.join(os.path.dirname(__file__))

fname = "ut_constants"
data = loadbunch(os.path.join(datadir, f"{fname}.mat"), masked=False)
np.savez(f"{fname}.npz", **data)

fname = "FUV0"
data = loadbunch(os.path.join(datadir, f"{fname}.mat"), masked=False)
np.savez(f"{fname}.npz", **data)
