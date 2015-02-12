import os

base_dir = os.path.join(os.path.dirname(__file__), 'data')
ut_constants = os.path.join(base_dir, 'ut_constants.mat')


from ut_solv import ut_solv
from ut_reconstr import ut_reconstr


__version__ = '0.1b0.dev0'
