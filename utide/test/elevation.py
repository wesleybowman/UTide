from __future__ import division
import numpy as np
import cPickle as pickle
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv

#time = np.arange(734776, 734795+1, 1)
#time = np.arange(735604, 735964, 1/(24*10))
time = np.arange(735604, 735964+1e-8, 1/(24*10))

y = 15*np.cos(time)
x = y*y

lat = 44.86061
Rayleigh = np.array([1])
coef = ut_solv(time, x, np.array([]), lat,
               'auto', Rayleigh[0], 'NoTrend', 'Rmin', 'OLS',
               'NoDiagn', 'LinCI')

print coef
pickle.dump(coef, open("elevCoef.p", "wb"))
