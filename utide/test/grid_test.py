from __future__ import division
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import cPickle as pickle
import netCDF4 as nc
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv, ut_reconstr

def mjd2num(x):

    y = x + 678942

    return y


filename = '/home/wesley/ncfiles/smallcape_force_0001.nc'

data = nc.Dataset(filename, 'r')
x = data.variables['x'][:]
y = data.variables['y'][:]
lonc = data.variables['lonc'][:]
latc = data.variables['latc'][:]
lon = data.variables['lon'][:]
lat = data.variables['lat'][:]
ua = data.variables['ua']
va = data.variables['va']
time = data.variables['time'][:]
trinodes = data.variables['nv'][:]

time = mjd2num(time)

Rayleigh = np.array([1])


order = ['M2','S2','N2','K2','K1','O1','P1','Q1']
#order = ['M2  ','S2  ','N2  ','K2  ','K1  ','O1  ','P1  ','Q1  ']

#coef = ut_solv(time, time_series, time_series, lat, cnstit=order,
#               notrend=True, rmin=0.95, method='ols',
#               nodiagn=True, linci=True, conf_int=True, ordercnstit=order)
#
#speedcoef = ut_solv(time, time_series, time_series, lat, cnstit=order, gwchnone=True,
#               nodsatnone=True, notrend=True, rmin=0.95, method='ols',
#               nodiagn=True, linci=True, conf_int=True, ordercnstit=order)

speedcoef = ut_solv(time, ua[:, 0], va[:, 0], latc[0], cnstit='auto',
               notrend=True, rmin=0.95, method='ols',
               nodiagn=True, linci=True, conf_int=True)

#
pickle.dump(speedcoef, open("speedpythoncoef.p", "wb"))

u,v = ut_reconstr(time, speedcoef)

pickle.dump([u,v], open("speedpythonrecon.p", "wb"))

#coef = ut_solv(time, time_series, [], lat, cnstit=order, gwchnone=True,
#               nodsatnone=True, notrend=True, rmin=0.95, method='ols',
#               nodiagn=True, linci=True, conf_int=True, ordercnstit=order)

#coef = ut_solv(time, ua[:,0], [], lat[0], cnstit='auto', gwchnone=True,
#               nodsatnone=True, notrend=True, rmin=0.95, method='ols',
#               nodiagn=True, linci=True, conf_int=True)

coef = ut_solv(time, ua[:,0], [], lat[0], cnstit='auto',
               notrend=True, rmin=0.95, method='ols',
               nodiagn=True, linci=True, conf_int=True)

#coef = ut_solv(time, time_series, [], lat, cnstit='auto',
#               notrend=True, rmin=0.95, method='ols',
#               nodiagn=True, linci=True, conf_int=True)

pickle.dump(coef, open("pythoncoef.p", "wb"))

ts_recon, _ = ut_reconstr(time, coef)
pickle.dump(ts_recon, open("pythonrecon.p", "wb"))
