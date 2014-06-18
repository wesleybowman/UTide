from __future__ import division
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import cPickle as pickle
import sys
sys.path.append('/home/wesley/github/UTide/')
from utide import ut_solv, ut_reconstr


ts = 735604
duration = 35
dt = 1/24.0

time = np.linspace(ts, ts+duration, 841)
time_origin = 6.939615e5

mat_contents = sio.loadmat('ut_constants.mat', struct_as_record=False, squeeze_me=True)

shallow = mat_contents['shallow']
const = mat_contents['const']

amp = 1.0
phase = 53
lat = 45.5
period = 1 / const.freq * 3600

jj=48-1

time_series = amp * np.cos((((time-time_origin) * (2*np.pi/period[jj]) *
                             (24*3600)) - 2 * np.pi * phase / 360))

order = ['M2','S2','N2','K2','K1','O1','P1','Q1']
#order = ['M2  ','S2  ','N2  ','K2  ','K1  ','O1  ','P1  ','Q1  ']

#coef = ut_solv(time, time_series, time_series, lat, cnstit=order,
#               notrend=True, rmin=0.95, method='ols',
#               nodiagn=True, linci=True, conf_int=True, ordercnstit=order)
#
#coef = ut_solv(time, time_series, time_series, lat, cnstit=order, gwchnone=True,
#               nodsatnone=True, notrend=True, rmin=0.95, method='ols',
#               nodiagn=True, linci=True, conf_int=True, ordercnstit=order)

#speedcoef = ut_solv(time, time_series, time_series, lat, cnstit='auto',
#               notrend=True, rmin=0.95, method='ols',
#               nodiagn=True, linci=True, conf_int=True)
#
#
#pickle.dump(speedcoef, open("speedpythoncoef.p", "wb"))

#coef = ut_solv(time, time_series, [], lat, cnstit=order, gwchnone=True,
#               nodsatnone=True, notrend=True, rmin=0.95, method='ols',
#               nodiagn=True, linci=True, conf_int=True, ordercnstit=order)

coef = ut_solv(time, time_series, [], lat, cnstit='auto', gwchnone=True,
               nodsatnone=True, notrend=True, rmin=0.95, method='ols',
               nodiagn=True, linci=True, conf_int=True)

coef = ut_solv(time, time_series, [], lat, cnstit='auto',
               notrend=True, rmin=0.95, method='ols',
               nodiagn=True, linci=True, conf_int=True)


pickle.dump(coef, open("pythoncoef.p", "wb"))

amp_err = amp - coef['A'][0]

phase_err = phase - coef['g'][0]

ts_recon, _ = ut_reconstr(time, coef)

err = np.sqrt(np.mean((time_series-ts_recon[0])**2))

ts_fvcom=coef['A'][0]*np.cos(2*np.pi*((time-np.mean(time))/(period[jj]/(24*3600))-coef['g'][0]/360))

PLOT = False
if PLOT:
    plt.plot(time, ts_recon)
    plt.plot(time, ts_fvcom)
    plt.show()
