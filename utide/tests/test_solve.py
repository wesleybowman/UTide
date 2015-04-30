"""
Tests for periodogram module.
"""

from __future__ import division

import numpy as np

from utide import ut_constants
from utide import solve
from utide import reconstruct

def test_roundtrip():
    # Minimal conversion from simple_utide_test
    ts = 735604
    duration = 35

    time = np.linspace(ts, ts+duration, 842)
    tref = (time[-1] + time[0]) / 2

    const = ut_constants.const

    amp = 1.0
    phase = 53
    lat = 45.5

    freq_cpd = 24 * const.freq

    jj = 48-1  # Python index for M2

    arg = 2 * np.pi * (time - tref) * freq_cpd[jj] - np.deg2rad(phase)
    time_series = amp * np.cos(arg)

    # Add a tiny bit of noise to prevent division by zero? Didn't do it.
    #time_series += 1e-10 * np.random.randn(len(time_series))

#    speed_coef = solve(time, time_series, time_series, lat=lat, cnstit='auto',
#                         notrend=True, rmin=0.95, method='ols',
#                         nodiagn=True, linci=True, conf_int=True)

#    elev_coef = solve(time, time_series, lat=lat, cnstit='auto',
#                        gwchnone=True, nodsatnone=True, notrend=True,
#                        rmin=0.95, method='ols', nodiagn=True, linci=True,
#                        conf_int=True)

    opts = dict(constit='auto',
                phase='raw',
                nodal=False,
                trend=False,
                method='ols',
                conf_int='linear',
                Rayleigh_min=0.95,
                )

    speed_coef = solve(time, time_series, time_series, lat=lat, **opts)
    elev_coef = solve(time, time_series, lat=lat, **opts)

    amp_err = amp - elev_coef['A'][0]
    phase_err = phase - elev_coef['g'][0]
    ts_recon = reconstruct(time, elev_coef).h

    vel = reconstruct(time, speed_coef)

    err = np.sqrt(np.mean((time_series-ts_recon)**2))

    print(amp_err, phase_err, err)
    print(elev_coef['aux']['reftime'], tref)
    print(elev_coef['aux']['opt'])

    np.testing.assert_almost_equal(amp_err, 0)
    np.testing.assert_almost_equal(phase_err, 0)
    np.testing.assert_almost_equal(err, 0)

