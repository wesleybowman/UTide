"""
Smoke testing--just see if the system runs.

"""

# These tests are quick and crude.
# TODO: extend the tests by cycling through various combinations
#       of configuration and data input.

import numpy as np
import pytest

from utide import reconstruct, solve
from utide._ut_constants import ut_constants
from utide.utilities import Bunch


ts = 735604
duration = 35

time = np.linspace(ts, ts + duration, 842)
tref = (time[-1] + time[0]) / 2

const = ut_constants.const

amp = 1.0
phase = 53
lat = 45.5

freq_cpd = 24 * const.freq

jj = 48 - 1  # Python index for M2

arg = 2 * np.pi * (time - tref) * freq_cpd[jj] - np.deg2rad(phase)
tide = amp * np.cos(arg)

np.random.seed(123454321)
noise = 1e-5 * np.random.randn(len(time))

time_series = tide + noise

# We omit the 'MC' case for now because with this test data, it
# fails with 'numpy.linalg.LinAlgError: SVD did not converge'.
@pytest.mark.parametrize("conf_int", ["linear", "none"])
def test_roundtrip(conf_int):
    """Minimal conversion from simple_utide_test."""

    opts = {
        "constit": "auto",
        "phase": "raw",
        "nodal": False,
        "trend": False,
        "method": "ols",
        "conf_int": conf_int,
        "Rayleigh_min": 0.95,
    }

    speed_coef = solve(time, time_series, time_series, lat=lat, **opts)
    elev_coef = solve(time, time_series, lat=lat, **opts)

    amp_err = amp - elev_coef["A"][0]
    phase_err = phase - elev_coef["g"][0]
    ts_recon = reconstruct(time, elev_coef).h

    # pure smoke testing of reconstruct
    vel = reconstruct(time, speed_coef)
    vel = reconstruct(time, speed_coef, constit=("M2", "S2"))
    htmp = reconstruct(time, elev_coef, constit=("M2", "S2"))
    vel = reconstruct(time, speed_coef, min_SNR=3)
    htmp = reconstruct(time, elev_coef, min_SNR=3)
    vel = reconstruct(time, speed_coef, min_PE=10)
    htmp = reconstruct(time, elev_coef, min_PE=10)
    vel = reconstruct(time, speed_coef, min_SNR=0)
    htmp = reconstruct(time, elev_coef, min_SNR=0)
    assert isinstance(vel, Bunch)
    assert isinstance(htmp, Bunch)

    # Now the round-trip check, just for the elevation.
    err = np.sqrt(np.mean((tide - ts_recon) ** 2))

    print(amp_err, phase_err, err)
    print(elev_coef["aux"]["reftime"], tref)
    print(elev_coef["aux"]["opt"])

    np.testing.assert_almost_equal(amp_err, 0, decimal=4)
    np.testing.assert_almost_equal(phase_err, 0, decimal=4)
    np.testing.assert_almost_equal(err, 0, decimal=4)


def test_masked_input():
    """Masked values in time and/or time series."""

    opts = {
        "constit": "auto",
        "phase": "raw",
        "nodal": False,
        "trend": False,
        "method": "ols",
        "conf_int": "linear",
        "Rayleigh_min": 0.95,
    }

    t = np.ma.array(time)
    t[[10, 15, 20, 21]] = np.ma.masked

    series = np.ma.array(time_series)
    series[[11, 17, 22, 25]] = np.ma.masked

    speed_coef = solve(t, series, series, lat=lat, **opts)
    elev_coef = solve(t, series, lat=lat, **opts)

    amp_err = amp - elev_coef["A"][0]
    phase_err = phase - elev_coef["g"][0]
    ts_recon = reconstruct(time, elev_coef).h
    assert isinstance(ts_recon, np.ndarray)

    # pure smoke testing of reconstruct
    vel = reconstruct(time, speed_coef)
    assert isinstance(vel, Bunch)

    elev = reconstruct(time, elev_coef)
    assert isinstance(elev, Bunch)

    np.testing.assert_almost_equal(amp_err, 0, decimal=4)
    np.testing.assert_almost_equal(phase_err, 0, decimal=4)


def test_robust():
    """
    Quick check that method='robust' works; no real checking
    of results, other than by using "py.test -s" and noting that
    the results are reasonable, and the weights for the outliers
    are very small.
    Minimal conversion from simple_utide_test

    """

    # Add noise
    np.random.seed(13579)
    noisy = tide + 0.01 * np.random.randn(len(time))

    # Add wild points
    noisy[:5] = 10
    noisy[-5:] = -10

    opts = {
        "constit": "auto",
        "phase": "raw",
        "nodal": False,
        "trend": False,
        "method": "robust",
        "conf_int": "linear",
        "Rayleigh_min": 0.95,
    }

    speed_coef = solve(time, noisy, noisy, lat=lat, **opts)
    elev_coef = solve(time, noisy, lat=lat, **opts)

    print(speed_coef.weights, elev_coef.weights)
    print(speed_coef.rf, elev_coef.rf)

    ts_recon = reconstruct(time, elev_coef).h
    err = np.std(tide - ts_recon)
    np.testing.assert_almost_equal(err, 0, decimal=2)


def test_MC():
    # Add noise
    np.random.seed(1)
    noisy = tide + 0.01 * np.random.randn(len(time))

    opts = {
        "constit": "auto",
        "phase": "raw",
        "nodal": False,
        "trend": False,
        "method": "ols",
        "conf_int": "MC",
        "white": False,
        "Rayleigh_min": 0.95,
    }

    speed_coef = solve(time, noisy, noisy, lat=lat, **opts)
    elev_coef = solve(time, noisy, lat=lat, **opts)

    for name, AA, AA_ci, gg, gg_ci in zip(
        elev_coef.name, elev_coef.A, elev_coef.A_ci, elev_coef.g, elev_coef.g_ci
    ):
        print("%5s %10.4g %10.4g  %10.4g %10.4g" % (name, AA, AA_ci, gg, gg_ci))

    for (name, Lsmaj, Lsmaj_ci, Lsmin, Lsmin_ci, theta, theta_ci, gg, gg_ci) in zip(
        speed_coef.name,
        speed_coef.Lsmaj,
        speed_coef.Lsmaj_ci,
        speed_coef.Lsmin,
        speed_coef.Lsmin_ci,
        speed_coef.theta,
        speed_coef.theta_ci,
        speed_coef.g,
        speed_coef.g_ci,
    ):
        print(
            "%5s %10.4g %10.4g  %10.4g %10.4g  %10.4g %10.4g  %10.4g %10.4g"
            % (name, Lsmaj, Lsmaj_ci, Lsmin, Lsmin_ci, theta, theta_ci, gg, gg_ci)
        )


def test_ordercnstit():
    # Add noise
    np.random.seed(1)
    noisy = tide + 0.01 * np.random.randn(len(time))

    opts = {
        "constit": "auto",
        "phase": "raw",
        "nodal": False,
        "trend": False,
        "method": "ols",
        "conf_int": "MC",
        "white": False,
        "Rayleigh_min": 0.95,
        "ordercnstit": "frq",
    }

    speed_coef = solve(time, noisy, noisy, lat=lat, **opts)
    elev_coef = solve(time, noisy, lat=lat, **opts)

    for name, AA, AA_ci, gg, gg_ci in zip(
        elev_coef.name, elev_coef.A, elev_coef.A_ci, elev_coef.g, elev_coef.g_ci
    ):
        print("%5s %10.4g %10.4g  %10.4g %10.4g" % (name, AA, AA_ci, gg, gg_ci))

    for (name, Lsmaj, Lsmaj_ci, Lsmin, Lsmin_ci, theta, theta_ci, gg, gg_ci) in zip(
        speed_coef.name,
        speed_coef.Lsmaj,
        speed_coef.Lsmaj_ci,
        speed_coef.Lsmin,
        speed_coef.Lsmin_ci,
        speed_coef.theta,
        speed_coef.theta_ci,
        speed_coef.g,
        speed_coef.g_ci,
    ):
        print(
            "%5s %10.4g %10.4g  %10.4g %10.4g  %10.4g %10.4g  %10.4g %10.4g"
            % (name, Lsmaj, Lsmaj_ci, Lsmin, Lsmin_ci, theta, theta_ci, gg, gg_ci)
        )
