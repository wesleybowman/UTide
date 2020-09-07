"""
Tests for periodogram module.
"""

import numpy as np

import utide.periodogram as pgram


def random_ts(ndays, dt_hours, is_complex=True):
    """Returns t (time in days) and x (random series)."""
    np.random.seed(1)
    npts = int(ndays * 24 / dt_hours)
    if npts % 2:
        npts -= 1
    t = np.arange(npts, dtype=float) * dt_hours / 24
    if is_complex:
        x = np.random.randn(npts) + 1j * np.random.randn(npts)
    else:
        x = np.random.randn(npts)
    return t, x


def test_fft_ls_consistency():
    """
    With the exception of the frequency band including the
    Nyquist, band_psd should yield identical results with
    equi=True (fft method) and False (Lomb-Scargle method).

    """
    t, x = random_ts(20, 0.5)
    y_fft = pgram.band_psd(t, x, [1 / 12.42], equi=True)
    y_ls = pgram.band_psd(t, x, [1 / 12.42], equi=False)
    # skip the last frequency bin because y_fft includes nyquist
    sl = slice(0, -1)
    for key in ["Puu", "Pvv", "Puv"]:
        np.testing.assert_array_almost_equal(y_fft[key][sl], y_ls[key][sl])


def test_uv_consistency():
    t, x = random_ts(20, 0.5, is_complex=False)
    xx = x + 1j * x
    y_fft = pgram.band_psd(t, xx, [1 / 12.42], equi=True)
    y_ls = pgram.band_psd(t, xx, [1 / 12.42], equi=False)
    for y in y_fft, y_ls:
        np.testing.assert_array_almost_equal(y.Puu, y.Pvv)
        np.testing.assert_array_almost_equal(y.Puu, y.Puv)


""" TODO
def test_lmb_uneven():
    n = 200
    frequency = 10
    ofac = 1
    w = None
    t = np.linspace(0, 1, n, endpoint=False)
    x = np.array(np.exp(1j*2*np.pi*frequency*t))
    x_keep = np.ones(x.shape, dtype=bool)
    x_keep[2:5] = False
    tc = np.compress(x_keep, t)
    xc = np.compress(x_keep, x)
    pxx, f = lmbscga(xc.real, tc, w, ofac)

def test_lmbscga():
    # Input parameters.
    n = 200
    frequency = 10
    ofac = 1
    t = np.linspace(0, 1, n, endpoint=False)
    x = np.array(np.exp(1j*2*np.pi*frequency*t))
    ampl = 1
    w = None
    # Calculate Lomb-Scargle periodogram.
    pxx, f = lmbscga(x.real, t, w, ofac)
    # assert_approx_equal(np.max(pxx), ampl, significant=2)
    np.testing.assert_almost_equal(pxx[frequency],
                                    n / (4 * np.pi), decimal=6)
    pxx[frequency] = 0
    np.testing.assert_array_almost_equal(pxx, np.zeros_like(f), decimal=6)

def test_lmbscgc():
    # Input parameters
    n = 200
    frequency = 10
    ofac = 1
    w = None
    t = np.linspace(0, 1, n, endpoint=False)
    x = np.array(np.exp(1j*2*np.pi*frequency*t))
    y = np.array(np.exp(1j*2*np.pi*frequency*t))

    # Calculate Lomb-Scargle periodogram.
    pxy, f = lmbscgc(x.real, y.real, t, w, ofac)
    # np.testing.assert_almost_equal(pxy[frequency], n/ (4 * np.pi), decimal=6)
    pxy[frequency] = 0
    np.testing.assert_array_almost_equal(pxy, np.zeros_like(f), decimal=6)

    assert_array_equal(Pxx, Pxy)
    assert_array_equal(freqsxx, freqsxy)

evenly spaced, odd number of points
evenly spaced, missing values
"""
