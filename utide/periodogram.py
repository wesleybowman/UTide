"""
Band-averaged spectral estimation for use in tidal error analysis.

There is only one function in this module that is imported by
other modules: `band_psd` (``ut_pdgm`` in the Matlab version), which
uses the residuals from the tidal fit to calculate power spectral
density averaged in bands near the constituent clusters.

Because the functions here are intended for internal use only,
they are stripped down to their essentials, with little or
no argument checking.  Some argument flexibility is included
to facilitate testing.
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np

from scipy import signal

from utide.utilities import Bunch


# __all__ = ['freq_bands', 'band_psd']

# Frequency bands in cycles per hour for
# estimating background noise levels.
freq_bands = np.array([[.00010, .00417],   # fortnightly etc.
                       [.03192, .04859],   # diurnal
                       [.07218, .08884],   # semidiurnal
                       [.11243, .12910],   # mk3, 2mk3
                       [.15269, .16936],   # about 4 cpd
                       [.19295, .20961],
                       [.23320, .25100],   # s6
                       [.26000, .29000],
                       [.30000, .50000]])  # m8 is 0.322; highest is 0.487


def fbndavg(P, freq, cfreq=None, fbands=None):
    """
    Band-average periodogram, excluding constituent frequencies

    Parameters
    ----------
    P : ndarray (nfreq,)
        periodogram in units^2 per CPH
    freq: ndarray (nfreq,)
        frequencies in CPH
    cfreq: None or array_like
        constituent frequencies (CPH) to be excluded
    fbands: array_like (nbands, 2)
        limits of frequency bands; defaults to 9 bands in `freq_bands`

    Returns
    -------
    ndarray (nbands,)
        mean of `P` within `fbands`

    Notes
    -----
    From ut_fbndavg, which was based on residual_spectrum.m of
    t_tide, Pawlowicz et al. (2002).
    """

    if fbands is None:
        fbands = freq_bands

    nfreq = len(freq)

    if cfreq is not None:
        # Exclude periodogram values that are nearest neighbors of
        # input constituent frequencies.
        P = np.ma.array(P, copy=True)
        indices = np.arange(nfreq)
        i_constit = np.round(np.interp(cfreq, freq, indices)).astype(int)
        i_constit = np.clip(i_constit, 0, nfreq)
        P[i_constit] = np.ma.masked

    nbands = fbands.shape[0]
    avP = np.ma.zeros((nbands,), dtype=P.dtype)

    for k, (f0, f1) in enumerate(fbands):
        i0, i1 = np.searchsorted(freq, [f0, f1])
        i1 = min(nfreq, i1+1)
        inband = slice(i0, i1)
        avP[k] = P[inband].mean()

    # Better to leave it as a masked array?
    avP = avP.filled(np.nan)

    return avP


def _lomb_freqs(t, fbands=None, ofac=1, max_per_band=500):
    """
    Calculate uniformly spaced frequencies in each band, with
    no more than `max_per_band` in any band.
    """

    if fbands is None:
        fbands = freq_bands

    # Estimated record length as n delta-t, where delta-t
    # is the *average* time per sample.
    n = len(t)
    delta_t = (t[-1] - t[0]) / (n-1)
    reclen = n * delta_t

    nf = n * ofac  # number of "Fourier frequencies" based on oversampling
    # Simplify by ignoring the 0 and Nyquist frequencies.
    # Divide by reclen to convert cycles/record to cycles/time unit.
    freq = np.arange(1, nf//2) / reclen
    nfreq = len(freq)

    freqs = []
    for k, (f0, f1) in enumerate(fbands):
        i0, i1 = np.searchsorted(freq, [f0, f1])
        i1 = min(nfreq, i1+1)
        inband = slice(i0, i1)
        if i1 - i0 > max_per_band:
            band = np.linspace(freq[i0], freq[i1-1], max_per_band)
        else:
            band = freq[inband]
        freqs.append(band)

    return np.hstack(freqs)


def _psd_lomb(t, x, window=None, freq=None, ofac=1):
    """
    Periodogram estimate for irregular sampling, Lomb-Scargle method

    Parameters
    ----------
    t : ndarray, (n,)
        time, monotonic but possibly irregularly spaced
    x : ndarray, (n,)
        signal, real or complex
    window : None or ndarray, (n,)
        if not None, a uniformly-sampled data window
    freq : ndarray (nfreq,)
        evaluation frequencies in cycles per unit time
    ofac : integer
        oversampling factor; defaults to 1

    Returns
    -------
    P : Bunch

        - P.F: Frequencies (units: cycles per time unit) of the Pxx estimates.
        - P.Pxx: One-sided auto-spectral density estimate for real(`x`).

        if `x` is complex:

        - P.Pyy: as above, for imag(`x`)
        - P.Pxy: complex cross-spectrum between real(`x`) and imag(`x`)

    Notes
    -----
    If `freq` is None, `P.F` is calculated to coincide with the
    Fourier frequencies for a series of n uniformly distributed
    times from min(`t`) to max(`t`). The mean and Nyquist are omitted
    because they are irrelevant in this context.

    PSD units are [`x`-units^2 per cycle per unit time]

    """
    out = Bunch()

    # copy inputs
    x = np.array(x)
    t = np.array(t, dtype=float)

    # remove mean
    x -= x.mean()

    n = len(x)

    if window is None:
        w = np.ones(t.shape, dtype=float)
    else:
        # interpolate window from uniform grid to nonuniform t
        t_uniform = np.linspace(np.min(t), np.max(t), n)
        w = np.interp(t, t_uniform, window)

        x *= w

    # Estimated record length as n delta-t, where delta-t
    # is the *average* time per sample.
    delta_t = (t[-1] - t[0]) / (n-1)
    reclen = n * delta_t

    if freq is None:
        ofac = int(round(ofac))
        nf = n * ofac  # number of "Fourier frequencies" based on oversampling
        # Simplify by ignoring the 0 and Nyquist frequencies.
        # Divide by reclen to convert cycles/record to cycles/time unit.
        freq = np.arange(1, nf//2) / reclen

    out.F = freq

    xr = np.real(x)

    # signal.lombscargle returns "(A**2) * N/4 for a harmonic signal
    # with amplitude A for sufficiently large N."
    # It takes *angular* frequencies as 3rd argument.
    freq_radian = freq * 2 * np.pi
    psdnorm = 2 * delta_t * n / (w**2).sum()
    out.Pxx = psdnorm * signal.lombscargle(t, xr, freq_radian)

    if x.dtype.kind == 'f':
        return out

    out.Pyy = psdnorm * signal.lombscargle(t, x.imag, freq_radian)

    # If we need to limit memory usage and don't want to use
    # Cython, we can segment the frequencies and loop over the
    # segments.  The speed penalty will be minimal.
    out.Pxy = psdnorm * _ls_cross(t, x, freq_radian)

    return out


def _ls_cross(t, x, fr):
    """
    Cross spectral counterpart of signal.lombscargle.
    """
    # The python code can be replaced with cython later.

    # Initial argument for phase shift calculation: double frequency
    arg = 2 * fr * t[:, np.newaxis]

    # One phase shift per frequency:
    tau = 0.5 * np.arctan2(np.sin(arg).sum(axis=0), np.cos(arg).sum(axis=0))

    # Phase-shifted argument of sine, cosine:
    arg *= 0.5   # from double frequency back to original frequency
    arg -= tau

    xr, xi = x.real, x.imag

    # Calculate unnormalized, 2-sided, mean-removed cross-periodogram

    tmpx = np.empty(fr.shape, dtype=complex)
    tmpy = np.empty(fr.shape, dtype=complex)

    c = np.cos(arg)
    a = 1 / np.sqrt(sum(c**2))
    tmpx.real = a * (xr[:, np.newaxis] * c).sum(axis=0)
    tmpy.real = a * (xi[:, np.newaxis] * c).sum(axis=0)

    del c

    s = np.sin(arg)
    b = 1 / np.sqrt(sum(s**2))
    tmpx.imag = b * (xr[:, np.newaxis] * s).sum(axis=0)
    tmpy.imag = -b * (xi[:, np.newaxis] * s).sum(axis=0)

    f0 = np.exp(1j * tau)
    tmpx *= f0
    tmpy *= np.conj(f0)
    pxy = 0.5 * (tmpx * tmpy)

    return pxy


def _psd(e, window, fs):
    """
    Autospectrum if `e` is real, or cross spectrum if complex.

    `e` must be uniformly spaced, with sampling frequency `fs` in
    cycles per time unit.  Number of points is assumed to be even.

    `window` is a set of data windowing weights, the length of `e`

    Returns PSD in cycles per time unit.
    """

    # iny = nt//2 +1
    # Detrending: just remove the mean.
    e = e - e.mean()
    e *= window
    if e.dtype.kind == 'c':
        cs = np.conj(np.fft.rfft(e.real)) * np.fft.rfft(e.imag)
    else:
        cs = np.abs(np.fft.rfft(e.real))**2

    # cs = cs[:iny]
    cs[1:-1] *= 2
    psdnorm = (1/fs) * (1/(window**2).sum())  # dt / sum of win squared.
    return cs * psdnorm


def band_psd(t, e, cfrq, equi=True, frqosamp=1):
    """
    Band-averaged power spectral density in vicinity of tidal constituents

    Parameters
    ----------
    t : ndarray
        time [days] relative to an arbitrary origin
    e : ndarray
        residual from tidal fit, complex/real
    cfrq : ndarray
        frequencies of NR & R constituents [CPH]
    equi : boolean
        True (default) if sample times are uniformly spaced
    frqosamp : integer
        Lomb-Scargle frequency oversampling factor, if equi = False

    Returns
    -------
    P : Bunch

        - P.fbnd : edges of freq bands [CPH]
        - P.Puu : 1-sided auto-sp dens of u err (=real(e)) [units^2/CPH]
        - P.Pvv : (complex case only)= as P.Puu but for v err (=imag(e))
        - P.Puv (complex case only)= cross-sp dens between u and v

    Note
    ----
    This is the only function in this module that is called from
    other modules.

    """

    P = Bunch(fbnd=freq_bands)

    nt = len(e)
    if nt % 2:
        e = e[:-1]
        t = t[:-1]
        nt -= 1

    hn = signal.windows.hann(nt, sym=False)

    # on real component
    if equi:  # If even sampling, FFT.
        # sample interval in hours; t is in days
        dt = 24 * (t[1] - t[0])

        fs = 1/dt  # sampling frequency: cycles (samples) per hour
        Puu1s = _psd(np.real(e), hn, fs)
        allfrq = np.arange(nt//2 + 1) / (nt * dt)

    else:  # If uneven, Lomb-scargle.
        # time in hours, for output in CPH and x^2 per CPH
        lfreq = _lomb_freqs(t * 24, fbands=freq_bands, ofac=frqosamp)
        ls_spec = _psd_lomb(t * 24, e, window=hn, freq=lfreq)
        Puu1s = ls_spec.Pxx
        allfrq = ls_spec.F

    P.Puu = fbndavg(Puu1s, allfrq, cfrq)

    # If e is complex, handle imaginary part.
    if e.dtype.kind == 'c':

        if equi:  # If even sampling, Welch.
            Pvv1s = _psd(e.imag, hn, fs)
            Puv1s = _psd(e, hn, fs)       # complex cross-periodogram

        else:  # If uneven, lomscargle.
            Pvv1s = ls_spec.Pyy
            Puv1s = ls_spec.Pxy

        P.Pvv = fbndavg(Pvv1s, allfrq, cfrq)

        # Take co-spectrum only; ignore the quadrature (imag) part.
        P.Puv = fbndavg(Puv1s, allfrq, cfrq).real

        # I don't think we want to throw away the sign of the
        # co-spectrum, which determines the quadrant in which the
        # semi-major variance ellipse lies.
        # If it does matter, we need to check the sign convention.
        # P.Puv = np.abs(P.Puv)

    return P
