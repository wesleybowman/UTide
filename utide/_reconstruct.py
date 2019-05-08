from __future__ import (absolute_import, division, print_function)

import numpy as np
from .harmonics import ut_E
from .utilities import Bunch
from ._time_conversion import _normalize_time


def reconstruct(t, coef,
                epoch='python',
                verbose=True,
                constit=None,
                min_SNR=2,
                min_PE=0):
    """
    Reconstruct a tidal signal.

    Parameters
    ----------
    t : array_like
        Time in days since ``epoch``.
    coef : `Bunch`
        Data structure returned by `utide.solve`.
    epoch : {string, `datetime.date`, `datetime.datetime`}, optional
        Valid strings are 'python' (default); 'matlab' if `t` is
        an array of Matlab datenums; or an arbitrary date in the
        form 'YYYY-MM-DD'.  The default corresponds to the Python
        standard library `datetime` proleptic Gregorian calendar,
        starting with 1 on January 1 of year 1.
    verbose : {True, False}, optional
        True to enable output message (default). False turns off all
        messages.
    constit : {None, array_like}, optional
        List of strings with standard letter abbreviations of
        tidal constituents, to be used in reconstruction if present
        in coef; alternative to the SNR and PE criteria.
    min_SNR : float, optional, default 2
        Include only the constituents with signal-to-noise SNR >= min_SNR,
        where SNR is based on the constituent confidence intervals in
        ``coef``.
    min_PE : float, optional, default 0
        Include only the constituents with percent energy PE >= min_PE,
        where PE is based on the amplitudes in ``coef``.

    Returns
    -------
    tide : `Bunch`
        Scalar time series is returned as `tide.h`; a vector
        series as `tide.u`, `tide.v`.  Each is an ndarray with
        ``np.nan`` as the missing value.
        Most input kwargs are included: 'epoch', 'constit',
        'min_SNR', and 'min_PE'.
        The input time array is included as 't_in', and 't_mpl';
        the former is the original input time argument, and the
        latter is the time as a matplotlib datenum.  If 'epoch'
        is 'python', these will be identical, and the names will
        point to the same array.

    """

    out = Bunch(t_in=t, epoch=epoch, constit=constit, min_SNR=min_SNR,
                min_PE=min_PE)
    t = np.atleast_1d(t)
    if t.ndim != 1:
        raise ValueError("t must be a 1-D array")
    t = _normalize_time(t, epoch)
    if epoch == 'python':
        out.t_mpl = out.t_in
    else:
        out.t_mpl = t
    t = np.ma.masked_invalid(t)
    goodmask = ~np.ma.getmaskarray(t)
    t = t.compressed()

    u, v = _reconstruct(t, goodmask, coef,
                        verbose=verbose,
                        constit=constit,
                        min_SNR=min_SNR,
                        min_PE=min_PE)

    if v is not None:
        out.u, out.v = u, v
    else:
        out.h = u
    return out


def _reconstruct(t, goodmask, coef, verbose, constit, min_SNR, min_PE):

    twodim = coef['aux']['opt']['twodim']

    # Determine constituents to include.
    if constit is not None:
        ind = [i for i, c in enumerate(coef['name']) if c in constit]
    elif min_SNR == 0 and min_PE == 0:
        ind = slice(None)
    else:
        if twodim:
            E = coef['Lsmaj']**2 + coef['Lsmin']**2
            N = (coef['Lsmaj_ci']/1.96)**2 + (coef['Lsmin_ci']/1.96)**2
        else:
            E = coef['A']**2
            N = (coef['A_ci']/1.96)**2
        SNR = E / N
        PE = 100 * E / E.sum()
        with np.errstate(invalid='ignore'):
            ind = np.logical_and(SNR >= min_SNR, PE >= min_PE)

    # Complex coefficients.
    rpd = np.pi/180
    if twodim:
        ap = 0.5 * ((coef['Lsmaj'][ind] + coef['Lsmin'][ind]) *
                    np.exp(1j*(coef['theta'][ind] - coef['g'][ind]) * rpd))
        am = 0.5 * ((coef['Lsmaj'][ind] - coef['Lsmin'][ind]) *
                    np.exp(1j*(coef['theta'][ind] + coef['g'][ind]) * rpd))
    else:
        ap = 0.5 * coef['A'][ind] * np.exp(-1j*coef['g'][ind] * rpd)
        am = np.conj(ap)

    ngflgs = [coef['aux']['opt']['nodsatlint'],
              coef['aux']['opt']['nodsatnone'],
              coef['aux']['opt']['gwchlint'],
              coef['aux']['opt']['gwchnone']]

    if verbose:
        print('prep/calcs ... ', end='')

    E = ut_E(t,
             coef['aux']['reftime'], coef['aux']['frq'][ind],
             coef['aux']['lind'][ind], coef['aux']['lat'], ngflgs,
             coef['aux']['opt']['prefilt'])

    fit = np.dot(E, ap) + np.dot(np.conj(E), am)

    # Mean (& trend).
    u = np.empty(goodmask.shape, dtype=float)
    u.fill(np.nan)
    trend = not coef['aux']['opt']['notrend']

    if twodim:
        v = u.copy()
        u[goodmask] = np.real(fit) + coef['umean']
        v[goodmask] = np.imag(fit) + coef['vmean']
        if trend:
            u[goodmask] += coef['uslope'] * (t - coef['aux']['reftime'])
            v[goodmask] += coef['vslope'] * (t - coef['aux']['reftime'])

    else:
        u[goodmask] = np.real(fit) + coef['mean']
        if trend:
            u[goodmask] += coef['slope'] * (t - coef['aux']['reftime'])
        v = None

    if verbose:
        print('done.')

    return u, v
