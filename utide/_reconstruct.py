from __future__ import (absolute_import, division, print_function)

import numpy as np
from .harmonics import ut_E
from .utilities import Bunch
from ._time_conversion import _normalize_time


def reconstruct(t, coef, epoch='python', verbose=True, **opts):
    """
    Reconstruct a tidal signal.

    Parameters
    ----------
    t : array_like
        Time in days since `epoch`.
    coef : `Bunch`
        Data structure returned by `utide.solve`
    epoch : {string, `datetime.date`, `datetime.datetime`}, optional
        Valid strings are 'python' (default); 'matlab' if `t` is
        an array of Matlab datenums; or an arbitrary date in the
        form 'YYYY-MM-DD'.  The default corresponds to the Python
        standard library `datetime` proleptic Gregorian calendar,
        starting with 1 on January 1 of year 1.
    verbose : {True, False}, optional
        True to enable output message (default). False turns off all
        messages.

    Returns
    -------
    tide : `Bunch`
        Scalar time series is returned as `tide.h`; a vector
        series as `tide.u`, `tide.v`.

    """

    out = Bunch()
    u, v = _reconstr1(t, coef, epoch=epoch, verbose=verbose, **opts)
    if coef['aux']['opt']['twodim']:
        out.u, out.v = u, v
    else:
        out.h = u
    return out


def _reconstr1(tin, coef, **opts):

    # Parse inputs and options.
    t, goodmask, opt = _rcninit(tin, **opts)

    if opt['RunTimeDisp']:
        print('reconstruct:', end='')

    # Determine constituents to include.
    # if ~isempty(opt.cnstit)
    # if not np.empty(opt['cnstit']):
    if opt['cnstit']:

        # [~,ind] = ismember(cellstr(opt.cnstit),coef.name);
        # opt['cnstit'] in coef['name']
        ind = np.where(opt['cnstit'] == coef['name'])

#        if ~isequal(length(ind),length(cellstr(opt.cnstit)))
#            error(['reconstruct: one or more of input constituents Cnstit '...
#                'not found in coef.name']);
    else:

        ind = np.arange(len(coef['aux']['frq']))
        if coef['aux']['opt']['twodim']:
            SNR = ((coef['Lsmaj']**2 + coef['Lsmin']**2) /
                   ((coef['Lsmaj_ci']/1.96)**2 + (coef['Lsmin_ci']/1.96)**2))

            PE = sum(coef['Lsmaj']**2 + coef['Lsmin']**2)
            PE = 100*(coef['Lsmaj']**2 + coef['Lsmin']**2)/PE
        else:
            SNR = (coef['A']**2)/((coef['A_ci']/1.96)**2)
            PE = 100*coef['A']**2/sum(coef['A']**2)

        # ind = ind[SNR[ind]>=opt['minsnr'] & PE[ind]>=opt['minpe']]
        ind = np.where(np.logical_and(SNR[ind] >= opt['minsnr'],
                                      PE[ind] >= opt['minpe']))[0]

    # Complex coefficients.
    rpd = np.pi/180
    if coef['aux']['opt']['twodim']:
        ap = 0.5 * ((coef['Lsmaj'][ind] + coef['Lsmin'][ind]) *
                    np.exp(1j*(coef['theta'][ind] - coef['g'][ind]) * rpd))
        am = 0.5 * ((coef['Lsmaj'][ind] - coef['Lsmin'][ind]) *
                    np.exp(1j*(coef['theta'][ind] + coef['g'][ind]) * rpd))
    else:
        ap = 0.5 * coef['A'][ind] * np.exp(-1j*coef['g'][ind] * rpd)
        am = np.conj(ap)

    # Exponentials.

    ngflgs = [coef['aux']['opt']['nodsatlint'],
              coef['aux']['opt']['nodsatnone'],
              coef['aux']['opt']['gwchlint'],
              coef['aux']['opt']['gwchnone']]

    if opt['RunTimeDisp']:
        print('prep/calcs ... ', end='')

    E = ut_E(t,
             coef['aux']['reftime'], coef['aux']['frq'][ind],
             coef['aux']['lind'][ind], coef['aux']['lat'], ngflgs,
             coef['aux']['opt']['prefilt'])

    # Fit.
    # fit = E*ap + np.conj(E)*am
    fit = np.dot(E, ap) + np.dot(np.conj(E), am)

    # Mean (& trend).
    u = np.nan * np.ones(tin.shape)
    whr = goodmask
    if coef['aux']['opt']['twodim']:
        v = np.nan * np.ones(tin.shape)
        if coef['aux']['opt']['notrend']:
            u[whr] = np.real(fit) + coef['umean']
            v[whr] = np.imag(fit) + coef['vmean']
        else:
            u[whr] = np.real(fit) + coef['umean']
            u[whr] += coef['uslope'] * (t-coef['aux']['reftime'])
            v[whr] = np.imag(fit) + coef['vmean']
            v[whr] += coef['vslope'] * (t-coef['aux']['reftime'])

    else:
        if coef['aux']['opt']['notrend']:
            u[whr] = np.real(fit) + coef['mean']
        else:
            u[whr] = np.real(fit) + coef['mean']
            u[whr] += coef['slope'] * (t-coef['aux']['reftime'])

        v = None

    if opt['RunTimeDisp']:
        print('done.')

    return u, v


def _rcninit(tin, **opts):

    t = tin[:]

    # Supporting only 1-D arrays for now; we can add "group"
    # support later.
    if tin.ndim != 1:
        raise ValueError("t must be a 1-D array")

    # Step 0: apply epoch to time.
    t = _normalize_time(tin, opts['epoch'])

    # Step 1: remove invalid times from tin
    t = np.ma.masked_invalid(t)
    goodmask = ~np.ma.getmaskarray(t)
    t = t.compressed()

    opt = {}

    opt['cnstit'] = False
    opt['minsnr'] = 2
    opt['minpe'] = 0

    for key, item in opts.items():
        # Be backward compatible with the MATLAB package syntax.
        if key == 'verbose':
            opt['RunTimeDisp'] = item

        try:
            opt[key] = item
        except KeyError:
            print('reconstruct: unrecognized input: {0}'.format(key))

    return t, goodmask, opt
