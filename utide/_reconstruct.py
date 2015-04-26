from __future__ import absolute_import, division

import numpy as np
from .harmonics import ut_E
from .utilities import Bunch

def reconstruct(t, coef, **opts):
    """
    Reconstruct a tidal signal.

    Parameters
    ----------
    t : array_like
        Time in days since `epoch`.
    coef : `Bunch`
        Data structure returned by `utide.solve`
    epoch : {string, int, `datetime.datetime`}, optional
        Not implemented yet. It will default to the epoch
        used for `coef`.

    Returns
    -------
    tide : `Bunch`
        Scalar time series is returned as `tide.h`; a vector
        series as `tide.u`, `tide.v`.

    """

    out = Bunch()
    u, v = _reconstr1(t, coef, **opts)
    if coef['aux']['opt']['twodim']:
        out.u, out.v = u, v
    else:
        out.h = u
    return out


def _reconstr1(tin, coef, **opts):

    print('reconstruct:')

    # parse inputs and options
    t, opt = _rcninit(tin, **opts)

    # determine constituents to include
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

    # exponentials

    ngflgs = [coef['aux']['opt']['nodsatlint'],
              coef['aux']['opt']['nodsatnone'],
              coef['aux']['opt']['gwchlint'],
              coef['aux']['opt']['gwchnone']]

    print('prep/calcs...')

    E = ut_E(t,
             coef['aux']['reftime'], coef['aux']['frq'][ind],
             coef['aux']['lind'][ind], coef['aux']['lat'], ngflgs,
             coef['aux']['opt']['prefilt'])

    # fit
    # fit = E*ap + np.conj(E)*am
    fit = np.dot(E, ap) + np.dot(np.conj(E), am)

    # mean (& trend)
    u = np.nan * np.ones(tin.shape)
    whr = ~np.isnan(tin)
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

        v = []

    print('Done.\n')

    return u, v

def _rcninit(tin, **opts):

    t = tin[:]

    t[np.isnan(t)] = []
    # t(isnan(t)) = []
    opt = {}

    opt['cnstit'] = False
    opt['minsnr'] = 2
    opt['minpe'] = 0

    for key, item in opts.items():
        try:
            opt[key] = item
        except KeyError:
            print('reconstruct: unrecognized input: {0}'.format(key))

    # args = list(args)
    # args = [string.lower() for string in args]

    # Need an example of the args

#    while ~isempty(args)
#        switch(lower(args{1}))
#            case 'cnstit'
#                opt.cnstit = args{2};
#                args(1:2) = [];
#            case 'minsnr'
#                opt.minsnr = args{2};
#                args(1:2) = [];
#            case 'minpe'
#                opt.minpe = args{2};
#                args(1:2) = [];
#            otherwise
#                error(['reconstruct: unrecognized input: ' args{1}]);
#        end
#    end

    return t, opt


