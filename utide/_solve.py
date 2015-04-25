"""
Central module for calculating the tidal amplitudes, phases, etc.
"""

from __future__ import absolute_import, division

import numpy as np

from .harmonics import ut_E
from .diagnostics import ut_diagn
from .ellipse_params import ut_cs2cep
from .constituent_selection import ut_cnstitsel
from .confidence import _confidence
from . import constit_index_dict

def solve(tin, uin, vin=None, lat=None, **opts):
    '''
    Need to put in docstring and figure a good way to put in all the optional
    parameters

    Keyword Arguments:
    conf_int=True
    cnstit='auto'
    notrend=0
    prefilt=[]
    nodsatlint=0
    nodsatnone=0
    gwchlint=0
    gwchnone=0
    infer=[]
    inferaprx=0
    rmin=1
    method='cauchy'
    tunrdn=1
    linci=0
    white=0
    nrlzn=200
    lsfrqosmp=1
    nodiagn=0
    diagnplots=0
    diagnminsnr=2
    ordercnstit=None
    runtimedisp='yyy'
    '''

    coef = _solv1(tin, uin, vin, lat, **opts)

    return coef


def _solv1(tin, uin, vin, lat, **opts):

    print('solve: ')

    # The following returns a possibly modified copy of tin (ndarray).
    # t, u, v are fully edited ndarrays (unless v is None)
    packed = _slvinit(tin, uin, vin, lat, **opts)
    tin, t, u, v, tref, lor, elor, opt = packed
    nt = len(t)

    # opt['cnstit'] = cnstit
    nNR, nR, nI, cnstit, coef = ut_cnstitsel(tref, opt['rmin']/(24*lor),
                                             opt['cnstit'], opt['infer'])

    # a function we don't need
    # coef.aux.rundescr = ut_rundescr(opt,nNR,nR,nI,t,tgd,uvgd,lat)

    coef['aux']['opt'] = opt
    coef['aux']['lat'] = lat

    print('matrix prep ... ')

    ngflgs = [opt['nodsatlint'], opt['nodsatnone'],
              opt['gwchlint'], opt['gwchnone']]

    # Make the model array, starting with the harmonics.
    E = ut_E(t, tref, cnstit['NR']['frq'], cnstit['NR']['lind'],
             lat, ngflgs, opt['prefilt'])

    # Positive and negative frequencies, and the mean.
    B = np.hstack((E, E.conj(), np.ones((nt, 1))))

    if not opt['notrend']:
        B = np.hstack((B, ((t-tref)/lor)[:, np.newaxis]))

    nm = B.shape[1]  #  2*(nNR + nR) + 1, plus 1 if trend is included

    print('Solution ...')

    if opt['twodim']:
        xraw = u + 1j*v
    else:
        xraw = u

    if opt['method'] == 'ols':
        # m = B\xraw;
        m = np.linalg.lstsq(B, xraw)[0]   # model coefficients
        # W = sparse(1:nt,1:nt,1);
        W = np.ones(nt)  # Uniform weighting; we could use a scalar 1, or None
    else:
        raise NotImplementedError("Only method 'ols' has been implemented")
#        [m,solnstats] = robustfit(B,ctranspose(xraw),...
#            opt.method,opt.tunconst,'off');
#        if isequal(lastwarn,'Iteration limit reached.')
#            # nan-fill, create coef.results, reorder coef fields,
#            % do runtime display
#            coef = ut_finish(coef,nNR,nR,nI,elor,cnstit);
#            % abort remainder of calcs
#            return;
#        W = sparse(1:nt,1:nt,solnstats.w);

    xmod = np.dot(B, m)   # model fit

    if not opt['twodim']:
        xmod = np.real(xmod)

    e = W*(xraw-xmod)  # Weighted residuals

    nc = nNR + nR

    ap = np.hstack((m[:nNR], m[2*nNR:2*nNR+nR]))
    i0 = 2*nNR + nR
    am = np.hstack((m[nNR:2*nNR], m[i0:i0+nR]))

    Xu = np.real(ap + am)
    Yu = -np.imag(ap - am)

    if not opt['twodim']:
        coef['A'], _, _, coef['g'] = ut_cs2cep(Xu, Yu)
        Xv = []
        Yv = []

    else:
        Xv = np.imag(ap+am)
        Yv = np.real(ap-am)
        packed = ut_cs2cep(Xu, Yu, Xv, Yv)
        coef['Lsmaj'], coef['Lsmin'], coef['theta'], coef['g'] = packed

    # mean and trend
    if opt['twodim']:
        if opt['notrend']:
            coef['umean'] = np.real(m[-1])
            coef['vmean'] = np.imag(m[-1])
        else:
            coef['umean'] = np.real(m[-2])
            coef['vmean'] = np.imag(m[-2])
            coef['uslope'] = np.real(m[-1])/lor
            coef['vslope'] = np.imag(m[-1])/lor
    else:
        if opt['notrend']:
            coef['mean'] = np.real(m[-1])
        else:
            coef['mean'] = np.real(m[-2])
            coef['slope'] = np.real(m[-1])/lor

    if opt['conf_int'] is True:
        coef = _confidence(coef, opt, t, e, tin, elor, xraw, xmod,
                           W, m, B, nc, Xu, Yu, Xv, Yv)

    # diagnostics
    if not opt['nodiagn']:
        coef, indPE = ut_diagn(coef, opt)

    # re-order constituents
    if opt['ordercnstit'] is not None:

        if opt['ordercnstit'] == 'frq':
            ind = coef['aux']['frq'].argsort()

        elif opt['ordercnstit'] == 'snr':
            if not opt['nodiagn']:
                ind = coef['diagn']['SNR'].argsort()[::-1]
            else:
                if opt['twodim']:
                    SNR = (coef['Lsmaj']**2 + coef['Lsmin']**2) / (
                        (coef['Lsmaj_ci']/1.96)**2 +
                        (coef['Lsmin_ci']/1.96)**2)

                else:
                    SNR = (coef['A']**2) / (coef['A_ci']/1.96)**2

                ind = SNR.argsort()[::-1]

        else:
            ilist = [constit_index_dict[name] for name in opt['ordercnstit']]
            ind = np.array(ilist, dtype=int)

    else:    # any other string: order by decreasing energy
        if not opt['nodiagn']:
            ind = indPE

        else:
            if opt['twodim']:
                PE = np.sum(coef['Lsmaj']**2 + coef['Lsmin']**2)
                PE = 100 * (coef['Lsmaj']**2 + coef['Lsmin']**2) / PE
            else:
                PE = 100 * coef['A']**2 / np.sum(coef['A']**2)

            ind = PE.argsort()[::-1]

    reorderlist = ['g', 'name']
    if opt['twodim']:
        reorderlist += ['Lsmaj', 'Lsmin', 'theta']
        if opt['conf_int']:
            reorderlist += ['Lsmaj_ci', 'Lsmin_ci', 'theta_ci', 'g_ci']
    else:
        reorderlist += ['A']
        if opt['conf_int']:
            reorderlist += ['A_ci']

    for key in reorderlist:
        coef[key] = coef[key][ind]

    coef['aux']['frq'] = coef['aux']['frq'][ind]
    coef['aux']['lind'] = coef['aux']['lind'][ind]

    print("Done.\n")

    return coef


def _slvinit(tin, uin, vin, lat, **opts):

    if lat is None:
        raise ValueError("Latitude must be supplied")

    # Supporting only 1-D arrays for now; we can add "group"
    # support later.
    if tin.shape != uin.shape or tin.ndim != 1 or uin.ndim != 1:
        raise ValueError("t and u must be 1-D arrays")

    if vin is not None and vin.shape != uin.shape:
        raise ValueError("v must have the same shape as u")

    opt = dict(twodim=(vin is not None))

    # Step 1: remove invalid times from tin, uin, vin
    tin = np.ma.masked_invalid(tin)
    uin = np.ma.masked_invalid(uin)
    if vin is not None:
        vin = np.ma.masked_invalid(vin)
    if np.ma.is_masked(tin):
        mask = np.ma.getmaskarray(tin)
        uin = uin.compress(mask)
        if vin is not None:
            vin = vin.compress(mask)

    tin = tin.compressed() # no longer masked

    # Step 2: generate t, u, v from edited tin, uin, vin
    v = None
    if np.ma.is_masked(uin) or np.ma.is_masked(vin):
        mask = np.ma.getmaskarray(uin)
        if vin is not None:
            mask = np.ma.mask_or(np.ma.getmaskarray(vin), mask)
        t = tin.compress(mask)
        u = uin.compress(mask).filled()
        if vin is not None:
            v = vin.compress(mask).filled()
    else:
        t = tin
        u = uin.filled()
        if vin is not None:
            v = vin.filled()

    # Now t, u, v, tin are clean ndarrays; uin and vin are masked,
    # but don't necessarily have masked values.

    # Are the times equally spaced?
    eps = np.finfo(np.float64).eps
    if np.var(np.unique(np.diff(tin))) < eps:
        opt['equi'] = True  # based on times; u/v can still have nans ("gappy")
        lor = np.ptp(tin)
        ntgood = len(tin)
        elor = lor*ntgood / (ntgood-1)
        tref = 0.5*(tin[0] + tin[-1])
    else:
        opt['equi'] = False
        lor = np.ptp(t)
        elor = lor*nt / (nt-1)
        tref = 0.5*(t[0]+t[-1])

    # Options.
    opt['conf_int'] = True
    opt['cnstit'] = 'auto'
    opt['notrend'] = 0
    opt['prefilt'] = []
    opt['nodsatlint'] = 0
    opt['nodsatnone'] = 0
    opt['gwchlint'] = 0
    opt['gwchnone'] = 0
    opt['infer'] = []
    opt['inferaprx'] = 0
    opt['rmin'] = 1
    # opt['method'] = 'cauchy'
    opt['method'] = 'ols'
    opt['tunrdn'] = 1
    opt['linci'] = 0
    opt['white'] = 0
    opt['nrlzn'] = 200
    opt['lsfrqosmp'] = 1
    opt['nodiagn'] = 0
    opt['diagnplots'] = 0
    opt['diagnminsnr'] = 2
    opt['ordercnstit'] = None
    opt['runtimedisp'] = 'yyy'

    # Update the default opt dictionary with the kwargs,
    # ensuring that every kwarg key matches a key in opt.
    for key, item in opts.items():
        try:
            opt[key] = item
        except KeyError:
            print('solve: unrecognized input: {0}'.format(key))

    allmethods = ['ols', 'andrews', 'bisquare', 'fair', 'huber',
                  'logistic', 'talwar', 'welsch']

    if opt['method'] != 'cauchy':
        ind = np.argwhere(opt['method'] in allmethods)[0][0]
        allconst = [np.nan, 1.339, 4.685, 1.400, 1.345, 1.205, 2.795, 2.985]
        opt['tunconst'] = allconst[ind]
    else:
        opt['tunconst'] = 2.385

    opt['tunconst'] = opt['tunconst'] / opt['tunrdn']

    return tin, t, u, v, tref, lor, elor, opt


