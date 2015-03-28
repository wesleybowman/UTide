"""
Central module for calculating the tidal amplitudes, phases, etc.
"""

from __future__ import absolute_import, division

import numpy as np
import scipy       # This will go away; see FIXME below

from .harmonics import ut_E
from .ut_diagn import ut_diagn
from .ut_cs2cep import ut_cs2cep
from .ut_cnstitsel import ut_cnstitsel
from .ut_confidence import ut_confidence

def solve(tin, uin, vin, lat, **opts):
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
    ordercnstit=[]
    runtimedisp='yyy'
    '''

    coef = _solv1(tin, uin, vin, lat, **opts)

    return coef


def _solv1(tin, uin, vin, lat, **opts):

    print('solve: ')
    packed = _slvinit(tin, uin, vin, **opts)
    nt, t, u, v, tref, lor, elor, opt, tgd, uvgd = packed

    # opt['cnstit'] = cnstit
    [nNR, nR, nI, cnstit, coef] = ut_cnstitsel(tref, opt['rmin']/(24*lor),
                                               opt['cnstit'], opt['infer'])

    # a function we don't need
    # coef.aux.rundescr = ut_rundescr(opt,nNR,nR,nI,t,tgd,uvgd,lat)

    coef['aux']['opt'] = opt
    coef['aux']['lat'] = lat

    print('matrix prep ... ')

    ngflgs = [opt['nodsatlint'], opt['nodsatnone'],
              opt['gwchlint'], opt['gwchnone']]

    E = ut_E(t, tref, cnstit['NR']['frq'], cnstit['NR']['lind'],
             lat, ngflgs, opt['prefilt'])

    B = np.hstack((E, E.conj()))

    # more infer stuff

    if opt['notrend']:
        B = np.hstack((B, np.ones((nt, 1))))
        nm = 2 * (nNR + nR) + 1
    else:
        B = np.hstack((B, np.ones((nt, 1)), ((t-tref)/lor)[:, np.newaxis]))
        nm = 2*(nNR + nR) + 2

    print('Solution ...')

    xraw = u

    if opt['twodim']:
        # xraw = complex(u,v);
        xraw = u+v*1j

    if opt['method'] == 'ols':
        # m = B\xraw;
        m = np.linalg.lstsq(B, xraw)[0]
        # W = sparse(1:nt,1:nt,1);
        W = scipy.sparse.identity(nt)    ## FIXME: we shouldn't need this
#    else:
#        lastwarn('');
#        [m,solnstats] = robustfit(B,ctranspose(xraw),...
#            opt.method,opt.tunconst,'off');
#        if isequal(lastwarn,'Iteration limit reached.')
#            # nan-fill, create coef.results, reorder coef fields,
#            % do runtime display
#            coef = ut_finish(coef,nNR,nR,nI,elor,cnstit);
#            % abort remainder of calcs
#            return;
#        W = sparse(1:nt,1:nt,solnstats.w);

    xmod = np.dot(B, m)

    if not opt['twodim']:
        xmod = np.real(xmod)

    e = W*(xraw-xmod)

    nc = nNR+nR
    ap = m[np.hstack((np.arange(nNR), 2*nNR+np.arange(nR)))]
    am = m[np.hstack((nNR+np.arange(nNR), 2*nNR+nR+np.arange(nR)))]

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
            coef['umean'] = np.real(m[-1-1])
            coef['vmean'] = np.imag(m[-1-1])
            coef['uslope'] = np.real(m[-1])/lor
            coef['vslope'] = np.imag(m[-1])/lor
    else:
        if opt['notrend']:
            coef['mean'] = np.real(m[-1])
        else:
            coef['mean'] = np.real(m[-1-1])
            coef['slope'] = np.real(m[-1])/lor

    if opt['conf_int'] is True:
        coef = ut_confidence(coef, opt, t, e, tin, tgd, uvgd, elor, xraw, xmod,
                             W, m, B, nm, nt, nc, Xu, Yu, Xv, Yv)

    # diagnostics
    if not opt['nodiagn']:
        coef, indPE = ut_diagn(coef, opt)

    # re-order constituents
    if len(opt['ordercnstit']) != 0:

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
            '''There has to be a better way to do this.'''
            ind = np.zeros((len(opt['ordercnstit'])))
            for j, v in enumerate(opt['ordercnstit']):
                temp = np.core.defchararray.replace(coef['name'], " ", "")
                v = np.core.defchararray.replace(v, " ", "")
                lind1 = np.where(temp == v)[0][0]
                ind[j] = lind1
            ind = ind.astype(int)

    else:
        if not opt['nodiagn']:
            ind = indPE

        else:
            if opt['twodim']:
                PE = np.sum(coef['Lsmaj']**2 + coef['Lsmin']**2)
                PE = 100 * (coef['Lsmaj']**2 + coef['Lsmin']**2) / PE
            else:
                PE = 100 * coef['A']**2 / np.sum(coef['A']**2)

            ind = PE.argsort()[::-1]

    coef['g'] = coef['g'][ind]
    coef['name'] = coef['name'][ind]
    if opt['twodim']:
        coef['Lsmaj'] = coef['Lsmaj'][ind]
        coef['Lsmin'] = coef['Lsmin'][ind]
        coef['theta'] = coef['theta'][ind]
        if opt['conf_int'] is True:
            coef['Lsmaj_ci'] = coef['Lsmaj_ci'][ind]
            coef['Lsmin_ci'] = coef['Lsmin_ci'][ind]
            coef['theta_ci'] = coef['theta_ci'][ind]
            coef['g_ci'] = coef['g_ci'][ind]

    else:
        coef['A'] = coef['A'][ind]
        if opt['conf_int'] is True:
            coef['A_ci'] = coef['A_ci'][ind]
            coef['g_ci'] = coef['g_ci'][ind]

    coef['aux']['frq'] = coef['aux']['frq'][ind]
    coef['aux']['lind'] = coef['aux']['lind'][ind]

    print("Done.\n")

    return coef


def _slvinit(tin, uin, vin, **opts):

    if len(tin) != len(uin):
        raise ValueError('''solve: vectors of input times and
               input values must be same size.''')

    opt = {}

    tgd = ~np.isnan(tin)
    uin = uin[tgd]
    tin = tin[tgd]

    if len(vin) == 0:
        opt['twodim'] = False
        v = np.array([])
    else:
        if len(tin) != len(vin):
            raise('''solve: vectors of input times and
                  input values must be same size.''')

        opt['twodim'] = True
        vin = vin[tgd]

    if opt['twodim']:
        uvgd = ~np.isnan(uin) & ~np.isnan(vin)
        v = vin[uvgd]
    else:
        uvgd = ~np.isnan(uin)

    t = tin[uvgd]
    nt = len(t)
    u = uin[uvgd]
    eps = np.finfo(np.float64).eps

    if np.var(np.unique(np.diff(tin))) < eps:
        opt['equi'] = 1  # based on times; u/v can still have nans ("gappy")
        lor = (np.max(tin)-np.min(tin))
        elor = lor*len(tin)/(len(tin)-1)
        tref = 0.5*(tin[0]+tin[-1])
    else:
        opt['equi'] = 1  # based on times; u/v can still have nans ("gappy")
        lor = (np.max(tin)-np.min(tin))
        elor = lor*len(tin)/(len(tin)-1)
        tref = 0.5*(tin[0]+tin[-1])

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
    opt['ordercnstit'] = []
    opt['runtimedisp'] = 'yyy'

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

    return nt, t, u, v, tref, lor, elor, opt, tgd, uvgd


