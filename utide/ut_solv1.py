import numpy as np
import scipy
from ut_slvinit import ut_slvinit
from ut_E import ut_E
from ut_cnstitsel import ut_cnstitsel
from ut_cs2cep import ut_cs2cep
from ut_confidence import ut_confidence

def ut_solv1(tin,uin,vin,lat,cnstit,Rayleigh,varargin):

    print 'ut_solv: '
    nt,t,u,v,tref,lor,elor,opt,tgd,uvgd = ut_slvinit(tin,uin,vin,cnstit,Rayleigh,varargin)

    opt['cnstit'] = cnstit
    [nNR,nR,nI,cnstit,coef] = ut_cnstitsel(tref, opt['rmin']/(24*lor),
                                           opt['cnstit'], opt['infer'])

    # a function we don't need
    # coef.aux.rundescr = ut_rundescr(opt,nNR,nR,nI,t,tgd,uvgd,lat)

    coef['aux']['opt'] = opt
    coef['aux']['lat'] = lat

    print 'matrix prep ... '

    ngflgs = [opt['nodsatlint'], opt['nodsatnone'],
              opt['gwchlint'], opt['gwchnone']]

    #ngflgs = [opt.nodsatlint, opt.nodsatnone opt.gwchlint opt.gwchnone];

    E = ut_E(t,tref,cnstit['NR']['frq'],cnstit['NR']['lind'],lat,ngflgs,opt['prefilt'])

    B = np.hstack((E,E.conj()))

    # more infer stuff

    if opt['notrend']:
        B = np.hstack((B,np.ones((nt,1))))
        nm = 2 * (nNR + nR) + 1
    else:
        B = np.hstack((B, np.ones((nt,1)), (t-tref)/lor))
        nm = 2*(nNR + nR) + 2

    print 'Solution ...'

    xraw = u

    if opt['twodim']:
        #xraw = complex(u,v);
        xraw = u+v*1j

    if opt['method']=='ols':
        #m = B\xraw;
        m = np.linalg.lstsq(B, xraw)[0]
        #W = sparse(1:nt,1:nt,1);
        W = scipy.sparse.identity(nt)
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

    xmod = B*m
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
        #XY = np.hstack((Xu, Yu))
        coef['A'], _ , _ ,coef['g '] = ut_cs2cep(Xu, Yu)
        #coef['A'], _ , _ ,coef['g '] = ut_cs2cep(XY)

    else:
        Xv = np.imag(ap+am)
        Yv = np.real(ap-am)
        #XY = np.vstack((Xu, Yu, Xv, Yv))
        coef['Lsmaj'], coef['Lsmin'], coef['theta'], coef['g'] = ut_cs2cep(Xu, Yu, Xv, Yv)
        #coef['Lsmaj'], coef['Lsmin'], coef['theta'], coef['g '] = ut_cs2cep(XY)

    ## mean and trend
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

    coef = ut_confidence(coef, opt, t, e, tin, tgd, uvgd, elor, xraw, xmod, W, m, B,
                   nm, nt, nc, Xu, Yu, Xv, Yv)

    if opt['twodim']:
        PE = np.sum(coef['Lsmaj']**2+coef['Lsmin']**2)
        PE = 100* (coef['Lsmaj']**2+coef['Lsmin']**2)/PE

    ind = PE.argsort()[::-1]
    coef['Lsmaj'] = coef['Lsmaj'][ind]
    coef['Lsmin'] = coef['Lsmin'][ind]
    coef['theta'] = coef['theta'][ind]
    coef['g'] = coef['g'][ind]
    coef['name'] = coef['name'][ind]

    coef['aux']['frq'] = coef['aux']['frq'][ind]
    coef['aux']['lind'] = coef['aux']['lind'][ind]

    return coef
