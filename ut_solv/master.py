import numpy as np
import scipy.io as sio
import scipy.sparse

def ut_solv(tin, uin, vin, lat, cnstit, Rayleigh, *varargin):

    coef = ut_solv1(tin,uin,vin,lat,cnstit,Rayleigh,varargin)

    return coef


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


#def ut_cs2cep(XY):
def ut_cs2cep(Xu, Yu, Xv=np.array([False]), Yv=np.array([False])):

    #Xu = XY[:, 0]
    #Yu = XY[:, 1]

    if not Xv.all():
        Xv = np.zeros(Xu.shape)
        Yv = np.zeros(Yu.shape)

#    if XY.shape[-1] > 2:
#        Xv = XY[:, 3]
#        Yv = XY[:, 4]
#    else:
#        Xv = np.zeros(Xu.shape)
#        Yv = np.zeros(Yu.shape)

    ap = ((Xu+Yv)+1j*(Xv-Yu))/2
    am = ((Xu-Yv)+1j*(Xv+Yu))/2
    Ap = np.abs(ap)
    Am = np.abs(am)
    Lsmaj = Ap+Am
    Lsmin = Ap-Am
    epsp = np.angle(ap)*180/np.pi
    epsm = np.angle(am)*180/np.pi

    theta = ((epsp+epsm)/2) % 180
    g = (-epsp+theta) % 360

    return Lsmaj,Lsmin,theta,g


def ut_E(t,tref,frq,lind,lat,ngflgs,prefilt):

    nt = len(t)
    nc = len(lind)
    if ngflgs[1] and ngflgs[3]:
        F = np.ones((nt,nc))
        U = np.zeros((nt,nc))
        V = 24*(t-tref)*frq
    else:
        F, U, V = ut_FUV(t,tref,lind,lat,ngflgs);

    E = F * np.exp(1j*(U+V)*2*np.pi)

    #if ~isempty(prefilt)
#    if len(prefilt)!=0:
#        P=interp1(prefilt.frq,prefilt.P,frq).T
#        P( P>max(prefilt.rng) | P<min(prefilt.rng) | isnan(P) )=1;
#        E = E*P(ones(nt,1),:);

    return E


def ut_FUV(t, tref, lind, lat, ngflgs):

    nt = len(t)
    nc = len(lind)
    ## nodsat

    if ngflgs[1]:
        F = np.ones((nt,nc))
        U = np.zeros((nt,nc))
    else:
        if ngflgs[0]:
            tt = tref
        else:
            tt = t

        ntt = len(tt)

        mat_contents = sio.loadmat('ut_constants.mat', struct_as_record=False, squeeze_me=True)
        sat = mat_contents['sat']
        const = mat_contents['const']
        shallow = mat_contents['shallow']

        astro, ader = ut_astron(tt)

        if abs(lat) < 5:
            lat=np.sign(lat)*5;

        slat = np.sin(np.pi * lat/180)
        rr = sat.amprat
        j = np.where(sat.ilatfac==1)[0]

        rr[j] = rr[j]*0.36309*(1.0-5.0*slat*slat)/slat;

        j = np.where(sat.ilatfac==2)

        rr[j]=rr[j]*2.59808*slat;

        uu = np.dot(sat.deldood, astro[3:6, :]) + sat.phcorr[:,None]*np.ones((1,ntt)) % 1

        nfreq=len(const.isat)
        mat = rr[:,None]*np.ones((1,ntt)) * np.exp(1j*2*np.pi*uu)

        F = np.ones((nfreq, ntt)) + 0j
        ind = np.unique(sat.iconst)

        for i in xrange(len(ind)):
            F[ind[i]-1, :] = 1+np.sum(mat[sat.iconst==ind[i],:], axis=0)

        #U = imag(log(F))/(2*pi); % faster than angle(F)
        U = np.imag(np.log(F)) / (2*np.pi)
        F = np.abs(F)

        for k in np.where(np.isfinite(const.ishallow))[0]:
            ik=const.ishallow[k]+np.arange(const.nshallow[k])
            ik = ik.astype(int)
            j = shallow.iname[ik-1]
            exp1 = shallow.coef[ik-1]
            exp2 = np.abs(exp1)
            temp1 = exp1*np.ones((ntt,1))
            temp2 = exp2*np.ones((ntt,1))
            temp1 = temp1.T
            temp2 = temp2.T
            F[k,:]=np.prod(F[j-1,:]**temp2,axis=0)
            U[k,:]=np.sum(U[j-1,:]*temp1,axis=0)

        F=F[lind,:].T
        U=U[lind,:].T

        if ngflgs[1]: # nodal/satellite with linearized times
            F = F[np.ones((nt,1)),:]
            U = U[np.ones((nt,1)),:]

    ## gwch (astron arg)
    if ngflgs[3]: # none (raw phase lags not greenwich phase lags)
#        if ~exist('const','var'):
#            load('ut_constants.mat','const');
#        [~,ader] = ut_astron(tref);
#        ii=isfinite(const.ishallow);
#        const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
#        for k=find(ii)'
#            ik=const.ishallow(k)+(0:const.nshallow(k)-1);
#            const.freq(k)=sum(const.freq(shallow.iname(ik)).*shallow.coef(ik));
        V = 24*(t-tref)*const.freq(lind).T
    else:
        if ngflgs[3]: # linearized times
            tt = tref
        else:
            tt = t # exact times

        ntt = len(tt)
#        if exist('astro','var')
#            if ~isequal(size(astro,2),ntt)
#                [astro,~]=ut_astron(tt');
#            end
#        else
#            [astro,~]=ut_astron(tt');
#
#        if ~exist('const','var')
#            load('ut_constants.mat')

        mat_contents = sio.loadmat('ut_constants.mat',
                                   struct_as_record=False, squeeze_me=True)
        sat = mat_contents['sat']
        const = mat_contents['const']
        shallow = mat_contents['shallow']
        astro, ader = ut_astron(tt)

        #V = np.dot(const.doodson, astro) + const.semi[:,None]*np.ones((1,ntt)) % 1
        V = np.dot(const.doodson, astro) + const.semi[:,None]*np.ones((1,ntt))
        #V = V % 1

        for k in np.where(np.isfinite(const.ishallow))[0]:
            ik=const.ishallow[k]+np.arange(const.nshallow[k])
            ik = ik.astype(int)
            j = shallow.iname[ik-1]
            exp1 = shallow.coef[ik-1]
            temp1 = exp1[:]*np.ones((ntt,1))
            temp1 = temp1.T

            V[k,:] = np.sum(V[j-1,:]*temp1,axis=0)

        V=V[lind,:].T

#        if ngflgs(3) % linearized times
#            [~,ader] = ut_astron(tref);
#            ii=isfinite(const.ishallow);
#            const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
#            for k=find(ii)'
#                ik=const.ishallow(k)+(0:const.nshallow(k)-1);
#                const.freq(k)=sum( const.freq(shallow.iname(ik)).* ...
#                    shallow.coef(ik) );
#            end
#            V = V(ones(1,nt),:) + 24*(t-tref)*const.freq(lind)';
#        end

    return F, U, V


def ut_cnstitsel(tref,minres,incnstit,infer):

    mat_contents = sio.loadmat('ut_constants.mat', struct_as_record=False, squeeze_me=True)
    shallow = mat_contents['shallow']
    const = mat_contents['const']

    cnstit = {}
    coef = {}

    astro, ader = ut_astron(tref)

    ii = np.isfinite(const.ishallow)
    const.freq[~ii] = np.dot(const.doodson[~ii, :], ader) / 24

    for k in ii.nonzero()[0]:
        ik = const.ishallow[k]+np.arange(const.nshallow[k])
        ik = ik.astype(int)-1
        const.freq[k] = np.sum(const.freq[shallow.iname[ik]]*shallow.coef[ik])

    ## cnstit.NR
    cnstit['NR'] = {}
    if incnstit.lower()=='auto':
        cnstit['NR']['lind'] = np.where(const.df >= minres)[0]
    else:
        pass

    # skipped some stuff here cause they involve infer

    cnstit['NR']['frq'] = const.freq[cnstit['NR']['lind']]
    cnstit['NR']['name'] = const.name[cnstit['NR']['lind']]
    nNR = len(cnstit['NR']['frq'])

    ## cnstit.R
    nR = 0
    nI = 0
    cnstit['R'] = []

    nallc = nNR+nR+nI

    coef['name'] = cnstit['NR']['name']
    coef['aux'] = {}
    coef['aux']['frq'] = cnstit['NR']['frq']
    coef['aux']['lind'] = cnstit['NR']['lind']

    # another infer if statement

    coef['aux']['reftime'] = tref


    return nNR, nR, nI, cnstit, coef


def ut_slvinit(tin,uin,vin,cnstit,Rayleigh,args):

    opt = {}
    args = list(args)
    tgd = ~np.isnan(tin)
    uin = uin[tgd]
    tin = tin[tgd]

    if vin.shape[0] == 0:
        opt['twodim'] = False
        #twodim = False
        v = np.array([])
    else:
        opt['twodim'] = True
        #twodim = True
        vin = vin[tgd]

    #if twodim:
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
        opt['equi'] = 1 # based on times; u/v can still have nans ("gappy")
        #equi = 1 # based on times; u/v can still have nans ("gappy")
        lor = (np.max(tin)-np.min(tin))
        elor = lor*len(tin)/(len(tin)-1)
        tref = 0.5*(tin[0]+tin[-1])
    else:
        opt['equi'] = 0
        #equi = 0;
        lor = (np.max(t) - np.min(t))
        elor = lor*nt/(nt-1)
        tref = 0.5*(t[0]+t[-1])

    ## options
    opt['notrend'] = 0
    opt['prefilt'] = []
    opt['nodsatlint'] = 0
    opt['nodsatnone'] = 0
    opt['gwchlint'] = 0
    opt['gwchnone'] = 0
    opt['infer'] = []
    opt['inferaprx'] = 0
    opt['rmin'] = 1
    opt['method'] = 'cauchy'
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

    #methnotset = 1
    allmethods = ['ols', 'andrews', 'bisquare', 'fair', 'huber',
                  'logistic', 'talwar', 'welsch']

    args = [string.lower() for string in args]

    if 'notrend' in args:
        #opt['notrend'] = 1
        opt['notrend'] = True

    if 'rmin' in args:
        opt['rmin'] = Rayleigh

    if 'nodiagn' in args:
        #opt['nodiagn']=1
        opt['nodiagn'] = True

    if 'linci' in args:
        #opt['linci'] = 1
        opt['linci'] = True

    if allmethods:
        methods = [i for i in allmethods if i in args]
        if len(methods) > 1:
            print 'ut_solv: Only one "method" option allowed.'
        else:
            opt['method'] = methods[0]

    if opt['method'] != 'cauchy':
        ind = np.argwhere(opt['method'] in allmethods)[0][0]
        allconst = [np.nan, 1.339, 4.685, 1.400, 1.345, 1.205, 2.795, 2.985]
        opt['tunconst'] = allconst[ind]
    else:
        opt['tunconst'] = 2.385

    opt['tunconst'] = opt['tunconst'] /opt['tunrdn']

# only needed if we sort the options
#    nf = len(opt)

    return nt, t, u, v, tref, lor, elor, opt, tgd, uvgd


def ut_astron(jd):
    '''
    UT_ASTRON()
    calculate astronomical constants
    input
    jd = time [datenum UTC] (1 x nt)
    outputs
    astro = matrix [tau s h p np pp]T, units are [cycles] (6 x nt)
    ader = matrix of derivatives of astro [cycles/day] (6 x nt)
    UTide v1p0 9/2011 d.codiga@gso.uri.edu
    (copy of t_astron.m from t_tide, Pawlowicz et al 2002)
    '''

    jd = np.array([jd])
    # datenum(1899,12,31,12,0,0)
    daten = 693961.500000000
    d = jd[:] - daten
    D = d / 10000

    #args = np.array([[np.ones(jd.shape)],[d],[D*D],[D**3]]).flatten()[:,None]
    args = np.vstack((np.ones(jd.shape), d, D*D, D**3))
    sc = np.array([270.434164, 13.1763965268, -0.0000850, 0.000000039])
    hc= np.array([ 279.696678, 0.9856473354, 0.00002267,0.000000000])
    pc= np.array([ 334.329556, 0.1114040803,-0.0007739,-0.00000026])
    npc= np.array([-259.183275, 0.0529539222,-0.0001557,-0.000000050])
    ppc= np.array([281.220844, 0.0000470684, 0.0000339, 0.000000070])

    astro = np.dot(np.vstack((sc, hc, pc, npc, ppc)), args) / 360 % 1
    tau = jd%1 + astro[1,:] - astro[0,:]
    astro = np.vstack((tau,astro))

#    dargs = np.array([[np.zeros(jd.shape[0])], [np.ones(jd.shape[0])],
#                      [2.0e-4*D], [3.0e-4*D*D]]).flatten()[:,None]

    dargs = np.vstack((np.zeros(jd.shape), np.ones(jd.shape),
                      2.0e-4*D, 3.0e-4*D*D))

    ader = np.dot(np.vstack((sc,hc,pc,npc,ppc)), dargs)/360.0
    dtau = 1.0 + ader[1,:] - ader[0, :]
    ader = np.vstack((dtau, ader))

    # might need to take out depending on shape of jd
    #astro = astro.flatten()
    #ader = ader.flatten()

    return astro,ader


def loadMAT(filename):

    mat_contents = sio.loadmat('ut_constants.mat', struct_as_record=False, squeeze_me=True)

    items = []
    items = {}

    for i in mat_contents:
        name = '{0}'.format(i)
        items[name] = mat_contents[name]
        #i = mat_contents[name]
        #items.append(i)

    return items
