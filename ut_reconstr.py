from ut_solv import ut_E, ut_FUV, ut_astron
import numpy as np

def ut_reconstr(tin,coef,varargin):

    u, v = ut_reconstr1(tin, coef, varargin)

    return u, v


def ut_reconstr1(tin, coef, varargin):

    print 'ut_reconstr: '

    # parse inputs and options
    t, opt = ut_rcninit(tin,varargin)

    # determine constituents to include
    #if ~isempty(opt.cnstit)
    if not np.empty(opt['cnstit']):

        #[~,ind] = ismember(cellstr(opt.cnstit),coef.name);
        #opt['cnstit'] in coef['name']
        ind = np.where(opt['cnstit'] == coef['name'])

#        if ~isequal(length(ind),length(cellstr(opt.cnstit)))
#            error(['ut_reconstr: one or more of input constituents Cnstit '...
#                'not found in coef.name']);
    else:
        #ind = 1:length(coef.aux.frq) # not needed
        if coef['aux']['opt']['twodim']:
        #if coef.aux.opt.twodim:
            SNR = (coef['Lsmaj']**2 +coef['Lsmin']**2)/((coef['Lsmaj_ci']/1.96)**2 + (coef['Lsmin_ci']/1.96)**2)
            PE = sum(coef['Lsmaj']**2 + coef['Lsmin']**2)
            PE = 100*(coef['Lsmaj']**2 + coef['Lsmin']**2)/PE
        else:
            SNR = (coef['A']**2)/((coef['A_ci']/1.96)**2)
            PE = 100*coef['A']**2/sum(coef['A']**2)

        ind = ind[SNR[ind]>=opt['minsnr ']& PE[ind]>=opt['minpe']]


    # complex coefficients
    rpd = np.pi/180
    #if coef.aux.opt.twodim
    if coef['aux']['opt']['twodim']:
        ap = 0.5*(coef['Lsmaj'][ind] + coef['Lsmin'][ind]) * np.exp(1j*(coef['theta'][ind] - coef['g'][ind])*rpd)
        am = 0.5*(coef['Lsmaj'][ind] - coef['Lsmin'][ind]) * np.exp(1j*(coef['theta'][ind] + coef['g'][ind])*rpd)
    else:
        ap = 0.5*coef['A'][ind]*np.exp(-1j*coef['g'][ind]*rpd)
        am = np.conj(ap)

    # exponentials
#    ngflgs = [coef.aux.opt.nodsatlint coef.aux.opt.nodsatnone ...
#        coef.aux.opt.gwchlint coef.aux.opt.gwchnone];

    ngflgs = [coef['aux']['opt']['nodsatlint'],coef['aux']['opt']['nodsatnone'],
              coef['aux']['opt']['gwchlint'],coef['aux']['opt']['gwchnone']]

    print 'prep/calcs ... '

#    E = ut_E(t,coef.aux.reftime,coef.aux.frq(ind),coef.aux.lind(ind),...
#        coef.aux.lat,ngflgs,coef.aux.opt.prefilt);

    E = ut_E(t,
             coef['aux']['reftime'],coef['aux']['frq'][ind],
             coef['aux']['lind'][ind],coef['aux']['lat'], ngflgs,
             coef['aux']['opt']['prefilt'] )

    # fit
    fit = E*ap + np.conj(E)*am

    # mean (& trend)
    u = np.nan*np.ones(tin.shape)
    whr = ~np.isnan(tin)
    if coef['aux']['opt']['twodim']:
        v = u
        if coef['aux']['opt']['notrend']:
        #if coef.aux.opt.notrend
            u[whr] = np.real(fit) + coef['umean']
            v[whr] = np.imag(fit) + coef['vmean']
        else:
            u[whr] = np.real(fit) + coef['umean ']+ coef['uslope']*(t-coef['aux']['reftime'])
            v[whr] = np.imag(fit) + coef['vmean ']+ coef['vslope']*(t-coef['aux']['reftime'])

    else:
        if coef['aux']['opt']['notrend']:
        #if coef.aux.opt.notrend
            u[whr] = np.real(fit) + coef['mean']
        else:
            u[whr] = np.real(fit) + coef['mean ']+ coef['slope']*(t-coef['aux']['reftime'])

        v = []


    print 'Done.\n'

    return u, v


def ut_rcninit(tin,args):

    t = tin[:]

    t[np.isnan(t)] = []
    #t(isnan(t)) = []
    opt = {}

    opt['cnstit'] = []
    opt['minsnr'] = 2
    opt['minpe'] = 0

    #args = list(args)
    #args = [string.lower() for string in args]

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
#                error(['ut_reconstr: unrecognized input: ' args{1}]);
#        end
#    end

    return t, opt
