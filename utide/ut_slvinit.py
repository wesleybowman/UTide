import numpy as np

#def ut_slvinit(tin,uin,vin,cnstit,Rayleigh,args):
#
#    opt = {}
#    args = list(args)
#    tgd = ~np.isnan(tin)
#    uin = uin[tgd]
#    tin = tin[tgd]
#
#    #if vin.shape[0] == 0:
#    if len(vin) == 0:
#        opt['twodim'] = False
#        v = np.array([])
#    else:
#        opt['twodim'] = True
#        #twodim = True
#        vin = vin[tgd]
#
#    #if twodim:
#    if opt['twodim']:
#        uvgd = ~np.isnan(uin) & ~np.isnan(vin)
#        v = vin[uvgd]
#    else:
#        uvgd = ~np.isnan(uin)
#
#    t = tin[uvgd]
#    nt = len(t)
#    u = uin[uvgd]
#    eps = np.finfo(np.float64).eps
#
#    if np.var(np.unique(np.diff(tin))) < eps:
#        opt['equi'] = 1 # based on times; u/v can still have nans ("gappy")
#        #equi = 1 # based on times; u/v can still have nans ("gappy")
#        lor = (np.max(tin)-np.min(tin))
#        elor = lor*len(tin)/(len(tin)-1)
#        tref = 0.5*(tin[0]+tin[-1])
#    else:
#        opt['equi'] = 0
#        #equi = 0;
#        lor = (np.max(t) - np.min(t))
#        elor = lor*nt/(nt-1)
#        tref = 0.5*(t[0]+t[-1])
#
#    ## options
#    opt['notrend'] = 0
#    opt['prefilt'] = []
#    opt['nodsatlint'] = 0
#    opt['nodsatnone'] = 0
#    opt['gwchlint'] = 0
#    opt['gwchnone'] = 0
#    opt['infer'] = []
#    opt['inferaprx'] = 0
#    opt['rmin'] = 1
#    opt['method'] = 'cauchy'
#    opt['tunrdn'] = 1
#    opt['linci'] = 0
#    opt['white'] = 0
#    opt['nrlzn'] = 200
#    opt['lsfrqosmp'] = 1
#    opt['nodiagn'] = 0
#    opt['diagnplots'] = 0
#    opt['diagnminsnr'] = 2
#    opt['ordercnstit'] = []
#    opt['runtimedisp'] = 'yyy'
#
#    #methnotset = 1
#    allmethods = ['ols', 'andrews', 'bisquare', 'fair', 'huber',
#                  'logistic', 'talwar', 'welsch']
#
#    args = [string.lower() for string in args]
#
#    if 'notrend' in args:
#        #opt['notrend'] = 1
#        opt['notrend'] = True
#
#    if 'rmin' in args:
#        opt['rmin'] = Rayleigh
#
#    if 'nodiagn' in args:
#        #opt['nodiagn']=1
#        opt['nodiagn'] = True
#
#    if 'linci' in args:
#        #opt['linci'] = 1
#        opt['linci'] = True
#
#    if allmethods:
#        methods = [i for i in allmethods if i in args]
#        if len(methods) > 1:
#            print 'ut_solv: Only one "method" option allowed.'
#        else:
#            opt['method'] = methods[0]
#
#    if opt['method'] != 'cauchy':
#        ind = np.argwhere(opt['method'] in allmethods)[0][0]
#        allconst = [np.nan, 1.339, 4.685, 1.400, 1.345, 1.205, 2.795, 2.985]
#        opt['tunconst'] = allconst[ind]
#    else:
#        opt['tunconst'] = 2.385
#
#    opt['tunconst'] = opt['tunconst'] /opt['tunrdn']
#
## only needed if we sort the options
##    nf = len(opt)
#
#    return nt, t, u, v, tref, lor, elor, opt, tgd, uvgd


def ut_slvinit(tin, uin, vin, **opts):

    opt = {}

    tgd = ~np.isnan(tin)
    uin = uin[tgd]
    tin = tin[tgd]

    #if vin.shape[0] == 0:
    if len(vin) == 0:
        opt['twodim'] = False
        v = np.array([])
    else:
        opt['twodim'] = True
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
        lor = (np.max(tin)-np.min(tin))
        elor = lor*len(tin)/(len(tin)-1)
        tref = 0.5*(tin[0]+tin[-1])
    else:
        opt['equi'] = 0
        lor = (np.max(t) - np.min(t))
        elor = lor*nt/(nt-1)
        tref = 0.5*(t[0]+t[-1])

    ## options
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

    for key, item in opts.items():
        opt[key] = item

    #methnotset = 1
    allmethods = ['ols', 'andrews', 'bisquare', 'fair', 'huber',
                  'logistic', 'talwar', 'welsch']
#
#    args = [string.lower() for string in args]
#
#    if 'notrend' in args:
#        #opt['notrend'] = 1
#        opt['notrend'] = True
#
#    if 'rmin' in args:
#        opt['rmin'] = Rayleigh
#
#    if 'nodiagn' in args:
#        #opt['nodiagn']=1
#        opt['nodiagn'] = True
#
#    if 'linci' in args:
#        #opt['linci'] = 1
#        opt['linci'] = True
#
#    if allmethods:
#        methods = [i for i in allmethods if i in args]
#        if len(methods) > 1:
#            print 'ut_solv: Only one "method" option allowed.'
#        else:
#            opt['method'] = methods[0]

    if opt['method'] != 'cauchy':
        ind = np.argwhere(opt['method'] in allmethods)[0][0]
        allconst = [np.nan, 1.339, 4.685, 1.400, 1.345, 1.205, 2.795, 2.985]
        opt['tunconst'] = allconst[ind]
    else:
        opt['tunconst'] = 2.385

    opt['tunconst'] = opt['tunconst'] /opt['tunrdn']

    return nt, t, u, v, tref, lor, elor, opt, tgd, uvgd


#cnstit='auto',
#notrend=0,
#prefilt=[],
#nodsatlint=0,
#nodsatnone=0,
#gwchlint=0,
#gwchnone=0,
#infer=[],
#inferaprx=0,
#rmin=1,
#method='cauchy',
#tunrdn=1,
#linci=0,
#white=0,
#nrlzn=200,
#lsfrqosmp=1,
#nodiagn=0,
#diagnplots=0,
#diagnminsnr=2,
#ordercnstit=[],
#runtimedisp='yyy'
