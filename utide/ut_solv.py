from __future__ import absolute_import, division

from .ut_solv1 import ut_solv1

# def ut_solv(tin, uin, vin, lat, cnstit, Rayleigh, *varargin):
#
#    coef = ut_solv1(tin,uin,vin,lat,cnstit,Rayleigh,varargin)
#
#    return coef


# def ut_solv(tin, uin, vin, lat, cnstit='auto', notrend=0, prefilt=[],
#            nodsatlint=0, nodsatnone=0, gwchlint=0, gwchnone=0, infer=[],
#            inferaprx=0, rmin=1, method='cauchy', tunrdn=1, linci=0, white=0,
#            nrlzn=200, lsfrqosmp=1, nodiagn=0, diagnplots=0, diagnminsnr=2,
#            ordercnstit=[], runtimedisp='yyy'):
#
#    coef = ut_solv1(tin, uin, vin, lat, **opts)
#
#    return coef


def ut_solv(tin, uin, vin, lat, **opts):
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

    coef = ut_solv1(tin, uin, vin, lat, **opts)

    return coef
