from ut_solv1 import ut_solv1

def ut_solv(tin, uin, vin, lat, cnstit, Rayleigh, *varargin):

    coef = ut_solv1(tin,uin,vin,lat,cnstit,Rayleigh,varargin)

    return coef
