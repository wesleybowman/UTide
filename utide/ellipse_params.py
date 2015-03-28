"""
Single small function for ellipse parameters.
It is used only by solve().
This will probably end up being consolidated in another
module.
"""

from __future__ import absolute_import, division

import numpy as np


def ut_cs2cep(Xu, Yu, Xv=np.array([False]), Yv=np.array([False])):
"""
% UT_CS2CEP()
% compute current ellipse parameters from cosine-sine coefficients
% inputs
%   two-dim case: XY = [Xu Yu Xv Yv] 4-column matrix
%   one-dim case: XY = [Xu Yu] 2-column matrix
%                      OR 4-column w/ Xv=Yv=zeros(size(Xu))
%      where: Xu,Yu are cosine, sine coeffs of u, & same for v
% outputs
%   two-dim case:
%     Lsmaj, Lsmin = column vectors [units of XY] (size of Xu)
%     theta = column vector [deg. ccw rel. +x-axis, 0-180] (size of Xu)
%     g = column vector [degrees, 0-360] (size of Xu)
%   one-dim case: same, where Lsmaj = A, and Lsmin and theta are zeros
% UTide v1p0 9/2011 d.codiga@gso.uri.edu
"""
    if not Xv.all():
        Xv = np.zeros(Xu.shape)
        Yv = np.zeros(Yu.shape)

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

    return Lsmaj, Lsmin, theta, g
