from __future__ import absolute_import, division

import numpy as np


def ut_cs2cep(Xu, Yu, Xv=np.array([False]), Yv=np.array([False])):

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
