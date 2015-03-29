"""
Single small function for ellipse parameters.
It is used only by solve().
This will probably end up being consolidated in another
module.
"""

from __future__ import absolute_import, division

import numpy as np


def ut_cs2cep(Xu, Yu, Xv=None, Yv=None):
    """
    Compute ellipse parameters from cosine and sine coefficients.

    For the 2-D case (currents), Xu, Yu are the cosine and sine
    coefficients of the zonal component, and Xv and Yv of the
    meridional component.

    For the 1-D case (height), Xv and Yv are None (default).

    Returns:
        Lsmaj, Lsmin: semi-major and semi-minor axes
        theta: major axis orientation, degrees ccw from x-axis, 0-180
        g: phase, degrees, 0-360

    1-D case: Lsmin and theta are zero

    Note: Following the matlab, Lsmin can be negative.

    % UTide v1p0 9/2011 d.codiga@gso.uri.edu
    """

    if Xv is None:
        ap = (Xu - 1j*Yu)
        Lsmaj = np.abs(ap)
        Lsmin = np.zeros_like(Xu)
        theta = np.zeros_like(Xu)
        g = -np.angle(ap, deg=True) % 360
        return Lsmaj, Lsmin, theta, g

    # 2-D case
    ap = ((Xu+Yv) + 1j*(Xv-Yu))/2
    am = ((Xu-Yv) + 1j*(Xv+Yu))/2
    Ap = np.abs(ap)
    Am = np.abs(am)
    Lsmaj = Ap+Am
    Lsmin = Ap-Am
    epsp = np.angle(ap, deg=True)
    epsm = np.angle(am, deg=True)

    theta = ((epsp+epsm)/2) % 180
    g = (-epsp+theta) % 360

    return Lsmaj, Lsmin, theta, g

