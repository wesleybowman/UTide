"""
Single small function for ellipse parameters.
It is used only by solve().
This will probably end up being consolidated in another
module.
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np


def ut_cs2cep(Xu, Yu=None, Xv=None, Yv=None):
    """
    Compute ellipse parameters from cosine and sine coefficients.

    Parameters
    ----------
    Xu : ndarray
        1-D array with cosine coefficients; or 2-D, (n,2) or (n,4)
        with columns being cosine and sine coefficients of height,
        or of velocity components.

    Yu, Xv, Yv : ndarray, 1-D, optional
        If Xu is 1-D, these will have the sine coefficients of
        the first variable (height or U velocity component), and
        potentially the cosine and sin coefficients of a second
        variable (V velocity component).

    Returns
    -------
        Lsmaj, Lsmin : ndarray
             semi-major and semi-minor axes
        theta : ndarray
             major axis orientation, degrees ccw from x-axis, 0-180
        g : ndarray
            phase, degrees, 0-360

    In the single-variable case (height): Lsmin and theta are zero

    Notes
    -----
    Following the matlab version, Lsmin can be negative.

    Based on UTide v1p0 9/2011 d.codiga@gso.uri.edu
    """

    if Yu is None:
        ndim = Xu.shape[-1]
        if not (ndim == 2 or ndim == 4):
            raise ValueError("invalid arguments")
        Yu = Xu[:, 1]
        if ndim == 4:
            Xv, Yv = Xu[:, 2:].T
        Xu = Xu[:, 0]

    if Xv is None:
        ap = Xu - 1j*Yu
        Lsmaj = np.abs(ap)
        Lsmin = np.zeros_like(Xu)
        theta = np.zeros_like(Xu)
        g = -np.angle(ap, deg=True) % 360
        return Lsmaj, Lsmin, theta, g

    # 2-D case
    ap = ((Xu+Yv) + 1j * (Xv-Yu))/2
    am = ((Xu-Yv) + 1j * (Xv+Yu))/2
    Ap = np.abs(ap)
    Am = np.abs(am)
    Lsmaj = Ap + Am
    Lsmin = Ap - Am
    epsp = np.angle(ap, deg=True)
    epsm = np.angle(am, deg=True)

    theta = ((epsp+epsm)/2) % 180
    g = (-epsp+theta) % 360

    return Lsmaj, Lsmin, theta, g
