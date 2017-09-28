"""
Function ut_E() returns complex exponential basis functions
for a given set of frequencies.
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np

from .astronomy import ut_astron
from ._ut_constants import ut_constants


sat = ut_constants.sat
const = ut_constants.const
shallow = ut_constants.shallow

nshallow = np.ma.masked_invalid(const.nshallow).astype(int)
ishallow = np.ma.masked_invalid(const.ishallow).astype(int) - 1
not_shallow = ishallow.mask  # True where it was masked.
nshallow = nshallow.compressed()
ishallow = ishallow.compressed()
kshallow = np.nonzero(~not_shallow)[0]


def linearized_freqs(tref):
    astro, ader = ut_astron(tref)
    freq = const.freq.copy()
    selected = np.dot(const.doodson[not_shallow, :], ader) / 24
    freq[not_shallow] = selected.squeeze()
    for i0, nshal, k in zip(ishallow, nshallow, kshallow):
        ik = i0 + np.arange(nshal)
        freq[k] = (freq[shallow.iname[ik] - 1] *
                   shallow.coef[ik]).sum()
    return freq


def ut_E(t, tref, frq, lind, lat, ngflgs, prefilt):
    """
    Compute complex exponential basis function.

    Parameters
    ----------
    t : array_like or float (nt,)
        time in days
    tref : float
        reference time in days
    frq : array_like or float (nc,)
        frequencies in cph
    lind : array_like or int (nc,)
        indices of constituents
    lat : float
        latitude, degrees N
    nflgs : array_like, bool
        [NodsatLint NodsatNone GwchLint GwchNone]
    prefilt: Bunch
        not implemented

    Returns
    -------
    E : array (nt, nc)
        complex exponential basis function; always returned as 2-D array
    """

    t = np.atleast_1d(t)
    frq = np.atleast_1d(frq)
    lind = np.atleast_1d(lind)
    nt = len(t)
    nc = len(frq)

    if ngflgs[1] and ngflgs[3]:
        F = np.ones((nt, nc))
        U = np.zeros((nt, nc))
        V = np.dot(24*(t-tref)[:, None], frq[:, None].T)
    else:
        F, U, V = FUV(t, tref, lind, lat, ngflgs)

    E = F * np.exp(1j*(U+V)*2*np.pi)

    # if ~isempty(prefilt)
    # if len(prefilt)!=0:
    #     P=interp1(prefilt.frq,prefilt.P,frq).T
    #     P( P>max(prefilt.rng) | P<min(prefilt.rng) | isnan(P) )=1;
    #     E = E*P(ones(nt,1),:);

    return E


def FUV(t, tref, lind, lat, ngflgs):
    """
    UT_FUV()
    compute nodal/satellite correction factors and astronomical argument
    inputs
      t = times [datenum UTC] (nt x 1)
      tref = reference time [datenum UTC] (1 x 1)
      lind = list indices of constituents in ut_constants.mat (nc x 1)
      lat = latitude [deg N] (1 x 1)
      ngflgs = [NodsatLint NodsatNone GwchLint GwchNone] each 0/1
    output
      F = real nodsat correction to amplitude [unitless] (nt x nc)
      U = nodsat correction to phase [cycles] (nt x nc)
      V = astronomical argument [cycles] (nt x nc)
    UTide v1p0 9/2011 d.codiga@gso.uri.edu
    (uses parts of t_vuf.m from t_tide, Pawlowicz et al 2002)
    """

    t = np.atleast_1d(t).flatten()
    nt = len(t)
    nc = len(lind)

    # nodsat

    if ngflgs[1]:
        F = np.ones((nt, nc))
        U = np.zeros((nt, nc))
    else:
        if ngflgs[0]:
            tt = np.array([tref])
        else:
            tt = t
        ntt = len(tt)

        astro, ader = ut_astron(tt)

        if abs(lat) < 5:
            lat = np.sign(lat) * 5

        slat = np.sin(np.deg2rad(lat))
        rr = sat.amprat.copy()

        j = sat.ilatfac == 1
        rr[j] *= 0.36309 * (1.0 - 5.0 * slat**2)/slat

        j = sat.ilatfac == 2
        rr[j] *= 2.59808 * slat

        # sat.deldood is (162, 3); all other sat vars are (162,)
        uu = np.dot(sat.deldood, astro[3:6, :]) + sat.phcorr[:, None]
        np.fmod(uu, 1, out=uu)  # fmod is matlab rem; differs from % op
        mat = rr[:, None] * np.exp(1j * 2 * np.pi * uu)

        nfreq = len(const.isat)  # 162
        F = np.ones((nfreq, ntt), dtype=complex)

        iconst = sat.iconst - 1
        ind = np.unique(iconst)
        for ii in ind:
            F[ii, :] = 1 + np.sum(mat[iconst == ii], axis=0)

        U = np.angle(F) / (2 * np.pi)  # cycles
        F = np.abs(F)

        for i0, nshal, k in zip(ishallow, nshallow, kshallow):
            ik = i0 + np.arange(nshal)
            j = shallow.iname[ik] - 1
            exp1 = shallow.coef[ik, None]
            exp2 = np.abs(exp1)
            F[k, :] = np.prod(F[j, :]**exp2, axis=0)
            U[k, :] = np.sum(U[j, :] * exp1, axis=0)

        F = F[lind, :].T
        U = U[lind, :].T

        # if ngflgs[0]:  # Nodal/satellite with linearized times.
        #    F = F[np.ones((nt, 1)), :]
        #    U = U[np.ones((nt, 1)), :]
        # Let's try letting broadcasting take care of it.

    # gwch (astron arg)
    if ngflgs[3]:  # None (raw phase lags not Greenwich phase lags).
        freq = linearized_freqs(tref)
        V = 24 * (t[:, np.newaxis] - tref) * freq[lind]

    else:
        if ngflgs[2]:  # Linearized times.
            tt = np.array([tref])
        else:
            tt = t  # Exact times.
        ntt = len(tt)

        astro, ader = ut_astron(tt)

        V = np.dot(const.doodson, astro) + const.semi[:, None]
        np.fmod(V, 1, out=V)

        for i0, nshal, k in zip(ishallow, nshallow, kshallow):
            ik = i0 + np.arange(nshal)
            j = shallow.iname[ik] - 1
            exp1 = shallow.coef[ik, None]
            V[k, :] = np.sum(V[j, :] * exp1, axis=0)

        V = V[lind, :].T

        if ngflgs[2]:  # linearized times
            freq = linearized_freqs(tref)
            V = V + 24*(t[:, None] - tref) * freq[None, lind]

    return F, U, V
