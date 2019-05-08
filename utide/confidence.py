"""
Confidence interval calculations for solve().

Includes a strictly private _confidence() function and
a more general ut_linci() function for linearized estimates
of ellipse parameter uncertainties.
"""

from __future__ import (absolute_import, division, print_function)

import warnings

import numpy as np

from utide.ellipse_params import ut_cs2cep
from utide.periodogram import band_psd
from utide.utilities import complex_interp


def band_averaged_psd_by_constit(tin, t, e, elor, coef, opt):
    # Band-averaged (ba) spectral densities at each constituent freq.
    constits = coef.aux.frq
    if opt.equi:
        e_ = e
        if len(tin) > len(t):
            e_ = complex_interp(tin, t, e)
        ba = band_psd(tin, e_, constits, equi=True)

    else:
        ba = band_psd(t, e, constits,
                      equi=False, frqosamp=opt.lsfrqosmp)

    # power [ (e units)^2 ] from spectral density [ (e units)^2 / cph ]
    df = 1 / (elor * 24)  # inverse of record length in hours
    ba.Puu *= df
    Puu = np.zeros_like(constits)
    Pvv = Puv = None

    if opt.twodim:
        ba.Pvv *= df
        ba.Puv *= df
        Pvv = np.zeros_like(constits)
        Puv = np.zeros_like(constits)

    for i, (lo, hi) in enumerate(ba.fbnd):
        inside = (constits >= lo) & (constits <= hi)
        Puu[inside] = ba.Puu[i]
        if opt.twodim:
            Pvv[inside] = ba.Pvv[i]
            Puv[inside] = ba.Puv[i]
    return Puu, Pvv, Puv


def cluster(x, ang=360):
    """
    Wrapping the values of x to +- ang/2 of x[0].
    """

    x = np.array(x)
    ha = ang/2
    ofs = - x[0] + ha
    y = (x+ofs) % ang - ofs
    return y


def _is_PD(A):
    """
    Helper for nearestSPD.  Testing PD via the cholesky call is
    much faster than testing for negative eigenvalues.
    """
    try:
        np.linalg.cholesky(A)
        return True
    except np.linalg.LinAlgError:
        return False


def nearestSPD(A):
    """
    Nearest Symmetric Positive Definite matrix to A.

    The Frobenius norm is used: the rms difference of the elements.

    Parameters
    ----------
    A : ndarray, 2-D
        Matrix that should be SPD, but might not be, perhaps because
        of floating point arithmetic, or limitations in the data
        available to estimate its elements.

    Returns
    -------
    Ahat : ndarray, 2-D
        (Almost) nearest positive definite matrix to A

    Notes
    -----
    From Higham: "The nearest symmetric positive semidefinite matrix in the
    Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
    where H is the symmetric polar factor of B=(A + A')/2."
    http://www.sciencedirect.com/science/article/pii/0024379588902236

    Code and docstring are based on the Matlab m-file by  John D'Errico:
    http://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

    The tweaking method was changed to be more robust in the common
    pathological case, where more than one eigenvalue is near zero, and
    when one eigenvalue is exactly zero.  To be conservative, we are
    always nudging the diagonal to larger values.
    """

    # Ensure symmetry:
    B = (A + A.T) / 2

    # Symmetric polar factor, H: (in numpy svd, B = U S V, not U S V').
    U, S, V = np.linalg.svd(B)
    H = np.dot(V.T * S, V)

    Ahat = (B + H) / 2

    Ahat = (Ahat + Ahat.T) / 2

    # At this point, given floating point errors and differences
    # among algorithms, Ahat might be on the PD boundary, or too
    # close to it for some numerical operations. Adjust it:

    n = A.shape[0]
    k = 0
    # The "k == 0" is included to avoid a warning from
    # np.random.multivariate_normal, which seems to have a PD-detection
    # algorithm that occasionally fails on matrices that pass the
    # cholesky test.
    while k == 0 or not _is_PD(Ahat):
        k += 1
        # Tweaking strategy differs from D'Errico version.  It
        # is still a very small adjustment, but much larger than
        # his.
        # Eigvals are or can be complex dtype, so take abs().
        maxeig = np.abs(np.linalg.eigvals(Ahat)).max()
        Ahat[np.diag_indices(n)] += np.spacing(maxeig)
        # Normally no more than one adjustment will be needed.
        if k > 100:
            warnings.warn('adjustment in nearestSPD did not converge; '
                          'returning diagonal')
            return np.diag(np.diag(A))
    return Ahat


def _confidence(coef, cnstit, opt, t, e, tin, elor, xraw, xmod, W, m, B,
                Xu, Yu, Xv, Yv):
    """
    This confidence interval calculation does not correspond
    to a single ut_ matlab function, but is based on code
    in ut_solv1 starting from the line

    '%% spectral power (Puu, Pvv, Puv) of residual'.

    It returns its first argument, the coef dictionary,
    with additional entries.
    """

    # Confidence Intervals.
    if not opt['white']:
        Puu, Pvv, Puv = band_averaged_psd_by_constit(tin, t, e, elor, coef, opt)  # noqa

    # Make temporaries for quantities needed more than once.
    _Wx = W * xraw
    _WB = W[:, np.newaxis] * B

    # Matlab is recalculating xmod here; but we already have it.
    # In the 1-D case xmod is only the real part, but in that
    # case _Wx is real, and we are taking the real part in the
    # end, so the imaginary part would not contribute anything.
    nt = len(xraw)
    nm = B.shape[1]

    # Mean Square Misfit: Eq. 52 (mean squared error).
    varMSM = np.real(np.dot(xraw.conj(), _Wx) -
                     np.dot(xmod.conj(), _Wx)) / (nt-nm)

    # Gamma_C: covariance Eq. 54.
    gamC = np.linalg.inv(np.dot(B.conj().T, _WB)) * varMSM

    # Gamma_P: pseudo-covariance Eq. 54.
    gamP = np.linalg.inv(np.dot(B.T, _WB))
    gamP *= (np.dot(xraw, _Wx) - np.dot(xmod, _Wx)) / (nt - nm)

    del _Wx, _WB

    # Eq. 55; convenient intermediate variables; see Eq. 51.
    Gall = gamC + gamP
    Hall = gamC - gamP

    nc = len(Xu)
    coef.g_ci = np.nan * np.ones_like(coef.g)
    if opt['twodim']:
        coef.Lsmaj_ci = coef.g_ci.copy()
        coef.Lsmin_ci = coef.g_ci.copy()
        coef.theta_ci = coef.g_ci.copy()
        varcov_mCw = np.nan * np.ones((nc, 4, 4))
    else:
        coef.A_ci = coef.g_ci.copy()
        varcov_mCw = np.nan * np.ones((nc, 2, 2))

    if not opt['white']:
        varcov_mCc = np.copy(varcov_mCw)

    for c in range(nc):
        G = np.array([[Gall[c, c], Gall[c, c+nc]],
                      [Gall[c + nc, c], Gall[c+nc, c+nc]]])
        H = np.array([[Hall[c, c], Hall[c, c+nc]],
                      [Hall[c + nc, c], Hall[c+nc, c+nc]]])
        varXu = np.real(G[0, 0] + G[1, 1] + 2 * G[0, 1]) / 2
        varYu = np.real(H[0, 0] + H[1, 1] - 2 * H[0, 1]) / 2

        if opt['twodim']:
            varXv = np.real(H[0, 0] + H[1, 1] + 2 * H[0, 1]) / 2
            varYv = np.real(G[0, 0] + G[1, 1] - 2 * G[0, 1]) / 2

        if opt['linci']:  # Linearized.
            if not opt['twodim']:
                varcov_mCw[c, :, :] = np.diag(np.array([varXu, varYu]))
                if not opt['white']:
                    den = varXu + varYu
                    varXu = Puu[c] * varXu / den
                    varYu = Puu[c] * varYu / den
                    varcov_mCc[c, :, :] = np.diag(np.array([varXu, varYu]))
                sig1, sig2 = ut_linci(Xu[c], Yu[c], np.sqrt(varXu),
                                      np.sqrt(varYu))
                coef['A_ci'][c] = 1.96 * sig1
                coef['g_ci'][c] = 1.96 * sig2
            else:
                varcov_mCw[c, :, :] = np.diag(np.array([varXu, varYu,
                                                        varXv, varYv]))
                if not opt['white']:
                    den = varXv + varYv
                    varXv = Pvv[c] * varXv / den
                    varYv = Pvv[c] * varYv / den
                    varcov_mCc[c, :, :] = np.diag(np.array([varXu, varYu,
                                                            varXv, varYv]))
                with np.errstate(invalid='ignore'):
                    sig1, sig2 = ut_linci(Xu[c] + 1j * Xv[c],
                                      Yu[c] + 1j * Yv[c],
                                      np.sqrt(varXu) + 1j * np.sqrt(varXv),
                                      np.sqrt(varYu) + 1j * np.sqrt(varYv))
                coef['Lsmaj_ci'][c] = 1.96 * np.real(sig1)
                coef['Lsmin_ci'][c] = 1.96 * np.imag(sig1)
                coef['g_ci'][c] = 1.96 * np.real(sig2)
                coef['theta_ci'][c] = 1.96 * np.imag(sig2)

        else:  # Monte Carlo.
            covXuYu = np.imag(H[0, 0] - H[0, 1] + H[1, 0] - H[1, 1]) / 2
            Duu = np.array([[varXu, covXuYu], [covXuYu, varYu]])
            varcov_mCw[c, :2, :2] = Duu

            if not opt.white:
                Duu = Puu[c] * Duu / np.trace(Duu)
                varcov_mCc[c, :2, :2] = Duu

            if not opt.twodim:
                if not opt.white:
                    varcov_mCc[c, :, :] = nearestSPD(varcov_mCc[c, :, :])
                    mCall = np.random.multivariate_normal((Xu[c], Yu[c]),
                                                          varcov_mCc[c],
                                                          opt.nrlzn)
                else:
                    mCall = np.random.multivariate_normal((Xu[c], Yu[c]),
                                                          varcov_mCw[c],
                                                          opt.nrlzn)
                A, _, _, g = ut_cs2cep(mCall)
                coef.A_ci[c] = 1.96 * np.median(np.abs(A - np.median(A))) / 0.6745  # noqa
                g[0] = coef.g[c]
                g = cluster(g, 360)
                coef.g_ci[c] = 1.96 * np.median(np.abs(g - np.median(g))) / 0.6745  # noqa
            else:
                covXvYv = np.imag(G[0, 0] - G[0, 1] + G[1, 0] - G[1, 1]) / 2
                Dvv = np.array([[varXv, covXvYv], [covXvYv, varYv]])
                varcov_mCw[c, 2:, 2:] = Dvv

                if not opt.white:
                    Dvv = Pvv[c] * Dvv / np.trace(Dvv)
                    varcov_mCc[c, 2:, 2:] = Dvv
                covXuXv = np.imag(- H[0, 0] - H[0, 1] - H[1, 0] - H[1, 1]) / 2
                covXuYv = np.real(G[0, 0] - G[1, 1]) / 2
                covYuXv = np.real(- H[0, 0] + H[1, 1]) / 2
                covYuYv = np.imag(- G[0, 0] + G[0, 1] + G[1, 0] - G[1, 1]) / 2
                Duv = np.array([[covXuXv, covXuYv], [covYuXv, covYuYv]])
                varcov_mCw[c, :2, 2:] = Duv
                varcov_mCw[c, 2:, :2] = Duv.T

                if not opt.white:
                    if np.abs(Duv).sum() > 0:
                        Duv = Puv[c] * Duv / np.abs(Duv).sum()
                        varcov_mCc[c, :2, 2:] = Duv
                        varcov_mCc[c, 2:, :2] = Duv.T
                    else:
                        varcov_mCc[c, :2, 2:] = 0
                        varcov_mCc[c, 2:, :2] = 0

                    varcov_mCc[c] = nearestSPD(varcov_mCc[c])
                    mCall = np.random.multivariate_normal((Xu[c], Yu[c], Xv[c], Yv[c]), varcov_mCc[c], opt.nrlzn)  # noqa
                else:
                    mCall = np.random.multivariate_normal((Xu[c], Yu[c], Xv[c], Yv[c]), varcov_mCw[c], opt.nrlzn)  # noqa
                Lsmaj, Lsmin, theta, g = ut_cs2cep(mCall)
                coef.Lsmaj_ci[c] = 1.96 * np.median(np.abs(Lsmaj - np.median(Lsmaj))) / 0.6745  # noqa
                coef.Lsmin_ci[c] = 1.96 * np.median(np.abs(Lsmin - np.median(Lsmin))) / 0.6745  # noqa
                theta[0] = coef.theta[c]
                theta = cluster(theta, 360)
                coef.theta_ci[c] = 1.96 * np.median(np.abs(theta - np.median(theta))) / 0.6745  # noqa
                g[0] = coef.g[c]
                g = cluster(g, 360)
                coef.g_ci[c] = 1.96 * np.median(np.abs(g - np.median(g))) / 0.6745  # noqa

    nNR = coef.nNR

    if opt.infer:
        ind = nc
        if opt.linci:
            for k, ref in enumerate(cnstit.R):
                varcov = varcov_mCw if opt.white else varcov_mCc
                varReap = 0.25 * varcov[nNR + k, 0, 0]
                varImap = 0.25 * varcov[nNR + k, 1, 1]
                if opt.twodim:
                    varReap += 0.25 * varcov[nNR + k, 3, 3]
                    varImap += 0.25 * varcov[nNR + k, 2, 2]
                rp = ref.I.Rp
                rm = ref.I.Rm
                varXuHH = ((rp.real**2 + rm.real**2) * varReap +
                           (rp.imag**2 + rm.imag**2) * varImap)
                varYuHH = ((rp.real**2 + rm.real**2) * varImap +
                           (rp.imag**2 + rm.imag**2) * varReap)
                for varX, varY in zip(varXuHH, varYuHH):
                    if not opt.twodim:
                        sig1, sig2 = ut_linci(Xu[nNR + k], Yu[nNR + k],
                                              np.sqrt(varX), np.sqrt(varY))
                        coef.A_ci[ind] = 1.96 * sig1
                        coef.g_ci[ind] = 1.96 * sig2
                    else:
                        sig1, sig2 = ut_linci(Xu[nNR + k] + 1j * Xv[nNR + k],
                                              Yu[nNR + k] + 1j * Yv[nNR + k],
                                              np.sqrt(varX) + 1j *
                                              np.sqrt(varY),
                                              np.sqrt(varY) + 1j *
                                              np.sqrt(varX))
                        coef.Lsmaj_ci[ind] = 1.96 * sig1.real
                        coef.Lsmin_ci[ind] = 1.96 * sig1.imag
                        coef.g_ci[ind] = 1.96 * sig2.real
                        coef.theta_ci[ind] = 1.96 * sig2.imag
                    ind += 1
        else:
            raise NotImplementedError('Monte Carlo inference not implemented')

    return coef


def ut_linci(X, Y, sigX, sigY):
    # UT_LINCI()
    # current ellipse parameter uncertainties from cosine/sine coefficient
    # uncertainties, by linearized relations w/ correlations presumed zero
    # inputs: (two-dim case complex, one-dim case real)
    # X = Xu + i*Xv
    # Y = Yu + i*Yv
    # for Xu =real(X) = u cosine coeff; Yu =real(Y) = u sine coeff
    # Xv =imag(X) = v cosine coeff; Yv =imag(Y) = v sine coeff
    # sigX = sigXu + i*sigXv
    # sigY = sigYu + i*sigYv
    # for sigXu =real(sigX) =stddev(Xu); sigYu =real(sigY) =stddev(Yu)
    # sigXv =imag(sigX) =stddev(Xv); sigYv =imag(sigY) =stddev(Yv)
    # outputs:
    # two-dim case, complex
    # sig1 = sig_Lsmaj +1i*sig_Lsmin [same units as inputs]
    # sig2 = sig_g + 1i*sig_theta [degrees]
    # one-dim case, real
    # sig1 = sig_A [same units as inputs]
    # sig2 = sig_g [degrees]
    # UTide v1p0 9/2011 d.codiga@gso.uri.edu
    # (adapted from errell.m of t_tide, Pawlowicz et al 2002)

    X = np.array([X])
    Y = np.array([Y])
    sigX = np.array([sigX])
    sigY = np.array([sigY])
    Xu = np.real(X[:])
    sigXu = np.real(sigX)
    Yu = np.real(Y[:])
    sigYu = np.real(sigY)

    Xv = np.imag(X[:])
    sigXv = np.imag(sigX[:])
    Yv = np.imag(Y[:])
    sigYv = np.imag(sigY[:])

    rp = 0.5 * np.sqrt((Xu+Yv)**2 + (Xv-Yu)**2)
    rm = 0.5 * np.sqrt((Xu-Yv)**2 + (Xv+Yu)**2)
    sigXu2 = sigXu**2
    sigYu2 = sigYu**2
    sigXv2 = sigXv**2
    sigYv2 = sigYv**2

    ex = (Xu+Yv) / rp
    fx = (Xu-Yv) / rm
    gx = (Yu-Xv) / rp
    hx = (Yu+Xv) / rm

    # major axis
    dXu2 = (0.25*(ex+fx))**2
    dYu2 = (0.25*(gx+hx))**2
    dXv2 = (0.25*(hx-gx))**2
    dYv2 = (0.25*(ex-fx))**2
    sig1 = np.sqrt(dXu2 * sigXu2 + dYu2 * sigYu2 +
                   dXv2 * sigXv2 + dYv2 * sigYv2)

    # phase
    rn = 2 * (Xu * Yu + Xv * Yv)
    rd = Xu**2 - Yu**2 + Xv**2 - Yv**2
    den = rn**2 + rd**2
    dXu2 = ((rd*Yu - rn*Xu) / den)**2
    dYu2 = ((rd*Xu + rn*Yu) / den)**2
    dXv2 = ((rd*Yv - rn*Xv) / den)**2
    dYv2 = ((rd*Xv + rn*Yv) / den)**2
    sig2 = (180/np.pi) * np.sqrt(dXu2 * sigXu2 + dYu2 * sigYu2 +
                                 dXv2 * sigXv2 + dYv2 * sigYv2)

    # if ~isreal(X)
    if not np.isreal(X):
        # Minor axis.
        dXu2 = (0.25 * (ex-fx))**2
        dYu2 = (0.25 * (gx-hx))**2
        dXv2 = (0.25 * (hx+gx))**2
        dYv2 = (0.25 * (ex+fx))**2
        sig1 = sig1 + 1j*np.sqrt(dXu2 * sigXu2 + dYu2 * sigYu2 +
                                 dXv2 * sigXv2 + dYv2 * sigYv2)

        # Orientation.
        rn = 2.0 * (Xu * Xv + Yu * Yv)
        rd = Xu**2 + Yu**2 - (Xv**2 + Yv**2)
        den = rn**2 + rd**2
        dXu2 = ((rd*Xv - rn*Xu) / den)**2
        dYu2 = ((rd*Yv - rn*Yu) / den)**2
        dXv2 = ((rd*Xu + rn*Xv) / den)**2
        dYv2 = ((rd*Yu + rn*Yv) / den)**2
        sig2 = sig2 + 1j*(180/np.pi) * np.sqrt(dXu2 * sigXu2 + dYu2 * sigYu2 +
                                               dXv2 * sigXv2 + dYv2 * sigYv2)

    return sig1, sig2
