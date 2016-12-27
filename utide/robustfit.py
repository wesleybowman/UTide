"""
Robust MLR via iteratively reweighted least squares.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np

from utide.utilities import Bunch


# Weighting functions:
def andrews(r):
    r = np.abs(r)
    r = max(np.sqrt(np.spacing(1)), r)
    w = (r < np.pi) * np.sin(r) / r
    return w


def bisquare(r):
    r = np.abs(r)
    w = (r < 1) * (1 - r**2)**2
    return w


def cauchy(r):
    r = np.abs(r)
    w = 1 / (1 + r**2)
    return w


def fair(r):
    w = 1 / (1 + np.abs(r))
    return w


def huber(r):
    w = 1 / max(1, np.abs(r))
    return w


def logistic(r):
    r = np.abs(r)
    r = max(np.sqrt(np.single(1)), r)
    w = np.tanh(r) / r
    return w


def ols(r):
    w = np.ones(len(r))
    return w


def talwar(r):
    w = (np.abs(r) < 1).astype(float)
    return w


def welsch(r):
    r = np.abs(r)
    w = np.exp(-(r**2))
    return w


wfuncdict = dict(andrews=andrews,
                 bisquare=bisquare,
                 cauchy=cauchy,
                 fair=fair,
                 huber=huber,
                 logistic=logistic,
                 ols=ols,
                 talwar=talwar,
                 welsch=welsch)

tune_defaults = {'andrews': 1.339,
                 'bisquare': 4.685,
                 'cauchy': 2.385,
                 'fair': 1.400,
                 'huber': 1.345,
                 'logistic': 1.205,
                 'ols': 1,
                 'talwar': 2.795,
                 'welsch': 2.985}


def sigma_hat(x):
    """
    Robust estimate of standard deviation based on medians.
    """
    # The center could be based on the mean or some other function.
    return np.median(np.abs(x - np.median(x))) / 0.6745


def leverage(x):
    """
    Calculate leverage as the diagonal of the "Hat" matrix of the
    model matrix, x.
    """

    # The Hat is x times its pseudo-inverse.
    # In einum, the diagonal is calculated for each row of x
    # and column of pinv as the dot product of column j of x.T
    # and column j of pinv; hence the 'j' in the output means
    # *don't* sum over j.

    hdiag = np.einsum('ij, ij -> j', x.T, np.linalg.pinv(x))
    # This should be real and positive, but with floating point
    # arithmetic the imaginary part is not exactly zero.
    return np.abs(hdiag)


def r_normed(R, rfac):
    """
    Normalized residuals from raw residuals and a multiplicative factor.
    """
    return rfac * R / sigma_hat(R)


def robustfit(X, y, weight_function='bisquare', tune=None,
              rcond=1, tol=0.001, maxit=50):
    """
    Multiple linear regression via iteratively reweighted least squares.

    Parameters
    ----------
    X : ndarray (n, p)
        MLR model with `p` parameters (independent variables) at `n` times
    y : ndarray (n,)
        dependent variable
    weight_function : string, optional
        name of weighting function
    tune : None or float, optional
        Tuning parameter for normalizing residuals in weight calculation;
        larger numbers *decrease* the sensitivity to outliers. If `None`,
        a default will be provided based on the `weight_function`.
    rcond : float, optional
        minimum condition number parameter for `np.linalg.lstsq`
    tol : float, optional
        When the fractional reduction in mean squared weighted residuals
        is less than `tol`, the iteration stops.
    maxit : integer, optional
        Maximum number of iterations.

    Returns
    -------
    rf : `utide.utilities.Bunch`

        - rf.b: model coefficients of the solution
        - rf.w: weights used for the solution
        - rf.s: singular values for each model component
        - rf.rms_resid: rms residuals (unweighted) from the fit
        - rf.leverage: sensitivity of the OLS estimate to each point in `y`
        - rf.ols_b: OLS model coefficients
        - rf.ols_rms_resid: rms residuals from the OLS fit
        - rf.iterations: number of iterations completed

    """

    if tune is None:
        tune = tune_defaults[weight_function]

    _wfunc = wfuncdict[weight_function]

    if X.ndim == 1:
        X = X.reshape((x.size, 1))
    n, p = X.shape

    lev = leverage(X)

    out = Bunch(weight_function=weight_function,
                tune=tune,
                rcond=rcond,
                tol=tol,
                maxit=maxit,
                leverage=lev)

    # LJ2009 has an incorrect expression for leverage in the
    # appendix, and an incorrect version of the following
    # multiplicative factor for scaling the residuals.

    rfac = 1 / (tune * np.sqrt(1 - lev))

    # We probably only need to keep track of the rmeansq, but
    # it's cheap to carry along rsumsq until we are positive.
    oldrsumsq = None
    oldrmeansq = None
    oldlstsq = None
    oldw = None
    iterations = 0  # 1-based iteration exit number
    w = np.ones(y.shape)

    for i in range(maxit):
        wX = w[:, np.newaxis] * X
        wy = w * y
        b, rsumsq, rank, sing = np.linalg.lstsq(wX, wy, rcond)
        rsumsq = rsumsq[0]
        if i == 0:
            rms_resid = np.sqrt(rsumsq / n)
            out.update(dict(ols_b=b,
                            ols_rms_resid=rms_resid))

        # Weighted mean of squared weighted residuals:
        rmeansq = rsumsq / w.sum()

        if oldrsumsq is not None:
            # improvement = (oldrsumsq - rsumsq) / oldrsumsq
            improvement = (oldrmeansq - rmeansq) / oldrmeansq
            # print("improvement:", improvement)

            if improvement < 0:
                b, rsumsq, rank, sing = oldlstsq
                w = oldw
                iterations = i
                break

            if improvement < tol:
                iterations = i + 1
                break

        # Save these values in case the next iteration
        # makes things worse.
        oldlstsq = b, rsumsq, rank, sing
        oldw = w
        oldrsumsq = rsumsq
        oldrmeansq = rmeansq

        # Residuals (unweighted) from latest fit:
        resid = y - np.dot(X, b)

        # Update weights based on these residuals.
        w = _wfunc(r_normed(resid, rfac))

    if iterations == 0:
        iterations = maxit  # Did not converge.

    rms_resid = np.sqrt(np.mean(np.abs(resid)**2))

    out.update(dict(iterations=iterations,
                    b=b,
                    s=sing,
                    w=w,
                    rank=rank,
                    rms_resid=rms_resid,
                    ))

    return out


# Some simple test cases; this probably will be removed.

if __name__ == '__main__':

    np.random.seed(1)
    n = 10000
    x = np.arange(n)
    x0 = np.ones_like(x)
    x1 = np.exp(1j * x/9)
    x2 = np.exp(1j * x/7)
    y = (1 + 1j) * x1 + (2 - 1j) * x2 + (0.1 * np.random.randn(n) +
                                         0.1 * 1j * np.random.randn(n))
    y[::10] = (np.random.randn(n) + 1j * np.random.randn(n))[::10]

    y[10] = 3
    y[20] = 2 * 1j
    y[30] = -2 - 3 * 1j

    A = np.vstack((x0, x1, x2)).T
    c = np.linalg.lstsq(A, y)
    print('OLS:', c[0])

    rf1 = robustfit(A, y)

    print('robust:', rf1.b)
    print('another test: a very short real series')

    x = np.arange(1, 21, dtype=float)
    x0 = np.ones_like(x)
    xx = np.vstack((x0, x)).T

    # Signal for the model: linear trend.
    y = 2 * x

    # Some outliers.
    y[0] = 1.5
    y[2] = -2
    y[4] = 9.6

    # Use a sine as the "noise" component; not part of the model.
    y = y + 0.1 * np.sin(x)

    rf2 = robustfit(xx, y)
    print(np.linalg.lstsq(xx, y)[0])
    print(rf2.b)
