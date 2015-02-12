from __future__ import division
import numpy as np


def robustfit(X, y, wfun='bisquare', tune=[], const=True,
              priorw=[], dowarn=True):

    tuning = {'andrews': 1.339,
              'bisquare': 4.685,
              'cauchy': 2.385,
              'fair': 1.400,
              'huber': 1.345,
              'logistic': 1.205,
              'ols': 1,
              'talwar': 2.795,
              'welsch': 2.985}

    tune = tuning[wfun]

    if not priorw:
        priorw = np.ones(y.shape)

    # statrobustfit(X, y, wfun, tune, wasnan, const, priorw, dowarn)
    b, stats = statrobustfit(X, y, wfun, tune, const, priorw, dowarn)


# def statrobustfit(X, y, wfun, tune, wasnan, const, priorw, dowarn):
def statrobustfit(X, y, wfun, tune, const, priorw, dowarn):

    n, p = X.shape

    if const:
        X = np.hstack((np.ones((n, 1)), X))
        p += 1

    if not np.all(priorw == 1):
        sw = np.sqrt(priorw)
        X = X * sw
        y = y * sw
    else:
        sw = 1

    Q, R = np.linalg.qr(X)
    tol = abs(R[0]) * max(n, p) * np.spacing(1)
    xrank = sum(abs(np.diag(R)) > tol)

    if xrank == p:
        b, resid, rank, s = np.linalg.lstsq(R, np.dot(Q.T, y))

    # FIXME: 'b0' is assigned to but never used.
    b0 = np.zeros(b.shape)

    stats = {}

    return b, stats


def andrews(r):
    r = max(np.sqrt(np.spacing(1)), abs(r))
    w = (abs(r) < np.pi) * np.sin(r) / r
    return w


def bisquare(r):
    w = (abs(r) < 1) * (1 - r**2)**2
    return w


def cauchy(r):
    w = 1 / (1 + r**2)
    return w


def fair(r):
    w = 1 / (1 + abs(r))
    return w


def huber(r):
    w = 1 / max(1, abs(r))
    return w


def logistic(r):
    r = max(np.sqrt(np.single(1)), abs(r))
    w = np.tanh(r) / r
    return w


def ols(r):
    w = np.ones(len(r))
    return w


def talwar(r):
    w = 1 * (abs(r) < 1)
    return w


def welsch(r):
    w = np.exp(-(r**2))
    return w


def test(r):
    print(andrews(r))
    print(bisquare(r))
    print(cauchy(r))
    print(fair(r))
    print(huber(r))
    print(logistic(r))
    print(talwar(r))
    print(welsch(r))

#        case 'andrews'
#            wfun = @andrews;
#            t = 1.339;
#        case 'bisquare'
#            wfun = @bisquare;
#            t = 4.685;
#        case 'cauchy'
#            wfun = @cauchy;
#            t= 2.385;
#        case 'fair'
#            wfun = @fair;
#            t = 1.400;
#        case 'huber'
#            wfun = @huber;
#            t = 1.345;
#        case 'logistic'
#            wfun = @logistic;
#            t = 1.205;
#        case 'ols'
#            wfun = @ols;
#            t = 1;
#        case 'talwar'
#            wfun = @talwar;
#            t = 2.795;
#        case 'welsch'
#            wfun = @welsch;
#            t = 2.985
