import numpy as np


def ut_fbndavg(P, allfrq, cfrq):
    # UT_FBNDAVG()
    # line-decimate and band-average spectra
    # inputs
    # P = periodogram to treat [e units^2/cph]
    # allfrq = frequency values of (equispaced) P estimates [cph]
    # cfrq = frequencies of constituents [cph] (nc x 1)
    # outputs
    # avP = line-decimated and band-averaged spectrum [e units^2/cph] (9 x 1)
    # fbnd = frequencies [cph] at edges of averaged bands (9 x 2)
    # UTide v1p0 9/2011 d.codiga@gso.uri.edu
    # (based on residual_spectrum.m of t_tide, Pawlowicz et al 2002)

    df = allfrq[2] - allfrq[1]
    # P[np.round(cfrq/df).astype(int)+1] = np.nan
    np.round(cfrq/df).astype(int)
    P[np.round(cfrq/df).astype(int)] = np.nan

    fbnd = np.array([[.00010, .00417],
                     [.03192, .04859],
                     [.07218, .08884],
                     [.11243, .12910],
                     [.15269, .16936],
                     [.19295, .20961],
                     [.23320, .25100],
                     [.26000, .29000],
                     [.30000, .50000]])

    # nfbnd=size(fbnd, 1)
    nfbnd = fbnd.shape[0]
    avP = np.zeros((nfbnd, 1))

    # import pdb; pdb.set_trace()

    # for k=nfbnd:-1:1,
    for k in np.arange(nfbnd-1, -1, -1):
        b1 = np.where(allfrq >= fbnd[k, 0])[0]
        b2 = np.where(allfrq <= fbnd[k, 1])[0]
        b3 = np.where(np.isfinite(P))[0]
        jbnd = np.intersect1d(np.intersect1d(b1, b2), b3)
        # Issue is here, taking the mean of Pvv1s when Pvv1s is already off.
        if jbnd.any():
            # avP[k]=np.mean(P[jbnd-1])
            avP[k] = np.mean(P[jbnd])
        elif k < nfbnd:
            avP[k] = P[k+1]

    return avP, fbnd
