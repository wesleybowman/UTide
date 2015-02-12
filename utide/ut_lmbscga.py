# FIXME: Incomplete.
import numpy as np

def ut_lmbscga(x, t, w, ofac):

    x = x[:].T
    t = t[:].T
    w = w[:].T

    ofac = round(ofac)
    n = len(x)

    dt = (max(t) - min(t)) / (n-1)
    w = interp1(min(t), max(t), w, t)
    xw = x * w

    return Pxx, F
