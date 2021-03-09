import numpy as np


def _PE(coef):
    """
    Return the energy percentage for each constituent.
    """
    if "Lsmaj" in coef:
        E = coef["Lsmaj"] ** 2 + coef["Lsmin"] ** 2
        PE = (100 / np.sum(E)) * E
    else:
        PE = 100 * coef["A"] ** 2 / np.sum(coef["A"] ** 2)
    return PE


def _SNR(coef):
    """
    Return the signal-to-noise ratio for each constituent.
    """
    if "Lsmaj" in coef:
        SNR = (coef["Lsmaj"] ** 2 + coef["Lsmin"] ** 2) / (
            (coef["Lsmaj_ci"] / 1.96) ** 2 + (coef["Lsmin_ci"] / 1.96) ** 2
        )
    else:
        SNR = (coef["A"] ** 2) / (coef["A_ci"] / 1.96) ** 2
    return SNR


def ut_diagn(coef):
    """
    Add to coef the names, PE, and SNR, *always* sorted by energy.

    To be eliminated...
    """
    coef["diagn"] = {}
    PE = _PE(coef)
    SNR = _SNR(coef)
    indPE = PE.argsort()[::-1]

    coef["diagn"]["name"] = coef["name"][indPE]
    coef["diagn"]["PE"] = PE[indPE]
    coef["diagn"]["SNR"] = SNR[indPE]

    return coef
