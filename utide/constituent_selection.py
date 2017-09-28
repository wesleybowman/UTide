from __future__ import (absolute_import, division, print_function)

from collections import OrderedDict

import numpy as np

from .astronomy import ut_astron
from ._ut_constants import ut_constants, constit_index_dict
from .utilities import Bunch


def ut_cnstitsel(tref, minres, incnstit, infer):
    """
    Select constituents and organize constituent data.

    inputs
      tref = reference time (UTC, days relative to Python datetime epoch)
      minres = freq separation (cph) used in decision tree
      incnstit = 'cnstit' input to ut_solv
      infer = 'opt.infer' input to ut_solv
    outputs
      cnstit.NR.name = list of 4-char names of NR constits
      cnstit.NR.frq = frequencies (cph) of NR constits
      cnstit.NR.lind = list indices (in ut_constants.mat) of NR constits
      cnstit.R = empty if no inference; otherwise, for each (i'th) R constit:
          cnstit.R[i].name, .frq, .lind = as above, but for R constits
          cnstit.R[i].I[j].name, .frq, .lind = as above for j'th I constit
      coef.nNR, coef.nR, coef.nI = number non-reference, reference,
          inferred constituents
      coef.name = list of names of all constituents (NR, R, and I)
      coef.aux.frq = frequencies (cph) of all constituents
      coef.aux.lind = list indices of all constituents
      coef.aux.reftime = tref
    """

    shallow = ut_constants.shallow
    const = ut_constants.const

    cnstit = Bunch()
    coef = Bunch()

    astro, ader = ut_astron(tref)

    ii = np.isfinite(const.ishallow)
    const.freq[~ii] = np.dot(const.doodson[~ii, :], ader[:, 0]) / 24

    for k in ii.nonzero()[0]:
        ik = const.ishallow[k] + np.arange(const.nshallow[k])
        ik = ik.astype(int) - 1
        const.freq[k] = np.sum(const.freq[shallow.iname[ik] - 1] *
                               shallow.coef[ik])

    # cnstit.NR
    cnstit['NR'] = Bunch()

    # if incnstit.lower() == 'auto':
    if incnstit == 'auto':
        cnstit.NR.lind = np.where(const.df >= minres)[0]
    else:
        cnstit.NR.lind = [constit_index_dict[n] for n in incnstit]

    # Remove from NR any R and I constituents.
    if infer is not None:
        RIset = set(infer.inferred_names) | set(infer.reference_names)
        RI_index_set = {constit_index_dict[n] for n in RIset}
        cnstit.NR.lind = [ind for ind in cnstit.NR.lind
                          if ind not in RI_index_set]

    cnstit.NR.frq = const.freq[cnstit.NR.lind]
    cnstit.NR.name = const.name[cnstit.NR.lind]
    nNR = len(cnstit.NR.frq)

    # cnstit.R
    nR = 0
    nI = 0
    cnstit.R = []

    if infer is not None:
        nI = len(infer.inferred_names)
        # Find unique reference names
        _r = infer.reference_names
        allrefs = list(OrderedDict(zip(_r, [1]*len(_r))).keys())
        nR = len(allrefs)
        for k, name in enumerate(allrefs):
            refstruct = Bunch(name=name)
            refstruct.lind = constit_index_dict[name]
            refstruct.frq = const.freq[refstruct.lind]
            ind = [i for i, rname in enumerate(infer.reference_names)
                   if name == rname]
            refstruct.nI = len(ind)
            refstruct.I = Bunch(Rp=[], Rm=[], name=[], lind=[], frq=[])
            for lk, ilk in enumerate(ind):
                refstruct.I.Rp.append(infer.amp_ratios[ilk] *
                                      np.exp(1j * infer.phase_offsets[ilk] *
                                      np.pi/180))
                if len(infer.amp_ratios) > nI:
                    refstruct.I.Rm.append(infer.amp_ratios[ilk + nI] *
                                          np.exp(-1j *
                                          infer.phase_offsets[ilk + nI] *
                                          np.pi / 180))
                else:
                    refstruct.I.Rm.append(np.conj(refstruct.I.Rp[lk]))

                iname = infer.inferred_names[ilk]
                refstruct.I.name.append(iname)
                lind = constit_index_dict[iname]
                refstruct.I.lind.append(lind)
                refstruct.I.frq.append(const.freq[lind])

            refstruct.I.Rp = np.array(refstruct.I.Rp)
            refstruct.I.Rm = np.array(refstruct.I.Rm)
            cnstit.R.append(refstruct)

    coef.name = list(cnstit.NR.name[:])
    coef.aux = Bunch(frq=list(cnstit.NR.frq[:]),
                     lind=list(cnstit.NR.lind[:]),
                     reftime=tref)

    if infer is not None:
        # Append reference values, and then inferred values, to the lists.
        coef.name.extend(allrefs)
        coef.aux.frq.extend([_ref.frq for _ref in cnstit.R])
        coef.aux.lind.extend([_ref.lind for _ref in cnstit.R])
        for ref in cnstit.R:
            coef.name.extend(ref.I.name)
            coef.aux.frq.extend(ref.I.frq)
            coef.aux.lind.extend(ref.I.lind)

    coef.name = np.array(coef.name, dtype=object)
    coef.aux.frq = np.array(coef.aux.frq, dtype=float)
    coef.aux.lind = np.array(coef.aux.lind, dtype=int)

    coef.nR = nR
    coef.nNR = nNR
    coef.nI = nI

    return cnstit, coef
