import numpy as np
import scipy.io as sio
from ut_astron import ut_astron

def ut_cnstitsel(tref,minres,incnstit,infer):

    #mat_contents = sio.loadmat('./ut_constants.mat', struct_as_record=False, squeeze_me=True)
    mat_contents = sio.loadmat(ut_constants, struct_as_record=False, squeeze_me=True)
    shallow = mat_contents['shallow']
    const = mat_contents['const']

    cnstit = {}
    coef = {}

    astro, ader = ut_astron(tref)

    ii = np.isfinite(const.ishallow)
    const.freq[~ii] = np.dot(const.doodson[~ii, :], ader) / 24

    for k in ii.nonzero()[0]:
        ik = const.ishallow[k]+np.arange(const.nshallow[k])
        ik = ik.astype(int)-1
        const.freq[k] = np.sum(const.freq[shallow.iname[ik]]*shallow.coef[ik])

    ## cnstit.NR
    cnstit['NR'] = {}
    if incnstit.lower()=='auto':
        cnstit['NR']['lind'] = np.where(const.df >= minres)[0]
    else:
        pass

    # skipped some stuff here cause they involve infer

    cnstit['NR']['frq'] = const.freq[cnstit['NR']['lind']]
    cnstit['NR']['name'] = const.name[cnstit['NR']['lind']]
    nNR = len(cnstit['NR']['frq'])

    ## cnstit.R
    nR = 0
    nI = 0
    cnstit['R'] = []

    nallc = nNR+nR+nI

    coef['name'] = cnstit['NR']['name']
    coef['aux'] = {}
    coef['aux']['frq'] = cnstit['NR']['frq']
    coef['aux']['lind'] = cnstit['NR']['lind']

    # another infer if statement

    coef['aux']['reftime'] = tref


    return nNR, nR, nI, cnstit, coef
