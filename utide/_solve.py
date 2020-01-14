"""
Central module for calculating the tidal amplitudes, phases, etc.
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np

from .harmonics import ut_E
from .diagnostics import ut_diagn
from .ellipse_params import ut_cs2cep
from .constituent_selection import ut_cnstitsel
from .confidence import _confidence
from .utilities import Bunch
from ._ut_constants import constit_index_dict
from .robustfit import robustfit
from ._time_conversion import _normalize_time

default_opts = dict(constit='auto',
                    conf_int='linear',
                    method='ols',
                    trend=True,
                    phase='Greenwich',
                    nodal=True,
                    infer=None,
                    MC_n=200,
                    Rayleigh_min=1,
                    robust_kw=dict(weight_function='cauchy'),
                    white=False,
                    verbose=True,
                    epoch='python',
                    )


def _process_opts(opts, is_2D):
    newopts = Bunch(default_opts)
    newopts.update_values(strict=True, **opts)
    # TODO: add more validations.
    newopts.infer = validate_infer(newopts.infer, is_2D)

    compat_opts = _translate_opts(newopts)

    return compat_opts


def _translate_opts(opts):
    # Temporary shim between new-style options and Matlab heritage.
    # Here or elsewhere, proper validation remains to be added.
    oldopts = Bunch()
    oldopts.cnstit = opts.constit
    oldopts.infer = opts.infer     # we will not use the matlab names, though

    oldopts.conf_int = True
    if opts.conf_int == 'linear':
        oldopts.linci = True
    elif opts.conf_int == 'MC':
        oldopts.linci = False
    elif opts.conf_int == 'none':
        oldopts.conf_int = False
        oldopts.nodiagn = 1
    else:
        raise ValueError("'conf_int' must be 'linear', 'MC', or 'none'")

    oldopts.notrend = not opts.trend
    oldopts['nodesatlint'] = False
    oldopts['nodesatnone'] = False
    oldopts['gwchlint'] = False
    oldopts['gwchnone'] = False
    if opts.nodal == 'linear_time':
        oldopts['nodsatlint'] = True
    elif not opts.nodal:
        oldopts['nodsatnone'] = True
    if opts.phase == 'linear_time':
        oldopts['gwchlint'] = True
    elif opts.phase == 'raw':
        oldopts['gwchnone'] = True
    # Otherwise it should be default, 'Greenwich.'
    oldopts.rmin = opts.Rayleigh_min
    oldopts.white = opts.white
    oldopts.newopts = opts  # So we can access new opts via the single "opt."
    oldopts['RunTimeDisp'] = opts.verbose
    oldopts.epoch = opts.epoch
    return oldopts


def validate_infer(infer, is_2D):
    if infer is None or infer == 'none':
        return None
    required_keys = {'inferred_names', 'reference_names', 'amp_ratios',
                     'phase_offsets'}
    keys = set(infer.keys())
    if keys < required_keys:
        raise ValueError("infer option must include %s" % required_keys)
    nI = len(infer.inferred_names)
    if len(infer.reference_names) != nI:
        raise ValueError("inferred_names must be same"
                         "  length as reference_names")
    nratios = 2 * nI if is_2D else nI
    if (len(infer.amp_ratios) != nratios or
            len(infer.phase_offsets) != nratios):
        raise ValueError("ratios and offsets need to have length %d" %
                         nratios)
    if 'approximate' not in infer:
        infer.approximate = False
    return infer


def solve(t, u, v=None, lat=None, **opts):
    """
    Calculate amplitude, phase, confidence intervals of tidal constituents.

    Parameters
    ----------
    t : array_like
        Time in days since `epoch`.
    u : array_like
        Sea-surface height, velocity component, etc.
    v : {None, array_like}, optional
        If `u` is a velocity component, `v` is the orthogonal component.
    lat : float, required
        Latitude in degrees.
    epoch : {string, `datetime.date`, `datetime.datetime`}, optional
        Valid strings are 'python' (default); 'matlab' if `t` is
        an array of Matlab datenums; or an arbitrary date in the
        form 'YYYY-MM-DD'.  The default corresponds to the Python
        standard library `datetime` proleptic Gregorian calendar,
        starting with 1 at 00:00 on January 1 of year 1; this is
        the 'datenum' used by Matplotlib.
    constit : {'auto', array_like}, optional
        List of strings with standard letter abbreviations of
        tidal constituents; or 'auto' to let the list be determined
        based on the time span.
    conf_int : {'linear', 'MC', 'none'}, optional
        If not 'none' (string), calculate linearized confidence
        intervals, or use a Monte-Carlo simulation.
    method : {'ols', 'robust'}, optional
        Solve with ordinary least squares, or with a robust algorithm.
    trend : bool, optional
        True (default) to include a linear trend in the model.
    phase : {'Greenwich', 'linear_time', 'raw'}, optional
        Give Greenwich-referenced phase lags, an approximation
        using linearized times, or raw lags.
    nodal : {True, False, 'linear_time'}, optional
        True (default) to include nodal/satellite corrections;
        'linear_time' to use the linearized time approximation;
        False to omit nodal corrections.

    Returns
    -------
    coef : Bunch
        Data container with all configuration and solution information:

    Other Parameters
    ----------------
    infer : {None, dict or Bunch}, optional; default is None.
        If not None, the items are:

        **inferred_names** : {sequence of N strings}
            inferred constituent names
        **reference_names** : {sequence of N strings}
            reference constituent names
        **amp_ratios** : {sequence, N or 2N floats}
            amplitude ratios (unitless)
        **phase_offsets** : {sequence, N or 2N floats}
            phase offsets (degrees)
        **approximate** : {bool, optional (default is False)}
            use approximate method

        amp_ratios and phase_offsets have length N for a scalar
        time series, or 2N for a vector series.

    MC_n : integer, optional
        Not yet implemented.
    robust_kw : dict, optional
        Keyword arguments for `robustfit`, if `method` is 'robust'.
    Rayleigh_min : float
        Minimum conventional Rayleigh criterion for automatic
        constituent selection; default is 1.
    white : bool
        If False (default), use band-averaged spectra from the
        residuals in the confidence limit estimates; if True,
        assume a white background spectrum.
    verbose : {True, False}, optional
        True (default) turns on verbose output. False emits no messages.

    Note
    ----
    `utide.reconstruct` requires the calculation of confidence intervals.

    Notes
    -----

    To be added: much additional explanation.

    There will also be more "Other Parameters".

    """

    compat_opts = _process_opts(opts, v is not None)

    coef = _solv1(t, u, v, lat, **compat_opts)

    return coef


def _solv1(tin, uin, vin, lat, **opts):

    # The following returns a possibly modified copy of tin (ndarray).
    # t, u, v are fully edited ndarrays (unless v is None).
    packed = _slvinit(tin, uin, vin, lat, **opts)
    tin, t, u, v, tref, lor, elor, opt = packed
    nt = len(t)
    if opt['RunTimeDisp']:
        print('solve: ', end='')

    # opt['cnstit'] = cnstit
    cnstit, coef = ut_cnstitsel(tref, opt['rmin']/(24*lor),
                                opt['cnstit'], opt['infer'])

    # a function we don't need
    # coef.aux.rundescr = ut_rundescr(opt,nNR,nR,nI,t,tgd,uvgd,lat)

    coef.aux.opt = opt
    coef.aux.lat = lat

    if opt['RunTimeDisp']:
        print('matrix prep ... ', end='')

    ngflgs = [opt['nodsatlint'], opt['nodsatnone'],
              opt['gwchlint'], opt['gwchnone']]

    E_args = (lat, ngflgs, opt.prefilt)

    # Make the model array, starting with the harmonics.
    E = ut_E(t, tref, cnstit.NR.frq, cnstit.NR.lind, *E_args)

    # Positive and negative frequencies
    B = np.hstack((E, E.conj()))

    if opt.infer is not None:

        Etilp = np.empty((nt, coef.nR), dtype=complex)
        Etilm = np.empty((nt, coef.nR), dtype=complex)

        if not opt.infer.approximate:
            for k, ref in enumerate(cnstit.R):
                E = ut_E(t, tref, ref.frq, ref.lind, *E_args)
                # (nt,1)
                Q = ut_E(t, tref, ref.I.frq, ref.I.lind, *E_args) / E
                # (nt,ni)
                Qsum_p = (Q * ref.I.Rp).sum(axis=1)
                Etilp[:, k] = E[:, 0] * (1 + Qsum_p)
                Qsum_m = (Q * np.conj(ref.I.Rm)).sum(axis=1)
                Etilm[:, k] = E[:, 0] * (1 + Qsum_m)

        else:
            # Approximate inference.
            Q = np.empty((coef.nR,), dtype=float)
            beta = np.empty((coef.nR,), dtype=float)

            for k, ref in enumerate(cnstit.R):
                E = ut_E(t, tref, ref.frq, ref.lind, *E_args)[:, 0]
                Etilp[:, k] = E
                Etilm[:, k] = E
                num = ut_E(tref, tref, ref.I.frq, ref.I.lind, *E_args).real
                den = ut_E(tref, tref, ref.frq, ref.lind, *E_args).real
                Q[k] = ((num/den))[0, 0]
                arg = np.pi*lor*24*(ref.I.frq - ref.frq)*(nt+1) / nt
                beta[k] = np.sin(arg) / arg

        B = np.hstack((B, Etilp, np.conj(Etilm)))

    # add the mean
    B = np.hstack((B, np.ones((nt, 1))))

    if not opt['notrend']:
        B = np.hstack((B, ((t-tref)/lor)[:, np.newaxis]))

    # nm = B.shape[1]  # 2*(nNR + nR) + 1, plus 1 if trend is included.

    if opt['RunTimeDisp']:
        print('solution ... ', end='')

    if opt['twodim']:
        xraw = u + 1j*v
    else:
        xraw = u

    if opt.newopts.method == 'ols':
        # Model coefficients.
        try:
            m = np.linalg.lstsq(B, xraw, rcond=None)[0]
        except TypeError:
            m = np.linalg.lstsq(B, xraw)[0]
        W = np.ones(nt)  # Uniform weighting; we could use a scalar 1, or None.
    else:
        rf = robustfit(B, xraw, **opt.newopts.robust_kw)
        m = rf.b
        W = rf.w
        coef.rf = rf
    coef.weights = W

    xmod = np.dot(B, m)  # Model fit.

    if not opt['twodim']:
        xmod = np.real(xmod)

    e = W*(xraw-xmod)  # Weighted residuals.

    nI, nR, nNR = coef.nI, coef.nR, coef.nNR

    ap = np.hstack((m[:nNR], m[2*nNR:2*nNR+nR]))
    i0 = 2*nNR + nR
    am = np.hstack((m[nNR:2*nNR], m[i0:i0+nR]))

    Xu = np.real(ap + am)
    Yu = -np.imag(ap - am)

    if not opt['twodim']:
        coef['A'], _, _, coef['g'] = ut_cs2cep(Xu, Yu)
        Xv = []
        Yv = []

    else:
        Xv = np.imag(ap+am)
        Yv = np.real(ap-am)
        packed = ut_cs2cep(Xu, Yu, Xv, Yv)
        coef['Lsmaj'], coef['Lsmin'], coef['theta'], coef['g'] = packed

    # Mean and trend.
    if opt['twodim']:
        if opt['notrend']:
            coef['umean'] = np.real(m[-1])
            coef['vmean'] = np.imag(m[-1])
        else:
            coef['umean'] = np.real(m[-2])
            coef['vmean'] = np.imag(m[-2])
            coef['uslope'] = np.real(m[-1])/lor
            coef['vslope'] = np.imag(m[-1])/lor
    else:
        if opt['notrend']:
            coef['mean'] = np.real(m[-1])
        else:
            coef['mean'] = np.real(m[-2])
            coef['slope'] = np.real(m[-1])/lor

    if opt.infer:
        # complex coefficients
        apI = np.empty((nI,), dtype=complex)
        amI = np.empty((nI,), dtype=complex)
        ind = 0

        for k, ref in enumerate(cnstit.R):
            apI[ind:ind + ref.nI] = ref.I.Rp * ap[nNR + k]
            amI[ind:ind + ref.nI] = ref.I.Rm * am[nNR + k]
            ind += ref.nI

        XuI = (apI + amI).real
        YuI = -(apI - amI).imag

        if not opt.twodim:
            A, _, _, g = ut_cs2cep(XuI, YuI)
            coef.A = np.hstack((coef.A, A))
            coef.g = np.hstack((coef.g, g))
        else:
            XvI = (apI + amI).imag
            YvI = (apI - amI).real
            Lsmaj, Lsmin, theta, g = ut_cs2cep(XuI, YuI, XvI, YvI)
            coef.Lsmaj = np.hstack((coef.Lsmaj, Lsmaj))
            coef.Lsmin = np.hstack((coef.Lsmin, Lsmin))
            coef.theta = np.hstack((coef.theta, theta))
            coef.g = np.hstack((coef.g, g))

    if opt['conf_int'] is True:
        coef = _confidence(coef, cnstit, opt, t, e, tin, elor, xraw, xmod,
                           W, m, B, Xu, Yu, Xv, Yv)

    # Diagnostics.
    if not opt['nodiagn']:
        coef, indPE = ut_diagn(coef, opt)

    # Re-order constituents.
    if opt['ordercnstit'] is not None:

        if opt['ordercnstit'] == 'frq':
            ind = coef['aux']['frq'].argsort()

        elif opt['ordercnstit'] == 'snr':
            if not opt['nodiagn']:
                ind = coef['diagn']['SNR'].argsort()[::-1]
            else:
                if opt['twodim']:
                    SNR = (coef['Lsmaj']**2 + coef['Lsmin']**2) / (
                        (coef['Lsmaj_ci']/1.96)**2 +
                        (coef['Lsmin_ci']/1.96)**2)

                else:
                    SNR = (coef['A']**2) / (coef['A_ci']/1.96)**2

                ind = SNR.argsort()[::-1]

        else:
            ilist = [constit_index_dict[name] for name in opt['ordercnstit']]
            ind = np.array(ilist, dtype=int)

    else:  # Default: order by decreasing energy.
        if not opt['nodiagn']:
            ind = indPE
        else:
            if opt['twodim']:
                PE = np.sum(coef['Lsmaj']**2 + coef['Lsmin']**2)
                PE = 100 * (coef['Lsmaj']**2 + coef['Lsmin']**2) / PE
            else:
                PE = 100 * coef['A']**2 / np.sum(coef['A']**2)

            ind = PE.argsort()[::-1]

    reorderlist = ['g', 'name']
    if opt.twodim:
        reorderlist.extend(['Lsmaj', 'Lsmin', 'theta'])
        if opt.conf_int:
            reorderlist.extend(['Lsmaj_ci', 'Lsmin_ci', 'theta_ci', 'g_ci'])
    else:
        reorderlist.append('A')
        if opt.conf_int:
            reorderlist.extend(['A_ci', 'g_ci'])

    for key in reorderlist:
        coef[key] = coef[key][ind]

    coef['aux']['frq'] = coef['aux']['frq'][ind]
    coef['aux']['lind'] = coef['aux']['lind'][ind]

    if opt['RunTimeDisp']:
        print("done.")

    return coef


def _slvinit(tin, uin, vin, lat, **opts):

    if lat is None:
        raise ValueError("Latitude must be supplied")

    # Supporting only 1-D arrays for now; we can add "group"
    # support later.
    if tin.shape != uin.shape or tin.ndim != 1 or uin.ndim != 1:
        raise ValueError("t and u must be 1-D arrays")

    if vin is not None and vin.shape != uin.shape:
        raise ValueError("v must have the same shape as u")

    opt = Bunch(twodim=(vin is not None))

    # Step 0: apply epoch to time.
    tin = _normalize_time(tin, opts['epoch'])

    # Step 1: remove invalid times from tin, uin, vin
    tin = np.ma.masked_invalid(tin)
    uin = np.ma.masked_invalid(uin)
    if vin is not None:
        vin = np.ma.masked_invalid(vin)
    if np.ma.is_masked(tin):
        goodmask = ~np.ma.getmaskarray(tin)
        uin = uin.compress(goodmask)
        if vin is not None:
            vin = vin.compress(goodmask)

    tin = tin.compressed()  # No longer masked.

    # Step 2: generate t, u, v from edited tin, uin, vin.
    v = None
    if np.ma.is_masked(uin) or np.ma.is_masked(vin):
        mask = np.ma.getmaskarray(uin)
        if vin is not None:
            mask = np.ma.mask_or(np.ma.getmaskarray(vin), mask)
        goodmask = ~mask
        t = tin.compress(goodmask)
        u = uin.compress(goodmask).filled()
        if vin is not None:
            v = vin.compress(goodmask).filled()
    else:
        t = tin
        u = uin.filled()
        if vin is not None:
            v = vin.filled()

    # Now t, u, v, tin are clean ndarrays; uin and vin are masked,
    # but don't necessarily have masked values.

    # Are the times equally spaced?
    eps = np.finfo(np.float64).eps
    if np.var(np.unique(np.diff(tin))) < eps:
        opt['equi'] = True  # based on times; u/v can still have nans ("gappy")
        lor = np.ptp(tin)
        ntgood = len(tin)
        elor = lor*ntgood / (ntgood-1)
        tref = 0.5*(tin[0] + tin[-1])
    else:
        opt['equi'] = False
        lor = np.ptp(t)
        nt = len(t)
        elor = lor*nt / (nt-1)
        tref = 0.5*(t[0]+t[-1])

    # Options.
    opt['conf_int'] = True
    opt['cnstit'] = 'auto'
    opt['notrend'] = 0
    opt['prefilt'] = []
    opt['nodsatlint'] = 0
    opt['nodsatnone'] = 0
    opt['gwchlint'] = 0
    opt['gwchnone'] = 0
    opt['infer'] = None
    opt['inferaprx'] = 0
    opt['rmin'] = 1
    opt['method'] = 'ols'
    opt['tunrdn'] = 1
    opt['linci'] = 0
    opt['white'] = 0
    opt['nrlzn'] = 200
    opt['lsfrqosmp'] = 1
    opt['nodiagn'] = 0
    opt['diagnplots'] = 0
    opt['diagnminsnr'] = 2
    opt['ordercnstit'] = None
    opt['runtimedisp'] = 'yyy'

    # Update the default opt dictionary with the kwargs,
    # ensuring that every kwarg key matches a key in opt.
    for key, item in opts.items():
        try:
            opt[key] = item
        except KeyError:
            print('solve: unrecognized input: {0}'.format(key))

    return tin, t, u, v, tref, lor, elor, opt
