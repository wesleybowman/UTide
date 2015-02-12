from __future__ import division
import numpy as np
import scipy.interpolate as sip
from ut_pdgm import ut_pdgm
from ut_linci import ut_linci


def ut_confidence(coef, opt, t, e, tin, tgd, uvgd, elor, xraw, xmod, W, m, B,
                  nm, nt, nc, Xu, Yu, Xv, Yv):

    # Confidence Intervals.
    print('conf. int vls... ')

    if not opt['white']:
        # band-averaged (ba) spectral densities
        if opt['equi']:
            if np.sum(tgd) > np.sum(uvgd):
                # efill = np.interp1(t, e, tin(tgd))
                # efill = np.interp(t, e, tin[tgd])
                eTemp = sip.interp1d(t, e)
                efill = eTemp(tin[tgd])
                # Fill start&/end nans w/ nearest good.
                if np.any(np.isnan(efill)):
                    ind = np.where(np.isnan(efill))[0]
                    # ind2 = ind(ind<find(~isnan(efill),1,'first'))
                    ind2 = ind[ind < np.where(~np.isnan(efill), 1, 'first')]
                    efill[ind2] = efill[np.max(ind2) + 1]
                    ind2 = ind[ind > np.where(~np.isnan(efill), 1, 'last')]
                    efill[ind2] = efill[np.min(ind2) - 1]

                ba = ut_pdgm(tin[tgd], efill, coef['aux']['frq'], 1, 0)
            else:
                ba = ut_pdgm(tin[tgd], e, coef['aux']['frq'], 1, 0)

        else:
            ba = ut_pdgm(t, e, coef['aux']['frq'], 0, opt['lsfrqosmp'])

        # import pdb; pdb.set_trace()
        # power [ (e units)^2 ] from spectral density [ (e units)^2 / cph ]
        df = 1 / (elor * 24)
        ba['Puu'] = ba['Puu'] * df

        if opt['twodim']:
            ba['Pvv'] = ba['Pvv'] * df
            ba['Puv'] = ba['Puv'] * df

        # Assign band-avg power values to NR & R freqs.
        Puu = np.zeros(coef['aux']['frq'].shape)
        if opt['twodim']:
            Pvv = np.zeros(coef['aux']['frq'].shape)
            Puv = np.zeros(coef['aux']['frq'].shape)
            # This was creating a copy of Puu and not a new array, so Puu was
            # getting overridden
            # Pvv = Puu
            # Puv = Pvv

        # import pdb; pdb.set_trace()
        for i in range(ba['Puu'].shape[0]):

            ind = np.logical_and(coef['aux']['frq'] >= ba['fbnd'][i, 0],
                                 coef['aux']['frq'] <= ba['fbnd'][i, 1])
            ind = np.where(ind[0])
            Puu[ind] = ba['Puu'][i]

            if opt['twodim']:
                Pvv[ind] = ba['Pvv'][i]
                Puv[ind] = ba['Puv'][i]

        # import pdb; pdb.set_trace()
    # varMSM = real((ctranspose(xraw)*W*xraw -
    #                ctranspose(m)*ctranspose(B)*W*xraw)/(nt-nm))
    # varMSM = np.real((np.conj(xraw).T * W * xraw -
    #                  np.conj(m).T[:,None] * np.conj(B).T * W * xraw)/(nt-nm))

    varMSM = np.real((np.dot(np.conj(xraw[:, None]).T * W, xraw[:, None]) -
                      np.dot(np.dot(np.conj(m[:, None]).T, np.conj(B).T) * W,
                             xraw[:, None]))/(nt-nm))

    # gamC = inv(ctranspose(B)*W*B)*varMSM
    gamC = np.linalg.inv(np.dot(np.conj(B).T * W, B)) * varMSM
    # gamP = inv(transpose(B)*W*B)*((transpose(xraw)*W*xraw -
    #            transpose(m)*transpose(B)*W*xraw)/(nt-nm))
    # gamP = np.dot(np.linalg.inv(np.dot(B.T * W, B)),
    #               ((xraw.T * W * xraw - m.T[:,None] * B.T * W *
    #                xraw) / (nt-nm)))

    gamP = (np.linalg.inv(np.dot(B.T * W, B)) *
            (np.dot(xraw[:, None].T * W, xraw[:, None]) -
             np.dot(np.dot(m[:, None].T, B.T) * W, xraw[:, None])) / (nt-nm))

    Gall = gamC + gamP
    Hall = gamC - gamP

    coef['g_ci'] = np.nan*np.ones(coef['g'].shape)
    # import pdb; pdb.set_trace()
    if opt['twodim']:
        # FIXME: change to np.ones_like.
        coef['Lsmaj_ci'] = np.nan * np.ones(coef['g'].shape)
        coef['Lsmin_ci'] = np.nan * np.ones(coef['g'].shape)
        coef['theta_ci'] = np.nan * np.ones(coef['g'].shape)
        # same issue with copying
        # coef['Lsmaj_ci']= coef['g_ci']
        # coef['Lsmin_ci']= coef['g_ci']
        # coef['theta_ci']= coef['g_ci']
        varcov_mCw = np.nan*np.ones((nc, 4, 4))
    else:
        coef['A_ci'] = np.nan*np.ones(coef['g'].shape)
        # coef['A_ci'] = coef['g_ci']
        varcov_mCw = np.nan * np.ones((nc, 2, 2))

    if not opt['white']:
        varcov_mCc = np.copy(varcov_mCw)
        # varcov_mCc = varcov_mCw

    # for c=1:nc
    for c in np.arange(nc):
        # G = [Gall(c,c) Gall(c,c+nc); Gall(c+nc,c) Gall(c+nc,c+nc);];
        G = np.array([[Gall[c, c], Gall[c, c+nc]],
                      [Gall[c+nc, c], Gall[c+nc, c+nc]]])
        H = np.array([[Hall[c, c], Hall[c, c+nc]],
                      [Hall[c+nc, c], Hall[c+nc, c+nc]]])
        # H = [Hall(c,c) Hall(c,c+nc); Hall(c+nc,c) Hall(c+nc,c+nc);];
        varXu = np.real(G[0, 0] + G[1, 1] + 2 * G[0, 1]) / 2
        varYu = np.real(H[0, 0] + H[1, 1] - 2 * H[0, 1]) / 2

        if opt['twodim']:
            varXv = np.real(H[0, 0] + H[1, 1] + 2 * H[0, 1]) / 2
            varYv = np.real(G[0, 0] + G[1, 1] - 2 * G[0, 1]) / 2
            # varXv = real(H(1,1)+H(2,2)+2*H(1,2))/2;
            # varYv = real(G(1,1)+G(2,2)-2*G(1,2))/2;

        if opt['linci']:  # Linearized.
            if not opt['twodim']:
                varcov_mCw[c, :, :] = np.diag(np.array([varXu, varYu]))
                if not opt['white']:
                    den = varXu + varYu
                    varXu = Puu[c]*varXu/den
                    varYu = Puu[c]*varYu/den
                    varcov_mCc[c, :, :] = np.diag(np.array([varXu, varYu]))
                sig1, sig2 = ut_linci(Xu[c], Yu[c], np.sqrt(varXu),
                                      np.sqrt(varYu))
                coef['A_ci'][c] = 1.96*sig1
                coef['g_ci'][c] = 1.96*sig2
            else:
                varcov_mCw[c, :, :] = np.diag(np.array([varXu, varYu,
                                                        varXv, varYv]))
                if not opt['white']:
                    den = varXv + varYv
                    varXv = Pvv[c] * varXv / den
                    varYv = Pvv[c] * varYv / den
                    varcov_mCc[c, :, :] = np.diag(np.array([varXu, varYu,
                                                            varXv, varYv]))
                sig1, sig2 = ut_linci(Xu[c] + 1j * Xv[c], Yu[c] + 1j * Yv[c],
                                      np.sqrt(varXu) + 1j * np.sqrt(varXv),
                                      np.sqrt(varYu) + 1j * np.sqrt(varYv))
                coef['Lsmaj_ci'][c] = 1.96*np.real(sig1)
                coef['Lsmin_ci'][c] = 1.96*np.imag(sig1)
                coef['g_ci'][c] = 1.96*np.real(sig2)
                coef['theta_ci'][c] = 1.96*np.imag(sig2)
                # import pdb; pdb.set_trace()

        else:  # TODO: Monte Carlo.
            pass

    return coef
