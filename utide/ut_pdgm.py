from __future__ import absolute_import, division

import numpy as np
import scipy.signal
import matplotlib.mlab as mlab

from .ut_fbndavg import ut_fbndavg
from .ut_lmbscga import ut_lmbscga


def ut_pdgm(t, e, cfrq, equi, frqosmp):

    P = {}
    nt = len(e)
    # hn = np.hanning(nt)
    # Matches matlab hanning.
    hn = np.hanning(nt+2)
    hn = hn[1:-1]

    if equi:
        # matlab pwelch
        # pwelch(x,window,noverlap,nfft)
        # [Puu1s,allfrq] = pwelch(real(e),hn,0,nt);
        # Puu1s, allfrq = scipy.signal.welch(np.real(e), window='hanning',
        #                                    noverlap=0, nfft=nt, fs=2*np.pi)
        # allfrq, Puu1s = scipy.signal.welch(np.real(e), window='hanning',
        #                                    noverlap=0, nfft=nt, fs=2*np.pi,
        #                                    detrend='constant',
        #                                    scaling='density')

        # allfrq, Puu1s = scipy.signal.periodogram(np.real(e),
        #                                    window='hanning',
        #                                    nfft=nt, fs=2*np.pi,
        #                                    detrend='constant',
        #                                    scaling='density')
        allfrq, Puu1s = scipy.signal.welch(np.real(e), window=hn, noverlap=0,
                                           nfft=nt, fs=2*np.pi)
        # hn = mlab.window_hanning(t)
        # Puu1s, allfrq = mlab.psd(np.real(e), window=hn, noverlap=0, NFFT=nt,
        #                          Fs=2*np.pi)

    else:
        Puu1s, allfrq = ut_lmbscga(np.real(e), t, hn, frqosmp)

    # import pdb; pdb.set_trace()

    fac = (nt-1)/(2*np.pi*(t[-1]-t[0])*24)  # conv fac: rad/sample to cph
    allfrq = allfrq*fac  # to [cycle/hour] from [rad/samp]
    Puu1s = Puu1s / fac  # to [e units^2/cph] from [e units^2/(rad/samp)]

    # import pdb; pdb.set_trace()

    P['Puu'], P['fbnd'] = ut_fbndavg(Puu1s, allfrq, cfrq)

    if not np.isreal(e).all():

        if equi:
            # Pvv1s, _ = pwelch(np.imag(e), hn, 0, nt)
            temp, Pvv1s = scipy.signal.welch(np.imag(e), window=hn,
                                             noverlap=0, nfft=nt, fs=2*np.pi)
            # temp, Pvv1s = scipy.signal.welch(np.imag(e), window=hn,
            #                                  noverlap=0, nfft=nt, fs=2*np.pi)

            # Should be able to use mlab.csd.
            # Puv1s, _ = cpsd(np.real(e), np.imag(e), hn, 0, nt)
            # Pvv1s, temp = mlab.psd(np.imag(e), window=hn, noverlap=0,
            #                        NFFT=nt, Fs=2*np.pi, sides='default')

            Puv1s, temp = mlab.csd(np.real(e), np.imag(e), noverlap=0,
                                   NFFT=nt, window=hn, Fs=2*np.pi)

        else:
            Pvv1s, _ = ut_lmbscga(np.imag(e), t, hn, frqosmp)
            # FIXME: Undefined function 'ut_lmbscgc'.
            Puv1s, _ = ut_lmbscgc(np.real(e), np.imag(e), t, hn, frqosmp)
            pass

        Pvv1s = Pvv1s / fac
        P['Pvv'], _ = ut_fbndavg(Pvv1s, allfrq, cfrq)
        Puv1s = np.real(Puv1s) / fac
        P['Puv'], _ = ut_fbndavg(Puv1s, allfrq, cfrq)
        P['Puv'] = np.abs(P['Puv'])

    return P
