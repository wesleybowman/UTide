import numpy as np


def ut_astron(jd):
    '''
    UT_ASTRON()
    calculate astronomical constants
    input
    jd = time [datenum UTC] (1 x nt)
    outputs
    astro = matrix [tau s h p np pp]T, units are [cycles] (6 x nt)
    ader = matrix of derivatives of astro [cycles/day] (6 x nt)
    UTide v1p0 9/2011 d.codiga@gso.uri.edu
    (copy of t_astron.m from t_tide, Pawlowicz et al 2002)
    '''

    jd = np.array([jd])
    # datenum(1899,12,31,12,0,0)
    daten = 693961.500000000
    d = jd[:] - daten
    D = d / 10000

    #args = np.array([[np.ones(jd.shape)],[d],[D*D],[D**3]]).flatten()[:,None]
    args = np.vstack((np.ones(jd.shape), d, D*D, D**3))
    sc = np.array([270.434164, 13.1763965268, -0.0000850, 0.000000039])
    hc= np.array([ 279.696678, 0.9856473354, 0.00002267,0.000000000])
    pc= np.array([ 334.329556, 0.1114040803,-0.0007739,-0.00000026])
    npc= np.array([-259.183275, 0.0529539222,-0.0001557,-0.000000050])
    ppc= np.array([281.220844, 0.0000470684, 0.0000339, 0.000000070])

    astro = np.dot(np.vstack((sc, hc, pc, npc, ppc)), args) / 360 % 1
    tau = jd%1 + astro[1,:] - astro[0,:]
    astro = np.vstack((tau,astro))

#    dargs = np.array([[np.zeros(jd.shape[0])], [np.ones(jd.shape[0])],
#                      [2.0e-4*D], [3.0e-4*D*D]]).flatten()[:,None]

    dargs = np.vstack((np.zeros(jd.shape), np.ones(jd.shape),
                      2.0e-4*D, 3.0e-4*D*D))

    ader = np.dot(np.vstack((sc,hc,pc,npc,ppc)), dargs)/360.0
    dtau = 1.0 + ader[1,:] - ader[0, :]
    ader = np.vstack((dtau, ader))

    # might need to take out depending on shape of jd
    #astro = astro.flatten()
    #ader = ader.flatten()

    return astro,ader


