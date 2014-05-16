import numpy as np


def ut_linci(X,Y,sigX,sigY):
    # UT_LINCI()
    # current ellipse parameter uncertainties from cosine/sine coefficient
    # uncertainties, by linearized relations w/ correlations presumed zero
    # inputs: (two-dim case complex, one-dim case real)
    # X = Xu + i*Xv
    # Y = Yu + i*Yv
    # for Xu =real(X) = u cosine coeff; Yu =real(Y) = u sine coeff
    # Xv =imag(X) = v cosine coeff; Yv =imag(Y) = v sine coeff
    # sigX = sigXu + i*sigXv
    # sigY = sigYu + i*sigYv
    # for sigXu =real(sigX) =stddev(Xu); sigYu =real(sigY) =stddev(Yu)
    # sigXv =imag(sigX) =stddev(Xv); sigYv =imag(sigY) =stddev(Yv)
    # outputs:
    # two-dim case, complex
    # sig1 = sig_Lsmaj +1i*sig_Lsmin [same units as inputs]
    # sig2 = sig_g + 1i*sig_theta [degrees]
    # one-dim case, real
    # sig1 = sig_A [same units as inputs]
    # sig2 = sig_g [degrees]
    # UTide v1p0 9/2011 d.codiga@gso.uri.edu
    # (adapted from errell.m of t_tide, Pawlowicz et al 2002)

    X = np.array([X])
    Y = np.array([Y])
    sigX = np.array([sigX])
    sigY = np.array([sigY])
    Xu = np.real(X[:])
    sigXu = np.real(sigX)
    Yu = np.real(Y[:])
    sigYu = np.real(sigY)

    Xv = np.imag(X[:])
    sigXv = np.imag(sigX[:])
    Yv = np.imag(Y[:])
    sigYv = np.imag(sigY[:])

    rp=.5*np.sqrt((Xu+Yv)**2+(Xv-Yu)**2)
    rm=.5*np.sqrt((Xu-Yv)**2+(Xv+Yu)**2)
    sigXu2=sigXu**2
    sigYu2=sigYu**2
    sigXv2=sigXv**2
    sigYv2=sigYv**2

    ex=(Xu+Yv)/rp
    fx=(Xu-Yv)/rm
    gx=(Yu-Xv)/rp
    hx=(Yu+Xv)/rm

    # major axis
    dXu2=(.25*(ex+fx))**2
    dYu2=(.25*(gx+hx))**2
    dXv2=(.25*(hx-gx))**2
    dYv2=(.25*(ex-fx))**2
    sig1 = np.sqrt(dXu2*sigXu2+dYu2*sigYu2+dXv2*sigXv2+dYv2*sigYv2)

    # phase
    rn=2*(Xu*Yu+Xv*Yv)
    rd=Xu**2-Yu**2+Xv**2-Yv**2
    den=rn**2+rd**2
    dXu2=((rd*Yu-rn*Xu)/den)**2
    dYu2=((rd*Xu+rn*Yu)/den)**2
    dXv2=((rd*Yv-rn*Xv)/den)**2
    dYv2=((rd*Xv+rn*Yv)/den)**2
    sig2 = (180/np.pi)*np.sqrt(dXu2*sigXu2+dYu2*sigYu2+dXv2*sigXv2+dYv2*sigYv2)

    #if ~isreal(X)
    if not np.isreal(X):
        # minor axis
        dXu2=(.25*(ex-fx))**2
        dYu2=(.25*(gx-hx))**2
        dXv2=(.25*(hx+gx))**2
        dYv2=(.25*(ex+fx))**2
        sig1 = sig1 + 1j*np.sqrt(dXu2*sigXu2+dYu2*sigYu2+dXv2*sigXv2+dYv2*sigYv2)

        # orientation
        rn=2.*(Xu*Xv+Yu*Yv)
        rd=Xu**2+Yu**2-(Xv**2+Yv**2)
        den=rn**2+rd**2
        dXu2=((rd*Xv-rn*Xu)/den)**2
        dYu2=((rd*Yv-rn*Yu)/den)**2
        dXv2=((rd*Xu+rn*Xv)/den)**2
        dYv2=((rd*Yu+rn*Yv)/den)**2
        sig2 = sig2 + 1j*(180/np.pi)*np.sqrt(dXu2*sigXu2+dYu2*sigYu2 + dXv2*sigXv2+dYv2*sigYv2)

    return sig1, sig2
