def ut_E(t,tref,frq,lind,lat,ngflgs,prefilt):

    nt = len(t)
    nc = len(lind)
    if ngflgs[1] and ngflgs[3]:
        F = np.ones((nt,nc))
        U = np.zeros((nt,nc))
        V = 24*(t-tref)*frq
    else:
        F, U, V = ut_FUV(t,tref,lind,lat,ngflgs);

    E = F * np.exp(1j*(U+V)*2*np.pi)

    #if ~isempty(prefilt)
#    if len(prefilt)!=0:
#        P=interp1(prefilt.frq,prefilt.P,frq).T
#        P( P>max(prefilt.rng) | P<min(prefilt.rng) | isnan(P) )=1;
#        E = E*P(ones(nt,1),:);

    return E
