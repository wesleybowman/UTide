import scipy.io as sio
import numpy as np
from ut_astron import ut_astron


def ut_FUV(t, tref, lind, lat, ngflgs):

    nt = len(t)
    nc = len(lind)
    ## nodsat

    if ngflgs[1]:
        F = np.ones((nt,nc))
        U = np.zeros((nt,nc))
    else:
        if ngflgs[0]:
            tt = tref
        else:
            tt = t

        ntt = len(tt)

        mat_contents = sio.loadmat('ut_constants.mat', struct_as_record=False, squeeze_me=True)
        sat = mat_contents['sat']
        const = mat_contents['const']
        shallow = mat_contents['shallow']

        astro, ader = ut_astron(tt)

        if abs(lat) < 5:
            lat=np.sign(lat)*5;

        slat = np.sin(np.pi * lat/180)
        rr = sat.amprat
        j = np.where(sat.ilatfac==1)[0]

        rr[j] = rr[j]*0.36309*(1.0-5.0*slat*slat)/slat;

        j = np.where(sat.ilatfac==2)

        rr[j]=rr[j]*2.59808*slat;

        uu = np.dot(sat.deldood, astro[3:6, :]) + sat.phcorr[:,None]*np.ones((1,ntt)) % 1

        nfreq=len(const.isat)
        mat = rr[:,None]*np.ones((1,ntt)) * np.exp(1j*2*np.pi*uu)

        F = np.ones((nfreq, ntt)) + 0j
        ind = np.unique(sat.iconst)

        for i in xrange(len(ind)):
            F[ind[i]-1, :] = 1+np.sum(mat[sat.iconst==ind[i],:], axis=0)

        #U = imag(log(F))/(2*pi); % faster than angle(F)
        U = np.imag(np.log(F)) / (2*np.pi)
        F = np.abs(F)

        for k in np.where(np.isfinite(const.ishallow))[0]:
            ik=const.ishallow[k]+np.arange(const.nshallow[k])
            ik = ik.astype(int)
            j = shallow.iname[ik-1]
            exp1 = shallow.coef[ik-1]
            exp2 = np.abs(exp1)
            temp1 = exp1*np.ones((ntt,1))
            temp2 = exp2*np.ones((ntt,1))
            temp1 = temp1.T
            temp2 = temp2.T
            F[k,:]=np.prod(F[j-1,:]**temp2,axis=0)
            U[k,:]=np.sum(U[j-1,:]*temp1,axis=0)

        F=F[lind,:].T
        U=U[lind,:].T

        if ngflgs[1]: # nodal/satellite with linearized times
            F = F[np.ones((nt,1)),:]
            U = U[np.ones((nt,1)),:]

    ## gwch (astron arg)
    if ngflgs[3]: # none (raw phase lags not greenwich phase lags)
#        if ~exist('const','var'):
#            load('ut_constants.mat','const');
#        [~,ader] = ut_astron(tref);
#        ii=isfinite(const.ishallow);
#        const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
#        for k=find(ii)'
#            ik=const.ishallow(k)+(0:const.nshallow(k)-1);
#            const.freq(k)=sum(const.freq(shallow.iname(ik)).*shallow.coef(ik));
        V = 24*(t-tref)*const.freq(lind).T
    else:
        if ngflgs[3]: # linearized times
            tt = tref
        else:
            tt = t # exact times

        ntt = len(tt)
#        if exist('astro','var')
#            if ~isequal(size(astro,2),ntt)
#                [astro,~]=ut_astron(tt');
#            end
#        else
#            [astro,~]=ut_astron(tt');
#
#        if ~exist('const','var')
#            load('ut_constants.mat')

        mat_contents = sio.loadmat('ut_constants.mat',
                                   struct_as_record=False, squeeze_me=True)
        sat = mat_contents['sat']
        const = mat_contents['const']
        shallow = mat_contents['shallow']
        astro, ader = ut_astron(tt)

        #V = np.dot(const.doodson, astro) + const.semi[:,None]*np.ones((1,ntt)) % 1
        V = np.dot(const.doodson, astro) + const.semi[:,None]*np.ones((1,ntt))
        #V = V % 1

        for k in np.where(np.isfinite(const.ishallow))[0]:
            ik=const.ishallow[k]+np.arange(const.nshallow[k])
            ik = ik.astype(int)
            j = shallow.iname[ik-1]
            exp1 = shallow.coef[ik-1]
            temp1 = exp1[:]*np.ones((ntt,1))
            temp1 = temp1.T

            V[k,:] = np.sum(V[j-1,:]*temp1,axis=0)

        V=V[lind,:].T

#        if ngflgs(3) % linearized times
#            [~,ader] = ut_astron(tref);
#            ii=isfinite(const.ishallow);
#            const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
#            for k=find(ii)'
#                ik=const.ishallow(k)+(0:const.nshallow(k)-1);
#                const.freq(k)=sum( const.freq(shallow.iname(ik)).* ...
#                    shallow.coef(ik) );
#            end
#            V = V(ones(1,nt),:) + 24*(t-tref)*const.freq(lind)';
#        end

    return F, U, V


