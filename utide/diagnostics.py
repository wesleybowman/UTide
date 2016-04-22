from __future__ import (absolute_import, division, print_function)

import numpy as np


def ut_diagn(coef, opt):

    if opt['RunTimeDisp']:
        print('diagnostics ... ', end='')
    coef['diagn'] = {}

    if opt['twodim']:
        PE = np.sum(coef['Lsmaj']**2 + coef['Lsmin']**2)
        PE = 100 * (coef['Lsmaj']**2 + coef['Lsmin']**2) / PE

        SNR = (coef['Lsmaj']**2 + coef['Lsmin']**2) / (
            (coef['Lsmaj_ci']/1.96)**2 +
            (coef['Lsmin_ci']/1.96)**2)

    else:
        PE = 100 * coef['A']**2 / np.sum(coef['A']**2)
        SNR = (coef['A']**2) / (coef['A_ci']/1.96)**2

    indPE = PE.argsort()[::-1]

    coef['diagn']['name'] = coef['name'][indPE]
    coef['diagn']['PE'] = PE[indPE]
    coef['diagn']['SNR'] = SNR[indPE]

    return coef, indPE


#    [~,indPE] = sort(PE,'descend');
#    coef.diagn.name = coef.name(indPE);
#    coef.diagn.PE = PE(indPE);
#    coef.diagn.SNR = SNR; % used in ut_diagntable; ordered by PE there
#    if opt.twodim
#        [coef.diagn,usnrc,vsnrc] = ut_diagntable(coef,cnstit,...
#            t,u,v,xmod,m,B,W,varMSM,Gall,Hall,elor,varcov_mCw,indPE);
#    else
#        [coef.diagn,usnrc,~] = ut_diagntable(coef,cnstit,...
#            t,u,[],xmod,m,B,W,varMSM,Gall,Hall,elor,varcov_mCw,indPE);
#    end
#    if opt.diagnplots
#        tmp = nan*ones(size(uin));
#        tmp(uvgd) = usnrc;
#        usnrc = tmp;
#        tmp = nan*ones(size(uin));
#        tmp(uvgd) = e;
#        e = tmp;
#        if opt.twodim
#            tmp = nan*ones(size(uin));
#            tmp(uvgd) = vsnrc;
#            vsnrc = tmp;
#            ut_diagnfigs(coef,indPE,tin,uin,vin,usnrc,vsnrc,e);
#        else
#            ut_diagnfigs(coef,indPE,tin,uin,[],usnrc,[],e);
#        end
#    end
# end
