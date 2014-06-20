import cPickle as pickle
import scipy.io as sio
import numpy as np


def compareRec(coefName, matName, twodim=False):

    matRec = sio.loadmat(matName, struct_as_record=False, squeeze_me=True)

    rec = pickle.load(open(coefName, 'rb'))

    if twodim:
        U = matRec['u']
        V = matRec['v']
        print 'U'
        print np.allclose(U,rec[0])
        print 'V'
        print np.allclose(V,rec[1])
        print rec[0]-rec[1]
    else:
        ts_recon = matRec['ts_recon']
        print np.allclose(ts_recon, rec)



#    for k in coef.keys():
#        mat = eval('mcoef.' + k)
#        pyth = coef[k]
#        try:
#            if not np.allclose(mat, pyth).all():
#                print '\n' + k
#                print mat
#                print pyth
#            #print np.allclose(mat, pyth)
#        except TypeError:
#            if k == 'name':
#                names = []
#                for i in coef[k]:
#                    names.append(i.replace(' ', ''))
#
#                mat = eval('mcoef.' + k)
#                ind = np.where( mat == names)[0].shape[0]
#                size = len(names)
#                #if not mat.all() == names:
#                if not ind == size:
#                    print mat
#                    print pyth
#                    print mat == pyth
#
#            else:
#                for i in coef[k]:
#
#                    mat = eval('mcoef.' + k + '.' + i)
#                    pyth = coef[k][i]
#                    try:
#                        if not np.allclose(mat, pyth).all():
#                            print '\n' + k + ':' + i
#                            print mat
#                            print pyth
#                        #print np.allclose(mat, pyth)
#                    except TypeError:
#                        for j in coef[k][i]:
#                            try:
#                                mat = eval('mcoef.' + k + '.' + i + '.' + j)
#                                pyth = coef[k][i][j]
#                                #print '\n' + k + ':' + i + ':' + j
#                                #print mat
#                                #print pyth
#                                if type(pyth) == str:
#                                    #print mat == pyth
#                                    if not mat == pyth:
#                                        print mat
#                                        print pyth
#
#                                else:
#                                    if not np.allclose(mat, pyth).all():
#                                        print '\n' + k + ':' + i + ':' + j
#                                        print mat
#                                        print pyth
#                                    #print np.allclose(mat, pyth)
#                            except AttributeError:
#                                pass

if __name__ == '__main__':

    print 'Elev'
    compareRec('pythonrecon.p', 'matlabrecon.mat')

    print '\nVelocity'
    compareRec('speedpythonrecon.p', 'speedmatlabrecon.mat', True)
