"""
Astronomy module test.

Test is based in t_tide t_astron which is functionally identical
to u_tide ut_astron::

    octave:10> dns = [datenum(1900, 5, 5, 0, 0, 0) ...
                      datenum(2015, 5, 10, 13, 15, 30)]
    octave:11> [a, ad] = t_astron(dns)

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np

from utide.astronomy import ut_astron


def test_astron():
    # Switch from Matlab epoch to python epoch.
    dns = [t - 366 for t in [694086.000000000, 736094.552430556]]

    a_expected = np.array([[-0.190238223090586, -0.181487296022524],
                           [0.308043143259513, 0.867328798490917],
                           [0.117804920168928, 0.133410946894728],
                           [0.967220455214981, 0.966971540511103],
                           [-0.701640310306205, 0.477567604831163],
                           [0.781185288947374, 0.786679406642858]])
    ad_expected = np.array([[9.66136807802411e-01, 9.66136808053116e-01],
                            [3.66011014627454e-02, 3.66011012649486e-02],
                            [2.73790926515680e-03, 2.73790931806424e-03],
                            [3.09455773258269e-04, 3.09453963277675e-04],
                            [1.47094227256402e-04, 1.47093863142766e-04],
                            [1.30745790039597e-07, 1.30825941667529e-07]])
    a, ad = ut_astron(dns)
    np.testing.assert_array_almost_equal(a, a_expected)
    np.testing.assert_array_almost_equal(ad, ad_expected)
