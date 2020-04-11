import os

import numpy as np
import pytest

from numpy.testing import assert_array_almost_equal, assert_array_equal

import utide

from utide.utilities import Bunch, convert_unicode_arrays, loadbunch


thisdir = os.path.dirname(__file__)
datapath = os.path.join(thisdir, "data")
HNL_path = os.path.join(datapath, "HNL2010.asc")

t, h = np.loadtxt(HNL_path, unpack=True)

epoch = "1700-01-01"

solve_kw = {
    "epoch": epoch,
    "lat": 21,
    "constit": "auto",
    "method": "ols",
    "conf_int": "linear",
}

coef_all = utide.solve(t, h, **solve_kw)

sl_month = slice(24 * 31)
coef_month = utide.solve(t[sl_month], h[sl_month], **solve_kw)

sl_week = slice(24 * 7)
coef_week = utide.solve(t[sl_week], h[sl_week], **solve_kw)

month_index_dict = dict({(name.strip(), i) for i, name in enumerate(coef_month.name)})
infer = Bunch()
infer.inferred_names = "S2", "N2", "O1"
infer.reference_names = "M2", "M2", "K1"
infer.amp_ratios, infer.phase_offsets = [], []
for ref, inf in zip(infer.reference_names, infer.inferred_names):
    iref = month_index_dict[ref]
    iinf = month_index_dict[inf]
    infer.amp_ratios.append(coef_month.A[iinf] / coef_month.A[iref])
    infer.phase_offsets.append(coef_month.g[iref] - coef_month.g[iinf])

coef_week_inf = utide.solve(t[sl_week], h[sl_week], infer=infer, **solve_kw)

runs = [
    (coef_all, "HNL2010.mat"),
    (coef_month, "HNL2010_Jan.mat"),
    (coef_week, "HNL2010_Jan_week1.mat"),
    (coef_week_inf, "HNL2010_Jan_week1_infer_S2.mat"),
]


@pytest.mark.parametrize("coef,matfile", runs)
def test_name_A_g(coef, matfile):
    coef_mat = loadbunch(os.path.join(datapath, matfile)).coef
    coef_mat = convert_unicode_arrays(coef_mat)
    assert_array_equal(coef.name, coef_mat.name)
    assert_array_almost_equal(coef.A, coef_mat.A)
    assert_array_almost_equal(coef.g, coef_mat.g)
