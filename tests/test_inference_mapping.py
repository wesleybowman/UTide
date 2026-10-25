from copy import deepcopy

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from utide import reconstruct, solve
from utide.utilities import Bunch


@pytest.mark.parametrize("vector", [False, True])
def test_inference_mapping(vector):
    t = np.arange(24 * 35) / 24
    u = np.cos(2 * np.pi * t * 24 / 12.42) + 0.1 * np.sin(t)
    infer = {
        "inferred_names": ["S2"],
        "reference_names": ["M2"],
        "amp_ratios": [0.5, 0.25] if vector else [0.5],
        "phase_offsets": [20, 30] if vector else [20],
    }
    original = deepcopy(infer)
    v = 0.5 * np.sin(t) if vector else None
    opts = {
        "epoch": "2020-01-01",
        "lat": 21,
        "constit": ["M2", "K1"],
        "conf_int": "none",
        "verbose": False,
    }
    expected = solve(t, u, v, infer=Bunch(infer), **opts)
    actual = solve(t, u, v, infer=infer, **opts)
    assert infer == original
    assert_array_equal(actual.name, expected.name)
    keys = ["Lsmaj", "Lsmin", "theta", "g"] if vector else ["A", "g"]
    for key in keys:
        assert_allclose(actual[key], expected[key])
    result = reconstruct(t, actual, epoch="2020-01-01", verbose=False)
    reference = reconstruct(t, expected, epoch="2020-01-01", verbose=False)
    for key in (["u", "v"] if vector else ["h"]):
        assert_allclose(result[key], reference[key])
