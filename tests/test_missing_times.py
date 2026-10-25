import numpy as np
import pytest
from numpy.testing import assert_allclose

from utide import reconstruct, solve


@pytest.mark.parametrize("kind", ["python", "matlab", "custom", "masked_date", "nat"])
@pytest.mark.parametrize("vector", [False, True])
def test_missing_times(kind, vector):
    hours = np.arange(24 * 35)
    dates = np.datetime64("2020-01-01", "h") + hours.astype("timedelta64[h]")
    missing = (hours % 19) == 7
    data = np.cos(2 * np.pi * hours / 12.42) + 0.1 * np.sin(hours)
    data[missing] = 100
    if kind in ("masked_date", "nat"):
        times = dates.copy()
        epoch = None
        if kind == "nat":
            times[missing] = np.datetime64("NaT", "h")
        else:
            times = np.ma.array(times, mask=missing)
    else:
        epoch = {"python": "python", "matlab": "matlab", "custom": "2020-01-01"}[kind]
        times = hours / 24 + {"python": 737425, "matlab": 737791, "custom": 0}[kind]
        times = np.ma.array(times, mask=missing)

    opts = {
        "epoch": epoch,
        "lat": 21,
        "constit": ["M2"],
        "conf_int": "none",
        "verbose": False,
    }
    v = 0.5 * data if vector else None
    reference = solve(
        times[~missing],
        data[~missing],
        None if v is None else v[~missing],
        **opts,
    )
    actual = solve(times, data, v, **opts)
    keys = (
        ["Lsmaj", "Lsmin", "theta", "g", "umean", "vmean"]
        if vector
        else ["A", "g", "mean"]
    )
    for key in keys:
        assert_allclose(actual[key], reference[key], rtol=1e-12, atol=1e-10)

    result = reconstruct(times, actual, epoch=epoch, verbose=False)
    expected = reconstruct(times[~missing], reference, epoch=epoch, verbose=False)
    for key in (["u", "v"] if vector else ["h"]):
        assert np.isnan(result[key][missing]).all()
        assert_allclose(result[key][~missing], expected[key], rtol=1e-12, atol=1e-10)
