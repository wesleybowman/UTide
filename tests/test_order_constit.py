import numpy as np
import pytest

from utide import reconstruct, solve
from utide._ut_constants import constit_index_dict, ut_constants

ts = 735604
duration = 35

time = np.linspace(ts, ts + duration, 842)
tref = (time[-1] + time[0]) / 2

const = ut_constants.const

amps = [1.0, 0.5, 0.6, 0.1]
names = ["M2", "S2", "K1", "O1"]
cpds = [24 * const.freq[constit_index_dict[name]] for name in names]
sinusoids = []
for amp, cpd in zip(amps, cpds):
    arg = 2 * np.pi * (time - tref) * cpd
    sinusoids.append(amp * np.cos(arg))
tide = np.hstack(tuple(sinusoids)).sum(axis=0)

np.random.seed(1234)
noise = 1e-3 * np.random.randn(len(time))

time_series = tide + noise

opts0 = {
    "constit": ["K1", "M2", "O1", "S2"],
    "order_constit": "frequency",
    "phase": "raw",
    "nodal": False,
    "trend": False,
    "method": "ols",
    "conf_int": "MC",
    "Rayleigh_min": 0.95,
    "epoch": "python",
}


@pytest.mark.parametrize("conf_int", ["none", "linear", "MC"])
def test_order(conf_int):
    orders = [None, "PE", "frequency", opts0["constit"]]
    if conf_int != "none":
        orders.append("SNR")
    elevs = []
    ts_elevs = []
    vels = []
    ts_vels = []
    for order in orders:
        opts = opts0.copy()
        opts["order_constit"] = order
        opts["conf_int"] = conf_int
        elevs.append(solve(time, time_series, lat=45, **opts))
        vels.append(solve(time, time_series, time_series, lat=45, **opts))
        ts_elevs.append(reconstruct(time, elevs[-1], min_SNR=0))
        ts_vels.append(reconstruct(time, vels[-1], min_SNR=0))

    # Are the reconstructions all the same?
    for i in range(1, len(elevs)):
        assert (ts_elevs[i].h == ts_elevs[0].h).all()
        assert (ts_vels[i].u == ts_vels[0].u).all()
        assert (ts_vels[i].v == ts_vels[0].v).all()

    # Is None equivalent to "PE"? (Just a spot check.)
    assert (elevs[0].name == elevs[1].name).all()
    assert (elevs[0].A == elevs[1].A).all()


def test_invalid_snr():
    opts = opts0.copy()
    opts["conf_int"] = "none"
    opts["order_constit"] = "SNR"
    with pytest.raises(ValueError):
        solve(time, time_series, lat=45, **opts)
