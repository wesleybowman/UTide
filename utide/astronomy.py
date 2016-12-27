from __future__ import (absolute_import, division, print_function)

import numpy as np

# (comments based on t_tide)
# Coefficients of the formulas in the Explan. Suppl.
_sc = np.array([270.434164, 13.1763965268, -0.0000850, 0.000000039])
_hc = np.array([279.696678, 0.9856473354, 0.00002267, 0.000000000])
_pc = np.array([334.329556, 0.1114040803, -0.0007739, -0.00000026])
_npc = np.array([-259.183275, 0.0529539222, -0.0001557, -0.000000050])
# First coeff was 281.220833 in Foreman but Expl. Suppl. has 44.
_ppc = np.array([281.220844, 0.0000470684, 0.0000339, 0.000000070])

_coefs = np.vstack((_sc, _hc, _pc, _npc, _ppc))


def ut_astron(jd):
    """
    Compute the astronomical variables and their time derivatives.

    Parameters
    ----------
    jd : float, scalar or sequence
        Time (UTC) in days starting with 1 on 1 Jan. of the year 1
        in the proleptic Gregorian calendar as in
        `datetime.date.toordinal`.

    Returns
    -------
    astro : array, (6, nt)
        rows are tau, s, h, p, np, pp (cycles)
    ader : array, (6, nt)
        time derivatives of the above (cycles/day)

    Notes
    -----
    2-D arrays are always returned.

    Variables are:

    ===  ====================================================
    tau  lunar time
    s    mean longitude of the moon
    h    mean longitude of the sun
    p    mean longitude of the lunar perigee
    np   negative of the longitude of the mean ascending node
    pp   mean longitude of the perihelion (solar perigee)
    ===  ====================================================


    Based on UTide v1p0 9/2011 d.codiga@gso.uri.edu, which
    in turn came from t_tide's t_astron.m, Pawlowicz et al 2002

    For more background information from t_tide, see the t_tide_doc
    string variable in this module.

    """

    jd = np.atleast_1d(jd).flatten()

    # Shift epoch to 1899-12-31 at noon:
    # daten = 693961.500000000  Matlab datenum version

    daten = 693595.5  # Python epoch is 366 days later than Matlab's

    d = jd - daten
    D = d / 10000

    args = np.vstack((np.ones(jd.shape), d, D*D, D**3))

    astro = np.fmod((np.dot(_coefs, args) / 360), 1)

    # lunar time: fractional part of solar day
    #             plus hour angle to longitude of sun
    #             minus longitude of moon
    tau = jd % 1 + astro[1, :] - astro[0, :]
    astro = np.vstack((tau, astro))

    # derivatives (polynomial)
    dargs = np.vstack((np.zeros(jd.shape), np.ones(jd.shape),
                      2.0e-4*D, 3.0e-4*D*D))

    ader = np.dot(_coefs, dargs)/360.0
    dtau = 1.0 + ader[1, :] - ader[0, :]
    ader = np.vstack((dtau, ader))

    return astro, ader


t_tide_doc = """
The following is taken verbatim from the t_tide t_astron.m file, with
permission.

The formulae for calculating these ephemerides (other than tau)
were taken from pages 98 and 107 of the Explanatory Supplement to
the Astronomical Ephemeris and the American Ephemeris and Nautical
Almanac (1961). They require EPHEMERIS TIME (ET), now TERRESTRIAL
TIME (TT) and are based on observations made in the 1700/1800s.
In a bizarre twist, the current definition of time is derived
by reducing observations of planetary motions using these formulas.

The current world master clock is INTERNATIONAL ATOMIC TIME (TAI).
The length of the second is based on inverting the actual
locations of the planets over the period 1956-65 into "time"
using these formulas, and an offset added to keep the scale
continuous with previous defns. Thus

                 TT = TAI + 32.184 seconds.

Universal Time UT is a time scale that is 00:00 at midnight (i.e.,
based on the earth's rotation rather than on planetary motions).
Coordinated Universal Time (UTC) is kept by atomic clocks, the
length of the second is the same as for TAI but leap seconds are
inserted at intervals so that it provides UT to within 1 second.
This is necessary because the period of the earth's rotation is
slowly increasing (the day was exactly 86400 seconds around 1820,
it is now about 2 ms longer). 22 leap seconds have been added in
the last 27 years.

As of 1/1/99,    TAI = UTC + 32 seconds.

Thus,             TT = UTC + 62.184 seconds

GPS time was synchronized with UTC 6/1/1980 ( = TAI - 19 secs),
but is NOT adjusted for leap seconds. Your receiver might do this
automatically...or it might not.

Does any of this matter? The moon longitude is the fastest changing
parameter at 13 deg/day. A time error of one minute implies a
position error of less than 0.01 deg. This would almost always be
unimportant for tidal work.

The lunar time (tau) calculation requires UT as a base.  UTC is
close enough - an error of 1 second, the biggest difference that
can occur between UT and UTC, implies a Greenwich phase error of
0.01 deg.  In Doodson's definition (Proc R. Soc. A, vol 100,
reprinted in International Hydrographic Review, Appendix to
Circular Letter 4-H, 1954) mean lunar time is taken to begin at
"lunar midnight".

B. Beardsley  12/29/98, 1/11/98
R. Pawlowicz  9/1/01
Version 1.0
"""
