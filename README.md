UTide
=====

Python implementation of the MatLab package UTide

Still in heavy development--everything is subject to change.

% For more information see:
% Codiga, D.L., 2011. Unified Tidal Analysis and Prediction Using the
% UTide Matlab Functions. Technical Report 2011-01. Graduate School
% of Oceanography, University of Rhode Island, Narragansett, RI.
% 59pp. ftp://www.po.gso.uri.edu/pub/downloads/codiga/pubs/
% 2011Codiga-UTide-Report.pdf
%
% UTide v1p0 9/2011 d.codiga@gso.uri.edu
% http://www.po.gso.uri.edu/~codiga/utide/utide.htm

***Added Experimental branch, where all code for a pull request should be submitted for testing.***

Installation
============

Download the zip file and unzip it.

```
python setup.py install
```
or
```
python setup.py install --user
```
if the user doesn't have access to all files.

If you want to work on developing the package, then
```
python setup.py develop
```
is the way to go. The public functions can then be imported using
```
from utide import solve, reconstruct
```

To test and make sure that the package has been installed and imported correctly, run:
```
from utide import simple_utide_test
simple_utide_test.simple_utide_test()
```


**Under Construction**
----
The only method that is currently implemented is 'ols.' For the rest to work,
we need to supply a suitable robust fit routine.

Diagnostics is still under work (and a lot of functions within it).
It will not run diagntable or the diagnplots.

Functions that aren't finished (there may be more that I overlooked):
ut_finish
ut_diagnfigs
ut_diagnrcn
ut_diagntable
ut_cluster
ut_nearposdef
ut_lmbscgc (could use scipy.signal.lombscargle)
ut_lmbscga
ut_rundescr

A sample call would be
```
from utide import solve
coef = solve(time, time_series_u, time_series_v, lat, cnstit='auto',
               notrend=True, rmin=0.95, method='ols',
               nodiagn=True, linci=True, conf_int=True)
```


**Optional Keywords**
----
These can be supplied to **solve**, to change the default values, which are
indicated.

    conf_int=True
    cnstit='auto'
    notrend=0
    prefilt=[]
    nodsatlint=0
    nodsatnone=0
    gwchlint=0
    gwchnone=0
    infer=[]
    inferaprx=0
    rmin=1
    method='cauchy'
    tunrdn=1
    linci=0
    white=0
    nrlzn=200
    lsfrqosmp=1
    nodiagn=0
    diagnplots=0
    diagnminsnr=2
    ordercnstit=[]
    runtimedisp='yyy'

These can be supplied to **reconstruct** to change the
default values, which are indicated.

    cnstit = []
    minsnr = 2
    minpe = 0


**File Structure for locating functions**
----
When changes are made to file structure, please update.

- _solve.py: solve, _solve1, _slvinit
- _reconstruct.py: reconstruct, _reconstr1, _rcinit
- astronomy.py: ut_astron
- band_average.py: ut_fbndavg
- confidence.py: _confidence, ut_linci
- constituent_selection.py: ut_cnstitsel
- diagnostics.py: ut_diagn
- ellipse_params.py: ut_cs2cep
- harmonics.py: ut_E, ut_FUV
- periodogram.py: ut_pdgm
- utilities.py: Bunch, showmatbunch, loadmatbunch
- simple_utide_test.py: simple_utide_test
- data/ut_constants.mat

