UTide
=====

Python distribution of the MatLab package UTide

Still in Development.

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
=====

Download the zip file and unzip it.

```
python setup.py install
```
or
```
python setup.py install --user
```
if the user doesn't have access to all files.
(untested, may not work)

**Optinal Keywords**
----
These can be inputted into ut_solv, to change the default values, which are
indicated.

    conf_int='True'
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


**File Structure for locating functions**
----
If changes are made to file structure, please update.

- ut_solv
- ut_solv1
- ut_slvinit
- ut_reconstr
- ut_reconstr1
- ut_rcinit
- ut_constants.mat
- ut_cs2cep
- ut_E
- ut_FUV
- ut_cnstitsel
- ut_astron
