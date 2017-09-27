UTide
=====

.. image:: https://travis-ci.org/wesleybowman/UTide.svg?branch=master
   :target: https://travis-ci.org/wesleybowman/UTide

Python re-implementation of the Matlab package UTide.

Still in heavy development--everything is subject to change!

Note: the user interface differs from the Matlab version, so
consult the Python function docstrings to see how to specify
parameters. Some functionality from the Matlab version is
not yet available.

::

    % For more information see:
    % Codiga, D.L., 2011. Unified Tidal Analysis and Prediction Using the
    % UTide Matlab Functions. Technical Report 2011-01. Graduate School
    % of Oceanography, University of Rhode Island, Narragansett, RI.
    % 59pp. ftp://www.po.gso.uri.edu/pub/downloads/codiga/pubs/
    % 2011Codiga-UTide-Report.pdf
    %
    % UTide v1p0 9/2011 d.codiga@gso.uri.edu
    % http://www.po.gso.uri.edu/~codiga/utide/utide.htm

Installation
============

.. code:: shell

    pip install utide

If you are using conda,

.. code:: shell

    conda install utide -c conda-forge


The public functions can be imported using

.. code:: python

    from utide import solve, reconstruct

To test and make sure that the package has been installed
and imported correctly, run:

.. code:: python

    from utide.tests.test_solve import test_roundtrip
    test_roundtrip()

If you have the pytest package installed, you can execute ``python setup.py test``
from within the source directory to run all tests.

A sample call would be

.. code:: python

    from utide import solve

    coef = solve(time, time_series_u, time_series_v,
                 lat=30,
                 nodal=False,
                 trend=False,
                 method='ols',
                 conf_int='linear',
                 Rayleigh_min=0.95,)


For more examples see the
`notebooks <https://nbviewer.jupyter.org/github/wesleybowman/UTide/tree/master/notebooks/>`__
folder.
