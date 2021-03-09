################
Usage
################

.. toctree::
   :maxdepth: 1

solve
-----
.. autofunction:: utide.solve

reconstruct
-----------
.. autofunction:: utide.reconstruct

Bunch
-----
The data structure used internally and to hold the output from
`solve` is a hybrid; it is a dictionary subclass that provides
attribute access to its data.  For example, ``coef['aux']`` and
``coef.aux`` can be used interchangeably.

.. autofunction:: utide.utilities.Bunch
