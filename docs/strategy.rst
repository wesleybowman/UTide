Strategy
--------

In translating the algorithms from Matlab to Python, the
first step is generating a set of Python functions that
mimic their Matlab counterparts.  Gradually, however, the
code evolves to become closer to what might have been done
if the original implementation had been in Python.  This
evolution will include extensive renaming of variables and
functions to improve readability.

Array dimension ordering
^^^^^^^^^^^^^^^^^^^^^^^^
One of the basic differences between Matlab and Python
(including numpy) is that the former is built around
matrices and linear algebra operations, while the core array
data structure in numpy is the ndarray--an n-dimensional
array, not a matrix.  Arithmetic operations are
element-by-element.  A feature of ndarray operations is
their ability to act as if dimensions were added to the left
as needed to make the dimensions of two operands match.
This is called 'broadcasting', and happens automatically.
Dimensions can also be added to the right by explicit
indexing.  The default broadcasting to the left matches
another feature of numpy: by default, array storage is in C
order, with the right-most index (the column index for a 2-D
array) varying fastest.  This is
in contrast to Matlab's Fortran order, with the left index
(the row index) varying fastest.

To take full advantage of numpy broadcasting, we will need to
assume that the time index is on the right for arrays with
more than one dimension.  This is the opposite of the Matlab
UTide case.

Public interface
^^^^^^^^^^^^^^^^
The package is called `utide`, and presently has a very
simple interface, with two function: `solve` and
`reconstruct`.  These are simply English spellings of their
Matlab counterparts.  Everything else
should be considered to be private, regardless of whether it
has a leading underscore.

There is an overwhelming number of options, and many of the
Matlab names are rather cryptic.  We have made some changes
in the way options are specified, but the process is not
complete.

Options are being held internally in a `Bunch` so as to
provide both dictionary and attribute access syntax.

Time
^^^^
Time inputs are arrays of time in days relative to the epoch
given in the `epoch` keyword argument.  In the Matlab version
of utide these would be Matlab datenums; in the python version,
using Matlab datenums requires `epoch = 'matlab'`.  The default is
`epoch = 'python'`, corresponding to the `matplotlib` date
numbers.  Any other epoch can be specified using either a
string, like `'2015-01-01'`, or a Python standard library
`datetime.datetime` or `datetime.date` instance.  The numpy
`datetime64` dtype is not yet supported, nor are any Pandas
constructs.

Missing values
^^^^^^^^^^^^^^^
The `t`, `u`, `v` inputs to `solve` and the `t` input to `
reconstruct` now support any combination
of nans and masked array inputs to indicate missing values.

The degree to which masked arrays will be used internally is
unclear, but most likely their use will be highly localized.
