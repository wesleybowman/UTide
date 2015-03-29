Strategy
--------

In translating the algorithms from Matlab to Python, the
first step is generating a set of Python functions that
mimic their Matlab counterparts.  Gradually, however, the
code evolves to become closer to what might have been done
if the original implementation had been in Python.

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
slightly shortened Matlab counterparts.  Everything else
should be considered to be private, regardless of whether it
has a leading underscore.

There is an overwhelming number of options; we might be able
to find ways of making this interface friendlier.

Options are being held internally in a `dict`.  We might
switch to using a `Bunch` or similar more flexible
structure.

Time
^^^^
Presently, time inputs are assumed to be Matlab `datenum`
arrays.  We need to make this more flexible, at the very
least including the ability to handle time in days since a
specified epoch. An array of Python datetime objects could
be supported, but this is not top priority. At some point
we will presumably handle the numpy datetime64 dtype, but we
can wait until it has been reworked and is no longer in a
semi-broken experimental state.  We will also need to
investigate handling whatever Pandas produces.

Missing values
^^^^^^^^^^^^^^^
We will add support for masked array inputs to the public
functions; the degree to
which masked arrays will be used internally is open to
discussion, but most likely their internal use will be at
most highly localized.


