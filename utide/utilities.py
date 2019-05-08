# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function)

import numpy as np
from scipy.io import loadmat

# This module began as an excerpt from the one in python-gsw.

# Based on Robert Kern's Bunch; taken from
# http://currents.soest.hawaii.edu/hgstage/pycurrents/
# pycurrents/system/utilities.py


def complex_interp(x, xp, fp, **kw):
    if fp.dtype.kind == 'c':
        fr = np.interp(x, xp, fp.real, **kw)
        fi = np.interp(x, xp, fp.imag, **kw)
        return fr + 1j*fi
    return np.interp(x, xp, fp, **kw)


class Bunch(dict):
    """
    A dictionary that also provides access via attributes.

    Additional methods update_values and update_None provide
    control over whether new keys are added to the dictionary
    when updating, and whether an attempt to add a new key is
    ignored or raises a KeyError.

    The Bunch also prints differently than a normal
    dictionary, using str() instead of repr() for its
    keys and values, and in key-sorted order.  The printing
    format can be customized by subclassing with a different
    str_ftm class attribute.  Do not assign directly to this
    class attribute, because that would substitute an instance
    attribute which would then become part of the Bunch, and
    would be reported as such by the keys() method.

    To output a string representation with
    a particular format, without subclassing, use the
    formatted() method.
    """

    str_fmt = "{0!s:<{klen}} : {1!s:>{vlen}}\n"

    def __init__(self, *args, **kwargs):
        """
        *args* can be dictionaries, bunches, or sequences of
        key,value tuples.  *kwargs* can be used to initialize
        or add key, value pairs.
        """
        dict.__init__(self)
        for arg in args:
            self.update(arg)
        self.update(kwargs)

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError("'Bunch' object has no attribute '%s'" % name)

    def __setattr__(self, name, value):
        self[name] = value

    def __str__(self):
        return self.formatted()

    def formatted(self, fmt=None, types=False):
        """
        Return a string with keys and/or values or types.

        *fmt* is a format string as used in the str.format() method.

        The str.format() method is called with key, value as positional
        arguments, and klen, vlen as kwargs.  The latter are the maxima
        of the string lengths for the keys and values, respectively,
        up to respective maxima of 20 and 40.
        """
        if fmt is None:
            fmt = self.str_fmt

        items = list(self.items())
        items.sort()

        klens = []
        vlens = []
        for i, (k, v) in enumerate(items):
            lenk = len(str(k))
            if types:
                v = type(v).__name__
            lenv = len(str(v))
            items[i] = (k, v)
            klens.append(lenk)
            vlens.append(lenv)

        klen = min(20, max(klens))
        vlen = min(40, max(vlens))
        slist = [fmt.format(key, value, klen=klen, vlen=vlen) for
                 key, value in items]
        return ''.join(slist)

    def from_pyfile(self, filename):
        """
        Read in variables from a python code file.
        """
        # We can't simply exec the code directly, because in
        # Python 3 the scoping for list comprehensions would
        # lead to a NameError.  Wrapping the code in a function
        # fixes this.
        d = dict()
        lines = ["def _temp_func():\n"]
        with open(filename) as f:
            lines.extend(["    " + line for line in f])
        lines.extend(["\n    return(locals())\n",
                      "_temp_out = _temp_func()\n",
                      "del(_temp_func)\n"])
        codetext = "".join(lines)
        code = compile(codetext, filename, 'exec')
        exec(code, globals(), d)
        self.update(d["_temp_out"])
        return self

    def update_values(self, *args, **kw):
        """
        arguments are dictionary-like; if present, they act as
        additional sources of kwargs, with the actual kwargs
        taking precedence.

        One reserved optional kwarg is "strict".  If present and
        True, then any attempt to update with keys that are not
        already in the Bunch instance will raise a KeyError.
        """
        strict = kw.pop("strict", False)
        newkw = dict()
        for d in args:
            newkw.update(d)
        newkw.update(kw)
        self._check_strict(strict, newkw)
        dsub = dict([(k, v) for (k, v) in newkw.items() if k in self])
        self.update(dsub)

    def update_None(self, *args, **kw):
        """
        Similar to update_values, except that an existing value
        will be updated only if it is None.
        """
        strict = kw.pop("strict", False)
        newkw = dict()
        for d in args:
            newkw.update(d)
        newkw.update(kw)
        self._check_strict(strict, newkw)
        dsub = dict([(k, v) for (k, v) in newkw.items()
                     if k in self and self[k] is None])
        self.update(dsub)

    def _check_strict(self, strict, kw):
        if strict:
            bad = set(kw.keys()) - set(self.keys())
            if bad:
                bk = list(bad)
                bk.sort()
                ek = list(self.keys())
                ek.sort()
                raise KeyError(
                    "Update keys %s don't match existing keys %s" % (bk, ek))


# The following functions ending with loadbunch() and showmatbunch()
# are taken from the repo
#     http://currents.soest.hawaii.edu/hgstage/pycurrents/,
# pycurrents/file/matfile.py.

def _crunch(arr, masked=True):
    """
    Handle all arrays that are not Matlab structures.
    """
    if arr.size == 1:
        arr = arr.item()  # Returns the contents.
        return arr

    # The following squeeze is discarding some information;
    # we might want to make it optional.
    arr = arr.squeeze()

    if masked and arr.dtype.kind == 'f':  # Check for complex also.
        arrm = np.ma.masked_invalid(arr)
        if arrm.count() < arrm.size:
            arr = arrm
        else:
            arr = np.array(arr)  # Copy to force a read.
    else:
        arr = np.array(arr)
    return arr


def _structured_to_bunch(arr, masked=True):
    """
    Recursively move through the structure tree, creating
    a Bunch for each structure.  When a non-structure is
    encountered, process it with crunch().
    """

    if isinstance(arr, str):
        return arr

    # A single "void" object comes from a Matlab structure.
    # Each Matlab structure field corresponds to a field in
    # a numpy structured dtype.

    if arr.dtype.kind == 'V' and arr.shape == (1, 1):
        b = Bunch()
        x = arr[0, 0]
        for name in x.dtype.names:
            b[name] = _structured_to_bunch(x[name], masked=masked)
        return b

    return _crunch(arr, masked=masked)


def _showmatbunch(b, elements=None, origin=None):
    if elements is None:
        elements = []
    if origin is None:
        origin = ''
    items = list(b.items())
    for k, v in items:
        _origin = "%s.%s" % (origin, k)
        if isinstance(v, Bunch):
            _showmatbunch(v, elements, _origin)
        else:
            if isinstance(v, str):
                slen = len(v)
                if slen < 50:
                    entry = v
                else:
                    entry = 'string, %d characters' % slen
            elif isinstance(v, np.ndarray):
                if np.ma.isMA(v):
                    entry = 'masked array, shape %s, dtype %s' % (v.shape, v.dtype)
                else:
                    entry = 'ndarray, shape %s, dtype %s' % (v.shape, v.dtype)
            else:
                entry = '%s %s' % (type(v).__name__, v)
            elements.append((_origin, entry))
    elements.sort()
    return elements


def showmatbunch(b):
    """
    Show the contents of a matfile as it has been, or would be, loaded
    by loadbunch.

    *b* can be either the name of a matfile or the output of loadbunch.

    Returns a multi-line string suitable for printing.
    """
    if isinstance(b, str):
        b = loadbunch(b)
    elist = _showmatbunch(b)
    names = [n for n, v in elist]
    namelen = min(40, max([len(n) for n in names]))
    str_fmt = "{0!s:<{namelen}} : {1!s}\n"
    strlist = [str_fmt.format(n[1:], v, namelen=namelen) for (n, v) in elist]
    return ''.join(strlist)


def loadbunch(fname, masked=True):
    """
    Wrapper for loadmat that dereferences (1,1) object arrays,
    converts floating point arrays to masked arrays, and uses
    nested Bunch objects in place of the matlab structures.
    """
    out = Bunch()
    if fname.endswith('.mat'):
        with open(fname, 'rb') as fobj:
            xx = loadmat(fobj, chars_as_strings=True)
    elif fname.endswith('.npz'):
        xx = np.load(fname, encoding='latin1', allow_pickle=True)
    else:
        raise ValueError('Unrecognized file {}'.format(fname))
    keys = [k for k in xx.keys() if not k.startswith("__")]
    for k in keys:
        out[k] = _structured_to_bunch(xx[k], masked=masked)
    return out


def convert_unicode_arrays(b):
    """
    Given a dict-like, *b*, find ndarrays of dtype unicode and
    convert them to object arrays of strings, or to a single
    string if there is only one.  The strings have trailing
    whitespace removed. A new Bunch is returned.
    """
    out = Bunch()
    for key, val in b.items():
        if isinstance(val, np.ndarray):
            if val.dtype.kind == 'O':
                newval = np.empty(shape=val.shape, dtype=val.dtype)
                for k, x in enumerate(val):
                    if isinstance(x, np.ndarray) and x.dtype.kind == 'U' and x.size == 1:
                        newval[k] = x.item()
                    else:
                        newval[k] = x
            elif val.dtype.kind == 'U' and val.ndim == 1:
                newval = np.array([s.rstrip() for s in val], dtype=object)
            else:
                newval = val
        elif isinstance(val, dict):
            newval = convert_unicode_arrays(val)
        else:
            newval = val
        out[key] = newval
    return out
