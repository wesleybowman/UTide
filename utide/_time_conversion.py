"""
Utility for allowing flexible time input.
"""

import numpy as np

# to be added to get 1 on 1st January of year 1 from unix epoch '1970-01-01'
_DAY_TO_GREGORIAN_EPOCH = 719163

# millisecond in a day
_MS_PER_DAY = 1000 * 86400


def _date2num(date, epoch="1970-01-01 00:00:00.000"):
    """
    Numpy based date to datenum calculator.

    `date` and `epoch` can be anything parsable by np.datetime64 - string, datetime, datetime64,
    pandas datetime.

    Default `epoch` is the unix epoch 1970-01-01 00:00:00.
    """
    date = np.asarray(date)

    try:
        date = date.astype("datetime64[ms]")
    except ValueError as err:
        raise ValueError(
            f"Cannot convert date argument ({date}) to scalar or array of numpy datetime64 dtype.",
        ) from err

    try:
        epoch = np.datetime64(epoch, "ms")
    except ValueError as err:
        raise ValueError(
            f"Cannot convert epoch argument ({epoch}) to numpy datetime64 dtype.",
        ) from err

    # datenum calculation
    datenum = (date - epoch).astype(float) / _MS_PER_DAY
    return datenum


def _python_gregorian_datenum(date):
    """
    Number of days since 0000-12-31.

    Python gregorian time is 1 on 1st day of 1st year. Essentially, it means,
    the epoch for python gregorian time is 0000-12-31. With _date2num() defined
    above, this amounts to 719163 days from the unix-epoch 1970-01-01 00:00:00.
    To avoid repetitive calculation, this is defined as _DAY_TO_GREGORIAN_EPOCH.
    """
    return _date2num(date) + _DAY_TO_GREGORIAN_EPOCH


def _normalize_time(t, epoch=None):
    """
    Convert datetime or datenum array to proper input datenum array with an
    epoch from '0000-12-31' - 1st Jan of 1st year is 1.

    `t` input time or datenum array
    `epoch` either 'python', 'matlab', or np.datetime64 compatible value
    """
    t = np.asarray(t)

    if epoch is None:
        # default datetime, datetime64, or datetime array
        return _python_gregorian_datenum(t)

    if t.dtype.kind in ("if"):
        if epoch == "python":
            return t
        elif epoch == "matlab":
            return t - 366
        else:
            try:
                ofs = _python_gregorian_datenum(epoch)
            except ValueError as err:
                raise ValueError(
                    "Cannot parse epoch as string or date or datetime",
                ) from err
            else:
                return t + ofs
    else:
        raise ValueError("Can not process time array as timestamp or datenum.")
