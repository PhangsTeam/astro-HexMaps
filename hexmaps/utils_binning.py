"""
Binned trends with uncertainties.

Bin one quantity (``y``), or the ratio of two (``y1 / y2``), in bins of
another (``x``), and propagate the per-point uncertainties to an uncertainty
of each bin's mean or median. Used by
:meth:`hexmaps.HexMapsAnalysis.bin_trend` and
:meth:`hexmaps.HexMapsAnalysis.bin_ratio`, but works on any arrays::

    from hexmaps.utils_binning import bin_trend, bin_ratio_trend
    tab = bin_trend(radius, intensity, y_err=intensity_err, nbins=8, oversample=5.2)

Conventions
-----------
Bins
    Bin *i* contains points with ``edges[i] < x <= edges[i+1]``; the first bin
    also includes ``x == edges[0]``, so every point between the outer edges
    is counted exactly once. Without explicit ``bins``, ``nbins`` bins span
    the range of the valid ``x`` values, linearly or (``logbins=True``)
    logarithmically over the positive ``x`` values.

Bin value
    The ``"median"`` (default) or ``"mean"`` of the points in the bin. For a
    ratio, the bin value is ``stat(y1) / stat(y2)``, i.e. the ratio of the
    binned numerator and denominator ("stacked" ratio), which stays defined
    where individual sightlines are not detected.

Uncertainty of a bin
    With per-point errors ``e`` of the ``n`` points in a bin::

        error of the mean   = sqrt(mean(e**2)) / sqrt(n / oversample)
        error of the median = 1.25 * error of the mean

    ``oversample`` is the number of data points per independent resolution
    element (beam area / area per point); neighbouring points of an
    oversampled map are correlated, so only ``n / oversample`` count as
    independent. The factor 1.25 (~ sqrt(pi/2)) is the larger scatter of the
    median of Gaussian noise. Ratio errors follow from Gaussian propagation
    of the numerator and denominator bin errors.

Scatter
    ``y_p16`` and ``y_p84`` are the 15.9th and 84.1th percentiles of the
    individual values in the bin (for a ratio: of the per-point ratios), i.e.
    the +-1 sigma range of the distribution, not an uncertainty of the bin.

Result
    An :class:`astropy.table.Table` with one row per non-empty bin (all bins
    with ``keep_empty=True``); columns are listed in :func:`bin_trend` and
    :func:`bin_ratio_trend`. Units of the inputs are carried over. The bin
    edges, method and oversampling factor are stored in ``table.meta``.

Ported from ``bindata.py`` (``get_bins_1d`` and ``get_bins_ratio_1d``) by
L. Neumann; the numerical results are unchanged.
"""

__author__ = "L. Neumann"

import numpy as np
import astropy.units as au
from astropy.table import Table

#: Percentiles of the in-bin distribution reported as y_p16 / y_p84 (+-1 sigma).
PERCENTILES = (50 - 34.1, 50 + 34.1)

#: Ratio of the uncertainty of the median to that of the mean (~ sqrt(pi/2)).
MEDIAN_ERROR_FACTOR = 1.25

_STATISTICS = {"median": np.nanmedian, "mean": np.nanmean}


# ----------------------------------------------------------------------------
# Building blocks
# ----------------------------------------------------------------------------


def _values_and_unit(data):
    """Return (flattened float array, unit or None) for an array, Quantity or Column."""
    if data is None:
        return None, None
    unit = getattr(data, "unit", None)
    if isinstance(data, au.Quantity):
        data = data.value
    values = np.asarray(data, dtype=float).ravel()
    unit = None if unit is None or unit == au.dimensionless_unscaled else au.Unit(unit)
    return values, unit


def make_bin_edges(x, bins=None, nbins=10, logbins=False):
    """
    Return the bin edges.

    Parameters
    ----------
    x       : array — values to bin (only finite values are used)
    bins    : array, optional — explicit, increasing bin edges; overrides
              ``nbins`` and ``logbins``
    nbins   : int — number of bins spanning the range of *x*
    logbins : bool — logarithmically spaced bins over the positive *x*
    """
    if bins is not None:
        edges = np.asarray(bins, dtype=float).ravel()
        if edges.size < 2 or np.any(np.diff(edges) <= 0):
            raise ValueError("bins must be at least two strictly increasing edges.")
        return edges
    if int(nbins) < 1:
        raise ValueError(f"nbins must be >= 1, got {nbins}.")
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if logbins:
        x = x[x > 0]
    if x.size == 0:
        raise ValueError("No valid x values to define the bins.")
    if logbins:
        return np.logspace(np.log10(x.min()), np.log10(x.max()), int(nbins) + 1)
    return np.linspace(x.min(), x.max(), int(nbins) + 1)


def assign_bins(x, edges):
    """
    Return the bin index of each value of *x*, or -1 if outside all bins.

    Bin *i* is ``edges[i] < x <= edges[i+1]``; the first bin also includes
    ``edges[0]``.
    """
    x = np.asarray(x, dtype=float)
    idx = np.searchsorted(edges, x, side="left") - 1
    idx[x == edges[0]] = 0
    idx[(idx < 0) | (idx >= len(edges) - 1) | ~np.isfinite(x)] = -1
    return idx


def bin_uncertainty(errors, oversample=1.0, method="median"):
    """
    Uncertainty of the mean or median of a bin from the per-point *errors*.

    ``sqrt(mean(errors**2)) / sqrt(n / oversample)``, times 1.25 for the
    median. NaN errors are ignored; returns NaN if none are finite.
    """
    errors = np.asarray(errors, dtype=float)
    errors = errors[np.isfinite(errors)]
    if errors.size == 0:
        return np.nan
    error_of_mean = np.sqrt(np.mean(errors**2)) / np.sqrt(errors.size / oversample)
    return error_of_mean * (MEDIAN_ERROR_FACTOR if method == "median" else 1.0)


def _check_method(method):
    if method not in _STATISTICS:
        raise ValueError(f"method must be 'median' or 'mean', got {method!r}.")
    return _STATISTICS[method]


def _prepare_extra(extra, n):
    """Normalise ``extra`` to {name: (values, errors or None, unit)}."""
    out = {}
    for name, item in (extra or {}).items():
        values, errors = item if isinstance(item, tuple) else (item, None)
        values, unit = _values_and_unit(values)
        errors, _ = _values_and_unit(errors)
        if values.size != n or (errors is not None and errors.size != n):
            raise ValueError(f"extra quantity {name!r} does not match the size of x.")
        out[name] = (values, errors, unit)
    return out


class _Bins:
    """The bins of one binning call: selection masks and per-bin statistics."""

    def __init__(self, x, valid, bins, nbins, logbins, method, oversample, keep_empty):
        self.edges = make_bin_edges(x[valid], bins, nbins, logbins)
        idx = np.where(valid, assign_bins(x, edges=self.edges), -1)
        self.index = [
            i for i in range(len(self.edges) - 1) if keep_empty or np.any(idx == i)
        ]
        self.masks = [idx == i for i in self.index]
        self.method = method
        self.oversample = oversample
        self._stat = _STATISTICS[method]

    def apply(self, func, values):
        """``func(values[in bin])`` for every bin; NaN for empty bins."""
        return np.array([func(values[m]) if m.any() else np.nan for m in self.masks])

    def value(self, values):
        """Mean or median of *values* in every bin."""
        return self.apply(self._stat, values)

    def error(self, errors):
        """Uncertainty of the bin value from per-point *errors* (NaN if None)."""
        if errors is None:
            return np.full(len(self.masks), np.nan)
        return self.apply(lambda e: bin_uncertainty(e, self.oversample, self.method), errors)

    def percentiles(self, values):
        """(16th, 84th) percentile arrays of *values* in every bin."""
        pct = [
            np.nanpercentile(values[m], PERCENTILES) if m.any() else (np.nan, np.nan)
            for m in self.masks
        ]
        return np.array([p[0] for p in pct]), np.array([p[1] for p in pct])

    def count(self, condition):
        """Number of points per bin where boolean array *condition* is True."""
        return np.array([int(np.sum(condition[m])) for m in self.masks], dtype=int)

    def table(self, x, x_err, x_unit, extra):
        """Start the result table: edges, counts, binned x and extra quantities."""
        tab = Table()
        tab["bin_lo"] = self.edges[self.index]
        tab["bin_hi"] = self.edges[np.add(self.index, 1)]
        tab["n"] = np.array([int(m.sum()) for m in self.masks], dtype=int)
        tab["x"], tab["x_err"] = self.value(x), self.error(x_err)
        for col in ("bin_lo", "bin_hi", "x", "x_err"):
            tab[col].unit = x_unit
        for name, (values, errors, unit) in extra.items():
            tab[name], tab[f"{name}_err"] = self.value(values), self.error(errors)
            tab[name].unit = tab[f"{name}_err"].unit = unit
        tab.meta.update(
            bin_edges=self.edges, method=self.method, oversample=float(self.oversample)
        )
        return tab


def _check_sizes(n, **arrays):
    for name, arr in arrays.items():
        if arr is not None and arr.size != n:
            raise ValueError(f"{name} has {arr.size} values but x has {n}.")


# ----------------------------------------------------------------------------
# Public functions
# ----------------------------------------------------------------------------


def bin_trend(
    x,
    y,
    y_err=None,
    x_err=None,
    bins=None,
    nbins=10,
    logbins=False,
    method="median",
    oversample=1.0,
    extra=None,
    keep_empty=False,
):
    """
    Bin *y* in bins of *x*, with uncertainties.

    Parameters
    ----------
    x, y        : array-like (any shape; flattened) — Quantities/Columns keep
                  their units in the result
    y_err       : array-like, optional — per-point uncertainty of *y*
    x_err       : array-like, optional — per-point uncertainty of *x*
    bins        : array-like, optional — explicit bin edges
    nbins       : int — number of bins if *bins* is not given (default 10)
    logbins     : bool — logarithmic bins if *bins* is not given
    method      : ``"median"`` (default) or ``"mean"``
    oversample  : float — data points per independent resolution element
                  (default 1, i.e. independent points)
    extra       : dict, optional — further quantities to average in the same
                  bins, ``{name: values}`` or ``{name: (values, errors)}``
    keep_empty  : bool — keep rows for empty bins (filled with NaN)

    Returns
    -------
    astropy.table.Table, one row per bin, with columns

    ``bin_lo``, ``bin_hi``
        bin edges
    ``n``
        number of points in the bin
    ``x``, ``x_err``
        binned *x* and its uncertainty (NaN without *x_err*)
    ``y``, ``y_err``
        binned *y* and its uncertainty (NaN without *y_err*)
    ``snr``
        ``y / y_err``
    ``y_p16``, ``y_p84``
        15.9th / 84.1th percentile of the *y* values in the bin
    ``n_pos``, ``n_neg``
        number of points with ``y > 0`` and ``y <= 0``
    ``<name>``, ``<name>_err``
        for each *extra* quantity

    Points with a non-finite *x*, *y* or (if given) *y_err* are ignored.
    """
    _check_method(method)
    x, x_unit = _values_and_unit(x)
    y, y_unit = _values_and_unit(y)
    y_err, _ = _values_and_unit(y_err)
    x_err, _ = _values_and_unit(x_err)
    _check_sizes(x.size, y=y, y_err=y_err, x_err=x_err)
    extra = _prepare_extra(extra, x.size)

    valid = np.isfinite(x) & np.isfinite(y)
    if y_err is not None:
        valid &= np.isfinite(y_err)
    b = _Bins(x, valid, bins, nbins, logbins, method, oversample, keep_empty)

    tab = b.table(x, x_err, x_unit, extra)
    tab["y"], tab["y_err"] = b.value(y), b.error(y_err)
    with np.errstate(divide="ignore", invalid="ignore"):
        tab["snr"] = tab["y"] / tab["y_err"]
    tab["y_p16"], tab["y_p84"] = b.percentiles(y)
    tab["n_pos"] = b.count(y > 0)
    tab["n_neg"] = tab["n"] - tab["n_pos"]
    for col in ("y", "y_err", "y_p16", "y_p84"):
        tab[col].unit = y_unit
    return tab


def bin_ratio_trend(
    x,
    y1,
    y2,
    y1_err,
    y2_err,
    x_err=None,
    bins=None,
    nbins=10,
    logbins=False,
    method="median",
    oversample=1.0,
    snr_limit=3.0,
    extra=None,
    keep_empty=False,
):
    """
    Bin the ratio ``y1 / y2`` in bins of *x*, with uncertainties and limits.

    Numerator and denominator are binned separately and the bin ratio is
    ``stat(y1) / stat(y2)``; its uncertainty follows from Gaussian error
    propagation of the two bin uncertainties.

    Parameters
    ----------
    x, y1, y2        : array-like (any shape; flattened)
    y1_err, y2_err   : array-like — per-point uncertainties of *y1*, *y2*
    snr_limit        : float — S/N threshold for upper/lower limits (default 3)
    x_err, bins, nbins, logbins, method, oversample, extra, keep_empty
                     : as in :func:`bin_trend`

    Returns
    -------
    astropy.table.Table with the columns of :func:`bin_trend` (where ``y`` is
    the bin ratio and ``y_p16``/``y_p84`` are percentiles of the per-point
    ratios), plus

    ``y1``, ``y1_err``, ``y2``, ``y2_err``
        binned numerator and denominator with uncertainties
    ``upper_limit``
        ``snr_limit * y1_err / y2`` where the numerator is not detected
        (``y1 / y1_err < snr_limit``) and the ratio is not significant
        (``snr < snr_limit``); NaN otherwise
    ``lower_limit``
        ``y1 / (snr_limit * y2_err)`` where the denominator is not detected
        and the ratio is not significant; NaN otherwise

    ``n_pos`` counts points with both ``y1 > 0`` and ``y2 > 0``. Points with
    a non-finite *x*, *y1*, *y2*, *y1_err* or *y2_err* are ignored.
    """
    _check_method(method)
    x, x_unit = _values_and_unit(x)
    y1, y1_unit = _values_and_unit(y1)
    y2, y2_unit = _values_and_unit(y2)
    y1_err, _ = _values_and_unit(y1_err)
    y2_err, _ = _values_and_unit(y2_err)
    x_err, _ = _values_and_unit(x_err)
    if y1_err is None or y2_err is None:
        raise ValueError("bin_ratio_trend needs the uncertainties y1_err and y2_err.")
    _check_sizes(x.size, y1=y1, y2=y2, y1_err=y1_err, y2_err=y2_err, x_err=x_err)
    extra = _prepare_extra(extra, x.size)

    valid = (
        np.isfinite(x) & np.isfinite(y1) & np.isfinite(y2)
        & np.isfinite(y1_err) & np.isfinite(y2_err)
    )
    b = _Bins(x, valid, bins, nbins, logbins, method, oversample, keep_empty)

    tab = b.table(x, x_err, x_unit, extra)
    b1, e1 = b.value(y1), b.error(y1_err)
    b2, e2 = b.value(y2), b.error(y2_err)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = b1 / b2
        ratio_err = np.sqrt((e1 / b2) ** 2 + (b1 / b2**2 * e2) ** 2)
        snr = ratio / ratio_err
        not_significant = snr < snr_limit
        upper = np.where((b1 / e1 < snr_limit) & not_significant, snr_limit * e1 / b2, np.nan)
        lower = np.where((b2 / e2 < snr_limit) & not_significant, b1 / (snr_limit * e2), np.nan)
        p16, p84 = b.percentiles(y1 / y2)

    tab["y"], tab["y_err"], tab["snr"] = ratio, ratio_err, snr
    tab["y_p16"], tab["y_p84"] = p16, p84
    tab["n_pos"] = b.count((y1 > 0) & (y2 > 0))
    tab["n_neg"] = tab["n"] - tab["n_pos"]
    tab["y1"], tab["y1_err"], tab["y2"], tab["y2_err"] = b1, e1, b2, e2
    tab["upper_limit"], tab["lower_limit"] = upper, lower

    ratio_unit = None
    if y1_unit is not None or y2_unit is not None:
        ratio_unit = (y1_unit or au.dimensionless_unscaled) / (y2_unit or au.dimensionless_unscaled)
        ratio_unit = None if ratio_unit == au.dimensionless_unscaled else ratio_unit
    for col in ("y", "y_err", "y_p16", "y_p84", "upper_limit", "lower_limit"):
        tab[col].unit = ratio_unit
    tab["y1"].unit = tab["y1_err"].unit = y1_unit
    tab["y2"].unit = tab["y2_err"].unit = y2_unit
    tab.meta["snr_limit"] = float(snr_limit)
    return tab
