HexMapsAnalysis Functions
=========================

The ``HexMapsAnalysis`` class is part of the installed ``hexmaps`` package.
Instantiate it with the path to an output ``.ecsv`` file:

.. code-block:: python

   from hexmaps import HexMapsAnalysis

   db = HexMapsAnalysis("output/ngc5194_hexmaps_27p0as_2025_01_01.ecsv")

The underlying Astropy table is available at ``db.struct``, and the list of
spectral lines at ``db.lines``.

----

Coordinate Extraction
---------------------

.. function:: db.get_coordinates(center=None)

   Return the sky coordinates of all sightlines.

   Works with any coordinate system (equatorial RA/Dec, galactic GLON/GLAT,
   ecliptic, etc.) by reading the coordinate column names from the table.

   If *center* is provided, coordinates are returned as offsets in arcseconds
   relative to that position.

   :param center: Reference coordinate string. For equatorial coordinates
                  use sexagesimal form ``"13:29:52.7 47:11:43"``. For
                  galactic coordinates use decimal degrees ``"202.47 47.19"``.
                  If ``None``, returns absolute decimal degree coordinates.
   :type center: str or None
   :return: ``axis1``, ``axis2`` — two 1D arrays (degrees or arcsec offsets).
   :rtype: numpy.ndarray

----

Quicklook Plots
---------------

.. function:: db.quickplot_map(line, quantity="MOM0", s=100, cmap="viridis", stretch="lin", center=None, ax=None)

   Scatter plot of a moment map for *line* on the hexagonal grid.

   Coordinate axes and x-axis inversion are set automatically based on
   the coordinate system of the overlay cube.

   :param line: Line name as it appears in the database, e.g. ``"12CO21"``.
   :type line: str
   :param quantity: Column prefix: ``"MOM0"``, ``"MOM1"``, ``"MOM2"``,
                    ``"TPEAK"``, ``"RMS"``, or ``"MAP"`` for a 2D map.
   :type quantity: str
   :param s: Marker size.
   :type s: int
   :param cmap: Matplotlib colormap name. Default: ``"viridis"``.
   :type cmap: str
   :param stretch: Colour stretch: ``"lin"``, ``"log"``, or ``"symlog"``.
   :type stretch: str
   :param center: If given, plot offset coordinates in arcseconds relative
                  to this sky position.
   :type center: str or None
   :param ax: Existing Axes to plot into. If ``None``, a new figure is created.
   :type ax: matplotlib.axes.Axes or None

.. function:: db.quickplot_spectrum(line, idx=None, show_mask=True, show_rms=True, ax=None)

   Plot the spectrum at a single sightline on the native velocity grid.

   :param line: Line name, e.g. ``"12CO21"``.
   :type line: str
   :param idx: Sightline index. If ``None``, the sightline closest to the
               target centre is used.
   :type idx: int or None
   :param show_mask: Shade the integration mask region.
   :type show_mask: bool
   :param show_rms: Overlay a horizontal line at the RMS level.
   :type show_rms: bool
   :param ax: Existing Axes to plot into.
   :type ax: matplotlib.axes.Axes or None

.. function:: db.quickplot_shuffled_spectrum(line, idx=None, ax=None)

   Plot the velocity-shuffled spectrum at a single sightline.

   :param line: Line name.
   :type line: str
   :param idx: Sightline index. Defaults to the sightline closest to the
               target centre.
   :type idx: int or None
   :param ax: Existing Axes to plot into.
   :type ax: matplotlib.axes.Axes or None

----

Binned Trends
-------------

Bin a quantity, or the ratio of two, in bins of another quantity and get the
uncertainty of every bin, e.g. a radial profile of the CO intensity or the
CO(2-1)/CO(1-0) ratio versus radius:

.. code-block:: python

   profile = db.bin_trend("RGAL_KPC", "12CO21", nbins=8)
   db.quickplot_binned(profile, logy=True)

   ratio = db.bin_ratio("RGAL_KPC", "12CO21", "12CO10", nbins=8)
   db.quickplot_binned(ratio)

Quantities are given as column names (``"RGAL_KPC"``, ``"TPEAK_12CO21"``,
``"MAP_SPIRE250"``, ...), as line names, which select the integrated
intensity (``"12CO21"`` → ``MOM0_12CO21``), or as arrays with one value per
sightline. The per-point uncertainties are looked up automatically
(``EMOM0_*`` for ``MOM0_*``, ``EMAP_*`` for ``MAP_*``, ``RMS_*`` for
``TPEAK_*``). The result is an Astropy Table with one row per bin, which can
be printed, plotted or saved with ``.write()``.

**How the bins are computed**

* **Bins:** bin *i* holds the points with ``edges[i] < x <= edges[i+1]``
  (the first bin also includes the lower edge). By default ``nbins`` bins span
  the range of *x*, linearly or, with ``logbins=True``, logarithmically.
* **Bin value:** the median (default) or mean of the points in the bin. For a
  ratio, numerator and denominator are binned separately and divided
  (a "stacked" ratio), which stays meaningful where individual sightlines are
  too faint for a ratio of their own.
* **Uncertainty:** from the per-point errors *e* of the *n* points in a bin,
  ``sqrt(mean(e**2)) / sqrt(n / oversample)``, and 1.25 times that for the
  median. Ratio uncertainties follow from Gaussian error propagation.
* **Oversampling:** neighbouring hexagonal points share a beam, so a bin of *n*
  points holds only ``n / oversample`` independent measurements. By default
  ``oversample`` is computed from the beam and the grid spacing
  (:func:`oversampling_factor`, ≈ 5.2 for half-beam sampling). Ignoring it
  (``oversample=1``) underestimates the bin uncertainties by √5.2 ≈ 2.3.
* **Scatter:** ``y_p16`` and ``y_p84`` are the 16th and 84th percentiles of
  the individual values in the bin, i.e. the spread of the data, not an
  uncertainty.
* **Limits (ratios):** where the binned numerator (denominator) is below
  ``snr_limit`` times its uncertainty and the ratio itself is not significant,
  an upper (lower) limit is given instead.

The underlying functions also work on plain arrays without a database:
``from hexmaps.utils_binning import bin_trend, bin_ratio_trend``.

.. function:: db.bin_trend(x, y, bins=None, nbins=10, logbins=False, method="median", oversample="auto", x_err="auto", y_err="auto", extra=None, keep_empty=False)

   Bin *y* in bins of *x*.

   :param x: Quantity defining the bins: column name, line name or array.
   :type x: str or array-like
   :param y: Quantity to bin: column name, line name or array.
   :type y: str or array-like
   :param bins: Explicit bin edges. Overrides *nbins* and *logbins*.
   :type bins: array-like or None
   :param nbins: Number of bins spanning the range of *x*. Default: 10.
   :type nbins: int
   :param logbins: Logarithmically spaced bins (over the positive *x*).
   :type logbins: bool
   :param method: Bin statistic, ``"median"`` (default) or ``"mean"``.
   :type method: str
   :param oversample: Sampling points per independent beam, or ``"auto"``
                      to compute it with :func:`oversampling_factor`.
   :type oversample: float or str
   :param x_err: Per-point uncertainties of *x*: ``"auto"`` (matching error
                 column, if any), ``None``, a column name or an array.
   :param y_err: Per-point uncertainties of *y*, as for *x_err*.
   :param extra: Further quantities to average in the same bins, as a list of
                 column/line names (errors looked up automatically) or a dict
                 ``{name: values}`` / ``{name: (values, errors)}``.
   :type extra: list or dict or None
   :param keep_empty: Keep rows (filled with NaN) for empty bins.
   :type keep_empty: bool
   :return: One row per bin with the columns ``bin_lo``, ``bin_hi``, ``n``,
            ``x``, ``x_err``, ``y``, ``y_err``, ``snr`` (``y / y_err``),
            ``y_p16``, ``y_p84``, ``n_pos``, ``n_neg`` (points with ``y > 0``
            and ``y <= 0``), and ``<name>``, ``<name>_err`` for each *extra*
            quantity. ``meta`` holds the bin edges, method and oversampling.
   :rtype: astropy.table.Table

.. function:: db.bin_ratio(x, numerator, denominator, bins=None, nbins=10, logbins=False, method="median", oversample="auto", snr_limit=3.0, x_err="auto", numerator_err="auto", denominator_err="auto", extra=None, keep_empty=False)

   Bin the ratio ``numerator / denominator`` in bins of *x*.

   :param numerator: Numerator: line name (→ ``MOM0``), column name or array.
   :param denominator: Denominator, as *numerator*.
   :param snr_limit: S/N threshold for upper and lower limits. Default: 3.
   :type snr_limit: float
   :param numerator_err: Per-point uncertainties of the numerator
                         (required; ``"auto"`` looks them up).
   :param denominator_err: Per-point uncertainties of the denominator.
   :return: The columns of :func:`bin_trend` (``y`` is the binned ratio,
            ``y_p16``/``y_p84`` are percentiles of the per-point ratios,
            ``n_pos`` counts points where both quantities are positive), plus
            ``y1``, ``y1_err``, ``y2``, ``y2_err`` (binned numerator and
            denominator), ``upper_limit`` and ``lower_limit`` (NaN where the
            ratio is significant).
   :rtype: astropy.table.Table

   The other parameters are as in :func:`bin_trend`.

.. function:: db.quickplot_binned(binned, ax=None, show_points=True, show_scatter=True, logx=False, logy=False, color="C0", label=None)

   Plot the result of :func:`bin_trend` or :func:`bin_ratio`: bin values with
   error bars, the 16th–84th percentile range as a shaded band, upper (▼) and
   lower (▲) limits, and the individual sightlines in grey (for ratios only
   those where both lines exceed ``snr_limit``). Bins that cannot be shown on
   a logarithmic axis (e.g. a median of 0 for undetected sightlines) are left
   out. Pass the same *ax* with different *color* and *label* to overplot
   several trends.

   :param binned: Table returned by a binning method.
   :type binned: astropy.table.Table
   :param ax: Existing Axes to plot into. If ``None``, a new figure is created.
   :type ax: matplotlib.axes.Axes or None
   :param show_points: Show the individual sightlines (when *x* and *y* were
                       given as column or line names).
   :type show_points: bool
   :param show_scatter: Shade the 16th–84th percentile range.
   :type show_scatter: bool
   :param logx: Logarithmic x axis.
   :type logx: bool
   :param logy: Logarithmic y axis.
   :type logy: bool
   :return: The Axes.
   :rtype: matplotlib.axes.Axes

.. function:: db.oversampling_factor()

   Number of sampling points per independent beam: the beam area
   (``beam_as`` in the table metadata) divided by the area per hexagonal
   sampling point, with the spacing measured from the sightline coordinates.
   About 5.2 for the default ``pixels_per_beam = 2``. Returns 1 (with a
   warning) if the beam size is not stored in the table.

   :rtype: float

----

Data Extraction
---------------

.. function:: db.get_mom0(line)

   Return the moment-0 array for *line*.

   :param line: Line name.
   :type line: str
   :return: 1D array of integrated intensities.
   :rtype: numpy.ndarray

.. function:: db.get_ratio(line1, line2, sn=5.0)

   Compute the line ratio ``line1 / line2``, with upper and lower limits
   for sightlines where only one line is detected.

   :param line1: Numerator line name.
   :type line1: str
   :param line2: Denominator line name.
   :type line2: str
   :param sn: S/N threshold for detections.
   :type sn: float
   :return: Dictionary with keys ``"ratio"``, ``"uc"``, ``"ulimit"``,
            ``"llimit"``.
   :rtype: dict

.. function:: db.get_2D_database(fname=None, save=False)

   Return a copy of the table with all ``SPEC_*`` columns removed, suitable
   for compact storage or sharing.

   :param fname: Output filename. Defaults to ``<target>_hexmaps_2D.ecsv``.
   :type fname: str or None
   :param save: If ``True``, write the table to *fname*.
   :type save: bool
   :return: Astropy Table with spectral columns removed.
   :rtype: astropy.table.Table

----

Provenance Recovery
-------------------

.. function:: db.get_config(save_to=None)

   Return the full content of the ``config.txt`` that was used to produce
   this database, as a plain-text string.

   :param save_to: If given, write the config content to this file path.
   :type save_to: str or None
   :return: The original ``config.txt`` content.
   :rtype: str

   Example::

      print(db.get_config())
      db.get_config(save_to="recovered_config.txt")

.. function:: db.get_log(save_to=None)

   Return the full pipeline log (Loading + Regrid + Products stages) embedded
   in the database at run time.

   :param save_to: If given, write the log to this file path.
   :type save_to: str or None
   :return: The pipeline log as a plain-text string.
   :rtype: str

.. function:: db.list_input_headers()

   List the labels of all raw FITS headers embedded in this database.

   Labels correspond to the input files: e.g. ``"12CO21"`` for the
   ``SPEC_12CO21`` cube, ``"SPIRE250"`` for a 2D map, ``"OVERLAY"`` for the
   overlay cube.

   :return: Sorted list of label strings.
   :rtype: list of str

   Example::

      db.list_input_headers()
      # ['12CO10', '12CO21', 'OVERLAY', 'SPIRE250']

.. function:: db.get_input_header(label)

   Return the raw FITS header of the input file identified by *label*, exactly
   as it was on disk before any pipeline processing.

   :param label: Label as returned by :func:`list_input_headers`.
                 The full metadata key (e.g. ``"input_header_12CO21"``) is
                 also accepted.
   :type label: str
   :return: The original FITS header.
   :rtype: astropy.io.fits.Header
   :raises KeyError: If *label* is not found.

   Example::

      hdr = db.get_input_header("12CO21")
      print(f"Native beam: {hdr['BMAJ'] * 3600:.1f} arcsec")
