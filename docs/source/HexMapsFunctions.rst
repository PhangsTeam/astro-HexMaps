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

.. function:: db.quickplot_radial_profile(line, quantity="MOM0", nbins=10, ax=None)

   Plot a binned median radial profile of a moment map.

   Requires galaxy geometry (``RGAL_KPC``) to be present in the table.
   Raises ``RuntimeError`` for non-galaxy targets.

   :param line: Line name.
   :type line: str
   :param quantity: Column prefix (``"MOM0"``, ``"TPEAK"``, etc.).
   :type quantity: str
   :param nbins: Number of radial bins.
   :type nbins: int
   :param ax: Existing Axes to plot into.
   :type ax: matplotlib.axes.Axes or None

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
