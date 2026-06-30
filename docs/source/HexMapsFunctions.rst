HexMapsAnalysis Functions
=========================

The ``HexMapsAnalysis`` class (in ``analysis/hexmaps_analysis.py``) provides
a set of methods for loading and working with HexMaps ``.ecsv`` databases.
Instantiate it with the path to an output file:

.. code-block:: python

   import sys
   sys.path.append("analysis/")
   from hexmaps_analysis import HexMapsAnalysis

   db = HexMapsAnalysis("output/ngc5194_hexmaps_27p0as_2025_01_01.ecsv")

The underlying Astropy table is available at ``db.struct``, and the list of
spectral lines at ``db.lines``.

----

Coordinate Extraction
---------------------

.. function:: db.get_coordinates(center=None)

   Return the RA and Dec coordinates of all sightlines.

   If *center* is provided, coordinates are returned as offsets in arcseconds
   relative to that position.

   :param center: Reference coordinate string, e.g. ``"13:29:52.7 47:11:43"``.
                  If ``None``, returns absolute decimal degree coordinates.
   :type center: str or None
   :return: ``ra``, ``dec`` — two 1D arrays (degrees or arcsec offsets).
   :rtype: numpy.ndarray

----

Quicklook Plots
---------------

.. function:: db.quickplot_map(line, s=50, cmap=None, ax=None)

   Scatter plot of the moment-0 (integrated intensity) map for *line*.

   :param line: Line name as it appears in the database, e.g. ``"12CO21"``.
   :type line: str
   :param s: Marker size. Adjust if hexagons overlap or leave gaps.
   :type s: int
   :param cmap: Matplotlib colormap. Default: ``"RdYlBu_r"``.
   :type cmap: str or None
   :param ax: Existing Axes to plot into. If ``None``, a new figure is created.
   :return: 2D scatter plot of integrated intensities.

.. function:: db.quickplot_spectrum(line, idx=None, ax=None)

   Plot the spectrum at a single sightline.

   :param line: Line name, e.g. ``"12CO21"``.
   :type line: str
   :param idx: Sightline index. If ``None``, the brightest sightline is used.
   :type idx: int or None
   :param ax: Existing Axes to plot into.

.. function:: db.quickplot_shuffled_spectrum(line, idx=None, ax=None)

   Plot the velocity-shuffled spectrum at a single sightline.

   :param line: Line name.
   :type line: str
   :param idx: Sightline index. Defaults to the brightest sightline.
   :type idx: int or None
   :param ax: Existing Axes to plot into.

.. function:: db.quickplot_radial_profile(line, ax=None)

   Plot the azimuthally averaged moment-0 radial profile.

   :param line: Line name.
   :type line: str
   :param ax: Existing Axes to plot into.

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

   Compute the line ratio ``line1 / line2``, masking sightlines below
   the S/N threshold in either line.

   :param line1: Numerator line name.
   :type line1: str
   :param line2: Denominator line name.
   :type line2: str
   :param sn: S/N threshold for sigma clipping.
   :type sn: float
   :return: Dictionary with key ``"ratio"`` containing the 1D ratio array.
   :rtype: dict

.. function:: db.get_2D_database(fname=None, save=False)

   Return a copy of the table with all ``SPEC_*`` columns removed, suitable
   for compact storage or sharing.

   :param fname: Output filename. Defaults to ``<source>_hexmaps_2D.ecsv``.
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

   The config is embedded in the ``.ecsv`` metadata at run time, making the
   database fully self-documenting.

   :param save_to: If given, write the config content to this file path.
   :type save_to: str or None
   :return: The original ``config.txt`` content.
   :rtype: str

   Example::

      print(db.get_config())
      db.get_config(save_to="recovered_config.txt")

.. function:: db.list_input_headers()

   List the labels of all raw FITS headers embedded in this database.

   Labels correspond to the table column keys: e.g. ``"12CO21"`` for
   ``SPEC_12CO21``, ``"SPIRE250"`` for ``MAP_SPIRE250``, ``"OVERLAY"`` for
   the overlay cube.

   :return: Sorted list of label strings.
   :rtype: list of str

   Example::

      db.list_input_headers()
      # ['12CO10', '12CO21', 'OVERLAY', 'SPIRE250']

.. function:: db.get_input_header(label)

   Return the raw FITS header of the input file identified by *label*, exactly
   as it was on disk before any pipeline processing.

   :param label: Label as returned by :func:`list_input_headers`, e.g.
                 ``"12CO21"``, ``"OVERLAY"``, ``"SPIRE250"``.
                 The full metadata key (``"input_header_12CO21"``) is also
                 accepted.
   :type label: str
   :return: The original FITS header.
   :rtype: astropy.io.fits.Header
   :raises KeyError: If *label* is not found. The error message lists all
                     available labels.

   Example::

      hdr = db.get_input_header("12CO21")
      print(f"Native beam: {hdr['BMAJ'] * 3600:.1f} arcsec")
      print(repr(hdr))    # print all header cards
