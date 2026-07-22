Target Definitions File
=======================

The ``keys/`` directory contains the two reference files that describe your
targets and spectral lines.  This page covers ``target_definitions.txt``;
see :doc:`hfs_lines` for the companion hyperfine structure file.

The ``keys/target_definitions.txt`` file stores the geometric and physical
parameters for every target that may ever be processed by HexMaps.  The
targets actually processed on a given run are selected in ``config.txt``
under ``[targets]``; all others are ignored.

Keeping all targets in a single shared file means you only need to update it
when you observe a new object, not for every new pipeline run.


File Format
-----------

The file is a comma-separated text file with no header row.  Lines beginning
with ``#`` are comments and are ignored.  Whitespace around commas is ignored,
so columns can be aligned freely for readability.

.. code-block:: text

   # target,   x_ctr,    y_ctr,   dist_mpc, e_dist_mpc, incl_deg, e_incl_deg, posang_deg, e_posang_deg, r25,   e_r25
   ngc5194,    202.4696, 47.1952, 8.58,     0.10,       22.0,     3.0,        173.0,      3.0,          3.54,  0.05


Column Descriptions
-------------------

Required columns
~~~~~~~~~~~~~~~~

``target``
   Target name.  Must match the FITS filename prefix exactly
   (e.g. ``ngc5194`` matches ``ngc5194_12co21.fits``) and must be listed in
   ``config.txt [targets]`` to be processed.  The comparison is
   case-insensitive.

``x_ctr``
   X-coordinate (longitude) of the target centre in degrees.  For equatorial
   overlays this is Right Ascension; for galactic overlays it is Galactic
   Longitude.  The coordinate system is inferred from the ``CTYPE1`` keyword
   of the overlay FITS header, so this column should always contain degrees
   in whichever system the overlay uses.

``y_ctr``
   Y-coordinate (latitude) of the target centre in degrees.  For equatorial
   overlays: Declination; for galactic overlays: Galactic Latitude.

``dist_mpc``
   Distance to the target in megaparsecs (Mpc).  Used to convert between
   angular and physical scales when ``resolution = physical`` is set, and
   stored in the output ``.ecsv`` metadata.

``e_dist_mpc``
   Uncertainty on the distance in Mpc.  Stored in the metadata but not
   currently used in any computation.  Set to ``0`` or ``NaN`` if unknown.

Optional columns
~~~~~~~~~~~~~~~~

The following columns are needed only for galaxy targets where you want
deprojected galactocentric coordinates (``RGAL_KPC``, ``RGAL_R25``,
``THETA_RAD``) in the output database.  If any of them are missing or blank
for a row, the pipeline skips deprojection for that target and prints a
warning.

``incl_deg``
   Inclination of the galaxy disk with respect to the plane of the sky, in
   degrees (0° = face-on, 90° = edge-on).

``e_incl_deg``
   Uncertainty on the inclination in degrees.

``posang_deg``
   Position angle of the galaxy's major axis, in degrees East of North.

``e_posang_deg``
   Uncertainty on the position angle in degrees.

``r25``
   Optical radius at the 25 mag arcsec⁻² isophote, in arcminutes.  Used to
   compute deprojected radii in units of ``r25`` (``RGAL_R25`` column).

``e_r25``
   Uncertainty on ``r25`` in arcminutes.


Using the Table
---------------

**Non-galaxy targets** (e.g. Galactic molecular clouds):
   Include only the four required columns.  The optional geometry columns
   can be omitted or left blank.  HexMaps will skip deprojection and will
   not produce ``RGAL_*`` output columns.

**Multiple targets in one run:**
   List all targets in this file, then select the subset to process via
   ``targets = ngc5194, ngc5457`` in ``config.txt``.

**Uncertainties:**
   Uncertainty columns are stored in the output metadata but are not
   propagated into moment error estimates at this time.


Example
-------

.. code-block:: text

   # =============================================================================
   # HexMaps target_definitions.txt
   # =============================================================================
   # target,       x_ctr,      y_ctr,      dist_mpc,  e_dist_mpc,  incl_deg,  e_incl_deg,  posang_deg,  e_posang_deg,  r25,    e_r25
   ngc5194,        202.4696,   47.1952,    8.58,      0.10,        22.0,      3.0,         173.0,       3.0,           3.54,   0.05
   ngc0628,        24.1739,    15.7836,    9.84,      0.61,        8.9,       12.2,        20.7,        1.0,           4.94,   1.28
   # Milky Way cloud (no galaxy geometry)
   g305,           196.0,      -62.0,      3.80,      0.50
