Configuration Guide
===================

``config.txt`` is the single file that controls every HexMaps run. It uses
a simple INI format with named sections (``[section]``); lines starting with
``#`` are comments. Relative paths are resolved relative to the location of
``config.txt`` itself.

Every key has a sensible default — you only need to set what differs from
those defaults. The sections below document each key individually.

For advanced options (resolution modes, mask combinations, spectral smoothing,
HFS correction, FOV erosion, and database fill mode) see the
:ref:`AdvancedConfig` page.

----


[meta]
------

Metadata that is stored in the header of the output ``.ecsv`` table for
provenance. Neither key affects pipeline behaviour.

.. code-block:: ini

   [meta]
   user     = Dr. Blocksberg
   comments = Example HexMaps run

``user``
   Your name or identifier. Stored as-is in the table metadata.

``comments``
   A free-form description of the run (dataset, goal, special settings, …).
   Stored as-is in the table metadata.

----


[paths]
-------

All file and directory locations. Relative paths are resolved relative to
``config.txt``.

.. code-block:: ini

   [paths]
   data_dir = data/
   out_dir = output/
   folder_savefits = ./saved_fits_files/
   # geom_file = keys/target_definitions.txt
   # hfs_file = keys/hfs_lines.txt

``data_dir``
   Directory that contains the input FITS files. The pipeline prepends the
   target name and appends the file extensions from the map and cube tables
   to build full file paths.

   *Default:* ``data/``

``out_dir``
   Directory for the output ``.ecsv`` database files.

   *Default:* ``output/``

``folder_savefits``
   Directory for FITS moment maps, cubes, and mask files written by the
   optional *fits* stage.

   *Default:* ``./saved_fits_files/``

``geom_file`` *(optional)*
   Path to the target geometry table. Must be a comma-separated file in the
   :doc:`target_definitions` format.

   *Default:* ``keys/target_definitions.txt`` (relative to ``config.txt``).
   Uncomment and set only if your file lives elsewhere.

``hfs_file`` *(optional)*
   Path to the hyperfine structure line definitions file. Only read when
   ``use_hfs_lines = true`` in ``[masking]``. Must be a comma-separated file
   in the :doc:`hfs_lines` format.

   *Default:* ``keys/hfs_lines.txt`` (relative to ``config.txt``).
   Uncomment and set only if your file lives elsewhere.

----


[targets]
---------

.. code-block:: ini

   [targets]
   targets = ngc5194, ngc5457

``targets``
   Comma-separated list of target names to process in this run. Each name
   must match the FITS filename prefix (e.g. ``ngc5194`` matches
   ``ngc5194_12co21.fits``) and must have an entry in
   ``keys/target_definitions.txt``.

   You can list many targets in ``target_definitions.txt`` and process only
   a subset here without modifying the geometry file.

----


[overlay]
---------

.. code-block:: ini

   [overlay]
   overlay_file = _12co21.fits

``overlay_file``
   File extension of the 3D spectral FITS cube used as the spatial and
   spectral reference for the hexagonal grid. The full filename is assembled
   as ``<target><overlay_file>`` (e.g. ``ngc5194_12co21.fits``).

   The overlay cube defines:

   * the WCS reference frame and pixel scale of all output products
   * the spectral axis (velocity range, channel width, rest frequency)
   * the spatial footprint within which sightlines are placed

   It must be a 3D cube. In most cases this is the brightest or most
   commonly observed line (e.g. ¹²CO(2–1)).

----

Maps (2D)
~~~~~~~~~

2D maps are listed as comma-separated rows immediately
after the ``# ---- maps ----`` comment markers.
These markers must be present; the rows between them define which datasets
are processed.

.. code-block:: ini

   # ---- maps ----
   # col 1: name         short identifier used as database key  → MAP_<name>
   # col 2: description  human-readable label
   # col 3: unit         physical unit (astropy-readable)
   # col 4: file_ext     file extension; full path = data_dir/<target><file_ext>
   # col 5: directory    directory containing the file (overrides data_dir)
   # col 6: uc_ext       uncertainty file extension (optional; leave blank if none)

   spire250,  SPIRE 250 um,  MJy/sr,  _spire250_gauss21.fits,  data/,  _spire250_gauss21_unc.fits

``name``
   Short identifier. The output column in the ``.ecsv`` table is named
   ``MAP_<name>`` (e.g. ``MAP_SPIRE250``).

``description``
   Human-readable label stored in the column description.

``unit``
   Physical unit string, readable by ``astropy.units``. Written to the
   column unit metadata and to FITS ``BUNIT`` headers.

``file_ext``
   File extension appended to the target name to build the full input path:
   ``<data_dir>/<target><file_ext>``.

``directory``
   Override directory for this file. Useful when maps live in a different
   folder from ``data_dir``.

``uc_ext`` *(optional)*
   Extension of the corresponding uncertainty map. When provided, the
   pipeline samples it at the same sightlines and stores it as
   ``MAP_<name>_UC``. Leave blank if no uncertainty map is available.

Cubes (3D)
~~~~~~~~~~

3D spectral cubes (ppv) are listed as comma-separated rows immediately
after the ``# ---- cubes ----`` comment markers.
These markers must be present; the rows between them define which datasets
are processed.

.. code-block:: ini

   # ---- cubes ----
   # col 1: name         short identifier  → SPEC_<name>, MOM0_<name>, …
   # col 2: description  human-readable label
   # col 3: unit         brightness temperature unit (K, Jy/beam, …)
   # col 4: file_ext     file extension
   # col 5: directory    directory containing the file
   # col 6: map_ext      pre-computed moment-0 map extension (optional)
   # col 7: map_uc_ext   uncertainty of the moment-0 map (optional)

   12co21,  12CO(2-1),  K,  _12co21.fits,  data/
   12co10,  12CO(1-0),  K,  _12co10.fits,  data/

``name``
   Short identifier. The pipeline generates a family of output columns:
   ``SPEC_<name>`` (spectrum), ``MOM0_<name>``, ``MOM1_<name>``,
   ``MOM2_<name>``, ``TPEAK_<name>``, ``RMS_<name>``, ``EW_<name>``, and
   their error columns.

``description``
   Human-readable label stored in the column description.

``unit``
   Brightness temperature (or flux density) unit. Used in moment computations
   and written to output FITS files.

``file_ext``
   File extension appended to the target name.

``directory``
   Override directory for this cube.

``map_ext`` *(optional)*
   Extension of a pre-computed integrated intensity (moment-0) map. When
   provided, it is sampled onto the hex grid and stored alongside the
   pipeline-computed moment maps.

``map_uc_ext`` *(optional)*
   Extension of the uncertainty map corresponding to ``map_ext``.

.. IMPORTANT::

   The **first cube in the list** is used as the default reference line for
   mask construction (``ref_line = first``). Place your brightest,
   highest-SNR spectral line first. See the :ref:`AdvancedConfig` page for
   advanced line-selection options.

----


[resolution]
------------

.. code-block:: ini

   [resolution]
   target_res = 27.0
   resolution = angular

``resolution``
   Controls how ``target_res`` is interpreted and whether the input data
   are convolved before sampling.

   * ``angular`` *(default)* — ``target_res`` is in **arcseconds** and is
     applied uniformly to every target.
   * ``physical`` — ``target_res`` is in **parsecs**.  HexMaps converts to
     arcseconds per target using ``dist_mpc`` from
     ``keys/target_definitions.txt``.  Targets without a valid distance are
     skipped with a warning.
   * ``native`` — no convolution is performed.  The overlay cube's native
     beam is used as the effective resolution.  ``target_res`` is ignored.

   *Default:* ``angular``

``target_res``
   Numeric value of the target resolution, in the units implied by
   ``resolution``:

   * arcseconds when ``resolution = angular``
   * parsecs when ``resolution = physical``
   * ignored when ``resolution = native``

   *Default:* ``27.0``

----


[masking]
---------

.. code-block:: ini

   [masking]
   ref_line = first

``ref_line``
   Comma-separated list of tokens that controls which masks are built and
   how they are combined. See the :ref:`AdvancedConfig` page for the full
   token reference and examples.

   Quick summary:

   * Line-selection: ``first`` *(default)*, ``all``, ``<n>``,
     ``<LINE_NAME>``, ``individual``

   *Default:* ``first``

----