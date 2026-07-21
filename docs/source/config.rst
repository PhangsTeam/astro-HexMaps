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
   data_dir        = data/
   out_dir         = output/
   folder_savefits = ./saved_fits_files/
   # geom_file     = keys/target_definitions.txt
   # hfs_file      = keys/hfs_lines.txt

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


Map and Cube Tables
-------------------

2D maps and 3D spectral cubes are listed as comma-separated rows immediately
after the ``# ---- maps ----`` and ``# ---- cubes ----`` comment markers.
These markers must be present; the rows between them define which datasets
are processed.

Maps (2D)
~~~~~~~~~

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


Mask Table
----------

External masks are defined as rows after the ``# ---- mask ----`` comment
marker. Three types of entry are supported, each distinguished by its key:

.. code-block:: ini

   # ---- mask ----

   # File mask — a pre-computed binary FITS mask sampled onto the hex grid:
   # input_mask = name, description, file_ext, directory
   # input_mask = co_mask, CO signal mask, _co_mask.fits, data/

   # Velocity-window mask — channels within a fixed velocity range:
   # window_mask = name, description, v_start, v_end, unit
   # window_mask = win, Fixed velocity window, 400, 600, km/s

   # Noise velocity windows — line-free channels for RMS estimation:
   # noise_mask = name, description, v_start, v_end, unit
   # noise_mask = noise_b, Noise blue, -300, -150, km/s
   # noise_mask = noise_r, Noise red,   150,  300, km/s

``input_mask``
   An external binary FITS file that is sampled onto the hexagonal grid at
   the regrid stage and stored as a ``SPEC_<name>`` column. To use it as
   part of the signal mask, add the ``input`` token to ``ref_line``.

   *Columns:* ``name``, ``description``, ``file_ext``, ``directory``

``window_mask``
   A fixed velocity window defined by start and end velocities. All channels
   within the window are set to 1; all others to 0. To use it as part of the
   signal mask, add the ``window`` token to ``ref_line``.

   *Columns:* ``name``, ``description``, ``v_start``, ``v_end``, ``unit``
   (``unit`` must be an ``astropy.units``-readable velocity unit, e.g.
   ``km/s``).

``noise_mask``
   One or more line-free velocity windows used for per-sightline RMS
   estimation. Multiple rows are OR-combined into a single noise channel mask.
   Enabled by ``use_fixed_noise_mask = true`` in ``[masking]``.

   *Columns:* ``name``, ``description``, ``v_start``, ``v_end``, ``unit``

   The pipeline automatically excludes any noise channels that overlap with
   the signal integration mask, preventing signal contamination of the RMS
   estimate.

----


[resolution]
------------

.. code-block:: ini

   [resolution]
   target_res      = 27.0
   resolution      = angular
   pixels_per_beam = 2
   max_rad         = auto
   NAXIS_shuff     = 200
   CDELT_SHUFF     = 4000.0

``target_res``
   Target angular resolution. Interpreted according to ``resolution``:

   * ``angular`` → value in **arcseconds**
   * ``physical`` → value in **parsecs** (converted per target using
     ``dist_mpc`` from ``target_definitions.txt``)
   * ``native`` → ignored; the overlay beam is used as-is

   *Default:* ``27.0``

``resolution``
   Controls how ``target_res`` is interpreted:

   * ``angular`` — use ``target_res`` directly in arcseconds *(default)*
   * ``physical`` — convert ``target_res`` (parsecs) to arcseconds using each
     target's distance
   * ``native`` — use the native resolution of the overlay cube; no
     convolution is performed

   *Default:* ``angular``

``pixels_per_beam``
   Number of hexagonal sampling points per beam diameter. The hex-grid
   spacing is ``target_res / pixels_per_beam``.

   * ``2`` → half-beam spacing *(default, recommended)*
   * ``1`` → one sightline per beam (coarser grid, fewer sightlines)
   * ``3`` → denser grid (sightlines overlap substantially)

   *Default:* ``2``

``max_rad``
   Maximum map radius in degrees, measured from the target centre. Set to
   ``auto`` to derive the radius from the overlap between the overlay
   footprint and the data coverage.

   *Default:* ``auto``

``NAXIS_shuff``
   Number of channels in the shuffled spectrum output. The shuffled spectra
   are centred on the systemic velocity of each sightline and extend
   ±(NAXIS_shuff/2 × CDELT_SHUFF) m/s.

   *Default:* ``200``

``CDELT_SHUFF``
   Channel width of the shuffled spectrum in m/s.

   *Default:* ``4000.0``

----


[masking]
---------

.. code-block:: ini

   [masking]
   ref_line             = first
   SN_processing        = 2, 4
   strict_mask          = false
   use_fixed_noise_mask = false
   use_hfs_lines        = false
   fov_erosion_beams    = 0.5
   mom_thresh           = 5
   conseq_channels      = 3
   mom2_method          = fwhm

``ref_line``
   Comma-separated list of tokens that controls which masks are built and
   how they are combined. See the :ref:`AdvancedConfig` page for the full
   token reference and examples.

   Quick summary:

   * Line-selection: ``first`` *(default)*, ``all``, ``<n>``,
     ``<LINE_NAME>``, ``individual``
   * External-mask: ``input``, ``window``
   * Combinator: ``OR`` *(default)*, ``AND``

   *Default:* ``first``

``SN_processing``
   Two S/N thresholds ``low, high`` for the two-level mask construction.
   Channels above ``high`` seed the core mask; the core is grown into
   adjacent channels above ``low`` to capture line wings.

   *Default:* ``2, 4``

``strict_mask``
   Optional post-processing coherence filter applied after the signal mask
   is built. Options:

   * ``false`` — no additional filtering *(default)*
   * ``strict`` — remove spatially isolated detections smaller than
     approximately one beam area per channel
   * ``broad`` — re-derive the mask from a spatially smoothed cube with
     two-level S/N dilation

   See the :ref:`AdvancedConfig` page for details and caveats.

   *Default:* ``false``

``use_fixed_noise_mask``
   When ``true``, use the velocity windows defined by ``noise_mask`` rows in
   the mask table for per-sightline RMS estimation, instead of using channels
   outside the integration mask.

   Useful when the baseline contains emission from other lines that would
   otherwise bias the noise estimate.

   *Default:* ``false``

``use_hfs_lines``
   When ``true``, extend the signal mask to the hyperfine satellite
   frequencies of lines listed in ``hfs_file``. Requires ``hfs_file`` to be
   set in ``[paths]``.

   *Default:* ``false``

``fov_erosion_beams``
   Trim the effective field-of-view by this multiple of the beam FWHM. Pixels
   near the map edge where the convolution kernel extends beyond the observed
   area are biased; erosion removes them.

   * ``0`` — disable erosion; keep the full overlay footprint
   * ``0.5`` — trim by half a beam *(default, recommended minimum)*
   * ``1.0`` — conservative; trim by one full beam

   The same value is applied to the hex-grid footprint and to all FITS
   output maps so that they share a consistent effective FOV.

   *Default:* ``0.5``

``mom_thresh``
   S/N threshold for moment-1, moment-2, and equivalent-width computation.
   Sightlines with peak S/N below this value are excluded from those
   quantities (but moment-0 is still computed).

   *Default:* ``5``

``conseq_channels``
   Minimum number of consecutive channels above the S/N threshold for a
   detection to be considered valid. Isolated single-channel spikes are
   rejected even if they exceed ``SN_processing[1]``.

   *Default:* ``3``

``mom2_method``
   Definition used for the line-width (moment-2) output:

   * ``fwhm`` — convert the intensity-weighted second moment to FWHM
     (multiply by 2√(2 ln 2)) *(default)*
   * ``sqrt`` — return √(mom2), the intensity-weighted velocity dispersion
   * ``math`` — return the raw mathematical second moment

   *Default:* ``fwhm``

----


[spectral]
----------

.. code-block:: ini

   [spectral]
   spec_smooth        = default
   spec_smooth_method = binned

``spec_smooth``
   Spectral smoothing applied to each cube before sampling:

   * ``default`` — no smoothing *(default)*
   * ``overlay`` — smooth to the spectral resolution of the overlay cube
   * ``<float>`` — convolve to this resolution in km/s (e.g. ``5.0``)

   *Default:* ``default``

``spec_smooth_method``
   Algorithm used when ``spec_smooth`` is not ``default``:

   * ``binned`` — bin channels to the nearest integer ratio *(default,
     recommended)*
   * ``gauss`` — Gaussian kernel convolution. Note: this can underestimate
     RMS by 10–15 % in low-S/N regions; use with caution for science.
   * ``combined`` — bin first, then apply a Gaussian to handle the
     fractional remainder

   *Default:* ``binned``

----


[output]
--------

.. code-block:: ini

   [output]
   save_cubes    = false
   save_mom_maps = true
   save_maps     = true
   save_mask     = false

``save_cubes``
   When ``true``, write each convolved PPV cube to a FITS file in
   ``folder_savefits``. Only applies to the *fits* stage
   (``--stages all``).

   *Default:* ``false``

``save_mom_maps``
   When ``true``, write moment maps (mom0, mom1, mom2, rms, Tpeak, EW, and
   their error maps) to FITS files in ``folder_savefits``.

   *Default:* ``true``

``save_maps``
   When ``true``, write the convolved and reprojected 2D band maps to FITS
   files in ``folder_savefits``.

   *Default:* ``true``

``save_mask``
   When ``true``, write the velocity-integration mask(s) to FITS cubes in
   ``folder_savefits``. One file is written for the combined mask; additional
   files are written per line when HFS masks or per-line masks are active.

   *Default:* ``false``

----


[structure]
-----------

.. code-block:: ini

   [structure]
   structure_creation = default
   # fname_fill = ngc5194_hexmaps_27p0as_2025_01_01.ecsv

``structure_creation``
   Controls how the pipeline handles existing output ``.ecsv`` files:

   * ``default`` — create or overwrite the output file each run *(default)*
   * ``fill`` — open an existing file and add only the maps/cubes that are
     not yet present. Useful for incrementally building up a database without
     re-processing everything.
   * ``archive`` — append a timestamp to the filename and create a new
     versioned copy each run. The original file is never overwritten.

   *Default:* ``default``

``fname_fill`` *(only used when* ``structure_creation = fill``\ *)*
   Explicit filename of the existing ``.ecsv`` to open in fill mode. If not
   set, the pipeline searches ``out_dir`` for the most recent matching file.
