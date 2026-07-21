Configuration Guide
===================

This page walks through ``config.txt`` section by section. Every key has a
sensible default; you only need to set what differs from those defaults.

``config.txt`` combines what used to be three separate files in the old
PyStructure (``master_key.txt``, ``data_key.txt``, ``config_key.txt``) into
a single INI-style file with named sections.


[meta]
------

Metadata stored in the output ``.ecsv`` table header for provenance.

.. code-block:: ini

   [meta]
   user     = Dr. Blocksberg
   comments = Example HexMaps run


[paths]
-------

All file and directory paths. Relative paths are resolved relative to the
location of ``config.txt``.

.. code-block:: ini

   [paths]
   data_dir        = data/
   out_dir         = output/
   folder_savefits = ./saved_fits_files/

   # Defaults to keys/target_definitions.txt next to config.txt
   # geom_file = keys/target_definitions.txt

   # Only needed if use_hfs_lines = true
   # hfs_file  = keys/hfs_lines.txt


[targets] and [overlay]
------------------------

.. code-block:: ini

   [targets]
   targets = ngc5194

   [overlay]
   overlay_file = _12co21.fits

``targets`` is a comma-separated list of target names. Each name is
prepended to the file extensions in the map and cube tables to form full
filenames. ``overlay_file`` defines the 3D spectral cube that sets the
spatial extent and spectral axis of the hexagonal grid.


Map and Cube Tables
-------------------

Maps (2D) and cubes (3D) are defined as comma-separated table rows
immediately after their comment markers.

.. code-block:: ini

   # ---- maps ----
   # name,  description,  unit,  file_extension,  directory,  [uc_extension]
   spire250,  SPIRE 250 um,  MJy/sr,  _spire250_gauss21.fits,  data/

   # ---- cubes ----
   # name,  description,  unit,  file_extension,  directory,  [map_ext],  [map_uc_ext]
   12co21,  12CO(2-1),  K,  _12co21.fits,  data/
   12co10,  12CO(1-0),  K,  _12co10.fits,  data/

.. IMPORTANT::

   By default, the **first cube in the list is used as the reference line**
   for mask construction. Put your brightest, highest-SNR line first.


Mask Table
----------

An optional mask section supports three types of entries, each defined by an
explicit key:

.. code-block:: ini

   # ---- mask ----

   # External FITS mask (sampled onto the hex grid at the regrid stage):
   # input_mask = name, description, file_extension, directory
   # input_mask = co_mask, CO signal mask, _co_mask.fits, data/

   # Fixed velocity window (signal integration range):
   # window_mask = name, description, v_start, v_end, unit
   # window_mask = window_mask, Fixed velocity window, 400, 600, km/s

   # Noise velocity windows (for RMS estimation):
   # noise_mask = name, description, v_start, v_end, unit
   # noise_mask = noise_mask, Noise blue, -300, -150, km/s
   # noise_mask = noise_mask, Noise red,   150,  300, km/s

To include these masks in the signal mask, add the corresponding tokens to
``ref_line`` (see the :ref:`masking section <step_masking>` below).


[resolution]
------------

.. code-block:: ini

   [resolution]
   target_res      = 27.0      # arcsec (for resolution = angular)
   resolution      = angular   # angular | physical | native
   pixels_per_beam = 2         # hex grid spacing in units of beam FWHM
   max_rad         = auto      # map radius in degrees, or "auto"
   NAXIS_shuff     = 200       # channels in shuffled spectrum
   CDELT_SHUFF     = 4000.0    # channel width of shuffled spectrum [m/s]

``resolution`` controls how ``target_res`` is interpreted:

* ``angular`` — use ``target_res`` directly in arcseconds
* ``physical`` — convert ``target_res`` (in pc) to arcseconds using the
  target distance from ``target_definitions.txt``
* ``native`` — use the native beam of the overlay cube


.. _step_masking:

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

``ref_line`` is a comma-separated list of tokens that controls which masks
are built and how they are combined.

**Line-selection tokens** (choose exactly one):

* ``first`` — first cube in the list (default)
* ``<LINE_NAME>`` — a specific named line, e.g. ``12co21``
* ``all`` — OR-combine all cubes
* ``n`` — OR-combine the first *n* cubes (integer ≥ 1)
* ``individual`` — build and apply one independent mask per line

**External-mask tokens** (optional, one or both):

* ``input`` — include the ``input_mask`` defined in the mask table
* ``window`` — include the ``window_mask`` defined in the mask table

**Combinator token** (optional, default ``OR``):

* ``OR`` — a sightline is masked if it passes *any* of the selected masks
* ``AND`` — a sightline is masked only if it passes *all* selected masks

Examples:

.. code-block:: ini

   ref_line = first                  # S/N mask from first cube only
   ref_line = 12co21                 # S/N mask from 12co21
   ref_line = first, input           # OR of first-cube mask and input mask
   ref_line = 12co21, input, AND     # 12co21 mask AND input mask
   ref_line = first, window, AND     # first-cube mask AND velocity window
   ref_line = all, input, window     # OR of all-cube + input + window masks
   ref_line = individual             # one independent mask per line


[spectral]
----------

.. code-block:: ini

   [spectral]
   spec_smooth        = default   # default | overlay | <float km/s>
   spec_smooth_method = binned    # binned | gauss | combined

.. WARNING::

   Gaussian spectral smoothing (``spec_smooth_method = gauss``) can
   underestimate RMS by 10–15%. Use ``binned`` or ``combined`` instead.


[output]
--------

.. code-block:: ini

   [output]
   save_cubes    = false   # save convolved PPV cubes (fits stage)
   save_mom_maps = true    # save moment map FITS files
   save_maps     = true    # save 2D band map FITS files
   save_mask     = false   # save the PPV mask as a FITS cube


[structure]
-----------

.. code-block:: ini

   [structure]
   structure_creation = default   # default | fill | archive
   # fname_fill = ngc5194_hexmaps_27p0as_2025_01_01.ecsv

* ``default`` — create or overwrite the ``.ecsv`` file each run
* ``fill`` — open an existing file and add only missing maps/cubes
* ``archive`` — create a new versioned file each run


.. _geomFile:

target_definitions.txt
-----------------------

The ``keys/target_definitions.txt`` file lists geometry for all targets
that may ever be processed. Add targets here once; only those listed in
``config.txt [targets]`` will be processed on any given run.

.. code-block:: text

   # target, x_ctr, y_ctr, dist_mpc, e_dist_mpc,
   #         incl_deg, e_incl_deg, posang_deg, e_posang_deg, r25, e_r25
   ngc5194, 202.4696, 47.1952, 8.58, 0.10, 22.0, 3.0, 173.0, 3.0, 3.54, 0.05

The coordinate columns (``x_ctr``, ``y_ctr``) contain the sky coordinates
of the target centre; their units and axis names are derived from the overlay
FITS header, so they work for both equatorial (RA/Dec) and galactic
(GLON/GLAT) data.

Galaxy geometry columns (``incl_deg``, ``posang_deg``, ``r25``) are
optional. Leave them blank or omit them for non-galaxy targets such as
Milky Way molecular clouds — the pipeline will skip deprojection and
print a warning.
