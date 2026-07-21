Configuration Guide
===================

This page walks through ``config.txt`` section by section. Every key has a
sensible default; you only need to set what differs from those defaults.


[meta]
------

Metadata stored in the output ``.ecsv`` table header for provenance.

.. code-block:: ini

   [meta]
   user = Dr. Blocksberg
   comments = Example HexMaps run


[paths]
-------

All file and directory paths. Relative paths are resolved relative to the
location of ``config.txt``. The ``geom_file`` and ``hfs_file`` keys are optional; 
if not set, the pipeline will look for ``keys/target_definitions.txt`` and ``keys/hfs_lines.txt``.

.. code-block:: ini

   [paths]
   data_dir = data/
   out_dir = output/
   geom_file = keys/target_definitions.txt
   hfs_file = keys/hfs_lines.txt
   folder_savefits = ./saved_fits_files/


[targets]
---------

.. code-block:: ini

   [targets]
   targets = ngc5194

``targets`` is a comma-separated list of target names. Each name is
prepended to the file extensions in the map and cube tables to form full
filenames.


[overlay]
---------

.. code-block:: ini

   [overlay]
   overlay_file = _12co21.fits

``overlay_file`` defines the 3D spectral cube that sets the
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
   for mask construction. Put your brightest, highest-SNR line first. For 
   more advanced line-selection options, see the ``ref_line`` key in the 
   :ref:`AdvancedConfig`.


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