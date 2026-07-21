Configuration Guide
===================

This page walks through ``config.txt`` section by section. Every key has a
sensible default; you only need to set what differs from those defaults.


Meta data
---------

Metadata stored in the output ``.ecsv`` table header for provenance.

.. code-block:: ini

   [meta]
   user = Dr. Blocksberg
   comments = Example HexMaps run


Directory paths
---------------

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


Target list
-----------

.. code-block:: ini

   [targets]
   targets = ngc5194

``targets`` is a comma-separated list of target names. Each name is
prepended to the file extensions in the map and cube tables to form full
filenames.


Overlay file
------------

.. code-block:: ini

   [overlay]
   overlay_file = _12co21.fits

``overlay_file`` defines the 3D spectral cube that sets the
spatial extent and spectral axis of the hexagonal grid.


Inpute maps and cubes
---------------------

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