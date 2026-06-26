Quick Start
===========

This page walks you through running HexMaps on the bundled NGC 5194 example
data that ships with the repository.


Step 1 — Initialise a Working Directory
----------------------------------------

After installing HexMaps, create a working directory populated with template
files:

.. code-block:: console

   $ hexmaps --init --workdir ~/hexmaps_example
   $ cd ~/hexmaps_example

You will find:

* ``config.txt`` — the main configuration file
* ``keys/target_definitions.txt`` — source geometry table
* ``keys/hfs_lines.txt`` — hyperfine structure definitions (not needed for this example)
* ``run_hexmaps.py`` — a ready-to-run Python script


Step 2 — Copy the Example Data
--------------------------------

Copy the ``data/`` folder from the cloned repository into your working
directory (or point ``data_dir`` in ``config.txt`` at it directly):

.. code-block:: console

   $ cp -r /path/to/PyStructure/data ./data

The example data contains NGC 5194 CO(2–1), CO(1–0), and SPIRE 250 µm FITS
files.


.. _run_example:

Step 3 — Run the Pipeline
--------------------------

Open ``config.txt`` and verify the paths are correct, then run:

.. code-block:: console

   $ hexmaps --conf config.txt

Or equivalently:

.. code-block:: console

   $ python run_hexmaps.py

The pipeline will print progress to the terminal as it runs. The result is
an ``.ecsv`` file in the ``output/`` folder.

To also produce FITS moment maps and band images, add the fits stage:

.. code-block:: console

   $ hexmaps --conf config.txt --stages all

.. NOTE::

   All input FITS filenames must follow the convention::

      <source_name><file_extension>

   For example: ``ngc5194_12co21.fits``, where ``ngc5194`` is the source
   name and ``_12co21.fits`` is the file extension defined in ``config.txt``.
   The source name in ``config.txt`` must match both the filename prefix and
   the entry in ``keys/target_definitions.txt``.


Step 4 — Inspect the Output
-----------------------------

For information on how to open and analyse the output, see :ref:`Analysis`.
