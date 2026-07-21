Quick Start
===========

This page walks you through running HexMaps on the NGC 5194 example data.


Step 1 — Initialise a Working Directory
-----------------------------------------

.. code-block:: console

   $ hexmaps --init --workdir ~/hexmaps_example
   $ cd ~/hexmaps_example

This creates ``config.txt``, ``keys/``, and ``run_hexmaps.py`` in the
working directory.


Step 2 — Download the Example Data
-------------------------------------

.. code-block:: console

   $ hexmaps --download-example --workdir ~/hexmaps_example

This fetches ~46 MB of NGC 5194 FITS files into ``~/hexmaps_example/data/``.
The bundled ``config.txt`` is already configured to use these files.

Use ``--force`` to re-download files that already exist.

You can also download the example Jupyter notebook:

.. code-block:: console

   $ hexmaps --download-notebook --workdir ~/hexmaps_example


.. _run_example:

Step 3 — Run the Pipeline
--------------------------

.. code-block:: console

   $ hexmaps --conf config.txt

Or equivalently:

.. code-block:: console

   $ python run_hexmaps.py

The pipeline prints progress to the terminal and writes the output to
``output/``. To also produce FITS moment maps and band images:

.. code-block:: console

   $ hexmaps --conf config.txt --stages all

.. NOTE::

   All input FITS filenames must follow the convention::

      <target_name><file_extension>

   For example: ``ngc5194_12co21.fits``, where ``ngc5194`` is the target
   name and ``_12co21.fits`` is the file extension defined in ``config.txt``.
   The target name must also appear in ``keys/target_definitions.txt``.


Step 4 — Inspect the Output
-----------------------------

See :ref:`Analysis` for how to open and explore the output database.
