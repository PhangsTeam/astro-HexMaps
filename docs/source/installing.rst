Installation
============

Prerequisites
-------------

HexMaps requires **Python ≥ 3.9**. All dependencies are installed
automatically by pip:

* astropy ≥ 5.0
* numpy ≥ 1.22
* pandas ≥ 1.4
* scipy ≥ 1.7
* matplotlib ≥ 3.4
* reproject ≥ 0.9
* radio_beam ≥ 0.3.4
* spectral_cube ≥ 0.6
* scikit-image ≥ 0.19

It is strongly recommended to work inside a dedicated conda or virtual
environment to avoid dependency conflicts.

Installing from GitHub
----------------------

The development version lives on the ``rename/hexmaps`` branch of:

https://github.com/lukas-neumann-astro/PyStructure

Clone and install in editable mode:

.. code-block:: console

   $ git clone -b rename/hexmaps https://github.com/lukas-neumann-astro/PyStructure.git
   $ cd PyStructure
   $ pip install -e ".[dev]"

Or install directly without cloning:

.. code-block:: console

   $ pip install git+https://github.com/lukas-neumann-astro/PyStructure.git@rename/hexmaps

Installing from PyPI
--------------------

Once the package is published on PyPI:

.. code-block:: console

   $ pip install hexmaps

Verifying the Installation
---------------------------

After installation, verify everything is working:

.. code-block:: console

   $ hexmaps --help

You should see the HexMaps command-line help message. To run the built-in
test suite:

.. code-block:: console

   $ python -m pytest hexmaps/test_hexmaps.py -q

Migrating from PyStructure
---------------------------

If you have existing PyStructure configuration files, the
``conversion_from_pystructure/`` folder contains three standalone migration
scripts that require no additional dependencies:

.. code-block:: console

   $ python conversion_from_pystructure/config_conversion.py \
         PyStructure.conf config.txt

   $ python conversion_from_pystructure/target_definitions_conversion.py \
         List_Files/geometry.txt keys/target_definitions.txt

   $ python conversion_from_pystructure/hfs_lines_conversion.py \
         List_Files/hfs_lines.txt keys/hfs_lines.txt

Each script takes the old file as its first argument and the desired output
path as its second.
