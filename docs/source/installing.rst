Installation
============

Prerequisites
-------------

HexMaps requires **Python ≥ 3.10**. All dependencies are installed
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
environment.

Installing from GitHub
----------------------

.. code-block:: console

   $ git clone https://github.com/lukas-neumann-astro/astro-HexMaps.git
   $ cd astro-HexMaps
   $ pip install -e ".[dev]"

Or install directly without cloning:

.. code-block:: console

   $ pip install git+https://github.com/lukas-neumann-astro/astro-HexMaps.git

Installing from PyPI
--------------------

Once published on PyPI:

.. code-block:: console

   $ pip install astro-hexmaps

Verifying the Installation
---------------------------

.. code-block:: console

   $ hexmaps --help

To run the built-in test suite:

.. code-block:: console

   $ python -m pytest hexmaps/test_hexmaps.py -q

Migrating from PyStructure
---------------------------

If you have existing PyStructure v4.x configuration files, three standalone
migration scripts convert them to the new format:

.. code-block:: console

   $ python conversion_from_pystructure/config_conversion.py \
         PyStructure.conf config.txt

   $ python conversion_from_pystructure/target_definitions_conversion.py \
         List_Files/geometry.txt keys/target_definitions.txt

   $ python conversion_from_pystructure/hfs_lines_conversion.py \
         List_Files/hfs_lines.txt keys/hfs_lines.txt

If your band and cube definitions live in separate list files, pass them
explicitly:

.. code-block:: console

   $ python conversion_from_pystructure/config_conversion.py \
         PyStructure.conf config.txt \
         --band-list List_Files/band_list.txt \
         --cube-list List_Files/cube_list.txt
