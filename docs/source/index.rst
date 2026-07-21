.. HexMaps documentation master file

HexMaps Documentation
=====================

**HexMaps** — *Hexagonal-grid Multi-data Analysis and Processing Software* —
is a Python package for homogenizing and analysing multi-wavelength
astronomical datasets on hexagonal grids. It is the successor to
`PyStructure (PhangsTeam) <https://github.com/PhangsTeam/PyStructure>`_,
fully rewritten as a pip-installable package with a clean CLI and a modular
stage architecture.

.. important::

   This code is actively developed. Check the
   `GitHub repository <https://github.com/lukas-neumann-astro/astro-HexMaps>`_
   for the latest changes before upgrading between versions.


Acknowledgements
----------------

HexMaps builds on the original PyStructure IDL scripts developed within the
PHANGS collaboration. The routines have been updated, extended, and fully
rewritten in Python.


List of Papers
--------------

The code has been used in the following peer-reviewed publications:

* Stuber et al. (2025), A&A, 702A, 66S
* Galić et al. (2025), arXiv:2508.15901
* Zhang et al. (2025), ApJ, 982, 21Z
* den Brok et al. (2025), ApJ, 988, 162D
* Kovačić et al. (2025), A&A, 694A, 87K
* den Brok et al. (2023), MNRAS, 526, 6347
* Eibensteiner et al. (2023), A&A, 675, 37
* Neumann et al. (2023), MNRAS, 521, 3348
* den Brok et al. (2022), A&A, 662, 89
* Eibensteiner et al. (2022), A&A, 659, 173
* den Brok et al. (2021), MNRAS, 504, 3221


.. toctree::
   :caption: Getting Started
   :maxdepth: 2
   :hidden:

   installing
   code_structure
   quickstart

.. toctree::
   :caption: Configuration Guide
   :maxdepth: 4
   :hidden:

   config
   advanced
   target_definitions
   hfs_lines

.. toctree::
   :caption: Analysis
   :maxdepth: 2
   :hidden:

   reading
   HexMapsFunctions
