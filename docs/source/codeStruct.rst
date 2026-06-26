Code Structure
==============

Repository Layout
-----------------

When you clone or install HexMaps, the repository has the following structure:

.. code-block:: text

   PyStructure/                         ← git root (pip install this)
   ├── hexmaps/                         ← installable Python package
   │   ├── handler_keys.py              reads & validates config.txt and key files
   │   ├── handler_sources.py           source geometry lookups
   │   ├── handler_pipeline.py          PipelineHandler: stage orchestration
   │   ├── stage_regrid.py              hex grid + convolution + sampling → .ecsv
   │   ├── stage_products.py            spectral masking, moments, shuffled spectra
   │   ├── stage_fits.py                FITS moment maps / cubes / band images
   │   ├── utils_fits.py                FITS/WCS helpers (convolution, reprojection)
   │   ├── utils_table.py               table I/O, spectral shuffle, moments
   │   ├── hexmapsLogger.py             centralised stage-labelled logger
   │   ├── init_workdir.py              --init scaffolding
   │   ├── cli.py                       hexmaps console-script entry point
   │   └── test_hexmaps.py              unit and integration tests
   ├── config.txt                       ← example / template config file
   ├── keys/
   │   ├── target_definitions.txt       ← source geometry table
   │   └── hfs_lines.txt                ← hyperfine structure definitions
   ├── analysis/
   │   ├── hexmaps_analysis.py          HexMapsAnalysis class: quicklook plots
   │   └── hexmaps_example.ipynb        example analysis notebook
   ├── conversion_from_pystructure/     ← migration scripts from old PyStructure
   │   ├── config_conversion.py
   │   ├── target_definitions_conversion.py
   │   └── hfs_lines_conversion.py
   ├── data/                            ← example FITS input (NGC 5194)
   ├── docs/                            ← Sphinx / Read the Docs source
   ├── images/                          ← logo and screenshot
   └── run_hexmaps.py                   ← example run script


Your Working Directory
-----------------------

The installed package and your project data live completely separately.
A typical project directory looks like:

.. code-block:: text

   ~/my_survey/
   ├── config.txt               ← edit this every run
   ├── keys/
   │   ├── target_definitions.txt   ← edit once, reuse across projects
   │   └── hfs_lines.txt
   ├── data/                    ← your FITS files
   ├── output/                  ← .ecsv database written here
   ├── saved_fits_files/        ← FITS moment maps and images written here
   └── run_hexmaps.py           ← copy and edit this

Create this layout with a single command:

.. code-block:: console

   $ hexmaps --init --workdir ~/my_survey


How HexMaps Works
-----------------

The pipeline has three stages:

1. **Regrid**
   Based on a user-defined overlay cube and target angular resolution, all
   input maps and cubes are convolved to a common beam and sampled onto a
   hexagonal grid. The grid spacing is ``target_res / pixels_per_beam``
   (default: half-beam spacing). The result is an Astropy ``.ecsv`` table
   with one row per hexagonal sightline.

2. **Products**
   For each spectral cube, a two-level S/N mask is constructed from a
   user-chosen reference line. Moment maps (integrated intensity, mean
   velocity, line width, peak temperature, rms, equivalent width) are
   computed for every line. Spectra are also shuffled by the line-of-sight
   velocity to enable spectral stacking.

3. **FITS** *(optional)*
   Convolved cubes are reconstructed from the raw inputs and moment maps are
   computed directly in position–position–velocity (PPV) space, then written
   as FITS images. This stage is independent of the ``.ecsv`` database and
   can be run separately.


Design Philosophy
-----------------

HexMaps is designed around the principle that the **installed package is
never modified by the user**. All project-specific files (config, keys,
data, outputs) live in a working directory that the user controls, completely
separate from the package installation. This means:

* Multiple projects can share a single HexMaps installation.
* Upgrading HexMaps does not affect existing project files.
* The working directory is fully self-contained and portable.

The ``.ecsv`` output format (Astropy Enhanced CSV) is human-readable,
stores units and metadata in a comment header, and can be opened with
``astropy.table.Table.read()`` without any HexMaps-specific code.
