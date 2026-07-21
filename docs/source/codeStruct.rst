Code Structure
==============

Repository Layout
-----------------

.. code-block:: text

   astro-HexMaps/                       ← git root (pip install this)
   ├── hexmaps/                         ← installable Python package
   │   ├── handler_keys.py              reads & validates config.txt and key files
   │   ├── handler_targets.py           target geometry lookups
   │   ├── handler_pipeline.py          PipelineHandler: stage orchestration
   │   ├── stage_regrid.py              hex grid + convolution + sampling → .ecsv
   │   ├── stage_products.py            spectral masking, moments, shuffled spectra
   │   ├── stage_fits.py                FITS moment maps / cubes / band images
   │   ├── utils_fits.py                FITS/WCS helpers (convolution, reprojection)
   │   ├── utils_table.py               table I/O, spectral shuffle, moments
   │   ├── hexmaps_analysis.py          HexMapsAnalysis class (installed with package)
   │   ├── logger.py                    centralised stage-labelled logger
   │   ├── init_workdir.py              --init scaffolding
   │   ├── download_example.py          --download-example / --download-notebook
   │   ├── cli.py                       hexmaps console-script entry point
   │   ├── test_hexmaps.py              unit and integration tests
   │   └── templates/                   files copied by --init
   ├── config.txt                       ← example / template config file
   ├── keys/
   │   ├── target_definitions.txt       ← target geometry table (PHANGS sample)
   │   └── hfs_lines.txt                ← hyperfine structure definitions
   ├── analysis/
   │   └── hexmaps_example.ipynb        example notebook
   ├── conversion_from_pystructure/     ← migration scripts from PyStructure v4
   ├── data/                            ← example FITS input (NGC 5194)
   ├── docs/                            ← this documentation
   └── pyproject.toml


Your Working Directory
-----------------------

The installed package and your project data are completely separate.
A typical project directory looks like:

.. code-block:: text

   ~/my_survey/
   ├── config.txt               ← edit this for every run
   ├── keys/
   │   ├── target_definitions.txt   ← add your targets here
   │   └── hfs_lines.txt
   ├── data/                    ← your FITS files
   ├── output/                  ← .ecsv database written here
   ├── saved_fits_files/        ← FITS moment maps written here (optional)
   └── run_hexmaps.py

Create this layout in one command:

.. code-block:: console

   $ hexmaps --init --workdir ~/my_survey


How HexMaps Works
-----------------

The pipeline has three stages:

**1 — Regrid**
   Based on the overlay cube and the target angular resolution, all input
   maps and cubes are convolved to a common beam and sampled onto a
   hexagonal grid. The grid spacing is ``target_res / pixels_per_beam``
   (default: half-beam). The result is an Astropy ``.ecsv`` table with one
   row per hexagonal sightline.

**2 — Products**
   For each spectral cube, an S/N mask is constructed from the reference
   line(s) specified in ``ref_line``. Moment maps (mom0, mom1, mom2, Tpeak,
   rms, equivalent width) are computed for every line. Spectra are also
   shuffled by the line-of-sight velocity to enable spectral stacking.

**3 — FITS** *(optional)*
   Convolved cubes are reconstructed from the raw inputs. Moment maps are
   computed directly in PPV space and written as FITS images. This stage
   runs independently and can be combined with the regrid+products run or
   executed on its own.


Design Philosophy
-----------------

The installed package is never modified by the user. All project-specific
files (config, keys, data, outputs) live in a working directory that the
user controls, completely separate from the installation. Multiple projects
can share a single HexMaps installation, and upgrading does not affect
existing project files.

The ``.ecsv`` output format is human-readable, stores units and metadata
in a comment header, and embeds the full provenance of every run
(``config.txt``, raw FITS headers, pipeline log) directly in the file.
