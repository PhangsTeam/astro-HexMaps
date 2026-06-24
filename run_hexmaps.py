#!/usr/bin/env python3
"""
run_hexmaps.py — HexMaps run script for this project.

Copy this file into your working directory (or run `hexmaps --init`)
and edit the settings below.  Then execute from your working directory:

    python run_hexmaps.py

Or use the installed CLI directly:

    hexmaps --conf config.txt --stages regrid products --targets ngc5194
"""

import hexmaps as pys

# ---------------------------------------------------------------------------
# USER SETTINGS — edit these
# ---------------------------------------------------------------------------

# Path to your configuration file
CONF_PATH = "config.txt"

# Stages to run. Choose any subset (in order) of:
#   "regrid"    – generate the hexagonal sampling grid, convolve and sample
#                 bands / cubes, write the .ecsv table
#   "products"  – process spectra, compute moments, write shuffled spectra
#   "fits"      – (optional) write FITS moment maps and band images
#
# Default: run regrid + products only. The fits stage is optional — it
# produces convenient FITS images but the primary deliverable is the .ecsv
# database. Set STAGES = ["regrid", "products", "fits"] (or STAGES = ["all"])
# to include it.
STAGES = None  # runs regrid + products (default); set ["all"] for every stage

# Sources to process. Must match entries in keys/target_definitions.txt.
# Set to None to process all sources defined in config.txt [sources].
TARGETS = None  # e.g. ["ngc5194", "ngc5457"]

# ---------------------------------------------------------------------------
# RUN — no need to edit below this line
# ---------------------------------------------------------------------------

handler = pys.PipelineHandler(conf_path=CONF_PATH)

if STAGES is None:
    # run_all() executes the default stages (regrid + products).
    # To include FITS output: set STAGES = ["regrid", "products", "fits"]
    handler.run_all(targets=TARGETS)
else:
    handler.run_stages(stages=STAGES, targets=TARGETS)
