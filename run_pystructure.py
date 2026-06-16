#!/usr/bin/env python3
"""
run_pystructure.py — PyStructure run script for this project.

Copy this file into your working directory (or run `pystructure --init`)
and edit the settings below.  Then execute from your working directory:

    python run_pystructure.py

Or use the installed CLI directly:

    pystructure --conf config.txt --stages regrid products --targets ngc5194
"""

import pystructurePipeline as pys

# ---------------------------------------------------------------------------
# USER SETTINGS — edit these
# ---------------------------------------------------------------------------

# Path to your configuration file
CONF_PATH = "config.txt"

# Stages to run. Choose any subset (in order) of:
#   "regrid"    – generate the hexagonal sampling grid, convolve and sample
#                 bands / cubes, write the .ecsv table
#   "products"  – process spectra, compute moments, write shuffled spectra
#   "fits"      – write FITS moment maps and band maps
# Set to None to run ALL stages.
STAGES = None  # e.g. ["regrid", "products"]

# Sources to process. Must match entries in keys/target_definitions.txt.
# Set to None to process all sources defined in config.txt [sources].
TARGETS = None  # e.g. ["ngc5194", "ngc5457"]

# ---------------------------------------------------------------------------
# RUN — no need to edit below this line
# ---------------------------------------------------------------------------

handler = pys.PipelineHandler(conf_path=CONF_PATH)

if STAGES is None:
    handler.run_all(targets=TARGETS)
else:
    handler.run_stages(stages=STAGES, targets=TARGETS)
