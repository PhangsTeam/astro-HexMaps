#!/usr/bin/env python3
"""
run_pystructure.py — PyStructure run script for this project.

Copy this file into your working directory (or run `pystructure --init`)
and edit the settings below.  Then execute from your working directory:

    python run_pystructure.py

Or use the installed CLI directly:

    pystructure --key_dir keys/ --stages regrid spectra --targets ngc5194
"""

import pystructure as pys

# ---------------------------------------------------------------------------
# USER SETTINGS — edit these
# ---------------------------------------------------------------------------

# Path to the directory containing your key files
KEY_DIR = "keys/"

# Stages to run. Choose any subset (in order) of:
#   "sampling"  – generate hexagonal grid
#   "regrid"    – convolve and sample bands / cubes
#   "spectra"   – process spectra, compute moments
#   "output"    – write FITS moment and band maps
# Set to None to run ALL stages.
STAGES = None  # e.g. ["sampling", "regrid"]

# Sources to process. Must match entries in keys/target_definitions.txt.
# Set to None to process all sources defined in keys/imaging_key.txt.
TARGETS = None  # e.g. ["ngc5194", "ngc5457"]

# ---------------------------------------------------------------------------
# RUN — no need to edit below this line
# ---------------------------------------------------------------------------

handler = pys.PipelineHandler(key_dir=KEY_DIR)

if STAGES is None:
    handler.run_all(targets=TARGETS)
else:
    handler.run_stages(stages=STAGES, targets=TARGETS)
