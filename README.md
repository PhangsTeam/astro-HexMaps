# HexMaps

[![Contributors](https://img.shields.io/github/contributors/jdenbrok/HexMaps.svg?style=for-the-badge)](https://github.com/jdenbrok/HexMaps/graphs/contributors)
[![MIT License](https://img.shields.io/github/license/jdenbrok/HexMaps.svg?style=for-the-badge)](LICENSE)

A Python package for homogenizing and analyzing multi-wavelength astronomical
datasets on hexagonal grids. Samples 2D images (maps) and 3D spectral cubes
at a common resolution and grid, producing Astropy Table output (`.ecsv`)
along with optional FITS moment and map images.

---

## Installation

```bash
# From PyPI (once published)
pip install hexmaps

# From GitHub (latest)
pip install git+https://github.com/lukas-neumann-astro/HexMaps.git

# Editable / development install (from inside the cloned repo)
git clone https://github.com/lukas-neumann-astro/HexMaps.git
cd HexMaps
pip install -e ".[dev]"
```

The installed package lives in your Python environment (e.g. your conda env).
Your data, key files, and outputs live separately in a **working directory**
that you control.

---

## Quick start

### 1 — Set up a working directory

```bash
# Creates config.txt, keys/, and run_hexmaps.py in the current folder
hexmaps --init

# Or choose a destination
hexmaps --init --workdir ~/projects/my_galaxy_survey
cd ~/projects/my_galaxy_survey
```

This copies a config file, a `keys/` subfolder, and a ready-to-run script
into your working directory. The installed package is never modified.

### 2 — Edit your configuration

Configuration is split into two parts, based on how often each one changes:

| File | What to edit | How often |
|------|-------------|-----------|
| `config.txt` | `data_dir`/`out_dir`, your name, source list, overlay cube, map/cube file definitions, target resolution, masking thresholds, output flags | every run |
| `keys/target_definitions.txt` | RA, Dec, distance, inclination for each source | rarely — set up once, reuse across projects |
| `keys/hfs_lines.txt` *(optional)* | Hyperfine structure line definitions | rarely |

`config.txt` is the single file you'll typically touch: it combines what
used to be three separate files (`master_key.txt`, `data_key.txt`,
`config_key.txt`) into one. `keys/target_definitions.txt` and
`keys/hfs_lines.txt` live in a fixed `keys/` subfolder next to `config.txt`
and are usually shared across many projects, so they're kept separate.

### 3 — Run

```bash
# Edit and run the script in your working directory
python run_hexmaps.py

# Or use the CLI directly (runs regrid + products by default)
hexmaps --conf config.txt

# Include the optional fits stage (FITS moment maps and band images)
hexmaps --conf config.txt --stages regrid products fits

# Single source
hexmaps --conf config.txt --targets ngc5194

# Write a log file in addition to stdout
hexmaps --conf config.txt --log_file hexmaps_run.log
```

### 4 — Use from Python

```python
import hexmaps as hm

handler = hm.PipelineHandler(conf_path="config.txt")
handler.run_all()  # default: regrid + products only

# Include the optional fits stage
handler.run_stages(["regrid", "products", "fits"])

# Or a specific subset
handler.run_stages(["regrid", "products"], targets=["ngc5194"])
```

> **Migrating from an older version?** TBD

---

## Repository layout

```
HexMaps/                            <- git repo - install this with pip
|-- hexmaps/                        <- installable package
|   |-- handler_keys.py             reads & validates all key files
|   |-- handler_sources.py          source geometry lookups
|   |-- handler_pipeline.py         PipelineHandler: stage orchestration
|   |-- stage_regrid.py             hex grid generation + convolution + sampling
|   |-- stage_products.py           spectral masking, moments, shuffled spectra
|   |-- stage_fits.py               FITS moment-map / map-image writing
|   |-- utils_fits.py               FITS/WCS helpers (convolution, deprojection, ...)
|   |-- utils_table.py              table I/O, spectral shuffle, moment computation
|   |-- hexmapsLogger.py            centralized [HexMaps] [Stage] [LEVEL] logger
|   |-- init_workdir.py             copies key-file templates (--init)
|   |-- cli.py                      `hexmaps` console-script entry point
|   `-- test_hexmaps.py
|-- config.txt                      <- example / template config file
|-- keys/                           <- example / template reference tables
|   |-- target_definitions.txt
|   `-- hfs_lines.txt
|-- analysis/                       <- post-processing & plotting helpers
|   |-- hexmapsAnalysis.py          load .ecsv, quicklook maps/spectra
|   `-- hexmaps_example.ipynb       example analysis notebook
|-- data/                           <- example FITS input (NGC 5194)
|-- run_hexmaps.py                  <- example run script
|-- pyproject.toml
`-- README.md

~/my_project/                       <- your working directory (anywhere on disk)
|-- config.txt                      <- edit this every run
|-- keys/
|   |-- target_definitions.txt      <- edit once, reuse across projects
|   `-- hfs_lines.txt
|-- data/                           <- your FITS files
|-- output/                         <- pipeline writes .ecsv tables here
|-- saved_FITS_files/                  FITS moment/map images land here
`-- run_hexmaps.py                  <- edit and run this
```

---

## Pipeline stages

The pipeline is organized into three stages, always executed in this order
regardless of the order you list them in:

| Stage | Module | Default | Description |
|-------|--------|---------|-------------|
| `regrid` | `stage_regrid.py` | ✓ | Generate the hexagonal sampling grid from the overlay cube, then convolve and sample maps & cubes onto it; write the HexMaps `.ecsv` |
| `products` | `stage_products.py` | ✓ | Build the S/N mask, compute moment maps (mom0/1/2, Tpeak, rms, EW), and shuffled spectra for every line |
| `fits` | `stage_fits.py` | — optional | Compute PPV-native moment maps directly on the convolved cubes and write FITS files (moment maps, band images, mask cubes) |

The default run (`run_all()` / `hexmaps --conf config.txt`) executes only
**regrid** and **products** — the primary pipeline deliverable is the `.ecsv`
database. The `fits` stage is an optional bonus that produces convenient FITS
images; enable it explicitly with `--stages regrid products fits` or
`run_stages(["regrid", "products", "fits"])`.

The hexagonal sampling grid generation, previously a separate `sampling`
stage, is now an internal step of `regrid`.

---

## Logging

All pipeline output goes through a centralized logger
(`hexmaps.hexmapsLogger`), giving every message a consistent,
column-aligned format:

```
YYYY-MM-DD HH:MM:SS[<Stage>][<LEVEL>] <message>
```

Stages used during a run:

| Stage label | When it's used |
|-------------|-----------------|
| `Loading`   | Reading key files, validating configuration, per-source setup |
| `Regrid`    | Hex grid generation, convolution, reprojection, sampling |
| `Products`  | Mask construction, moment maps, shuffled spectra |
| `FITS`      | Writing FITS moment maps and 2D map images |
| `Return`    | Per-source error reporting and the final run summary |

Example:

```
YYYY-MM-DD HH:MM:SS [Loading]  [INFO]   Loading key files...
YYYY-MM-DD HH:MM:SS [Loading]  [INFO]   Loaded 1 source(s): ['ngc5194']
YYYY-MM-DD HH:MM:SS [Regrid]   [INFO]   Hexagonal grid generated: 1060 sampling points (spacing = 13.5 arcsec).
YYYY-MM-DD HH:MM:SS [Regrid]   [INFO]   Cube 12co21 sampled successfully.
YYYY-MM-DD HH:MM:SS [Products] [INFO]   Mask complete. Computing moments.
YYYY-MM-DD HH:MM:SS [FITS]     [INFO]   Moment map FITS files written to: ./saved_FITS_files/
YYYY-MM-DD HH:MM:SS [Return]   [INFO]   --- Run summary ---
```

Pass `--log_file run.log` (CLI) or `PipelineHandler(..., log_file="run.log")`
(Python) to additionally stream every message to a file. The full log history
can also be written at any time with `handler.save_log("run.log")`, or
inspected programmatically via `hexmapsLogger.logger.get_records(...)`.

---

## Reading the output

```python
from hexmaps.utils_table import load_hexmaps

table = load_hexmaps("output/ngc5194_data_struct_27as_2025_01_01.ecsv")

import numpy as np, matplotlib.pyplot as plt
mom0 = table["MOM0_12CO21"]
plt.scatter(table["ra_deg"], table["dec_deg"], c=mom0, marker="h")
plt.show()
```

For richer quicklook plots (maps, spectra, shuffled spectra, radial
profiles), see `analysis/hexmapsAnalysis.py` and the accompanying
`analysis/hexmaps_example.ipynb` notebook.

---

## License

Distributed under the MIT License — see [LICENSE](LICENSE) for details.

## Contact

* Dr. Jakob den Brok — jadenbrok@mpia.de
* Dr. Lukas Neumann — lukas.neumann@eso.org