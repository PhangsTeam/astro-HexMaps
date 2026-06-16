# PyStructure

[![Contributors](https://img.shields.io/github/contributors/jdenbrok/PyStructure.svg?style=for-the-badge)](https://github.com/jdenbrok/PyStructure/graphs/contributors)
[![MIT License](https://img.shields.io/github/license/jdenbrok/PyStructure.svg?style=for-the-badge)](LICENSE)

A Python package for homogenizing and analyzing multi-wavelength astronomical
datasets on hexagonal grids. Samples 2D images (maps) and 3D spectral cubes
at a common resolution and grid, producing Astropy Table output (`.ecsv`)
along with optional FITS moment and map images.

---

## Installation

```bash
# From PyPI (once published)
pip install pystructurePipeline

# From GitHub (latest)
pip install git+https://github.com/lukas-neumann-astro/PyStructure.git

# Editable / development install (from inside the cloned repo)
git clone https://github.com/lukas-neumann-astro/PyStructure.git
cd PyStructure
pip install -e ".[dev]"
```

The installed package lives in your Python environment (e.g. your conda env).
Your data, key files, and outputs live separately in a **working directory**
that you control.

---

## Quick start

### 1 — Set up a working directory

```bash
# Creates keys/ and run_pystructure.py in the current folder
pystructure --init

# Or choose a destination
pystructure --init --workdir ~/projects/my_galaxy_survey
cd ~/projects/my_galaxy_survey
```

This copies the editable key files and a ready-to-run script into your
working directory. The installed package is never modified.

### 2 — Edit the key files

| File | What to edit |
|------|-------------|
| `keys/master_key.txt` | `data_dir`, `out_dir`, your name |
| `keys/target_definitions.txt` | RA, Dec, distance, inclination for each source |
| `keys/data_key.txt` | Source list, overlay cube, map and cube file definitions |
| `keys/config_key.txt` | Target resolution, masking thresholds, output flags |
| `keys/hfs_lines.txt` *(optional)* | Hyperfine structure line definitions |

### 3 — Run

```bash
# Edit and run the script in your working directory
python run_pystructure.py

# Or use the CLI directly
pystructure --key_dir keys/

# Specific stages only
pystructure --key_dir keys/ --stages regrid products

# Single source
pystructure --key_dir keys/ --targets ngc5194

# Write a log file in addition to stdout
pystructure --key_dir keys/ --log_file pystructure_run.log
```

### 4 — Use from Python

```python
import pystructurePipeline as pys

handler = pys.PipelineHandler(key_dir="keys/")
handler.run_all()

# Or selectively
handler.run_stages(["regrid", "products"], targets=["ngc5194"])
```

---

## Repository layout

```
PyStructure/                      <- git repo - install this with pip
|-- pystructurePipeline/          <- installable package
|   |-- handler_keys.py               reads & validates all key files
|   |-- handler_sources.py            source geometry lookups
|   |-- handler_pipeline.py           PipelineHandler: stage orchestration
|   |-- stage_regrid.py               hex grid generation + convolution + sampling
|   |-- stage_products.py             spectral masking, moments, shuffled spectra
|   |-- stage_fits.py                 FITS moment-map / map-image writing
|   |-- utils_fits.py                 FITS/WCS helpers (convolution, deprojection, ...)
|   |-- utils_table.py                table I/O, spectral shuffle, moment computation
|   |-- pystructureLogger.py          centralized [pyStructure] [Stage] [LEVEL] logger
|   |-- init_workdir.py               copies key-file templates (--init)
|   |-- cli.py                        `pystructure` console-script entry point
|   `-- test_pystructure.py
|-- keys/                          <- example / template key files
|   |-- master_key.txt
|   |-- target_definitions.txt
|   |-- data_key.txt
|   |-- config_key.txt
|   `-- hfs_lines.txt
|-- analysis/                      <- post-processing & plotting helpers
|   |-- pystructureAnalysis.py        load .ecsv, quicklook maps/spectra
|   `-- pystructure_example.ipynb     example analysis notebook
|-- data/                          <- example FITS input (NGC 5194)
|-- run_pystructure.py             <- example run script
|-- pyproject.toml
`-- README.md

~/my_project/                      <- your working directory (anywhere on disk)
|-- keys/
|   |-- master_key.txt              <- edit these
|   |-- target_definitions.txt
|   |-- data_key.txt
|   |-- config_key.txt
|   `-- hfs_lines.txt
|-- data/                           <- your FITS files
|-- output/                         <- pipeline writes .ecsv tables here
|-- saved_FITS_files/                  FITS moment/map images land here
`-- run_pystructure.py              <- edit and run this
```

---

## Pipeline stages

The pipeline is organized into three stages, always executed in this order
regardless of the order you list them in:

| Stage | Module | Description |
|-------|--------|-------------|
| `regrid` | `stage_regrid.py` | Generate the hexagonal sampling grid from the overlay cube, then convolve and sample maps & cubes onto it; write the PyStructure `.ecsv` |
| `products` | `stage_products.py` | Build the S/N mask, compute moment maps (mom0/1/2, Tpeak, rms, EW), and shuffled spectra for every line |
| `fits` | `stage_fits.py` | Regrid the moment maps and 2D maps back onto a rectangular pixel grid and write FITS files |

Run a subset of stages with `--stages` (CLI) or `run_stages([...])` (Python).
The hexagonal sampling grid generation, previously a separate `sampling`
stage, is now an internal step of `regrid`.

---

## Logging

All pipeline output goes through a centralized logger
(`pystructurePipeline.pystructureLogger`), giving every message a consistent,
column-aligned format:

```
[pyStructure] [<Stage>]   [<LEVEL>]   <message>
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
[pyStructure] [Loading]   [INFO]     Loading key files...
[pyStructure] [Loading]   [INFO]     Loaded 1 source(s): ['ngc5194']
[pyStructure] [Regrid]    [INFO]     Hexagonal grid generated: 1060 sampling points (spacing = 13.5 arcsec).
[pyStructure] [Regrid]    [INFO]     Cube 12co21 sampled successfully.
[pyStructure] [Products]  [INFO]     Mask complete. Computing moments.
[pyStructure] [FITS]      [INFO]     Moment map FITS files written to: ./saved_FITS_files/
[pyStructure] [Return]    [INFO]     --- Run summary ---
```

Pass `--log_file run.log` (CLI) or `PipelineHandler(..., log_file="run.log")`
(Python) to additionally stream every message to a file. The full log history
can also be written at any time with `handler.save_log("run.log")`, or
inspected programmatically via `pystructureLogger.logger.get_records(...)`.

---

## Reading the output

```python
from pystructurePipeline.utils_table import load_pystructure

table = load_pystructure("output/ngc5194_data_struct_27as_2025_01_01.ecsv")

import numpy as np, matplotlib.pyplot as plt
mom0 = table["MOM0_12CO21"]
plt.scatter(table["ra_deg"], table["dec_deg"], c=mom0, marker="h")
plt.show()
```

For richer quicklook plots (maps, spectra, shuffled spectra, radial
profiles), see `analysis/pystructureAnalysis.py` and the accompanying
`analysis/pystructure_example.ipynb` notebook.

---

## License

Distributed under the MIT License — see [LICENSE](LICENSE) for details.

## Contact

Dr. Jakob den Brok — jadenbrok@mpia.de
Dr. Lukas Neumann — lukas.neumann@eso.org
