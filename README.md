# PyStructure

[![Contributors](https://img.shields.io/github/contributors/jdenbrok/PyStructure.svg?style=for-the-badge)](https://github.com/jdenbrok/PyStructure/graphs/contributors)
[![MIT License](https://img.shields.io/github/license/jdenbrok/PyStructure.svg?style=for-the-badge)](LICENSE)

A Python package for homogenizing and analyzing multi-wavelength astronomical
datasets on hexagonal grids.  Samples 2D images (bands) and 3D spectral cubes
at a common resolution and grid, producing Astropy Table output (`.ecsv`) along
with optional FITS moment and band maps.

---

## Installation

```bash
# From PyPI (once published)
pip install pystructure

# From GitHub (latest)
pip install git+https://github.com/PhangsTeam/PyStructure.git

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

This copies four editable key files and a ready-to-run script into your
working directory. The installed package is never modified.

### 2 — Edit the key files

| File | What to edit |
|------|-------------|
| `keys/master_key.txt` | `data_dir`, `out_dir`, your name |
| `keys/target_definitions.txt` | RA, Dec, distance, inclination for each galaxy |
| `keys/imaging_key.txt` | Source list, overlay image, band and cube file definitions |
| `keys/config_key.txt` | Target resolution, masking thresholds, output flags |

### 3 — Run

```bash
# Edit and run the script in your working directory
python run_pystructure.py

# Or use the CLI directly
pystructure --key_dir keys/

# Specific stages only
pystructure --key_dir keys/ --stages sampling regrid

# Single source
pystructure --key_dir keys/ --targets ngc5194
```

### 4 — Use from Python

```python
import pystructure as pys

handler = pys.PipelineHandler(key_dir="keys/")
handler.run_all()

# Or selectively
handler.run_stages(["regrid", "spectra"], targets=["ngc5194"])
```

---

## Repository layout

```
PyStructure/                 ← git repo — install this with pip
├── pystructure/             ← installable package
│   ├── handlers/
│   │   ├── key_handler.py       reads all key files
│   │   ├── target_handler.py    galaxy geometry lookups
│   │   └── pipeline_handler.py  stage orchestration
│   ├── stages/
│   │   ├── stage_sampling.py    hexagonal grid generation
│   │   ├── stage_regrid.py      convolution and sampling
│   │   ├── stage_spectra.py     spectral processing and moments
│   │   └── stage_output.py      FITS map writing
│   ├── utils/
│   │   ├── fits_utils.py
│   │   └── table_utils.py
│   ├── templates/           ← bundled templates (copied by --init)
│   │   ├── keys/
│   │   │   ├── master_key.txt
│   │   │   ├── target_definitions.txt
│   │   │   ├── imaging_key.txt
│   │   │   └── config_key.txt
│   │   └── run_pystructure.py
│   ├── init_workdir.py
│   └── cli.py
├── tests/
├── pyproject.toml
└── README.md

~/my_project/                ← your working directory (anywhere on disk)
├── keys/
│   ├── master_key.txt       ← edit these
│   ├── target_definitions.txt
│   ├── imaging_key.txt
│   └── config_key.txt
├── data/                    ← your FITS files
├── Output/                  ← pipeline writes .ecsv tables here
├── saved_FITS_files/        ← FITS moment/band maps land here
└── run_pystructure.py       ← edit and run this
```

---

## Pipeline stages

| Stage | Description |
|-------|-------------|
| `sampling` | Generate hexagonal sampling grid from the overlay cube |
| `regrid` | Convolve and sample bands & cubes; write the PyStructure `.ecsv` |
| `spectra` | Mask spectra, compute moments (mom0/1/2, EW), shuffle |
| `output` | Write FITS moment maps and band maps |

---

## Reading the output

```python
from pystructure.utils import load_pystructure

table = load_pystructure("Output/ngc5194_data_struct_27as_2025_01_01.ecsv")

import numpy as np, matplotlib.pyplot as plt
mom0 = np.nansum(table["SPEC_12CO21"], axis=1)
plt.scatter(table["ra_deg"], table["dec_deg"], c=mom0, marker="h")
plt.show()
```

---

## License

Distributed under the MIT License — see [LICENSE](LICENSE) for details.

## Contact

Dr. Jakob den Brok — jadenbrok@mpia.de  
Dr. Lukas Neumann — lukas.neumann@eso.org
