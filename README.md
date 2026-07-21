<!-- back-to-top anchor -->
<a name="readme-top"></a>

<br />
<div align="center">
  <a href="https://github.com/lukas-neumann-astro/astro-HexMaps">
    <img src="images/logo.png" alt="Logo" width="100" height="100">
  </a>

  <h3 align="center">HexMaps</h3>

  <p align="center">
    Hexagonal-grid Multi-data Analysis and Processing Software
    <br />
    <a href="https://astro-hexmaps.readthedocs.io/en/latest/"><strong>Explore the docs »</strong></a>
    <br /><br />
    <a href="https://astro-hexmaps.readthedocs.io/en/latest/quickstart.html">View Demo</a>
    ·
    <a href="https://github.com/lukas-neumann-astro/astro-HexMaps/issues">Report Bug</a>
    ·
    <a href="https://github.com/lukas-neumann-astro/astro-HexMaps/issues">Request Feature</a>
  </p>
</div>

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]

---

## About the Project

[![HexMaps screenshot][product-screenshot]](https://astro-hexmaps.readthedocs.io/en/latest/)

**HexMaps** is a Python package for homogenizing and analysing
multi-wavelength astronomical datasets on hexagonal grids. It is ideally
suited for combining a heterogeneous set of 2D maps and 3D spectral cubes
observed at different angular resolutions and pixel scales into a single,
science-ready table for sightline-by-sightline analysis.

Given a set of input data, HexMaps:

1. **Regrids** all maps and cubes onto a common hexagonal sampling grid at a
   user-specified angular resolution, convolving each dataset to a common
   beam.
2. **Processes** each spectral cube: builds an S/N mask from a reference line,
   then computes moment maps (integrated intensity, mean velocity, line width,
   peak temperature, rms, equivalent width) and shuffled spectra for every
   line.
3. **Writes FITS** output: PPV-native moment maps, optionally convolved cubes,
   and 2D band images — all on the same pixel grid as the overlay cube.

The primary deliverable is an Astropy `.ecsv` table with one row per
hexagonal sightline, containing all spectra, moment maps, and 2D map values
side by side. This makes line-ratio analysis, radial profile extraction, and
spectral stacking straightforward with standard Python tools.

HexMaps is the successor to
[PyStructure (PhangsTeam)](https://github.com/PhangsTeam/PyStructure),
refactored into a pip-installable package with a clean CLI, an INI-style
single configuration file, and a modular stage architecture.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Getting Started

### Prerequisites

HexMaps requires **Python ≥ 3.10**. All Python dependencies are installed
automatically by pip:

```
astropy  numpy  pandas  scipy  matplotlib
reproject  radio_beam  spectral_cube  scikit-image
```

### Installation

```bash
# From PyPI (once published)
pip install astro-hexmaps

# From GitHub — latest version
pip install git+https://github.com/lukas-neumann-astro/astro-HexMaps.git

# Editable / development install
git clone https://github.com/lukas-neumann-astro/astro-HexMaps.git
cd astro-HexMaps
pip install -e ".[dev]"
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Usage

### 1 — Initialise a working directory

```bash
# Creates config.txt, keys/, and run_hexmaps.py in the current folder
hexmaps --init

# Or target a specific directory
hexmaps --init --workdir ~/my_project
cd ~/my_project
```

This copies template configuration files into your working directory.
The installed package is never modified.

### 2 — Download the example dataset *(optional)*

To try the pipeline immediately with real data, download the bundled
NGC 5194 example dataset (~46 MB):

```bash
hexmaps --download-example --workdir ~/my_project
```

This fetches the input FITS files (CO(2–1) cube, CO(1–0) cube, SPIRE 250 µm
map and associated files) into `~/my_project/data/`. The bundled `config.txt`
is already configured to use these files. Use `--force` to re-download
existing files.

You can also download the example Jupyter notebook separately:

```bash
hexmaps --download-notebook --workdir ~/my_project
```

### 3 — Edit your configuration

| File | What to configure | How often |
|------|-------------------|-----------|
| `config.txt` | data directory, target list, overlay cube, maps/cubes, target resolution, masking, output flags | every run |
| `keys/target_definitions.txt` | sky coordinates, distance, inclination per target | only when adding new targets |
| `keys/hfs_lines.txt` *(optional)* | Hyperfine structure line definitions | rarely |

`config.txt` replaces the `PyStructure.conf` file of the old PyStructure.

> **Migrating from PyStructure?** Use the conversion scripts:
>
> ```bash
> python conversion_from_pystructure/config_conversion.py PyStructure.conf config.txt
> python conversion_from_pystructure/target_definitions_conversion.py geometry.txt keys/target_definitions.txt
> python conversion_from_pystructure/hfs_lines_conversion.py hfs_lines.txt keys/hfs_lines.txt
> ```

### 4 — Run

```bash
# Default: regrid + products (writes .ecsv database)
hexmaps --conf config.txt

# All stages including FITS output
hexmaps --conf config.txt --stages all

# Single target
hexmaps --conf config.txt --targets ngc5194

# Also write a log file
hexmaps --conf config.txt --log_file run.log
```

Or from Python:

```python
import hexmaps as hm

handler = hm.PipelineHandler(conf_path="config.txt")
handler.run_all()                                      # regrid + products
handler.run_stages(["regrid", "products", "fits"])     # include FITS output
handler.run_stages(["fits"], targets=["ngc5194"])      # re-run one stage only
```

*For more examples, please refer to the [Documentation](https://astro-hexmaps.readthedocs.io/en/latest/).*

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Repository Layout

```
astro-HexMaps/                       ← git repository root (pip install this)
├── hexmaps/                         ← installable package
│   ├── handler_keys.py              reads & validates config and key files
│   ├── handler_targets.py           target geometry lookups
│   ├── handler_pipeline.py          PipelineHandler: stage orchestration
│   ├── stage_regrid.py              hex grid + convolution + sampling → .ecsv
│   ├── stage_products.py            spectral masking, moments, shuffled spectra
│   ├── stage_fits.py                FITS moment maps / cubes / band images
│   ├── utils_fits.py                FITS/WCS helpers (convolution, reprojection)
│   ├── utils_table.py               table I/O, spectral shuffle, moments
│   ├── hexmaps_analysis.py          HexMapsAnalysis class (importable as part of package)
│   ├── logger.py                    centralised stage-labelled logger
│   ├── init_workdir.py              --init scaffolding
│   ├── download_example.py          --download-example and --download-notebook
│   ├── cli.py                       hexmaps console-script entry point
│   ├── test_hexmaps.py              unit and integration tests
│   └── templates/                   template files copied by --init
│       ├── config.txt
│       ├── run_hexmaps.py
│       └── keys/
│           ├── target_definitions.txt
│           └── hfs_lines.txt
├── config.txt                       ← example / template config file
├── keys/
│   ├── target_definitions.txt       ← target geometry table (PHANGS example)
│   └── hfs_lines.txt                ← hyperfine structure definitions
├── analysis/
│   └── hexmaps_example.ipynb        example analysis notebook
├── conversion_from_pystructure/     ← migration scripts from old PyStructure
│   ├── config_conversion.py
│   ├── target_definitions_conversion.py
│   └── hfs_lines_conversion.py
├── data/                            ← example FITS input (NGC 5194)
├── docs/                            ← Sphinx / Read the Docs source
├── images/                          ← README images (logo, screenshot)
├── run_hexmaps.py                   ← example run script
└── pyproject.toml
```

---

## Pipeline Stages

The pipeline runs three stages always in this order:

| Stage | Module | Default | Description |
|-------|--------|---------|-------------|
| `regrid` | `stage_regrid.py` | ✓ | Generate the hexagonal sampling grid; convolve and sample all maps & cubes onto it; write the `.ecsv` database |
| `products` | `stage_products.py` | ✓ | Build the S/N mask; compute moment maps (mom0/1/2, Tpeak, rms, EW) and shuffled spectra |
| `fits` | `stage_fits.py` | optional | Compute PPV-native moment maps on the convolved cubes; write FITS images |

The default run (`hexmaps --conf config.txt`) executes **regrid + products** only.
Add `--stages all` to also run the fits stage.

---

## Reading the Output

```python
from hexmaps.utils_table import load_hexmaps
import matplotlib.pyplot as plt

table = load_hexmaps("output/ngc5194_hexmaps_27p0as_2025_01_01.ecsv")

plt.figure(figsize=(5, 5))
plt.scatter(table["RA"], table["DEC"],
            c=table["MOM0_12CO21"], marker="h", s=60, cmap="inferno")
plt.gca().invert_xaxis()
plt.xlabel("R.A. [deg]"); plt.ylabel("Dec. [deg]")
plt.show()
```

For richer quicklook plots use the `HexMapsAnalysis` class, which is now
part of the installed package:

```python
from hexmaps import HexMapsAnalysis

db = HexMapsAnalysis("output/ngc5194_hexmaps_27p0as_2025_01_01.ecsv")
db.quickplot_map("12CO21")
db.quickplot_spectrum("12CO21")

# Recover provenance information embedded at run time
print(db.get_config())               # recover config.txt used for this run
print(db.get_log())                  # recover the full pipeline log
hdr = db.get_input_header("12CO21")  # recover raw FITS header of input cube
print(db.list_input_headers())       # list all embedded headers
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Roadmap

- [ ] PyPI release
- [ ] Expanded documentation and tutorials
- [ ] Additional analysis utilities in `hexmaps_analysis.py`
- [ ] Various feature updates (e.g. "island-method" masking, chunking for large datasets)

See the [open issues](https://github.com/lukas-neumann-astro/astro-HexMaps/issues)
for a full list of proposed features and known issues.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Contributing

Contributions are greatly appreciated! If you have a suggestion that would
make this better, please fork the repository and create a pull request. You
can also open an issue with the tag "enhancement".

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## License

Distributed under the MIT License — see [LICENSE](LICENSE) for details.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Contact

Dr. Jakob den Brok — jadenbrok@mpia.de

Dr. Lukas Neumann — lukas.neumann@eso.org

Project Link: [https://github.com/lukas-neumann-astro/astro-HexMaps](https://github.com/lukas-neumann-astro/astro-HexMaps)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

---

## Acknowledgements

HexMaps builds on the original PyStructure IDL scripts developed within the
PHANGS collaboration. The routines have been updated, improved, and fully
rewritten in Python.

The code has been employed in several peer-reviewed publications including
den Brok et al. (2021, 2022, 2023, 2025), Eibensteiner et al. (2022, 2023),
Neumann et al. (2023), Stuber et al. (2025), and others. See the
[documentation](https://astro-hexmaps.readthedocs.io) for the full list.

* [PHANGS collaboration](https://sites.google.com/view/phangs/home)
* [Prof. Bigiel's research group](https://www.astro.uni-bonn.de/~bigiel/)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- MARKDOWN LINKS & IMAGES -->
[contributors-shield]: https://img.shields.io/github/contributors/lukas-neumann-astro/astro-HexMaps.svg?style=for-the-badge
[contributors-url]: https://github.com/lukas-neumann-astro/astro-HexMaps/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/lukas-neumann-astro/astro-HexMaps.svg?style=for-the-badge
[forks-url]: https://github.com/lukas-neumann-astro/astro-HexMaps/network/members
[stars-shield]: https://img.shields.io/github/stars/lukas-neumann-astro/astro-HexMaps.svg?style=for-the-badge
[stars-url]: https://github.com/lukas-neumann-astro/astro-HexMaps/stargazers
[issues-shield]: https://img.shields.io/github/issues/lukas-neumann-astro/astro-HexMaps.svg?style=for-the-badge
[issues-url]: https://github.com/lukas-neumann-astro/astro-HexMaps/issues
[license-shield]: https://img.shields.io/github/license/lukas-neumann-astro/astro-HexMaps.svg?style=for-the-badge
[license-url]: https://github.com/lukas-neumann-astro/astro-HexMaps/blob/master/LICENSE
[product-screenshot]: images/screenshot.png
