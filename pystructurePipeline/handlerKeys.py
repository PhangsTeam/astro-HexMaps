"""
KeyHandler: reads and validates all PyStructure key files.

Replaces the monolithic PyStructure.conf parser. The configuration is split
into separate key files by concern:

  master_key.txt         - paths, global settings, run metadata
  target_definitions.txt - list of sources (replaces geom_file reference)
  data_key.txt        - map and cube definitions
  config_key.txt         - resolution, masking, output settings
"""

import os
import configparser
import pandas as pd
from pathlib import Path


# ---------------------------------------------------------------------------
# Column definitions for the tabular key files
# ---------------------------------------------------------------------------

MAP_COLUMNS = ["map_name", "map_desc", "map_unit", "map_ext", "map_dir", "map_uc"]
CUBE_COLUMNS = ["line_name", "line_desc", "line_unit", "line_ext", "line_dir", "map_ext", "map_uc"]
MASK_COLUMNS_VEL  = ["mask_name", "mask_desc", "mask_start", "mask_end", "mask_unit"]
MASK_COLUMNS_FILE = ["mask_name", "mask_desc", "mask_ext", "mask_dir"]
TARGET_COLUMNS = [
    "source", "ra_ctr", "dec_ctr", "dist_mpc", "e_dist_mpc",
    "incl_deg", "e_incl_deg", "posang_deg", "e_posang_deg",
    "r25", "e_r25",
]
HFS_COLUMNS = ["hfs_name", "hfs_ref_freq", "hfs_freq", "unit"]


class KeyHandler:
    """
    Reads and validates all PyStructure key files from a key directory.

    Parameters
    ----------
    key_dir : str or Path
        Directory containing the key files.

    Example
    -------
    >>> kh = KeyHandler("./keys/")
    >>> print(kh.sources)
    >>> print(kh.maps)
    """

    def __init__(self, key_dir: str):
        self.key_dir = Path(key_dir)
        self._validate_key_dir()

        # Parsed data (populated by load())
        self.meta = {}          # global settings from master_key.txt / config_key.txt
        self.sources = []       # list of source name strings
        self.source_table = None  # full geometry DataFrame
        self.maps = None       # DataFrame of map (2D image) definitions
        self.cubes = None       # DataFrame of cube definitions
        self.input_mask = None  # DataFrame of mask definitions
        self.hfs_data = None    # DataFrame of hyperfine structure lines (optional)

        self.load()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def load(self):
        """Load and validate all key files."""
        self._load_master_key()
        self._load_config_key()
        self._load_target_definitions()
        self._load_data_key()
        self._load_hfs_key()

    def get_sources(self) -> list:
        """Return the list of source names to be processed."""
        return list(self.sources)

    def get_maps(self) -> pd.DataFrame:
        return self.maps

    def get_cubes(self) -> pd.DataFrame:
        return self.cubes

    def get_input_mask(self) -> pd.DataFrame:
        return self.input_mask

    def get_source_table(self) -> pd.DataFrame:
        return self.source_table

    def get_hfs_data(self):
        return self.hfs_data

    # ------------------------------------------------------------------
    # Private loaders
    # ------------------------------------------------------------------

    def _validate_key_dir(self):
        if not self.key_dir.is_dir():
            raise FileNotFoundError(f"Key directory not found: {self.key_dir}")

    def _load_master_key(self):
        """
        Load master_key.txt.

        Expected ini-style file with a [paths] section and a [meta] section.
        Example:

            [paths]
            data_dir    = data/
            out_dir     = Output/
            geom_file   = keys/target_definitions.txt
            data_key = keys/data_key.txt
            config_key  = keys/config_key.txt

            [meta]
            user     = Dr. Hubble Trouble
            comments = Example run
        """
        master_path = self.key_dir / "master_key.txt"
        if not master_path.exists():
            raise FileNotFoundError(f"master_key.txt not found in {self.key_dir}")

        cfg = configparser.ConfigParser(inline_comment_prefixes=("#",))
        cfg.read(master_path)

        paths = dict(cfg["paths"]) if "paths" in cfg else {}
        meta  = dict(cfg["meta"])  if "meta"  in cfg else {}

        # Store paths — resolve relative to key_dir parent so the package
        # works regardless of the working directory.
        base = self.key_dir.resolve().parent
        self.meta["data_dir"]    = str(base / paths.get("data_dir",    "data/"))
        self.meta["out_dir"]     = str(base / paths.get("out_dir",     "Output/"))
        self.meta["geom_file"]   = str(base / paths.get("geom_file",   "keys/target_definitions.txt"))
        self.meta["data_key"] = str(base / paths.get("data_key", "keys/data_key.txt"))
        self.meta["config_key"]  = str(base / paths.get("config_key",  "keys/config_key.txt"))
        self.meta["hfs_file"]    = str(base / paths.get("hfs_file",    "")) if paths.get("hfs_file") else None

        self.meta["user"]     = meta.get("user",     "Unknown user")
        self.meta["comments"] = meta.get("comments", "")
        self.meta["_base"]    = str(base)  # absolute repo root for resolving relative paths

    def _load_config_key(self):
        """
        Load config_key.txt — resolution, masking, output flags, spectral settings.

        Expected ini-style file with a [resolution], [masking], [output],
        [spectral] and [structure] section.
        """
        config_path = Path(self.meta["config_key"])
        if not config_path.exists():
            raise FileNotFoundError(f"config_key not found: {config_path}")

        cfg = configparser.ConfigParser(inline_comment_prefixes=("#",))
        cfg.read(config_path)

        def _get(section, key, fallback):
            return cfg.get(section, key, fallback=str(fallback))

        # Resolution settings
        self.meta["target_res"]       = float(_get("resolution", "target_res",       27.0))
        self.meta["resolution"]       = _get("resolution", "resolution",       "angular")
        self.meta["spacing_per_beam"] = float(_get("resolution", "spacing_per_beam", 2.0))
        self.meta["max_rad"]          = _get("resolution", "max_rad",          "auto")
        self.meta["NAXIS_shuff"]      = int(float(_get("resolution", "NAXIS_shuff",  200)))
        self.meta["CDELT_SHUFF"]      = float(_get("resolution", "CDELT_SHUFF",      4000.0))

        # Masking settings
        self.meta["ref_line"]          = _get("masking", "ref_line",          "first")
        self.meta["SN_processing"]     = [
            float(x) for x in _get("masking", "SN_processing", "2,4").split(",")
        ]
        self.meta["strict_mask"]       = _get("masking", "strict_mask",       "false").lower() == "true"
        self.meta["use_input_mask"]    = _get("masking", "use_input_mask",    "false").lower() == "true"
        self.meta["use_fixed_vel_mask"]= _get("masking", "use_fixed_vel_mask","false").lower() == "true"
        self.meta["use_hfs_lines"]     = _get("masking", "use_hfs_lines",     "false").lower() == "true"
        self.meta["mom_thresh"]        = float(_get("masking", "mom_thresh",        5.0))
        self.meta["conseq_channels"]   = int(float(_get("masking", "conseq_channels",   3)))
        self.meta["mom2_method"]       = _get("masking", "mom2_method",       "fwhm")

        # Output settings
        self.meta["save_fits"]         = _get("output", "save_fits",          "false").lower() == "true"
        self.meta["save_mom_maps"]     = _get("output", "save_mom_maps",      "true").lower()  == "true"
        self.meta["save_maps"]    = _get("output", "save_maps",     "true").lower()  == "true"
        self.meta["folder_savefits"]   = _get("output", "folder_savefits",    "./saved_FITS_files/")

        # Spectral smoothing
        self.meta["spec_smooth"]        = _get("spectral", "spec_smooth",        "default")
        self.meta["spec_smooth_method"] = _get("spectral", "spec_smooth_method", "binned")

        # Structure creation mode
        self.meta["structure_creation"] = _get("structure", "structure_creation", "default")
        self.meta["fname_fill"]          = _get("structure", "fname_fill",          "")

    def _load_target_definitions(self):
        """
        Load target_definitions.txt (tab-separated, '#' comments).

        Columns: source  ra_ctr  dec_ctr  dist_mpc  e_dist_mpc
                 incl_deg  e_incl_deg  posang_deg  e_posang_deg  r25  e_r25

        The sources to process may be a subset defined in data_key.txt
        [sources]. If not specified there, all rows are used.
        """
        geom_path = Path(self.meta["geom_file"])
        if not geom_path.exists():
            raise FileNotFoundError(f"target_definitions not found: {geom_path}")

        self.source_table = pd.read_csv(
            geom_path, sep="\t", names=TARGET_COLUMNS, comment="#"
        )

    def _load_data_key(self):
        """
        Load data_key.txt.

        Expected ini-style header section [sources] and [overlay], then
        freeform tabular sections for bands, cubes, and masks separated by
        labelled comment lines (mirrors the original .conf format).

        [sources]
        sources = ngc5194, ngc5457

        [overlay]
        overlay_file = _12co21.fits

        # ---- bands ----
        spire250, SPIRE250, MJy/sr, _spire250_gauss21.fits, ./data/

        # ---- cubes ----
        12co21, 12CO2-1, K, _12co21.fits, data/

        # ---- mask ----
        (optional)
        """
        imaging_path = Path(self.meta["data_key"])
        if not imaging_path.exists():
            raise FileNotFoundError(f"data_key not found: {imaging_path}")

        map_rows = []
        cube_rows = []
        mask_rows = []
        section = None

        # --- Pass 1: read ini-style sections with configparser ---
        # We only read the file up to the first tabular section so that
        # comma-separated data rows do not confuse configparser.
        ini_lines = []
        with open(imaging_path, "r") as f:
            for raw_line in f:
                stripped = raw_line.strip()
                low = stripped.lower()
                # Stop feeding ini lines once we hit a tabular section header
                if stripped.startswith("#") and any(
                    kw in low for kw in ("---- map", "---- cube", "---- mask")
                ):
                    break
                ini_lines.append(raw_line)

        cfg = configparser.ConfigParser(inline_comment_prefixes=("#",))
        cfg.read_string("".join(ini_lines))

        # Parse sources
        if "sources" in cfg:
            raw = cfg["sources"].get("sources", "")
            self.sources = [s.strip() for s in raw.split(",") if s.strip()]
        else:
            self.sources = list(self.source_table["source"])

        # Parse overlay
        if "overlay" in cfg:
            self.meta["overlay_file"] = cfg["overlay"].get("overlay_file", "")
        else:
            self.meta["overlay_file"] = ""

        # --- Pass 2: parse the tabular sections from free-form lines ---
        with open(imaging_path, "r") as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue
                low = line.lower()
                if "---- map" in low and line.startswith("#"):
                    section = "bands"
                    continue
                if "---- cube" in low and line.startswith("#"):
                    section = "cubes"
                    continue
                if "---- mask" in low and line.startswith("#"):
                    section = "mask"
                    continue
                if line.startswith("#") or line.startswith("["):
                    continue

                # skip ini-style key = value lines that have no comma
                if "=" in line and "," not in line:
                    continue

                parts = [p.strip() for p in line.split(",")]
                if section == "bands" and len(parts) >= 4:
                    while len(parts) < len(MAP_COLUMNS):
                        parts.append("")
                    map_rows.append(parts[:len(MAP_COLUMNS)])
                elif section == "cubes" and len(parts) >= 4:
                    while len(parts) < len(CUBE_COLUMNS):
                        parts.append("")
                    cube_rows.append(parts[:len(CUBE_COLUMNS)])
                elif section == "mask" and len(parts) >= 3:
                    mask_rows.append(parts)

        self.maps = pd.DataFrame(map_rows, columns=MAP_COLUMNS)
        self.cubes = pd.DataFrame(cube_rows, columns=CUBE_COLUMNS)

        # Resolve relative map_dir and line_dir to absolute paths so the
        # pipeline works regardless of the current working directory.
        _base = Path(self.meta.get("_base", "."))
        if len(self.maps) > 0:
            self.maps["map_dir"] = self.maps["map_dir"].apply(
                lambda d: str((_base / d.strip()).resolve()) if d.strip() else d
            )
        if len(self.cubes) > 0:
            self.cubes["line_dir"] = self.cubes["line_dir"].apply(
                lambda d: str((_base / d.strip()).resolve()) if d.strip() else d
            )

        if self.meta.get("use_fixed_vel_mask"):
            cols = MASK_COLUMNS_VEL
        else:
            cols = MASK_COLUMNS_FILE

        if mask_rows:
            padded = [r + [""] * max(0, len(cols) - len(r)) for r in mask_rows]
            self.input_mask = pd.DataFrame(
                [r[:len(cols)] for r in padded], columns=cols
            )
        else:
            self.input_mask = pd.DataFrame(columns=cols)

    def _load_hfs_key(self):
        """Load optional hyperfine structure file."""
        hfs_path = self.meta.get("hfs_file")
        if not hfs_path:
            self.hfs_data = None
            return
        hfs_path = Path(hfs_path)
        if not hfs_path.exists():
            self.hfs_data = None
            return
        self.hfs_data = pd.read_csv(
            hfs_path, sep="\t", names=HFS_COLUMNS, comment="#"
        )

    # ------------------------------------------------------------------
    # Validation helpers
    # ------------------------------------------------------------------

    def validate(self):
        """
        Basic validation of loaded keys.  Prints warnings for common issues.
        """
        issues = []

        if self.maps is None or len(self.maps) == 0:
            issues.append("No bands defined in data_key.")
        if self.cubes is None or len(self.cubes) == 0:
            issues.append("No cubes defined in data_key.")
        if not self.sources:
            issues.append("No sources defined.")
        if not self.meta.get("overlay_file"):
            issues.append("No overlay_file defined in data_key.")

        for issue in issues:
            print(f"[WARNING] KeyHandler: {issue}")

        return len(issues) == 0

    def __repr__(self):
        n_maps = len(self.maps) if self.maps is not None else 0
        n_cubes = len(self.cubes) if self.cubes is not None else 0
        return (
            f"KeyHandler(key_dir='{self.key_dir}', "
            f"sources={self.sources}, "
            f"n_maps={n_maps}, n_cubes={n_cubes})"
        )
