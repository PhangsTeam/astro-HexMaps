"""
handler_keys.py — KeyHandler: reads and validates all PyStructure key files.

The pipeline configuration is intentionally split across four separate key files
rather than a single monolithic config file. This makes it easy to share source
geometry tables across projects, swap imaging configurations, and version-control
each concern independently.

Key files
---------
master_key.txt
    Paths to all other key files, the data and output directories, and
    free-form metadata (user name, comments). This is the only file the
    PipelineHandler needs to know about directly.

target_definitions.txt
    Tab-separated table of source geometric parameters: RA/Dec centre,
    distance, inclination, position angle, and optical radius. One row per
    source. All sources that may ever be processed should be listed here;
    the subset to actually run is controlled by data_key.txt [sources].

data_key.txt
    Defines which sources to process, the overlay FITS cube (used to define
    the sampling grid and spectral axis), and the input 2D maps and spectral
    cubes. Also contains the optional external mask definition.

config_key.txt
    All numerical and boolean pipeline settings: target resolution, masking
    thresholds, spectral smoothing, output flags, and structure-creation mode.
"""

import os
import re
import configparser
import pandas as pd
from pathlib import Path

from pystructurePipeline.pystructureLogger import get_logger

LOG = get_logger("Loading")


# ---------------------------------------------------------------------------
# Column name definitions for the tabular sections of data_key.txt
#
# MAP_COLUMNS:  columns expected in the "---- maps ----" section
# CUBE_COLUMNS: columns expected in the "---- cubes ----" section
# MASK_COLUMNS_VEL:  columns for a fixed-velocity-window mask
# MASK_COLUMNS_FILE: columns for an external FITS mask file
# TARGET_COLUMNS: columns in target_definitions.txt
# HFS_COLUMNS:  columns in the optional hyperfine-structure file
# ---------------------------------------------------------------------------

MAP_COLUMNS        = ["map_name", "map_desc", "map_unit", "map_ext", "map_dir", "map_uc"]
CUBE_COLUMNS       = ["line_name", "line_desc", "line_unit", "line_ext", "line_dir", "map_ext", "map_uc"]
MASK_COLUMNS_VEL   = ["mask_name", "mask_desc", "mask_start", "mask_end", "mask_unit"]
MASK_COLUMNS_FILE  = ["mask_name", "mask_desc", "mask_ext", "mask_dir"]
TARGET_COLUMNS     = [
    "source", "ra_ctr", "dec_ctr", "dist_mpc", "e_dist_mpc",
    "incl_deg", "e_incl_deg", "posang_deg", "e_posang_deg",
    "r25", "e_r25",
]
HFS_COLUMNS = ["hfs_name", "hfs_ref_freq", "hfs_freq", "unit"]


class KeyHandler:
    """
    Reads and validates all PyStructure key files from a key directory.

    The handler is the single source of truth for all pipeline configuration.
    Every other pipeline module receives either ``meta`` (a plain dict of
    scalar settings) or one of the DataFrames (``maps``, ``cubes``, etc.)
    returned by the getter methods below.

    Parameters
    ----------
    key_dir : str or Path
        Directory containing master_key.txt (and, by default, the other key
        files unless they are located elsewhere as specified in master_key).

    Attributes
    ----------
    meta         : dict   — scalar settings from master_key + config_key
    sources      : list   — source names to process (from data_key [sources])
    source_table : pd.DataFrame — full geometry table from target_definitions
    maps         : pd.DataFrame — 2D map definitions (from data_key)
    cubes        : pd.DataFrame — spectral cube definitions (from data_key)
    input_mask   : pd.DataFrame — mask definition (from data_key, may be empty)
    hfs_data     : pd.DataFrame or None — hyperfine structure data (optional)

    Example
    -------
    >>> kh = KeyHandler("./keys/")
    >>> print(kh.sources)
    >>> print(kh.maps)
    """

    def __init__(self, key_dir: str):
        self.key_dir = Path(key_dir)
        self._validate_key_dir()

        # All parsed data; populated in load()
        self.meta         = {}
        self.sources      = []
        self.source_table = None
        self.maps         = None
        self.cubes        = None
        self.input_mask   = None
        self.hfs_data     = None

        self.load()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def load(self):
        """
        Load all key files in dependency order.

        master_key is loaded first because it provides the paths to all
        other files.  config_key is loaded before data_key so that masking
        flags (use_fixed_vel_mask etc.) are available when parsing data_key.
        """
        self._load_master_key()
        self._load_config_key()
        self._load_target_definitions()
        self._load_data_key()
        self._load_hfs_key()

    def get_sources(self) -> list:
        """Return the ordered list of source names to be processed."""
        return list(self.sources)

    def get_maps(self) -> pd.DataFrame:
        """Return the DataFrame of 2D map definitions."""
        return self.maps

    def get_cubes(self) -> pd.DataFrame:
        """Return the DataFrame of spectral cube definitions."""
        return self.cubes

    def get_input_mask(self) -> pd.DataFrame:
        """Return the DataFrame of mask definitions (may be empty)."""
        return self.input_mask

    def get_source_table(self) -> pd.DataFrame:
        """Return the full source geometry table."""
        return self.source_table

    def get_hfs_data(self):
        """Return the hyperfine structure DataFrame, or None if not configured."""
        return self.hfs_data

    # ------------------------------------------------------------------
    # Private loaders
    # ------------------------------------------------------------------

    def _validate_key_dir(self):
        """Raise FileNotFoundError if key_dir does not exist."""
        if not self.key_dir.is_dir():
            LOG.error(
                f"Key directory not found: {self.key_dir}"
            )
            raise FileNotFoundError(
                f"Key directory not found: {self.key_dir}"
            )

    def _load_master_key(self):
        """
        Parse master_key.txt.

        Reads the [paths] and [meta] sections using configparser.  All path
        values are resolved to absolute paths relative to the *parent* of
        key_dir (i.e. the project root), so the pipeline works regardless of
        the current working directory at runtime.

        Stores the resolved absolute base path in ``self.meta["_base"]`` for
        use by later loaders that need to resolve relative directories found
        inside data_key.txt.

        Expected format::

            [paths]
            data_dir   = data/
            out_dir    = output/
            geom_file  = keys/target_definitions.txt
            data_key   = keys/data_key.txt
            config_key = keys/config_key.txt

            [meta]
            user     = Your Name
            comments = Free-form description of this run
        """
        master_path = self.key_dir / "master_key.txt"
        if not master_path.exists():
            LOG.error(
                f"master_key.txt not found in {self.key_dir}"
            )
            raise FileNotFoundError(
                f"master_key.txt not found in {self.key_dir}"
            )

        cfg = configparser.ConfigParser(inline_comment_prefixes=("#",))
        cfg.read(master_path)

        paths = dict(cfg["paths"]) if "paths" in cfg else {}
        meta  = dict(cfg["meta"])  if "meta"  in cfg else {}

        # Resolve all paths relative to the project root (key_dir's parent).
        # Using .resolve() converts key_dir to an absolute path first, so
        # this is safe even when key_dir itself is given as a relative path.
        base = self.key_dir.resolve().parent

        self.meta["data_dir"]    = str(base / paths.get("data_dir",    "data/"))
        self.meta["out_dir"]     = str(base / paths.get("out_dir",     "output/"))
        self.meta["geom_file"]   = str(base / paths.get("geom_file",   "keys/target_definitions.txt"))
        self.meta["data_key"]    = str(base / paths.get("data_key",    "keys/data_key.txt"))
        self.meta["config_key"]  = str(base / paths.get("config_key",  "keys/config_key.txt"))
        self.meta["hfs_file"]    = str(base / paths.get("hfs_file", "")) if paths.get("hfs_file") else None

        self.meta["user"]        = meta.get("user",     "Unknown user")
        self.meta["comments"]    = meta.get("comments", "")
        # Store the absolute project root so _load_data_key can resolve
        # relative map_dir / line_dir entries to absolute paths.
        self.meta["_base"]       = str(base)

    def _load_config_key(self):
        """
        Parse config_key.txt.

        Reads five ini-style sections — [resolution], [masking], [spectral],
        [output], [structure] — and stores all values in ``self.meta``.

        All settings have sensible defaults so a minimal config file with
        only the settings you want to change is perfectly valid.

        Resolution settings
        -------------------
        target_res       : float  — target beam FWHM (arcsec for angular mode,
                                    pc for physical mode)
        resolution       : str    — "angular" | "physical" | "native"
        spacing_per_beam : float  — number of sampling points per beam diameter
        max_rad          : float | "auto"  — maximum map radius in degrees
        NAXIS_shuff      : int    — number of channels in the shuffled spectrum
        CDELT_SHUFF      : float  — channel width of the shuffled spectrum (m/s)

        Masking settings
        ----------------
        ref_line          : str   — which line to use for mask construction
        SN_processing     : list  — [low_SN, high_SN] thresholds
        strict_mask       : bool  — apply spatial connectivity filter
        use_input_mask    : bool  — use an external FITS mask from data_key
        use_fixed_vel_mask: bool  — use a fixed velocity-window mask
        use_hfs_lines     : bool  — apply HFS correction (requires hfs_file)
        mom_thresh        : float — S/N threshold for moment computation
        conseq_channels   : int   — minimum consecutive channels for valid mask
        mom2_method       : str   — "fwhm" | "sqrt" | "math"

        Output settings
        ---------------
        save_fits      : bool — save convolved intermediate FITS cubes
        save_mom_maps  : bool — save moment maps as FITS files
        save_maps      : bool — save 2D map FITS files
        folder_savefits: str  — output folder for FITS maps

        Spectral smoothing
        ------------------
        spec_smooth        : "default" | float (target resolution in km/s)
        spec_smooth_method : "binned" | "gauss" | "combined"

        Structure creation
        ------------------
        structure_creation : "default" | "fill" | "archive"
        fname_fill         : str — pin a specific output filename for fill mode
        """
        config_path = Path(self.meta["config_key"])
        if not config_path.exists():
            LOG.error(
                f"config_key not found: {config_path}"
            )
            raise FileNotFoundError(
                f"config_key not found: {config_path}"
            )

        cfg = configparser.ConfigParser(inline_comment_prefixes=("#",))
        cfg.read(config_path)

        def _get(section, key, fallback):
            return cfg.get(section, key, fallback=str(fallback))

        # Resolution
        self.meta["target_res"]       = float(_get("resolution", "target_res",       27.0))
        self.meta["resolution"]       =       _get("resolution", "resolution",       "angular")
        self.meta["spacing_per_beam"] = float(_get("resolution", "spacing_per_beam", 2.0))
        self.meta["max_rad"]          =       _get("resolution", "max_rad",          "auto")
        self.meta["NAXIS_shuff"]      = int(float(_get("resolution", "NAXIS_shuff",  200)))
        self.meta["CDELT_SHUFF"]      = float(_get("resolution", "CDELT_SHUFF",      4000.0))

        # Masking
        self.meta["ref_line"]           =       _get("masking", "ref_line",           "first")
        self.meta["SN_processing"]      = [float(x) for x in _get("masking", "SN_processing", "2,4").split(",")]
        self.meta["strict_mask"]        =       _get("masking", "strict_mask",        "false").lower() == "true"
        self.meta["use_input_mask"]     =       _get("masking", "use_input_mask",     "false").lower() == "true"
        self.meta["use_fixed_vel_mask"] =       _get("masking", "use_fixed_vel_mask", "false").lower() == "true"
        self.meta["use_hfs_lines"]      =       _get("masking", "use_hfs_lines",      "false").lower() == "true"
        self.meta["mom_thresh"]         = float(_get("masking", "mom_thresh",         5.0))
        self.meta["conseq_channels"]    = int(float(_get("masking", "conseq_channels", 3)))
        self.meta["mom2_method"]        =       _get("masking", "mom2_method",        "fwhm")

        # Output
        self.meta["save_fits"]        =       _get("output", "save_fits",        "false").lower() == "true"
        self.meta["save_mom_maps"]    =       _get("output", "save_mom_maps",    "true").lower()  == "true"
        self.meta["save_maps"]        =       _get("output", "save_maps",        "true").lower()  == "true"
        self.meta["folder_savefits"]  =       _get("output", "folder_savefits",  "./saved_FITS_files/")

        # Spectral smoothing
        self.meta["spec_smooth"]        = _get("spectral", "spec_smooth",        "default")
        self.meta["spec_smooth_method"] = _get("spectral", "spec_smooth_method", "binned")

        # Structure creation
        self.meta["structure_creation"] = _get("structure", "structure_creation", "default")
        self.meta["fname_fill"]         = _get("structure", "fname_fill",         "")

    def _load_target_definitions(self):
        """
        Parse target_definitions.txt.

        The file is a tab-separated table with no header row.  Comment lines
        beginning with '#' are ignored.  Columns must appear in the order
        defined by TARGET_COLUMNS.

        The full table is stored in ``self.source_table`` (a DataFrame). The
        subset of sources to actually process is determined later when
        data_key.txt is parsed; at this stage we load everything.
        """
        geom_path = Path(self.meta["geom_file"])
        if not geom_path.exists():
            LOG.error(
                f"target_definitions not found: {geom_path}"
            )
            raise FileNotFoundError(
                f"target_definitions not found: {geom_path}"
            )
        self.source_table = pd.read_csv(
            geom_path, sep="\t", names=TARGET_COLUMNS, comment="#"
        )

    def _load_data_key(self):
        """
        Parse data_key.txt.

        The file has a hybrid format: an ini-style header (parsed by
        configparser) followed by free-form comma-separated tabular sections
        for maps, cubes, and an optional mask.

        Parsing strategy
        ----------------
        **Pass 1** — configparser reads only the lines before the first
        tabular section divider.  This safely extracts the [sources] and
        [overlay] sections without configparser choking on the comma-separated
        data rows that follow.

        The stop condition uses a precise regex (``^#\\s*----\\s*(map|cube|mask)``)
        that matches only the section-divider comment lines, not any comment
        prose in the file header that might mention those keywords.

        **Pass 2** — a simple line-by-line parser reads the tabular sections.
        Lines are routed to the correct section based on the most recently
        seen divider comment.  Each comma-separated row is padded or trimmed to
        match the expected column count for that section.

        Path resolution
        ---------------
        After building the DataFrames, ``map_dir`` and ``line_dir`` values are
        resolved to absolute paths using the project root stored in
        ``self.meta["_base"]``.  This ensures that file paths constructed
        later in stage_regrid.py are valid regardless of the current working
        directory.

        Expected format::

            [sources]
            sources = ngc5194, ngc5457

            [overlay]
            overlay_file = _12co21.fits

            # ---- maps ----
            spire250, SPIRE 250 um, MJy/sr, _spire250_gauss27.fits, data/

            # ---- cubes ----
            12co21, 12CO(2-1), K, _12co21.fits, data/
            12co10, 12CO(1-0), K, _12co10.fits, data/

            # ---- mask ----
            # (leave empty if no external mask is used)
        """
        imaging_path = Path(self.meta["data_key"])
        if not imaging_path.exists():
            LOG.error(
                f"data_key not found: {imaging_path}"
            )
            raise FileNotFoundError(
                f"data_key not found: {imaging_path}"
            )

        map_rows, cube_rows, mask_rows = [], [], []
        section = None

        # ------------------------------------------------------------------
        # Pass 1: feed only the ini-style header into configparser
        # ------------------------------------------------------------------
        ini_lines = []
        with open(imaging_path, "r") as f:
            for raw_line in f:
                stripped = raw_line.strip()
                # Stop as soon as we hit an exact section-divider comment.
                # The regex matches lines like "# ---- maps ----" but NOT
                # prose like "# MAP TABLE (after '# ---- maps ----')".
                if re.match(r'^#\s*----\s*(map|cube|mask)', stripped, re.IGNORECASE):
                    break
                ini_lines.append(raw_line)

        cfg = configparser.ConfigParser(inline_comment_prefixes=("#",))
        cfg.read_string("".join(ini_lines))

        # Source list: explicit [sources] section, or fall back to all targets
        if "sources" in cfg:
            raw = cfg["sources"].get("sources", "")
            self.sources = [s.strip() for s in raw.split(",") if s.strip()]
        else:
            self.sources = list(self.source_table["source"])

        # Overlay file extension (prepended with source name at runtime)
        if "overlay" in cfg:
            self.meta["overlay_file"] = cfg["overlay"].get("overlay_file", "")
        else:
            self.meta["overlay_file"] = ""

        # ------------------------------------------------------------------
        # Pass 2: parse the tabular sections line by line
        # ------------------------------------------------------------------
        with open(imaging_path, "r") as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue
                low = line.lower()

                # Update the current section based on divider comments
                if "---- map"  in low and line.startswith("#"):
                    section = "maps";  continue
                if "---- cube" in low and line.startswith("#"):
                    section = "cubes"; continue
                if "---- mask" in low and line.startswith("#"):
                    section = "mask";  continue

                # Skip all other comments and ini-style [section]/key=value lines
                if line.startswith("#") or line.startswith("["):
                    continue
                if "=" in line and "," not in line:
                    continue

                parts = [p.strip() for p in line.split(",")]

                if section == "maps" and len(parts) >= 4:
                    while len(parts) < len(MAP_COLUMNS):
                        parts.append("")
                    map_rows.append(parts[:len(MAP_COLUMNS)])

                elif section == "cubes" and len(parts) >= 4:
                    while len(parts) < len(CUBE_COLUMNS):
                        parts.append("")
                    cube_rows.append(parts[:len(CUBE_COLUMNS)])

                elif section == "mask" and len(parts) >= 3:
                    mask_rows.append(parts)

        # Build DataFrames
        self.maps  = pd.DataFrame(map_rows,  columns=MAP_COLUMNS)
        self.cubes = pd.DataFrame(cube_rows, columns=CUBE_COLUMNS)

        # ------------------------------------------------------------------
        # Resolve relative map_dir / line_dir to absolute paths.
        # This is essential so that stage_regrid can construct valid file
        # paths regardless of where the user runs the pipeline from.
        # ------------------------------------------------------------------
        _base = Path(self.meta.get("_base", "."))
        if len(self.maps) > 0:
            self.maps["map_dir"] = self.maps["map_dir"].apply(
                lambda d: str((_base / d.strip()).resolve()) if d.strip() else d
            )
        if len(self.cubes) > 0:
            self.cubes["line_dir"] = self.cubes["line_dir"].apply(
                lambda d: str((_base / d.strip()).resolve()) if d.strip() else d
            )

        # Build mask DataFrame with the appropriate column set
        cols = MASK_COLUMNS_VEL if self.meta.get("use_fixed_vel_mask") else MASK_COLUMNS_FILE
        if mask_rows:
            padded = [r + [""] * max(0, len(cols) - len(r)) for r in mask_rows]
            self.input_mask = pd.DataFrame([r[:len(cols)] for r in padded], columns=cols)
        else:
            self.input_mask = pd.DataFrame(columns=cols)

    def _load_hfs_key(self):
        """
        Load the optional hyperfine structure file.

        The file is tab-separated with columns: hfs_name, hfs_ref_freq,
        hfs_freq, unit.  If the hfs_file path is not set in master_key.txt
        or the file does not exist, ``self.hfs_data`` is set to None and no
        error is raised — HFS correction is simply not applied.
        """
        hfs_path = self.meta.get("hfs_file")
        if not hfs_path:
            self.hfs_data = None
            return
        hfs_path = Path(hfs_path)
        if not hfs_path.exists():
            self.hfs_data = None
            return
        self.hfs_data = pd.read_csv(hfs_path, sep="\t", names=HFS_COLUMNS, comment="#")

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate(self) -> bool:
        """
        Run basic sanity checks on the loaded configuration.

        Prints a [WARNING] for each problem found but does not raise.
        Returns True if all checks pass, False otherwise.

        Checks performed
        ----------------
        - At least one map defined
        - At least one cube defined
        - At least one source defined
        - overlay_file is set
        """
        issues = []
        if self.maps  is None or len(self.maps)  == 0:
            issues.append("No maps defined in data_key.")
        if self.cubes is None or len(self.cubes) == 0:
            issues.append("No cubes defined in data_key.")
        if not self.sources:
            issues.append("No sources defined.")
        if not self.meta.get("overlay_file"):
            issues.append("No overlay_file defined in data_key.")

        for issue in issues:
            LOG.warning(f"{issue}")

        return len(issues) == 0

    def __repr__(self):
        n_maps  = len(self.maps)  if self.maps  is not None else 0
        n_cubes = len(self.cubes) if self.cubes is not None else 0
        return (
            f"KeyHandler(key_dir='{self.key_dir}', "
            f"sources={self.sources}, n_maps={n_maps}, n_cubes={n_cubes})"
        )
