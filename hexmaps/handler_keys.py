"""
handler_keys.py — KeyHandler: reads and validates all HexMaps configuration.

HexMaps configuration lives in two places:

config.txt
    A single file in the working directory containing everything needed to
    run the pipeline: paths/metadata (formerly master_key.txt), the source
    list/overlay/maps/cubes/mask tables (formerly data_key.txt), and all
    numerical/boolean pipeline settings (formerly config_key.txt). This is
    the file you pass to the CLI via ``hexmaps --conf config.txt`` and
    the one you're expected to edit on every run.

keys/ subfolder (next to config.txt)
    target_definitions.txt
        Tab-separated table of source geometric parameters: RA/Dec centre,
        distance, inclination, position angle, and optical radius. One row
        per source. All sources that may ever be processed should be listed
        here; the subset to actually run is controlled by config.txt
        [sources]. Kept separate because this table is normally shared
        across many projects and changes rarely. Its path is set via
        [paths] geom_file in config.txt (default: keys/target_definitions.txt)
        and is REQUIRED — KeyHandler raises if it is missing.
    hfs_lines.txt (optional)
        Hyperfine structure line definitions. Also normally shared and
        rarely changed. Its path is set via [paths] hfs_file in config.txt
        (default: keys/hfs_lines.txt) and is OPTIONAL — if not found, HFS
        correction is simply unavailable and no error is raised.

This keeps the file you edit constantly (config.txt) separate from the
reference tables you set up once and reuse (keys/).
"""

import os
import re
import configparser
import numpy as np
import pandas as pd
from pathlib import Path

from hexmaps.logger import get_logger

LOG = get_logger("Loading")


# ---------------------------------------------------------------------------
# Column name definitions for the tabular sections of config.txt
#
# MAP_COLUMNS:  columns expected in the "---- maps ----" section
# CUBE_COLUMNS: columns expected in the "---- cubes ----" section
# MASK_COLUMNS_VEL:  columns for a fixed-velocity-window mask
# MASK_COLUMNS_FILE: columns for an external FITS mask file
# TARGET_COLUMNS: columns in keys/target_definitions.txt
# HFS_COLUMNS:  columns in the optional keys/hfs_lines.txt
# ---------------------------------------------------------------------------

MAP_COLUMNS = ["map_name", "map_desc", "map_unit", "map_ext", "map_dir", "map_uc"]
CUBE_COLUMNS = [
    "line_name",
    "line_desc",
    "line_unit",
    "line_ext",
    "line_dir",
    "map_ext",
    "map_uc",
]
MASK_COLUMNS_VEL = ["mask_name", "mask_desc", "mask_start", "mask_end", "mask_unit"]
MASK_COLUMNS_FILE = ["mask_name", "mask_desc", "mask_ext", "mask_dir"]
# Columns for noise-estimation velocity windows (same layout as MASK_COLUMNS_VEL)
NOISE_MASK_COLUMNS = ["mask_name", "mask_desc", "mask_start", "mask_end", "mask_unit"]
TARGET_COLUMNS = [
    "source",
    "ra_ctr",
    "dec_ctr",
    "dist_mpc",
    "e_dist_mpc",
    "incl_deg",
    "e_incl_deg",
    "posang_deg",
    "e_posang_deg",
    "r25",
    "e_r25",
]
HFS_COLUMNS = ["hfs_name", "hfs_ref_freq", "hfs_freq", "unit"]


class KeyHandler:
    """
    Reads and validates all HexMaps configuration from config.txt.

    The handler is the single source of truth for all pipeline configuration.
    Every other pipeline module receives either ``meta`` (a plain dict of
    scalar settings) or one of the DataFrames (``maps``, ``cubes``, etc.)
    returned by the getter methods below.

    Parameters
    ----------
    conf_path : str or Path
        Path to config.txt. The geometry table (required; default
        keys/target_definitions.txt) and optional HFS file (default
        keys/hfs_lines.txt) both default to a `keys/` subfolder next to this
        file, but either can be pointed elsewhere via [paths] geom_file /
        [paths] hfs_file in config.txt.

    Attributes
    ----------
    meta         : dict   — scalar settings from [paths]/[meta]/[resolution]/
                            [masking]/[spectral]/[output]/[structure]
    sources      : list   — source names to process (from config.txt [sources])
    source_table : pd.DataFrame — full geometry table from target_definitions
    maps         : pd.DataFrame — 2D map definitions (from config.txt)
    cubes        : pd.DataFrame — spectral cube definitions (from config.txt)
    input_mask   : pd.DataFrame — mask definition (from config.txt, may be empty)
    hfs_data     : pd.DataFrame or None — hyperfine structure data (optional)

    Example
    -------
    >>> kh = KeyHandler("./config.txt")
    >>> print(kh.sources)
    >>> print(kh.maps)
    """

    def __init__(self, conf_path: str):
        self.conf_path = Path(conf_path)
        self._validate_conf_path()

        # All parsed data; populated in load()
        self.meta = {}
        self.sources = []
        self.source_table = None
        self.maps = None
        self.cubes = None
        self.input_mask = None
        self.noise_mask = None
        self.hfs_data = None

        self.load()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def load(self):
        """
        Load the full configuration in dependency order.

        [resolution]/[masking]/[spectral]/[output]/[structure] are parsed
        before the [sources]/maps/cubes/mask tables so that masking flags
        (use_fixed_vel_mask etc.) are available when parsing the mask table.
        _resolve_resolution runs last because it needs both overlay_file
        (set by _load_sources_and_tables) and source distances (set by
        _load_target_definitions).
        """
        self._load_paths_and_meta()
        self._load_settings()
        self._load_target_definitions()
        self._load_sources_and_tables()
        self._resolve_resolution()
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

    def get_noise_mask(self) -> pd.DataFrame:
        """Return the DataFrame of noise velocity windows (may be empty)."""
        return self.noise_mask

    def get_source_table(self) -> pd.DataFrame:
        """Return the full source geometry table."""
        return self.source_table

    def get_hfs_data(self):
        """Return the hyperfine structure DataFrame, or None if not configured."""
        return self.hfs_data

    # ------------------------------------------------------------------
    # Private loaders
    # ------------------------------------------------------------------

    def _validate_conf_path(self):
        """Raise FileNotFoundError if config.txt does not exist."""
        if not self.conf_path.is_file():
            LOG.error(f"Config file not found: {self.conf_path}")
            raise FileNotFoundError(f"Config file not found: {self.conf_path}")

    def _ini_lines_for_configparser(self, path: Path) -> str:
        """
        Return the portion of *path* that is safe to feed to configparser.

        config.txt interleaves standard ``[section]`` blocks with free-form
        comma-separated tables (the maps/cubes/mask rows). configparser can
        only handle the former, so this strips out exactly the table rows —
        i.e. every line between a ``# ---- maps/cubes/mask ----`` divider and
        the next ``[section]`` header or another divider — while keeping all
        ``[section]`` blocks, including ones that appear *after* the tables
        (such as [resolution], [masking], [spectral], [output], [structure]
        in the unified config.txt).

        The divider-detection regex matches only the exact divider comment
        lines, not comment prose elsewhere in the file that happens to
        mention "maps"/"cubes"/"mask".
        """
        ini_lines = []
        in_table = False
        with open(path, "r") as f:
            for raw_line in f:
                stripped = raw_line.strip()

                if re.match(r"^#\s*----\s*(map|cube|mask)", stripped, re.IGNORECASE):
                    in_table = True
                    continue

                if in_table and stripped.startswith("["):
                    # A new [section] header ends the table region, even
                    # without an explicit divider comment.
                    in_table = False

                if in_table:
                    continue

                ini_lines.append(raw_line)
        return "".join(ini_lines)

    def _load_paths_and_meta(self):
        """
        Parse the [paths] and [meta] sections of config.txt.

        All path values are resolved to absolute paths relative to the
        directory containing config.txt, so the pipeline works regardless of
        the current working directory at runtime.

        Stores the resolved absolute base directory in ``self.meta["_base"]``
        for use by later loaders that need to resolve relative directories
        found in the maps/cubes tables.

        Expected format::

            [paths]
            data_dir   = data/
            out_dir    = output/
            # geom_file = keys/target_definitions.txt   (required; this is the default)
            # hfs_file  = keys/hfs_lines.txt             (optional)

            [meta]
            user     = Your Name
            comments = Free-form description of this run
        """
        # Only the [paths]/[meta] header is safe to feed to configparser
        # directly (the file also contains comma-separated tables further
        # down), so reuse the same "stop at first table divider" logic.
        ini_text = self._ini_lines_for_configparser(self.conf_path)
        cfg = configparser.ConfigParser(inline_comment_prefixes=("#",))
        cfg.read_string(ini_text)

        paths = dict(cfg["paths"]) if "paths" in cfg else {}
        meta = dict(cfg["meta"]) if "meta" in cfg else {}

        def _get_path(key, fallback):
            if key in paths:
                return paths[key]
            LOG.warning(
                f"[paths] {key} not set in config.txt; using default: {fallback}"
            )
            return fallback

        def _get_meta(key, fallback):
            if key in meta:
                return meta[key]
            LOG.warning(
                f"[meta] {key} not set in config.txt; using default: {fallback!r}"
            )
            return fallback

        # Resolve all paths relative to the directory containing config.txt.
        # Using .resolve() converts conf_path to an absolute path first, so
        # this is safe even when conf_path itself is given as a relative path.
        base = self.conf_path.resolve().parent

        self.meta["data_dir"] = str(base / _get_path("data_dir", "data/"))
        self.meta["out_dir"] = str(base / _get_path("out_dir", "output/"))
        self.meta["folder_savefits"] = str(
            base / _get_path("folder_savefits", "./saved_fits_files/")
        )
        # geom_file and hfs_file intentionally fall back to a default path
        # next to config.txt without a warning: this is documented behavior
        # (geom_file is required regardless and will raise its own error if
        # missing; hfs_file is optional and its "fallback" already depends
        # on whether the default file happens to exist, not just on whether
        # the key was set).
        self.meta["geom_file"] = str(
            base / paths.get("geom_file", "keys/target_definitions.txt")
        )
        self.meta["hfs_file"] = (
            str(base / paths.get("hfs_file", "")) if paths.get("hfs_file") else None
        )

        self.meta["user"] = _get_meta("user", "Unknown user")
        self.meta["comments"] = _get_meta("comments", "")

        # Store the absolute config file path so downstream code (e.g.
        # run_regrid) can read and embed the config content in the .ecsv.
        self.meta["conf_path"] = str(self.conf_path.resolve())

        # Store the absolute project root so _load_sources_and_tables can
        # resolve relative map_dir / line_dir entries to absolute paths.
        self.meta["_base"] = str(base)

        # geom_file (target_definitions.txt) is REQUIRED — handled like any
        # other [paths] entry, with a default pointing at keys/ next to
        # config.txt, but _load_target_definitions raises if it's missing.
        #
        # hfs_file (hfs_lines.txt) is OPTIONAL — if not explicitly set in
        # [paths], fall back to keys/hfs_lines.txt next to config.txt, but
        # only use it if that file actually exists; otherwise HFS correction
        # is simply unavailable (no error).
        if not self.meta["hfs_file"]:
            default_hfs = base / "keys" / "hfs_lines.txt"
            self.meta["hfs_file"] = str(default_hfs) if default_hfs.exists() else None

    def _load_settings(self):
        """
        Parse the [resolution], [masking], [spectral], [output], and
        [structure] sections of config.txt.

        All values have sensible defaults, so a minimal config.txt with only
        the settings you want to change is perfectly valid.

        Resolution settings
        -------------------
        target_res       : float  — target beam FWHM (arcsec for angular mode,
                                    pc for physical mode). MANDATORY.
        resolution       : str    — "angular" | "physical" | "native". MANDATORY.
        pixels_per_beam : float  — number of sampling points per beam diameter
        max_rad          : float | "auto"  — maximum map radius in degrees
        NAXIS_shuff      : int    — number of channels in the shuffled spectrum
        CDELT_SHUFF      : float  — channel width of the shuffled spectrum (m/s)

        Masking settings
        ----------------
        ref_line          : str   — which line to use for mask construction. MANDATORY.
        SN_processing     : list  — [low_SN, high_SN] thresholds
        strict_mask       : bool  — apply spatial connectivity filter
        use_input_mask    : bool  — use an external FITS mask from the [mask] table
        use_fixed_vel_mask: bool  — use a fixed velocity-window mask
        use_fixed_noise_mask: bool — use explicit velocity windows for noise estimation
        use_hfs_lines     : bool  — apply HFS correction (requires hfs_file)
        fov_erosion_beams : float — FOV erosion in units of the beam FWHM (default 0.5);
                                    set to 0 to disable erosion entirely
        mom_thresh        : float — S/N threshold for moment computation
        conseq_channels   : int   — minimum consecutive channels for valid mask
        mom2_method       : str   — "fwhm" | "sqrt" | "math"

        Output settings
        ---------------
        save_cubes     : bool — save convolved PPV cubes as FITS files in the fits stage
        save_mom_maps  : bool — save moment maps as FITS files
        save_maps      : bool — save 2D map FITS files
        save_mask      : bool — save the velocity-integration mask(s) as a
                                3D FITS cube (default False)

        Spectral smoothing
        ------------------
        spec_smooth        : "default" | float (target resolution in km/s)
        spec_smooth_method : "binned" | "gauss" | "combined"

        Structure creation
        ------------------
        structure_creation : "default" | "fill" | "archive"
        fname_fill         : str — pin a specific output filename for fill mode
                                   (rarely used; no fallback warning is logged
                                   when left unset)
        """
        ini_text = self._ini_lines_for_configparser(self.conf_path)
        cfg = configparser.ConfigParser(inline_comment_prefixes=("#",))
        cfg.read_string(ini_text)

        def _get(section, key, fallback, warn=True):
            if cfg.has_option(section, key):
                return cfg.get(section, key)
            if warn:
                LOG.warning(
                    f"[{section}] {key} not set in config.txt; "
                    f"using default: {fallback}"
                )
            return str(fallback)

        def _require(section, key):
            """Read *key* from *section*; raise ConfigError if absent."""
            if cfg.has_option(section, key):
                return cfg.get(section, key)
            LOG.error(
                f"Mandatory key '{key}' missing from [{section}] in config.txt. "
                f"This key must be explicitly set — there is no default value."
            )
            raise KeyError(
                f"Mandatory key '{key}' missing from [{section}] in config.txt. "
                f"This key must be explicitly set — there is no default value."
            )

        # Resolution
        self.meta["target_res"] = float(_require("resolution", "target_res"))
        self.meta["resolution"] = _require("resolution", "resolution")
        self.meta["pixels_per_beam"] = float(_get("resolution", "pixels_per_beam", 2.0))
        self.meta["max_rad"] = _get("resolution", "max_rad", "auto")
        self.meta["NAXIS_shuff"] = int(float(_get("resolution", "NAXIS_shuff", 200)))
        self.meta["CDELT_SHUFF"] = float(_get("resolution", "CDELT_SHUFF", 4000.0))

        # Masking
        self.meta["ref_line"] = _require("masking", "ref_line")
        self.meta["SN_processing"] = [
            float(x) for x in _get("masking", "SN_processing", "2,4").split(",")
        ]
        self.meta["strict_mask"] = (
            _get("masking", "strict_mask", "false").lower() == "true"
        )
        self.meta["use_input_mask"] = (
            _get("masking", "use_input_mask", "false").lower() == "true"
        )
        self.meta["use_fixed_vel_mask"] = (
            _get("masking", "use_fixed_vel_mask", "false").lower() == "true"
        )
        self.meta["use_fixed_noise_mask"] = (
            _get("masking", "use_fixed_noise_mask", "false").lower() == "true"
        )
        self.meta["use_hfs_lines"] = (
            _get("masking", "use_hfs_lines", "false").lower() == "true"
        )
        self.meta["fov_erosion_beams"] = float(
            _get("masking", "fov_erosion_beams", 0.5)
        )
        self.meta["mom_thresh"] = float(_get("masking", "mom_thresh", 5.0))
        self.meta["conseq_channels"] = int(float(_get("masking", "conseq_channels", 3)))
        self.meta["mom2_method"] = _get("masking", "mom2_method", "fwhm")

        # Output
        self.meta["save_cubes"] = _get("output", "save_cubes", "true").lower() == "true"
        self.meta["save_mom_maps"] = (
            _get("output", "save_mom_maps", "true").lower() == "true"
        )
        self.meta["save_maps"] = _get("output", "save_maps", "true").lower() == "true"
        self.meta["save_mask"] = _get("output", "save_mask", "true").lower() == "true"

        # Spectral smoothing
        self.meta["spec_smooth"] = _get("spectral", "spec_smooth", "default")
        self.meta["spec_smooth_method"] = _get(
            "spectral", "spec_smooth_method", "binned"
        )

        # Structure creation
        self.meta["structure_creation"] = _get(
            "structure", "structure_creation", "default"
        )
        # fname_fill is an optional, rarely-used override (only relevant
        # when structure_creation = "fill"), so its fallback to "" is
        # expected for most users and not worth a warning.
        self.meta["fname_fill"] = _get("structure", "fname_fill", "", warn=False)

    def _load_target_definitions(self):
        """
        Parse the geometry table at [paths] geom_file in config.txt
        (default: keys/target_definitions.txt next to config.txt).

        This file is REQUIRED — unlike hfs_file, there is no "just don't use
        it" fallback, since every source needs geometry. If the resolved
        path does not exist, this raises FileNotFoundError.

        The file is a comma-separated table with no header row.  Any spaces
        or tabs surrounding a comma are ignored, so columns can be aligned
        with extra whitespace for readability.  Comment lines beginning with
        '#' are ignored.  Columns must appear in the order defined by
        TARGET_COLUMNS.

        The full table is stored in ``self.source_table`` (a DataFrame). The
        subset of sources to actually process is determined later when the
        [sources] section of config.txt is parsed; at this stage we load
        everything.
        """
        geom_path = Path(self.meta["geom_file"])
        if not geom_path.exists():
            LOG.error(f"target_definitions not found: {geom_path}")
            raise FileNotFoundError(f"target_definitions not found: {geom_path}")
        # Comma-separated, but tolerant of stray spaces/tabs around each
        # field (e.g. "ngc5194,  202.4696,\t47.1952"). A regex separator
        # consumes the comma plus any surrounding whitespace in one go.
        #
        # Galaxy geometry columns (incl_deg, posang_deg, r25 and their
        # uncertainties) are OPTIONAL: rows may omit them entirely, or supply
        # empty/blank fields. Missing columns are filled with NaN so the rest
        # of the pipeline can detect absence and skip galaxy-specific
        # computations with an appropriate warning.
        self.source_table = pd.read_csv(
            geom_path,
            sep=r"\s*,\s*",
            engine="python",
            names=TARGET_COLUMNS,
            comment="#",
        )
        # Coerce all numeric columns; blank / "nan" / empty strings → NaN
        for col in TARGET_COLUMNS[1:]:   # skip "source"
            if col in self.source_table.columns:
                self.source_table[col] = pd.to_numeric(
                    self.source_table[col], errors="coerce"
                )
        # Fill any columns that are entirely absent (fewer columns in file)
        for col in TARGET_COLUMNS:
            if col not in self.source_table.columns:
                self.source_table[col] = float("nan")

    def _load_sources_and_tables(self):
        """
        Parse the [sources], [overlay], and maps/cubes/mask tables of config.txt.

        config.txt has a hybrid format: an ini-style header (parsed by
        configparser, shared with _load_paths_and_meta / _load_settings)
        followed by free-form comma-separated tabular sections for maps,
        cubes, and an optional mask.

        Parsing strategy
        ----------------
        **Pass 1** — configparser reads only the lines before the first
        tabular section divider, giving access to [sources] and [overlay]
        without configparser choking on the comma-separated data rows that
        follow.

        **Pass 2** — a simple line-by-line parser reads the tabular sections.
        Lines are routed to the correct section based on the most recently
        seen divider comment (``# ---- maps ----`` etc). Each comma-separated
        row is padded or trimmed to match the expected column count for that
        section.

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
        ini_text = self._ini_lines_for_configparser(self.conf_path)
        cfg = configparser.ConfigParser(inline_comment_prefixes=("#",))
        cfg.read_string(ini_text)

        # Source list — mandatory
        if "sources" not in cfg or not cfg["sources"].get("sources", "").strip():
            LOG.error(
                "Mandatory key 'sources' missing from [sources] in config.txt. "
                "This key must be explicitly set — there is no default value."
            )
            raise KeyError(
                "Mandatory key 'sources' missing from [sources] in config.txt. "
                "This key must be explicitly set — there is no default value."
            )
        self.sources = [
            s.strip() for s in cfg["sources"]["sources"].split(",") if s.strip()
        ]

        # Overlay file extension — mandatory
        if "overlay" not in cfg or not cfg["overlay"].get("overlay_file", "").strip():
            LOG.error(
                "Mandatory key 'overlay_file' missing from [overlay] in config.txt. "
                "This key must be explicitly set — there is no default value."
            )
            raise KeyError(
                "Mandatory key 'overlay_file' missing from [overlay] in config.txt. "
                "This key must be explicitly set — there is no default value."
            )
        self.meta["overlay_file"] = cfg["overlay"]["overlay_file"]

        # ------------------------------------------------------------------
        # Pass 2: parse the tabular sections line by line
        # ------------------------------------------------------------------
        map_rows, cube_rows, mask_rows, noise_mask_rows = [], [], [], []
        section = None

        with open(self.conf_path, "r") as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line:
                    continue
                low = line.lower()

                # Update the current section based on divider comments
                if "---- map" in low and line.startswith("#"):
                    section = "maps"
                    continue
                if "---- cube" in low and line.startswith("#"):
                    section = "cubes"
                    continue
                if "---- mask" in low and line.startswith("#"):
                    section = "mask"
                    continue

                # Skip all other comments and ini-style [section]/key=value lines
                if line.startswith("#") or line.startswith("["):
                    continue
                if "=" in line and "," not in line:
                    continue

                parts = [p.strip() for p in line.split(",")]

                if section == "maps" and len(parts) >= 4:
                    while len(parts) < len(MAP_COLUMNS):
                        parts.append("")
                    map_rows.append(parts[: len(MAP_COLUMNS)])

                elif section == "cubes" and len(parts) >= 4:
                    while len(parts) < len(CUBE_COLUMNS):
                        parts.append("")
                    cube_rows.append(parts[: len(CUBE_COLUMNS)])

                elif section == "mask" and len(parts) >= 3:
                    # Route by first field: "noise_mask" rows go to noise_mask_rows
                    if parts[0].lower() == "noise_mask":
                        noise_mask_rows.append(parts)
                    else:
                        mask_rows.append(parts)

        # Build DataFrames
        self.maps = pd.DataFrame(map_rows, columns=MAP_COLUMNS)
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
        cols = (
            MASK_COLUMNS_VEL
            if self.meta.get("use_fixed_vel_mask")
            else MASK_COLUMNS_FILE
        )
        if mask_rows:
            padded = [r + [""] * max(0, len(cols) - len(r)) for r in mask_rows]
            self.input_mask = pd.DataFrame(
                [r[: len(cols)] for r in padded], columns=cols
            )
        else:
            self.input_mask = pd.DataFrame(columns=cols)

        # Build noise_mask DataFrame (always velocity-window format)
        if noise_mask_rows:
            padded = [
                r + [""] * max(0, len(NOISE_MASK_COLUMNS) - len(r))
                for r in noise_mask_rows
            ]
            self.noise_mask = pd.DataFrame(
                [r[: len(NOISE_MASK_COLUMNS)] for r in padded],
                columns=NOISE_MASK_COLUMNS,
            )
            LOG.info(
                f"Loaded {len(self.noise_mask)} noise velocity window(s) "
                "from the [mask] table."
            )
        else:
            self.noise_mask = pd.DataFrame(columns=NOISE_MASK_COLUMNS)

    def _resolve_resolution(self):
        """
        Resolve the target resolution to arcseconds and parsecs, and compute
        the filename suffix. Writes three keys into ``self.meta``:

        ``target_res``
            Always **arcseconds**. For angular/native mode this is the primary
            working resolution; for physical mode it is converted from parsecs
            using the first source's distance. For native mode a placeholder
            value is stored here; ``run_sampling`` overwrites it with the
            exact beam size read from the overlay FITS header once the source
            loop begins.

        ``target_res_pc``
            Always **parsecs**. Derived from ``target_res`` and the first
            source's distance (dist_mpc).

        ``res_suffix``
            Single string used as the resolution part of all output filenames,
            e.g. ``"27p0as"``, ``"12p8as"``, ``"100pc"``.

        Called last in ``load()`` so that both ``overlay_file``
        (from _load_sources_and_tables) and ``source_table`` (from
        _load_target_definitions) are available.

        No log messages are emitted here — resolution logging is deferred
        to ``run_sampling`` where the source context is known and the overlay
        header is actually opened.
        """
        resolution = self.meta.get("resolution", "angular")
        # target_res as written in config — arcsec for angular/native, pc for physical
        target_res_config = float(self.meta.get("target_res", 27.0))

        # First source's distance for physical mode and pc conversion.
        if self.source_table is not None and len(self.source_table) > 0:
            dist_mpc = float(self.source_table["dist_mpc"].iloc[0])
        else:
            dist_mpc = 1.0

        if resolution == "physical":
            # Convert pc → arcsec using the first source's distance as a
            # placeholder.  run_sampling will recompute per-source correctly.
            target_res_as = (
                3600.0 * 180.0 / np.pi * 1e-6 * target_res_config / dist_mpc
            )
            # Store the original pc value separately so run_sampling can
            # re-convert it correctly for each source's distance.
            self.meta["target_res_config"] = target_res_config

        elif resolution == "native":
            # Placeholder; run_sampling will overwrite with the overlay beam.
            target_res_as = target_res_config

        else:
            # Angular: config value is already in arcseconds
            target_res_as = target_res_config

        # Parsec equivalent (placeholder; overwritten per-source by run_sampling)
        target_res_pc = target_res_as / 3600.0 * np.pi / 180.0 * dist_mpc * 1e6

        # Filename suffix
        if resolution == "physical":
            res_suffix = str(int(round(target_res_config))) + "pc"
        else:
            res_suffix = str(np.round(target_res_as, 1)).replace(".", "p") + "as"

        # meta["target_res"] always holds arcseconds from here on.
        # For physical mode the original pc value is preserved in
        # meta["target_res_config"] so run_sampling can re-convert per source.
        self.meta["target_res"]    = target_res_as
        self.meta["target_res_pc"] = target_res_pc
        self.meta["res_suffix"]    = res_suffix
        
    def _load_hfs_key(self):
        """
        Load the optional keys/hfs_lines.txt file.

        The file is comma-separated (spaces/tabs around each comma are
        ignored) with columns: hfs_name, hfs_ref_freq,
        hfs_freq, unit.  If no hfs_file is configured (explicitly via
        [paths] hfs_file, or implicitly via keys/hfs_lines.txt existing) or
        the file does not exist, ``self.hfs_data`` is set to None and no
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
        # Comma-separated, but tolerant of stray spaces/tabs around each
        # field — see the matching comment in _load_target_definitions.
        self.hfs_data = pd.read_csv(
            hfs_path,
            sep=r"\s*,\s*",
            engine="python",
            names=HFS_COLUMNS,
            comment="#",
        )

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
        if self.maps is None or len(self.maps) == 0:
            issues.append("No maps defined in config.txt.")
        if self.cubes is None or len(self.cubes) == 0:
            issues.append("No cubes defined in config.txt.")
        if not self.sources:
            issues.append("No sources defined.")
        if not self.meta.get("overlay_file"):
            issues.append("No overlay_file defined in config.txt.")

        for issue in issues:
            LOG.warning(f"{issue}")

        return len(issues) == 0

    def __repr__(self):
        n_maps = len(self.maps) if self.maps is not None else 0
        n_cubes = len(self.cubes) if self.cubes is not None else 0
        return (
            f"KeyHandler(conf_path='{self.conf_path}', "
            f"sources={self.sources}, n_maps={n_maps}, n_cubes={n_cubes})"
        )
