#!/usr/bin/env python3
"""
config_conversion.py — Convert a PyStructure.conf file (old format,
PhangsTeam/PyStructure) to a HexMaps config.txt (new format,
lukas-neumann-astro/PyStructure rename/hexmaps branch).

Usage:
    python config_conversion.py PyStructure.conf config.txt

The script reads every recognised key from the old flat config file and
maps it to the corresponding section and key in the new INI-style format.
Unrecognised lines are collected and appended as comments at the end so
nothing is silently dropped.

Old format reference : https://github.com/PhangsTeam/PyStructure
New format reference : https://github.com/lukas-neumann-astro/PyStructure
                       (branch rename/hexmaps)
"""

import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _strip_val(raw: str) -> str:
    """Remove surrounding quotes, whitespace and trailing comments."""
    v = raw.strip()
    # strip inline comment
    v = re.sub(r"\s*#.*$", "", v)
    # strip surrounding quotes (single or double)
    v = v.strip("'\"")
    return v.strip()


def _parse_old_config(path: Path) -> dict:
    """
    Parse the old flat PyStructure.conf into a dict of {key: raw_value}.
    Table rows (bands, cubes, masks) are collected separately.
    Lines that cannot be parsed as key=value are kept as 'unknown'.
    """
    data = {
        "kv": {},          # key -> stripped value
        "bands": [],       # raw band table lines
        "cubes": [],       # raw cube table lines
        "masks": [],       # raw mask table lines
        "unknown": [],     # unrecognised non-empty, non-comment lines
    }

    # Section tracking for table rows
    section = None
    # A line is a table row if it starts with a word followed by commas
    table_re = re.compile(r"^\s*\w[\w\d]*\s*[,\t]")

    # Keys that appear in the old config for which there is no equivalent
    skip_keys = {"save_fits", "save_band_maps"}

    with open(path, encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            stripped = line.strip()

            if not stripped or stripped.startswith("#"):
                # Track section from comment headers
                low = stripped.lower()
                if "step 4" in low or "band" in low:
                    section = "bands"
                elif "step 5" in low or "cube" in low:
                    section = "cubes"
                elif "step 6" in low or "mask" in low:
                    section = "masks"
                continue

            # key = value lines
            if "=" in stripped and not table_re.match(stripped):
                key, _, val_raw = stripped.partition("=")
                key = key.strip().lower()
                val = _strip_val(val_raw)
                if key not in skip_keys:
                    data["kv"][key] = val
                continue

            # Table rows (comma- or tab-separated)
            if section == "bands":
                data["bands"].append(stripped)
            elif section == "cubes":
                data["cubes"].append(stripped)
            elif section == "masks":
                data["masks"].append(stripped)
            else:
                data["unknown"].append(stripped)

    return data


def _convert_mask_rows(rows: list) -> tuple:
    """
    Convert old mask rows to new format rows.
    Returns (input_mask_rows, vel_mask_rows, noise_mask_rows,
             use_input_mask, use_fixed_vel_mask).
    """
    input_mask_rows = []
    vel_mask_rows = []
    noise_mask_rows = []
    use_input_mask = False
    use_fixed_vel_mask = False

    for row in rows:
        parts = [p.strip() for p in re.split(r"[\t,]+", row) if p.strip()]
        if len(parts) < 3:
            continue
        name = parts[0].lower()

        # Noise velocity ranges (old name: noise_vel)
        if name in ("noise_vel", "noise_mask"):
            noise_mask_rows.append(
                f"noise_mask, {', '.join(parts[1:])}"
            )
            continue

        # Detect by number of columns:
        # 4 cols → file mask (name, desc, ext, dir)
        # 5 cols → velocity window (name, desc, start, end, unit)
        if len(parts) == 5:
            # velocity window
            vel_mask_rows.append(f"vel_mask, {', '.join(parts[1:])}")
            use_fixed_vel_mask = True
        elif len(parts) >= 4:
            # file mask
            input_mask_rows.append(
                f"{parts[0]}, {', '.join(parts[1:])}"
            )
            use_input_mask = True

    return (input_mask_rows, vel_mask_rows, noise_mask_rows,
            use_input_mask, use_fixed_vel_mask)


def _convert_sn(val: str) -> str:
    """Convert '[2,4]' or '2,4' or '[2, 4]' to '2, 4'."""
    v = val.strip("[]() ")
    parts = [p.strip() for p in v.split(",")]
    if len(parts) == 2:
        return f"{parts[0]}, {parts[1]}"
    return val


def convert(old_path: Path, new_path: Path):
    d = _parse_old_config(old_path)
    kv = d["kv"]

    # -- resolve old→new key renames -----------------------------------------
    user          = kv.get("user", "")
    comments      = kv.get("comments", "")
    data_dir      = kv.get("data_dir", "data/")
    geom_file     = kv.get("geom_file", "")
    hfs_file      = kv.get("hfs_file", "")
    overlay_file  = kv.get("overlay_file", "")
    out_dir       = kv.get("out_dic", kv.get("out_dir", "output/"))
    folder_savefits = kv.get("folder_savefits", "./saved_fits_files/")
    sources       = kv.get("sources", "")
    target_res    = kv.get("target_res", "27.0")
    resolution    = kv.get("resolution", "angular")
    pixels_per_beam = kv.get("spacing_per_beam",
                              kv.get("pixels_per_beam", "2"))
    max_rad       = kv.get("max_rad", "auto")
    naxis_shuff   = kv.get("naxis_shuff", "200")
    cdelt_shuff   = kv.get("cdelt_shuff", "4000.0")
    ref_line      = kv.get("ref_line", "first")
    sn_processing = _convert_sn(kv.get("sn_processing", "2, 4"))
    strict_mask   = kv.get("strict_mask", "false").lower()
    use_input_mask     = kv.get("use_input_mask", "false").lower()
    use_fixed_vel_mask = kv.get("use_fixed_vel_mask", "false").lower()
    use_noise_vel      = kv.get("use_noise_vel_ranges",
                                 kv.get("use_fixed_noise_mask", "false")).lower()
    use_hfs_lines = kv.get("use_hfs_lines", "false").lower()
    mom_thresh    = kv.get("mom_thresh", "5")
    conseq_ch     = kv.get("conseq_channels", "3")
    mom2_method   = kv.get("mom2_method", "fwhm")
    spec_smooth   = kv.get("spec_smooth", "default")
    spec_smooth_method = kv.get("spec_smooth_method", "binned")
    save_mom_maps = kv.get("save_mom_maps", "true").lower()
    save_maps     = kv.get("save_band_maps",   # old name
                            kv.get("save_maps", "true")).lower()
    structure_creation = kv.get("structure_creation", "default")
    fname_fill    = kv.get("fname_fill", "")

    # -- handle mask table ----------------------------------------------------
    (input_mask_rows, vel_mask_rows, noise_mask_rows,
     has_input_mask, has_vel_mask) = _convert_mask_rows(d["masks"])

    # Override flags if mask rows were found
    if has_input_mask:
        use_input_mask = "true"
    if has_vel_mask:
        use_fixed_vel_mask = "true"
    if noise_mask_rows:
        use_noise_vel = "true"

    # -- build map table (new format) -----------------------------------------
    map_lines = []
    for row in d["bands"]:
        parts = [p.strip() for p in re.split(r"[\t,]+", row) if p.strip()]
        while len(parts) < 6:
            parts.append("")
        # old: name, desc, unit, ext, dir, uc_ext
        map_lines.append(
            f"{parts[0]},  {parts[1]},  {parts[2]},  {parts[3]},  "
            f"{parts[4]},  {parts[5]}"
        )

    # -- build cube table (new format) ----------------------------------------
    cube_lines = []
    for row in d["cubes"]:
        parts = [p.strip() for p in re.split(r"[\t,]+", row) if p.strip()]
        while len(parts) < 7:
            parts.append("")
        # old: name, desc, unit, ext, dir, map_ext, map_uc_ext
        cube_lines.append(
            f"{parts[0]},  {parts[1]},  {parts[2]},  {parts[3]},  "
            f"{parts[4]},  {parts[5]},  {parts[6]}"
        )

    # -- optional geom_file / hfs_file lines ----------------------------------
    geom_line = (f"geom_file = {geom_file}"
                 if geom_file and geom_file not in ("", "keys/target_definitions.txt")
                 else "# geom_file = keys/target_definitions.txt")
    hfs_line  = (f"hfs_file = {hfs_file}"
                 if hfs_file and hfs_file not in ("", "keys/hfs_lines.txt")
                 else "# hfs_file  = keys/hfs_lines.txt")

    fname_fill_line = (f"fname_fill = {fname_fill}"
                       if fname_fill else
                       "# fname_fill = <filename>.ecsv")

    # -- assemble new config --------------------------------------------------
    sections = []

    sections.append(f"""\
# =============================================================================
# HexMaps config.txt  (converted from {old_path.name})
# =============================================================================
# Converted by config_conversion.py
# Old format: PhangsTeam/PyStructure (PyStructure.conf)
# New format: lukas-neumann-astro/PyStructure  rename/hexmaps branch
# =============================================================================

[meta]
user = {user}
comments = {comments}

[paths]
data_dir         = {data_dir}
out_dir          = {out_dir}
{geom_line}
{hfs_line}
folder_savefits  = {folder_savefits}

[sources]
sources = {sources}

[overlay]
overlay_file = {overlay_file}

# ---- maps ----""")

    for ml in map_lines:
        sections.append(ml)
    if not map_lines:
        sections.append("# (no band/map entries found in old config)")

    sections.append("\n# ---- cubes ----")
    for cl in cube_lines:
        sections.append(cl)
    if not cube_lines:
        sections.append("# (no cube entries found in old config)")

    sections.append("\n# ---- mask ----")
    for row in input_mask_rows:
        sections.append(row)
    for row in vel_mask_rows:
        sections.append(row)
    for row in noise_mask_rows:
        sections.append(row)
    if not (input_mask_rows or vel_mask_rows or noise_mask_rows):
        sections.append("# (no mask entries found in old config)")

    sections.append(f"""

[resolution]
target_res      = {target_res}
resolution      = {resolution}
pixels_per_beam = {pixels_per_beam}
max_rad         = {max_rad}
NAXIS_shuff     = {naxis_shuff}
CDELT_SHUFF     = {cdelt_shuff}

[masking]
ref_line              = {ref_line}
SN_processing         = {sn_processing}
strict_mask           = {strict_mask}
use_input_mask        = {use_input_mask}
use_fixed_vel_mask    = {use_fixed_vel_mask}
use_fixed_noise_mask  = {use_noise_vel}
use_hfs_lines         = {use_hfs_lines}
fov_erosion_beams     = 0.5
mom_thresh            = {mom_thresh}
conseq_channels       = {conseq_ch}
mom2_method           = {mom2_method}

[spectral]
spec_smooth        = {spec_smooth}
spec_smooth_method = {spec_smooth_method}

[output]
save_cubes    = false
save_mom_maps = {save_mom_maps}
save_maps     = {save_maps}
save_mask     = false

[structure]
structure_creation = {structure_creation}
{fname_fill_line}""")

    # -- unknown / unrecognised lines -----------------------------------------
    if d["unknown"]:
        sections.append(
            "\n# ---- unrecognised lines from old config (review manually) ----"
        )
        for u in d["unknown"]:
            sections.append(f"# {u}")

    output = "\n".join(sections) + "\n"
    new_path.write_text(output, encoding="utf-8")
    print(f"[OK] Written: {new_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print("Usage: python config_conversion.py <old_PyStructure.conf> <new_config.txt>")
        sys.exit(1)
    old_path = Path(sys.argv[1])
    new_path = Path(sys.argv[2])
    if not old_path.exists():
        print(f"[ERROR] Input file not found: {old_path}")
        sys.exit(1)
    convert(old_path, new_path)


if __name__ == "__main__":
    main()
