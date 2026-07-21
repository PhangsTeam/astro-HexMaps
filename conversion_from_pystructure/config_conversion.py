#!/usr/bin/env python3
"""
config_conversion.py — Convert a PyStructure v4.x config file
(PyStructure.conf, old flat key=value format) to a HexMaps config.txt
(new INI-style format, lukas-neumann-astro/astro-HexMaps main branch).

Usage:
    python config_conversion.py PyStructure.conf config.txt

    # If band/cube definitions live in separate list files:
    python config_conversion.py PyStructure.conf config.txt \\
        --band-list List_Files/band_list.txt \\
        --cube-list List_Files/cube_list.txt

Key changes from PyStructure v4 → HexMaps v5
---------------------------------------------
- Config format: flat key=value → INI sections
- out_dic → out_dir
- save_fits dropped (replaced by save_cubes in [output])
- save_band_maps → save_maps
- use_input_mask / use_fixed_vel_mask booleans → ref_line tokens (input/window)
- use_noise_vel_ranges → use_fixed_noise_mask
- strict_mask bool → false | strict | broad
- ref_line = 'ref+HI' is removed (warning issued)
- Mask table rows now use explicit keys: input_mask, window_mask, noise_mask
"""

import re
import sys
import argparse
from pathlib import Path

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _strip_val(raw):
    v = raw.strip()
    v = re.sub(r"\s*#.*$", "", v)
    v = v.strip("'\"").strip()
    return v


def _convert_sn(val):
    v = val.strip("[]() ")
    parts = [p.strip() for p in v.split(",")]
    return f"{parts[0]}, {parts[1]}" if len(parts) == 2 else val


def _parse_table_rows(lines):
    result = []
    for row in lines:
        row = re.sub(r"\s*#.*$", "", row).strip()
        if not row:
            continue
        parts = [p.strip() for p in re.split(r"[\t,]+", row) if p.strip()]
        if parts:
            result.append(parts)
    return result


def _read_list_file(path):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            rows.append(line)
    return rows


# ---------------------------------------------------------------------------
# Old-config parser
# ---------------------------------------------------------------------------


def _parse_old_config(path):
    data = {"kv": {}, "bands": [], "cubes": [], "masks": [], "unknown": []}
    _DROPPED = {"save_fits", "save_band_maps"}
    table_row_re = re.compile(r"^[^=]+[\t,]")
    section = None

    with open(path, encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                low = stripped.lower()
                if "step 4" in low or ("band" in low and "step" in low):
                    section = "bands"
                elif "step 5" in low or ("cube" in low and "step" in low):
                    section = "cubes"
                elif "step 6" in low or ("mask" in low and "step" in low):
                    section = "masks"
                continue
            if "=" in stripped and not table_row_re.match(stripped):
                key, _, val_raw = stripped.partition("=")
                key = key.strip().lower()
                val = _strip_val(val_raw)
                if key not in _DROPPED:
                    data["kv"][key] = val
                continue
            if section == "bands":
                data["bands"].append(stripped)
            elif section == "cubes":
                data["cubes"].append(stripped)
            elif section == "masks":
                data["masks"].append(stripped)
            else:
                data["unknown"].append(stripped)

    return data


# ---------------------------------------------------------------------------
# Mask row converter
# ---------------------------------------------------------------------------


def _convert_mask_rows(rows):
    input_lines, window_lines, noise_lines = [], [], []

    for row in rows:
        row = re.sub(r"\s*#.*$", "", row).strip()
        if not row:
            continue
        parts = [p.strip() for p in re.split(r"[\t,]+", row) if p.strip()]
        if len(parts) < 3:
            continue

        first = parts[0].lower()

        # Noise velocity ranges
        if first in ("noise_vel", "noise_mask", "noise"):
            if len(parts) >= 5:
                noise_lines.append(f"noise_mask = {parts[0]}, {', '.join(parts[1:5])}")
            continue

        # Distinguish file mask (4 parts) vs velocity window (5 parts)
        if len(parts) >= 5:
            window_lines.append(f"window_mask = {parts[0]}, {', '.join(parts[1:5])}")
        elif len(parts) == 4:
            input_lines.append(f"input_mask = {parts[0]}, {', '.join(parts[1:4])}")

    return {
        "input_mask_lines": input_lines,
        "window_mask_lines": window_lines,
        "noise_mask_lines": noise_lines,
        "has_input": bool(input_lines),
        "has_window": bool(window_lines),
        "has_noise": bool(noise_lines),
    }


# ---------------------------------------------------------------------------
# ref_line token builder
# ---------------------------------------------------------------------------


def _build_ref_line(raw_ref, use_input, use_window):
    warnings = []
    ref = raw_ref.strip().strip("'\"").strip()

    if ref.lower() in ("ref+hi", "ref + hi"):
        warnings.append(
            "ref_line = 'ref+HI' is no longer supported and has been replaced "
            "with 'first'. Review your masking configuration."
        )
        ref = "first"

    tokens = [ref]
    if use_input:
        tokens.append("input")
    if use_window:
        tokens.append("window")

    return ", ".join(tokens), warnings


# ---------------------------------------------------------------------------
# Main converter
# ---------------------------------------------------------------------------


def convert(old_path, new_path, band_list_path=None, cube_list_path=None):
    warnings_out = []
    d = _parse_old_config(old_path)
    kv = d["kv"]

    # Scalar settings
    user = kv.get("user", "")
    comments = kv.get("comments", "")
    data_dir = kv.get("data_dir", "data/")
    out_dir = kv.get("out_dic", kv.get("out_dir", "output/"))
    geom_file = kv.get("geom_file", "")
    hfs_file = kv.get("hfs_file", "")
    overlay_file = kv.get("overlay_file", "")
    folder_savefits = kv.get("folder_savefits", "./saved_fits_files/")
    targets = kv.get("targets", "")
    target_res = kv.get("target_res", "27.0")
    resolution = kv.get("resolution", "angular").strip("'\"")
    pixels_per_beam = kv.get("spacing_per_beam", kv.get("pixels_per_beam", "2"))
    max_rad = kv.get("max_rad", "auto").strip("'\"")
    naxis_shuff = kv.get("naxis_shuff", "200")
    cdelt_shuff = kv.get("cdelt_shuff", "4000.0")
    raw_ref = kv.get("ref_line", "first")
    sn_processing = _convert_sn(kv.get("sn_processing", "2, 4"))
    mom_thresh = kv.get("mom_thresh", "5")
    conseq_ch = kv.get("conseq_channels", "3")
    mom2_method = kv.get("mom2_method", "fwhm").strip("'\"")
    spec_smooth = kv.get("spec_smooth", "default").strip("'\"")
    spec_smooth_method = kv.get("spec_smooth_method", "binned").strip("'\"")
    save_mom_maps = kv.get("save_mom_maps", "true").lower().strip("'\"")
    save_maps = (
        kv.get("save_band_maps", kv.get("save_maps", "true")).lower().strip("'\"")
    )
    structure_creation = kv.get("structure_creation", "default").strip("'\"")
    fname_fill = kv.get("fname_fill", "").strip("'\"")

    # strict_mask: bool → false/strict
    old_strict = kv.get("strict_mask", "false").lower().strip("'\"")
    if old_strict in ("true", "1", "yes"):
        strict_mask = "strict"
        warnings_out.append(
            "strict_mask = True converted to strict_mask = strict. "
            "New options: false | strict | broad."
        )
    else:
        strict_mask = "false"

    # Old boolean flags → ref_line tokens
    def _is_true(v):
        return v.lower().strip("'\"") in ("true", "1", "yes")

    use_input_mask = _is_true(kv.get("use_input_mask", "false"))
    use_fixed_vel = _is_true(kv.get("use_fixed_vel_mask", "false"))
    use_noise = _is_true(
        kv.get("use_noise_vel_ranges", kv.get("use_fixed_noise_mask", "false"))
    )
    use_hfs_lines = kv.get("use_hfs_lines", "false").lower().strip("'\"")

    if "save_fits" in kv:
        warnings_out.append(
            "save_fits is no longer supported. "
            "Use save_cubes = true in [output] to save convolved cube FITS files."
        )

    # Build new ref_line
    ref_line, ref_warns = _build_ref_line(raw_ref, use_input_mask, use_fixed_vel)
    warnings_out.extend(ref_warns)

    # Mask table
    masks = _convert_mask_rows(d["masks"])
    if use_input_mask and not masks["has_input"]:
        warnings_out.append(
            "use_input_mask = True but no file mask row found in Step 6. "
            "Add an input_mask = ... row to the [mask] table."
        )
    if use_fixed_vel and not masks["has_window"]:
        warnings_out.append(
            "use_fixed_vel_mask = True but no velocity-window row found in Step 6. "
            "Add a window_mask = ... row to the [mask] table."
        )
    if use_noise and not masks["has_noise"]:
        warnings_out.append(
            "use_noise_vel_ranges = True but no noise_vel row found in Step 6. "
            "Add noise_mask = ... row(s) to the [mask] table."
        )

    # Band and cube tables
    band_rows_raw = _read_list_file(band_list_path) if band_list_path else d["bands"]
    cube_rows_raw = _read_list_file(cube_list_path) if cube_list_path else d["cubes"]

    map_lines = []
    for row in _parse_table_rows(band_rows_raw):
        while len(row) < 6:
            row.append("")
        map_lines.append(
            f"{row[0]},  {row[1]},  {row[2]},  {row[3]},  {row[4]},  {row[5]}"
        )

    cube_lines = []
    for row in _parse_table_rows(cube_rows_raw):
        while len(row) < 7:
            row.append("")
        cube_lines.append(
            f"{row[0]},  {row[1]},  {row[2]},  {row[3]},  {row[4]},  {row[5]},  {row[6]}"
        )

    # Optional path lines
    _default_geom = ("keys/target_definitions.txt", "List_Files/geometry.txt", "")
    geom_line = (
        f"geom_file        = {geom_file}"
        if geom_file and geom_file not in _default_geom
        else "# geom_file        = keys/target_definitions.txt  # (default)"
    )
    _default_hfs = ("keys/hfs_lines.txt", "List_Files/hfs_lines.txt", "")
    hfs_line = (
        f"hfs_file         = {hfs_file}"
        if hfs_file and hfs_file not in _default_hfs
        else "# hfs_file         = keys/hfs_lines.txt  # (uncomment if needed)"
    )
    fname_fill_line = (
        f"fname_fill       = {fname_fill}"
        if fname_fill
        else "# fname_fill       = <filename>.ecsv"
    )

    warning_block = ""
    if warnings_out:
        warning_block = (
            "\n# ---- CONVERSION WARNINGS (review before running) ----\n"
            + "".join(f"# [WARN] {w}\n" for w in warnings_out)
            + "# -------------------------------------------------------\n"
        )

    # Assemble output
    sections = []

    sections.append(
        f"# =============================================================================\n"
        f"# HexMaps config.txt  (converted from {old_path.name})\n"
        f"# =============================================================================\n"
        f"# Source: PyStructure v4.x  (PhangsTeam/astro-HexMaps PyStructure_v4p2)\n"
        f"# Target: HexMaps v5+       (lukas-neumann-astro/astro-HexMaps main)\n"
        f"# =============================================================================\n"
        f"{warning_block}\n"
        f"[meta]\n"
        f"user             = {user}\n"
        f"comments         = {comments}\n"
        f"\n"
        f"[paths]\n"
        f"data_dir         = {data_dir}\n"
        f"out_dir          = {out_dir}\n"
        f"{geom_line}\n"
        f"{hfs_line}\n"
        f"folder_savefits  = {folder_savefits}\n"
        f"\n"
        f"[targets]\n"
        f"targets          = {targets}\n"
        f"\n"
        f"[overlay]\n"
        f"overlay_file     = {overlay_file}\n"
        f"\n"
        f"# ---- maps ----"
    )

    for ml in map_lines:
        sections.append(ml)
    if not map_lines:
        sections.append("# (no band/map entries found — add rows here)")

    sections.append("\n# ---- cubes ----")
    for cl in cube_lines:
        sections.append(cl)
    if not cube_lines:
        sections.append("# (no cube entries found — add rows here)")

    sections.append("\n# ---- mask ----")
    for row in masks["input_mask_lines"]:
        sections.append(row)
    for row in masks["window_mask_lines"]:
        sections.append(row)
    for row in masks["noise_mask_lines"]:
        sections.append(row)
    if not (
        masks["input_mask_lines"]
        or masks["window_mask_lines"]
        or masks["noise_mask_lines"]
    ):
        sections.append(
            "# (no mask rows found — add input_mask, window_mask, or noise_mask rows here)"
        )

    sections.append(
        f"\n\n[resolution]\n"
        f"target_res       = {target_res}\n"
        f"resolution       = {resolution}\n"
        f"pixels_per_beam  = {pixels_per_beam}\n"
        f"max_rad          = {max_rad}\n"
        f"NAXIS_shuff      = {naxis_shuff}\n"
        f"CDELT_SHUFF      = {cdelt_shuff}\n"
        f"\n"
        f"[masking]\n"
        f"# ref_line tokens: first | all | <n> | <LINE_NAME> | individual\n"
        f"#   optional: input, window  (include external masks)\n"
        f"#   optional: AND | OR       (combinator; default OR)\n"
        f"ref_line              = {ref_line}\n"
        f"SN_processing         = {sn_processing}\n"
        f"# strict_mask: false | strict | broad\n"
        f"strict_mask           = {strict_mask}\n"
        f"use_fixed_noise_mask  = {'true' if use_noise else 'false'}\n"
        f"use_hfs_lines         = {use_hfs_lines}\n"
        f"fov_erosion_beams     = 0.5\n"
        f"mom_thresh            = {mom_thresh}\n"
        f"conseq_channels       = {conseq_ch}\n"
        f"mom2_method           = {mom2_method}\n"
        f"\n"
        f"[spectral]\n"
        f"spec_smooth           = {spec_smooth}\n"
        f"spec_smooth_method    = {spec_smooth_method}\n"
        f"\n"
        f"[output]\n"
        f"save_cubes            = false\n"
        f"save_mom_maps         = {save_mom_maps}\n"
        f"save_maps             = {save_maps}\n"
        f"save_mask             = false\n"
        f"\n"
        f"[structure]\n"
        f"structure_creation    = {structure_creation}\n"
        f"{fname_fill_line}"
    )

    if d["unknown"]:
        sections.append(
            "\n# ---- unrecognised lines from old config (review manually) ----"
        )
        for u in d["unknown"]:
            sections.append(f"# {u}")

    new_path.write_text("\n".join(sections) + "\n", encoding="utf-8")
    print(f"[OK] Written: {new_path}")
    if warnings_out:
        print(
            f"\n[WARN] {len(warnings_out)} conversion warning(s) — "
            "see comment block at the top of the output file."
        )
        for w in warnings_out:
            print(f"  • {w}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Convert a PyStructure v4.x PyStructure.conf to a HexMaps config.txt."
    )
    parser.add_argument("old_conf", help="Path to the old PyStructure.conf")
    parser.add_argument("new_conf", help="Path for the new HexMaps config.txt")
    parser.add_argument(
        "--band-list",
        metavar="FILE",
        help="Path to an external band_list.txt (overrides inline Step 4 rows)",
    )
    parser.add_argument(
        "--cube-list",
        metavar="FILE",
        help="Path to an external cube_list.txt (overrides inline Step 5 rows)",
    )
    args = parser.parse_args()

    old_path = Path(args.old_conf)
    new_path = Path(args.new_conf)
    if not old_path.exists():
        print(f"[ERROR] Input file not found: {old_path}")
        sys.exit(1)

    convert(
        old_path,
        new_path,
        Path(args.band_list) if args.band_list else None,
        Path(args.cube_list) if args.cube_list else None,
    )


if __name__ == "__main__":
    main()
