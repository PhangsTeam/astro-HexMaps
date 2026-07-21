#!/usr/bin/env python3
"""
target_definitions_conversion.py — Convert the old PyStructure geometry.txt
(tab-separated, 11 columns + optional extras) to the new HexMaps
target_definitions.txt (comma-separated, same 11 columns, with a header).

Usage:
    python target_definitions_conversion.py geometry.txt target_definitions.txt

Old format (PhangsTeam/PyStructure  List_Files/geometry.txt):
    Tab-separated, no header row, columns:
      1  source name
      2  RA (deg)
      3  Dec (deg)
      4  dist_mpc
      5  e_dist_mpc
      6  incl_deg
      7  e_incl_deg
      8  posang_deg
      9  e_posang_deg
      10 r25 (arcmin)
      11 e_r25

New format (keys/target_definitions.txt):
    Comma-separated with a comment-header, same 11 columns.
    Spaces and tabs around commas are ignored by the parser.
"""

import re
import sys
from pathlib import Path

HEADER = """\
# =============================================================================
# HexMaps target_definitions.txt  (converted from {src})
# =============================================================================
# Comma-separated table of target parameters.
# Spaces and tabs around each comma are ignored; so feel free to align
# columns with extra whitespace for readability.
# All targets that may ever be used should be listed here.
# The targets actually processed are defined in config.txt [targets].
#
# Columns (comma-separated; no header row):
# [required]
#   target       - target name (must match the FITS filename prefix)
#   x_ctr        - x-coordinate of target centre (e.g. RA in degrees; read from FITS header of overlay file)
#   y_ctr        - y-coordinate of target centre (e.g. Dec in degrees; read from FITS header of overlay file)
#   dist_mpc     - Distance (Mpc)
#   e_dist_mpc   - Uncertainty of distance (Mpc)
# [optional]
#   incl_deg     - Inclination (degrees)
#   e_incl_deg   - Uncertainty of inclination (degrees)
#   posang_deg   - Position angle (degrees; East of North)
#   e_posang_deg - Uncertainty of position angle (degrees)
#   r25          - Optical radius r25 (arcmin)
#   e_r25        - Uncertainty of r25 (arcmin)
# =============================================================================
# target,       x_ctr,      y_ctr,      dist_mpc,   e_dist_mpc,   incl_deg,   e_incl_deg,   posang_deg,   e_posang_deg,   r25,     e_r25
"""


def convert(old_path: Path, new_path: Path):
    rows = []
    skipped = []

    with open(old_path, encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n").strip()
            if not line or line.startswith("#"):
                continue

            # Split on any mix of tabs and spaces (old files use single tab,
            # but some have extra trailing spaces)
            parts = re.split(r"\t+", line)
            # Also handle space-only separated files
            if len(parts) < 4:
                parts = line.split()

            if len(parts) < 11:
                skipped.append(line)
                continue

            # Pad to exactly 11 columns, clean each value
            cols = [p.strip() for p in parts[:11]]
            while len(cols) < 11:
                cols.append("NaN")

            # Column 12 (if present) is a literature reference — not used by
            # HexMaps but we note it in a trailing comment so it is not lost.
            ref_comment = ""
            if len(parts) >= 12:
                ref_comment = "  # ref: " + " ".join(p.strip() for p in parts[11:])

            # Format: align target name left, numbers right-padded for readability
            target = cols[0]
            nums = cols[1:]
            row = f"{target:<12}, " + ", ".join(f"{v:>12}" for v in nums) + ref_comment
            rows.append(row)

    header = HEADER.format(src=old_path.name)
    body = "\n".join(rows) + "\n"

    new_path.write_text(header + body, encoding="utf-8")
    print(f"[OK] {len(rows)} target(s) written to: {new_path}")

    if skipped:
        print(f"[WARN] {len(skipped)} line(s) skipped (fewer than 11 columns):")
        for s in skipped:
            print(f"       {s}")


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: python target_definitions_conversion.py "
            "<old_geometry.txt> <new_target_definitions.txt>"
        )
        sys.exit(1)
    old_path = Path(sys.argv[1])
    new_path = Path(sys.argv[2])
    if not old_path.exists():
        print(f"[ERROR] Input file not found: {old_path}")
        sys.exit(1)
    convert(old_path, new_path)


if __name__ == "__main__":
    main()
