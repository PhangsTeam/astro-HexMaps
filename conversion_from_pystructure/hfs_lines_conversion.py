#!/usr/bin/env python3
"""
hfs_lines_conversion.py — Convert the old PyStructure hfs_lines.txt
(tab-separated) to the new HexMaps hfs_lines.txt (comma-separated).

Usage:
    python hfs_lines_conversion.py old_hfs_lines.txt new_hfs_lines.txt

Old format (PhangsTeam/PyStructure  List_Files/hfs_lines.txt):
    Tab-separated, no header row, columns:
      1  line name  (must match database entry)
      2  reference frequency
      3  hyperfine transition frequency
      4  unit (astropy.units readable)

New format (keys/hfs_lines.txt):
    Comma-separated with a comment-header, same 4 columns.
    Spaces and tabs around commas are ignored by the parser.
"""

import re
import sys
from pathlib import Path

HEADER = """\
# =============================================================================
# HexMaps hfs_lines.txt  (converted from {src})
# =============================================================================
# Comma-separated hyperfine structure line definitions.
# Spaces and tabs around each comma are ignored by the parser.
#
# Columns (comma-separated; no header row):
#   line      - line name (must match the cube short name in config.txt)
#   ref_freq  - reference (main) frequency
#   hfs_freq  - hyperfine transition frequency
#   unit      - frequency unit (astropy.units readable, e.g. GHz)
# =============================================================================
# line,      ref_freq,         hfs_freq,         unit
"""


def convert(old_path: Path, new_path: Path):
    rows = []
    skipped = []

    with open(old_path, encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n").strip()
            if not line or line.startswith("#"):
                continue

            # Split on tabs or commas (handle already-comma-separated files)
            parts = re.split(r"[\t,]+", line)
            parts = [p.strip() for p in parts if p.strip()]

            if len(parts) < 4:
                skipped.append(line)
                continue

            line_name, ref_freq, hfs_freq, unit = (
                parts[0], parts[1], parts[2], parts[3]
            )
            row = (
                f"{line_name:<10}, "
                f"{ref_freq:>16}, "
                f"{hfs_freq:>16}, "
                f"    {unit}"
            )
            rows.append(row)

    header = HEADER.format(src=old_path.name)
    body = "\n".join(rows) + "\n"

    new_path.write_text(header + body, encoding="utf-8")
    print(f"[OK] {len(rows)} row(s) written to: {new_path}")

    if skipped:
        print(f"[WARN] {len(skipped)} line(s) skipped (fewer than 4 columns):")
        for s in skipped:
            print(f"       {s}")


def main():
    if len(sys.argv) != 3:
        print("Usage: python hfs_lines_conversion.py "
              "<old_hfs_lines.txt> <new_hfs_lines.txt>")
        sys.exit(1)
    old_path = Path(sys.argv[1])
    new_path = Path(sys.argv[2])
    if not old_path.exists():
        print(f"[ERROR] Input file not found: {old_path}")
        sys.exit(1)
    convert(old_path, new_path)


if __name__ == "__main__":
    main()
