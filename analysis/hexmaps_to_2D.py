#!/usr/bin/env python
"""
hexmaps_to_2D.py: convert a HexMaps .ecsv database to its 2D-only version
(i.e. with all SPEC_ columns/metadata removed).

Usage
-----
    python hexmaps_to_2D.py <fname.ecsv>

This loads the given .ecsv file with HexMapsAnalysis (from the hexmaps
package), strips out the 3D/spectral (SPEC_) columns via
`get_2D_database`, and saves the result as `<fname>_2D.ecsv`.
"""

__author__ = "J. den Brok & L. Neumann"

import argparse
import os

from hexmaps import HexMapsAnalysis


def main():
    parser = argparse.ArgumentParser(
        description="Convert a HexMaps .ecsv database to its 2D-only version."
    )
    parser.add_argument(
        "fname",
        type=str,
        help="Path to the input .ecsv file produced by the HexMaps pipeline.",
    )
    args = parser.parse_args()

    fname = args.fname
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"Input file not found: {fname}")

    hm = HexMapsAnalysis(fname)

    base, ext = os.path.splitext(fname)
    out_fname = f"{base}_2D{ext}"

    hm.get_2D_database(fname=out_fname, save=True)


if __name__ == "__main__":
    main()