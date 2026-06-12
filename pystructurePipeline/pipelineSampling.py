"""
stage_sampling: generate hexagonal sampling grid for a source.

All logic is contained here — no imports from the legacy scripts/ directory.
"""

import numpy as np
from astropy.io import fits
from os import path

from pystructurePipeline.utilsFits import twod_head, make_sampling_points


def run_sampling(source: str, params: dict, meta: dict) -> dict:
    """
    Generate hexagonal sampling grid for *source*.

    Parameters
    ----------
    source : str
    params : dict  — from SourceHandler.get_source_params()
    meta   : dict  — from KeyHandler.meta

    Returns
    -------
    dict with keys:
        samp_ra, samp_dec, ov_hdr, mask_hdr, target_res_as
    """
    data_dir     = meta.get("data_dir", "data/")
    overlay_file = meta.get("overlay_file", "")

    overlay_fname = (path.join(data_dir, overlay_file)
                     if source in overlay_file
                     else path.join(data_dir, source + overlay_file))

    if not path.exists(overlay_fname):
        raise FileNotFoundError(
            f"Overlay file not found for {source}: {overlay_fname}")

    ov_cube, ov_hdr = fits.getdata(overlay_fname, header=True)

    if ov_hdr["NAXIS"] == 4:
        raise ValueError(f"4D overlay cube for {source}. Need 3D.")

    # Resolve target resolution
    resolution = meta.get("resolution", "angular")
    target_res = meta.get("target_res", 27.0)

    if resolution == "native":
        target_res_as = max(ov_hdr.get("BMIN", 0), ov_hdr.get("BMAJ", 0)) * 3600.0
    elif resolution == "physical":
        dist_mpc = params.get("dist_mpc", 1.0)
        target_res_as = 3600.0 * 180.0 / np.pi * 1e-6 * float(target_res) / dist_mpc
    else:
        target_res_as = float(target_res)

    mask     = np.sum(np.isfinite(ov_cube), axis=0) >= 1
    mask_hdr = twod_head(ov_hdr)

    spacing_per_beam = meta.get("spacing_per_beam", 2.0)
    max_rad          = meta.get("max_rad", "auto")
    spacing          = target_res_as / 3600.0 / float(spacing_per_beam)

    samp_ra, samp_dec = make_sampling_points(
        ra_ctr     = params["ra_ctr"],
        dec_ctr    = params["dec_ctr"],
        max_rad    = max_rad,
        spacing    = spacing,
        mask       = mask,
        hdr_mask   = mask_hdr,
        overlay_in = overlay_fname,
        show       = False,
    )

    print(f'{"[INFO]":<10} Hexagonal grid: {len(samp_ra)} sampling points.')

    return dict(
        samp_ra       = samp_ra,
        samp_dec      = samp_dec,
        ov_hdr        = ov_hdr,
        mask_hdr      = mask_hdr,
        target_res_as = target_res_as,
    )
