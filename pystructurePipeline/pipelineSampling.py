"""
pipelineSampling.py — generate the hexagonal sampling grid for a source.

The sampling grid is the spatial backbone of the PyStructure database.
Every subsequent stage (regrid, spectra, output) works at these positions.

Design
------
The grid is hexagonal (close-packed circles) rather than rectangular because
it provides more uniform spatial coverage with the least number of beams and
minimises the number of correlated positions for a given beam spacing.

The grid is generated in RA/Dec space, clipped to the footprint of the overlay
cube (pixels with at least one finite channel), and the spacing is derived from
the target resolution and the spacing_per_beam parameter in config_key.txt.
"""

import numpy as np
from astropy.io import fits
from os import path

from pystructurePipeline.utilsFits import twod_head, make_sampling_points

from pystructurePipeline.pystructureLogger import get_logger

_log = get_logger("Sampling")



def run_sampling(source: str, params: dict, meta: dict) -> dict:
    """
    Generate the hexagonal sampling grid for *source*.

    Steps
    -----
    1. Load the overlay FITS cube and check it is 3-D.
    2. Determine the target resolution in arcseconds (angular, physical, or
       native mode).
    3. Collapse the cube along the spectral axis to create a binary footprint
       mask (True where at least one channel is finite).
    4. Build the hexagonal grid centred on the source and clip it to the mask.

    Parameters
    ----------
    source : str
        Source name; used to construct the overlay filename as
        ``{data_dir}/{source}{overlay_file}``.
    params : dict
        Source geometric parameters from SourceHandler.get_source_params().
        Required keys: ra_ctr, dec_ctr, dist_mpc.
    meta : dict
        Pipeline settings from KeyHandler.meta.
        Used keys: data_dir, overlay_file, resolution, target_res,
        spacing_per_beam, max_rad.

    Returns
    -------
    dict with keys:
        samp_ra       : np.ndarray — RA of each sampling point (degrees)
        samp_dec      : np.ndarray — Dec of each sampling point (degrees)
        ov_hdr        : astropy FITS Header — header of the overlay cube
        mask_hdr      : astropy FITS Header — 2-D version of ov_hdr
        target_res_as : float — target resolution in arcseconds

    Raises
    ------
    FileNotFoundError if the overlay FITS file does not exist.
    ValueError if the overlay cube is 4-D.
    """
    data_dir     = meta.get("data_dir", "data/")
    overlay_file = meta.get("overlay_file", "")

    # Construct the overlay filename: if the source name is already embedded in
    # overlay_file use it as-is, otherwise prepend the source name.
    overlay_fname = (path.join(data_dir, overlay_file)
                     if source in overlay_file
                     else path.join(data_dir, source + overlay_file))

    if not path.exists(overlay_fname):
        _log.error(
            f"Overlay file not found for {source}: {overlay_fname}"
        )
        raise FileNotFoundError(
            f"Overlay file not found for {source}: {overlay_fname}"
        )

    ov_cube, ov_hdr = fits.getdata(overlay_fname, header=True)

    if ov_hdr["NAXIS"] == 4:
        _log.error(
            f"4D overlay cube for {source}. "
            "Please provide a 3D cube."
        )
        raise ValueError(
            f"4D overlay cube for {source}. "
            "Please provide a 3D cube."
        )

    # ------------------------------------------------------------------
    # Determine target resolution in arcseconds
    # ------------------------------------------------------------------
    resolution = meta.get("resolution", "angular")
    target_res = meta.get("target_res", 27.0)

    if resolution == "native":
        # Use the native beam of the overlay cube
        target_res_as = max(ov_hdr.get("BMIN", 0), ov_hdr.get("BMAJ", 0)) * 3600.0
        _log.info(f"Native resolution: {target_res_as:.1f} arcsec.")
    elif resolution == "physical":
        # Convert target_res (parsecs) to arcseconds using the source distance
        dist_mpc = params.get("dist_mpc", 1.0)
        target_res_as = 3600.0 * 180.0 / np.pi * 1e-6 * float(target_res) / dist_mpc
        _log.info(f"Physical resolution: {target_res} pc "
                  f"= {target_res_as:.1f} arcsec at {dist_mpc} Mpc.")
    else:
        # Angular: use target_res directly in arcseconds
        target_res_as = float(target_res)
        _log.info(f"Angular resolution: {target_res_as:.1f} arcsec.")

    # ------------------------------------------------------------------
    # Build the footprint mask and generate the hex grid
    # ------------------------------------------------------------------
    # The mask is True wherever at least one spectral channel is finite.
    # This clips the grid to the mapped area without relying on a separate mask file.
    mask     = np.sum(np.isfinite(ov_cube), axis=0) >= 1
    mask_hdr = twod_head(ov_hdr)

    spacing_per_beam = meta.get("spacing_per_beam", 2.0)
    max_rad          = meta.get("max_rad", "auto")
    # Spacing in degrees: one beam FWHM divided by spacing_per_beam
    spacing = target_res_as / 3600.0 / float(spacing_per_beam)

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

    _log.info(f"Hexagonal grid generated: "
              f"{len(samp_ra)} sampling points "
              f"(spacing = {spacing * 3600:.1f} arcsec).")

    return dict(
        samp_ra       = samp_ra,
        samp_dec      = samp_dec,
        ov_hdr        = ov_hdr,
        mask_hdr      = mask_hdr,
        target_res_as = target_res_as,
    )
