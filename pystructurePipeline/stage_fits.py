"""
stage_fits.py — write FITS moment maps and 2D map images for a source.

This stage converts the hex-sampled data stored in the .ecsv table back into
regular FITS images on a rectangular pixel grid for use in standard image
viewers (DS9, CARTA) and for sharing with collaborators who do not use PyStructure.

Workflow
--------
1. Load the overlay cube to obtain a reference WCS and footprint mask.
2. For each moment/quantity column, regrid the hex-sampled values onto a
   rectangular grid using nearest-neighbour interpolation via scipy.griddata.
3. Multiply by the footprint mask so that pixels outside the mapped area
   are set to NaN.
4. Optionally resample to a coarser pixel scale matching the beam (1/3 beam
   per pixel) if the overlay has a finer native pixel scale.
5. Write one FITS file per quantity.

Pixel scale of output FITS files
---------------------------------
The output pixel scale is set to target_res / 3 arcsec/pixel, giving three
pixels per beam FWHM.  This is sufficient for a clean image while keeping
file sizes manageable.  If the overlay WCS already has a coarser pixel scale,
no resampling is done.

Output filename convention
--------------------------
{source}_{line}_{quantity}.fits

e.g.  ngc5194_12CO21_mom0.fits
      ngc5194_SPIRE250_map.fits
"""

import os
import copy
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from scipy.interpolate import griddata
from reproject import reproject_interp

from pystructurePipeline.utils_fits import twod_head

from pystructurePipeline.pystructureLogger import get_logger

LOG = get_logger("FITS")



# ============================================================================
# Grid helpers
# ============================================================================

def sample_to_hdr(in_data, ra_samp, dec_samp, in_hdr):
    """
    Regrid hex-sampled 1-D data onto a 2-D rectangular pixel grid.

    Uses scipy.griddata with nearest-neighbour interpolation.  The sampling
    points are first converted to pixel coordinates using the input WCS, and
    the output grid covers the full NAXIS1 × NAXIS2 extent of *in_hdr*.

    Parameters
    ----------
    in_data   : array-like (n_pts,)  — values at the hex-grid positions
    ra_samp   : array-like (n_pts,)  — RA of hex-grid points (degrees)
    dec_samp  : array-like (n_pts,)  — Dec of hex-grid points (degrees)
    in_hdr    : FITS Header           — 2-D WCS defining the output grid

    Returns
    -------
    gridded : np.ndarray (NAXIS2, NAXIS1) — regridded map
    """
    x_axis = np.arange(in_hdr["NAXIS1"])
    y_axis = np.arange(in_hdr["NAXIS2"])
    grid_x, grid_y = np.meshgrid(x_axis, y_axis)

    wcs          = WCS(in_hdr)
    pixel_coords = wcs.all_world2pix(np.column_stack((ra_samp, dec_samp)), 0)
    return griddata(pixel_coords, in_data, (grid_x, grid_y), method="nearest")


def resample_hdr(hdr_ov, target_res):
    """
    Build a new 2-D WCS header with a pixel scale of target_res / 3 arcsec/px.

    The reference pixel is placed at the bottom-left corner of the overlay
    footprint.  The output grid size is computed so that it covers the same
    sky area as the input overlay at the new (coarser) pixel scale.

    Parameters
    ----------
    hdr_ov     : FITS Header — 2-D overlay header
    target_res : float       — target beam FWHM in arcseconds

    Returns
    -------
    hdr_new : FITS Header — 2-D WCS header for the output FITS image
    """
    wcs_new = WCS(naxis=2)
    wcs_new.wcs.crpix = [1, 1]
    wcs_ov = WCS(hdr_ov)
    ra_ref, dec_ref = wcs_ov.all_pix2world(0, 0, 0)
    wcs_new.wcs.crval = [ra_ref, dec_ref]
    wcs_new.wcs.cunit = ["deg", "deg"]
    wcs_new.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    delta_px = target_res / 3600.0 / 3.0   # 3 pixels per beam FWHM
    wcs_new.wcs.cdelt = [-delta_px, delta_px]

    xaxis_n = int(np.round(hdr_ov["NAXIS1"] * abs(hdr_ov["CDELT1"]) / delta_px))
    yaxis_n = int(np.round(hdr_ov["NAXIS2"] * abs(hdr_ov["CDELT2"]) / delta_px))
    wcs_new.array_shape = [xaxis_n, yaxis_n]

    hdr_new          = wcs_new.to_header()
    hdr_new["NAXIS"]  = 2
    hdr_new["NAXIS1"] = xaxis_n
    hdr_new["NAXIS2"] = yaxis_n
    hdr_new["BMAJ"]   = target_res / 3600.0
    hdr_new["BMIN"]   = target_res / 3600.0
    hdr_new["BPA"]    = 0.0
    return hdr_new


def save_to_fits(ra, dec, hdr_in, ov_slice, key, filename,
                 this_source, this_data, line, folder, target_res):
    """
    Regrid one table column and write it to a FITS file.

    Silently skips if the requested column does not exist in *this_data*
    (e.g. when EMOM0 has not been computed because all pixels are undetected).

    Parameters
    ----------
    ra, dec    : arrays       — hex-grid RA/Dec (degrees)
    hdr_in     : FITS Header  — 2-D overlay header (pixel grid definition)
    ov_slice   : np.ndarray   — footprint mask (1.0 inside mapped area, NaN outside)
    key        : str          — column prefix, e.g. "MOM0", "TPEAK", "MAP_"
    filename   : str          — quantity label used in the output filename, e.g. "mom0"
    this_source: str          — source name
    this_data  : Table        — the PyStructure table
    line       : str          — line/map name, e.g. "12CO21" or "SPIRE250"
    folder     : str          — output directory
    target_res : float        — target beam FWHM in arcseconds (for header)
    """
    col_name = f"{key}_{line.upper()}"
    if col_name not in this_data.colnames:
        return

    data_in  = copy.deepcopy(this_data[col_name])
    map_cart = sample_to_hdr(data_in, ra, dec, hdr_in)

    # Apply footprint mask (NaN outside the observed area)
    map_cart = ov_slice * map_cart

    # Resample to a coarser grid if the overlay pixel scale is finer than the beam
    native_beam_as = 3600.0 * min(hdr_in.get("BMAJ", 1e6), hdr_in.get("BMIN", 1e6))
    if native_beam_as < 0.99 * target_res:
        hdr_repr = resample_hdr(hdr_in, target_res)
        map_cart, _ = reproject_interp((map_cart, hdr_in), hdr_repr)
        hdr_in = hdr_repr

    fname_fits = os.path.join(folder, f"{this_source}_{line}_{filename}.fits")
    fits.writeto(fname_fits, data=map_cart, header=hdr_in, overwrite=True)


# ============================================================================
# Stage entry point
# ============================================================================

def run_fits(source, fname, meta, maps, cubes, params):
    """
    Write FITS moment maps and 2D map images for *source*.

    This is the entry point for the "fits" pipeline stage.

    Reads the processed .ecsv table and writes one FITS file per moment quantity
    per line (if save_mom_maps is True) and one FITS file per 2D map (if
    save_maps is True).

    Parameters
    ----------
    source : str
    fname  : str          — path to the processed .ecsv file
    meta   : dict         — from KeyHandler.meta
    maps   : pd.DataFrame — map definitions from KeyHandler
    cubes  : pd.DataFrame — cube definitions from KeyHandler
    params : dict         — source geometry from SourceHandler
    """
    save_mom_maps    = meta.get("save_mom_maps",    True)
    save_maps        = meta.get("save_maps",        True)
    folder           = meta.get("folder_savefits",  "./saved_FITS_files/")
    target_res_as    = _resolve_target_res(params, meta)
    spacing_per_beam = meta.get("spacing_per_beam", 2.0)

    if not (save_mom_maps or save_maps):
        LOG.info(f"Output writing disabled for {source}; skipping.")
        return

    # Warn if spacing is too coarse for good image reconstruction
    if float(spacing_per_beam) < 4:
        LOG.warning(f"Spacing_per_beam = {spacing_per_beam} < 4. "
                    "The output FITS images may show hexagonal grid artefacts. "
                    "Consider using spacing_per_beam ≥ 4 for publication-quality output maps.")

    os.makedirs(folder, exist_ok=True)

    # ------------------------------------------------------------------
    # Load overlay cube to get the reference WCS and footprint mask
    # ------------------------------------------------------------------
    data_dir     = meta.get("data_dir", "data/")
    overlay_file = meta.get("overlay_file", "")
    from os import path as _path
    overlay_fname = (_path.join(data_dir, overlay_file) if source in overlay_file
                     else _path.join(data_dir, source + overlay_file))

    ov_cube, ov_hdr = fits.getdata(overlay_fname, header=True)

    # Build a binary footprint mask from the middle channel of the overlay cube.
    # Pixels outside the observed area are set to NaN after gridding.
    ov_slice = ov_cube[ov_hdr["NAXIS3"] // 2, :, :].copy()
    ov_slice[np.isfinite(ov_slice)] = 1.0
    ov_hdr_2d = twod_head(ov_hdr)

    this_data = Table.read(fname)
    ra_deg    = this_data["ra_deg"]
    dec_deg   = this_data["dec_deg"]

    # ------------------------------------------------------------------
    # Moment maps
    # ------------------------------------------------------------------
    if save_mom_maps:
        for line in cubes["line_name"]:
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "MOM0",  "mom0",  source, this_data, line, folder, target_res_as)
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "EMOM0", "emom0", source, this_data, line, folder, target_res_as)
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "MOM1",  "mom1",  source, this_data, line, folder, target_res_as)
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "EMOM1", "emom1", source, this_data, line, folder, target_res_as)
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "MOM2",  "mom2",  source, this_data, line, folder, target_res_as)
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "EMOM2", "emom2", source, this_data, line, folder, target_res_as)
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "TPEAK", "tpeak", source, this_data, line, folder, target_res_as)
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "RMS",   "rms",   source, this_data, line, folder, target_res_as)
        LOG.info(f"Moment map FITS files written to: {folder}")

    # ------------------------------------------------------------------
    # 2D map images
    # ------------------------------------------------------------------
    if save_maps:
        for map_name in maps["map_name"]:
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "MAP_",  "map",  source, this_data, map_name, folder, target_res_as)
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "EMAP_", "emap", source, this_data, map_name, folder, target_res_as)
        LOG.info(f"2D map FITS files written to: {folder}")


def _resolve_target_res(params, meta):
    """
    Convert the configured target resolution to arcseconds.

    In physical mode, converts the target resolution from parsecs to
    arcseconds using the source distance stored in *params*.
    """
    resolution = meta.get("resolution", "angular")
    target_res = float(meta.get("target_res", 27.0))
    if resolution == "physical":
        dist_mpc = params.get("dist_mpc", 1.0)
        return 3600.0 * 180.0 / np.pi * 1e-6 * target_res / dist_mpc
    return target_res
