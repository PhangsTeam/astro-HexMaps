"""
stage_output: write FITS moment maps and 2D map FITS files for a source.

All logic is contained here — no imports from the legacy scripts/ directory.
"""

import os
import copy
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from scipy.interpolate import griddata
from reproject import reproject_interp

from pystructurePipeline.utilsFits import twod_head


# ============================================================================
# Grid helpers
# ============================================================================

def sample_to_hdr(in_data, ra_samp, dec_samp, in_hdr):
    """
    Regrid hexagonal-sampled data onto a rectangular pixel grid.

    Port of save_moment_maps.sample_to_hdr (L. Neumann / J. den Brok).
    """
    x_axis    = np.arange(in_hdr["NAXIS1"])
    y_axis    = np.arange(in_hdr["NAXIS2"])
    grid_x, grid_y = np.meshgrid(x_axis, y_axis)

    wcs = WCS(in_hdr)
    pixel_coords = wcs.all_world2pix(np.column_stack((ra_samp, dec_samp)), 0)
    return griddata(pixel_coords, in_data, (grid_x, grid_y), method="nearest")


def resample_hdr(hdr_ov, target_res):
    """
    Build a new WCS header at 1/3 beam pixel scale for FITS output.

    Port of save_moment_maps.resample_hdr.
    """
    wcs_new = WCS(naxis=2)
    wcs_new.wcs.crpix = [1, 1]
    wcs_ov = WCS(hdr_ov)
    ra_ref, dec_ref = wcs_ov.all_pix2world(0, 0, 0)
    wcs_new.wcs.crval  = [ra_ref, dec_ref]
    wcs_new.wcs.cunit  = ["deg", "deg"]
    wcs_new.wcs.ctype  = ["RA---TAN", "DEC--TAN"]

    delta_px = target_res / 3600.0 / 3.0
    wcs_new.wcs.cdelt  = [-delta_px, delta_px]

    xaxis_n = int(np.round(hdr_ov["NAXIS1"] * abs(hdr_ov["CDELT1"]) / delta_px))
    yaxis_n = int(np.round(hdr_ov["NAXIS2"] * abs(hdr_ov["CDELT2"]) / delta_px))
    wcs_new.array_shape = [xaxis_n, yaxis_n]

    hdr_new = wcs_new.to_header()
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
    Regrid one map column and write to FITS.

    Port of save_moment_maps.save_to_fits.
    """
    col_name = f"{key}_{line.upper()}"
    if col_name not in this_data.colnames:
        return

    data_in = copy.deepcopy(this_data[col_name])
    map_cart = sample_to_hdr(data_in, ra, dec, hdr_in)

    # Mask to data footprint
    map_cart = ov_slice * map_cart

    # Resample if overlay resolution is finer than target
    if 3600.0 * min(hdr_in.get("BMAJ", 1e6), hdr_in.get("BMIN", 1e6)) < 0.99 * target_res:
        hdr_repr = resample_hdr(hdr_in, target_res)
        map_cart, _ = reproject_interp((map_cart, hdr_in), hdr_repr)
        hdr_in = hdr_repr

    fname_fits = os.path.join(folder, f"{this_source}_{line}_{filename}.fits")
    fits.writeto(fname_fits, data=map_cart, header=hdr_in, overwrite=True)


# ============================================================================
# Stage entry point
# ============================================================================

def run_output(source, fname, meta, maps, cubes, params):
    """
    Write FITS moment maps and 2D map FITS files for *source*.

    Parameters
    ----------
    source  : str
    fname   : str — path to the processed .ecsv file
    meta    : dict from KeyHandler
    maps    : pd.DataFrame
    cubes   : pd.DataFrame
    params  : dict from SourceHandler
    """
    save_mom_maps   = meta.get("save_mom_maps",  True)
    save_maps   = meta.get("save_maps", True)
    folder          = meta.get("folder_savefits", "./saved_FITS_files/")
    target_res_as   = _resolve_target_res(params, meta)
    spacing_per_beam = meta.get("spacing_per_beam", 2.0)

    if not (save_mom_maps or save_maps):
        print(f'{"[INFO]":<10}', f'Output writing disabled for {source}; skipping.')
        return

    if float(spacing_per_beam) < 4:
        print(f'{"[WARNING]":<10}',
              f'spacing_per_beam < 4 ({spacing_per_beam}); expect image artefacts.')

    os.makedirs(folder, exist_ok=True)

    # Load overlay header and build footprint mask
    data_dir     = meta.get("data_dir", "data/")
    overlay_file = meta.get("overlay_file", "")
    from os import path as _path
    overlay_fname = (_path.join(data_dir, overlay_file) if source in overlay_file
                     else _path.join(data_dir, source + overlay_file))

    ov_cube, ov_hdr = fits.getdata(overlay_fname, header=True)
    ov_slice = ov_cube[ov_hdr["NAXIS3"] // 2, :, :].copy()
    ov_slice[np.isfinite(ov_slice)] = 1.0
    ov_hdr_2d = twod_head(ov_hdr)

    this_data = Table.read(fname)
    ra_deg    = this_data["ra_deg"]
    dec_deg   = this_data["dec_deg"]

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
        print(f'{"[INFO]":<10}', f'Moment maps written to {folder}')

    if save_maps:
        for map_name in maps["map_name"]:
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "MAP_",  "map",  source, this_data, map_name, folder, target_res_as)
            save_to_fits(ra_deg, dec_deg, ov_hdr_2d, ov_slice, "EMAP_", "emap", source, this_data, map_name, folder, target_res_as)
        print(f'{"[INFO]":<10}', f'Map FITS files written to {folder}')


def _resolve_target_res(params, meta):
    import numpy as np
    resolution = meta.get("resolution", "angular")
    target_res = float(meta.get("target_res", 27.0))
    if resolution == "physical":
        dist_mpc = params.get("dist_mpc", 1.0)
        return 3600.0 * 180.0 / np.pi * 1e-6 * target_res / dist_mpc
    return target_res
