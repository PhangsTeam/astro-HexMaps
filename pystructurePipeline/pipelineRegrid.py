"""
pipelineRegrid.py — convolve and sample all maps and cubes onto the hex grid.

This stage is the core data-ingestion step of the pipeline.  It:

1. Calls pipelinePrepare to (re-)generate the hex grid and overlay header.
2. Initialises the output Astropy Table with source metadata and
   deprojected galactocentric coordinates.
3. For each 2D map:   convolves to the target beam → samples at hex points.
4. For each cube:     convolves to the target beam → reprojects onto overlay
                      WCS → samples at hex points.
5. Optionally samples an external mask cube/image.
6. Writes the table to disk as a .ecsv file.

The .ecsv file written here is the primary output format of PyStructure.
Subsequent stages (pipelineProducts, pipelineFITS) read and enrich it.

Column naming convention
------------------------
MAP_<NAME>   : sampled 2D map intensity
EMAP_<NAME>  : uncertainty on the 2D map
SPEC_<NAME>  : sampled spectral cube  (n_pts × n_chan array column)
ESPEC_<NAME> : uncertainty cube
rgal_as      : deprojected galactocentric radius in arcseconds
rgal_kpc     : deprojected galactocentric radius in kpc
rgal_r25     : deprojected galactocentric radius in units of r25
theta_rad    : deprojected polar angle in radians
"""

import os
import copy
import warnings
import numpy as np
import pandas as pd
from os import path
from datetime import date
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column
from astropy import units as au
from reproject import reproject_interp

from pystructurePipeline.utilsFits import twod_head, conv_with_gauss, deproject
from pystructurePipeline.pipelineSampling import run_sampling

warnings.filterwarnings("ignore")

_LOG_PREFIX = "[pyStructure] [Regrid]  "


# ============================================================================
# Velocity-axis helpers
# ============================================================================

def _get_vaxis(hdr):
    """
    Reconstruct the velocity axis from a FITS header.

    Uses the standard FITS WCS keywords CRVAL3, CDELT3, CRPIX3, NAXIS3.
    Returns an array of length NAXIS3 in whatever units CUNIT3 specifies.
    """
    v = np.arange(hdr["NAXIS3"])
    return (v - (hdr["CRPIX3"] - 1)) * hdr["CDELT3"] + hdr["CRVAL3"]


def _ensure_ms(hdr, data=None):
    """
    Ensure the velocity axis is in m/s and monotonically increasing.

    Some FITS cubes store velocities in km/s (CDELT3 < 200) or have a
    decreasing velocity axis (CDELT3 < 0).  Both need to be normalised so
    that the reprojection and shuffle steps work correctly.

    Modifies the header in place.  If *data* is provided and the axis needs
    flipping, the data array is also flipped along axis 0 and returned.
    """
    # Convert km/s → m/s
    if abs(hdr["CDELT3"]) < 200:
        hdr["CDELT3"] *= 1000
        hdr["CRVAL3"] *= 1000
        hdr["CUNIT3"]  = "m/s"

    # Flip decreasing axis
    if data is not None and hdr["CDELT3"] < 0:
        vaxis_inv    = _get_vaxis(hdr)
        hdr["CDELT3"] = abs(hdr["CDELT3"])
        hdr["CRPIX3"] = 1
        hdr["CRVAL3"] = vaxis_inv[-1]
        data = np.flip(data, axis=0)

    return hdr, data


# ============================================================================
# Spectral smoothing
# ============================================================================

def _spectral_smooth(data, hdr_out, spec_smooth):
    """
    Optionally smooth a data cube along the spectral axis.

    Controlled by config_key.txt settings spec_smooth and spec_smooth_method.

    Parameters
    ----------
    data        : np.ndarray (n_chan × n_y × n_x)
    hdr_out     : FITS header with spectral WCS
    spec_smooth : list [mode, method]
        mode   — "default" (no smoothing) | float (target resolution in km/s)
        method — "binned" | "gauss" | "combined"

    Returns
    -------
    data, hdr_out  — modified in place if smoothing was applied

    Smoothing methods
    -----------------
    gauss
        Convolve each spectrum with a Gaussian kernel whose width is chosen
        to bring the native channel width up to the target resolution.
    binned
        Average consecutive channel groups so that the new channel width
        equals the target resolution (rounded to the nearest integer ratio).
        Fast and preserves the noise properties of the data.
    combined
        Apply binned first, then Gaussian to reach the exact target resolution.
        Use this when the integer-ratio approximation of binned is not accurate
        enough.
    """
    from astropy.convolution import Gaussian1DKernel, convolve

    mode, method = spec_smooth[0], spec_smooth[1]

    # No smoothing requested
    if mode == "default":
        return data, hdr_out

    # mode must be a number (target resolution in km/s)
    if not isinstance(mode, (int, float)):
        return data, hdr_out

    spec_res    = abs(hdr_out["CDELT3"]) / 1000.0   # current channel width in km/s
    fwhm_factor = np.sqrt(8 * np.log(2))
    dim_data    = np.shape(data)

    if spec_res >= mode:
        print(f"{_LOG_PREFIX} [INFO]  No spectral smoothing; already at target resolution.")
        return data, hdr_out

    print(f"{_LOG_PREFIX} [INFO]  Spectral smoothing to {round(mode, 3)} km/s ({method}).")

    if method == "gauss":
        # Convolve with a Gaussian whose width bridges native → target resolution
        pix    = ((mode**2 - spec_res**2)**0.5 / spec_res) / fwhm_factor
        kernel = Gaussian1DKernel(pix)
        for s in range(dim_data[1] * dim_data[2]):
            y, x = s % dim_data[1], s // dim_data[1]
            data[:, y, x] = convolve(data[:, y, x], kernel, nan_treatment="fill")

    elif method in ("binned", "combined"):
        vaxis   = _get_vaxis(hdr_out)
        n_ratio = int(mode / spec_res)
        if (mode / spec_res - n_ratio) > 0.9:
            n_ratio += 1
        new_len = len(vaxis) // n_ratio

        if n_ratio > 1:
            new_vaxis = np.array([np.nanmean(vaxis[n_ratio*j:n_ratio*(j+1)])
                                   for j in range(new_len)])
            data      = np.array([np.nanmean(data[n_ratio*j:n_ratio*(j+1), :, :], axis=0)
                                   for j in range(new_len)])
            hdr_out["NAXIS3"] = new_len
            hdr_out["CDELT3"] = new_vaxis[1] - new_vaxis[0]
            hdr_out["CRVAL3"] = new_vaxis[0] + (hdr_out["CRPIX3"] - 1) * hdr_out["CDELT3"]

        # For "combined": apply a small Gaussian to reach the exact target
        if method == "combined" and n_ratio * spec_res < mode:
            pix    = ((mode**2 - (n_ratio * spec_res)**2)**0.5 / spec_res) / fwhm_factor
            kernel = Gaussian1DKernel(pix)
            for s in range(dim_data[1] * dim_data[2]):
                y, x = s % dim_data[1], s // dim_data[1]
                data[:, y, x] = convolve(data[:, y, x], kernel, nan_treatment="fill")

    return data, hdr_out


# ============================================================================
# Core sampling function
# ============================================================================

def sample_at_res(in_data, ra_samp, dec_samp, in_hdr=None,
                  target_res_as=None, target_hdr=None,
                  line_name="", source="", save_fits=False,
                  path_save_fits="", perbeam=False,
                  spec_smooth=("default", "binned"), unc=False):
    """
    Convolve *in_data* to *target_res_as* arcsec and sample at the hex points.

    This is the workhorse function of the regrid stage.  It handles both
    2D maps (NAXIS=2) and 3D spectral cubes (NAXIS=3) through the same code
    path, branching only where cube-specific steps are needed.

    Steps
    -----
    1. Load the FITS data if a path is given, or use the array directly.
    2. Convolve to the target beam using conv_with_gauss (skipped if the input
       beam is already larger or if no BMAJ header keyword is present).
    3. For cubes: ensure velocity axis is in m/s and monotonically increasing.
    4. Apply optional spectral smoothing.
    5. Reproject onto the target WCS using nearest-neighbour interpolation.
    6. Sample the reprojected data at the hex-grid pixel positions.

    Parameters
    ----------
    in_data       : str or np.ndarray — input FITS path or data array
    ra_samp       : np.ndarray        — RA of sampling points (degrees)
    dec_samp      : np.ndarray        — Dec of sampling points (degrees)
    in_hdr        : FITS Header       — required if in_data is an array
    target_res_as : float             — target beam FWHM in arcseconds
    target_hdr    : FITS Header       — WCS to reproject onto (overlay header)
    line_name     : str               — label for optional FITS output
    source        : str               — source name for optional FITS output
    save_fits     : bool              — write the convolved intermediate FITS
    path_save_fits: str               — directory for intermediate FITS output
    perbeam       : bool              — correct for beam area change (use for
                                        maps in Jy/beam or K units)
    spec_smooth   : (mode, method)    — spectral smoothing parameters
    unc           : bool              — treat as uncertainty map (square before
                                        convolving, sqrt after)

    Returns
    -------
    result  : np.ndarray — sampled values; shape (n_pts,) for 2D maps or
                           (n_pts, n_chan) for cubes.
    trg_hdr : FITS Header — the header onto which the data was projected.
    """
    if len(ra_samp) != len(dec_samp):
        print(f"{_LOG_PREFIX} [ERROR]  RA and Dec arrays must have the same length.")
        return ra_samp * np.nan, None

    # Load data
    if isinstance(in_data, str):
        if not path.exists(in_data):
            print(f"{_LOG_PREFIX} [ERROR]  File not found: {in_data}")
            return ra_samp * np.nan, None
        data, hdr = fits.getdata(in_data, header=True)
    else:
        data = copy.deepcopy(in_data)
        hdr  = in_hdr

    if target_res_as is None:
        target_res_as = 0

    dim_data = np.shape(data)
    is_cube  = (len(dim_data) == 3)
    trg_hdr  = copy.deepcopy(target_hdr) if target_hdr is not None else None

    # For 2D data, reduce the target header to 2D
    if not is_cube and trg_hdr is not None:
        trg_hdr = twod_head(trg_hdr)

    # ------------------------------------------------------------------
    # Spatial convolution
    # ------------------------------------------------------------------
    if "BMAJ" not in hdr:
        print(f"{_LOG_PREFIX} [WARNING]  No BMAJ in header of {in_data if isinstance(in_data, str) else 'array'}; skipping convolution.")
        hdr_out = copy.copy(hdr)
    elif hdr["BMAJ"] < 0.99 * target_res_as / 3600.0:
        print(f"{_LOG_PREFIX} [INFO]  Convolving {line_name} to {round(target_res_as, 2)} arcsec.")
        data, hdr_out = conv_with_gauss(
            in_data     = data,
            in_hdr      = hdr,
            target_beam = target_res_as * np.array([1.0, 1.0, 0.0]),
            quiet       = True,
            perbeam     = perbeam,
            unc         = unc,
        )
    else:
        print(f"{_LOG_PREFIX} [INFO]  {line_name} already at target resolution; skipping convolution.")
        hdr_out = copy.copy(hdr)

    # ------------------------------------------------------------------
    # Velocity-axis normalisation for cubes
    # ------------------------------------------------------------------
    if is_cube:
        if trg_hdr is not None:
            trg_hdr, _ = _ensure_ms(trg_hdr)
        hdr_out, data = _ensure_ms(hdr_out, data)

        # Flip the target header if its axis is also decreasing
        if trg_hdr is not None and trg_hdr["CDELT3"] < 0:
            vaxis_inv     = _get_vaxis(trg_hdr)
            trg_hdr["CDELT3"] = abs(trg_hdr["CDELT3"])
            trg_hdr["CRPIX3"] = 1
            trg_hdr["CRVAL3"] = vaxis_inv[-1]

        if hdr_out["CDELT3"] < 0:
            vaxis_inv     = _get_vaxis(hdr_out)
            hdr_out["CDELT3"] = abs(hdr_out["CDELT3"])
            hdr_out["CRPIX3"] = 1
            hdr_out["CRVAL3"] = vaxis_inv[-1]
            data = np.flip(data, axis=0)

    # ------------------------------------------------------------------
    # Spectral smoothing
    # ------------------------------------------------------------------
    data, hdr_out = _spectral_smooth(data, hdr_out, list(spec_smooth))

    # ------------------------------------------------------------------
    # Reprojection onto the overlay WCS
    # ------------------------------------------------------------------
    if trg_hdr is not None:
        # Remove REST frequency keyword which can confuse reproject
        for h in (hdr_out, trg_hdr):
            h.remove("RESTF*", ignore_missing=True)

        # Adjust target spectral axis if spectral smoothing changed the channel width
        if isinstance(spec_smooth[0], (int, float)) and spec_smooth[0] != "default":
            if spec_smooth[0] > trg_hdr.get("CDELT3", 0) / 1000.0:
                vaxis_ov  = _get_vaxis(trg_hdr)
                new_vaxis = np.arange(vaxis_ov[0], vaxis_ov[-1], spec_smooth[0] * 1000)
                trg_hdr["NAXIS3"] = len(new_vaxis)
                trg_hdr["CDELT3"] = spec_smooth[0] * 1000
                trg_hdr["CRVAL3"] = new_vaxis[0] + (trg_hdr["CRPIX3"] - 1) * trg_hdr["CDELT3"]

        data, _ = reproject_interp((data, hdr_out), trg_hdr, order="nearest-neighbor")

        if save_fits:
            out_hdr        = copy.copy(trg_hdr)
            out_hdr["BMAJ"] = target_res_as / 3600.0
            out_hdr["BMIN"] = target_res_as / 3600.0
            out_hdr["LINE"] = line_name
            fits.writeto(
                path.join(path_save_fits, f"{source}_{line_name}_{target_res_as}as.fits"),
                data=data, header=out_hdr, overwrite=True,
            )
    else:
        print(f"{_LOG_PREFIX} [INFO]  No target header supplied; skipping reprojection.")
        trg_hdr = hdr_out

    # ------------------------------------------------------------------
    # Sampling at the hex-grid positions
    # ------------------------------------------------------------------
    wcs_t = WCS(trg_hdr)
    if is_cube:
        pixel_coords = wcs_t.all_world2pix(
            np.column_stack((ra_samp, dec_samp, np.zeros(len(dec_samp)))), 0)
    else:
        pixel_coords = wcs_t.all_world2pix(np.column_stack((ra_samp, dec_samp)), 0)

    samp_x   = np.array(np.rint(pixel_coords[:, 0]), dtype=int)
    samp_y   = np.array(np.rint(pixel_coords[:, 1]), dtype=int)
    n_pts    = len(samp_x)
    dim_data = np.shape(data)

    result = np.full((n_pts, dim_data[0]), np.nan) if is_cube else np.full(n_pts, np.nan)

    if is_cube:
        in_bounds = np.where(
            (samp_x > 0) & (samp_x < dim_data[2]) &
            (samp_y > 0) & (samp_y < dim_data[1])
        )[0]
        for kk in in_bounds:
            result[kk, :] = data[:, samp_y[kk], samp_x[kk]]
    else:
        in_bounds = np.where(
            (samp_x > 0) & (samp_x < dim_data[1]) &
            (samp_y > 0) & (samp_y < dim_data[0])
        )[0]
        result[in_bounds] = data[samp_y[in_bounds], samp_x[in_bounds]]

    return result, trg_hdr


def sample_mask(in_data, ra_samp, dec_samp, in_hdr=None, target_hdr=None):
    """
    Reproject and sample a binary mask cube or image onto the hex-grid points.

    Identical to sample_at_res but skips convolution (masks should not be
    smoothed) and uses nearest-neighbour interpolation to preserve binary values.

    Returns
    -------
    result  : np.ndarray — sampled mask values (0 or 1)
    trg_hdr : FITS Header
    """
    if isinstance(in_data, str):
        if not path.exists(in_data):
            return ra_samp * np.nan, None
        data, hdr = fits.getdata(in_data, header=True)
    else:
        data = copy.deepcopy(in_data)
        hdr  = in_hdr

    dim_data = np.shape(data)
    is_cube  = (len(dim_data) == 3)
    trg_hdr  = copy.deepcopy(target_hdr)

    if not is_cube and trg_hdr is not None:
        trg_hdr = twod_head(trg_hdr)

    if trg_hdr is not None:
        hdr_out = copy.copy(hdr)
        for h in (hdr_out, trg_hdr):
            h.remove("RESTF*", ignore_missing=True)
        data, _ = reproject_interp((data, hdr_out), trg_hdr, order="nearest-neighbor")

    wcs_t = WCS(trg_hdr)
    if is_cube:
        pixel_coords = wcs_t.all_world2pix(
            np.column_stack((ra_samp, dec_samp, np.zeros(len(dec_samp)))), 0)
    else:
        pixel_coords = wcs_t.all_world2pix(np.column_stack((ra_samp, dec_samp)), 0)

    samp_x   = np.array(np.rint(pixel_coords[:, 0]), dtype=int)
    samp_y   = np.array(np.rint(pixel_coords[:, 1]), dtype=int)
    n_pts    = len(samp_x)
    dim_data = np.shape(data)

    result = np.full((n_pts, dim_data[0]), np.nan) if is_cube else np.full(n_pts, np.nan)

    if is_cube:
        in_bounds = np.where(
            (samp_x > 0) & (samp_x < dim_data[2]) &
            (samp_y > 0) & (samp_y < dim_data[1])
        )[0]
        for kk in in_bounds:
            result[kk, :] = data[:, samp_y[kk], samp_x[kk]]
    else:
        in_bounds = np.where(
            (samp_x > 0) & (samp_x < dim_data[1]) &
            (samp_y > 0) & (samp_y < dim_data[0])
        )[0]
        result[in_bounds] = data[samp_y[in_bounds], samp_x[in_bounds]]

    return result, trg_hdr


# ============================================================================
# Stage entry point
# ============================================================================

def run_regrid(source, params, meta, maps, cubes, input_mask):
    """
    Convolve and sample all maps and cubes for *source*.

    This function drives the full regrid stage:
    - calls run_sampling to get the hex grid and overlay header
    - initialises the output Astropy Table
    - loops over 2D maps and spectral cubes
    - optionally samples the external mask
    - writes the table to a .ecsv file

    Parameters
    ----------
    source     : str
    params     : dict  — from SourceHandler.get_source_params()
    meta       : dict  — from KeyHandler.meta
    maps       : pd.DataFrame — 2D map definitions from handlerKeys
    cubes      : pd.DataFrame — spectral cube definitions from handlerKeys
    input_mask : pd.DataFrame — mask definition from handlerKeys

    Returns
    -------
    fname : str — path of the written .ecsv file
    """
    from pystructurePipeline import __version__, __author__, __email__, __credits__

    # Generate sampling grid
    sampling      = run_sampling(source=source, params=params, meta=meta)
    samp_ra       = sampling["samp_ra"]
    samp_dec      = sampling["samp_dec"]
    ov_hdr        = sampling["ov_hdr"]
    target_res_as = sampling["target_res_as"]
    n_chan         = ov_hdr["NAXIS3"]
    n_pts          = len(samp_ra)

    fname              = _build_fname(source, meta)
    structure_creation = meta.get("structure_creation", "default")
    data_dir           = meta.get("data_dir", "data/")
    save_fits          = meta.get("save_fits", False)

    # Decide whether to create a fresh table or fill an existing one
    if "fill" in structure_creation and path.exists(fname):
        print(f"{_LOG_PREFIX} [INFO]  Fill mode: loading existing table from {fname}.")
        this_data, fill_maps, fill_cubes = _fill_checker(fname, samp_ra, samp_dec, maps, cubes)
    else:
        this_data  = _init_table(source, params, meta, samp_ra, samp_dec, ov_hdr,
                                 target_res_as, __version__, __author__, __email__, __credits__)
        fill_maps, fill_cubes = [], []

    # ------------------------------------------------------------------
    # Process 2D maps
    # ------------------------------------------------------------------
    for _, map_entry in maps.iterrows():
        if map_entry["map_name"] in fill_maps:
            print(f"{_LOG_PREFIX} [INFO]  Map {map_entry['map_name']} already present; skipping.")
            continue

        map_file = path.join(str(map_entry["map_dir"]), source + str(map_entry["map_ext"]))
        if not path.exists(map_file):
            print(f"{_LOG_PREFIX} [ERROR]  Map {map_entry['map_name']} not found: {map_file}")
            continue

        perbeam  = "/beam" in str(map_entry.get("map_unit", ""))
        this_int, _ = sample_at_res(
            map_file, samp_ra, samp_dec,
            target_res_as  = target_res_as,
            target_hdr     = ov_hdr,
            line_name      = map_entry["map_name"],
            source         = source,
            path_save_fits = data_dir,
            save_fits      = save_fits,
            perbeam        = perbeam,
        )
        this_data["MAP_" + map_entry["map_name"].upper()] = Column(
            this_int, unit=au.Unit(str(map_entry["map_unit"])),
            description=map_entry["map_desc"])

        # Optional uncertainty map
        if str(map_entry.get("map_uc", "")).strip():
            uc_file = path.join(str(map_entry["map_dir"]), source + str(map_entry["map_uc"]))
            if path.exists(uc_file):
                uc_int, _ = sample_at_res(
                    uc_file, samp_ra, samp_dec,
                    target_res_as = target_res_as, target_hdr = ov_hdr,
                    perbeam=perbeam, unc=True,
                )
                this_data["EMAP_" + map_entry["map_name"].upper()] = Column(
                    uc_int, unit=au.Unit(str(map_entry["map_unit"])),
                    description=f'Uncertainty: {map_entry["map_desc"]}')

        print(f"{_LOG_PREFIX} [INFO]  Map {map_entry['map_name']} sampled successfully.")

    # ------------------------------------------------------------------
    # Process spectral cubes
    # ------------------------------------------------------------------
    for _, cube in cubes.iterrows():
        if cube["line_name"] in fill_cubes:
            print(f"{_LOG_PREFIX} [INFO]  Cube {cube['line_name']} already present; skipping.")
            continue

        cube_file = path.join(str(cube["line_dir"]), source + str(cube["line_ext"]))
        if not path.exists(cube_file):
            print(f"{_LOG_PREFIX} [ERROR]  Cube {cube['line_name']} not found: {cube_file}")
            continue

        this_spec, _ = sample_at_res(
            cube_file, samp_ra, samp_dec,
            target_res_as  = target_res_as,
            target_hdr     = ov_hdr,
            line_name      = cube["line_name"],
            source         = source,
            path_save_fits = data_dir,
            save_fits      = save_fits,
        )
        this_data["SPEC_" + cube["line_name"].upper()] = Column(
            this_spec, unit=au.Unit(str(cube["line_unit"])),
            description=cube["line_desc"])

        # Optional 2D integrated-intensity map provided alongside the cube
        map_ext = str(cube.get("map_ext", "")).strip()
        if map_ext and map_ext not in ("nan", ""):
            b2d_file = path.join(str(cube["line_dir"]), source + map_ext)
            if path.exists(b2d_file):
                b2d, _ = sample_at_res(b2d_file, samp_ra, samp_dec,
                                       target_res_as=target_res_as, target_hdr=ov_hdr)
                this_data["MAP_" + cube["line_name"].upper()] = Column(
                    b2d, unit=au.Unit(str(cube["line_unit"])),
                    description=cube["line_desc"])

        # Optional 2D uncertainty map for the cube
        map_uc = str(cube.get("map_uc", "")).strip()
        if map_uc and map_uc not in ("nan", ""):
            uc_file = path.join(str(cube["line_dir"]), source + map_uc)
            if path.exists(uc_file):
                uc, _ = sample_at_res(uc_file, samp_ra, samp_dec,
                                      target_res_as=target_res_as, target_hdr=ov_hdr, unc=True)
                this_data["EMAP_" + cube["line_name"].upper()] = Column(
                    uc, unit=au.Unit(str(cube["line_unit"])),
                    description=f'Uncertainty: {cube["line_desc"]}')

        print(f"{_LOG_PREFIX} [INFO]  Cube {cube['line_name']} sampled successfully.")

    # ------------------------------------------------------------------
    # Optional external mask
    # ------------------------------------------------------------------
    if len(input_mask) > 0:
        use_fixed = meta.get("use_fixed_vel_mask", False)

        if use_fixed:
            # Build a binary mask from a fixed velocity window
            mask_unit  = input_mask["mask_unit"].iloc[0]
            mask_start = float(input_mask["mask_start"].iloc[0]) * au.Unit(mask_unit)
            mask_end   = float(input_mask["mask_end"].iloc[0])   * au.Unit(mask_unit)
            unit_v     = ov_hdr.get("CUNIT3", "m/s")
            v0, dv, crpix = ov_hdr["CRVAL3"], ov_hdr["CDELT3"], ov_hdr["CRPIX3"]
            vaxis     = (v0 + (np.arange(n_chan) - (crpix - 1)) * dv) * au.Unit(unit_v)
            vaxis     = vaxis.to(au.Unit(mask_unit))
            spec_mask = np.zeros((n_pts, n_chan))
            spec_mask[:, (vaxis >= mask_start) & (vaxis <= mask_end)] = 1.0
            print(f"{_LOG_PREFIX} [INFO]  Fixed velocity mask applied "
                  f"({mask_start} to {mask_end}).")
        else:
            # Sample an external FITS mask file
            mask_file = path.join(
                str(input_mask["mask_dir"].iloc[0]),
                source + str(input_mask["mask_ext"].iloc[0]),
            )
            if not path.exists(mask_file):
                print(f"{_LOG_PREFIX} [ERROR]  Mask file not found: {mask_file}")
                spec_mask = np.zeros((n_pts, n_chan))
            else:
                spec_mask, _ = sample_mask(mask_file, samp_ra, samp_dec, target_hdr=ov_hdr)
                print(f"{_LOG_PREFIX} [INFO]  External mask sampled.")

        tag = "SPEC_" + str(input_mask["mask_name"].iloc[0]).upper()
        this_data[tag] = Column(
            spec_mask, unit=au.dimensionless_unscaled,
            description=str(input_mask["mask_desc"].iloc[0]),
        )

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    os.makedirs(meta.get("out_dir", "output/"), exist_ok=True)
    this_data.write(fname, format="ascii.ecsv", overwrite=True)
    print(f"{_LOG_PREFIX} [INFO]  PyStructure table written to: {fname}")
    return fname


# ============================================================================
# Private helpers
# ============================================================================

def _build_fname(source, meta):
    """
    Construct the output .ecsv filename.

    Encodes source name, resolution (value + unit suffix), and today's date.
    """
    resolution = meta.get("resolution", "angular")
    target_res = meta.get("target_res", 27.0)
    out_dir    = meta.get("out_dir", "output/")
    suffix = (str(int(target_res)) + "as" if resolution == "angular"
              else str(int(target_res)) + "pc" if resolution == "physical"
              else "native")
    date_str = date.today().strftime("%Y_%m_%d")
    return os.path.join(out_dir, f"{source}_data_struct_{suffix}_{date_str}.ecsv")


def _init_table(source, params, meta, samp_ra, samp_dec, ov_hdr,
                target_res_as, version, author, email, credits_):
    """
    Create and populate an empty Astropy Table for *source*.

    Writes provenance metadata (version, author, date) and coordinate
    columns including deprojected galactocentric radius and polar angle.

    Parameters
    ----------
    source        : str
    params        : dict — source geometry from SourceHandler
    meta          : dict — pipeline settings from KeyHandler
    samp_ra/dec   : arrays — hex-grid positions
    ov_hdr        : FITS Header — spectral axis information
    target_res_as : float — beam FWHM in arcsec (written to metadata)
    version/author/email/credits_ : package metadata strings

    Returns
    -------
    this_data : astropy.table.Table
    """
    from datetime import date as _date
    this_data = Table()

    # Provenance metadata stored in the table header
    this_data.meta.update({
        "Name":     "PyStructure",
        "Version":  version,
        "Authors":  author,
        "Contacts": email,
        "Credits":  credits_,
        "User":     meta.get("user", ""),
        "Comments": meta.get("comments", ""),
        "Date":     _date.today().strftime("%Y_%m_%d"),
        "Source":   source,
    })

    # Sky coordinates
    this_data["ra_deg"]  = Column(samp_ra,  unit=au.deg, description="Right ascension (J2000)")
    this_data["dec_deg"] = Column(samp_dec, unit=au.deg, description="Declination (J2000)")

    # Source geometry metadata
    this_data.meta["dist_mpc"]   = params["dist_mpc"]   * au.Mpc
    this_data.meta["posang_deg"] = params["posang_deg"] * au.deg
    this_data.meta["incl_deg"]   = params["incl_deg"]   * au.deg
    this_data.meta["beam_as"]    = target_res_as         * au.arcsec

    # Spectral axis metadata (from the overlay cube header)
    unit_v = ov_hdr.get("CUNIT3", "m/s")
    this_data.meta["SPEC_VCHAN0"] = ov_hdr["CRVAL3"] * au.Unit(unit_v)
    this_data.meta["SPEC_DELTAV"] = ov_hdr["CDELT3"] * au.Unit(unit_v)
    this_data.meta["SPEC_CRPIX"]  = ov_hdr["CRPIX3"]
    this_data.meta["input_maps"]  = ""
    this_data.meta["input_cubes"] = ""

    # Deprojected galactocentric coordinates
    rgal_deg, theta_rad = deproject(
        samp_ra, samp_dec,
        [params["posang_deg"], params["incl_deg"], params["ra_ctr"], params["dec_ctr"]],
        vector=True,
    )
    dist_mpc = params["dist_mpc"]
    r25      = params["r25"]

    this_data["rgal_as"]   = Column(rgal_deg * 3600,
                                    unit=au.arcsec,
                                    description="Deprojected galactocentric radius")
    this_data["rgal_kpc"]  = Column(np.deg2rad(rgal_deg) * dist_mpc * 1e3,
                                    unit=au.kpc,
                                    description="Deprojected galactocentric radius")
    this_data["rgal_r25"]  = Column(rgal_deg / (r25 / 60.0),
                                    description="Deprojected galactocentric radius (r25 units)")
    this_data["theta_rad"] = Column(theta_rad,
                                    unit=au.rad,
                                    description="Deprojected polar angle")
    return this_data


def _fill_checker(fname, samp_ra, samp_dec, maps, cubes):
    """
    Load an existing .ecsv file and identify which maps/cubes still need filling.

    Validates that the coordinates in the existing file match the current
    sampling grid before adding new columns.  Raises ValueError if they differ
    (which would indicate that the key files have changed in a way that altered
    the grid, requiring a full re-run with structure_creation = "default").

    Returns
    -------
    this_data   : Table — the loaded existing table
    fill_maps   : list  — map names that are already present and can be skipped
    fill_cubes  : list  — cube names that are already present and can be skipped
    """
    this_data = Table.read(fname)
    diff = (abs(np.nansum(this_data["ra_deg"]  - samp_ra  * au.deg))
            + abs(np.nansum(this_data["dec_deg"] - samp_dec * au.deg)))
    if diff > 1e-12 * au.deg:
        raise ValueError(
            f"{_LOG_PREFIX} [ERROR]  Existing file coordinates do not match the "
            "current sampling grid.  Set structure_creation = 'default' to overwrite."
        )
    fill_maps  = [b for b in maps["map_name"]    if f"MAP_{b.upper()}"  in this_data.colnames]
    fill_cubes = [c for c in cubes["line_name"]  if f"MOM0_{c.upper()}" in this_data.colnames]
    return this_data, fill_maps, fill_cubes
