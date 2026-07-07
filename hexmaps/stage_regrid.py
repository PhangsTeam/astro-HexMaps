"""
stage_regrid.py — generate the hex sampling grid, then convolve and sample
all maps and cubes onto it.

This stage is the core data-ingestion step of the pipeline.  It:

1. Generates the hexagonal sampling grid from the overlay cube (run_sampling).
   The grid is hexagonal (close-packed circles) rather than rectangular
   because it provides more uniform spatial coverage with the fewest beams
   and minimises the number of correlated positions for a given beam spacing.
   The grid is clipped to the footprint of the overlay cube (pixels with at
   least one finite channel), with spacing derived from the target resolution
   and the pixels_per_beam parameter in config.txt.
2. Initialises the output Astropy Table with target metadata and
   deprojected galactocentric coordinates.
3. For each 2D map:   convolves to the target beam → samples at hex points.
4. For each cube:     convolves to the target beam → reprojects onto overlay
                      WCS → samples at hex points.
5. Optionally samples an external mask cube/image.
6. Writes the table to disk as a .ecsv file.

The .ecsv file written here is the primary output format of HexMaps.
Subsequent stages (stage_products, stage_fits) read and enrich it.

Column naming convention
------------------------
MAP_<NAME>   : sampled 2D map intensity
EMAP_<NAME>  : uncertainty on the 2D map
SPEC_<NAME>  : sampled spectral cube  (n_pts × n_chan array column)
ESPEC_<NAME> : uncertainty cube
RGAL_AS      : deprojected galactocentric radius in arcseconds
RGAL_KPC     : deprojected galactocentric radius in kpc
RGAL_R25     : deprojected galactocentric radius in units of r25
THETA_RAD    : deprojected polar angle in radians
"""

import os
import copy
import warnings
import numpy as np
import pandas as pd
from os import path
from pathlib import Path
from datetime import date
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column
from astropy import units as au
from reproject import reproject_interp
from scipy.ndimage import binary_erosion
from skimage.morphology import disk

from hexmaps.utils_fits import (
    twod_head,
    conv_with_gauss,
    deproject,
    make_sampling_points,
    reproject_cube,
)

from hexmaps.logger import get_logger

LOG = get_logger("Regrid")

warnings.filterwarnings("ignore")


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

    Modifies the header in place.  If *data* is provided and the axis needs
    flipping, the data array is also flipped along axis 0 and returned.
    """

    # Flip decreasing axis
    if data is not None and hdr["CDELT3"] < 0:
        vaxis_inv = _get_vaxis(hdr)
        hdr["CDELT3"] = abs(hdr["CDELT3"])
        hdr["CRPIX3"] = 1
        hdr["CRVAL3"] = vaxis_inv[-1]
        data = np.flip(data, axis=0)

    return hdr, data


# def _fix_gildas_header(hdr):
#     """Fix common GILDAS/CLASS FITS header issues for astropy compatibility."""
#     hdr = hdr.copy()

#     # Fix 1: GLS -> SFL (affects both cubes)
#     for key in ['CTYPE1', 'CTYPE2']:
#         if 'GLS' in hdr.get(key, ''):
#             hdr[key] = hdr[key].replace('GLS', 'SFL')

#     # Fix 2: CRPIX3 = 0 (only target cube, but check both)
#     if hdr.get('CRPIX3', 1) == 0.0:
#         hdr['CRVAL3'] = hdr['CRVAL3'] + hdr['CDELT3']
#         hdr['CRPIX3'] = 1.0

#     return hdr

# def _harmonize_restfreq(hdr_in, hdr_target):
#     """
#     Ensure both headers carry a consistent RESTFRQ before reprojection.

#     astropy's WCS spectral machinery needs RESTFRQ to be present to construct
#     the spectral coordinate transform — even when both axes are simple linear
#     velocity scales in m/s and the actual CTYPE3 convention (e.g. "VRAD" vs
#     "VELO-LSR") would not otherwise matter. If RESTFRQ is missing from either
#     header, reproject_interp silently returns an all-NaN array with an empty
#     footprint for the spectral axis, with no warning or error raised.

#     This previously was handled by removing RESTFRQ from both headers
#     entirely, which avoided an astropy "latest version" issue with stale
#     RESTFRQ values but reintroduces the all-NaN failure whenever the two
#     cubes use different CTYPE3 conventions (a common situation when combining
#     cubes from different surveys/pipelines).

#     The fix: if either header has a RESTFRQ, copy that value to both headers
#     (preferring the input cube's own value if both have one). The small
#     radio/optical velocity-convention difference this introduces is of order
#     (v/c)^2 — many orders of magnitude below the channel width — and is
#     negligible since HexMaps's own velocity axis (CRVAL3/CDELT3/CRPIX3 in
#     m/s) is used for all scientific calculations, not the WCS spectral
#     transform. Only if NEITHER header has a RESTFRQ is it removed from both,
#     as before.

#     Parameters
#     ----------
#     hdr_in     : FITS Header — input cube header (modified in place)
#     hdr_target : FITS Header — target/overlay header (modified in place)
#     """
#     hdr_in = hdr_in.copy()
#     hdr_target = hdr_target.copy()

#     restfreq = hdr_in.get("RESTFRQ") or hdr_target.get("RESTFRQ")
#     if restfreq:
#         hdr_in["RESTFRQ"] = restfreq
#         hdr_target["RESTFRQ"] = restfreq
#     else:
#         for h in (hdr_in, hdr_target):
#                 h.remove("RESTF*", ignore_missing=True)

#     return hdr_in, hdr_target

# ============================================================================
# Spectral smoothing
# ============================================================================


def _spectral_smooth(data, hdr_out, spec_smooth):
    """
    Optionally smooth a data cube along the spectral axis.

    Controlled by config.txt settings spec_smooth and spec_smooth_method.

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

    spec_res = abs(hdr_out["CDELT3"]) / 1000.0  # current channel width in km/s
    fwhm_factor = np.sqrt(8 * np.log(2))
    dim_data = np.shape(data)

    if spec_res >= mode:
        LOG.info(f"No spectral smoothing; already at target resolution.")
        return data, hdr_out

    LOG.info(f"Spectral smoothing to {round(mode, 3)} km/s ({method}).")

    if method == "gauss":
        # Convolve with a Gaussian whose width bridges native → target resolution
        pix = ((mode**2 - spec_res**2) ** 0.5 / spec_res) / fwhm_factor
        kernel = Gaussian1DKernel(pix)
        for s in range(dim_data[1] * dim_data[2]):
            y, x = s % dim_data[1], s // dim_data[1]
            data[:, y, x] = convolve(data[:, y, x], kernel, nan_treatment="fill")

    elif method in ("binned", "combined"):
        vaxis = _get_vaxis(hdr_out)
        n_ratio = int(mode / spec_res)
        if (mode / spec_res - n_ratio) > 0.9:
            n_ratio += 1
        new_len = len(vaxis) // n_ratio

        if n_ratio > 1:
            new_vaxis = np.array(
                [
                    np.nanmean(vaxis[n_ratio * j : n_ratio * (j + 1)])
                    for j in range(new_len)
                ]
            )
            data = np.array(
                [
                    np.nanmean(data[n_ratio * j : n_ratio * (j + 1), :, :], axis=0)
                    for j in range(new_len)
                ]
            )
            hdr_out["NAXIS3"] = new_len
            hdr_out["CDELT3"] = new_vaxis[1] - new_vaxis[0]
            hdr_out["CRVAL3"] = (
                new_vaxis[0] + (hdr_out["CRPIX3"] - 1) * hdr_out["CDELT3"]
            )

        # For "combined": apply a small Gaussian to reach the exact target
        if method == "combined" and n_ratio * spec_res < mode:
            pix = (
                (mode**2 - (n_ratio * spec_res) ** 2) ** 0.5 / spec_res
            ) / fwhm_factor
            kernel = Gaussian1DKernel(pix)
            for s in range(dim_data[1] * dim_data[2]):
                y, x = s % dim_data[1], s // dim_data[1]
                data[:, y, x] = convolve(data[:, y, x], kernel, nan_treatment="fill")

    return data, hdr_out


# ============================================================================
# Hexagonal sampling grid
# ============================================================================


def run_sampling(target: str, params: dict, meta: dict) -> dict:
    """
    Generate the hexagonal sampling grid for *target*.

    Steps
    -----
    1. Load the overlay FITS cube and check it is 3-D.
    2. Determine the target resolution in arcseconds for this target:

       * ``angular``  — use ``meta["target_res"]`` directly (arcseconds).
       * ``physical`` — convert ``meta["target_res_config"]`` (parsecs) to
                        arcseconds using this target's ``dist_mpc``.
       * ``native``   — read BMAJ/BMIN from this target's overlay header.

    3. Collapse the cube along the spectral axis to create a binary footprint
       mask (True where at least one channel is finite).
    4. Build the hexagonal grid centred on the target and clip it to the mask.

    Parameters
    ----------
    target : str
        Target name; used to construct the overlay filename as
        ``{data_dir}/{target}{overlay_file}``.
    params : dict
        Target geometric parameters from TargetHandler.get_target_params().
        Required keys: ra_ctr, dec_ctr, dist_mpc.
    meta : dict
        Pipeline settings from KeyHandler.meta.
        Used keys: data_dir, overlay_file, resolution, target_res,
        pixels_per_beam, max_rad.

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
    data_dir = meta.get("data_dir", "data/")
    overlay_file = meta.get("overlay_file", "")

    # Construct the overlay filename: if the target name is already embedded in
    # overlay_file use it as-is, otherwise prepend the target name.
    overlay_fname = (
        path.join(data_dir, overlay_file)
        if target in overlay_file
        else path.join(data_dir, target + overlay_file)
    )

    if not path.exists(overlay_fname):
        LOG.error(f"Overlay file not found for {target}: {overlay_fname}")
        raise FileNotFoundError(f"Overlay file not found for {target}: {overlay_fname}")

    LOG.info(f"Overlay file: {overlay_fname}")
    ov_cube, ov_hdr = fits.getdata(overlay_fname, header=True)

    if ov_hdr["NAXIS"] == 4:
        LOG.error(f"4D overlay cube for {target}. " "Please provide a 3D cube.")
        raise ValueError(f"4D overlay cube for {target}. " "Please provide a 3D cube.")

    # ------------------------------------------------------------------
    # Build the footprint mask and generate the hex grid
    # ------------------------------------------------------------------
    # The mask is True wherever at least one spectral channel is finite.
    ov_footprint = np.sum(np.isfinite(ov_cube), axis=0) >= 1
    mask_hdr = twod_head(ov_hdr)

    # Resolve per-target target resolution (target_res, target_res_pc,
    # res_suffix) — single authoritative implementation shared with the
    # fits stage so both stages always compute the same values.
    from hexmaps.utils_fits import resolve_meta_resolution
    resolve_meta_resolution(target, params, meta, ov_hdr=ov_hdr, log=LOG)
    target_res_as = meta["target_res"]
    # Apply configurable FOV erosion to the hex grid footprint so that the
    # hex grid and all FITS products share the same effective footprint.
    # fov_erosion_beams is read from meta (set by handler_keys); default 0.5
    # gives the conventional half-beam trim. Set to 0 to disable erosion.
    fov_erosion_beams = meta.get("fov_erosion_beams", 0.5)
    pix_scale_as = abs(ov_hdr["CDELT1"]) * 3600.0
    trim_radius_pix = int(np.floor(fov_erosion_beams * target_res_as / pix_scale_as))
    if trim_radius_pix > 0:
        mask = binary_erosion(ov_footprint, structure=disk(trim_radius_pix))
        trim_as = fov_erosion_beams * target_res_as
        LOG.info(
            f"Hex grid footprint eroded by {trim_radius_pix} px "
            f"({fov_erosion_beams} beam = {trim_as:.1f} arcsec); "
            f"{mask.sum()} of {ov_footprint.sum()} pixels retained."
        )
    else:
        mask = ov_footprint
        if fov_erosion_beams > 0:
            LOG.warning(
                "Edge trim radius is <= 0 pixels; no hex grid edge removal applied."
            )
        else:
            LOG.info(
                "FOV erosion disabled (fov_erosion_beams = 0); "
                "full footprint used for hex grid."
            )

    pixels_per_beam = meta.get("pixels_per_beam", 2.0)
    max_rad = meta.get("max_rad", "auto")
    # Spacing in degrees: one beam FWHM divided by pixels_per_beam
    spacing = target_res_as / 3600.0 / float(pixels_per_beam)

    samp_ra, samp_dec = make_sampling_points(
        ra_ctr=params["ra_ctr"],
        dec_ctr=params["dec_ctr"],
        max_rad=max_rad,
        spacing=spacing,
        mask=mask,
        hdr_mask=mask_hdr,
        overlay_in=overlay_fname,
        show=False,
        log=LOG,
    )

    LOG.info(
        f"Hexagonal grid generated: "
        f"{len(samp_ra)} sampling points "
        f"(spacing = {spacing * 3600:.1f} arcsec)."
    )

    return dict(
        samp_ra=samp_ra,
        samp_dec=samp_dec,
        ov_hdr=ov_hdr,
        mask_hdr=mask_hdr,
        target_res_as=target_res_as,
    )


# ============================================================================
# Core sampling function
# ============================================================================


def sample_at_res(
    in_data,
    ra_samp,
    dec_samp,
    meta,
    in_hdr=None,
    target_res_as=None,
    target_hdr=None,
    line_name="",
    target="",
    perbeam=False,
    spec_smooth=("default", "binned"),
    unc=False,
):
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
    line_name     : str               — label used in log messages
    target        : str               — target name used in log messages
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
        LOG.error(f"RA and Dec arrays must have the same length.")
        return ra_samp * np.nan, None

    # Load data
    if isinstance(in_data, str):
        if not path.exists(in_data):
            LOG.error(f"File not found: {in_data}")
            return ra_samp * np.nan, None
        data, hdr = fits.getdata(in_data, header=True)
    else:
        data = copy.deepcopy(in_data)
        hdr = in_hdr

    if target_res_as is None:
        target_res_as = 0

    dim_data = np.shape(data)
    is_cube = len(dim_data) == 3
    trg_hdr = copy.deepcopy(target_hdr) if target_hdr is not None else None

    # For 2D data, reduce the target header to 2D
    if not is_cube and trg_hdr is not None:
        trg_hdr = twod_head(trg_hdr)

    # ------------------------------------------------------------------
    # Spatial convolution
    # ------------------------------------------------------------------
    if "BMAJ" not in hdr:
        LOG.warning(
            f"No BMAJ in header of {in_data if isinstance(in_data, str) else 'array'}; skipping convolution."
        )
        hdr_out = copy.copy(hdr)
    elif hdr["BMAJ"] < 0.99 * target_res_as / 3600.0:
        LOG.info(f"Convolving {line_name} to {round(target_res_as, 1)} arcsec.")
        data, hdr_out = conv_with_gauss(
            in_data=data,
            in_hdr=hdr,
            target_beam=target_res_as * np.array([1.0, 1.0, 0.0]),
            quiet=True,
            perbeam=perbeam,
            unc=unc,
            log=LOG,
        )
    else:
        LOG.info(f"{line_name} already at target resolution; skipping convolution.")
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
            vaxis_inv = _get_vaxis(trg_hdr)
            trg_hdr["CDELT3"] = abs(trg_hdr["CDELT3"])
            trg_hdr["CRPIX3"] = 1
            trg_hdr["CRVAL3"] = vaxis_inv[-1]

        if hdr_out["CDELT3"] < 0:
            vaxis_inv = _get_vaxis(hdr_out)
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

        # # Fix common GILDAS/CLASS FITS header issues for astropy compatibility
        # hdr_out = _fix_gildas_header(hdr_out)
        # trg_hdr = _fix_gildas_header(trg_hdr)

        # # Harmonize RESTFRQ between the two headers (see _harmonize_restfreq
        # # for why this is needed rather than simply removing it from both).
        # hdr_out, trg_hdr = _harmonize_restfreq(hdr_out, trg_hdr)

        # Adjust target spectral axis if spectral smoothing changed the channel width
        if isinstance(spec_smooth[0], (int, float)) and spec_smooth[0] != "default":
            if spec_smooth[0] > trg_hdr.get("CDELT3", 0) / 1000.0:
                vaxis_ov = _get_vaxis(trg_hdr)
                new_vaxis = np.arange(vaxis_ov[0], vaxis_ov[-1], spec_smooth[0] * 1000)
                trg_hdr["NAXIS3"] = len(new_vaxis)
                trg_hdr["CDELT3"] = spec_smooth[0] * 1000
                trg_hdr["CRVAL3"] = (
                    new_vaxis[0] + (trg_hdr["CRPIX3"] - 1) * trg_hdr["CDELT3"]
                )

        # Reproject the data onto the target header using bilinear interpolation
        if is_cube:
            data, _ = reproject_cube((data, hdr_out), trg_hdr, order="bilinear")
        else:
            data, _ = reproject_interp((data, hdr_out), trg_hdr, order="bilinear")

    else:
        LOG.info(f"No target header supplied; skipping reprojection.")
        trg_hdr = hdr_out

    # ------------------------------------------------------------------
    # Sampling at the hex-grid positions
    # ------------------------------------------------------------------
    wcs_t = WCS(trg_hdr)
    if is_cube:
        pixel_coords = wcs_t.all_world2pix(
            np.column_stack((ra_samp, dec_samp, np.zeros(len(dec_samp)))), 0
        )
    else:
        pixel_coords = wcs_t.all_world2pix(np.column_stack((ra_samp, dec_samp)), 0)

    samp_x = np.array(np.rint(pixel_coords[:, 0]), dtype=int)
    samp_y = np.array(np.rint(pixel_coords[:, 1]), dtype=int)
    n_pts = len(samp_x)
    dim_data = np.shape(data)

    result = (
        np.full((n_pts, dim_data[0]), np.nan) if is_cube else np.full(n_pts, np.nan)
    )

    if is_cube:
        in_bounds = np.where(
            (samp_x > 0)
            & (samp_x < dim_data[2])
            & (samp_y > 0)
            & (samp_y < dim_data[1])
        )[0]
        for kk in in_bounds:
            result[kk, :] = data[:, samp_y[kk], samp_x[kk]]
    else:
        in_bounds = np.where(
            (samp_x > 0)
            & (samp_x < dim_data[1])
            & (samp_y > 0)
            & (samp_y < dim_data[0])
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
        hdr = in_hdr

    dim_data = np.shape(data)
    is_cube = len(dim_data) == 3
    trg_hdr = copy.deepcopy(target_hdr)

    if not is_cube and trg_hdr is not None:
        trg_hdr = twod_head(trg_hdr)

    if trg_hdr is not None:
        hdr_out = copy.copy(hdr)
        # Harmonize RESTFRQ between the two headers (see _harmonize_restfreq
        # for why this is needed rather than simply removing it from both).
        # hdr_out, trg_hdr = _harmonize_restfreq(hdr_out, trg_hdr)
        if is_cube:
            data, _ = reproject_cube((data, hdr_out), trg_hdr, order="nearest-neighbor")
        else:
            data, _ = reproject_interp(
                (data, hdr_out), trg_hdr, order="nearest-neighbor"
            )

    wcs_t = WCS(trg_hdr)
    if is_cube:
        pixel_coords = wcs_t.all_world2pix(
            np.column_stack((ra_samp, dec_samp, np.zeros(len(dec_samp)))), 0
        )
    else:
        pixel_coords = wcs_t.all_world2pix(np.column_stack((ra_samp, dec_samp)), 0)

    samp_x = np.array(np.rint(pixel_coords[:, 0]), dtype=int)
    samp_y = np.array(np.rint(pixel_coords[:, 1]), dtype=int)
    n_pts = len(samp_x)
    dim_data = np.shape(data)

    result = (
        np.full((n_pts, dim_data[0]), np.nan) if is_cube else np.full(n_pts, np.nan)
    )

    if is_cube:
        in_bounds = np.where(
            (samp_x > 0)
            & (samp_x < dim_data[2])
            & (samp_y > 0)
            & (samp_y < dim_data[1])
        )[0]
        for kk in in_bounds:
            result[kk, :] = data[:, samp_y[kk], samp_x[kk]]
    else:
        in_bounds = np.where(
            (samp_x > 0)
            & (samp_x < dim_data[1])
            & (samp_y > 0)
            & (samp_y < dim_data[0])
        )[0]
        result[in_bounds] = data[samp_y[in_bounds], samp_x[in_bounds]]

    return result, trg_hdr


# ============================================================================
# Stage entry point
# ============================================================================


def run_regrid(target, params, meta, maps, cubes, input_mask, window_mask=None):
    """
    Convolve and sample all maps and cubes for *target*.

    This function drives the full regrid stage:
    - calls run_sampling (defined in this module) to get the hex grid and overlay header
    - initialises the output Astropy Table
    - loops over 2D maps and spectral cubes
    - optionally samples the external mask
    - writes the table to a .ecsv file

    Parameters
    ----------
    target     : str
    params     : dict  — from TargetHandler.get_target_params()
    meta       : dict  — from KeyHandler.meta
    maps       : pd.DataFrame — 2D map definitions from handler_keys
    cubes      : pd.DataFrame — spectral cube definitions from handler_keys
    input_mask : pd.DataFrame — mask definition from handler_keys

    Returns
    -------
    fname : str — path of the written .ecsv file
    """
    from hexmaps import __version__, __author__, __email__, __credits__

    # Generate sampling grid
    sampling = run_sampling(target=target, params=params, meta=meta)
    samp_ra = sampling["samp_ra"]
    samp_dec = sampling["samp_dec"]
    ov_hdr = sampling["ov_hdr"]
    target_res_as = sampling["target_res_as"]
    n_chan = ov_hdr["NAXIS3"]
    n_pts = len(samp_ra)

    fname = _build_fname(target, meta)
    structure_creation = meta.get("structure_creation", "default")
    data_dir = meta.get("data_dir", "data/")

    # Decide whether to create a fresh table or fill an existing one
    if "fill" in structure_creation and path.exists(fname):
        LOG.info(f"Fill mode: loading existing table from {fname}.")
        this_data, fill_maps, fill_cubes = _fill_checker(
            fname, samp_ra, samp_dec, maps, cubes
        )
    else:
        this_data = _init_table(
            target,
            params,
            meta,
            samp_ra,
            samp_dec,
            ov_hdr,
            target_res_as,
            __version__,
            __author__,
            __email__,
            __credits__,
        )
        fill_maps, fill_cubes = [], []

    # ------------------------------------------------------------------
    # Collect raw FITS headers for all input files before any processing.
    # Populated during the maps/cubes loops below and written into the
    # table metadata at the end.
    # ------------------------------------------------------------------
    overlay_file = meta.get("overlay_file", "")
    overlay_fname = (
        path.join(data_dir, overlay_file)
        if target in overlay_file
        else path.join(data_dir, target + overlay_file)
    )
    input_headers = {}

    # Overlay header — raw, before any processing
    if path.exists(overlay_fname):
        input_headers["OVERLAY"] = fits.getheader(overlay_fname)

    # ------------------------------------------------------------------
    # Process 2D maps
    # ------------------------------------------------------------------
    for _, map_entry in maps.iterrows():
        if map_entry["map_name"] in fill_maps:
            LOG.info(f"Map {map_entry['map_name']} already present; skipping.")
            continue

        map_file = path.join(
            str(map_entry["map_dir"]), target + str(map_entry["map_ext"])
        )
        if not path.exists(map_file):
            LOG.error(f"Map {map_entry['map_name']} not found: {map_file}")
            continue

        # Capture the raw header before any processing
        input_headers[map_entry["map_name"].upper()] = fits.getheader(map_file)

        perbeam = "/beam" in str(map_entry.get("map_unit", ""))
        this_int, _ = sample_at_res(
            map_file,
            samp_ra,
            samp_dec,
            meta,
            target_res_as=target_res_as,
            target_hdr=ov_hdr,
            line_name=map_entry["map_name"],
            target=target,
            perbeam=perbeam,
        )
        this_data["MAP_" + map_entry["map_name"].upper()] = Column(
            this_int,
            unit=au.Unit(str(map_entry["map_unit"])),
            description=map_entry["map_desc"],
        )

        # Optional uncertainty map
        if str(map_entry.get("map_uc", "")).strip():
            uc_file = path.join(
                str(map_entry["map_dir"]), target + str(map_entry["map_uc"])
            )
            if path.exists(uc_file):
                input_headers["EMAP_" + map_entry["map_name"].upper()] = fits.getheader(
                    uc_file
                )
                uc_int, _ = sample_at_res(
                    uc_file,
                    samp_ra,
                    samp_dec,
                    meta,
                    target_res_as=target_res_as,
                    target_hdr=ov_hdr,
                    line_name=f"{map_entry['map_name']}_err",
                    perbeam=perbeam,
                    unc=True,
                )
                this_data["EMAP_" + map_entry["map_name"].upper()] = Column(
                    uc_int,
                    unit=au.Unit(str(map_entry["map_unit"])),
                    description=f'Uncertainty: {map_entry["map_desc"]}',
                )

        LOG.info(f"Map {map_entry['map_name']} sampled successfully.")

    # ------------------------------------------------------------------
    # Process spectral cubes
    # ------------------------------------------------------------------
    for _, cube in cubes.iterrows():
        if cube["line_name"] in fill_cubes:
            LOG.info(f"Cube {cube['line_name']} already present; skipping.")
            continue

        cube_file = path.join(str(cube["line_dir"]), target + str(cube["line_ext"]))
        if not path.exists(cube_file):
            LOG.error(f"Cube {cube['line_name']} not found: {cube_file}")
            continue

        # Capture the raw header before any processing
        input_headers[cube["line_name"].upper()] = fits.getheader(cube_file)

        this_spec, _ = sample_at_res(
            cube_file,
            samp_ra,
            samp_dec,
            meta,
            target_res_as=target_res_as,
            target_hdr=ov_hdr,
            line_name=cube["line_name"],
            target=target,
        )
        this_data["SPEC_" + cube["line_name"].upper()] = Column(
            this_spec,
            unit=au.Unit(str(cube["line_unit"])),
            description=cube["line_desc"],
        )

        # Optional 2D integrated-intensity map provided alongside the cube
        map_ext = str(cube.get("map_ext", "")).strip()
        if map_ext and map_ext not in ("nan", ""):
            b2d_file = path.join(str(cube["line_dir"]), target + map_ext)
            if path.exists(b2d_file):
                input_headers["MAP_" + cube["line_name"].upper()] = fits.getheader(
                    b2d_file
                )
                b2d, _ = sample_at_res(
                    b2d_file,
                    samp_ra,
                    samp_dec,
                    meta,
                    target_res_as=target_res_as,
                    target_hdr=ov_hdr,
                )
                this_data["MAP_" + cube["line_name"].upper()] = Column(
                    b2d,
                    unit=au.Unit(str(cube["line_unit"])),
                    description=cube["line_desc"],
                )

        # Optional 2D uncertainty map for the cube
        map_uc = str(cube.get("map_uc", "")).strip()
        if map_uc and map_uc not in ("nan", ""):
            uc_file = path.join(str(cube["line_dir"]), target + map_uc)
            if path.exists(uc_file):
                input_headers["EMAP_" + cube["line_name"].upper()] = fits.getheader(
                    uc_file
                )
                uc, _ = sample_at_res(
                    uc_file,
                    samp_ra,
                    samp_dec,
                    meta,
                    target_res_as=target_res_as,
                    target_hdr=ov_hdr,
                    unc=True,
                )
                this_data["EMAP_" + cube["line_name"].upper()] = Column(
                    uc,
                    unit=au.Unit(str(cube["line_unit"])),
                    description=f'Uncertainty: {cube["line_desc"]}',
                )

        LOG.info(f"Cube {cube['line_name']} sampled successfully.")

    # ------------------------------------------------------------------
    # Optional external masks — sampled onto the hex grid
    # ------------------------------------------------------------------
    # input_mask: external FITS file mask
    if input_mask is not None and len(input_mask) > 0:
        mask_file = path.join(
            str(input_mask["mask_dir"].iloc[0]),
            target + str(input_mask["mask_ext"].iloc[0]),
        )
        if not path.exists(mask_file):
            LOG.error(f"Input mask file not found: {mask_file}")
        else:
            _mask_label = "SPEC_" + str(input_mask["mask_name"].iloc[0]).upper()
            input_headers[_mask_label] = fits.getheader(mask_file)
            spec_mask, _ = sample_mask(
                mask_file, samp_ra, samp_dec, target_hdr=ov_hdr
            )
            tag = "SPEC_" + str(input_mask["mask_name"].iloc[0]).upper()
            this_data[tag] = Column(
                spec_mask,
                unit=au.dimensionless_unscaled,
                description=str(input_mask["mask_desc"].iloc[0]),
            )
            LOG.info(f"External input mask sampled ({tag}).")

    # window_mask: fixed velocity-window mask built from the spectral axis
    if window_mask is not None and len(window_mask) > 0:
        mask_unit  = window_mask["mask_unit"].iloc[0]
        mask_start = float(window_mask["mask_start"].iloc[0]) * au.Unit(mask_unit)
        mask_end   = float(window_mask["mask_end"].iloc[0])   * au.Unit(mask_unit)
        unit_v = ov_hdr.get("CUNIT3", "m/s")
        v0, dv, crpix = ov_hdr["CRVAL3"], ov_hdr["CDELT3"], ov_hdr["CRPIX3"]
        vaxis = (v0 + (np.arange(n_chan) - (crpix - 1)) * dv) * au.Unit(unit_v)
        vaxis = vaxis.to(au.Unit(mask_unit))
        spec_mask = np.zeros((n_pts, n_chan))
        spec_mask[:, (vaxis >= mask_start) & (vaxis <= mask_end)] = 1.0
        tag = "SPEC_" + str(window_mask["mask_name"].iloc[0]).upper()
        this_data[tag] = Column(
            spec_mask,
            unit=au.dimensionless_unscaled,
            description=str(window_mask["mask_desc"].iloc[0]),
        )
        LOG.info(
            f"Fixed velocity-window mask sampled "
            f"({mask_start} to {mask_end}, {tag})."
        )

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    # Embed all collected raw FITS headers into the table metadata.
    # Each header is serialised as a compact 80-char-per-card string
    # (FITS standard format, no newlines) via header.tostring(), so it
    # survives the ECSV round-trip without any encoding tricks.
    # Key: f"input_header_{label}" where label matches the table column key
    # (e.g. "input_header_12CO21", "input_header_SPIRE250", "input_header_OVERLAY").
    # Recover with: fits.Header.fromstring(table.meta["input_header_{label}"])
    for label, hdr in input_headers.items():
        try:
            this_data.meta[f"input_header_{label}"] = hdr.tostring()
        except Exception as e:
            LOG.warning(f"Could not serialise header for {label}: {e}")

    # Store the complete pipeline log up to this point in the table metadata.
    # Newlines are encoded as "\\n" (two characters) for the same ECSV
    # round-trip reason as config_file. Decode with:
    #   text = table.meta["pipeline_log"].replace("\\n", "\n")
    from hexmaps.logger import logger as _logger
    try:
        this_data.meta["pipeline_log"] = _logger.as_text().replace("\n", "\\n")
    except Exception as e:
        LOG.warning(f"Could not embed pipeline log in metadata: {e}")
        this_data.meta["pipeline_log"] = ""

    os.makedirs(meta.get("out_dir", "output/"), exist_ok=True)
    this_data.write(fname, format="ascii.ecsv", overwrite=True)
    LOG.info(f"Database written to: {fname}")

    return fname


# ============================================================================
# Private helpers
# ============================================================================


def _build_fname(target, meta):
    """
    Construct the output .ecsv filename.

    Encodes target name, resolution suffix, and today's date.
    The resolution suffix is read from ``meta["res_suffix"]`` which is set
    by handler_keys._resolve_resolution and kept current by run_sampling.
    """
    out_dir = meta.get("out_dir", "output/")
    res_suffix = meta.get("res_suffix", "27p0as")
    date_str = date.today().strftime("%Y_%m_%d")
    return os.path.join(out_dir, f"{target}_hexmaps_{res_suffix}_{date_str}.ecsv")


def _get_coord_names(ov_hdr):
    """
    Derive human-readable coordinate column names and descriptions from the
    spatial WCS axes of the overlay header.

    Reads CTYPE1 / CTYPE2 (e.g. "RA---TAN", "GLAT-CAR", "ELON-ZEA") and
    strips the projection code to produce clean axis labels.  CUNIT1 /
    CUNIT2 are inspected so non-degree units can be detected (they are still
    stored as degrees in the output columns after conversion).

    Returns
    -------
    col1_name  : str — column name for axis 1 (e.g. "RA", "GLON")
    col2_name  : str — column name for axis 2 (e.g. "DEC", "GLAT")
    col1_desc  : str — human-readable description for axis 1
    col2_desc  : str — human-readable description for axis 2
    """
    # Extract the axis type, stripping the projection suffix (e.g. "---TAN")
    raw1 = str(ov_hdr.get("CTYPE1", "RA")).split("-")[0].strip().upper()
    raw2 = str(ov_hdr.get("CTYPE2", "DEC")).split("-")[0].strip().upper()

    # Map known FITS type codes to friendly names
    _name_map = {
        "RA":   ("ra",   "Right ascension"),
        "DEC":  ("dec",  "Declination"),
        "GLON": ("glon", "Galactic longitude"),
        "GLAT": ("glat", "Galactic latitude"),
        "ELON": ("elon", "Ecliptic longitude"),
        "ELAT": ("elat", "Ecliptic latitude"),
        "HGLN": ("hgln", "Heliographic longitude"),
        "HGLT": ("hglt", "Heliographic latitude"),
    }
    short1, desc1 = _name_map.get(raw1, (raw1.upper(), raw1))
    short2, desc2 = _name_map.get(raw2, (raw2.upper(), raw2))

    col1_name = f"{short1.upper()}"
    col2_name = f"{short2.upper()}"
    col1_desc = f"{desc1} (degrees)"
    col2_desc = f"{desc2} (degrees)"

    return col1_name, col2_name, col1_desc, col2_desc


def _init_table(
    target,
    params,
    meta,
    samp_ra,
    samp_dec,
    ov_hdr,
    target_res_as,
    version,
    author,
    email,
    credits_,
):
    """
    Create and populate an empty Astropy Table for *target*.

    Writes provenance metadata (version, author, date) and coordinate
    columns including deprojected galactocentric radius and polar angle.

    Parameters
    ----------
    target        : str
    params        : dict — target geometry from TargetHandler
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
    this_data.meta.update(
        {
            "Name": "HexMaps",
            "Version": version,
            "Authors": author,
            "Contacts": email,
            "Credits": credits_,
            "User": meta.get("user", ""),
            "Comments": meta.get("comments", ""),
            "Date": _date.today().strftime("%Y_%m_%d"),
            "Target": target,
        }
    )

    # Store the full config file content so the .ecsv is self-contained.
    # Newlines are encoded as the two-character sequence "\\n" so the value
    # survives the ECSV round-trip (astropy collapses literal newlines in
    # metadata string values when reading back). Decode with:
    #   content = table.meta["config_file"].replace("\\n", "\n")
    conf_path = meta.get("conf_path", "")
    config_file_content = ""
    if conf_path and os.path.exists(conf_path):
        try:
            raw = Path(conf_path).read_text(encoding="utf-8")
            # Encode newlines so the string is single-line in the ECSV header
            config_file_content = raw.replace("\n", "\\n")
        except Exception as e:
            LOG.warning(f"Could not read config file for metadata: {e}")
    this_data.meta["config_file"] = config_file_content

    # ------------------------------------------------------------------
    # Sky coordinates — derive column names from the overlay WCS so the
    # table is correct for any coordinate system (RA/Dec, galactic, etc.)
    # ------------------------------------------------------------------
    col1_name, col2_name, col1_desc, col2_desc = _get_coord_names(ov_hdr)
    this_data[col1_name] = Column(
        samp_ra, unit=au.deg, description=col1_desc
    )
    this_data[col2_name] = Column(
        samp_dec, unit=au.deg, description=col2_desc
    )

    # Target geometry metadata — only the always-available values
    this_data.meta["dist_mpc"]  = params.get("dist_mpc", float("nan"))  * au.Mpc
    this_data.meta["beam_as"]   = target_res_as * au.arcsec

    # Optional galaxy geometry (PA, inclination) — stored only when valid
    incl_deg   = float(params.get("incl_deg",   float("nan")))
    posang_deg = float(params.get("posang_deg", float("nan")))
    r25        = float(params.get("r25",        float("nan")))
    import math
    has_geom = not any(math.isnan(v) for v in (incl_deg, posang_deg, r25))

    if has_geom:
        this_data.meta["posang_deg"] = posang_deg * au.deg
        this_data.meta["incl_deg"]   = incl_deg   * au.deg

    # Spectral axis metadata (from the overlay cube header)
    unit_v = ov_hdr.get("CUNIT3", "m/s")
    this_data.meta["SPEC_VCHAN0"] = ov_hdr["CRVAL3"] * au.Unit(unit_v)
    this_data.meta["SPEC_DELTAV"] = ov_hdr["CDELT3"] * au.Unit(unit_v)
    this_data.meta["SPEC_CRPIX"]  = ov_hdr["CRPIX3"]

    # ------------------------------------------------------------------
    # Deprojected galactocentric coordinates — galaxy targets only.
    # Requires incl_deg, posang_deg, and r25 to all be non-NaN in
    # target_definitions.txt.  Skipped with a warning for Milky Way clouds
    # or any other target type where these parameters are not meaningful.
    # ------------------------------------------------------------------
    if has_geom:
        ra_ctr  = float(params.get("ra_ctr",  0.0))
        dec_ctr = float(params.get("dec_ctr", 0.0))
        dist_mpc = float(params.get("dist_mpc", float("nan")))

        rgal_deg, theta_rad = deproject(
            samp_ra,
            samp_dec,
            [posang_deg, incl_deg, ra_ctr, dec_ctr],
            vector=True,
        )

        this_data["RGAL_AS"] = Column(
            rgal_deg * 3600,
            unit=au.arcsec,
            description="Deprojected galactocentric radius",
        )
        if not math.isnan(dist_mpc):
            this_data["RGAL_KPC"] = Column(
                np.deg2rad(rgal_deg) * dist_mpc * 1e3,
                unit=au.kpc,
                description="Deprojected galactocentric radius",
            )
        this_data["RGAL_R25"] = Column(
            rgal_deg / (r25 / 60.0),
            description="Deprojected galactocentric radius (r25 units)",
        )
        this_data["THETA_RAD"] = Column(
            theta_rad,
            unit=au.rad,
            description="Deprojected polar angle",
        )
    else:
        LOG.warning(
            f"Galaxy geometry columns (incl_deg, posang_deg, r25) are missing "
            f"or NaN for '{target}' in target_definitions.txt."
        )
        LOG.warning(
            f"Deprojected radius columns (RGAL_AS, RGAL_KPC, RGAL_R25, THETA_RAD) "
            f"will not be added to the database."
        )
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

    # Determine the coordinate column names from the existing table — they
    # depend on the WCS of the overlay cube and may be RA/DEC,
    # GLON/GLAT, etc.
    coord_cols = [c for c in this_data.colnames
                  if c.endswith("_deg") and c not in ("incl_deg", "posang_deg")]
    if len(coord_cols) >= 2:
        col1_name, col2_name = coord_cols[0], coord_cols[1]
    else:
        col1_name, col2_name = "RA", "DEC"

    diff = abs(np.nansum(this_data[col1_name] - samp_ra * au.deg)) + abs(
        np.nansum(this_data[col2_name] - samp_dec * au.deg)
    )
    if diff > 1e-12 * au.deg:
        LOG.error(
            f"Existing file coordinates do not match the "
            "current sampling grid.  Set structure_creation = 'default' to overwrite."
        )
        raise ValueError(
            f"Existing file coordinates do not match the "
            "current sampling grid.  Set structure_creation = 'default' to overwrite."
        )
    fill_maps = [
        b for b in maps["map_name"] if f"MAP_{b.upper()}" in this_data.colnames
    ]
    fill_cubes = [
        c for c in cubes["line_name"] if f"MOM0_{c.upper()}" in this_data.colnames
    ]
    return this_data, fill_maps, fill_cubes
