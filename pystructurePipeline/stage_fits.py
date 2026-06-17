"""
stage_fits.py — write FITS moment maps, 2D map images, and mask cube(s) for a source.

Moment maps: PPV-native computation (no hex grid)
---------------------------------------------------
Moment maps (MOM0/1/2, TPEAK, RMS, EW and their uncertainties) are computed
directly on the rectangular pixel-position-velocity (PPV) grid — never via
the hex-grid .ecsv table. This avoids the information loss and gridding
artefacts that come from sampling onto the irregular hex grid and then
regridding back onto a rectangular grid for FITS output.

For each cube, the pipeline:
  1. Obtains a convolved, overlay-WCS-aligned PPV cube for the line, either
     by reading it from disk (the file stage_regrid writes when save_fits is
     True) or, if that file is absent, by convolving and reprojecting the
     raw input cube itself (convolve_cube_to_target / reproject_cube_to_overlay).
  2. Builds (or reads) a PPV mask using exactly the same two-level S/N
     threshold + dilation algorithm as the hex-grid path (construct_mask_ppv
     mirrors stage_products.construct_mask channel-by-channel on the cube),
     including the same ref_line combination logic, the strict spatial
     filter (here implemented as 2-D connected-component labelling per
     channel, the rectangular-grid equivalent of the hex-grid distance-based
     filter), HFS mask extension, and the use_input_mask / use_fixed_vel_mask
     external-mask options.
  3. Computes moments on the masked PPV cube by reshaping it to
     (n_pix, n_chan) and calling utils_table.get_mom_maps — the exact same
     function used by the hex-grid path — then reshaping the results back to
     (ny, nx) maps.
  4. Writes one FITS file per moment quantity per line.

2D maps: unchanged; mask cube(s): now PPV-native too
---------------------------------------------------------
2D band/map columns (MAP_*/EMAP_*) have no PPV-cube equivalent (they are
already 2-D quantities in the .ecsv), so they are still regridded from the
hex-grid table via save_to_fits, exactly as before.

The velocity-integration mask(s), however, are now written PPV-native as
well (when save_mask is True): the same mask array built and used inside
run_moments_ppv (construct_mask_ppv / external_mask_ppv / etc.) is written
directly to FITS via save_ppv_mask_to_fits, with no hex-grid table involved
at any point. This requires save_mom_maps to also be True, since the mask
is only constructed while computing moments.

Output filename convention
--------------------------
{source}_{line}_{quantity}.fits        (moment maps, PPV-native)
{source}_{map}_{quantity}.fits         (2D maps, hex-grid regridded)
{source}_mask.fits / {source}_mask_<line>.fits   (mask cubes)

e.g.  ngc5194_12CO21_mom0.fits
      ngc5194_SPIRE250_map.fits
"""

import os
import copy
import numpy as np
import astropy.units as au
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from astropy.stats import median_absolute_deviation
from scipy.interpolate import griddata
from scipy.ndimage import label
from reproject import reproject_interp
from typing import Sequence, Union

from pystructurePipeline.utils_fits import twod_head, conv_with_gauss
from pystructurePipeline.stage_regrid import _harmonize_restfreq, _ensure_ms, _get_vaxis
from pystructurePipeline.utils_table import get_mom_maps

from pystructurePipeline.pystructureLogger import get_logger

LOG = get_logger("FITS")



# ============================================================================
# Grid helpers
# ============================================================================

def sample_to_hdr(in_data: Union[np.ndarray, Sequence[float]],
                  ra_samp: Union[np.ndarray, Sequence[float]],
                  dec_samp: Union[np.ndarray, Sequence[float]],
                  in_hdr: fits.Header) -> np.ndarray:
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


# ============================================================================
# PPV-native pipeline: convolution, reprojection, masking, moments
#
# Everything in this section operates on plain 2-D/3-D numpy arrays and
# FITS headers — never on the hex-grid .ecsv table. This mirrors the logic
# in stage_regrid.py (convolution, reprojection) and stage_products.py
# (mask construction) exactly, just applied directly to the rectangular
# PPV grid instead of sampling at hex-grid points first.
# ============================================================================

def get_convolved_ppv_cube(source, line_name, line_dir, line_ext, target_res_as,
                           ov_hdr, fits_dir, log=None):
    """
    Obtain a convolved PPV cube for *line_name*, aligned to the overlay WCS.

    Tries, in order:
      1. Read the cube stage_regrid already wrote when save_fits was True
         ({fits_dir}/{source}_{line_name}_{target_res_as}as.fits). This is
         already convolved to target_res_as and reprojected onto the overlay
         WCS, so it can be used directly.
      2. If that file is absent, fall back to reading the raw input cube
         from {line_dir}/{source}{line_ext}, convolving it to target_res_as
         (convolve_cube_to_target) and reprojecting it onto the overlay WCS
         (reproject_cube_to_overlay) from scratch.

    Parameters
    ----------
    source        : str
    line_name     : str   — cube name, e.g. "12co21"
    line_dir      : str   — directory containing the raw input cube
    line_ext      : str   — filename extension of the raw input cube
    target_res_as : float — target beam FWHM in arcseconds
    ov_hdr        : FITS Header — overlay header (3-D) defining the target WCS
    fits_dir      : str   — directory to look for the stage_regrid save_fits output
    log           : StageLogger, optional

    Returns
    -------
    data : np.ndarray (n_chan, ny, nx) — convolved PPV cube on the overlay grid
    hdr  : FITS Header — header matching *data*
    """
    log = log or LOG

    cached_path = os.path.join(fits_dir, f"{source}_{line_name}_{target_res_as}as.fits")
    if os.path.exists(cached_path):
        log.info(f"Using existing convolved cube for line: {line_name}: {cached_path}")
        return fits.getdata(cached_path, header=True)

    raw_path = os.path.join(line_dir, source + line_ext)
    log.info(f"No saved convolved cube found for line: {line_name}")
    log.info(f"Expected: {cached_path}")
    log.info(f"Convolving raw input from scratch: {raw_path}")
    if not os.path.exists(raw_path):
        log.error(f"Raw input cube not found for line: {line_name}: {raw_path}")
        raise FileNotFoundError(f"Raw input cube not found for line: {line_name}: {raw_path}")

    data, hdr = fits.getdata(raw_path, header=True)
    data, hdr = convolve_cube_to_target(data, hdr, target_res_as, log=log)
    data, hdr = reproject_cube_to_overlay(data, hdr, ov_hdr, log=log)
    return data, hdr


def convolve_cube_to_target(data, hdr, target_res_as, log=None):
    """
    Convolve a PPV cube to *target_res_as* arcsec, in place on its native grid.

    Thin wrapper around utils_fits.conv_with_gauss with the same calling
    convention stage_regrid.sample_at_res uses, so a from-scratch convolution
    here behaves identically to the one stage_regrid would have produced.

    Parameters
    ----------
    data          : np.ndarray (n_chan, ny, nx)
    hdr           : FITS Header — native header of *data*
    target_res_as : float — target beam FWHM in arcseconds
    log           : StageLogger, optional

    Returns
    -------
    data, hdr : convolved cube and updated header (BMAJ/BMIN updated)
    """
    log = log or LOG
    if "BMAJ" not in hdr:
        log.warning("No BMAJ in header; skipping convolution.")
        return data, hdr
    if hdr["BMAJ"] >= 0.99 * target_res_as / 3600.0:
        log.info("Cube already at or above target resolution; skipping convolution.")
        return data, hdr

    data, hdr_out = conv_with_gauss(
        in_data     = data,
        in_hdr      = hdr,
        target_beam = target_res_as * np.array([1.0, 1.0, 0.0]),
        quiet       = True,
        log         = log,
    )
    return data, hdr_out


def reproject_cube_to_overlay(data, hdr, ov_hdr, log=None):
    """
    Reproject a PPV cube onto the overlay's spatial+spectral WCS.

    Ensures the velocity axis is in m/s and monotonically increasing
    (matching the overlay convention), harmonizes RESTFRQ between the two
    headers (see stage_regrid._harmonize_restfreq), then reprojects with
    nearest-neighbour interpolation — the same approach stage_regrid.sample_at_res
    uses for cubes.

    Parameters
    ----------
    data   : np.ndarray (n_chan, ny, nx)
    hdr    : FITS Header — native header of *data*
    ov_hdr : FITS Header — overlay header (3-D) defining the target WCS
    log    : StageLogger, optional

    Returns
    -------
    data, hdr : reprojected cube and the (copied) overlay-aligned header
    """
    log = log or LOG
    trg_hdr = copy.deepcopy(ov_hdr)
    trg_hdr, _ = _ensure_ms(trg_hdr)
    hdr, data  = _ensure_ms(copy.copy(hdr), data)

    _harmonize_restfreq(hdr, trg_hdr)
    data, _ = reproject_interp((data, hdr), trg_hdr, order="nearest-neighbor")
    return data, trg_hdr


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


def construct_mask_ppv(ref_cube, SN_processing):
    """
    Build a two-level S/N velocity-integration mask from a PPV cube.

    Pixel-for-pixel identical algorithm to stage_products.construct_mask,
    just applied to a (n_chan, ny, nx) array instead of a (n_pts, n_chan)
    hex-grid table column. Each spatial pixel (y, x) is treated exactly like
    one hex-grid row: per-pixel noise via two-pass MAD, a high-S/N core mask
    requiring 3-of-3 consecutive channels, dilation into the low-S/N mask
    (5 passes), then a 2-channel edge grow.

    Parameters
    ----------
    ref_cube      : np.ndarray (n_chan, ny, nx) — reference line PPV cube
    SN_processing : list[float] — [low_SN_thresh, high_SN_thresh]

    Returns
    -------
    mask : np.ndarray (n_chan, ny, nx) — 0/1 integration mask
    """
    n_chan, ny, nx = ref_cube.shape

    # Two-pass global MAD to estimate the noise floor (identical to the
    # hex-grid version, just over the whole cube instead of the whole table)
    rms = median_absolute_deviation(ref_cube, axis=None, ignore_nan=True)
    rms = median_absolute_deviation(ref_cube[np.where(ref_cube < 3 * rms)], ignore_nan=True)

    # Per-pixel noise: MAD of channels below the global 3-sigma threshold.
    # axis=0 is the spectral axis here (vs axis=1 for the hex-grid table).
    mask_rough  = ref_cube < 3 * rms
    masked_cube = np.where(mask_rough, ref_cube, np.nan)
    med_mask    = np.nanmedian(masked_cube, axis=0)
    mad_mask    = np.nanmedian(np.abs(masked_cube - med_mask[None, :, :]), axis=0)

    low_thresh  = SN_processing[0] * mad_mask[None, :, :]
    high_thresh = SN_processing[1] * mad_mask[None, :, :]

    # Initial high-S/N mask: channel above high_thresh with adjacent support.
    # np.roll along axis=0 (spectral) replaces axis=1 from the hex-grid version.
    mask     = (ref_cube > high_thresh).astype(int)
    low_mask = (ref_cube > low_thresh).astype(int)
    mask = mask & (np.roll(mask, 1, 0) | np.roll(mask, -1, 0))

    # Require >= 3 of 3 consecutive channels to suppress single-channel spikes
    mask     = ((mask     + np.roll(mask,     1, 0) + np.roll(mask,     -1, 0)) >= 3).astype(int)
    low_mask = ((low_mask + np.roll(low_mask, 1, 0) + np.roll(low_mask, -1, 0)) >= 3).astype(int)

    # Dilate high-S/N core into low-S/N wings (5 passes)
    for _ in range(5):
        mask = (((mask + np.roll(mask, 1, 0) + np.roll(mask, -1, 0)) >= 1).astype(int)
                * low_mask)

    # Grow mask edge by 2 channels to ensure full line coverage
    for _ in range(2):
        mask = ((mask + np.roll(mask, 1, 0) + np.roll(mask, -1, 0)) >= 1).astype(int)

    return mask


def apply_strict_mask_ppv(mask, min_pixels=5):
    """
    Remove spatially isolated mask features, the PPV-grid equivalent of
    stage_products._apply_strict_mask.

    The hex-grid version labels spatially connected groups using a pairwise
    distance comparison between irregularly-spaced points — appropriate for
    a sparse hex grid, but both incorrect (the neighbour distance assumption
    doesn't hold) and prohibitively slow (O(n_pix^2) per channel) on a dense
    rectangular grid. The natural rectangular-grid equivalent is connected-
    component labelling on the regular pixel grid, which scipy.ndimage.label
    computes directly using 4-connectivity (matching the hex-grid filter's
    intent of "spatially adjacent" support).

    Mask features (connected components) smaller than *min_pixels* are
    removed, channel by channel.

    Parameters
    ----------
    mask       : np.ndarray (n_chan, ny, nx) — 0/1 mask array
    min_pixels : int — minimum connected-component size to keep (default 5,
                matching the hex-grid filter's hardcoded threshold)

    Returns
    -------
    mask : np.ndarray — filtered mask (same shape)
    """
    mask = mask.copy()
    n_chan = mask.shape[0]
    for ch in range(n_chan):
        labels, n_labels = label(mask[ch])
        if n_labels == 0:
            continue
        sizes = np.bincount(labels.ravel())
        small_labels = np.where(sizes < min_pixels)[0]
        small_labels = small_labels[small_labels != 0]   # never touch background
        if len(small_labels):
            mask[ch][np.isin(labels, small_labels)] = 0
    return mask


def build_hfs_mask_ppv(mask, line_name, hfs_data, delta_v_kms):
    """
    Extend a PPV mask to cover hyperfine satellite lines.

    Identical logic to stage_products._build_hfs_mask, applied along the
    spectral axis (axis=0) of a PPV cube instead of axis=1 of a hex-grid
    table column.

    Parameters
    ----------
    mask        : np.ndarray (n_chan, ny, nx) — existing 0/1 mask
    line_name   : str — name of the line to look up in hfs_data
    hfs_data    : pd.DataFrame — hyperfine structure table from handler_keys
    delta_v_kms : float — channel width in km/s

    Returns
    -------
    mask_hfs : np.ndarray (n_chan, ny, nx) — extended mask, or None if
               line_name is not in the HFS table.
    """
    lines_hfs = list(set(hfs_data["hfs_name"]))
    if line_name not in lines_hfs:
        return None

    idx_cols  = hfs_data["hfs_name"] == line_name
    restfreqs = [f * au.Unit(str(u)) for f, u in
                 zip(hfs_data["hfs_ref_freq"][idx_cols], hfs_data["unit"][idx_cols])]
    hfs_freqs = [f * au.Unit(str(u)) for f, u in
                 zip(hfs_data["hfs_freq"][idx_cols],     hfs_data["unit"][idx_cols])]

    mask_hfs = mask.copy()
    for freq, restfreq in zip(hfs_freqs, restfreqs):
        v_shift  = freq.to(au.km / au.s, equivalencies=au.doppler_radio(restfreq))
        shift_ch = int(np.rint(v_shift.value / delta_v_kms))

        mask_shift = np.zeros_like(mask, dtype=float)
        if shift_ch > 0:
            mask_shift[shift_ch:] = mask[:-shift_ch]
        elif shift_ch < 0:
            mask_shift[:shift_ch] = mask[-shift_ch:]
        else:
            mask_shift = mask.copy()

        mask_hfs[mask_shift == 1] = 1

    return mask_hfs


def fixed_velocity_mask_ppv(shape, ov_hdr, mask_start, mask_end, mask_unit):
    """
    Build a binary PPV mask from a fixed velocity window, the array-native
    equivalent of stage_regrid's use_fixed_vel_mask handling.

    Parameters
    ----------
    shape      : tuple (n_chan, ny, nx) — output mask shape
    ov_hdr     : FITS Header — overlay header (3-D), provides the velocity axis
    mask_start : astropy Quantity — start of the velocity window
    mask_end   : astropy Quantity — end of the velocity window
    mask_unit  : str — unit string for the window bounds

    Returns
    -------
    mask : np.ndarray (n_chan, ny, nx) — 0/1 mask, constant across all pixels
    """
    n_chan = shape[0]
    unit_v = ov_hdr.get("CUNIT3", "m/s")
    v0, dv, crpix = ov_hdr["CRVAL3"], ov_hdr["CDELT3"], ov_hdr["CRPIX3"]
    vaxis = (v0 + (np.arange(n_chan) - (crpix - 1)) * dv) * au.Unit(unit_v)
    vaxis = vaxis.to(au.Unit(mask_unit))

    chan_mask = (vaxis >= mask_start) & (vaxis <= mask_end)
    mask = np.zeros(shape)
    mask[chan_mask, :, :] = 1.0
    return mask


def external_mask_ppv(mask_file, ov_hdr, log=None):
    """
    Reproject an external FITS mask file onto the overlay grid, the
    array-native equivalent of stage_regrid.sample_mask's external-mask path
    (minus the final hex-grid sampling step, which doesn't apply here).

    Parameters
    ----------
    mask_file : str — path to the external mask FITS file (2-D or 3-D)
    ov_hdr    : FITS Header — overlay header defining the target WCS
    log       : StageLogger, optional

    Returns
    -------
    mask : np.ndarray (n_chan, ny, nx) — reprojected mask, broadcast across
           the spectral axis if the input mask was 2-D
    """
    log = log or LOG
    data, hdr = fits.getdata(mask_file, header=True)
    is_cube   = (data.ndim == 3)

    trg_hdr = copy.deepcopy(ov_hdr)
    if not is_cube:
        trg_hdr = twod_head(trg_hdr)

    hdr_out = copy.copy(hdr)
    _harmonize_restfreq(hdr_out, trg_hdr)
    data, _ = reproject_interp((data, hdr_out), trg_hdr, order="nearest-neighbor")

    if not is_cube:
        data = np.broadcast_to(data, (ov_hdr["NAXIS3"], *data.shape)).copy()

    return data


def get_mom_maps_ppv(cube, mask, vaxis, mom_calc):
    """
    Compute moment maps directly on a PPV cube, reusing utils_table.get_mom_maps
    exactly as-is.

    get_mom_maps expects a (n_pts, n_chan) array (one row per hex-grid
    point). A PPV cube is (n_chan, ny, nx). Rather than duplicate or modify
    get_mom_maps, this function reshapes the cube to (ny*nx, n_chan) — i.e.
    treating every pixel as one "point" — calls get_mom_maps unchanged, and
    reshapes the (ny*nx,) results back to (ny, nx) maps. This guarantees the
    PPV moments are computed with the literal same code as the hex-grid
    moments, just on a different (denser, regular) set of "points".

    Parameters
    ----------
    cube     : astropy Quantity (n_chan, ny, nx) — brightness temperature cube
    mask     : array-like (n_chan, ny, nx) — 0/1 integration mask
    vaxis    : astropy Quantity (n_chan,) — velocity axis
    mom_calc : tuple (SN_thresh, conseq_channels, mom2_method)

    Returns
    -------
    dict mapping str -> astropy Quantity (ny, nx): same keys as get_mom_maps
    (rms, tpeak, mom0, mom0_err, mom1, mom1_err, mom2, mom2_err, ew, ew_err)
    """
    n_chan, ny, nx = cube.shape

    # (n_chan, ny, nx) -> (ny, nx, n_chan) -> (ny*nx, n_chan): treat every
    # pixel as one "point", matching get_mom_maps' expected row-major layout.
    cube_pts = np.moveaxis(cube.value, 0, -1).reshape(ny * nx, n_chan) * cube.unit
    mask_pts = np.moveaxis(np.asarray(mask), 0, -1).reshape(ny * nx, n_chan)

    mom_maps_pts = get_mom_maps(cube_pts, mask_pts, vaxis, mom_calc)

    mom_maps = {}
    for key, val in mom_maps_pts.items():
        mom_maps[key] = val.reshape(ny, nx)
    return mom_maps


def save_ppv_mask_to_fits(mask, ov_hdr, source, filename, folder):
    """
    Write a PPV-native velocity-integration mask to a 3-D FITS cube.

    Unlike the hex-grid path's mask regridding (which used to reproject a
    SPEC_MASK* table column onto the overlay grid), the mask here is already a plain numpy array on
    the overlay's native PPV grid — produced directly by construct_mask_ppv
    / build_hfs_mask_ppv / fixed_velocity_mask_ppv / external_mask_ppv — so
    no resampling or re-binarization is needed; it is written out as-is.

    Parameters
    ----------
    mask     : np.ndarray (n_chan, ny, nx) — 0/1 mask array
    ov_hdr   : FITS Header (3-D) — overlay header; supplies both the spatial
              WCS and the spectral axis for the output cube
    source   : str — source name
    filename : str — quantity label used in the output filename, e.g. "mask"
              or "mask_12co21"
    folder   : str — output directory

    Output filename: {source}_{filename}.fits
    """
    hdr_out = copy.copy(ov_hdr)
    fname_fits = os.path.join(folder, f"{source}_{filename}.fits")
    fits.writeto(fname_fits, data=np.asarray(mask, dtype=float), header=hdr_out, overwrite=True)


def run_moments_ppv(source, meta, cubes, input_mask, hfs_data, params, folder,
                    save_mask=False):
    """
    Compute and write PPV-native moment maps for every cube of *source*.

    This function reproduces the mask-construction orchestration of
    stage_products.run_products (ref_line selection, ref_line combination
    modes, ref+HI, strict_mask, HFS extension, use_input_mask /
    use_fixed_vel_mask) exactly, but operates on convolved PPV cubes
    (get_convolved_ppv_cube) and computes moments with get_mom_maps_ppv
    instead of working through the hex-grid .ecsv table.

    Required inputs (raises FileNotFoundError if missing, per cube):
      - the convolved PPV cube, either the stage_regrid save_fits output or
        (as a fallback) the raw input cube to convolve from scratch — see
        get_convolved_ppv_cube.
      - the overlay cube (for the WCS / spectral axis / footprint).

    Parameters
    ----------
    source     : str
    meta       : dict — from KeyHandler.meta
    cubes      : pd.DataFrame — cube definitions from KeyHandler
    input_mask : pd.DataFrame — mask definition from KeyHandler
    hfs_data   : pd.DataFrame or None — hyperfine data from KeyHandler
    params     : dict — source geometry from SourceHandler
    folder     : str — output directory for the moment FITS files
    save_mask  : bool — if True, also write the PPV mask(s) used here to FITS
                (see save_ppv_mask_to_fits): the combined mask once as
                {source}_mask.fits, plus one {source}_mask_<line>.fits for
                every line whose HFS-extended mask actually differs from
                the combined mask.
    """
    use_input_mask     = meta.get("use_input_mask",     False)
    use_fixed_vel_mask = meta.get("use_fixed_vel_mask", False)
    use_mask           = use_input_mask or use_fixed_vel_mask
    if input_mask is None:
        input_mask = []
    use_hfs_lines      = meta.get("use_hfs_lines",      False)
    strict_mask        = meta.get("strict_mask",        False)
    ref_line_method    = meta.get("ref_line",           "first")
    SN_processing      = meta.get("SN_processing",      [2, 4])
    mom_calc           = [meta.get("mom_thresh",        5),
                          meta.get("conseq_channels",   3),
                          meta.get("mom2_method",       "fwhm")]
    target_res_as       = _resolve_target_res(params, meta)
    data_dir            = meta.get("data_dir", "data/")

    line_names = [str(l) for l in cubes["line_name"]]
    n_lines    = len(line_names)

    ref_line = (ref_line_method.upper()
                if ref_line_method in line_names
                else line_names[0].upper())

    # ------------------------------------------------------------------
    # Load the overlay cube to get the reference WCS/spectral axis.
    # ------------------------------------------------------------------
    overlay_file  = meta.get("overlay_file", "")
    overlay_fname = (os.path.join(data_dir, overlay_file) if source in overlay_file
                     else os.path.join(data_dir, source + overlay_file))
    if not os.path.exists(overlay_fname):
        LOG.error(f"Overlay file not found: {overlay_fname}")
        raise FileNotFoundError(f"Overlay file not found: {overlay_fname}")
    _, ov_hdr = fits.getdata(overlay_fname, header=True)
    ov_hdr, _ = _ensure_ms(copy.copy(ov_hdr))

    delta_v_kms = (ov_hdr["CDELT3"] * au.Unit(ov_hdr.get("CUNIT3", "m/s"))).to(au.km / au.s).value
    vaxis = (_get_vaxis(ov_hdr) * au.Unit(ov_hdr.get("CUNIT3", "m/s"))).to(au.km / au.s)

    # ------------------------------------------------------------------
    # Load every cube's convolved PPV data up front (needed both for mask
    # construction and the per-line moment computation below).
    # ------------------------------------------------------------------
    cube_data = {}
    for _, row in cubes.iterrows():
        name = str(row["line_name"])
        data, _ = get_convolved_ppv_cube(
            source, name, str(row["line_dir"]), str(row["line_ext"]),
            target_res_as, ov_hdr, folder, log=LOG,
        )
        cube_data[name.upper()] = data * au.Unit(str(row["line_unit"]))

    if ref_line.upper() not in cube_data:
        LOG.error(f"Reference line {ref_line} not found among loaded cubes for {source}.")
        raise FileNotFoundError(f"Reference line {ref_line} not found among loaded cubes for {source}.")

    # ------------------------------------------------------------------
    # Mask construction — mirrors stage_products.run_products exactly.
    # ------------------------------------------------------------------
    if use_mask:
        if len(input_mask) == 0:
            LOG.error("use_mask is True but no mask defined in config.txt.")
            raise ValueError("use_mask is True but no mask defined in config.txt.")

        if use_fixed_vel_mask:
            mask_unit  = input_mask["mask_unit"].iloc[0]
            mask_start = float(input_mask["mask_start"].iloc[0]) * au.Unit(mask_unit)
            mask_end   = float(input_mask["mask_end"].iloc[0])   * au.Unit(mask_unit)
            mask = fixed_velocity_mask_ppv(cube_data[ref_line].shape, ov_hdr,
                                           mask_start, mask_end, mask_unit)
            LOG.info(f"Fixed velocity mask applied ({mask_start} to {mask_end}).")
        else:
            mask_file = os.path.join(
                str(input_mask["mask_dir"].iloc[0]),
                source + str(input_mask["mask_ext"].iloc[0]),
            )
            if not os.path.exists(mask_file):
                LOG.error(f"Mask file not found: {mask_file}")
                raise FileNotFoundError(f"Mask file not found: {mask_file}")
            mask = external_mask_ppv(mask_file, ov_hdr, log=LOG)
            LOG.info("External mask sampled onto PPV grid.")
    else:
        LOG.info(f"Building PPV velocity mask from {ref_line}.")
        mask = construct_mask_ppv(cube_data[ref_line].value, SN_processing)

        if ref_line_method == "all":
            n_mask = n_lines
        elif isinstance(ref_line_method, int):
            n_mask = min(n_lines, ref_line_method)
        else:
            n_mask = 0   # "first": only the reference line

        for n_mask_i in range(1, n_mask + 1):
            line_i = line_names[n_mask_i].upper()
            if line_i not in cube_data:
                continue
            mask_i = construct_mask_ppv(cube_data[line_i].value, SN_processing)
            mask = ((mask.astype(int) | mask_i.astype(int)))
            LOG.info(f"Combined PPV mask includes {line_i}.")

        if ref_line_method == "ref+HI":
            if "HI" in cube_data:
                mask_hi = construct_mask_ppv(cube_data["HI"].value, SN_processing)
                mask = (mask.astype(int) | mask_hi.astype(int))
                LOG.info("ref+HI mask: combined reference line and HI masks.")
            else:
                LOG.warning("HI not found among loaded cubes; ignoring ref+HI option.")

        if strict_mask:
            LOG.info("Applying strict spatial mask filter (connected-component, PPV grid).")
            mask = apply_strict_mask_ppv(mask.astype(int))

    if save_mask:
        save_ppv_mask_to_fits(mask, ov_hdr, source, "mask", folder)
        LOG.info(f"PPV mask cube written to: {folder}")

    # ------------------------------------------------------------------
    # Compute and write moments for every line.
    # ------------------------------------------------------------------
    for jj, row in cubes.iterrows():
        line_name = str(row["line_name"])
        if line_name.upper() not in cube_data:
            continue

        active_mask = mask
        if use_hfs_lines and hfs_data is not None:
            mask_hfs = build_hfs_mask_ppv(mask, line_name, hfs_data, delta_v_kms)
            if mask_hfs is not None:
                active_mask = mask_hfs
                LOG.info(f"Using HFS-extended PPV mask for {line_name}.")
                if save_mask and not np.array_equal(mask_hfs, mask):
                    save_ppv_mask_to_fits(mask_hfs, ov_hdr, source,
                                         f"mask_{line_name.lower()}", folder)
                    LOG.info(f"PPV mask cube for {line_name} written to: {folder}")

        mom_maps = get_mom_maps_ppv(cube_data[line_name.upper()], active_mask, vaxis, mom_calc)

        ov_hdr_2d  = twod_head(ov_hdr)
        line_desc  = str(row["line_desc"])
        line_unit  = str(row["line_unit"])

        quantities = {
            "mom0":  (mom_maps["mom0"],     "K km/s" if line_unit == "K" else f"{line_unit} km/s"),
            "emom0": (mom_maps["mom0_err"], None),
            "mom1":  (mom_maps["mom1"],     "km/s"),
            "emom1": (mom_maps["mom1_err"], "km/s"),
            "mom2":  (mom_maps["mom2"],     "km/s"),
            "emom2": (mom_maps["mom2_err"], "km/s"),
            "tpeak": (mom_maps["tpeak"],    line_unit),
            "rms":   (mom_maps["rms"],      line_unit),
            "ew":    (mom_maps["ew"],       "km/s"),
            "eew":   (mom_maps["ew_err"],   "km/s"),
        }
        for quantity, (arr, bunit) in quantities.items():
            hdr_out = copy.copy(ov_hdr_2d)
            if bunit:
                hdr_out["BUNIT"] = bunit
            hdr_out["LINE"] = line_name
            fname_fits = os.path.join(folder, f"{source}_{line_name}_{quantity}.fits")
            data_out = arr.value if hasattr(arr, "value") else arr
            fits.writeto(fname_fits, data=data_out, header=hdr_out, overwrite=True)

        LOG.info(f"Compute moment maps and write to file for line: {line_name}.")


def run_fits(source, fname, meta, maps, cubes, params, input_mask=None, hfs_data=None):
    """
    Write FITS moment maps, 2D map images, and mask cube(s) for *source*.

    This is the entry point for the "fits" pipeline stage.

    Moment maps (if save_mom_maps is True) are computed PPV-native: directly
    on the convolved, overlay-aligned PPV cubes, never via the hex-grid
    .ecsv table — see run_moments_ppv and the module docstring. This
    requires either a save_fits cube on disk from stage_regrid, or the raw
    input cube to convolve from scratch as a fallback; it raises
    FileNotFoundError if neither is available for a given line.

    2D map images (if save_maps is True) are still regridded from the
    hex-grid .ecsv table via save_to_fits, since 2D map columns have no PPV
    cube equivalent. Mask cube(s) (if save_mask is True) are now written
    PPV-native too, as a byproduct of run_moments_ppv: the combined mask
    once as {source}_mask.fits, plus one {source}_mask_<line>.fits for
    every line whose HFS-extended mask differs from the combined mask. This
    means save_mask now requires save_mom_maps to also be True (a warning
    is logged if save_mask is requested without save_mom_maps).

    Parameters
    ----------
    source     : str
    fname      : str          — path to the processed .ecsv file
    meta       : dict         — from KeyHandler.meta
    maps       : pd.DataFrame — map definitions from KeyHandler
    cubes      : pd.DataFrame — cube definitions from KeyHandler
    params     : dict         — source geometry from SourceHandler
    input_mask : pd.DataFrame, optional — mask definition from KeyHandler
                (required if use_input_mask or use_fixed_vel_mask is set)
    hfs_data   : pd.DataFrame or None, optional — hyperfine data from KeyHandler
                (required if use_hfs_lines is set)
    """
    save_mom_maps    = meta.get("save_mom_maps",    True)
    save_maps        = meta.get("save_maps",        True)
    save_mask        = meta.get("save_mask",        False)
    folder           = meta.get("folder_savefits",  "./saved_fits_files/")
    target_res_as    = _resolve_target_res(params, meta)
    pixels_per_beam = meta.get("pixels_per_beam", 2.0)

    if not (save_mom_maps or save_maps or save_mask):
        LOG.info(f"Output writing disabled for {source}; skipping.")
        return

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
    # Moment maps — PPV-native, NOT from the hex-grid .ecsv table.
    # The PPV mask(s) used here are also written out (if save_mask is True)
    # as a byproduct of this same call, since the mask only exists as a
    # plain array inside run_moments_ppv.
    # ------------------------------------------------------------------
    if save_mom_maps:
        run_moments_ppv(source, meta, cubes, input_mask, hfs_data, params, folder,
                        save_mask=save_mask)
        LOG.info(f"Moment map FITS files written to: {folder}")
    elif save_mask:
        LOG.warning(f"save_mask is True but save_mom_maps is False for {source}; "
                   "the PPV mask is only built while computing moments, so no "
                   "mask FITS file(s) will be written. Set save_mom_maps = true "
                   "to enable mask output.")

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
