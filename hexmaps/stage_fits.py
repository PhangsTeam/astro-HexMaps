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
     by reading it from disk (a cube previously saved by the fits stage when save_cubes is
     True) or, if that file is absent, by convolving and reprojecting the
     raw input cube itself (convolve_cube_to_target / reproject_cube_to_overlay).
  2. Builds (or reads) a PPV mask using exactly the same two-level S/N
     threshold + dilation algorithm as the hex-grid path (construct_mask_ppv
     mirrors stage_products.construct_mask channel-by-channel on the cube),
     including the same ref_line combination logic, the strict spatial
     filter (here implemented as 2-D connected-component labelling per
     channel, the rectangular-grid equivalent of the hex-grid distance-based
     filter), HFS mask extension, and the use_input_mask / use_fixed_vel_mask
     external-mask options. Before masking, two additional steps constrain
     the valid pixel area:
       a. Overlay footprint masking: pixels where the overlay cube has no
          finite values along the entire velocity axis (outside the observed
          area) are set to NaN in all convolved cubes, so the footprint for
          edge erosion reflects actual data coverage, not just the
          reprojected grid extent.
       b. Edge erosion: the footprint is eroded by half a beam width
          (FWHM / 2) using a circular structuring element to remove
          convolution-artefact pixels at the map boundary.
  3. Computes moments on the masked PPV cube by reshaping it to
     (n_pix, n_chan) and calling utils_table.get_mom_maps — the exact same
     function used by the hex-grid path — then reshaping the results back to
     (ny, nx) maps.
  4. Writes one FITS file per moment quantity per line.

2D maps: also PPV-native; mask cube(s): PPV-native too
---------------------------------------------------------
2D band/map columns (MAP_*/EMAP_*) are now also processed PPV-native,
mirroring the cube path exactly: the raw input FITS file is read, convolved
to the target resolution on its native pixel grid with conv_with_gauss, and
reprojected onto the overlay's 2-D spatial WCS with bilinear interpolation
(get_convolved_map). This replaces the old nearest-neighbour hex-grid
regridding (save_to_fits/sample_to_hdr) which produced blocky artefacts.

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
from scipy.ndimage import label, binary_erosion
from skimage.morphology import disk
from datetime import date

from hexmaps.utils_fits import twod_head, conv_with_gauss, reproject_cube, resolve_meta_resolution
from hexmaps.stage_regrid import _ensure_ms, _get_vaxis
from hexmaps.utils_table import get_mom_maps, build_noise_mask

from hexmaps import __version__ as _HEXMAPS_VERSION
from hexmaps.logger import get_logger

LOG = get_logger("FITS")


# ============================================================================
# Grid helpers
# ============================================================================


def get_convolved_ppv_cube(
    source, line_name, line_dir, line_ext, meta, ov_hdr, log=None
):
    """
    Convolve the raw input cube for *line_name* to the target resolution and
    reproject it onto the overlay WCS.

    Always reads the raw input FITS file and performs the convolution and
    reprojection from scratch. There is no cache lookup — reproducibility is
    guaranteed by always starting from the original data.

    Parameters
    ----------
    source    : str — source name, prepended to line_ext to form the filename
    line_name : str — line label (used in log messages)
    line_dir  : str — directory containing the raw input FITS file
    line_ext  : str — filename extension of the raw input file
    meta      : dict — pipeline settings (target_res, res_suffix, etc.)
    ov_hdr    : FITS Header — overlay header defining the target WCS
    log       : StageLogger, optional

    Returns
    -------
    data : np.ndarray (n_chan, ny, nx) — convolved cube on the overlay grid
    hdr  : FITS Header — header matching the output array
    """
    log = log or LOG

    raw_path = os.path.join(line_dir, source + line_ext)
    log.info(f"Convolving {line_name} from: {raw_path}")
    if not os.path.exists(raw_path):
        log.error(f"Raw input cube not found for line: {line_name}: {raw_path}")
        raise FileNotFoundError(
            f"Raw input cube not found for line: {line_name}: {raw_path}"
        )

    data, hdr = fits.getdata(raw_path, header=True)
    data, hdr = convolve_cube_to_target(
        data, hdr, meta.get("target_res", 27.0), log=log
    )
    data, hdr = reproject_cube_to_overlay(data, hdr, ov_hdr, log=log)
    return data, hdr


def convolve_cube_to_target(data, hdr, target_res_as, log=None):
    """
    Convolve a PPV cube to *target_res_as* arcsec, in place on its native grid.
    """
    log = log or LOG
    if "BMAJ" not in hdr:
        log.warning("No BMAJ in header; skipping convolution.")
        return data, hdr
    if hdr["BMAJ"] >= 0.99 * target_res_as / 3600.0:
        log.info("Cube already at or above target resolution; skipping convolution.")
        return data, hdr

    data, hdr_out = conv_with_gauss(
        in_data=data,
        in_hdr=hdr,
        target_beam=target_res_as * np.array([1.0, 1.0, 0.0]),
        quiet=True,
        log=log,
    )
    return data, hdr_out


def reproject_cube_to_overlay(data, hdr, ov_hdr, log=None):
    """
    Reproject a PPV cube onto the overlay's spatial+spectral WCS.
    """
    log = log or LOG
    trg_hdr = copy.deepcopy(ov_hdr)
    trg_hdr, _ = _ensure_ms(trg_hdr)
    hdr, data = _ensure_ms(copy.copy(hdr), data)

    data, _ = reproject_cube((data, hdr), trg_hdr, order="bilinear")
    return data, trg_hdr


def reproject_map_to_overlay(data, hdr, ov_hdr, log=None):
    """
    Reproject a 2-D map onto the overlay's spatial WCS.

    Uses bilinear interpolation (same as the cube path) to place the map on
    the overlay pixel grid.  The overlay header is collapsed to 2-D first so
    spectral keywords don't confuse the reprojection.

    Parameters
    ----------
    data   : np.ndarray (ny, nx)
    hdr    : FITS Header — native 2-D header of *data*
    ov_hdr : FITS Header — overlay header (2-D or 3-D) defining the target WCS
    log    : StageLogger, optional

    Returns
    -------
    data, hdr : reprojected map and the 2-D overlay-aligned header
    """
    log = log or LOG
    from hexmaps.utils_fits import twod_head

    trg_hdr = twod_head(copy.deepcopy(ov_hdr))
    data, _ = reproject_cube((data, hdr), trg_hdr, order="bilinear")
    return data, trg_hdr


def get_convolved_map(
    source, map_name, map_dir, map_ext, target_res_as, ov_hdr, log=None
):
    """
    Obtain a convolved 2-D map for *map_name*, reprojected onto the overlay WCS.

    Mirrors get_convolved_ppv_cube but for 2-D maps: reads the raw input
    FITS file, convolves it to *target_res_as* using conv_with_gauss, and
    reprojects it onto the overlay's 2-D spatial WCS with bilinear
    interpolation.

    Unlike the cube path there is no on-disk cache for maps (stage_regrid
    does not write a save_cubes copy for 2-D maps), so the convolution and
    reprojection always run from the raw input file.

    Parameters
    ----------
    source        : str
    map_name      : str   — map name, e.g. "spire250"
    map_dir       : str   — directory containing the raw input FITS file
    map_ext       : str   — filename extension of the raw input file
    target_res_as : float — target beam FWHM in arcseconds
    ov_hdr        : FITS Header — overlay header (3-D or 2-D) defining the target WCS
    log           : StageLogger, optional

    Returns
    -------
    data : np.ndarray (ny, nx) — convolved map on the overlay spatial grid
    hdr  : FITS Header — 2-D header matching *data*
    """
    log = log or LOG

    raw_path = os.path.join(map_dir, source + map_ext)
    if not os.path.exists(raw_path):
        log.error(f"Raw input map not found for {map_name}: {raw_path}")
        raise FileNotFoundError(f"Raw input map not found for {map_name}: {raw_path}")

    log.info(f"Processing map {map_name} from: {raw_path}")
    data, hdr = fits.getdata(raw_path, header=True)

    # Convolve to target resolution on the native grid
    if "BMAJ" not in hdr:
        log.warning(f"No BMAJ in header for {map_name}; skipping convolution.")
    elif hdr["BMAJ"] >= 0.99 * target_res_as / 3600.0:
        log.info(
            f"Map {map_name} already at or above target resolution; skipping convolution."
        )
    else:
        data, hdr = conv_with_gauss(
            in_data=data,
            in_hdr=hdr,
            target_beam=target_res_as * np.array([1.0, 1.0, 0.0]),
            quiet=True,
            log=log,
        )

    # Reproject onto the overlay 2-D spatial grid
    data, hdr = reproject_map_to_overlay(data, hdr, ov_hdr, log=log)
    return data, hdr


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
    rms = median_absolute_deviation(
        ref_cube[np.where(ref_cube < 3 * rms)], ignore_nan=True
    )

    # Per-pixel noise: MAD of channels below the global 3-sigma threshold.
    # axis=0 is the spectral axis here (vs axis=1 for the hex-grid table).
    mask_rough = ref_cube < 3 * rms
    masked_cube = np.where(mask_rough, ref_cube, np.nan)
    med_mask = np.nanmedian(masked_cube, axis=0)
    mad_mask = np.nanmedian(np.abs(masked_cube - med_mask[None, :, :]), axis=0)

    low_thresh = SN_processing[0] * mad_mask[None, :, :]
    high_thresh = SN_processing[1] * mad_mask[None, :, :]

    # Initial high-S/N mask: channel above high_thresh with adjacent support.
    # np.roll along axis=0 (spectral) replaces axis=1 from the hex-grid version.
    mask = (ref_cube > high_thresh).astype(int)
    low_mask = (ref_cube > low_thresh).astype(int)
    mask = mask & (np.roll(mask, 1, 0) | np.roll(mask, -1, 0))

    # Require >= 3 of 3 consecutive channels to suppress single-channel spikes
    mask = ((mask + np.roll(mask, 1, 0) + np.roll(mask, -1, 0)) >= 3).astype(int)
    low_mask = (
        (low_mask + np.roll(low_mask, 1, 0) + np.roll(low_mask, -1, 0)) >= 3
    ).astype(int)

    # Dilate high-S/N core into low-S/N wings (5 passes)
    for _ in range(5):
        mask = ((mask + np.roll(mask, 1, 0) + np.roll(mask, -1, 0)) >= 1).astype(
            int
        ) * low_mask

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
        small_labels = small_labels[small_labels != 0]  # never touch background
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

    idx_cols = hfs_data["hfs_name"] == line_name
    restfreqs = [
        f * au.Unit(str(u))
        for f, u in zip(hfs_data["hfs_ref_freq"][idx_cols], hfs_data["unit"][idx_cols])
    ]
    hfs_freqs = [
        f * au.Unit(str(u))
        for f, u in zip(hfs_data["hfs_freq"][idx_cols], hfs_data["unit"][idx_cols])
    ]

    mask_hfs = mask.copy()
    for freq, restfreq in zip(hfs_freqs, restfreqs):
        v_shift = freq.to(au.km / au.s, equivalencies=au.doppler_radio(restfreq))
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
    is_cube = data.ndim == 3

    trg_hdr = copy.deepcopy(ov_hdr)
    if not is_cube:
        trg_hdr = twod_head(trg_hdr)

    hdr_out = copy.copy(hdr)
    data, _ = reproject_cube((data, hdr_out), trg_hdr, order="nearest-neighbor")

    if not is_cube:
        data = np.broadcast_to(data, (ov_hdr["NAXIS3"], *data.shape)).copy()

    return data


def get_mom_maps_ppv(cube, mask, vaxis, mom_calc, noise_mask=None):
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
    cube       : astropy Quantity (n_chan, ny, nx) — brightness temperature cube
    mask       : array-like (n_chan, ny, nx) — 0/1 integration mask
    vaxis      : astropy Quantity (n_chan,) — velocity axis
    mom_calc   : tuple (SN_thresh, conseq_channels, mom2_method)
    noise_mask : array-like (n_chan, ny, nx), optional — channels to use for
                 noise estimation (built by build_noise_mask; see get_mom_maps
                 for details)

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
    noise_pts = (
        np.moveaxis(np.asarray(noise_mask), 0, -1).reshape(ny * nx, n_chan)
        if noise_mask is not None
        else None
    )

    mom_maps_pts = get_mom_maps(
        cube_pts, mask_pts, vaxis, mom_calc, noise_mask=noise_pts
    )

    mom_maps = {}
    for key, val in mom_maps_pts.items():
        mom_maps[key] = val.reshape(ny, nx)
    return mom_maps


def make_clean_header(ov_hdr, is_cube=False, bunit=None, btype=None,
                      line_name=None, line_desc=None, object_name=None,
                      bmaj_as=None, restfrq=None, meta=None):
    """
    Build a clean FITS header containing only the mandatory keywords needed
    to locate and interpret the data — no residual keywords from the input
    overlay file.

    Starting from scratch (not from a copy of ov_hdr) guarantees that
    instrument-specific, history, comment, and pipeline-internal keywords
    from the original input do not bleed into the output.

    Spatial WCS keywords (CTYPE/CRVAL/CRPIX/CDELT/CUNIT for axes 1 and 2)
    are copied from *ov_hdr* unchanged.  If the header has a projection
    keyword EQUINOX or RADESYS those are copied too.

    For 3-D cubes (is_cube=True) the spectral axis keywords (axis 3 plus
    SPECSYS / VELREF) are also copied from *ov_hdr*.  RESTFRQ is taken from
    the *restfrq* argument when provided — which should be the value from the
    raw input cube header — rather than from the overlay, whose RESTFRQ would
    be wrong for any line other than the overlay line itself.

    Content and provenance keywords follow the order:
        OBJECT, LINE, BUNIT, BTYPE, ORIGIN, AUTHOR, DATE
    followed by a COMMENT block crediting the HexMaps package.

    Parameters
    ----------
    ov_hdr      : FITS Header — overlay header (2-D or 3-D) used as WCS source
    is_cube     : bool — True for PPV cubes, False for 2-D images
    bunit       : str or None — value for BUNIT
    btype       : str or None — value for BTYPE (moment type); written in upper
                  case (e.g. "MOM0" not "mom0")
    line_name   : str or None — short line identifier; only used as a fallback
                  for LINE when line_desc is not provided
    line_desc   : str or None — human-readable line description written to LINE
                  (preferred over line_name; e.g. "12CO(2-1)" instead of "12co21")
    object_name : str or None — source name for OBJECT keyword
    bmaj_as     : float or None — beam FWHM in arcsec; overrides ov_hdr BMAJ/BMIN.
                  When None the beam is copied from ov_hdr if present.
    restfrq     : float or None — rest frequency in Hz for the output cube.
                  Should be read from the raw input cube header, not the overlay.
                  When None, RESTFRQ is copied from ov_hdr if present (cubes only).
    meta        : dict or None — pipeline meta dict; used for AUTHOR

    Returns
    -------
    astropy.io.fits.Header — clean header ready to pass to fits.writeto
    """
    h = fits.Header()

    # -- Spatial WCS (axes 1 and 2) ------------------------------------------
    for ax in ("1", "2"):
        for kw in ("CTYPE", "CRVAL", "CRPIX", "CDELT", "CUNIT"):
            key = kw + ax
            if key in ov_hdr:
                h[key] = (ov_hdr[key], ov_hdr.comments[key])

    for kw in ("EQUINOX", "RADESYS", "LONPOLE", "LATPOLE"):
        if kw in ov_hdr:
            h[kw] = (ov_hdr[kw], ov_hdr.comments[kw])

    # -- Spectral axis (axis 3) — only for cubes -----------------------------
    if is_cube:
        for kw in ("CTYPE3", "CRVAL3", "CRPIX3", "CDELT3", "CUNIT3"):
            if kw in ov_hdr:
                h[kw] = (ov_hdr[kw], ov_hdr.comments[kw])
        for kw in ("SPECSYS", "VELREF"):
            if kw in ov_hdr:
                h[kw] = (ov_hdr[kw], ov_hdr.comments[kw])
        # RESTFRQ from the input cube, not the overlay.
        # Also check the older RESTFREQ keyword (no trailing Q) for
        # backwards compatibility with pre-FITS-4 files.
        if restfrq is not None:
            h["RESTFRQ"] = restfrq
        elif "RESTFRQ" in ov_hdr:
            h["RESTFRQ"] = ov_hdr["RESTFRQ"]
        elif "RESTFREQ" in ov_hdr:
            h["RESTFRQ"] = ov_hdr["RESTFREQ"]

    # -- Beam ----------------------------------------------------------------
    if bmaj_as is not None:
        h["BMAJ"] = bmaj_as / 3600.0
        h["BMIN"] = bmaj_as / 3600.0
        h["BPA"]  = 0.0
    else:
        for kw in ("BMAJ", "BMIN", "BPA"):
            if kw in ov_hdr:
                h[kw] = (ov_hdr[kw], ov_hdr.comments[kw])

    # -- Content and provenance — fixed order: OBJECT LINE BUNIT BTYPE -------
    #    ORIGIN AUTHOR DATE
    if object_name is not None:
        h["OBJECT"] = object_name.upper()
    # LINE stores the human-readable description (line_desc) when available,
    # falling back to the short identifier (line_name) otherwise.
    line_value = line_desc if line_desc else line_name
    if line_value is not None:
        h["LINE"] = line_value
    if bunit is not None:
        h["BUNIT"] = bunit
    if btype is not None:
        h["BTYPE"] = str(btype).upper()

    h["ORIGIN"] = "HexMaps"
    if meta is not None:
        author = meta.get("user", "")
        if author:
            h["AUTHOR"] = author
    h["DATE"] = date.today().isoformat()

    # -- Package credit comment ----------------------------------------------
    h["COMMENT"] = "Created with HexMaps (HEXagonal-grid Multi-data"
    h["COMMENT"] = f"Analysis and Processing Software) version {_HEXMAPS_VERSION}"
    h["COMMENT"] = "https://github.com/PhangsTeam/HexMaps"
    h["COMMENT"] = "Contact: Jakob den Brok <jadenbrok@mpia.de>"
    h["COMMENT"] = "         Lukas Neumann  <lukas.neumann@eso.org>"

    return h


def save_ppv_mask_to_fits(mask, ov_hdr, source, filename, folder,
                          out_nan_mask=None, meta=None):
    """
    Write a PPV-native velocity-integration mask to a 3-D FITS cube.

    Unlike the hex-grid path's mask regridding (which used to reproject a
    SPEC_MASK* table column onto the overlay grid), the mask here is already
    a plain numpy array on the overlay's native PPV grid — produced directly
    by construct_mask_ppv / build_hfs_mask_ppv / fixed_velocity_mask_ppv /
    external_mask_ppv — so no resampling or re-binarization is needed.

    Parameters
    ----------
    mask         : np.ndarray (n_chan, ny, nx) — 0/1 mask array
    ov_hdr       : FITS Header (3-D) — overlay header; supplies both the spatial
                  WCS and the spectral axis for the output cube
    source       : str — source name
    filename     : str — quantity label used in the output filename, e.g. "mask"
                  or "mask_12co21"
    folder       : str — output directory
    out_nan_mask : np.ndarray (ny, nx) bool, optional — if supplied, pixels
                  where out_nan_mask is True are set to NaN across all channels,
                  matching the NaN pattern of the moment maps and convolved cubes.

    Output filename: {source}_{filename}.fits
    """
    hdr_out = make_clean_header(
        ov_hdr, is_cube=True, bunit="", btype="mask",
        object_name=source, meta=None,
    )
    data_out = np.asarray(mask, dtype=float).copy()
    if out_nan_mask is not None:
        data_out[:, out_nan_mask] = np.nan
    fname_fits = os.path.join(folder, f"{source}_{filename}.fits")
    fits.writeto(fname_fits, data=data_out, header=hdr_out, overwrite=True)


def build_edge_mask(ov_footprint, ov_hdr, target_res_as, fov_erosion_beams=0.5):
    """
    Build a 2-D spatial edge mask by eroding the observed non-NaN footprint.

    Convolution near the edge of the observed area is unreliable: the beam
    extends beyond the data boundary and the convolved values are computed
    from only a partial beam footprint, which systematically under- or
    over-estimates the true brightness. Removing a border of
    ``fov_erosion_beams × beam_FWHM`` in pixels eliminates the worst-affected
    region. The default of 0.5 (half a beam) is the conventional minimum
    safe margin; larger values give a more conservative trim.

    The footprint to erode is the *non-NaN area* of the overlay cube —
    the irregular blob of pixels where the overlay actually has data — not
    the rectangular reprojected grid extent (which includes pixels filled in
    by interpolation beyond the true observed boundary). Eroding the true
    non-NaN blob means the edge-removal correctly follows the shape of the
    observed field, including any holes or concavities.

    The algorithm:
      1. Accept *ov_footprint*: a 2-D bool array, True where the overlay has
         at least one finite value along the velocity axis.
      2. Compute the trim radius in pixels:
         floor(fov_erosion_beams × target_res_as / pixel_scale).
         Uses floor to be conservative. If fov_erosion_beams is 0 (or
         the resulting radius is < 1 pixel), no erosion is applied and
         the full footprint is returned.
      3. Erode *ov_footprint* using a circular structuring element of that
         radius via scipy.ndimage.binary_erosion. A circular (rather than
         square) kernel gives isotropic edge removal matching the shape of
         the Gaussian PSF.
      4. Return the eroded footprint as a float array (1 inside, 0 on edges).

    Parameters
    ----------
    ov_footprint      : np.ndarray (ny, nx) bool — True where the overlay has
                        at least one finite value along the velocity axis.
                        Derived in the caller as
                        ``np.any(np.isfinite(ov_data), axis=0)``.
    ov_hdr            : FITS Header (3-D) — overlay header; provides CDELT1 for the
                        pixel scale used to convert target_res_as to pixels
    target_res_as     : float — target beam FWHM in arcseconds
    fov_erosion_beams : float — erosion radius in units of the beam FWHM
                        (default 0.5 = half beam). Set to 0 to disable erosion.

    Returns
    -------
    edge_mask : np.ndarray (ny, nx) — float 0/1 array; multiply into the
                spatial dimension of any mask to remove convolution-edge pixels
    """
    # Trim radius in pixels
    pix_scale_as = abs(ov_hdr["CDELT1"]) * 3600.0  # arcsec/pixel
    trim_radius_pix = int(np.floor(fov_erosion_beams * target_res_as / pix_scale_as))

    if trim_radius_pix <= 0:
        if fov_erosion_beams > 0:
            LOG.warning("Edge trim radius is <= 0 pixels; no edge removal applied.")
        else:
            LOG.info("FOV erosion disabled (fov_erosion_beams = 0).")
        return ov_footprint.astype(float)

    trim_as = fov_erosion_beams * target_res_as
    LOG.info(
        f"Removing {trim_radius_pix} pixel edge border "
        f"(= {fov_erosion_beams} beam = {trim_as:.1f} arcsec at "
        f"{pix_scale_as:.2f} arcsec/px)."
    )

    # Circular structuring element for isotropic erosion of the non-NaN blob
    structure = disk(trim_radius_pix)

    eroded = binary_erosion(ov_footprint, structure=structure)
    return eroded.astype(float)


def run_moments_ppv(
    source,
    meta,
    cubes,
    input_mask,
    hfs_data,
    params,
    folder,
    save_mask=False,
    noise_mask_df=None,
):
    """
    Compute and write PPV-native moment maps for every cube of *source*.

    This function reproduces the mask-construction orchestration of
    stage_products.run_products (ref_line selection, ref_line combination
    modes, ref+HI, strict_mask, HFS extension, use_input_mask /
    use_fixed_vel_mask) exactly, but operates on convolved PPV cubes
    (get_convolved_ppv_cube) and computes moments with get_mom_maps_ppv
    instead of working through the hex-grid .ecsv table.

    Required inputs (raises FileNotFoundError if missing, per cube):
      - the convolved PPV cube, either a previously saved save_cubes output or
        (as a fallback) the raw input cube to convolve from scratch — see
        get_convolved_ppv_cube.
      - the overlay cube (for the WCS / spectral axis / footprint).

    Parameters
    ----------
    source        : str
    meta          : dict — from KeyHandler.meta
    cubes         : pd.DataFrame — cube definitions from KeyHandler
    input_mask    : pd.DataFrame — mask definition from KeyHandler
    hfs_data      : pd.DataFrame or None — hyperfine data from KeyHandler
    params        : dict — source geometry from SourceHandler
    folder        : str — output directory for the moment FITS files
    save_mask     : bool — if True, also write the PPV mask(s) used here to FITS
    noise_mask_df : pd.DataFrame or None — noise velocity window table from
                   KeyHandler.get_noise_mask(). When use_fixed_noise_mask is
                   True and this DataFrame is non-empty, the RMS noise is
                   estimated from channels within these windows rather than
                   from channels outside the integration mask.
    """
    use_input_mask = meta.get("use_input_mask", False)
    use_fixed_vel_mask = meta.get("use_fixed_vel_mask", False)
    use_mask = use_input_mask or use_fixed_vel_mask
    if input_mask is None:
        input_mask = []
    use_hfs_lines = meta.get("use_hfs_lines", False)
    strict_mask = meta.get("strict_mask", False)
    ref_line_method = meta.get("ref_line", "first")
    SN_processing = meta.get("SN_processing", [2, 4])
    mom_calc = [
        meta.get("mom_thresh", 5),
        meta.get("conseq_channels", 3),
        meta.get("mom2_method", "fwhm"),
    ]
    data_dir = meta.get("data_dir", "data/")

    line_names = [str(l) for l in cubes["line_name"]]
    n_lines = len(line_names)

    ref_line = (
        ref_line_method.upper()
        if ref_line_method in line_names
        else line_names[0].upper()
    )

    # ------------------------------------------------------------------
    # Load the overlay cube to get the reference WCS/spectral axis.
    # ------------------------------------------------------------------
    overlay_file = meta.get("overlay_file", "")
    overlay_fname = (
        os.path.join(data_dir, overlay_file)
        if source in overlay_file
        else os.path.join(data_dir, source + overlay_file)
    )
    if not os.path.exists(overlay_fname):
        LOG.error(f"Overlay file not found: {overlay_fname}")
        raise FileNotFoundError(f"Overlay file not found: {overlay_fname}")
    #LOG.info(f"Overlay file: {overlay_fname}")
    ov_data, ov_hdr = fits.getdata(overlay_fname, header=True)
    ov_hdr, _ = _ensure_ms(copy.copy(ov_hdr))

    target_res_as = meta["target_res"]

    delta_v_kms = (
        (ov_hdr["CDELT3"] * au.Unit(ov_hdr.get("CUNIT3", "m/s"))).to(au.km / au.s).value
    )
    vaxis = (_get_vaxis(ov_hdr) * au.Unit(ov_hdr.get("CUNIT3", "m/s"))).to(au.km / au.s)

    # ------------------------------------------------------------------
    # Build the noise channel mask (PPV path).
    # If use_fixed_noise_mask is True and noise_mask_df is non-empty,
    # build a (n_chan, ny, nx) boolean mask selecting the noise channels.
    # This is passed to get_mom_maps_ppv so RMS is estimated from those
    # channels rather than from the inverse of the integration mask.
    # ------------------------------------------------------------------
    use_fixed_noise_mask = meta.get("use_fixed_noise_mask", False)
    ppv_noise_mask = None
    if use_fixed_noise_mask:
        if noise_mask_df is not None and len(noise_mask_df) > 0:
            # Shape determined after cubes are loaded; use a placeholder for
            # now and rebuild with the correct shape inside the per-line loop.
            _noise_mask_df = noise_mask_df
            LOG.info(
                f"Noise RMS will be estimated from {len(noise_mask_df)} "
                "fixed velocity window(s)."
            )
        else:
            LOG.warning(
                "use_fixed_noise_mask is True but no noise_mask rows found "
                "in the [mask] table. Falling back to mask-inverted noise."
            )
            _noise_mask_df = None
    else:
        _noise_mask_df = None

    # ------------------------------------------------------------------
    # Load every cube's convolved PPV data up front (needed both for mask
    # construction and the per-line moment computation below).
    # Also read the raw input headers before convolution so we can copy the
    # correct per-line RESTFRQ into output FITS headers (reprojection
    # overwrites the cube header with the overlay header, losing RESTFRQ).
    # ------------------------------------------------------------------
    cube_data = {}
    cube_hdrs = {}
    for _, row in cubes.iterrows():
        name = str(row["line_name"])
        raw_path = os.path.join(str(row["line_dir"]), source + str(row["line_ext"]))
        cube_hdrs[name.upper()] = (
            fits.getheader(raw_path) if os.path.exists(raw_path) else {}
        )
        data, _ = get_convolved_ppv_cube(
            source,
            name,
            str(row["line_dir"]),
            str(row["line_ext"]),
            meta,
            ov_hdr,
            log=LOG,
        )
        cube_data[name.upper()] = data * au.Unit(str(row["line_unit"]))

    if ref_line.upper() not in cube_data:
        LOG.error(
            f"Reference line {ref_line} not found among loaded cubes for {source}."
        )
        raise FileNotFoundError(
            f"Reference line {ref_line} not found among loaded cubes for {source}."
        )

    # ------------------------------------------------------------------
    # Step 1: constrain every convolved cube to the overlay's valid-pixel
    # footprint. Pixels where the overlay has no finite values along the
    # entire velocity axis (i.e. outside the observed area) are set to NaN
    # in all convolved cubes. This ensures the footprint used for edge
    # erosion reflects the actual data coverage, not just the reprojected
    # grid extent, and propagates NaN correctly through to the moments.
    # ------------------------------------------------------------------
    ov_footprint = np.any(np.isfinite(ov_data), axis=0)  # (ny, nx) bool
    for name in cube_data:
        cube_val = cube_data[name].value.copy()
        cube_val[:, ~ov_footprint] = np.nan
        cube_data[name] = cube_val * cube_data[name].unit

    # ------------------------------------------------------------------
    # Step 2: build a 2-D edge mask by eroding the *overlay's* non-NaN
    # blob (ov_footprint) by half a beam width. Using the overlay footprint
    # directly — rather than deriving a footprint from the convolved cube —
    # ensures the erosion follows the true irregular observed boundary,
    # including any concavities, rather than the rectangular reprojected
    # grid extent.
    # ------------------------------------------------------------------
    edge_mask = build_edge_mask(
        ov_footprint,
        ov_hdr,
        target_res_as,
        fov_erosion_beams=meta.get("fov_erosion_beams", 0.5),
    )

    # ------------------------------------------------------------------
    # Mask construction — mirrors stage_products.run_products exactly.
    # ------------------------------------------------------------------
    if use_mask:
        if len(input_mask) == 0:
            LOG.error("use_mask is True but no mask defined in config.txt.")
            raise ValueError("use_mask is True but no mask defined in config.txt.")

        if use_fixed_vel_mask:
            mask_unit = input_mask["mask_unit"].iloc[0]
            mask_start = float(input_mask["mask_start"].iloc[0]) * au.Unit(mask_unit)
            mask_end = float(input_mask["mask_end"].iloc[0]) * au.Unit(mask_unit)
            mask = fixed_velocity_mask_ppv(
                cube_data[ref_line].shape, ov_hdr, mask_start, mask_end, mask_unit
            )
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
            n_mask = 0  # "first": only the reference line

        for n_mask_i in range(1, n_mask + 1):
            line_i = line_names[n_mask_i].upper()
            if line_i not in cube_data:
                continue
            mask_i = construct_mask_ppv(cube_data[line_i].value, SN_processing)
            mask = mask.astype(int) | mask_i.astype(int)
            LOG.info(f"Combined PPV mask includes {line_i}.")

        if ref_line_method == "ref+HI":
            if "HI" in cube_data:
                mask_hi = construct_mask_ppv(cube_data["HI"].value, SN_processing)
                mask = mask.astype(int) | mask_hi.astype(int)
                LOG.info("ref+HI mask: combined reference line and HI masks.")
            else:
                LOG.warning("HI not found among loaded cubes; ignoring ref+HI option.")

        if strict_mask:
            LOG.info(
                "Applying strict spatial mask filter (connected-component, PPV grid)."
            )
            mask = apply_strict_mask_ppv(mask.astype(int))

    # ------------------------------------------------------------------
    # Apply edge trimming: zero out the half-beam border of the footprint
    # across all channels. This removes convolution-edge artefacts from
    # both the moments and the saved mask FITS file.
    # edge_mask is (ny, nx), broadcast across all channels via None axis.
    # ------------------------------------------------------------------
    mask = (np.asarray(mask) * edge_mask[None, :, :]).astype(int)

    # Combined output NaN mask: True wherever the overlay has no data OR
    # the pixel lies in the eroded edge strip. Applied consistently to all
    # output files (moment maps, PPV mask cubes, and saved convolved cubes)
    # so every pipeline output shares the same NaN pattern.
    out_nan_mask = ~(ov_footprint & edge_mask.astype(bool))

    # Apply out_nan_mask to every in-memory cube BEFORE computing moments.
    # This is critical: without it, the edge-strip pixels (which have biased
    # values due to partial-beam convolution at the footprint boundary) would
    # enter the RMS and moment calculations even though they are masked out
    # in the integration mask. Setting them to NaN here ensures the moment
    # estimators see a clean cube with no partial-beam contamination.
    for name in cube_data:
        cube_val = cube_data[name].value.copy()
        cube_val[:, out_nan_mask] = np.nan
        cube_data[name] = cube_val * cube_data[name].unit

    if save_mask:
        save_ppv_mask_to_fits(
            mask, ov_hdr, source, "mask", folder, out_nan_mask=out_nan_mask
        )
        LOG.info(
            f"PPV mask cube written to: {os.path.join(folder, f'{source}_mask.fits')}"
        )

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
                    hfs_mask_name = f"mask_{line_name.lower()}"
                    save_ppv_mask_to_fits(
                        mask_hfs,
                        ov_hdr,
                        source,
                        hfs_mask_name,
                        folder,
                        out_nan_mask=out_nan_mask,
                    )
                    LOG.info(
                        f"PPV mask cube for {line_name} written to: {os.path.join(folder, f'{source}_{hfs_mask_name}.fits')}"
                    )

        mom_maps = get_mom_maps_ppv(
            cube_data[line_name.upper()],
            active_mask,
            vaxis,
            mom_calc,
            noise_mask=(
                build_noise_mask(
                    _noise_mask_df, vaxis, cube_data[line_name.upper()].shape
                )
                if _noise_mask_df is not None
                else None
            ),
        )

        ov_hdr_2d = twod_head(ov_hdr)
        line_desc = str(row["line_desc"])
        line_unit = str(row["line_unit"])

        quantities = {
            "mom0": (
                mom_maps["mom0"],
                "K km/s" if line_unit == "K" else f"{line_unit} km/s",
            ),
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

        # out_nan_mask was computed once above, from ov_footprint & edge_mask.
        _raw_hdr = cube_hdrs.get(line_name.upper(), {})
        _restfrq = (
            (_raw_hdr.get("RESTFRQ") or _raw_hdr.get("RESTFREQ"))
            if _raw_hdr else None
        )

        for quantity, (arr, bunit) in quantities.items():
            hdr_out = make_clean_header(
                ov_hdr_2d,
                is_cube=False,
                bunit=bunit,
                btype=quantity,
                line_name=line_name,
                line_desc=line_desc,
                object_name=source,
                bmaj_as=target_res_as,
                restfrq=None,   # 2D moment maps have no spectral axis
                meta=meta,
            )
            res_suffix = meta.get("res_suffix", "27p0as")
            fname_fits = os.path.join(
                folder, f"{source}_{line_name}_{res_suffix}_{quantity}.fits"
            )
            data_out = (
                arr.value.copy()
                if hasattr(arr, "value")
                else np.asarray(arr, dtype=float).copy()
            )
            data_out[out_nan_mask] = np.nan
            fits.writeto(fname_fits, data=data_out, header=hdr_out, overwrite=True)

        LOG.info(f"Compute moment maps and write to file for line: {line_name}.")


def run_fits(
    source,
    fname,
    meta,
    maps,
    cubes,
    params,
    input_mask=None,
    hfs_data=None,
    noise_mask_df=None,
):
    """
    Write FITS moment maps, 2D map images, and mask cube(s) for *source*.

    This is the entry point for the "fits" pipeline stage.

    Moment maps (if save_mom_maps is True) are computed PPV-native: directly
    on the convolved, overlay-aligned PPV cubes, never via the hex-grid
    .ecsv table — see run_moments_ppv and the module docstring. This
    requires either a save_cubes cube on disk from a prior fits stage run, or the raw
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
    save_mom_maps = meta.get("save_mom_maps", True)
    save_maps = meta.get("save_maps", True)
    save_mask = meta.get("save_mask", True)
    save_cubes = meta.get("save_cubes", True)
    folder = meta.get("folder_savefits", "./saved_fits_files/")

    if not (save_mom_maps or save_maps or save_mask or save_cubes):
        LOG.info(f"Output writing disabled for {source}; skipping.")
        return

    os.makedirs(folder, exist_ok=True)

    # ------------------------------------------------------------------
    # Load overlay cube to get the reference WCS and footprint mask
    # ------------------------------------------------------------------
    data_dir = meta.get("data_dir", "data/")
    overlay_file = meta.get("overlay_file", "")
    from os import path as _path

    overlay_fname = (
        _path.join(data_dir, overlay_file)
        if source in overlay_file
        else _path.join(data_dir, source + overlay_file)
    )

    LOG.info(f"Overlay file: {overlay_fname}")
    ov_cube, ov_hdr = fits.getdata(overlay_fname, header=True)

    # Resolve per-source target resolution (target_res, target_res_pc,
    # res_suffix) from the overlay header and source params.  This handles
    # all three modes (angular / physical / native) correctly without
    # requiring the regrid stage to have run first.
    resolve_meta_resolution(source, params, meta, ov_hdr=ov_hdr, log=LOG)
    target_res_as = meta["target_res"]

    # Build the overlay footprint: True wherever at least one spectral channel
    # is finite. Used as the authoritative NaN/valid mask for all output files.
    ov_footprint = np.any(np.isfinite(ov_cube), axis=0)  # (ny, nx) bool

    # ------------------------------------------------------------------
    # Moment maps — PPV-native, NOT from the hex-grid .ecsv table.
    # The PPV mask(s) used here are also written out (if save_mask is True)
    # as a byproduct of this same call, since the mask only exists as a
    # plain array inside run_moments_ppv.
    # ------------------------------------------------------------------
    if save_mom_maps:
        run_moments_ppv(
            source,
            meta,
            cubes,
            input_mask,
            hfs_data,
            params,
            folder,
            save_mask=save_mask,
            noise_mask_df=noise_mask_df,
        )
        LOG.info(f"Moment map FITS files written to: {folder}")
    elif save_mask:
        LOG.warning(
            f"save_mask is True but save_mom_maps is False for {source}; "
            "the PPV mask is only built while computing moments, so no "
            "mask FITS file(s) will be written. Set save_mom_maps = true "
            "to enable mask output."
        )

    # ------------------------------------------------------------------
    # 2D map images — PPV-native, analogous to the cube/moment path.
    #
    # Instead of regridding hex-grid MAP_*/EMAP_* columns via nearest-
    # neighbour interpolation (save_to_fits), we now read the raw input
    # map FITS files directly, convolve them to the target resolution on
    # their native pixel grid, and reproject onto the overlay WCS with
    # bilinear interpolation — exactly the same pipeline as for spectral
    # cubes (get_convolved_ppv_cube). This eliminates the blocky nearest-
    # neighbour artefacts that arose from interpolating sparse hex-grid
    # points back onto the dense overlay grid.
    # ------------------------------------------------------------------
    if save_maps:
        _edge = build_edge_mask(
            ov_footprint,
            ov_hdr,
            target_res_as,
            fov_erosion_beams=meta.get("fov_erosion_beams", 0.5),
        )
        _out_nan = ~(ov_footprint & _edge.astype(bool))
        res_suffix = meta.get("res_suffix", "27p0as")

        for _, map_row in maps.iterrows():
            map_name = str(map_row["map_name"])
            map_dir = str(map_row["map_dir"])
            map_ext = str(map_row["map_ext"])
            map_uc = str(map_row.get("map_uc", "")).strip()

            # --- primary map ---
            try:
                map_data, _ = get_convolved_map(
                    source, map_name, map_dir, map_ext,
                    target_res_as, ov_hdr, log=LOG,
                )
                map_data = np.where(_out_nan, np.nan, map_data)
                map_hdr = make_clean_header(
                    ov_hdr,
                    is_cube=False,
                    bunit=str(map_row.get("map_unit", "")),
                    btype=map_name,
                    object_name=source,
                    bmaj_as=target_res_as,
                    meta=meta,
                )
                fname_map = os.path.join(folder, f"{source}_{map_name}_{res_suffix}.fits")
                fits.writeto(fname_map, data=map_data, header=map_hdr, overwrite=True)
            except FileNotFoundError:
                LOG.warning(f"Skipping map {map_name}: raw input file not found.")
                continue

            # --- uncertainty map (optional; only if map_uc is defined) ---
            if map_uc and map_uc not in ("nan", ""):
                try:
                    unc_data, _ = get_convolved_map(
                        source, map_name, map_dir, map_uc,
                        target_res_as, ov_hdr, log=LOG,
                    )
                    unc_data = np.where(_out_nan, np.nan, unc_data)
                    unc_hdr = make_clean_header(
                        ov_hdr,
                        is_cube=False,
                        bunit=str(map_row.get("map_unit", "")),
                        btype=f"{map_name}_err",
                        object_name=source,
                        bmaj_as=target_res_as,
                        meta=meta,
                    )
                    fname_unc = os.path.join(
                        folder, f"{source}_{map_name}_{res_suffix}_err.fits"
                    )
                    fits.writeto(fname_unc, data=unc_data, header=unc_hdr, overwrite=True)
                except FileNotFoundError:
                    LOG.warning(
                        f"Skipping uncertainty map for {map_name}: "
                        f"raw input file {map_uc} not found."
                    )

        LOG.info(f"2D map FITS files written to: {folder}")

    # ------------------------------------------------------------------
    # Convolved cubes — independent of save_mom_maps so users can save
    # cubes without computing moment maps.
    # Cubes are convolved from raw input, reprojected onto the overlay WCS,
    # footprint-masked (ov_footprint) and edge-eroded (fov_erosion_beams)
    # exactly like the moment maps.
    # ------------------------------------------------------------------
    if save_cubes:
        _edge = build_edge_mask(
            ov_footprint,
            ov_hdr,
            target_res_as,
            fov_erosion_beams=meta.get("fov_erosion_beams", 0.5),
        )
        _out_nan = ~(ov_footprint & _edge.astype(bool))
        res_suffix = meta.get("res_suffix", "27p0as")

        for _, row in cubes.iterrows():
            try:
                # Read raw header before convolution to preserve the
                # per-cube rest frequency (RESTFRQ/RESTFREQ), which is
                # overwritten by the overlay header during reprojection.
                raw_cube_path = os.path.join(
                    str(row["line_dir"]), source + str(row["line_ext"])
                )
                _raw_hdr = fits.getheader(raw_cube_path) if os.path.exists(raw_cube_path) else {}
                _restfrq = _raw_hdr.get("RESTFRQ") or _raw_hdr.get("RESTFREQ")

                cube_data, cube_hdr = get_convolved_ppv_cube(
                    source,
                    str(row["line_name"]),
                    str(row["line_dir"]),
                    str(row["line_ext"]),
                    meta,
                    ov_hdr,
                    log=LOG,
                )
            except FileNotFoundError:
                LOG.warning(
                    f"Skipping cube {row['line_name']}: raw input file not found."
                )
                continue

            cube_data[:, _out_nan] = np.nan

            cube_fits_path = os.path.join(
                folder, f"{source}_{row['line_name']}_{res_suffix}.fits"
            )
            fits.writeto(
                cube_fits_path,
                data=cube_data,
                header=make_clean_header(
                    ov_hdr,
                    is_cube=True,
                    bunit=str(row.get("line_unit", "")),
                    line_name=str(row["line_name"]),
                    line_desc=str(row.get("line_desc", "")),
                    object_name=source,
                    bmaj_as=target_res_as,
                    restfrq=_restfrq,
                    meta=meta,
                ),
                overwrite=True,
            )
            LOG.info(f"Convolved cube written to: {cube_fits_path}")


def _resolve_target_res(params, meta, ov_hdr=None):
    """
    Return the target resolution in arcseconds from *meta*.

    After ``stage_regrid.run_sampling`` runs, ``meta["target_res"]`` always
    holds the resolved arcsecond value regardless of the original
    ``resolution`` mode (angular / physical / native).

    The corresponding parsec value is in ``meta["target_res_pc"]``.

    The *params* and *ov_hdr* arguments are kept for call-site compatibility
    but are no longer used here; all resolution conversion happens once in
    ``run_sampling`` and the arcsecond result is stored in ``meta["target_res"]``.
    """
    return float(meta.get("target_res", 27.0))
