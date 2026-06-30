"""
utils_fits.py — FITS and astronomy utility functions.

All functions are self-contained (no imports from the legacy scripts directory).
They are used by multiple pipeline stages and are designed to be importable
and usable independently of the pipeline infrastructure.

Contents
--------
Basic FITS I/O
    get_beam_arcsec  — extract beam size from a FITS header
    read_fits_cube   — load a cube, squeezing any degenerate 4th axis

Header utilities
    twod_head        — reduce a 3D/4D FITS header to 2D

Sampling grid
    hex_grid              — generate a hexagonal RA/Dec grid
    make_sampling_points  — generate hex grid points clipped to a mask

Deprojection
    deproject        — compute galactocentric radii and polar angles

Gaussian PSF
    gaussian_PSF_2D  — create a 2-D rotated Gaussian kernel

Beam deconvolution
    deconvolve_gauss — deconvolve one Gaussian from another (MIRIAD port)

Spatial convolution
    conv_with_gauss  — convolve a cube or map to a target Gaussian beam
"""

import copy
import warnings
import numpy as np
from pathlib import Path
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel
from astropy.utils.console import ProgressBar

from hexmaps.logger import get_logger

# Utility functions do not have their own pipeline "stage" — log messages
# should appear under whichever stage is calling them. Each function below
# accepts an optional `log` parameter (a StageLogger from get_logger); if not
# provided, falls back to the "Loading" stage (the typical context for
# standalone/analysis use of these functions).
_DEFAULT_LOG = get_logger("Loading")

warnings.filterwarnings("ignore")


# ============================================================================
# Basic FITS I/O
# ============================================================================


def get_beam_arcsec(fits_path: str, log=None) -> u.Quantity:
    """
    Return the beam major axis (BMAJ) in arcseconds from a FITS header.

    Parameters
    ----------
    fits_path : str or Path
    log : StageLogger, optional
        Logger to use for error messages (from get_logger()). Defaults to
        the "Loading" stage if not provided.

    Returns
    -------
    beam_as : astropy.units.Quantity in arcseconds

    Raises
    ------
    FileNotFoundError : if the file does not exist
    KeyError          : if BMAJ is absent from the header
    """
    log = log or _DEFAULT_LOG
    fits_path = Path(fits_path)
    if not fits_path.exists():
        log.error(f"FITS file not found: {fits_path}")
        raise FileNotFoundError(f"FITS file not found: {fits_path}")
    hdr = fits.getheader(fits_path)
    if "BMAJ" not in hdr:
        log.error(f"BMAJ not found in header of {fits_path}")
        raise KeyError(f"BMAJ not found in header of {fits_path}")
    return (hdr["BMAJ"] * u.deg).to(u.arcsec)


def read_fits_cube(fits_path: str, log=None):
    """
    Read a FITS file, squeezing any degenerate Stokes (4th) axis.

    Many radio FITS cubes have a 4th Stokes axis of length 1.  This function
    removes it so that the returned array is always 3-D (channels × y × x).

    Parameters
    ----------
    fits_path : str or Path
    log : StageLogger, optional
        Logger to use for error messages (from get_logger()). Defaults to
        the "Loading" stage if not provided.

    Returns
    -------
    data : np.ndarray (3D)
    hdr  : astropy.io.fits.Header (updated to reflect 3D)

    Raises
    ------
    FileNotFoundError : if the file does not exist
    """
    log = log or _DEFAULT_LOG
    fits_path = Path(fits_path)
    if not fits_path.exists():
        log.error(f"FITS file not found: {fits_path}")
        raise FileNotFoundError(f"FITS file not found: {fits_path}")
    data, hdr = fits.getdata(fits_path, header=True)
    if hdr["NAXIS"] == 4:
        data = np.squeeze(data, axis=0)
        hdr["NAXIS"] = 3
        for key in ["NAXIS4", "CTYPE4", "CRVAL4", "CDELT4", "CRPIX4", "CUNIT4"]:
            hdr.remove(key, ignore_missing=True)
    return data, hdr


# ============================================================================
# Header utilities
# ============================================================================


def twod_head(hdul_header):
    """
    Reduce a FITS header to 2-D by removing all axes beyond the second.

    This is used to create a 2-D WCS header from a 3-D cube header so that
    astropy.wcs.WCS can be used for spatial operations without the spectral axis.

    Port of IDL twod_head (A. Leroy, 2008) by J. den Brok (2019).

    Parameters
    ----------
    hdul_header : astropy.io.fits.Header — header with NAXIS ≥ 2

    Returns
    -------
    header_copy : Header — new header with NAXIS=2 and no axis-3+ keywords
    """
    header_copy = copy.copy(hdul_header)
    naxis = hdul_header["NAXIS"]
    header_copy["NAXIS"] = 2
    if "WCSAXES" in header_copy:
        header_copy["WCSAXES"] = 2
    if naxis > 2:
        header_copy["WCSAXES"] = 2
        for i in range(3, naxis + 1):
            del header_copy["*{}*".format(int(i))]
    return header_copy


# ============================================================================
# Sampling grid
# ============================================================================


def hex_grid(ctr_x, ctr_y, spacing, radec=False, r_limit=None, e_limit=None):
    """
    Generate a hexagonal close-packed grid centred on (ctr_x, ctr_y).

    The grid is constructed in Cartesian coordinates and then optionally
    corrected for the cos(Dec) foreshortening in RA.

    Parameters
    ----------
    ctr_x, ctr_y : float  — grid centre coordinates
    spacing      : float  — separation between adjacent grid points
                            (same units as r_limit or e_limit)
    radec        : bool   — if True, divide x offsets by cos(Dec) to correct
                            for RA foreshortening; set True when working in
                            RA/Dec degrees
    r_limit      : float  — keep only points within this circular radius
    e_limit      : float  — keep only points within a square of this half-extent
                            (one of r_limit or e_limit must be provided)

    Returns
    -------
    xout, yout : np.ndarrays — grid point coordinates
                 Returns (np.nan, np.nan) if no points survive the clipping.

    Raises
    ------
    TypeError if neither r_limit nor e_limit is provided.
    """
    x_spacing = spacing
    y_spacing = spacing * np.sin(np.deg2rad(60))  # row offset for hex packing

    if e_limit is None and r_limit is not None:
        scale = r_limit
    elif r_limit is None and e_limit is not None:
        scale = e_limit / 2
    else:
        raise TypeError("Provide exactly one of r_limit or e_limit to hex_grid.")

    half_ny = np.ceil(scale / y_spacing)
    half_nx = np.ceil(scale / x_spacing) + 1

    # Build 2-D coordinate arrays
    x = np.outer(np.ones(2 * int(half_ny) + 1), np.arange(2 * int(half_nx) + 1))
    y = np.outer(np.arange(2 * int(half_ny) + 1), np.ones(2 * int(half_nx) + 1))
    x -= half_nx
    y -= half_ny
    x *= x_spacing
    # Offset every other row by half a spacing to achieve hex packing
    x += 0.5 * x_spacing * (np.dot(abs(y) % 2 == 1, 1))
    y *= y_spacing

    r = np.sqrt(x**2 + y**2)
    keep = (
        np.where(r < r_limit)
        if r_limit is not None
        else np.where(np.logical_and(abs(x) < e_limit / 2, abs(y) < e_limit / 2))
    )

    if len(keep[0]) == 0:
        return np.nan, np.nan

    yout = y[keep] + ctr_y
    xout = (x[keep] / np.cos(np.deg2rad(yout)) + ctr_x) if radec else (x[keep] + ctr_x)
    return xout, yout


def make_sampling_points(
    ra_ctr,
    dec_ctr,
    max_rad,
    spacing,
    mask,
    hdr_mask,
    overlay_in=None,
    overlay_hdr_in=None,
    show=False,
    log=None,
):
    """
    Generate hexagonal sampling points clipped to a binary sky mask.

    Steps
    -----
    1. If *mask* is 3-D, collapse it along axis 0 to get a 2-D footprint.
    2. If max_rad is "auto", compute the half-diagonal of the mask array as
       the maximum radius.
    3. Generate a hex grid using hex_grid.
    4. Convert the grid RA/Dec to pixel coordinates using the mask WCS.
    5. Remove points outside the array boundary.
    6. Remove points where the mask is False (zero).

    Port of sampling.py (J. den Brok, 2019).

    Parameters
    ----------
    ra_ctr, dec_ctr : float — grid centre (degrees)
    max_rad         : float | "auto" — maximum radius in degrees
    spacing         : float — hex grid spacing in degrees
    mask            : np.ndarray (2D or 3D) — binary footprint mask
    hdr_mask        : FITS Header — WCS for the mask
    overlay_in      : str or array, optional — overlay for visualisation
    overlay_hdr_in  : FITS Header, optional — header for overlay
    show            : bool — if True, display sampling points on the overlay
    log : StageLogger, optional
        Logger to use for progress/error messages (from get_logger()).
        Defaults to the "Loading" stage if not provided.

    Returns
    -------
    samp_ra, samp_dec : np.ndarrays — coordinates of the surviving grid points.
                        Returns (np.nan, np.nan) if no points survive.
    """
    log = log or _DEFAULT_LOG

    # Collapse 3-D mask to 2-D
    if len(np.shape(mask)) == 3:
        log.info(f"Collapsing 3D mask to 2D footprint.")
        mask = np.sum(np.isfinite(mask), axis=0) >= 1
        hdr_mask = twod_head(hdr_mask)

    mask_dim = np.shape(mask)
    wcs = WCS(hdr_mask)

    # Auto-determine maximum radius from the mask array diagonal
    if max_rad == "auto":
        from astropy.coordinates import SkyCoord

        c1 = SkyCoord.from_pixel(0, 0, wcs)
        c2 = SkyCoord.from_pixel(mask_dim[1], mask_dim[0], wcs)
        max_rad = c1.separation(c2).value / 2
        log.info(f"Auto max_rad = {np.round(max_rad, 3)} deg.")

    samp_ra, samp_dec = hex_grid(ra_ctr, dec_ctr, spacing, radec=True, r_limit=max_rad)

    # Convert to pixel coordinates
    try:
        pixel_coords = wcs.all_world2pix(np.column_stack((samp_ra, samp_dec)), 0)
    except Exception:
        pixel_coords = wcs.all_world2pix(
            np.column_stack((samp_ra, samp_dec, np.zeros(len(samp_ra)))), 0
        )

    samp_x = np.array(np.rint(pixel_coords[:, 0]), dtype=int)
    samp_y = np.array(np.rint(pixel_coords[:, 1]), dtype=int)

    # Keep only points inside the array boundary
    keep = np.where(
        (samp_x >= 0) & (samp_y >= 0) & (samp_x < mask_dim[1]) & (samp_y < mask_dim[0])
    )[0]
    if len(keep) == 0:
        log.error(f"No sampling points inside mask bounds.")
        return np.nan, np.nan

    samp_ra, samp_dec = samp_ra[keep], samp_dec[keep]
    samp_x, samp_y = samp_x[keep], samp_y[keep]

    # Keep only points where the mask is True
    keep = np.where(mask[samp_y, samp_x])[0]
    if len(keep) == 0:
        log.error(f"No sampling points survive mask clipping.")
        return np.nan, np.nan

    samp_ra, samp_dec = samp_ra[keep], samp_dec[keep]

    if show:
        _show_sampling_points(
            samp_ra, samp_dec, mask, hdr_mask, overlay_in, overlay_hdr_in
        )

    return samp_ra, samp_dec


def _show_sampling_points(
    samp_ra, samp_dec, mask, hdr_mask, overlay_in, overlay_hdr_in
):
    """Visualise sampling points overlaid on the mask or overlay image."""
    import matplotlib.pyplot as plt

    if overlay_in is not None:
        if isinstance(overlay_in, str):
            overlay, overlay_hdr = fits.getdata(overlay_in, header=True)
        else:
            overlay = copy.deepcopy(overlay_in)
            overlay_hdr = overlay_hdr_in if overlay_hdr_in is not None else hdr_mask
        if len(np.shape(overlay)) == 3:
            overlay = np.nansum(overlay, 0)
            overlay_hdr = twod_head(overlay_hdr)
    else:
        overlay, overlay_hdr = copy.deepcopy(mask), hdr_mask

    wcs_ov = WCS(overlay_hdr)
    px = wcs_ov.all_world2pix(np.column_stack((samp_ra, samp_dec)), 0)
    plt.figure()
    plt.plot(px[:, 0], px[:, 1], "h", markersize=16)
    plt.show()


# ============================================================================
# Deprojection
# ============================================================================


def deproject(ra, dec, galpos, vector=False, gal=None):
    """
    Compute deprojected galactocentric radii and polar angles.

    Applies a rotation by the position angle and a stretching by 1/cos(incl)
    to convert observed (RA, Dec) offsets into intrinsic disk-plane coordinates.

    Port of IDL deproject (A. Leroy, 2001) by J. den Brok (2019).

    Parameters
    ----------
    ra, dec : array-like — observed coordinates (degrees, J2000)
    galpos  : list       — galaxy geometry:
                [pa_deg, inc_deg, ra_ctr, dec_ctr]  (4 elements, standard)
                [vlsr, pa_deg, inc_deg, ra_ctr, dec_ctr]  (5 elements)
    vector  : bool       — if True, ra/dec are already paired vectors of the
                           same length; if False, a 2-D grid is computed
    gal     : dict       — alternative to galpos; keys: posang_deg, incl_def,
                           RA, DEC

    Returns
    -------
    rgrid : np.ndarray — deprojected radius in degrees (same shape as ra/dec)
    tgrid : np.ndarray — deprojected polar angle in radians
    """
    np.seterr(divide="ignore", invalid="ignore")

    if gal is not None:
        pa = np.deg2rad(gal["posang_deg"])
        inc = np.deg2rad(gal["incl_def"])
        xctr = gal["RA"]
        yctr = gal["DEC"]
    elif len(galpos) == 5:
        pa = np.deg2rad(galpos[1])
        inc = np.deg2rad(galpos[2])
        xctr = galpos[3]
        yctr = galpos[4]
    else:
        pa = np.deg2rad(galpos[0])
        inc = np.deg2rad(galpos[1])
        xctr = galpos[2]
        yctr = galpos[3]

    ra_size = np.shape(ra)
    if ra_size[0] == 1 and not vector:
        rimg = np.outer(ra, np.ones(len(dec)))
        dimg = np.outer(np.ones(len(ra)), dec)
    else:
        rimg, dimg = ra, dec

    # Offset from centre, correcting RA for cos(Dec) projection
    xgrid = (rimg - xctr) * np.cos(np.deg2rad(yctr))
    ygrid = dimg - yctr

    # Rotate by (PA - 90°) to align with the major axis
    rotang = -(pa - np.pi / 2.0)
    deproj_x = xgrid * np.cos(rotang) + ygrid * np.sin(rotang)
    deproj_y = ygrid * np.cos(rotang) - xgrid * np.sin(rotang)

    # Stretch along the minor axis to correct for inclination
    deproj_y = deproj_y / np.cos(inc)

    rgrid = np.sqrt(deproj_x**2 + deproj_y**2)
    tgrid = np.arctan2(deproj_y, deproj_x)
    return rgrid, tgrid


# ============================================================================
# Gaussian PSF
# ============================================================================


def gaussian_PSF_2D(npix, a, center=False, normalize=False, log=None):
    """
    Create a 2-D rotated Gaussian PSF kernel array.

    Parameters
    ----------
    npix      : int or (int, int) — array size in pixels (square or rectangular)
    a         : list [offset, peak, fwhm_x, fwhm_y, cen_x, cen_y, rot_rad]
        offset  — additive baseline
        peak    — peak amplitude
        fwhm_x  — FWHM along the (rotated) x axis, in pixels
        fwhm_y  — FWHM along the (rotated) y axis, in pixels
        cen_x   — x centre pixel (ignored if center=True)
        cen_y   — y centre pixel (ignored if center=True)
        rot_rad — rotation angle in radians (CCW from x axis)
    center    : bool — if True, place the PSF at the centre of the array
    normalize : bool — if True, normalise so the kernel sums to 1
    log : StageLogger, optional
        Logger to use for error messages (from get_logger()). Defaults to
        the "Loading" stage if not provided.

    Returns
    -------
    output : np.ndarray (ny, nx) — the Gaussian kernel
    """
    log = log or _DEFAULT_LOG

    if isinstance(npix, (int, float)):
        nx = ny = int(npix)
    elif hasattr(npix, "__len__") and len(npix) == 2:
        nx, ny = int(npix[0]), int(npix[1])
    else:
        log.error(f"Invalid npix: {npix}")
        return None

    xarr = np.tile(np.arange(nx), ny).reshape(ny, nx).astype(float)
    yarr = np.repeat(np.arange(ny), nx).reshape(ny, nx).astype(float)

    cenx = (nx - 1) / 2 if center else a[4]
    ceny = (ny - 1) / 2 if center else a[5]

    fac = 2 * np.sqrt(2 * np.log(2))  # FWHM → sigma conversion
    ang = a[6]
    widthx = a[2] / fac
    widthy = a[3] / fac
    s, c = np.sin(ang), np.cos(ang)

    xarr -= cenx
    yarr -= ceny
    t = xarr * (c / widthx) + yarr * (s / widthx)
    yarr = xarr * (s / widthy) - yarr * (c / widthy)
    xarr = t

    output = a[0] + a[1] * np.exp(-0.5 * (xarr**2 + yarr**2))
    if normalize:
        output /= np.sum(output)
    return output


# ============================================================================
# Beam deconvolution
# ============================================================================


def deconvolve_gauss(
    meas_maj, beam_maj, meas_min=None, meas_pa=None, beam_min=None, beam_pa=None
):
    """
    Deconvolve a Gaussian beam from a measured Gaussian source size.

    Finds the intrinsic source size by subtracting the beam in quadrature
    (in 2-D, including position angle rotation).  Port of MIRIAD gaupar.for.

    This is used in conv_with_gauss to compute the convolution kernel needed
    to bring a native beam up to the target resolution.

    Parameters
    ----------
    meas_maj : float — measured major axis FWHM (arcsec)
    beam_maj : float — beam major axis FWHM (arcsec)
    meas_min : float, optional — measured minor axis FWHM (default: meas_maj)
    meas_pa  : float, optional — measured position angle (degrees, default: 0)
    beam_min : float, optional — beam minor axis FWHM (default: beam_maj)
    beam_pa  : float, optional — beam position angle (degrees, default: 0)

    Returns
    -------
    src_maj, src_min, src_pa : float — intrinsic source Gaussian parameters
    info : [worked, point_source]
        worked       — True if deconvolution succeeded
        point_source — True if the source is unresolved (within tolerance)
    """
    if beam_min is None:
        meas_min = meas_maj
    if meas_pa is None:
        meas_pa = 0.0
    if beam_pa is None:
        beam_pa = 0.0

    mt = np.deg2rad(meas_pa)
    bt = np.deg2rad(beam_pa)

    alpha = (
        (meas_maj * np.cos(mt)) ** 2
        + (meas_min * np.sin(mt)) ** 2
        - (beam_maj * np.cos(bt)) ** 2
        - (beam_min * np.sin(bt)) ** 2
    )
    beta = (
        (meas_maj * np.sin(mt)) ** 2
        + (meas_min * np.cos(mt)) ** 2
        - (beam_maj * np.sin(bt)) ** 2
        - (beam_min * np.cos(bt)) ** 2
    )
    gamma = 2 * (
        (meas_min**2 - meas_maj**2) * np.sin(mt) * np.cos(mt)
        - (beam_min**2 - beam_maj**2) * np.sin(bt) * np.cos(bt)
    )

    s = alpha + beta
    t = np.sqrt((alpha - beta) ** 2 + gamma**2)
    limit = (
        0.1 * min(meas_min or meas_maj, meas_maj, beam_maj, beam_min or beam_maj) ** 2
    )

    if alpha < 0 or beta < 0 or s < t:
        worked = False
        point = (0.5 * (s - t) < limit) and (alpha > -limit) and (beta > -limit)
        return 0.0, 0.0, 0.0, [worked, point]

    src_maj = np.sqrt(0.5 * (s + t))
    src_min = np.sqrt(0.5 * (s - t))
    src_pa = (
        0.0
        if (abs(gamma) + abs(alpha - beta)) == 0
        else np.rad2deg(0.5 * np.arctan(-gamma / (alpha - beta)))
    )
    return src_maj, src_min, src_pa, [True, False]


# ============================================================================
# Spatial convolution
# ============================================================================


def _round_sig(x, sig=2):
    """Round *x* to *sig* significant figures."""
    return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)


def _get_pixel_scale(hdr, tol=0.1, log=None):
    """
    Return the pixel scale in degrees from a FITS header.

    Issues a warning if the x and y pixel scales differ by more than *tol*
    arcseconds and returns the geometric mean in that case.

    Parameters
    ----------
    log : StageLogger, optional
        Logger to use for the warning (from get_logger()). Defaults to
        the "Loading" stage if not provided.
    """
    log = log or _DEFAULT_LOG
    w = WCS(hdr)
    scales = proj_plane_pixel_scales(w)
    px_dx = scales[0] * u.deg
    px_dy = scales[1] * u.deg
    if abs(px_dx - px_dy) > tol * u.arcsec:
        log.warning(
            f"Pixel scale differs in X and Y: "
            f"{px_dx.to(u.arcsec):.3f} vs {px_dy.to(u.arcsec):.3f}. "
            "Using geometric mean."
        )
        return np.sqrt(px_dx * px_dy).value
    return px_dx.value


def _convolve_func(data, kernel, method="fft"):
    """Dispatch to the appropriate astropy convolution function."""
    if method == "direct":
        return convolve(data, kernel, allow_huge=True)
    return convolve_fft(data, kernel, allow_huge=True)


def conv_with_gauss(
    in_data,
    in_hdr=None,
    start_beam=None,
    pix_deg=None,
    target_beam=None,
    no_ft=False,
    in_weight=None,
    out_weight_file=None,
    out_file=None,
    unc=False,
    perbeam=False,
    quiet=False,
    log=None,
):
    """
    Convolve a 2-D map or 3-D cube to a target Gaussian beam.

    Port of IDL conv_with_gauss (A. Leroy / J. den Brok, 2020).

    The convolution kernel is computed by deconvolving the current beam from
    the target beam.  The kernel is then applied plane-by-plane for cubes, or
    directly for 2-D maps.

    Unit corrections
    ----------------
    unc=True
        Treat *in_data* as an uncertainty map.  The data is squared before
        convolution (so that uncertainties add in quadrature) and the square
        root is taken after.  The beam-area correction is also applied.
    perbeam=True
        Correct for the change in beam solid angle when the data is in
        surface-brightness units per beam (e.g. Jy/beam or K).

    Parameters
    ----------
    in_data      : np.ndarray or str — input data array or FITS path
    in_hdr       : FITS Header       — required if in_data is an array
    start_beam   : float or list     — override the input beam size (arcsec)
    pix_deg      : float             — override the pixel scale (degrees)
    target_beam  : float or list     — target beam FWHM in arcseconds;
                                       list [maj, min, pa] for elliptical beams
    no_ft        : bool              — use direct convolution instead of FFT
    out_file     : str               — write the convolved data to this FITS path
    unc          : bool              — treat as uncertainty map (see above)
    perbeam      : bool              — apply per-beam correction (see above)
    quiet        : bool              — suppress progress output
    log : StageLogger, optional
        Logger to use for progress/error messages (from get_logger()).
        Defaults to the "Loading" stage if not provided.

    Returns
    -------
    data : np.ndarray — convolved data
    hdr  : FITS Header — updated with new BMAJ/BMIN/BPA keywords
    Returns (None, None) if the deconvolution fails.
    """
    log = log or _DEFAULT_LOG

    # Load data
    if isinstance(in_data, str):
        data, hdr = fits.getdata(in_data, header=True)
    else:
        data = copy.deepcopy(in_data)
        hdr = in_hdr

    # Normalise target_beam to a 3-element array [maj, min, pa]
    if target_beam is not None:
        if isinstance(target_beam, (float, int)):
            target_beam = np.array([float(target_beam), float(target_beam), 0.0])
        elif hasattr(target_beam, "__len__"):
            target_beam = np.array(list(target_beam) + [0.0] * (3 - len(target_beam)))

    flux_before = np.nansum(data)

    # Pixel scale
    if pix_deg is None:
        pix_deg = _get_pixel_scale(hdr, tol=0.1 * hdr.get("BMIN", 1.0), log=log)
    as_per_pix = pix_deg * 3600.0

    if unc:
        data = data**2  # square uncertainties before convolution

    # Identify the starting beam
    if start_beam is not None:
        if isinstance(start_beam, float):
            current_beam = [start_beam, start_beam, 0.0]
        elif hasattr(start_beam, "__len__"):
            current_beam = list(start_beam) + [0.0] * (3 - len(start_beam))
        else:
            log.error(f"Unknown start_beam format.")
            return None, None
    else:
        bmaj = hdr["BMAJ"] * 3600
        bmin = hdr["BMIN"] * 3600
        bpa = hdr.get("BPA", 0.0)
        current_beam = [bmaj, bmin, bpa]

    # Compute the convolution kernel by deconvolving the current beam
    src_maj, src_min, src_pa, info = deconvolve_gauss(
        target_beam[0],
        current_beam[0],
        target_beam[1],
        target_beam[2],
        current_beam[1],
        current_beam[2],
    )
    if not info[0]:
        log.error(
            f"Cannot compute convolution kernel: "
            f"target beam {target_beam} is smaller than current beam {current_beam}."
        )
        return None, None
    if info[1]:
        log.warning(
            f"Target and starting beam are nearly identical; "
            "kernel will be very small."
        )

    # Build the kernel array (6 × FWHM in pixels, capped to the data size)
    kern_size = int(6.0 * np.rint(src_maj / as_per_pix) + 1)
    dim_data = np.shape(data)
    dim_x = dim_data[1] if len(dim_data) == 3 else dim_data[0]
    dim_y = dim_data[2] if len(dim_data) == 3 else dim_data[1]
    if kern_size > dim_x or kern_size > dim_y:
        kern_size = int(np.floor(min(dim_x, dim_y) / 2 - 2) * 2 + 1)

    kernel = gaussian_PSF_2D(
        kern_size,
        [
            0.0,
            1.0,
            src_maj / as_per_pix,
            src_min / as_per_pix,
            0.0,
            0.0,
            np.pi / 2.0 + np.deg2rad(src_pa),
        ],
        center=True,
        normalize=True,
        log=log,
    )

    method = "direct" if no_ft else "fft"

    # Convolve: plane by plane for cubes, direct for 2-D maps
    if len(dim_data) == 3:
        new_data = copy.deepcopy(data)
        if not quiet:
            log.info(f"Convolving cube ({dim_data[0]} planes):")
        for plane in ProgressBar(range(dim_data[0])):
            new_data[plane, :, :] = _convolve_func(data[plane, :, :], kernel, method)
        data = new_data
    else:
        data = _convolve_func(data, kernel, method)

    # Beam-area correction for per-beam or uncertainty maps
    if unc or perbeam:
        cur_fwhm = np.sqrt(current_beam[0] * current_beam[1])
        ppbeam_start = (cur_fwhm / as_per_pix / 2) ** 2 / np.log(2) * np.pi
        tgt_fwhm = np.sqrt(target_beam[0] * target_beam[1])
        ppbeam_final = (tgt_fwhm / as_per_pix / 2) ** 2 / np.log(2) * np.pi

    if unc:
        data = np.sqrt(data) * np.sqrt(ppbeam_start / ppbeam_final)
    if perbeam:
        data *= ppbeam_final / ppbeam_start

    if not quiet:
        log.info(
            f"Pixel scale: {_round_sig(as_per_pix, 3)} arcsec/px  |  "
            f"Input beam: {[round(b, 1) for b in current_beam]} arcsec  |  "
            f"Target beam: {[round(b, 1) for b in target_beam]} arcsec  |  "
            f"Flux ratio: {_round_sig(np.nansum(data) / flux_before)}"
        )

    # Update header with the new beam
    hdr["BMAJ"] = (target_beam[0] / 3600.0, "FWHM BEAM IN DEGREES")
    hdr["BMIN"] = (target_beam[1] / 3600.0, "FWHM BEAM IN DEGREES")
    hdr["BPA"] = (target_beam[2], "POSITION ANGLE IN DEGREES")
    hdr["HISTORY"] = (
        f"conv_with_gauss: convolved with "
        f"[{src_maj:.2f}, {src_min:.2f}, {src_pa:.2f}] arcsec kernel"
    )

    if out_file is not None:
        fits.writeto(out_file, data, hdr, overwrite=True)

    return data, hdr


# ============================================================================
# Cube reprojection
# ============================================================================

"""
Drop-in replacement for reproject.reproject_interp that transparently fixes
failure modes produced by GILDAS/CLASS and other radio-astronomy FITS cubes:

  1. GLS projection  -- ``RA---GLS`` / ``DEC--GLS`` are non-standard GILDAS
     aliases for the Sanson-Flamsteed projection (``SFL``).  astropy's WCS
     does not recognise GLS and silently misplaces every coordinate, giving
     all-NaN output.  Fixed by renaming the CTYPE values before any WCS call.

  2. Different rest frequencies -- when the source and target cubes belong to
     different molecular lines (e.g. HCN and CO), a full 3-D WCS reprojection
     maps every source channel to a spectral index far outside the target array
     because the two rest frequencies encode different absolute frequencies for
     the same radial velocity.  Fixed by stripping the spectral axis from both
     WCS objects and reprojecting each spatial plane individually.

  3. Spectral resampling -- after fix 2, the spatially reprojected cube still
     carries the *source* velocity grid.  The cube is resampled onto the target
     velocity grid using 1-D interpolation along the spectral axis so that the
     output shape and WCS exactly match the requested target header.

  4. RESTFRQ vs RESTFREQ key name -- the FITS standard 8-character keyword
     limit means the rest frequency may be stored as either ``RESTFREQ`` or
     ``RESTFRQ`` depending on the writing software.  Both spellings are now
     recognised everywhere rest frequency is read.

  5. VELO-LSR spectral axis type -- a deprecated FITS keyword that astropy
     internally rewrites to ``VRAD`` via its ``spcfix`` mechanism.  This
     rewrite can produce all-NaN output when the two cubes have different
     ``CTYPE3`` variants (``VELO-LSR`` vs ``VRAD``) or when the spcfix
     transformation is applied inconsistently between source and target.
     Fixed by explicitly normalising ``CTYPE3`` to ``VRAD`` before any
     reprojection call.

Usage
-----
Replace every call to ``reproject_interp`` with ``reproject_interp_gildas``.
The function signature is identical, with one optional extra keyword:

    from reproject_gildas import reproject_interp_gildas

    array_out, footprint = reproject_interp_gildas(src_hdu, target_header)

    # Choose spectral interpolation order (default 'linear'):
    array_out, footprint = reproject_interp_gildas(
        src_hdu, target_header, spectral_order='cubic')
"""

import warnings

import numpy as np
import scipy.interpolate
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from reproject import reproject_interp

# ─── Internal helpers ────────────────────────────────────────────────────────


def _extract_data_header(input_data, hdu_in):
    """Return (ndarray, Header) from any input type accepted by reproject_interp."""
    if isinstance(input_data, (str,)) or hasattr(input_data, "__fspath__"):
        with fits.open(input_data) as hdul:
            hdu = hdul[hdu_in]
            return hdu.data.copy(), hdu.header.copy()
    if isinstance(input_data, fits.HDUList):
        hdu = input_data[hdu_in]
        return hdu.data.copy(), hdu.header.copy()
    if isinstance(input_data, (fits.PrimaryHDU, fits.ImageHDU, fits.CompImageHDU)):
        return input_data.data.copy(), input_data.header.copy()
    if isinstance(input_data, tuple) and len(input_data) == 2:
        arr, meta = input_data
        if isinstance(meta, fits.Header):
            return np.asarray(arr), meta.copy()
        # WCS object — convert to a minimal header
        hdr = meta.to_header(relax=True)
        hdr["NAXIS"] = arr.ndim
        for i, n in enumerate(np.asarray(arr).shape[::-1], 1):
            hdr[f"NAXIS{i}"] = n
        return np.asarray(arr), hdr
    # NDData
    if hasattr(input_data, "data") and hasattr(input_data, "wcs"):
        hdr = input_data.wcs.to_header(relax=True)
        hdr["NAXIS"] = input_data.data.ndim
        for i, n in enumerate(input_data.data.shape[::-1], 1):
            hdr[f"NAXIS{i}"] = n
        return np.asarray(input_data.data), hdr
    raise TypeError(f"Unsupported input_data type: {type(input_data)}")


def _extract_output_header(output_projection, shape_out):
    """Return a Header from any output_projection type."""
    if isinstance(output_projection, fits.Header):
        return output_projection.copy()
    hdr = output_projection.to_header(relax=True)
    if shape_out is not None:
        hdr["NAXIS"] = len(shape_out)
        for i, n in enumerate(shape_out[::-1], 1):
            hdr[f"NAXIS{i}"] = n
    return hdr


def _fix_gls(hdr):
    """
    Fix GILDAS-specific header issues (operates on a copy).

    Changes applied:
    - RA---GLS  → RA---SFL  (and likewise for DEC)
    - CRPIX3 = 0 → CRPIX3 = 1, CRVAL3 shifted by one CDELT3
      (FITS pixel coordinates are 1-based; CLASS writes 0)
    """
    for key in ("CTYPE1", "CTYPE2"):
        val = str(hdr.get(key, ""))
        if "GLS" in val:
            hdr[key] = val.strip().replace("GLS", "SFL")
    if "CRPIX3" in hdr and float(hdr["CRPIX3"]) == 0.0:
        hdr["CRVAL3"] = float(hdr["CRVAL3"]) + float(hdr["CDELT3"])
        hdr["CRPIX3"] = 1.0
    return hdr


def _restfreq(hdr):
    """
    Return the rest frequency (Hz) as float, or None if absent.

    Checks both ``RESTFREQ`` (used by GILDAS/CLASS and some other packages)
    and ``RESTFRQ`` (the 8-character FITS standard keyword) so that neither
    spelling is silently missed.
    """
    if "RESTFREQ" in hdr:
        return float(hdr["RESTFREQ"])
    if "RESTFRQ" in hdr:
        return float(hdr["RESTFRQ"])
    return None


def _normalize_spectral_ctype(hdr):
    """
    Normalise deprecated or non-standard CTYPE3 spectral axis names to
    their canonical FITS WCS equivalents so that both cubes always hand
    the same spectral type string to astropy.

    Mappings applied:
    - ``VELO-LSR``, ``VELO-HEL``, ``VELO-OBS``  → ``VRAD``
      (these are old AIPS/CLASS conventions; astropy's spcfix already does
      this internally, but doing it explicitly here prevents inconsistencies
      when one cube has been fixed and the other hasn't)
    - ``VELOCITY``  → ``VRAD``  (another informal alias)
    """
    _VELO_ALIASES = {"VELO-LSR", "VELO-HEL", "VELO-OBS", "VELOCITY"}
    ctype3 = str(hdr.get("CTYPE3", "")).strip().upper()
    if ctype3 in _VELO_ALIASES:
        hdr["CTYPE3"] = "VRAD"
    return hdr


def _velo_grid(hdr):
    """
    Return the 1-D velocity array (m/s) for the spectral axis of hdr.

    Uses the standard FITS linear WCS: v[i] = CRVAL3 + (i+1 - CRPIX3)*CDELT3
    (channels are 1-based in FITS).
    """
    n = int(hdr["NAXIS3"])
    crval = float(hdr["CRVAL3"])
    cdelt = float(hdr["CDELT3"])
    crpix = float(hdr["CRPIX3"])
    return crval + (np.arange(1, n + 1) - crpix) * cdelt


def _spatial_footprints_overlap(src_hdr, tgt_hdr):
    """Return True if at least one corner of the source falls inside the target."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs_s = WCS(src_hdr).celestial
            wcs_t = WCS(tgt_hdr).celestial
        ny = int(src_hdr["NAXIS2"])
        nx = int(src_hdr["NAXIS1"])
        corners = np.array([[0, 0], [nx, 0], [0, ny], [nx, ny]], float)
        sky = wcs_s.pixel_to_world(corners[:, 0], corners[:, 1])
        px, py = wcs_t.world_to_pixel(sky)
        inside = (
            np.isfinite(px)
            & np.isfinite(py)
            & (px >= 0)
            & (px <= tgt_hdr["NAXIS1"])
            & (py >= 0)
            & (py <= tgt_hdr["NAXIS2"])
        )
        return bool(inside.any())
    except Exception:
        return False


def _resample_spectral_axis(cube, v_in, v_out, kind="linear"):
    """
    Resample *cube* (shape: n_in × ny × nx) along axis 0 from the velocity
    grid *v_in* onto *v_out*.

    Parameters
    ----------
    cube : ndarray, shape (n_in, ny, nx)
    v_in : 1-D array, length n_in  — source velocity grid (m/s)
    v_out : 1-D array, length n_out — target velocity grid (m/s)
    kind : str
        Interpolation kind passed to ``scipy.interpolate.interp1d``:
        ``'linear'`` (default), ``'nearest'``, ``'cubic'``, or ``'quadratic'``.

    Returns
    -------
    out : ndarray, shape (n_out, ny, nx)
        Channels outside *v_in* are NaN.
    """
    n_out = len(v_out)
    ny, nx = cube.shape[1], cube.shape[2]
    out = np.full((n_out, ny, nx), np.nan, dtype=np.float64)

    # Only interpolate pixels that have at least one finite value
    has_data = np.any(np.isfinite(cube), axis=0)  # (ny, nx)
    ys, xs = np.where(has_data)
    if len(ys) == 0:
        return out

    # Replace NaN with 0 for the interpolator (NaN-aware interp1d is slow);
    # the footprint mask will handle the spatial NaN correctly afterwards.
    spectra = cube[:, ys, xs].copy()  # (n_in, n_valid)
    spectra[~np.isfinite(spectra)] = 0.0

    interp = scipy.interpolate.interp1d(
        v_in,
        spectra,
        axis=0,
        kind=kind,
        bounds_error=False,
        fill_value=np.nan,
        assume_sorted=(v_in[0] < v_in[-1]),
    )
    out[:, ys, xs] = interp(v_out)
    return out


def _reproject_2d_planes(
    data,
    src_hdr,
    tgt_hdr,
    order,
    roundtrip_coords,
    block_size,
    parallel,
    spectral_order,
):
    """
    Reproject a 3-D cube by:
      1. Reprojecting each spatial plane onto the target spatial grid.
      2. Resampling the resulting cube along the spectral axis onto the
         target velocity grid.

    Returns (out_arr, out_fp) with shape (n_tgt_chan, ny_tgt, nx_tgt).
    """
    # Build pure 2-D spatial headers (celestial only, no spectral axis)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        src_hdr_2d = WCS(src_hdr).celestial.to_header()
        tgt_hdr_2d = WCS(tgt_hdr).celestial.to_header()

    for h, ref in ((src_hdr_2d, src_hdr), (tgt_hdr_2d, tgt_hdr)):
        h["NAXIS"] = 2
        h["NAXIS1"] = ref["NAXIS1"]
        h["NAXIS2"] = ref["NAXIS2"]

    # ── Step 1: spatial reprojection ────────────────────────────────────────
    n_src_chan = data.shape[0]
    ny_out = int(tgt_hdr["NAXIS2"])
    nx_out = int(tgt_hdr["NAXIS1"])
    spatial = np.full((n_src_chan, ny_out, nx_out), np.nan, dtype=np.float64)
    footprint = np.zeros((n_src_chan, ny_out, nx_out), dtype=np.float64)

    for i in range(n_src_chan):
        plane_hdu = fits.PrimaryHDU(data=data[i], header=src_hdr_2d)
        plane_out, plane_fp = reproject_interp(
            plane_hdu,
            tgt_hdr_2d,
            order=order,
            roundtrip_coords=roundtrip_coords,
            block_size=block_size,
            parallel=parallel,
            return_footprint=True,
        )
        spatial[i] = plane_out
        footprint[i] = plane_fp

    # ── Step 2: spectral resampling onto target velocity grid ────────────────
    v_src = _velo_grid(src_hdr)
    v_tgt = _velo_grid(tgt_hdr)

    out_arr = _resample_spectral_axis(spatial, v_src, v_tgt, kind=spectral_order)

    # Resample footprint the same way (nearest is fine for a binary mask)
    out_fp = _resample_spectral_axis(footprint, v_src, v_tgt, kind="nearest")
    out_fp = np.where(np.isfinite(out_fp), out_fp, 0.0)

    return out_arr, out_fp


# ─── Public API ──────────────────────────────────────────────────────────────


def reproject_cube(
    input_data,
    output_projection,
    shape_out=None,
    hdu_in=0,
    order="bilinear",
    roundtrip_coords=True,
    output_array=None,
    output_footprint=None,
    return_footprint=True,
    block_size=None,
    parallel=False,
    return_type=None,
    # dask_method=None,
    spectral_order="linear",
):
    """
    Drop-in replacement for ``reproject.reproject_interp`` with automatic
    fixes for various historic FITS cube conventions (e.g. GILDAS/CLASS).

    Parameters
    ----------
    input_data : str, Path, HDUList, PrimaryHDU, ImageHDU, (array, header/WCS),
                 or NDData
        The input data to reproject.  Accepts every type that
        ``reproject_interp`` accepts.
    output_projection : astropy.io.fits.Header or WCS
        The target projection / coordinate system.
    shape_out : tuple, optional
        Shape of the output array.  Required when output_projection is a WCS
        instance without embedded shape information.
    hdu_in : int or str, optional
        HDU index or name when input_data is a FITS file or HDUList.
        Default is 0.
    order : int or str, optional
        Spatial interpolation order: ``'nearest-neighbor'``, ``'bilinear'``
        (default), ``'biquadratic'``, or ``'bicubic'``.
    roundtrip_coords : bool, optional
        Whether to verify that coordinate transformations are defined in both
        directions.  Passed unchanged to ``reproject_interp``.
    output_array : ndarray or None, optional
        Pre-allocated array in which to store the result.
    output_footprint : ndarray or None, optional
        Pre-allocated array in which to store the footprint.
    return_footprint : bool, optional
        Whether to return the footprint alongside the data array.
        Default is True.
    block_size : tuple or ``'auto'``, optional
        Block size for tiled reprojection.  Passed to ``reproject_interp``.
    parallel : bool or int or str, optional
        Parallelism control.  Passed to ``reproject_interp``.
    return_type : {'numpy', 'dask'}, optional
        Return type for the output array.
    #dask_method : {'memmap', 'none'}, optional
    #    Method to use when the input array is a dask array.
    spectral_order : str, optional
        Interpolation kind used when resampling the spectral axis onto the
        target velocity grid.  Any value accepted by
        ``scipy.interpolate.interp1d``: ``'linear'`` (default),
        ``'nearest'``, ``'quadratic'``, or ``'cubic'``.
        Only used when source and target RESTFREQ values differ.

    Returns
    -------
    array_new : np.ndarray
        The reprojected data array, with shape matching the target header
        (NAXIS3 × NAXIS2 × NAXIS1).
    footprint : np.ndarray
        Coverage footprint (0 = no data, 1 = valid data).
        Only returned when ``return_footprint=True``.

    Notes
    -----
    The following fixes are applied automatically to both headers:

    **GLS → SFL**
        ``RA---GLS`` / ``DEC--GLS`` CTYPE values are renamed to
        ``RA---SFL`` / ``DEC--SFL``.  A ``CRPIX3 = 0`` value is corrected
        to ``1`` with a matching shift of ``CRVAL3``.

    **VELO-LSR normalisation**
        Deprecated spectral axis keywords (``VELO-LSR``, ``VELO-HEL``,
        ``VELO-OBS``, ``VELOCITY``) are renamed to ``VRAD`` so that both
        cubes always present the same CTYPE3 string to astropy, preventing
        inconsistencies from astropy's internal ``spcfix`` rewrite.

    **Spectral-axis isolation**
        When source and target headers carry detectably different rest
        frequencies (checked under both ``RESTFREQ`` and ``RESTFRQ`` key
        names; relative tolerance 1 × 10⁻⁶), the spectral axis is stripped
        from both WCS objects and each spatial plane is reprojected
        individually, avoiding the all-NaN failure caused by incommensurable
        rest frequencies.

    **Spectral resampling**
        After spatial reprojection the cube is resampled along axis 0 from the
        source velocity grid onto the target velocity grid using 1-D
        interpolation (``spectral_order``).  Channels outside the source
        velocity range are set to NaN.  The output array shape therefore
        exactly matches the target header (NAXIS3 × NAXIS2 × NAXIS1).

        A safety-net fallback applies the same full treatment whenever a
        standard 3-D reprojection returns ≥ 99 % NaN despite the spatial
        footprints overlapping.
    """
    # -- Extract arrays and headers ------------------------------------------
    data, src_hdr = _extract_data_header(input_data, hdu_in)
    tgt_hdr = _extract_output_header(output_projection, shape_out)

    # -- Apply header fixes to both headers ----------------------------------
    src_hdr = _fix_gls(src_hdr)
    tgt_hdr = _fix_gls(tgt_hdr)
    src_hdr = _normalize_spectral_ctype(src_hdr)
    tgt_hdr = _normalize_spectral_ctype(tgt_hdr)

    # -- Decide whether we need the 2-D + spectral-resample path -------------
    is_3d = (
        data.ndim == 3
        and int(src_hdr.get("NAXIS", data.ndim)) == 3
        and int(tgt_hdr.get("NAXIS", 3)) == 3
    )

    use_2d_path = False
    if is_3d:
        rf_src = _restfreq(src_hdr)
        rf_tgt = _restfreq(tgt_hdr)
        if (
            rf_src is not None
            and rf_tgt is not None
            and not np.isclose(rf_src, rf_tgt, rtol=1e-6)
        ):
            use_2d_path = True

    # -- 2-D spatial + spectral resample path --------------------------------
    if use_2d_path:
        out_arr, out_fp = _reproject_2d_planes(
            data,
            src_hdr,
            tgt_hdr,
            order,
            roundtrip_coords,
            block_size,
            parallel,
            spectral_order,
        )
        return (out_arr, out_fp) if return_footprint else out_arr

    # -- Standard path (GLS already fixed) -----------------------------------
    fixed_hdu = fits.PrimaryHDU(data=data, header=src_hdr)
    result = reproject_interp(
        fixed_hdu,
        tgt_hdr,
        shape_out=shape_out,
        order=order,
        roundtrip_coords=roundtrip_coords,
        output_array=output_array,
        output_footprint=output_footprint,
        return_footprint=return_footprint,
        block_size=block_size,
        parallel=parallel,
        return_type=return_type,
        # dask_method=dask_method,
    )

    # -- Safety net: fall back if still all-NaN despite spatial overlap ------
    if is_3d:
        arr_check = result[0] if return_footprint else result
        if np.isnan(arr_check).mean() > 0.99 and _spatial_footprints_overlap(
            src_hdr, tgt_hdr
        ):
            out_arr, out_fp = _reproject_2d_planes(
                data,
                src_hdr,
                tgt_hdr,
                order,
                roundtrip_coords,
                block_size,
                parallel,
                spectral_order,
            )
            return (out_arr, out_fp) if return_footprint else out_arr

    return result


def resolve_meta_resolution(source, params, meta, ov_hdr=None, log=None):
    """
    Resolve the per-source target resolution and write the results back into
    *meta* in-place.

    This is the single authoritative implementation used by all pipeline
    stages (regrid, products, fits) to ensure ``meta["target_res"]``,
    ``meta["target_res_pc"]``, and ``meta["res_suffix"]`` are correct for
    *source* before any stage-specific processing begins.

    Parameters
    ----------
    source  : str — source name (used only in log messages)
    params  : dict — source geometry from SourceHandler.get_source_params()
    meta    : dict — pipeline settings; updated in-place
    ov_hdr  : FITS Header or None — overlay cube header; required for
              ``resolution == "native"`` to read BMAJ/BMIN.  When None
              and the resolution is native, a warning is logged and the
              existing placeholder value is kept.
    log     : logger, optional

    Updated meta keys
    -----------------
    target_res    : float — target beam FWHM in arcseconds
    target_res_pc : float — target beam FWHM in parsecs
    res_suffix    : str   — filename suffix, e.g. "27p0as" or "100pc"
    """
    import math as _math
    if log is None:
        log = get_logger("Loading")

    resolution = meta.get("resolution", "angular")
    dist_mpc   = float(params.get("dist_mpc", 1.0))

    if resolution == "native":
        if ov_hdr is not None:
            candidate = max(ov_hdr.get("BMIN", 0), ov_hdr.get("BMAJ", 0)) * 3600.0
            if candidate > 0:
                target_res_as = candidate
                log.info(
                    f"Native resolution: {target_res_as:.1f} arcsec "
                    f"(from overlay header)."
                )
            else:
                target_res_as = meta.get("target_res", 27.0)
                log.warning(
                    f"Native resolution: no BMAJ/BMIN in overlay header; "
                    f"using placeholder {target_res_as:.1f} arcsec."
                )
        else:
            target_res_as = meta.get("target_res", 27.0)
            log.warning(
                f"Native resolution: overlay header not available; "
                f"using placeholder {target_res_as:.1f} arcsec."
            )

    elif resolution == "physical":
        target_res_pc_config = meta.get("target_res_config",
                                        meta.get("target_res", 27.0))
        target_res_as = (
            3600.0 * 180.0 / _math.pi * 1e-6 * target_res_pc_config / dist_mpc
        )
        log.info(
            f"Physical resolution: "
            f"{target_res_pc_config:.1f} pc = {target_res_as:.1f} arcsec "
            f"(distance {dist_mpc:.2f} Mpc)."
        )

    else:
        target_res_as = meta.get("target_res", 27.0)
        log.info(f"Angular resolution: {target_res_as:.1f} arcsec.")

    meta["target_res"]    = target_res_as
    meta["target_res_pc"] = target_res_as / 3600.0 * _math.pi / 180.0 * dist_mpc * 1e6
    if resolution == "physical":
        meta["res_suffix"] = (
            str(int(round(meta.get("target_res_config", target_res_as)))) + "pc"
        )
    else:
        meta["res_suffix"] = (
            str(np.round(target_res_as, 1)).replace(".", "p") + "as"
        )
