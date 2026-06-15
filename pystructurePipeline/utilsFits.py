"""
utilsFits.py — FITS and astronomy utility functions.

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
from astropy import units as au
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel
from astropy.utils.console import ProgressBar

warnings.filterwarnings("ignore")

_LOG_PREFIX = "[pyStructure] [Utils]   "


# ============================================================================
# Basic FITS I/O
# ============================================================================

def get_beam_arcsec(fits_path: str) -> au.Quantity:
    """
    Return the beam major axis (BMAJ) in arcseconds from a FITS header.

    Parameters
    ----------
    fits_path : str or Path

    Returns
    -------
    beam_as : astropy.units.Quantity in arcseconds

    Raises
    ------
    FileNotFoundError : if the file does not exist
    KeyError          : if BMAJ is absent from the header
    """
    fits_path = Path(fits_path)
    if not fits_path.exists():
        raise FileNotFoundError(f"{_LOG_PREFIX} [ERROR]  FITS file not found: {fits_path}")
    hdr = fits.getheader(fits_path)
    if "BMAJ" not in hdr:
        raise KeyError(f"{_LOG_PREFIX} [ERROR]  BMAJ not found in header of {fits_path}")
    return (hdr["BMAJ"] * au.deg).to(au.arcsec)


def read_fits_cube(fits_path: str):
    """
    Read a FITS file, squeezing any degenerate Stokes (4th) axis.

    Many radio FITS cubes have a 4th Stokes axis of length 1.  This function
    removes it so that the returned array is always 3-D (channels × y × x).

    Parameters
    ----------
    fits_path : str or Path

    Returns
    -------
    data : np.ndarray (3D)
    hdr  : astropy.io.fits.Header (updated to reflect 3D)

    Raises
    ------
    FileNotFoundError : if the file does not exist
    """
    fits_path = Path(fits_path)
    if not fits_path.exists():
        raise FileNotFoundError(f"{_LOG_PREFIX} [ERROR]  FITS file not found: {fits_path}")
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

    r    = np.sqrt(x**2 + y**2)
    keep = (np.where(r < r_limit) if r_limit is not None
            else np.where(np.logical_and(abs(x) < e_limit / 2, abs(y) < e_limit / 2)))

    if len(keep[0]) == 0:
        return np.nan, np.nan

    yout = y[keep] + ctr_y
    xout = (x[keep] / np.cos(np.deg2rad(yout)) + ctr_x) if radec else (x[keep] + ctr_x)
    return xout, yout


def make_sampling_points(ra_ctr, dec_ctr, max_rad, spacing, mask, hdr_mask,
                          overlay_in=None, overlay_hdr_in=None, show=False):
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

    Returns
    -------
    samp_ra, samp_dec : np.ndarrays — coordinates of the surviving grid points.
                        Returns (np.nan, np.nan) if no points survive.
    """
    # Collapse 3-D mask to 2-D
    if len(np.shape(mask)) == 3:
        print(f"{_LOG_PREFIX} [INFO]  Collapsing 3D mask to 2D footprint.")
        mask     = np.sum(np.isfinite(mask), axis=0) >= 1
        hdr_mask = twod_head(hdr_mask)

    mask_dim = np.shape(mask)
    wcs      = WCS(hdr_mask)

    # Auto-determine maximum radius from the mask array diagonal
    if max_rad == "auto":
        from astropy.coordinates import SkyCoord
        c1  = SkyCoord.from_pixel(0,           0,           wcs)
        c2  = SkyCoord.from_pixel(mask_dim[1], mask_dim[0], wcs)
        max_rad = c1.separation(c2).value / 2
        print(f"{_LOG_PREFIX} [INFO]  Auto max_rad = {np.round(max_rad, 3)} deg.")

    samp_ra, samp_dec = hex_grid(ra_ctr, dec_ctr, spacing, radec=True, r_limit=max_rad)

    # Convert to pixel coordinates
    try:
        pixel_coords = wcs.all_world2pix(np.column_stack((samp_ra, samp_dec)), 0)
    except Exception:
        pixel_coords = wcs.all_world2pix(
            np.column_stack((samp_ra, samp_dec, np.zeros(len(samp_ra)))), 0)

    samp_x = np.array(np.rint(pixel_coords[:, 0]), dtype=int)
    samp_y = np.array(np.rint(pixel_coords[:, 1]), dtype=int)

    # Keep only points inside the array boundary
    keep = np.where(
        (samp_x >= 0) & (samp_y >= 0) &
        (samp_x < mask_dim[1]) & (samp_y < mask_dim[0])
    )[0]
    if len(keep) == 0:
        print(f"{_LOG_PREFIX} [ERROR]  No sampling points inside mask bounds.")
        return np.nan, np.nan

    samp_ra, samp_dec = samp_ra[keep], samp_dec[keep]
    samp_x, samp_y   = samp_x[keep],  samp_y[keep]

    # Keep only points where the mask is True
    keep = np.where(mask[samp_y, samp_x])[0]
    if len(keep) == 0:
        print(f"{_LOG_PREFIX} [ERROR]  No sampling points survive mask clipping.")
        return np.nan, np.nan

    samp_ra, samp_dec = samp_ra[keep], samp_dec[keep]

    if show:
        _show_sampling_points(samp_ra, samp_dec, mask, hdr_mask, overlay_in, overlay_hdr_in)

    return samp_ra, samp_dec


def _show_sampling_points(samp_ra, samp_dec, mask, hdr_mask, overlay_in, overlay_hdr_in):
    """Visualise sampling points overlaid on the mask or overlay image."""
    import matplotlib.pyplot as plt
    if overlay_in is not None:
        if isinstance(overlay_in, str):
            overlay, overlay_hdr = fits.getdata(overlay_in, header=True)
        else:
            overlay     = copy.deepcopy(overlay_in)
            overlay_hdr = overlay_hdr_in if overlay_hdr_in is not None else hdr_mask
        if len(np.shape(overlay)) == 3:
            overlay     = np.nansum(overlay, 0)
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
                           ra_deg, dec_deg

    Returns
    -------
    rgrid : np.ndarray — deprojected radius in degrees (same shape as ra/dec)
    tgrid : np.ndarray — deprojected polar angle in radians
    """
    np.seterr(divide="ignore", invalid="ignore")

    if gal is not None:
        pa   = np.deg2rad(gal["posang_deg"])
        inc  = np.deg2rad(gal["incl_def"])
        xctr = gal["ra_deg"]
        yctr = gal["dec_deg"]
    elif len(galpos) == 5:
        pa   = np.deg2rad(galpos[1])
        inc  = np.deg2rad(galpos[2])
        xctr = galpos[3]
        yctr = galpos[4]
    else:
        pa   = np.deg2rad(galpos[0])
        inc  = np.deg2rad(galpos[1])
        xctr = galpos[2]
        yctr = galpos[3]

    ra_size = np.shape(ra)
    if ra_size[0] == 1 and not vector:
        rimg = np.outer(ra, np.ones(len(dec)))
        dimg = np.outer(np.ones(len(ra)), dec)
    else:
        rimg, dimg = ra, dec

    # Offset from centre, correcting RA for cos(Dec) projection
    xgrid  = (rimg - xctr) * np.cos(np.deg2rad(yctr))
    ygrid  = dimg - yctr

    # Rotate by (PA - 90°) to align with the major axis
    rotang = -(pa - np.pi / 2.0)
    deproj_x =  xgrid * np.cos(rotang) + ygrid * np.sin(rotang)
    deproj_y =  ygrid * np.cos(rotang) - xgrid * np.sin(rotang)

    # Stretch along the minor axis to correct for inclination
    deproj_y = deproj_y / np.cos(inc)

    rgrid = np.sqrt(deproj_x**2 + deproj_y**2)
    tgrid = np.arctan2(deproj_y, deproj_x)
    return rgrid, tgrid


# ============================================================================
# Gaussian PSF
# ============================================================================

def gaussian_PSF_2D(npix, a, center=False, normalize=False):
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

    Returns
    -------
    output : np.ndarray (ny, nx) — the Gaussian kernel
    """
    if isinstance(npix, (int, float)):
        nx = ny = int(npix)
    elif hasattr(npix, "__len__") and len(npix) == 2:
        nx, ny = int(npix[0]), int(npix[1])
    else:
        print(f"{_LOG_PREFIX} [ERROR]  Invalid npix: {npix}")
        return None

    xarr = np.tile(np.arange(nx), ny).reshape(ny, nx).astype(float)
    yarr = np.repeat(np.arange(ny), nx).reshape(ny, nx).astype(float)

    cenx = (nx - 1) / 2 if center else a[4]
    ceny = (ny - 1) / 2 if center else a[5]

    fac    = 2 * np.sqrt(2 * np.log(2))   # FWHM → sigma conversion
    ang    = a[6]
    widthx = a[2] / fac
    widthy = a[3] / fac
    s, c   = np.sin(ang), np.cos(ang)

    xarr -= cenx
    yarr -= ceny
    t    = xarr * (c / widthx) + yarr * (s / widthx)
    yarr = xarr * (s / widthy) - yarr * (c / widthy)
    xarr = t

    output = a[0] + a[1] * np.exp(-0.5 * (xarr**2 + yarr**2))
    if normalize:
        output /= np.sum(output)
    return output


# ============================================================================
# Beam deconvolution
# ============================================================================

def deconvolve_gauss(meas_maj, beam_maj,
                     meas_min=None, meas_pa=None,
                     beam_min=None, beam_pa=None):
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
    if beam_min is None: meas_min = meas_maj
    if meas_pa  is None: meas_pa  = 0.0
    if beam_pa  is None: beam_pa  = 0.0

    mt = np.deg2rad(meas_pa)
    bt = np.deg2rad(beam_pa)

    alpha = ((meas_maj * np.cos(mt))**2 + (meas_min * np.sin(mt))**2
             - (beam_maj * np.cos(bt))**2 - (beam_min * np.sin(bt))**2)
    beta  = ((meas_maj * np.sin(mt))**2 + (meas_min * np.cos(mt))**2
             - (beam_maj * np.sin(bt))**2 - (beam_min * np.cos(bt))**2)
    gamma = 2 * ((meas_min**2 - meas_maj**2) * np.sin(mt) * np.cos(mt)
                 - (beam_min**2 - beam_maj**2) * np.sin(bt) * np.cos(bt))

    s     = alpha + beta
    t     = np.sqrt((alpha - beta)**2 + gamma**2)
    limit = 0.1 * min(meas_min or meas_maj, meas_maj, beam_maj, beam_min or beam_maj)**2

    if alpha < 0 or beta < 0 or s < t:
        worked = False
        point  = (0.5 * (s - t) < limit) and (alpha > -limit) and (beta > -limit)
        return 0.0, 0.0, 0.0, [worked, point]

    src_maj = np.sqrt(0.5 * (s + t))
    src_min = np.sqrt(0.5 * (s - t))
    src_pa  = (0.0 if (abs(gamma) + abs(alpha - beta)) == 0
               else np.rad2deg(0.5 * np.arctan(-gamma / (alpha - beta))))
    return src_maj, src_min, src_pa, [True, False]


# ============================================================================
# Spatial convolution
# ============================================================================

def _round_sig(x, sig=2):
    """Round *x* to *sig* significant figures."""
    return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)


def _get_pixel_scale(hdr, tol=0.1):
    """
    Return the pixel scale in degrees from a FITS header.

    Issues a warning if the x and y pixel scales differ by more than *tol*
    arcseconds and returns the geometric mean in that case.
    """
    w      = WCS(hdr)
    scales = proj_plane_pixel_scales(w)
    px_dx  = scales[0] * au.deg
    px_dy  = scales[1] * au.deg
    if abs(px_dx - px_dy) > tol * au.arcsec:
        print(f"{_LOG_PREFIX} [WARNING]  Pixel scale differs in X and Y: "
              f"{px_dx.to(au.arcsec):.3f} vs {px_dy.to(au.arcsec):.3f}. "
              "Using geometric mean.")
        return np.sqrt(px_dx * px_dy).value
    return px_dx.value


def _convolve_func(data, kernel, method="fft"):
    """Dispatch to the appropriate astropy convolution function."""
    if method == "direct":
        return convolve(data, kernel, allow_huge=True)
    return convolve_fft(data, kernel, allow_huge=True)


def conv_with_gauss(in_data, in_hdr=None, start_beam=None, pix_deg=None,
                    target_beam=None, no_ft=False, in_weight=None,
                    out_weight_file=None, out_file=None,
                    unc=False, perbeam=False, quiet=False):
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

    Returns
    -------
    data : np.ndarray — convolved data
    hdr  : FITS Header — updated with new BMAJ/BMIN/BPA keywords
    Returns (None, None) if the deconvolution fails.
    """
    # Load data
    if isinstance(in_data, str):
        data, hdr = fits.getdata(in_data, header=True)
    else:
        data = copy.deepcopy(in_data)
        hdr  = in_hdr

    # Normalise target_beam to a 3-element array [maj, min, pa]
    if target_beam is not None:
        if isinstance(target_beam, (float, int)):
            target_beam = np.array([float(target_beam), float(target_beam), 0.0])
        elif hasattr(target_beam, "__len__"):
            target_beam = np.array(list(target_beam) + [0.0] * (3 - len(target_beam)))

    flux_before = np.nansum(data)

    # Pixel scale
    if pix_deg is None:
        pix_deg = _get_pixel_scale(hdr, tol=0.1 * hdr.get("BMIN", 1.0))
    as_per_pix = pix_deg * 3600.0

    if unc:
        data = data**2   # square uncertainties before convolution

    # Identify the starting beam
    if start_beam is not None:
        if isinstance(start_beam, float):
            current_beam = [start_beam, start_beam, 0.0]
        elif hasattr(start_beam, "__len__"):
            current_beam = list(start_beam) + [0.0] * (3 - len(start_beam))
        else:
            print(f"{_LOG_PREFIX} [ERROR]  Unknown start_beam format.")
            return None, None
    else:
        bmaj = hdr["BMAJ"] * 3600
        bmin = hdr["BMIN"] * 3600
        bpa  = hdr.get("BPA", 0.0)
        current_beam = [bmaj, bmin, bpa]

    # Compute the convolution kernel by deconvolving the current beam
    src_maj, src_min, src_pa, info = deconvolve_gauss(
        target_beam[0], current_beam[0],
        target_beam[1], target_beam[2],
        current_beam[1], current_beam[2],
    )
    if not info[0]:
        print(f"{_LOG_PREFIX} [ERROR]  Cannot compute convolution kernel: "
              f"target beam {target_beam} is smaller than current beam {current_beam}.")
        return None, None
    if info[1]:
        print(f"{_LOG_PREFIX} [WARNING]  Target and starting beam are nearly identical; "
              "kernel will be very small.")

    # Build the kernel array (6 × FWHM in pixels, capped to the data size)
    kern_size = int(6.0 * np.rint(src_maj / as_per_pix) + 1)
    dim_data  = np.shape(data)
    dim_x = dim_data[1] if len(dim_data) == 3 else dim_data[0]
    dim_y = dim_data[2] if len(dim_data) == 3 else dim_data[1]
    if kern_size > dim_x or kern_size > dim_y:
        kern_size = int(np.floor(min(dim_x, dim_y) / 2 - 2) * 2 + 1)

    kernel = gaussian_PSF_2D(
        kern_size,
        [0., 1., src_maj / as_per_pix, src_min / as_per_pix, 0., 0.,
         np.pi / 2.0 + np.deg2rad(src_pa)],
        center=True, normalize=True,
    )

    method = "direct" if no_ft else "fft"

    # Convolve: plane by plane for cubes, direct for 2-D maps
    if len(dim_data) == 3:
        new_data = copy.deepcopy(data)
        if not quiet:
            print(f"{_LOG_PREFIX} [INFO]  Convolving cube ({dim_data[0]} planes):")
        for plane in ProgressBar(range(dim_data[0])):
            new_data[plane, :, :] = _convolve_func(data[plane, :, :], kernel, method)
        data = new_data
    else:
        data = _convolve_func(data, kernel, method)

    # Beam-area correction for per-beam or uncertainty maps
    if unc or perbeam:
        cur_fwhm      = np.sqrt(current_beam[0] * current_beam[1])
        ppbeam_start  = (cur_fwhm / as_per_pix / 2)**2 / np.log(2) * np.pi
        tgt_fwhm      = np.sqrt(target_beam[0] * target_beam[1])
        ppbeam_final  = (tgt_fwhm / as_per_pix / 2)**2 / np.log(2) * np.pi

    if unc:
        data = np.sqrt(data) * np.sqrt(ppbeam_start / ppbeam_final)
    if perbeam:
        data *= ppbeam_final / ppbeam_start

    if not quiet:
        print(f"{_LOG_PREFIX} [INFO]  Pixel scale: {_round_sig(as_per_pix, 3)} arcsec/px  |  "
              f"Input beam: {[round(b, 1) for b in current_beam]} arcsec  |  "
              f"Target beam: {[round(b, 1) for b in target_beam]} arcsec  |  "
              f"Flux ratio: {_round_sig(np.nansum(data) / flux_before)}")

    # Update header with the new beam
    hdr["BMAJ"]    = (target_beam[0] / 3600.0, "FWHM BEAM IN DEGREES")
    hdr["BMIN"]    = (target_beam[1] / 3600.0, "FWHM BEAM IN DEGREES")
    hdr["BPA"]     = (target_beam[2],           "POSITION ANGLE IN DEGREES")
    hdr["HISTORY"] = (f"conv_with_gauss: convolved with "
                      f"[{src_maj:.2f}, {src_min:.2f}, {src_pa:.2f}] arcsec kernel")

    if out_file is not None:
        fits.writeto(out_file, data, hdr, overwrite=True)

    return data, hdr
