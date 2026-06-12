"""
pystructurePipeline.utils.fits_utils
=====================================
FITS helper utilities used throughout the pipeline.

Contains (all inline, no legacy imports):
  - FITS read/write helpers
  - twod_head           — strip a FITS header down to 2-D
  - hex_grid            — generate a hexagonal RA/Dec grid
  - deproject           — galactocentric deprojection
  - gaussian_PSF_2D     — 2-D Gaussian PSF kernel
  - deconvolve_gauss    — Gaussian deconvolution (MIRIAD port)
  - conv_with_gauss     — convolve a cube/map to a target resolution
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
from astropy.stats import median_absolute_deviation
from astropy.utils.console import ProgressBar

warnings.filterwarnings("ignore")


# ============================================================================
# Basic FITS I/O
# ============================================================================

def get_beam_arcsec(fits_path: str) -> au.Quantity:
    """Return the beam major axis in arcseconds from a FITS header."""
    fits_path = Path(fits_path)
    if not fits_path.exists():
        raise FileNotFoundError(f"FITS file not found: {fits_path}")
    hdr = fits.getheader(fits_path)
    if "BMAJ" not in hdr:
        raise KeyError(f"BMAJ not found in header of {fits_path}")
    return (hdr["BMAJ"] * au.deg).to(au.arcsec)


def read_fits_cube(fits_path: str):
    """Read a FITS cube, squeezing any degenerate 4th axis."""
    fits_path = Path(fits_path)
    if not fits_path.exists():
        raise FileNotFoundError(f"FITS file not found: {fits_path}")
    data, hdr = fits.getdata(fits_path, header=True)
    if hdr["NAXIS"] == 4:
        data = np.squeeze(data, axis=0)
        hdr["NAXIS"] = 3
        for key in ["NAXIS4", "CTYPE4", "CRVAL4", "CDELT4", "CRPIX4", "CUNIT4"]:
            hdr.remove(key, ignore_missing=True)
    return data, hdr


# ============================================================================
# twod_head
# ============================================================================

def twod_head(hdul_header):
    """
    Reduce a FITS header to 2-D by removing axes 3+.
    Port of IDL twod_head (leroy@mpia.de, 2008) by J. den Brok (2019).
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
# Hexagonal grid
# ============================================================================

def hex_grid(ctr_x, ctr_y, spacing, radec=False, r_limit=None, e_limit=None):
    """
    Generate a hexagonal grid centred on (ctr_x, ctr_y).

    Parameters
    ----------
    ctr_x, ctr_y : float — grid centre
    spacing      : float — point spacing (same units as the limit)
    radec        : bool  — if True, correct x for cos(Dec) foreshortening
    r_limit      : float — keep points within this radius
    e_limit      : float — keep points within this half-extent
    """
    x_spacing = spacing
    y_spacing = spacing * np.sin(np.deg2rad(60))

    if e_limit is None and r_limit is not None:
        scale = r_limit
    elif r_limit is None and e_limit is not None:
        scale = e_limit / 2
    else:
        raise TypeError("Provide either r_limit or e_limit to hex_grid.")

    half_ny = np.ceil(scale / y_spacing)
    half_nx = np.ceil(scale / x_spacing) + 1

    x = np.outer(np.ones(2 * int(half_ny) + 1), np.arange(2 * int(half_nx) + 1))
    y = np.outer(np.arange(2 * int(half_ny) + 1), np.ones(2 * int(half_nx) + 1))
    x -= half_nx
    y -= half_ny
    x *= x_spacing
    x += 0.5 * x_spacing * (np.dot(abs(y) % 2 == 1, 1))
    y *= y_spacing

    r = np.sqrt(x ** 2 + y ** 2)
    keep = np.where(r < r_limit) if r_limit is not None else \
           np.where(np.logical_and(abs(x) < e_limit / 2, abs(y) < e_limit / 2))
    if len(keep[0]) == 0:
        return np.nan, np.nan

    yout = y[keep] + ctr_y
    xout = (x[keep] / np.cos(np.deg2rad(yout)) + ctr_x) if radec else (x[keep] + ctr_x)
    return xout, yout


def make_sampling_points(ra_ctr, dec_ctr, max_rad, spacing, mask, hdr_mask,
                          overlay_in=None, overlay_hdr_in=None, show=False):
    """
    Generate hexagonal sampling points clipped to a binary mask.

    Port of sampling.py (J. den Brok, 2019).
    """
    n_dim_mask = len(np.shape(mask))
    if n_dim_mask == 3:
        print(f'{"[INFO]":<10}', 'Collapsing mask to two dimensions.')
        mask = np.sum(np.isfinite(mask), axis=0) >= 1
        hdr_mask = twod_head(hdr_mask)

    mask_dim = np.shape(mask)
    wcs = WCS(hdr_mask)

    if max_rad in ["auto"]:
        from astropy.coordinates import SkyCoord
        dx, dy = mask_dim[1], mask_dim[0]
        c1 = SkyCoord.from_pixel(0, 0, wcs)
        c2 = SkyCoord.from_pixel(dx, dy, wcs)
        max_rad = c1.separation(c2).value / 2
        print(f'{"[INFO]":<10}', f'Overlay size set to {np.round(max_rad, 3)} deg.')

    samp_ra, samp_dec = hex_grid(ra_ctr, dec_ctr, spacing, radec=True, r_limit=max_rad)

    try:
        pixel_coords = wcs.all_world2pix(np.column_stack((samp_ra, samp_dec)), 0)
    except Exception:
        pixel_coords = wcs.all_world2pix(
            np.column_stack((samp_ra, samp_dec, np.zeros(len(samp_ra)))), 0)

    samp_x = np.array(np.rint(pixel_coords[:, 0]), dtype=int)
    samp_y = np.array(np.rint(pixel_coords[:, 1]), dtype=int)

    keep = np.where(
        (samp_x >= 0) & (samp_y >= 0) &
        (samp_x < mask_dim[1]) & (samp_y < mask_dim[0])
    )[0]
    if len(keep) == 0:
        print(f'{"[ERROR]":<10}', 'No sampling points inside mask bounds. Returning NaNs.')
        return np.nan, np.nan

    samp_ra, samp_dec = samp_ra[keep], samp_dec[keep]
    samp_x, samp_y   = samp_x[keep],  samp_y[keep]

    keep = np.where(mask[samp_y, samp_x])[0]
    if len(keep) == 0:
        print(f'{"[ERROR]":<10}', 'No sampling points survive mask. Returning NaNs.')
        return np.nan, np.nan

    samp_ra, samp_dec = samp_ra[keep], samp_dec[keep]

    if show:
        _show_sampling_points(samp_ra, samp_dec, mask, hdr_mask,
                              overlay_in, overlay_hdr_in)

    return samp_ra, samp_dec


def _show_sampling_points(samp_ra, samp_dec, mask, hdr_mask,
                           overlay_in, overlay_hdr_in):
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
    Compute deprojected galactocentric radii and angles.

    Port of IDL deproject (A. Leroy 2001, J. den Brok 2019).

    Parameters
    ----------
    ra, dec  : arrays of RA/Dec (degrees)
    galpos   : [pa_deg, inc_deg, ra_ctr, dec_ctr]  (4 elements)
               or [vlsr, pa, inc, xctr, yctr]       (5 elements)
    vector   : bool — if True, inputs are already matched vectors
    gal      : dict with keys posang_deg, incl_def, ra_deg, dec_deg (optional)

    Returns
    -------
    rgrid : array of deprojected radii (degrees)
    tgrid : array of polar angles (radians)
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

    xgrid  = (rimg - xctr) * np.cos(np.deg2rad(yctr))
    ygrid  = dimg - yctr
    rotang = -(pa - np.pi / 2.0)

    deproj_x = xgrid * np.cos(rotang) + ygrid * np.sin(rotang)
    deproj_y = ygrid * np.cos(rotang) - xgrid * np.sin(rotang)
    deproj_y = deproj_y / np.cos(inc)

    rgrid = np.sqrt(deproj_x ** 2 + deproj_y ** 2)
    tgrid = np.arctan2(deproj_y, deproj_x)
    return rgrid, tgrid


# ============================================================================
# Gaussian PSF
# ============================================================================

def gaussian_PSF_2D(npix, a, center=False, normalize=False):
    """
    Create a 2-D rotated Gaussian PSF kernel.

    Parameters
    ----------
    npix      : int or [nx, ny]
    a         : [offset, peak, fwhm_major, fwhm_minor, cen_x, cen_y, rot_rad]
    center    : bool — centre the PSF in the array
    normalize : bool — normalise so it sums to 1
    """
    if isinstance(npix, (int, float)):
        nx = ny = int(npix)
    elif hasattr(npix, "__len__") and len(npix) == 2:
        nx, ny = int(npix[0]), int(npix[1])
    else:
        print("[ERROR]    Invalid npix.")
        return None

    xarr = np.tile(np.arange(nx), ny).reshape(ny, nx).astype(float)
    yarr = np.repeat(np.arange(ny), nx).reshape(ny, nx).astype(float)

    cenx = (nx - 1) / 2 if center else a[4]
    ceny = (ny - 1) / 2 if center else a[5]

    fac    = 2 * np.sqrt(2 * np.log(2))
    ang    = a[6]
    widthx = a[2] / fac
    widthy = a[3] / fac
    s, c   = np.sin(ang), np.cos(ang)

    xarr -= cenx
    yarr -= ceny
    t    = xarr * (c / widthx) + yarr * (s / widthx)
    yarr = xarr * (s / widthy) - yarr * (c / widthy)
    xarr = t

    output = a[0] + a[1] * np.exp(-0.5 * (xarr ** 2 + yarr ** 2))
    if normalize:
        output /= np.sum(output)
    return output


# ============================================================================
# Gaussian deconvolution
# ============================================================================

def deconvolve_gauss(meas_maj, beam_maj,
                     meas_min=None, meas_pa=None,
                     beam_min=None, beam_pa=None):
    """
    Deconvolve one Gaussian from another (MIRIAD gaupar.for port).

    Returns
    -------
    src_maj, src_min, src_pa, [worked, point_source]
    """
    if beam_min is None:
        meas_min = meas_maj
    if meas_pa is None:
        meas_pa = 0.0
    if beam_pa is None:
        beam_pa = 0.0

    mt = np.deg2rad(meas_pa)
    bt = np.deg2rad(beam_pa)

    alpha = ((meas_maj * np.cos(mt)) ** 2 + (meas_min * np.sin(mt)) ** 2
             - (beam_maj * np.cos(bt)) ** 2 - (beam_min * np.sin(bt)) ** 2)
    beta  = ((meas_maj * np.sin(mt)) ** 2 + (meas_min * np.cos(mt)) ** 2
             - (beam_maj * np.sin(bt)) ** 2 - (beam_min * np.cos(bt)) ** 2)
    gamma = 2 * ((meas_min ** 2 - meas_maj ** 2) * np.sin(mt) * np.cos(mt)
                 - (beam_min ** 2 - beam_maj ** 2) * np.sin(bt) * np.cos(bt))

    s = alpha + beta
    t = np.sqrt((alpha - beta) ** 2 + gamma ** 2)
    limit = 0.1 * min(meas_min or meas_maj, meas_maj, beam_maj, beam_min or beam_maj) ** 2

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
# Convolution with Gaussian beam
# ============================================================================

def _round_sig(x, sig=2):
    return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)


def _get_pixel_scale(hdr, tol=0.1):
    """Return pixel scale in degrees."""
    w = WCS(hdr)
    scales = proj_plane_pixel_scales(w)
    px_dx = scales[0] * au.deg
    px_dy = scales[1] * au.deg
    if abs(px_dx - px_dy) > tol * au.arcsec:
        print(f'{"[WARNING]":<10}', 'Pixel scale differs in X and Y.')
        return np.sqrt(px_dx * px_dy).value
    return px_dx.value


def _convolve_func(data, kernel, method="fft"):
    if method == "direct":
        return convolve(data, kernel, allow_huge=True)
    return convolve_fft(data, kernel, allow_huge=True)


def conv_with_gauss(in_data, in_hdr=None, start_beam=None, pix_deg=None,
                    target_beam=None, no_ft=False, in_weight=None,
                    out_weight_file=None, out_file=None,
                    unc=False, perbeam=False, quiet=False):
    """
    Convolve a cube or map to a target Gaussian beam.

    Port of IDL conv_with_gauss (leroy@mpia.de, J. den Brok 2020).

    Parameters
    ----------
    in_data      : np.ndarray or str (FITS path)
    in_hdr       : FITS header (required if in_data is an array)
    target_beam  : [maj_as, min_as, pa_deg] — target beam in arcseconds
    unc          : bool — treat as uncertainty map (square before convolution)
    perbeam      : bool — correct per-beam units after convolution
    quiet        : bool — suppress reporting
    """
    if isinstance(in_data, str):
        data, hdr = fits.getdata(in_data, header=True)
    else:
        data = copy.deepcopy(in_data)
        hdr  = in_hdr

    if target_beam is not None:
        if isinstance(target_beam, (list, np.ndarray)):
            if len(target_beam) == 1:
                target_beam = np.array([target_beam[0], target_beam[0], 0.0])
            elif len(target_beam) == 2:
                target_beam = np.array([target_beam[0], target_beam[1], 0.0])
        elif isinstance(target_beam, (float, int)):
            target_beam = np.array([float(target_beam), float(target_beam), 0.0])

    flux_before = np.nansum(data)

    if pix_deg is None:
        pix_deg = _get_pixel_scale(hdr, tol=0.1 * hdr.get("BMIN", 1.0))
    as_per_pix = pix_deg * 3600.0

    if unc:
        data = data ** 2

    # --- Identify starting beam ---
    if start_beam is not None:
        if isinstance(start_beam, float):
            current_beam = [start_beam, start_beam, 0.0]
        elif hasattr(start_beam, "__len__"):
            current_beam = list(start_beam) + [0.0] * (3 - len(start_beam))
        else:
            print(f'{"[ERROR]":<10}', 'Unknown start_beam format.')
            return None, None
    else:
        bmaj = hdr["BMAJ"] * 3600
        bmin = hdr["BMIN"] * 3600
        bpa  = hdr.get("BPA", 0.0)
        current_beam = [bmaj, bmin, bpa]

    src_maj, src_min, src_pa, info = deconvolve_gauss(
        target_beam[0], current_beam[0],
        target_beam[1], target_beam[2],
        current_beam[1], current_beam[2],
    )
    if not info[0]:
        print(f'{"[ERROR]":<10}', 'Cannot compute convolution kernel for these beams.')
        return None, None
    if info[1]:
        print(f'{"[WARNING]":<10}', 'Target and starting beam are very close.')

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

    if len(dim_data) == 3:
        new_data = copy.deepcopy(data)
        print(f'{"[INFO]":<10}', 'Start cube convolution:')
        for plane in ProgressBar(range(dim_data[0])):
            new_data[plane, :, :] = _convolve_func(data[plane, :, :], kernel, method)
        data = new_data
    else:
        data = _convolve_func(data, kernel, method)

    # Beam-area correction terms
    if unc or perbeam:
        cur_fwhm    = np.sqrt(current_beam[0] * current_beam[1])
        ppbeam_start = (cur_fwhm / as_per_pix / 2) ** 2 / np.log(2) * np.pi
        tgt_fwhm    = np.sqrt(target_beam[0] * target_beam[1])
        ppbeam_final = (tgt_fwhm / as_per_pix / 2) ** 2 / np.log(2) * np.pi

    if unc:
        data = np.sqrt(data) * np.sqrt(ppbeam_start / ppbeam_final)
    if perbeam:
        data *= ppbeam_final / ppbeam_start

    if not quiet:
        print(f'{"[INFO]":<10}', f'Pixel scale [as] = {_round_sig(as_per_pix, 3)}')
        print(f'{"[INFO]":<10}', f'Starting beam [as] = {current_beam}')
        print(f'{"[INFO]":<10}', f'Target FWHM [as] = {target_beam}')
        print(f'{"[INFO]":<10}', f'Flux ratio = {_round_sig(np.nansum(data) / flux_before)}')

    hdr["BMAJ"] = (target_beam[0] / 3600.0, "FWHM BEAM IN DEGREES")
    hdr["BMIN"] = (target_beam[1] / 3600.0, "FWHM BEAM IN DEGREES")
    hdr["BPA"]  = (target_beam[2], "POSITION ANGLE IN DEGREES")
    hdr["HISTORY"] = f"conv_with_gauss: convolved with [{src_maj:.2f}, {src_min:.2f}, {src_pa:.2f}] arcsec Gaussian"

    if out_file is not None:
        fits.writeto(out_file, data, hdr, overwrite=True)

    return data, hdr
