"""
stage_regrid: convolve and sample all bands and cubes onto the hex grid.

All logic is contained here — no imports from the legacy scripts/ directory.
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
from astropy.stats import median_absolute_deviation
from reproject import reproject_interp

from pystructurePipeline.utilsFits import (
    twod_head, conv_with_gauss, deproject,
)
from pystructurePipeline.pipelineSampling import run_sampling

warnings.filterwarnings("ignore")


# ============================================================================
# Velocity-axis helpers
# ============================================================================

def _get_vaxis(hdr):
    v = np.arange(hdr["NAXIS3"])
    return (v - (hdr["CRPIX3"] - 1)) * hdr["CDELT3"] + hdr["CRVAL3"]


def _ensure_ms(hdr, data=None):
    """Convert velocity axis from km/s to m/s in-place if needed."""
    if abs(hdr["CDELT3"]) < 200:
        hdr["CDELT3"] *= 1000
        hdr["CRVAL3"] *= 1000
        hdr["CUNIT3"] = "m/s"
    if data is not None and hdr["CDELT3"] < 0:
        vaxis_inv = _get_vaxis(hdr)
        hdr["CDELT3"] = abs(hdr["CDELT3"])
        hdr["CRPIX3"] = 1
        hdr["CRVAL3"] = vaxis_inv[-1]
        data = np.flip(data, axis=0)
    return hdr, data


# ============================================================================
# Spectral smoothing
# ============================================================================

def _spectral_smooth(data, hdr_out, spec_smooth):
    """Apply spectral smoothing to *data* according to *spec_smooth*."""
    from astropy.convolution import Gaussian1DKernel, convolve

    mode, method = spec_smooth[0], spec_smooth[1]
    if mode in ["default"]:
        return data, hdr_out

    spec_res = abs(hdr_out["CDELT3"]) / 1000.0   # km/s
    fwhm_factor = np.sqrt(8 * np.log(2))
    dim_data = np.shape(data)

    if not (isinstance(mode, (int, float))):
        return data, hdr_out

    if spec_res >= mode:
        print(f'{"[INFO]":<10}', 'No spectral smoothing; already at target resolution.')
        return data, hdr_out

    print(f'{"[INFO]":<10}', f'Spectral smoothing to {round(mode, 3)} km/s ({method}).')

    if method == "gauss":
        pix = ((mode ** 2 - spec_res ** 2) ** 0.5 / spec_res) / fwhm_factor
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
            new_vaxis = np.array([np.nanmean(vaxis[n_ratio*j:n_ratio*(j+1)])
                                   for j in range(new_len)])
            data = np.array([np.nanmean(data[n_ratio*j:n_ratio*(j+1), :, :], axis=0)
                              for j in range(new_len)])
            hdr_out["NAXIS3"] = new_len
            hdr_out["CDELT3"] = new_vaxis[1] - new_vaxis[0]
            hdr_out["CRVAL3"] = new_vaxis[0] + (hdr_out["CRPIX3"] - 1) * hdr_out["CDELT3"]

        if method == "combined" and n_ratio * spec_res < mode:
            pix = ((mode ** 2 - (n_ratio * spec_res) ** 2) ** 0.5 / spec_res) / fwhm_factor
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
    Convolve *in_data* to *target_res_as* and sample at the hex-grid points.

    Port of sampling_at_resol.sample_at_res (J. den Brok / L. Neumann).
    """
    if len(ra_samp) != len(dec_samp):
        print(f'{"[ERROR]":<10}', 'RA and Dec arrays must have the same length.')
        return ra_samp * np.nan, None

    if isinstance(in_data, str):
        if not path.exists(in_data):
            print(f'{"[ERROR]":<10}', f'File not found: {in_data}')
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

    if not is_cube and trg_hdr is not None:
        trg_hdr = twod_head(trg_hdr)

    # --- Convolution ---
    if "BMAJ" not in hdr:
        print(f'{"[WARNING]":<10}', 'No BMAJ in header; skipping convolution.')
        hdr_out = copy.copy(hdr)
    elif hdr["BMAJ"] < 0.99 * target_res_as / 3600.0:
        print(f'{"[INFO]":<10}', f'Convolving to {round(target_res_as, 3)} arcsec.')
        data, hdr_out = conv_with_gauss(
            in_data=data, in_hdr=hdr,
            target_beam=target_res_as * np.array([1.0, 1.0, 0.0]),
            quiet=True, perbeam=perbeam, unc=unc,
        )
    else:
        print(f'{"[INFO]":<10}', 'Already at target resolution; skipping convolution.')
        hdr_out = copy.copy(hdr)

    # --- Velocity-axis checks for cubes ---
    if is_cube:
        if trg_hdr is not None:
            trg_hdr, _ = _ensure_ms(trg_hdr)
        hdr_out, data = _ensure_ms(hdr_out, data)

        # Flip inverted velocity axis
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

    # --- Spectral smoothing ---
    data, hdr_out = _spectral_smooth(data, hdr_out, list(spec_smooth))

    # --- Reproject ---
    if trg_hdr is not None:
        for key in ["RESTF*"]:
            hdr_out.remove(key, ignore_missing=True)
            trg_hdr.remove(key, ignore_missing=True)

        # Adjust target spectral axis for smoothing
        if not spec_smooth[0] in ("default",) and isinstance(spec_smooth[0], (int, float)):
            if spec_smooth[0] > trg_hdr.get("CDELT3", 0) / 1000.0:
                vaxis_ov  = _get_vaxis(trg_hdr)
                new_vaxis = np.arange(vaxis_ov[0], vaxis_ov[-1], spec_smooth[0] * 1000)
                trg_hdr["NAXIS3"] = len(new_vaxis)
                trg_hdr["CDELT3"] = spec_smooth[0] * 1000
                trg_hdr["CRVAL3"] = new_vaxis[0] + (trg_hdr["CRPIX3"] - 1) * trg_hdr["CDELT3"]

        data_out, _ = reproject_interp((data, hdr_out), trg_hdr, order="nearest-neighbor")
        data = data_out

        if save_fits:
            out_hdr = copy.copy(trg_hdr)
            out_hdr["BMAJ"] = target_res_as / 3600.0
            out_hdr["BMIN"] = target_res_as / 3600.0
            out_hdr["LINE"] = line_name
            fits.writeto(
                path.join(path_save_fits, f"{source}_{line_name}_{target_res_as}as.fits"),
                data=data, header=out_hdr, overwrite=True,
            )
    else:
        print(f'{"[INFO]":<10}', 'No alignment; no target header supplied.')
        trg_hdr = hdr_out

    # --- Sample ---
    wcs_t = WCS(trg_hdr)
    if is_cube:
        pixel_coords = wcs_t.all_world2pix(
            np.column_stack((ra_samp, dec_samp, np.zeros(len(dec_samp)))), 0)
    else:
        pixel_coords = wcs_t.all_world2pix(np.column_stack((ra_samp, dec_samp)), 0)

    samp_x = np.array(np.rint(pixel_coords[:, 0]), dtype=int)
    samp_y = np.array(np.rint(pixel_coords[:, 1]), dtype=int)
    n_pts  = len(samp_x)
    dim_data = np.shape(data)

    result = np.zeros((n_pts, dim_data[0])) * np.nan if is_cube else np.zeros(n_pts) * np.nan

    if is_cube:
        in_map = np.where(
            (samp_x > 0) & (samp_x < dim_data[2]) &
            (samp_y > 0) & (samp_y < dim_data[1])
        )[0]
        for kk in in_map:
            result[kk, :] = data[:, samp_y[kk], samp_x[kk]]
    else:
        in_map = np.where(
            (samp_x > 0) & (samp_x < dim_data[1]) &
            (samp_y > 0) & (samp_y < dim_data[0])
        )[0]
        result[in_map] = data[samp_y[in_map], samp_x[in_map]]

    return result, trg_hdr


def sample_mask(in_data, ra_samp, dec_samp, in_hdr=None, target_hdr=None):
    """Sample a binary mask cube/image onto the hex-grid points."""
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
        for key in ["RESTF*"]:
            hdr_out.remove(key, ignore_missing=True)
            trg_hdr.remove(key, ignore_missing=True)
        data_out, _ = reproject_interp((data, hdr_out), trg_hdr, order="nearest-neighbor")
        data = data_out

    wcs_t = WCS(trg_hdr)
    if is_cube:
        pixel_coords = wcs_t.all_world2pix(
            np.column_stack((ra_samp, dec_samp, np.zeros(len(dec_samp)))), 0)
    else:
        pixel_coords = wcs_t.all_world2pix(np.column_stack((ra_samp, dec_samp)), 0)

    samp_x = np.array(np.rint(pixel_coords[:, 0]), dtype=int)
    samp_y = np.array(np.rint(pixel_coords[:, 1]), dtype=int)
    n_pts  = len(samp_x)
    dim_data = np.shape(data)

    result = np.zeros((n_pts, dim_data[0])) * np.nan if is_cube else np.zeros(n_pts) * np.nan

    if is_cube:
        in_map = np.where(
            (samp_x > 0) & (samp_x < dim_data[2]) &
            (samp_y > 0) & (samp_y < dim_data[1])
        )[0]
        for kk in in_map:
            result[kk, :] = data[:, samp_y[kk], samp_x[kk]]
    else:
        in_map = np.where(
            (samp_x > 0) & (samp_x < dim_data[1]) &
            (samp_y > 0) & (samp_y < dim_data[0])
        )[0]
        result[in_map] = data[samp_y[in_map], samp_x[in_map]]

    return result, trg_hdr


# ============================================================================
# Stage entry point
# ============================================================================

def run_regrid(source, params, meta, maps, cubes, input_mask):
    """
    Convolve and sample all bands and cubes for *source*.

    Writes the resulting Astropy Table to disk as an .ecsv file and
    returns the file path.
    """
    from pystructurePipeline import __version__, __author__, __email__, __credits__

    sampling      = run_sampling(source=source, params=params, meta=meta)
    samp_ra       = sampling["samp_ra"]
    samp_dec      = sampling["samp_dec"]
    ov_hdr        = sampling["ov_hdr"]
    target_res_as = sampling["target_res_as"]
    n_chan         = ov_hdr["NAXIS3"]
    n_pts          = len(samp_ra)

    fname = _build_fname(source, meta)
    structure_creation = meta.get("structure_creation", "default")
    data_dir  = meta.get("data_dir", "data/")
    save_fits = meta.get("save_fits", False)

    if "fill" in structure_creation and path.exists(fname):
        this_data, fill_maps, fill_cubes = _fill_checker(fname, samp_ra, samp_dec, maps, cubes)
    else:
        this_data = _init_table(source, params, meta, samp_ra, samp_dec, ov_hdr, target_res_as,
                                __version__, __author__, __email__, __credits__)
        fill_maps, fill_cubes = [], []

    # --- Bands ---
    for _, map_entry in maps.iterrows():
        if map_entry["map_name"] in fill_maps:
            continue
        map_file = path.join(str(map_entry["map_dir"]), source + str(map_entry["map_ext"]))
        if not path.exists(map_file):
            print(f'{"[ERROR]":<10}', f'Map {map_entry["map_name"]} not found: {map_file}')
            continue
        perbeam = "/beam" in str(map_entry.get("map_unit", ""))
        this_int, _ = sample_at_res(
            map_file, samp_ra, samp_dec,
            target_res_as=target_res_as, target_hdr=ov_hdr,
            line_name=map_entry["map_name"], source=source,
            path_save_fits=data_dir, save_fits=save_fits, perbeam=perbeam,
        )
        this_data["MAP_" + map_entry["map_name"].upper()] = Column(
            this_int, unit=au.Unit(str(map_entry["map_unit"])), description=map_entry["map_desc"])

        if str(map_entry.get("map_uc", "")).strip():
            uc_file = path.join(str(map_entry["map_dir"]), source + str(map_entry["map_uc"]))
            if path.exists(uc_file):
                uc_int, _ = sample_at_res(
                    uc_file, samp_ra, samp_dec,
                    target_res_as=target_res_as, target_hdr=ov_hdr,
                    perbeam=perbeam, unc=True,
                )
                this_data["EMAP_" + map_entry["map_name"].upper()] = Column(
                    uc_int, unit=au.Unit(str(map_entry["map_unit"])),
                    description=f'Uncertainty: {map_entry["map_desc"]}')
        print(f'{"[INFO]":<10}', f'Map {map_entry["map_name"]} sampled.')

    # --- Cubes ---
    for _, cube in cubes.iterrows():
        if cube["line_name"] in fill_cubes:
            continue
        cube_file = path.join(str(cube["line_dir"]), source + str(cube["line_ext"]))
        if not path.exists(cube_file):
            print(f'{"[ERROR]":<10}', f'Cube {cube["line_name"]} not found: {cube_file}')
            continue
        this_spec, _ = sample_at_res(
            cube_file, samp_ra, samp_dec,
            target_res_as=target_res_as, target_hdr=ov_hdr,
            line_name=cube["line_name"], source=source,
            path_save_fits=data_dir, save_fits=save_fits,
        )
        this_data["SPEC_" + cube["line_name"].upper()] = Column(
            this_spec, unit=au.Unit(str(cube["line_unit"])), description=cube["line_desc"])

        # Optional 2D map
        map_ext = str(cube.get("map_ext", "")).strip()
        if map_ext and map_ext not in ("nan", ""):
            b2d_file = path.join(str(cube["line_dir"]), source + map_ext)
            if path.exists(b2d_file):
                b2d, _ = sample_at_res(b2d_file, samp_ra, samp_dec,
                                       target_res_as=target_res_as, target_hdr=ov_hdr)
                this_data["MAP_" + cube["line_name"].upper()] = Column(
                    b2d, unit=au.Unit(str(cube["line_unit"])), description=cube["line_desc"])

        # Optional 2D uncertainty
        map_uc = str(cube.get("map_uc", "")).strip()
        if map_uc and map_uc not in ("nan", ""):
            uc_file = path.join(str(cube["line_dir"]), source + map_uc)
            if path.exists(uc_file):
                uc, _ = sample_at_res(uc_file, samp_ra, samp_dec,
                                      target_res_as=target_res_as, target_hdr=ov_hdr, unc=True)
                this_data["EMAP_" + cube["line_name"].upper()] = Column(
                    uc, unit=au.Unit(str(cube["line_unit"])),
                    description=f'Uncertainty: {cube["line_desc"]}')

        print(f'{"[INFO]":<10}', f'Cube {cube["line_name"]} sampled.')

    # --- Mask ---
    if len(input_mask) > 0:
        use_fixed = meta.get("use_fixed_vel_mask", False)
        if use_fixed:
            mask_unit  = input_mask["mask_unit"].iloc[0]
            mask_start = float(input_mask["mask_start"].iloc[0]) * au.Unit(mask_unit)
            mask_end   = float(input_mask["mask_end"].iloc[0])   * au.Unit(mask_unit)
            unit_v     = ov_hdr.get("CUNIT3", "m/s")
            v0, dv, crpix = ov_hdr["CRVAL3"], ov_hdr["CDELT3"], ov_hdr["CRPIX3"]
            vaxis = (v0 + (np.arange(n_chan) - (crpix - 1)) * dv) * au.Unit(unit_v)
            vaxis = vaxis.to(au.Unit(mask_unit))
            spec_mask = np.zeros((n_pts, n_chan))
            spec_mask[:, (vaxis >= mask_start) & (vaxis <= mask_end)] = 1.0
        else:
            mask_file = path.join(
                str(input_mask["mask_dir"].iloc[0]),
                source + str(input_mask["mask_ext"].iloc[0]))
            if not path.exists(mask_file):
                print(f'{"[ERROR]":<10}', f'Mask file not found: {mask_file}')
                spec_mask = np.zeros((n_pts, n_chan))
            else:
                spec_mask, _ = sample_mask(mask_file, samp_ra, samp_dec, target_hdr=ov_hdr)

        tag = "SPEC_" + str(input_mask["mask_name"].iloc[0]).upper()
        this_data[tag] = Column(spec_mask, unit=au.dimensionless_unscaled,
                                description=str(input_mask["mask_desc"].iloc[0]))
        print(f'{"[INFO]":<10}', f'Mask sampled.')

    os.makedirs(meta.get("out_dir", "Output/"), exist_ok=True)
    this_data.write(fname, format="ascii.ecsv", overwrite=True)
    print(f'{"[INFO]":<10}', f'PyStructure written: {fname}')
    return fname


# ============================================================================
# Helpers
# ============================================================================

def _build_fname(source, meta):
    resolution = meta.get("resolution", "angular")
    target_res = meta.get("target_res", 27.0)
    out_dir    = meta.get("out_dir", "Output/")
    suffix = (str(int(target_res)) + "as" if resolution == "angular"
              else str(int(target_res)) + "pc" if resolution == "physical"
              else "native")
    date_str = date.today().strftime("%Y_%m_%d")
    return os.path.join(out_dir, f"{source}_data_struct_{suffix}_{date_str}.ecsv")


def _init_table(source, params, meta, samp_ra, samp_dec, ov_hdr, target_res_as,
                version, author, email, credits_):
    from datetime import date as _date
    this_data = Table()
    this_data.meta.update({
        "Name": "PyStructure", "Version": version, "Authors": author,
        "Contacts": email, "Credits": credits_,
        "User": meta.get("user", ""), "Comments": meta.get("comments", ""),
        "Date": _date.today().strftime("%Y_%m_%d"), "Source": source,
    })
    this_data["ra_deg"]  = Column(samp_ra,  unit=au.deg, description="Right ascension (J2000)")
    this_data["dec_deg"] = Column(samp_dec, unit=au.deg, description="Declination (J2000)")
    this_data.meta["dist_mpc"]   = params["dist_mpc"]   * au.Mpc
    this_data.meta["posang_deg"] = params["posang_deg"] * au.deg
    this_data.meta["incl_deg"]   = params["incl_deg"]   * au.deg
    this_data.meta["beam_as"]    = target_res_as         * au.arcsec

    unit_v = ov_hdr.get("CUNIT3", "m/s")
    this_data.meta["SPEC_VCHAN0"] = ov_hdr["CRVAL3"] * au.Unit(unit_v)
    this_data.meta["SPEC_DELTAV"] = ov_hdr["CDELT3"] * au.Unit(unit_v)
    this_data.meta["SPEC_CRPIX"]  = ov_hdr["CRPIX3"]
    this_data.meta["input_bands"] = ""
    this_data.meta["input_cubes"] = ""

    rgal_deg, theta_rad = deproject(
        samp_ra, samp_dec,
        [params["posang_deg"], params["incl_deg"], params["ra_ctr"], params["dec_ctr"]],
        vector=True,
    )
    dist_mpc = params["dist_mpc"]
    r25      = params["r25"]
    this_data["rgal_as"]   = Column(rgal_deg * 3600,                    unit=au.arcsec, description="Deprojected galactocentric radius")
    this_data["rgal_kpc"]  = Column(np.deg2rad(rgal_deg) * dist_mpc * 1e3, unit=au.kpc, description="Deprojected galactocentric radius")
    this_data["rgal_r25"]  = Column(rgal_deg / (r25 / 60.0),               description="Deprojected galactocentric radius (r25 units)")
    this_data["theta_rad"] = Column(theta_rad,                          unit=au.rad,    description="Deprojected polar angle")
    return this_data


def _fill_checker(fname, samp_ra, samp_dec, maps, cubes):
    this_data = Table.read(fname)
    diff = (abs(np.nansum(this_data["ra_deg"]  - samp_ra  * au.deg))
            + abs(np.nansum(this_data["dec_deg"] - samp_dec * au.deg)))
    if diff > 1e-12 * au.deg:
        raise ValueError("Existing file coordinates do not match. "
                         'Set structure_creation = "default" to overwrite.')
    fill_maps = [b for b in maps["map_name"] if f"BAND_{b.upper()}" in this_data.colnames]
    fill_cubes = [c for c in cubes["line_name"] if f"MOM0_{c.upper()}" in this_data.colnames]
    return this_data, fill_maps, fill_cubes
