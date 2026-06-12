"""
pystructurePipeline.utils.table_utils
======================================
Astropy Table I/O helpers, spectral shuffle, and moment computation.

Contains (all inline, no legacy imports):
  - load_pystructure / save_pystructure / find_latest_pystructure
  - shuffle         — remap a spectrum onto a new velocity axis
  - get_mom_maps    — compute moment-0/1/2, Tpeak, rms, EW
"""

import os
import copy
import glob
import numpy as np
from pathlib import Path
from astropy import units as au
from astropy.table import Table


# ============================================================================
# I/O helpers
# ============================================================================

def load_pystructure(fname: str) -> Table:
    """Load a PyStructure .ecsv file."""
    fname = Path(fname)
    if not fname.exists():
        raise FileNotFoundError(f"PyStructure file not found: {fname}")
    return Table.read(fname)


def save_pystructure(table: Table, fname: str, overwrite: bool = True) -> None:
    """Save a PyStructure Astropy Table to an .ecsv file."""
    fname = Path(fname)
    os.makedirs(fname.parent, exist_ok=True)
    table.write(str(fname), format="ascii.ecsv", overwrite=overwrite)


def find_latest_pystructure(out_dir: str, source: str) -> str:
    """Return the most recently dated PyStructure file for *source*."""
    pattern = os.path.join(out_dir, f"{source}_data_struct_*.ecsv")
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No PyStructure file found for '{source}' in '{out_dir}'"
        )
    return matches[-1]


def get_column_names(fname: str) -> list:
    return Table.read(fname).colnames


def get_spec_lines(fname: str) -> list:
    return [c[5:] for c in get_column_names(fname) if c.startswith("SPEC_")]


def get_map_names(fname: str) -> list:
    return [c[5:] for c in get_column_names(fname) if c.startswith("MAP_")]


# ============================================================================
# Spectral shuffle
# ============================================================================

def shuffle(spec, vaxis, zero=None, new_vaxis=None,
            new_naxis=None, new_crval=None, new_crpix=None, new_cdelt=None,
            interp=None, missing=None, quiet=False):
    """
    Remap a spectrum (or array of spectra) onto a new velocity axis.

    Port of IDL shuffle (cpropstoo, J. den Brok 2019).

    Parameters
    ----------
    spec     : 1-D, 2-D (n_spec × n_chan), or 3-D array
    vaxis    : original velocity axis
    zero     : scalar or array — velocity shift applied per spectrum
    new_vaxis: target velocity axis; or build from new_crval/crpix/cdelt/naxis
    interp   : 0 = nearest-neighbour, 1 = linear (default)
    missing  : fill value for out-of-range channels (default NaN)
    """
    if new_vaxis is None:
        if new_cdelt is None:
            new_cdelt = vaxis[1] - vaxis[0]
        if new_crval is None or new_crpix is None:
            new_crval, new_crpix = vaxis[0], 1
        if new_naxis is None:
            new_naxis = len(vaxis)
        new_vaxis = (np.arange(new_naxis) - (new_crpix - 1.0)) * new_cdelt + new_crval

    if len(new_vaxis) == len(vaxis) and np.sum(new_vaxis != vaxis) == 0:
        return spec

    n_chan   = len(new_vaxis)
    dim_spec = np.shape(spec)

    if len(dim_spec) == 2:
        shape, n_spec = "ARRAY", dim_spec[0]
    elif len(dim_spec) == 3:
        shape, n_spec = "CUBE", dim_spec[1] * dim_spec[2]
    else:
        shape, n_spec = "SPEC", 1

    if zero is None:
        zero = 0.0
    if missing is None:
        missing = np.nan
    if interp is None:
        interp = 1

    orig_nchan  = len(vaxis)
    orig_chan   = np.arange(orig_nchan)
    new_nchan   = len(new_vaxis)
    orig_deltav = vaxis[1] - vaxis[0]
    new_deltav  = new_vaxis[1] - new_vaxis[0]

    if len(dim_spec) == 1:
        output = np.full(n_chan, missing, dtype=float)
    elif len(dim_spec) == 2:
        output = np.full((dim_spec[0], n_chan), missing, dtype=float)
    else:
        output = np.full((dim_spec[0], dim_spec[1], n_chan), missing, dtype=float)

    no_overlap_ct = 0
    for ii in range(n_spec):
        if len(dim_spec) == 3:
            yy = ii // dim_spec[0]
            xx = ii % dim_spec[0]
            this_spec = copy.copy(spec[xx, yy, :])
            this_zero = zero[xx, yy] if hasattr(zero, "__len__") else zero
        elif len(dim_spec) == 2:
            this_spec = copy.copy(spec[ii, :])
            this_zero = zero[ii] if hasattr(zero, "__len__") else zero
        else:
            this_spec = copy.copy(spec)
            this_zero = zero

        this_vaxis = vaxis - this_zero

        if orig_deltav < 0 and (this_vaxis[1] - this_vaxis[0]) < 0:
            this_vaxis = np.flip(this_vaxis)
            this_spec  = np.flip(this_spec)
        if new_deltav < 0 and (new_vaxis[1] - new_vaxis[0]) < 0:
            new_vaxis  = np.flip(new_vaxis)

        channel_mapping = np.interp(new_vaxis, this_vaxis, orig_chan)
        overlap = np.where(
            (channel_mapping > 0.0) & (channel_mapping < orig_nchan - 1)
        )[0]
        if len(overlap) == 0:
            no_overlap_ct += 1
            continue

        new_spec = np.full(new_nchan, missing, dtype=float)
        if interp == 0:
            new_spec[overlap] = this_spec[
                np.array(np.rint(channel_mapping[overlap]), dtype=int)]
        else:
            new_spec[overlap] = np.interp(
                new_vaxis[overlap], this_vaxis, this_spec)

        if new_deltav < 0:
            new_spec = np.flip(new_spec)

        if len(dim_spec) == 3:
            output[xx, yy, :] = new_spec
        elif len(dim_spec) == 2:
            output[ii, :] = new_spec
        else:
            output = new_spec

    return output


# ============================================================================
# Moment computation
# ============================================================================

def get_mom_maps(spec_cube, mask, vaxis, mom_calc=(3, 3, "fwhm")):
    """
    Compute moment maps from a masked spectral cube.

    Port of mom_computer.py (J. den Brok / L. Neumann).

    Parameters
    ----------
    spec_cube : astropy Quantity (n_pts × n_chan)
    mask      : array-like (n_pts × n_chan)
    vaxis     : astropy Quantity (n_chan,)
    mom_calc  : [SN_thresh, conseq_channels, mom2_method]
                mom2_method ∈ {"fwhm", "sqrt", "math"}

    Returns
    -------
    dict of astropy Quantities keyed by:
      rms, tpeak, mom0, mom0_err, mom1, mom1_err, mom2, mom2_err, ew, ew_err
    """
    spec_vals  = spec_cube.value
    v_vals     = vaxis.value
    dv         = abs(v_vals[0] - v_vals[1])
    spec_unit  = spec_cube.unit
    v_unit     = vaxis.unit

    SNthresh        = mom_calc[0]
    conseq_channels = int(max(float(mom_calc[1]), 3))
    mom2_method     = mom_calc[2]
    fac_mom2        = np.sqrt(8 * np.log(2)) if mom2_method == "fwhm" else 1.0

    n_pts = spec_vals.shape[0]
    mom2_unit = v_unit if mom2_method == "fwhm" else v_unit ** 2

    mom_maps = {
        "rms":      np.full(n_pts, np.nan) * spec_unit,
        "tpeak":    np.full(n_pts, np.nan) * spec_unit,
        "mom0":     np.full(n_pts, np.nan) * spec_unit * v_unit,
        "mom0_err": np.full(n_pts, np.nan) * spec_unit * v_unit,
        "mom1":     np.full(n_pts, np.nan) * v_unit,
        "mom1_err": np.full(n_pts, np.nan) * v_unit,
        "mom2":     np.full(n_pts, np.nan) * mom2_unit,
        "mom2_err": np.full(n_pts, np.nan) * mom2_unit,
        "ew":       np.full(n_pts, np.nan) * v_unit,
        "ew_err":   np.full(n_pts, np.nan) * v_unit,
    }

    for m in range(n_pts):
        spectrum = spec_vals[m, :]
        mask_m   = np.array(mask[m, :], dtype=float)

        if np.nansum(spectrum != 0) < 1:
            continue

        # RMS
        rms = np.nanstd(spectrum[np.logical_and(mask_m == 0, spectrum != 0)])
        mom_maps["rms"][m] = rms * spec_unit

        # Tpeak
        tpeak = np.nanmax(spectrum * mask_m)
        mom_maps["tpeak"][m] = tpeak * spec_unit

        # Mom0
        mom0 = np.nansum(spectrum * mask_m) * dv
        mom_maps["mom0"][m]     = mom0     * spec_unit * v_unit
        mom_maps["mom0_err"][m] = np.sqrt(np.nansum(mask_m)) * rms * dv * spec_unit * v_unit

        # High-S/N mask for moments 1 and 2
        hsmask = (spectrum * mask_m > SNthresh * rms).astype(int)
        hsmask = ((hsmask + np.roll(hsmask, 1) + np.roll(hsmask, -1)) >= 3).astype(int)
        if np.nansum(hsmask) < conseq_channels - 2:
            continue
        for _ in range(5):
            hsmask = ((hsmask + np.roll(hsmask, 1) + np.roll(hsmask, -1)) >= 1).astype(int)

        den1 = np.nansum(spectrum * hsmask)

        # Mom1
        mom1 = np.nansum(spectrum * v_vals * hsmask) / den1
        mom_maps["mom1"][m] = mom1 * v_unit
        numer = rms ** 2 * np.nansum(hsmask * (v_vals - mom1) ** 2)
        mom_maps["mom1_err"][m] = np.sqrt(numer / den1 ** 2) * v_unit

        # Mom2
        mom2_math = np.nansum(spectrum * hsmask * (v_vals - mom1) ** 2) / den1
        numer     = rms ** 2 * np.nansum((hsmask * (v_vals - mom1) ** 2 - mom2_math) ** 2)
        mom2_err  = np.sqrt(numer / den1 ** 2)
        if mom2_method == "fwhm":
            mom_maps["mom2"][m]     = fac_mom2 * np.sqrt(mom2_math) * v_unit
            mom_maps["mom2_err"][m] = fac_mom2 * mom2_err / (2 * np.sqrt(mom2_math)) * v_unit
        else:
            mom_maps["mom2"][m]     = mom2_math * v_unit ** 2
            mom_maps["mom2_err"][m] = mom2_err  * v_unit ** 2

        # EW
        ew = np.nansum(spectrum * hsmask) * dv / tpeak / np.sqrt(2 * np.pi)
        mom_maps["ew"][m] = ew * v_unit
        term1 = rms ** 2 * np.nansum(hsmask) * dv ** 2 / (2 * np.pi * tpeak ** 2)
        term2 = ew ** 2 - ew * dv / np.sqrt(2 * np.pi)
        mom_maps["ew_err"][m] = np.sqrt(term1 + term2) * v_unit

    return mom_maps
