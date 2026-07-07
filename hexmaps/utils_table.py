"""
utils_table.py — Astropy Table I/O helpers, spectral shuffle, and moment maps.

All functions are self-contained and usable independently of the pipeline.

Contents
--------
I/O helpers
    load_hexmaps         — load a .ecsv database
    save_hexmaps         — save a Table to .ecsv
    find_latest_hexmaps  — find the most recently dated .ecsv for a target
    get_column_names         — return column names without loading all data
    get_spec_lines           — return spectral line names
    get_map_names            — return 2D map names

Spectral shuffle
    shuffle          — remap spectra onto a new (shifted) velocity axis

Moment computation
    get_mom_maps     — compute moment-0/1/2, Tpeak, rms, and EW maps
"""

import os
import copy
import glob
import numpy as np
from pathlib import Path
from astropy import units as u
from astropy.table import Table

from hexmaps.logger import get_logger

# Utility functions do not have their own pipeline "stage" — log messages
# should appear under whichever stage is calling them. Each function below
# accepts an optional `log` parameter (a StageLogger from get_logger); if not
# provided, falls back to the "Loading" stage (the typical context for
# standalone/analysis use of these functions).
_DEFAULT_LOG = get_logger("Loading")


# ============================================================================
# I/O helpers
# ============================================================================


def load_hexmaps(fname: str, log=None) -> Table:
    """
    Load a HexMaps .ecsv file into an Astropy Table.

    Parameters
    ----------
    fname : str or Path — path to the .ecsv file
    log : StageLogger, optional
        Logger to use for the error message (from get_logger()). Defaults to
        the "Loading" stage if not provided.

    Returns
    -------
    table : astropy.table.Table

    Raises
    ------
    FileNotFoundError if the file does not exist.
    """
    log = log or _DEFAULT_LOG
    fname = Path(fname)
    if not fname.exists():
        log.error(f"HexMaps file not found: {fname}")
        raise FileNotFoundError(f"HexMaps file not found: {fname}")
    return Table.read(fname)


def save_hexmaps(table: Table, fname: str, overwrite: bool = True) -> None:
    """
    Save an Astropy Table to a HexMaps .ecsv file.

    Creates the parent directory if it does not exist.

    Parameters
    ----------
    table     : astropy.table.Table
    fname     : str or Path
    overwrite : bool — if False, raise if *fname* already exists
    """
    fname = Path(fname)
    os.makedirs(fname.parent, exist_ok=True)
    table.write(str(fname), format="ascii.ecsv", overwrite=overwrite)


def find_latest_hexmaps(out_dir: str, target: str, log=None) -> str:
    """
    Find the most recently dated HexMaps .ecsv file for *target*.

    Files are matched by the glob pattern
    ``{out_dir}/{target}_hexmaps_*.ecsv`` and sorted lexicographically
    (which is equivalent to date-order for the YYYY_MM_DD filename convention).

    Parameters
    ----------
    out_dir : str — directory to search
    target  : str — target name
    log : StageLogger, optional
        Logger to use for the error message (from get_logger()). Defaults to
        the "Loading" stage if not provided.

    Returns
    -------
    path : str — path to the most recent matching file

    Raises
    ------
    FileNotFoundError if no matching file is found.
    """
    log = log or _DEFAULT_LOG
    pattern = os.path.join(out_dir, f"{target}_hexmaps_*.ecsv")
    matches = sorted(glob.glob(pattern))
    if not matches:
        log.error(f"No HexMaps file found for " f"'{target}' in '{out_dir}'")
        raise FileNotFoundError(
            f"No HexMaps file found for " f"'{target}' in '{out_dir}'"
        )
    return matches[-1]


def get_column_names(fname: str) -> list:
    """Return the column names of a HexMaps file."""
    return Table.read(fname).colnames


def get_spec_lines(fname: str) -> list:
    """Return the spectral line names stored in a HexMaps file (from SPEC_ columns)."""
    return [c[5:] for c in get_column_names(fname) if c.startswith("SPEC_")]


def get_map_names(fname: str) -> list:
    """Return the 2D map names stored in a HexMaps file (from MAP_ columns)."""
    return [c[4:] for c in get_column_names(fname) if c.startswith("MAP_")]


# ============================================================================
# Spectral shuffle
# ============================================================================


def shuffle(
    spec,
    vaxis,
    zero=None,
    new_vaxis=None,
    new_naxis=None,
    new_crval=None,
    new_crpix=None,
    new_cdelt=None,
    interp=None,
    missing=None,
    quiet=False,
):
    """
    Remap a spectrum (or array of spectra) onto a new velocity axis.

    The "shuffle" operation shifts each spectrum so that a reference velocity
    (``zero``) maps to v = 0 on the output axis.  This is used to stack spectra
    from different positions that have different systemic velocities, improving
    the sensitivity to faint emission.

    Port of IDL shuffle (cpropstoo, A. Leroy) by J. den Brok (2019).

    Parameters
    ----------
    spec      : np.ndarray — 1-D (single spectrum), 2-D (n_pts × n_chan),
                             or 3-D (nx × ny × n_chan) array
    vaxis     : array-like — original velocity axis (same units as zero/new_vaxis)
    zero      : scalar or array — velocity shift to apply per spectrum.
                Scalar: same shift for all spectra.
                Array: one shift per row (2-D input) or one per pixel (3-D).
    new_vaxis : array-like, optional — explicit output velocity axis.
                If not given, it is constructed from new_crval, new_crpix,
                new_cdelt, new_naxis (all defaulting to the input axis).
    interp    : int — 0 = nearest-neighbour (preserves noise statistics),
                       1 = linear (default, smoother but correlates noise)
    missing   : float — fill value for channels outside the valid range
                        (default: NaN)

    Returns
    -------
    output : np.ndarray — shuffled spectrum/array, same leading shape as *spec*
             but with n_chan = len(new_vaxis) along the last axis.

    Notes
    -----
    If new_vaxis is identical to vaxis and zero is 0, the function returns
    spec unchanged without any resampling.
    """
    # Build output velocity axis if not provided
    if new_vaxis is None:
        if new_cdelt is None:
            new_cdelt = vaxis[1] - vaxis[0]
        if new_crval is None or new_crpix is None:
            new_crval, new_crpix = vaxis[0], 1
        if new_naxis is None:
            new_naxis = len(vaxis)
        new_vaxis = (np.arange(new_naxis) - (new_crpix - 1.0)) * new_cdelt + new_crval

    # No-op check: same axis and no shift
    if len(new_vaxis) == len(vaxis) and np.sum(new_vaxis != vaxis) == 0:
        return spec

    n_chan = len(new_vaxis)
    dim_spec = np.shape(spec)

    if len(dim_spec) == 2:
        n_spec = dim_spec[0]
    elif len(dim_spec) == 3:
        n_spec = dim_spec[1] * dim_spec[2]
    else:
        n_spec = 1

    if zero is None:
        zero = 0.0
    if missing is None:
        missing = np.nan
    if interp is None:
        interp = 1

    orig_nchan = len(vaxis)
    orig_chan = np.arange(orig_nchan)
    new_nchan = len(new_vaxis)
    orig_deltav = vaxis[1] - vaxis[0]
    new_deltav = new_vaxis[1] - new_vaxis[0]

    # Pre-allocate output array
    if len(dim_spec) == 1:
        output = np.full(n_chan, missing, dtype=float)
    elif len(dim_spec) == 2:
        output = np.full((dim_spec[0], n_chan), missing, dtype=float)
    else:
        output = np.full((dim_spec[0], dim_spec[1], n_chan), missing, dtype=float)

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

        # Shift the original velocity axis by -this_zero so that emission
        # at this_zero maps to v=0 on the output axis
        this_vaxis = vaxis - this_zero

        # Ensure both axes are monotonically increasing
        if orig_deltav < 0 and (this_vaxis[1] - this_vaxis[0]) < 0:
            this_vaxis = np.flip(this_vaxis)
            this_spec = np.flip(this_spec)
        if new_deltav < 0 and (new_vaxis[1] - new_vaxis[0]) < 0:
            new_vaxis = np.flip(new_vaxis)

        # Map new velocity axis positions to channel indices in the shifted original axis
        channel_mapping = np.interp(new_vaxis, this_vaxis, orig_chan)
        overlap = np.where(
            (channel_mapping > 0.0) & (channel_mapping < orig_nchan - 1)
        )[0]
        if len(overlap) == 0:
            continue  # no valid overlap; output channels stay as missing

        new_spec = np.full(new_nchan, missing, dtype=float)
        if interp == 0:
            # Nearest-neighbour: integer channel lookup
            new_spec[overlap] = this_spec[
                np.array(np.rint(channel_mapping[overlap]), dtype=int)
            ]
        else:
            # Linear interpolation
            new_spec[overlap] = np.interp(new_vaxis[overlap], this_vaxis, this_spec)

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


def get_mom_maps(spec_cube, mask, vaxis, mom_calc=(3, 3, "fwhm"), noise_mask=None):
    """
    Compute integrated spectral properties from a masked spectral cube.

    For each sampling point, computes a set of quantities that characterise the
    emission line: integrated intensity (mom0), mean velocity (mom1), velocity
    dispersion (mom2), peak brightness (Tpeak), noise rms, and equivalent width
    (EW).  Uncertainties are propagated analytically.

    Moment definitions
    ------------------
    mom0 = ∑ T_i × dv  (integrated intensity, summed over masked channels)
    mom1 = ∑ T_i × v_i / ∑ T_i  (intensity-weighted mean velocity)
    mom2 = sqrt(∑ T_i × (v_i - mom1)² / ∑ T_i)  [math definition]
         → × sqrt(8 ln2)  to give FWHM  [mom2_method = "fwhm"]
    EW   = mom0 / Tpeak / sqrt(2π)  (equivalent width under a Gaussian profile)

    Moments 1 and 2 are computed using a high-S/N submask (pixels above
    SNthresh × rms with ≥ 3 consecutive channels) to reduce bias from low-S/N
    wings.

    Port of mom_computer.py (J. den Brok / L. Neumann).

    Parameters
    ----------
    spec_cube : astropy Quantity (n_pts × n_chan) — brightness temperature cube
    mask      : array-like (n_pts × n_chan)       — 0/1 integration mask
    vaxis     : astropy Quantity (n_chan,)          — velocity axis
    mom_calc  : tuple (SN_thresh, conseq_channels, mom2_method)
        SN_thresh       : float — S/N threshold for high-S/N submask
        conseq_channels : int   — min consecutive channels for submask
        mom2_method     : str   — "fwhm" | "sqrt" | "math"
    noise_mask : array-like (n_pts × n_chan), optional — 0/1 mask selecting
        the channels to use for noise (RMS) estimation. When supplied,
        the RMS is computed from channels where noise_mask == 1 (instead of
        the default channels-outside-the-integration-mask approach). Multiple
        velocity windows can be combined by OR-ing them into one mask before
        passing here.

    Returns
    -------
    dict mapping str → astropy Quantity (n_pts,):
        rms, tpeak, mom0, mom0_err, mom1, mom1_err, mom2, mom2_err, ew, ew_err
    """
    spec_vals = spec_cube.value
    v_vals = vaxis.value
    dv = abs(v_vals[0] - v_vals[1])
    spec_unit = spec_cube.unit
    v_unit = vaxis.unit

    SNthresh = mom_calc[0]
    conseq_channels = int(max(float(mom_calc[1]), 3))
    mom2_method = mom_calc[2]
    fac_mom2 = np.sqrt(8 * np.log(2)) if mom2_method == "fwhm" else 1.0

    n_pts = spec_vals.shape[0]
    mom2_unit = v_unit if mom2_method == "fwhm" else v_unit**2

    # Initialise all output arrays with NaN
    mom_maps = {
        "rms": np.full(n_pts, np.nan) * spec_unit,
        "tpeak": np.full(n_pts, np.nan) * spec_unit,
        "mom0": np.full(n_pts, np.nan) * spec_unit * v_unit,
        "mom0_err": np.full(n_pts, np.nan) * spec_unit * v_unit,
        "mom1": np.full(n_pts, np.nan) * v_unit,
        "mom1_err": np.full(n_pts, np.nan) * v_unit,
        "mom2": np.full(n_pts, np.nan) * mom2_unit,
        "mom2_err": np.full(n_pts, np.nan) * mom2_unit,
        "ew": np.full(n_pts, np.nan) * v_unit,
        "ew_err": np.full(n_pts, np.nan) * v_unit,
    }

    for m in range(n_pts):
        spectrum = spec_vals[m, :]
        mask_m = np.array(mask[m, :], dtype=float)

        # Skip points with no valid data
        if np.nansum(spectrum != 0) < 1:
            continue

        # RMS noise: use explicit noise channels if provided, otherwise
        # use channels outside the integration mask.
        if noise_mask is not None:
            noise_chans = np.asarray(noise_mask)[m].astype(bool)
            rms_vals = spectrum[noise_chans & (spectrum != 0)]
        else:
            rms_vals = spectrum[np.logical_and(mask_m == 0, spectrum != 0)]
        rms = np.nanstd(rms_vals)
        mom_maps["rms"][m] = rms * spec_unit

        # Peak brightness within the mask
        tpeak = np.nanmax(spectrum * mask_m)
        mom_maps["tpeak"][m] = tpeak * spec_unit

        # Moment 0: integrated intensity
        mom0 = np.nansum(spectrum * mask_m) * dv
        mom_maps["mom0"][m] = mom0 * spec_unit * v_unit
        # Uncertainty: noise per channel × sqrt(N_mask channels) × dv
        mom_maps["mom0_err"][m] = (
            np.sqrt(np.nansum(mask_m)) * rms * dv * spec_unit * v_unit
        )

        # High-S/N submask for moments 1 and 2
        # Requires SNthresh × rms AND ≥ 3 consecutive channels above threshold
        hsmask = (spectrum * mask_m > SNthresh * rms).astype(int)
        hsmask = ((hsmask + np.roll(hsmask, 1) + np.roll(hsmask, -1)) >= 3).astype(int)
        if np.nansum(hsmask) < conseq_channels - 2:
            continue  # insufficient high-S/N channels; skip moments 1 and 2
        # Dilate the high-S/N mask to include wings
        for _ in range(5):
            hsmask = ((hsmask + np.roll(hsmask, 1) + np.roll(hsmask, -1)) >= 1).astype(
                int
            )

        den1 = np.nansum(spectrum * hsmask)

        # Moment 1: intensity-weighted mean velocity
        mom1 = np.nansum(spectrum * v_vals * hsmask) / den1
        mom_maps["mom1"][m] = mom1 * v_unit
        numer = rms**2 * np.nansum(hsmask * (v_vals - mom1) ** 2)
        mom_maps["mom1_err"][m] = np.sqrt(numer / den1**2) * v_unit

        # Moment 2: velocity dispersion
        mom2_math = np.nansum(spectrum * hsmask * (v_vals - mom1) ** 2) / den1
        numer = rms**2 * np.nansum((hsmask * (v_vals - mom1) ** 2 - mom2_math) ** 2)
        mom2_err = np.sqrt(numer / den1**2)
        if mom2_method == "fwhm":
            mom_maps["mom2"][m] = fac_mom2 * np.sqrt(mom2_math) * v_unit
            mom_maps["mom2_err"][m] = (
                fac_mom2 * mom2_err / (2 * np.sqrt(mom2_math)) * v_unit
            )
        else:
            mom_maps["mom2"][m] = mom2_math * v_unit**2
            mom_maps["mom2_err"][m] = mom2_err * v_unit**2

        # Equivalent width
        ew = np.nansum(spectrum * hsmask) * dv / tpeak / np.sqrt(2 * np.pi)
        mom_maps["ew"][m] = ew * v_unit
        term1 = rms**2 * np.nansum(hsmask) * dv**2 / (2 * np.pi * tpeak**2)
        term2 = ew**2 - ew * dv / np.sqrt(2 * np.pi)
        mom_maps["ew_err"][m] = np.sqrt(term1 + term2) * v_unit

    return mom_maps


def build_noise_mask(noise_mask_df, vaxis, shape):
    """
    Build a noise-channel mask from a DataFrame of velocity windows.

    Combines all rows of *noise_mask_df* (each defining a velocity range
    [mask_start, mask_end] in *mask_unit*) into a single boolean mask by
    OR-ing them together, allowing multiple line-free windows to be specified
    in config and merged automatically.

    Parameters
    ----------
    noise_mask_df : pd.DataFrame with columns mask_start, mask_end, mask_unit
    vaxis         : astropy Quantity (n_chan,) — velocity axis in any unit
    shape         : tuple — output shape, either

                    * (n_pts, n_chan) for the hex-grid path, where the channel
                      axis is last and the mask is broadcast over all points, or
                    * (n_chan, ny, nx) for the PPV path, where the channel axis
                      is first and the mask is broadcast over all pixels.

    Returns
    -------
    noise_mask : np.ndarray of bool, same shape as *shape* — True in channels
                 that should be used for noise (RMS) estimation, or None if
                 no channels fell within any window (with a warning logged).
    """
    from astropy import units as _au

    LOG = get_logger("Loading")
    n_chan = shape[-1] if len(shape) == 2 else shape[0]
    chan_mask = np.zeros(n_chan, dtype=bool)

    for _, row in noise_mask_df.iterrows():
        try:
            mask_unit = str(row["mask_unit"]).strip()
            mask_start = float(row["mask_start"]) * _au.Unit(mask_unit)
            mask_end = float(row["mask_end"]) * _au.Unit(mask_unit)
            vaxis_conv = vaxis.to(_au.Unit(mask_unit))
            window = (vaxis_conv >= mask_start) & (vaxis_conv <= mask_end)
            chan_mask |= (
                window.value if hasattr(window, "value") else np.asarray(window)
            )
        except Exception as e:
            LOG.warning(
                f"Could not parse noise velocity window row " f"{row.to_dict()} — {e}"
            )

    if not chan_mask.any():
        LOG.warning(
            "use_fixed_noise_mask is True but no channels fall within the "
            "specified noise velocity windows. Falling back to mask-inverted noise."
        )
        return None

    # Broadcast to the requested output shape
    if len(shape) == 2:
        # hex-grid: (n_pts, n_chan) — repeat channel mask for all points
        return np.broadcast_to(chan_mask[None, :], shape).copy()
    else:
        # PPV: (n_chan, ny, nx) — channel axis is first
        return np.broadcast_to(chan_mask[:, None, None], shape).copy()


def parse_ref_line(ref_line_method, line_names):
    """
    Parse the ``ref_line`` config value into a structured specification.

    The value is a comma-separated list of tokens in three categories:

    **Combinator token** (optional, default ``OR``):
        ``AND`` — all masks must be True.
        ``OR``  — any mask being True is sufficient (default).

    **External-mask tokens** (zero or more):
        ``input``  — include the external FITS input mask.
        ``window`` — include the fixed velocity-window mask.

    **Line-selection token** (exactly one, mandatory unless only
    ``input``/``window`` are given):
        ``first``        — S/N mask from the first line in the cube list.
        ``all``          — S/N mask from all lines in the cube list.
        ``<n>``          — S/N mask from the first *n* lines (integer ≥ 1).
        ``individual``   — one independent S/N mask per line, each applied
                           only to that line's moments.  External masks are
                           still combined with each per-line mask.
        ``<LINE_NAME>``  — one or more named lines from the cube list,
                           matched case-insensitively.

    Parameters
    ----------
    ref_line_method : str or int
    line_names : list of str  — cube line names from the config.

    Returns
    -------
    mask_lines     : list of str
        Line names to build S/N masks from (original case).
        Empty list means no S/N mask requested (only external masks).
        For ``individual`` mode this is the full cube list;
        ``use_individual`` tells the caller how to apply them.
    use_individual : bool
        True  → build one S/N mask per line, apply it only to that line.
        False → combine all masks into a single master mask for every line.
    use_input  : bool — include the external input mask.
    use_window : bool — include the fixed velocity-window mask.
    combinator : str  — ``"AND"`` or ``"OR"``.

    Raises
    ------
    ValueError — empty input, unknown token, conflicting selection modes,
                 or a named line not found in the cube list.

    Examples
    --------
    ``first``                 master mask from the first line
    ``12co21``                master mask from 12co21
    ``12co21, 12co10``        master mask: OR of both S/N masks
    ``all``                   master mask: OR of all S/N masks
    ``2``                     master mask: OR of first-2 S/N masks
    ``individual``            per-line masks, applied per-line
    ``individual, input, AND``  per-line mask AND input mask
    ``first, input``          master mask: OR of first-line S/N + input
    ``12co21, input, AND``    master mask: 12co21 S/N AND input
    ``first, window, AND``    master mask: first-line S/N AND window
    ``input``                 only input mask (no S/N mask)
    """
    if not line_names:
        raise ValueError("parse_ref_line: line_names is empty.")

    # --- tokenise -------------------------------------------------------
    raw_tokens = [t.strip() for t in str(ref_line_method).split(",") if t.strip()]
    if not raw_tokens:
        raise ValueError(
            "parse_ref_line: ref_line is empty. "
            "Provide at least one token (e.g. 'first', a line name, "
            "'input', or 'window')."
        )

    # --- classify -------------------------------------------------------
    upper_to_orig = {ln.upper(): ln for ln in line_names}

    combinator     = "OR"
    use_input      = False
    use_window     = False
    use_individual = False
    line_tokens    = []   # resolved line-names or keyword strings

    for raw in raw_tokens:

        if raw in ("AND", "OR"):
            combinator = raw
            continue

        if raw == "input":
            use_input = True
            continue

        if raw == "window":
            use_window = True
            continue

        if raw == "individual":
            use_individual = True
            line_tokens.append("individual")
            continue

        if raw in ("first", "all"):
            line_tokens.append(raw)
            continue

        # Positive integer?
        try:
            n = int(raw)
        except ValueError:
            n = None
        if n is not None:
            if n < 1:
                raise ValueError(
                    f"parse_ref_line: integer token must be \u2265 1, got '{raw}'."
                )
            line_tokens.append(raw)   # store as string for uniform handling
            continue

        # Named line from the cube list?
        if raw in upper_to_orig:
            line_tokens.append(upper_to_orig[raw])   # preserve original case
            continue

        # Nothing matched
        raise ValueError(
            f"parse_ref_line: unrecognised token '{raw}'. "
            f"Expected one of: AND, OR, input, window, first, all, individual, "
            f"a positive integer, or a line name from {line_names}."
        )

    # --- resolve line-selection tokens → mask_lines ---------------------
    keyword_found = [t for t in line_tokens if t in ("first", "all", "individual")]
    named_found   = [t for t in line_tokens if t in line_names]
    int_found     = []
    for t in line_tokens:
        try:
            int_found.append(int(t))
        except ValueError:
            pass

    n_modes = (
        (1 if use_individual else 0)
        + (1 if "all"   in keyword_found else 0)
        + (1 if "first" in keyword_found else 0)
        + (1 if int_found else 0)
        + (1 if named_found else 0)
    )
    if n_modes > 1:
        raise ValueError(
            f"parse_ref_line: conflicting line-selection tokens in "
            f"'{ref_line_method}'. Provide exactly one of: first, all, "
            f"individual, an integer, or line name(s)."
        )

    if not line_tokens:
        # Only external-mask tokens — no S/N mask
        if not use_input and not use_window:
            raise ValueError(
                f"parse_ref_line: no valid token found in '{ref_line_method}'. "
                f"Provide at least one of: first, all, individual, an integer, "
                f"a line name, 'input', or 'window'."
            )
        mask_lines = []

    elif use_individual:
        mask_lines = list(line_names)   # all lines get their own mask

    elif "all" in keyword_found:
        mask_lines = list(line_names)

    elif "first" in keyword_found:
        mask_lines = [line_names[0]]

    elif int_found:
        n = max(1, min(int_found[0], len(line_names)))
        mask_lines = list(line_names[:n])

    else:
        mask_lines = named_found   # original case preserved

    return mask_lines, use_individual, use_input, use_window, combinator