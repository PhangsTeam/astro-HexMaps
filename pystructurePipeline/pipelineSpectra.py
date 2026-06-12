"""
stage_spectra: spectral processing, masking, and moment computation.

All logic is contained here — no imports from the legacy scripts/ directory.
"""

import numpy as np
import pandas as pd
from astropy import units as au
from astropy.stats import median_absolute_deviation
from astropy.table import Table, Column

from pystructurePipeline.utilsTable import shuffle, get_mom_maps


# ============================================================================
# Mask construction
# ============================================================================

def construct_mask(ref_line, this_data, SN_processing):
    """
    Build a 2-level S/N mask from the reference spectral line.

    Port of processing_spec.construct_mask (J. den Brok / L. Neumann).

    Returns
    -------
    mask         : Column (n_pts × n_chan)
    line_vmean   : array of mean velocities (km/s) per sampling point
    line_vaxis   : velocity axis (km/s)
    """
    ref_line_data = this_data["SPEC_" + ref_line]
    n_pts  = np.shape(ref_line_data)[0]
    n_chan = np.shape(ref_line_data)[1]

    line_vaxis = (this_data.meta["SPEC_VCHAN0"]
                  + (np.arange(n_chan) - (this_data.meta["SPEC_CRPIX"] - 1))
                  * this_data.meta["SPEC_DELTAV"])
    line_vaxis = line_vaxis.to(au.km / au.s)

    # Per-spectrum RMS using 2-pass MAD
    rms = median_absolute_deviation(ref_line_data, axis=None, ignore_nan=True)
    rms = median_absolute_deviation(
        ref_line_data[np.where(ref_line_data < 3 * rms)], ignore_nan=True)

    mask_rough  = ref_line_data < 3 * rms
    masked_cube = np.where(mask_rough, ref_line_data, np.nan)
    med_mask    = np.nanmedian(masked_cube, axis=1)
    mad_mask    = np.nanmedian(np.abs(masked_cube - med_mask[:, None]), axis=1)

    low_thresh  = SN_processing[0] * mad_mask[:, None]
    high_thresh = SN_processing[1] * mad_mask[:, None]

    mask     = (ref_line_data > high_thresh).astype(int)
    low_mask = (ref_line_data > low_thresh).astype(int)

    # Require adjacent channel agreement
    mask = mask & (np.roll(mask, 1, 1) | np.roll(mask, -1, 1))

    # Remove spectral spikes
    mask     = ((mask     + np.roll(mask,     1, 1) + np.roll(mask,     -1, 1)) >= 3).astype(int)
    low_mask = ((low_mask + np.roll(low_mask, 1, 1) + np.roll(low_mask, -1, 1)) >= 3).astype(int)

    # Grow high-S/N core to low-S/N wings
    for _ in range(5):
        mask = (((mask + np.roll(mask, 1, 1) + np.roll(mask, -1, 1)) >= 1).astype(int) * low_mask)

    # Grow mask edge by 2 channels
    for _ in range(2):
        mask = ((mask + np.roll(mask, 1, 1) + np.roll(mask, -1, 1)) >= 1).astype(int)

    # Mean velocity per point
    mask_q = mask * au.dimensionless_unscaled
    line_vmean = np.zeros(n_pts) * np.nan * au.km / au.s
    for jj in range(n_pts):
        denom = np.nansum(ref_line_data[jj, :] * mask_q[jj, :])
        if denom != 0:
            line_vmean[jj] = (np.nansum(line_vaxis * ref_line_data[jj, :] * mask_q[jj, :])
                              / denom)

    return mask_q, line_vmean, line_vaxis


# ============================================================================
# Strict spatial mask
# ============================================================================

def _apply_strict_mask(mask, this_data):
    """Remove spatially isolated mask features (connected-component filter)."""
    ra, dec = this_data["ra_deg"], this_data["dec_deg"]
    n_chan  = np.shape(mask)[1]
    sep     = this_data.meta["beam_as"] / 3600 / 2

    for jj in range(n_chan):
        mask_spec  = mask[:, jj]
        mask_labels = np.zeros_like(mask_spec)
        label = 1

        for n in range(len(mask_labels)):
            if mask_labels[n] != 0:
                continue
            if mask_spec[n] == 0:
                mask_labels[n] = -99
                continue
            dist_array = np.sqrt((ra - ra[n]) ** 2 + (dec - dec[n]) ** 2)
            idx_neigh  = np.where(
                abs(dist_array - sep) < 0.1 * this_data.meta["beam_as"].to(au.deg))
            labels_given = np.unique(mask_labels[idx_neigh])
            index = labels_given[labels_given > 0]
            if len(index) > 0:
                mask_labels[n] = index[0]
                for i in range(len(index) - 1):
                    mask_labels[mask_labels == index[i + 1]] = index[0]
            else:
                mask_labels[n] = label
                label += 1

        for lab in np.unique(mask_labels):
            if lab <= 0:
                continue
            if len(mask[:, jj][mask_labels == lab]) < 5:
                mask[:, jj][mask_labels == lab] = 0

    return mask


# ============================================================================
# Hyperfine structure mask
# ============================================================================

def _build_hfs_mask(mask, line_name, hfs_data, this_data):
    """Shift the mask to cover hyperfine satellite lines."""
    lines_hfs = list(set(hfs_data["hfs_name"]))
    if line_name not in lines_hfs:
        return None

    idx_line  = lines_hfs.index(line_name)
    idx_cols  = hfs_data["hfs_name"] == line_name
    restfreqs = [f * au.Unit(str(u)) for f, u in
                 zip(hfs_data["hfs_ref_freq"][idx_cols], hfs_data["unit"][idx_cols])]
    hfs_freqs = [f * au.Unit(str(u)) for f, u in
                 zip(hfs_data["hfs_freq"][idx_cols], hfs_data["unit"][idx_cols])]

    v_ch     = this_data.meta["SPEC_DELTAV"].to(au.km / au.s)
    mask_hfs = np.copy(mask)

    for freq, restfreq in zip(hfs_freqs, restfreqs):
        v_shift  = freq.to(au.km / au.s, equivalencies=au.doppler_radio(restfreq))
        shift_ch = int(np.rint(v_shift.value / v_ch.value))

        mask_shift = np.zeros_like(mask, dtype=float)
        if shift_ch > 0:
            mask_shift[:, shift_ch:] = mask[:, :-shift_ch]
        elif shift_ch < 0:
            mask_shift[:, :shift_ch] = mask[:, -shift_ch:]
        else:
            mask_shift = mask.copy()

        mask_hfs[mask_shift == 1] = 1

    return mask_hfs * au.dimensionless_unscaled


# ============================================================================
# Stage entry point
# ============================================================================

def run_spectra(source, fname, meta, cubes, input_mask, hfs_data):
    """
    Process spectra for *source*: build mask, compute moments, shuffle.

    Reads the .ecsv written by stage_regrid, enriches it with moment
    columns and shuffled spectra, then overwrites the file.
    """
    use_input_mask     = meta.get("use_input_mask", False)
    use_fixed_vel_mask = meta.get("use_fixed_vel_mask", False)
    use_mask           = use_input_mask or use_fixed_vel_mask
    use_hfs_lines      = meta.get("use_hfs_lines", False)
    strict_mask        = meta.get("strict_mask", False)

    ref_line_method = meta.get("ref_line", "first")
    SN_processing   = meta.get("SN_processing", [2, 4])
    mom_calc        = [meta.get("mom_thresh", 5),
                       meta.get("conseq_channels", 3),
                       meta.get("mom2_method", "fwhm")]
    shuff_axis      = [meta.get("NAXIS_shuff", 200), meta.get("CDELT_SHUFF", 4000.0)]

    this_data = Table.read(fname)
    n_lines   = len(cubes["line_name"])
    line_names = [str(l) for l in cubes["line_name"]]

    # Resolve reference line
    if ref_line_method in line_names:
        ref_line = ref_line_method.upper()
    else:
        ref_line = line_names[0].upper()

    n_chan = np.shape(this_data["SPEC_" + ref_line])[1]

    # ------------------------------------------------------------------
    # Build / load mask
    # ------------------------------------------------------------------
    if use_mask:
        if len(input_mask) == 0:
            print(f'{"[ERROR]":<10}', 'use_mask is True but no mask defined in imaging_key.')
        mask_tag = f'SPEC_{str(input_mask["mask_name"].iloc[0]).upper()}'
        mask = this_data[mask_tag]
        del this_data[mask_tag]
        _, ref_line_vmean, ref_line_vaxis = construct_mask(ref_line, this_data, SN_processing)
    else:
        print(f'{"[INFO]":<10}', 'Building mask from prior line(s).')
        mask, ref_line_vmean, ref_line_vaxis = construct_mask(
            ref_line, this_data, SN_processing)
        this_data["SPEC_MASK_" + ref_line] = Column(
            mask, unit=au.dimensionless_unscaled,
            description=f"Velocity-integration mask for {ref_line}")

        # Multi-line mask
        if ref_line_method == "first":
            n_mask = 0
        elif ref_line_method == "all":
            n_mask = n_lines
        elif isinstance(ref_line_method, int):
            n_mask = min(n_lines, ref_line_method)
        else:
            n_mask = 0

        for n_mask_i in range(1, n_mask + 1):
            line_i = line_names[n_mask_i].upper()
            mask_i, _, _ = construct_mask(line_i, this_data, SN_processing)
            this_data["SPEC_MASK_" + line_i] = Column(
                mask_i, unit=au.dimensionless_unscaled,
                description=f"Velocity-integration mask for {line_i}")
            mask = (mask.value.astype(int) | mask_i.value.astype(int)) * au.dimensionless_unscaled

        # HI combined mask
        if ref_line_method == "ref+HI":
            if "hi" in line_names:
                mask_hi, vmean_hi, _ = construct_mask("HI", this_data, SN_processing)
                mask = (mask.value.astype(int) | mask_hi.value.astype(int)) * au.dimensionless_unscaled
                rgal  = this_data["rgal_r25"]
                n_pts = len(rgal)
                vmean_comb = np.zeros(n_pts) * np.nan
                for jj in range(n_pts):
                    vmean_comb[jj] = (ref_line_vmean[jj].value if rgal[jj] < 0.23
                                      else vmean_hi[jj].value)
                ref_line_vmean = vmean_comb
            else:
                print(f'{"[WARNING]":<10}', 'HI not in PyStructure; skipping HI mask.')

        if strict_mask:
            mask_arr = _apply_strict_mask(mask.value.astype(int), this_data)
            mask = mask_arr * au.dimensionless_unscaled

    # Determine HFS lines
    lines_hfs = list(set(hfs_data["hfs_name"])) if (use_hfs_lines and hfs_data is not None) else []

    # ------------------------------------------------------------------
    # HFS masks
    # ------------------------------------------------------------------
    if use_hfs_lines and hfs_data is not None:
        for jj in range(n_lines):
            line_name = line_names[jj]
            if line_name in lines_hfs:
                print(f'{"[INFO]":<10}', f'HFS mask for {line_name}.')
                mask_hfs = _build_hfs_mask(mask.value, line_name, hfs_data, this_data)
                if mask_hfs is not None:
                    this_data[f"SPEC_MASK_{line_name.upper()}"] = Column(
                        mask_hfs, unit=au.dimensionless_unscaled,
                        description=f"HFS mask for {line_name.upper()}")

    # Store combined mask
    this_data["SPEC_MASK"] = Column(
        mask, unit=au.dimensionless_unscaled,
        description="Velocity-integration mask (used for integrated products)")

    print(f'{"[INFO]":<10}', 'Mask done. Computing moments.')

    # ------------------------------------------------------------------
    # Loop over lines: moments + shuffle
    # ------------------------------------------------------------------
    cdelt      = shuff_axis[1] * au.m / au.s
    naxis_shuff = int(shuff_axis[0])
    new_vaxis  = (cdelt * (np.arange(naxis_shuff) - naxis_shuff / 2)).to(au.km / au.s)

    for jj in range(n_lines):
        line_name = line_names[jj]
        tag_spec  = "SPEC_" + line_name.upper()

        if tag_spec not in this_data.keys():
            print(f'{"[ERROR]":<10}', f'{tag_spec} not found. Skipping.')
            continue
        this_spec = this_data[tag_spec]
        if np.nansum(this_spec, axis=None) == 0:
            print(f'{"[ERROR]":<10}', f'{line_name} appears empty. Skipping.')
            continue

        dim_sz   = np.shape(this_spec)
        n_pts_l  = dim_sz[0]
        n_chan_l = dim_sz[1]

        this_v0     = this_data.meta["SPEC_VCHAN0"]
        this_deltav = this_data.meta["SPEC_DELTAV"]
        this_crpix  = this_data.meta["SPEC_CRPIX"]
        this_vaxis  = (this_v0 + (np.arange(n_chan_l) - (this_crpix - 1)) * this_deltav).to(au.km / au.s)

        this_data["SPEC_VAXIS"] = Column(
            np.array([this_vaxis] * n_pts_l), unit=au.km / au.s,
            description="Velocity axis")

        # Choose mask for this line
        if use_hfs_lines and line_name in lines_hfs:
            hfs_tag = f"SPEC_MASK_{line_name.upper()}"
            active_mask = this_data[hfs_tag] * au.Unit(1) if hfs_tag in this_data.keys() else mask
        else:
            active_mask = mask

        # Moments
        mom_maps = get_mom_maps(this_spec, active_mask, this_vaxis, mom_calc)
        line_desc = str(cubes["line_desc"].iloc[jj])

        # Only store moments if no 2D map was supplied for this cube
        band_ext_val = str(cubes["map_ext"].iloc[jj]).strip()
        if band_ext_val in ("", "nan"):
            this_data["MOM0_"  + line_name.upper()] = Column(mom_maps["mom0"],     description=f"{line_desc} integrated intensity (mom0)")
            this_data["EMOM0_" + line_name.upper()] = Column(mom_maps["mom0_err"], description=f"Error: {line_desc} mom0")
            this_data["TPEAK_" + line_name.upper()] = Column(mom_maps["tpeak"],    description=f"{line_desc} peak brightness")
            this_data["RMS_"   + line_name.upper()] = Column(mom_maps["rms"],      description=f"{line_desc} rms noise")
            this_data["MOM1_"  + line_name.upper()] = Column(mom_maps["mom1"],     description=f"{line_desc} mean velocity (mom1)")
            this_data["EMOM1_" + line_name.upper()] = Column(mom_maps["mom1_err"], description=f"Error: {line_desc} mom1")
            this_data["MOM2_"  + line_name.upper()] = Column(mom_maps["mom2"],     description=f"{line_desc} velocity dispersion (mom2; {mom_calc[2]})")
            this_data["EMOM2_" + line_name.upper()] = Column(mom_maps["mom2_err"], description=f"Error: {line_desc} mom2")
            this_data["EW_"    + line_name.upper()] = Column(mom_maps["ew"],       description=f"{line_desc} equivalent width")
            this_data["EEW_"   + line_name.upper()] = Column(mom_maps["ew_err"],   description=f"Error: {line_desc} EW")
            print(f'{"[INFO]":<10}', f'Moments computed for {line_name}.')
        else:
            print(f'{"[INFO]":<10}', f'2D map provided for {line_name}; skipping moments.')

        # Shuffle
        shuffled = shuffle(
            spec      = this_spec,
            vaxis     = this_vaxis,
            zero      = ref_line_vmean,
            new_vaxis = new_vaxis,
            interp    = 0,
        )
        this_data["SPEC_SHUFF" + line_name.upper()] = Column(
            shuffled, unit=this_spec.unit,
            description=f"Shuffled {line_desc} brightness temperature")
        this_data.meta["SPEC_VCHAN0_SHUFF"] = new_vaxis[0]
        this_data.meta["SPEC_DELTAV_SHUFF"] = new_vaxis[1] - new_vaxis[0]
        this_data["SPEC_VAXISSHUFF"] = Column(
            np.array([new_vaxis] * n_pts_l), unit=au.km / au.s,
            description="Shuffled velocity axis")

    this_data.meta["SPEC_CRPIX_SHUFF"] = 1
    this_data.write(fname, format="ascii.ecsv", overwrite=True)
    print(f'{"[INFO]":<10}', f'Spectra processed for {source}.')
