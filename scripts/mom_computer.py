import numpy as np

def get_mom_maps(spec_cube, mask, vaxis, mom_calc=[3, 3, "fwhm"], noise_mask=None):
    """
    Function to compute moment maps.

    :param spec_cube: 2D array (n_pts, n_chan) of spectral data with astropy units.
    :param mask: 2D array (n_pts, n_chan) velocity-integration mask (1 = signal channels).
    :param vaxis: 1D velocity axis array with astropy units.
    :param mom_calc: list [SNthresh, conseq_channels, mom2_method].
    :param noise_mask: optional 2D boolean/int array (n_pts, n_chan) or 1D (n_chan,) boolean
                       array marking channels to use for RMS estimation.  When provided,
                       the RMS is computed from those channels only (instead of all
                       off-signal channels, i.e. mask==0).  A 1D array is broadcast
                       across all spatial points.  If None (default), the original
                       behaviour is preserved: RMS is estimated from all channels where
                       mask==0 and spectrum!=0.
    """

    # --- strip units ONCE ---
    spec_vals = spec_cube.value            # (n_pts, n_chan)
    v_vals    = vaxis.value                # (n_chan,)
    dv        = abs(v_vals[0] - v_vals[1])      # scalar
    spec_unit = spec_cube.unit
    v_unit    = vaxis.unit

    # --- pre-process noise_mask ---
    if noise_mask is not None:
        # strip astropy units if present
        if hasattr(noise_mask, 'value'):
            noise_mask_vals = noise_mask.value
        else:
            noise_mask_vals = np.asarray(noise_mask)
        # broadcast 1D mask to 2D (n_pts, n_chan)
        if noise_mask_vals.ndim == 1:
            noise_mask_vals = np.broadcast_to(noise_mask_vals, spec_vals.shape)
        noise_mask_vals = noise_mask_vals.astype(bool)
    else:
        noise_mask_vals = None

    # --- set up output maps WITH units ---
    mom_maps = {}

    dim_sz = np.shape(spec_vals)
    n_pts  = dim_sz[0]
    # n_chan = dim_sz[1]

    SNthresh         = mom_calc[0]
    conseq_channels  = int(np.nanmax((float(mom_calc[1]),3)))
    mom2_method      = mom_calc[2]
    fac_mom2         = np.sqrt(8*np.log(2)) if mom2_method == "fwhm" else 1.0

    # Output array templates
    mom_maps["rms"]       = np.full(n_pts, np.nan) * spec_unit
    mom_maps["tpeak"]     = np.full(n_pts, np.nan) * spec_unit

    mom_maps["mom0"]      = np.full(n_pts, np.nan) * spec_unit * v_unit
    mom_maps["mom0_err"]  = np.full(n_pts, np.nan) * spec_unit * v_unit

    mom_maps["mom1"]      = np.full(n_pts, np.nan) * v_unit
    mom_maps["mom1_err"]  = np.full(n_pts, np.nan) * v_unit

    mom2_unit = v_unit if mom2_method == "fwhm" else v_unit**2
    mom_maps["mom2"]      = np.full(n_pts, np.nan) * mom2_unit
    mom_maps["mom2_err"]  = np.full(n_pts, np.nan) * mom2_unit

    mom_maps["ew"]        = np.full(n_pts, np.nan) * v_unit
    mom_maps["ew_err"]    = np.full(n_pts, np.nan) * v_unit

    # -------------------------------
    #       MAIN LOOP (NO UNITS)
    # -------------------------------
    for m in range(n_pts):

        spectrum = spec_vals[m, :]      # 1D float array
        mask_m   = mask[m, :]

        # Skip empty spectra
        if np.nansum(spectrum != 0) < 1:
            continue

        # ---------------- RMS ----------------
        if noise_mask_vals is not None:
            # Use the explicitly defined noise velocity window(s)
            noise_sel = noise_mask_vals[m]
            rms_data = spectrum[np.logical_and(noise_sel, spectrum != 0)]
        else:
            # Default: use all off-signal (mask==0) channels
            rms_data = spectrum[np.logical_and(mask_m == 0, spectrum != 0)]
        rms = np.nanstd(rms_data) if len(rms_data) > 0 else np.nan
        mom_maps["rms"][m] = rms * spec_unit

        # ---------------- Tpeak ----------------
        tpeak = np.nanmax(spectrum * mask_m)
        mom_maps["tpeak"][m] = tpeak * spec_unit

        # ---------------- Mom0 ----------------
        mom0 = np.nansum(spectrum * mask_m) * dv
        mom_maps["mom0"][m] = mom0 * spec_unit * v_unit

        mom0_err = np.sqrt(np.nansum(mask_m)) * rms * dv
        mom_maps["mom0_err"][m] = mom0_err * spec_unit * v_unit

        # ---------------- Build high-signal mask ----------------
        masked = (spectrum * mask_m > SNthresh * rms).astype(int)
        masked = ((masked + np.roll(masked,1) + np.roll(masked,-1)) >= 3).astype(int)

        if np.nansum(masked) < conseq_channels - 2:
            continue

        for _ in range(5):
            masked = ((masked + np.roll(masked,1) + np.roll(masked,-1)) >= 1).astype(int)

        # ---------------- Mom1 ----------------
        num1 = np.nansum(spectrum * v_vals * masked)
        den1 = np.nansum(spectrum * masked)

        mom1 = num1 / den1
        mom_maps["mom1"][m] = mom1 * v_unit

        numer = rms**2 * np.nansum(masked * (v_vals - mom1)**2)
        mom1_err = np.sqrt(numer / den1**2)
        mom_maps["mom1_err"][m] = mom1_err * v_unit

        # ---------------- Mom2 ----------------
        mom2_math = np.nansum(spectrum * masked * (v_vals - mom1)**2) / den1

        numer = rms**2 * np.nansum((masked * (v_vals - mom1)**2 - mom2_math)**2)
        mom2_err = np.sqrt(numer / den1**2)

        if mom2_method == "fwhm":
            mom_maps["mom2"][m]     = fac_mom2 * np.sqrt(mom2_math) * v_unit
            mom_maps["mom2_err"][m] = fac_mom2 * mom2_err / (2 * np.sqrt(mom2_math)) * v_unit
        else:
            mom_maps["mom2"][m]     = mom2_math * v_unit**2
            mom_maps["mom2_err"][m] = mom2_err * v_unit**2

        # ---------------- EW ----------------
        ew = np.nansum(spectrum * masked) * dv / tpeak / np.sqrt(2*np.pi)
        mom_maps["ew"][m] = ew * v_unit

        term1 = rms**2 * np.nansum(masked) * dv**2 / (2*np.pi * tpeak**2)
        term2 = (ew**2 - ew * dv/np.sqrt(2*np.pi))
        ew_err = np.sqrt(term1 + term2)
        mom_maps["ew_err"][m] = ew_err * v_unit

    return mom_maps