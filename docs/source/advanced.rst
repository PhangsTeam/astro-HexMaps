Advanced Configuration
======================

This page documents options that are less commonly changed but useful for
specific use cases.


Resolution Modes
----------------

The ``resolution`` key in ``[resolution]`` controls how the target resolution
is interpreted:

.. code-block:: ini

   [resolution]
   resolution = angular    # use target_res in arcseconds (default)
   # resolution = physical # convert target_res (in pc) to arcsec via distance
   # resolution = native   # use the beam size of the overlay cube as-is

When ``resolution = physical``, set ``target_res`` in **parsecs**. HexMaps
reads the distance from ``keys/target_definitions.txt`` and converts
automatically.

When ``resolution = native``, the beam of the overlay cube is used unchanged
and no convolution is applied to the overlay itself.


Grid Parameters
---------------

.. code-block:: ini

   [resolution]
   pixels_per_beam = 2      # hex spacing = target_res / pixels_per_beam
   max_rad         = auto   # maximum map radius in degrees; "auto" derives it
                            # from the overlay cube footprint
   NAXIS_shuff     = 200    # number of channels in the shuffled spectrum
   CDELT_SHUFF     = 4000.0 # channel width of shuffled spectrum [m/s]

``pixels_per_beam = 2`` gives half-beam spacing (Nyquist sampling of the
beam), which is the standard for moment-map analysis. Increase for sparser
sampling (faster), decrease for denser (slower, more points).


FOV Edge Erosion
----------------

When convolving to a target beam, pixels near the map boundary are computed
from a partial kernel footprint and are therefore biased. HexMaps removes
these pixels by eroding the observed footprint by a configurable margin:

.. code-block:: ini

   [masking]
   fov_erosion_beams = 0.5   # trim by 0.5 × beam FWHM (default)
   # fov_erosion_beams = 0   # disable — keep full overlay footprint
   # fov_erosion_beams = 1.0 # conservative — trim one full beam

The same erosion is applied to the hexagonal sampling grid, the moment maps,
and any saved FITS cubes, so all outputs share an identical footprint.


Reference Line and Masking
---------------------------

The ``ref_line`` key controls which cube(s) define the signal mask:

.. code-block:: ini

   [masking]
   ref_line = first       # use first cube (default — recommended: use brightest line)
   # ref_line = 12co21    # use a specific named cube
   # ref_line = all       # OR-combine all cubes into one mask
   # ref_line = 2         # use first 2 cubes
   # ref_line = ref+HI    # use first cube combined with HI

The two-level S/N mask is controlled by:

.. code-block:: ini

   SN_processing = 2, 4   # [low_SN, high_SN]
   strict_mask   = false  # if true, apply spatial connectivity filter

Pixels above ``high_SN`` seed the mask; connected pixels above ``low_SN``
are included. Enabling ``strict_mask`` removes isolated detections that are
not spatially connected to the main signal region.


Velocity Windows
----------------

Instead of (or in addition to) the S/N mask, you can define explicit velocity
windows for signal integration and noise estimation in the ``# ---- mask ----``
table:

.. code-block:: ini

   # ---- mask ----
   # Fixed velocity window for signal integration:
   vel_mask,  Signal window,  -200,  200,  km/s

   # Noise estimation windows (multiple ranges are OR-combined):
   noise_mask,  Noise blue,  -300,  -150,  km/s
   noise_mask,  Noise red,    150,   300,  km/s

Enable these with:

.. code-block:: ini

   [masking]
   use_fixed_vel_mask   = true
   use_fixed_noise_mask = true

When ``use_fixed_noise_mask = true``, the noise (RMS) is estimated from
channels inside the noise windows rather than from the inverse of the signal
mask. This is useful when the spectral baseline is contaminated by other lines.


Spectral Smoothing
------------------

.. code-block:: ini

   [spectral]
   spec_smooth        = default  # no smoothing
   # spec_smooth      = overlay  # smooth to overlay spectral resolution
   # spec_smooth      = 5.0      # convolve to 5.0 km/s resolution

   spec_smooth_method = binned   # recommended
   # spec_smooth_method = gauss      # ±10-15% RMS bias — use with caution
   # spec_smooth_method = combined   # bin first, then Gaussian residual


Hyperfine Structure Correction
--------------------------------

For lines with hyperfine structure (e.g. HCN, N₂H⁺, CN, CCH), HexMaps can
extend the signal mask to cover satellite components:

.. code-block:: ini

   [paths]
   hfs_file = keys/hfs_lines.txt

   [masking]
   use_hfs_lines = true

The ``keys/hfs_lines.txt`` file lists, for each line, the reference frequency
and the frequencies of all hyperfine satellites. Add entries for any line in
your cube list that has known hyperfine structure.


Database Fill Mode
------------------

To add new maps or cubes to an existing ``.ecsv`` without re-running the full
pipeline:

.. code-block:: ini

   [structure]
   structure_creation = fill

HexMaps will open the existing file and add only the maps/cubes that are not
yet present, then save the updated file. Optionally pin the filename:

.. code-block:: ini

   fname_fill = ngc5194_hexmaps_27p0as_2025_01_01.ecsv
