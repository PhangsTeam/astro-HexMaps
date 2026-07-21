.. _AdvancedConfig:

Advanced Configuration
======================

This page documents options that are less commonly changed but useful for
specific use cases.


Resolution Modes
----------------

.. code-block:: ini

   [resolution]
   resolution = angular    # arcseconds (default)
   # resolution = physical # target_res in parsecs, converted via distance
   # resolution = native   # use the overlay beam as-is

For ``resolution = physical``, set ``target_res`` in **parsecs**. HexMaps
converts to arcseconds per target using its distance from
``target_definitions.txt``.


Grid Parameters
---------------

.. code-block:: ini

   [resolution]
   pixels_per_beam = 2       # spacing = target_res / pixels_per_beam
   max_rad         = auto    # "auto" derives radius from overlay footprint
   NAXIS_shuff     = 200     # channels in the shuffled spectrum
   CDELT_SHUFF     = 4000.0  # channel width of shuffled spectrum [m/s]


Reference Line and Mask Combinations
--------------------------------------

The ``ref_line`` key is a comma-separated list of tokens. All tokens are
optional except that at least one line-selection or external-mask token must
be present.

**Line-selection tokens** (choose exactly one):

.. code-block:: ini

   ref_line = first         # first cube (default)
   ref_line = 12co21        # specific named line
   ref_line = all           # OR-combine all lines
   ref_line = 2             # first 2 lines
   ref_line = individual    # one mask per line, applied independently

**External-mask tokens** (append to any line selection):

.. code-block:: ini

   input    # include the input_mask defined in the mask table
   window   # include the window_mask defined in the mask table

**Combinator token** (default ``OR``):

.. code-block:: ini

   OR       # a sightline passes if it passes ANY mask (default)
   AND      # a sightline passes only if it passes ALL masks

Examples:

.. code-block:: ini

   ref_line = 12co21, input, AND    # 12co21 S/N mask AND input mask
   ref_line = first, window, AND    # first-line mask AND velocity window
   ref_line = all, input, window    # OR of all-line + input + window masks
   ref_line = individual, input     # per-line S/N mask OR input mask

The external masks (``input_mask``, ``window_mask``) must be defined in the
``# ---- mask ----`` table of ``config.txt``.

**Two-level S/N mask:**

.. code-block:: ini

   SN_processing   = 2, 4   # [low_SN, high_SN] thresholds
   strict_mask     = false  # false | strict | broad
   conseq_channels = 3      # min consecutive channels for a valid detection


Strict and Broad Mask Modes
-----------------------------

.. warning::

   The ``strict`` and ``broad`` options for ``strict_mask`` are experimental
   and have not been thoroughly tested across a wide range of datasets.
   Use the default value (``false``) unless you have a specific reason to
   apply a coherence filter.  If you do use ``strict`` or ``broad``, treat
   the results with care and validate them against your expectations for the
   source structure and noise level.

The ``strict_mask`` key controls an optional post-processing coherence filter:

* ``false`` — no additional spatial filtering (default)
* ``strict`` — remove connected components smaller than approximately one
  beam area per channel. Analogous to the ACES beam-area pruning step;
  suppresses isolated noise spikes while preserving spatially extended emission.
* ``broad`` — re-derive the mask from a spatially smoothed cube using a
  two-level S/N dilation strategy (core mask at high S/N, grown into a
  wing mask at low S/N), following the PHANGS-ALMA broad-mask approach.
  More inclusive than the raw S/N mask; better suited for faint extended
  emission.

.. code-block:: ini

   strict_mask = false   # no filter
   strict_mask = strict  # beam-area pruning
   strict_mask = broad   # smoothed-cube two-level dilation


Velocity Windows
----------------

Explicit velocity windows for signal integration and noise estimation are
defined in the ``# ---- mask ----`` table:

.. code-block:: ini

   # ---- mask ----
   window_mask = window_mask, Fixed velocity window, -200, 200, km/s
   noise_mask  = noise_mask,  Noise blue,            -300, -150, km/s
   noise_mask  = noise_mask,  Noise red,              150,  300, km/s

To use the velocity window for signal integration, add ``window`` to
``ref_line``:

.. code-block:: ini

   ref_line = first, window, AND    # first-line S/N mask AND velocity window

To use fixed noise windows for RMS estimation instead of channels outside the
signal mask:

.. code-block:: ini

   [masking]
   use_fixed_noise_mask = true


Noise Mask Overlap Check
------------------------

When ``use_fixed_noise_mask = true``, the pipeline automatically removes
any channels from the noise window that overlap with the signal integration
mask before computing per-sightline RMS. This prevents signal contamination
of the noise estimate even when the noise windows are defined broadly.


Spectral Smoothing
------------------

.. code-block:: ini

   [spectral]
   spec_smooth        = default   # no smoothing
   # spec_smooth      = overlay   # smooth to overlay spectral resolution
   # spec_smooth      = 5.0       # convolve to 5.0 km/s

   spec_smooth_method = binned    # recommended
   # spec_smooth_method = gauss       # ±10-15% RMS bias; avoid for science
   # spec_smooth_method = combined    # bin first, then Gaussian residual


Hyperfine Structure Correction
--------------------------------

For lines with hyperfine structure (HCN, N₂H⁺, CN, CCH), HexMaps can
extend the signal mask to satellite components:

.. code-block:: ini

   [paths]
   hfs_file = keys/hfs_lines.txt

   [masking]
   use_hfs_lines = true

Add entries for each line with hyperfine structure to ``hfs_lines.txt``.

When ``use_hfs_lines = true``:

* In **combined mode**: the master mask is extended to each HFS-capable
  line's satellite frequencies after it is built, producing a per-line mask
  stored as ``SPEC_MASK_<LINE>``. Moments for each line are computed using
  its own extended mask.
* In **individual mode** (``ref_line = individual``): the per-line S/N mask
  already finds the satellite emission by construction, so no additional HFS
  extension is applied to it. External masks (``input``, ``window``) are
  still extended to the satellite frequencies before combining.


FOV Edge Erosion
----------------

Pixels near the map boundary are computed from a partial kernel and are
biased. HexMaps trims these by eroding the footprint:

.. code-block:: ini

   [masking]
   fov_erosion_beams = 0.5   # trim 0.5 × beam FWHM (default)
   # fov_erosion_beams = 0   # disable — keep full overlay footprint
   # fov_erosion_beams = 1.0 # conservative — trim one full beam

The same erosion is applied to the hex-grid, moment maps, and FITS outputs.


Database Fill Mode
------------------

To add new maps or cubes to an existing ``.ecsv`` without re-running the
full pipeline:

.. code-block:: ini

   [structure]
   structure_creation = fill
   fname_fill = ngc5194_hexmaps_27p0as_2025_01_01.ecsv

HexMaps opens the existing file and adds only the maps/cubes that are not
yet present.
