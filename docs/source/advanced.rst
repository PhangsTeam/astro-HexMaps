.. _AdvancedConfig:

Advanced Configuration
======================

This page documents every key in the ``[resolution]``, ``[masking]``, and
``[spectral]`` sections that is not covered on the main :doc:`config` page,
plus a dedicated section on each advanced masking strategy.  The keys are
discussed one by one.

----


[resolution]
------------

Grid Sampling
~~~~~~~~~~~~~

.. code-block:: ini

   [resolution]
   pixels_per_beam = 2
   max_rad = auto

``pixels_per_beam``
   Number of hex-grid sampling points per beam FWHM.  The grid spacing
   (centre-to-centre distance between adjacent hex points) is:

   .. math::

      \Delta\theta = \frac{\text{target_res}}{\text{pixels_per_beam}}

   * ``2`` *(default)* — half-beam spacing.  Adjacent sightlines share
     roughly 50 % of their beam area; a good balance between spatial
     sampling and independence.
   * ``1`` — one sightline per beam diameter; the coarsest sensible grid.
   * Values above ``2`` produce denser grids where sightlines overlap
     substantially; useful for spectral stacking but increases file size
     and processing time.

   *Default:* ``2``

``max_rad``
   Maximum radius of the hexagonal grid in degrees, measured from the
   target centre coordinate (``x_ctr``, ``y_ctr`` in
   ``target_definitions.txt``).

   * ``auto`` *(default)* — the radius is derived from the overlap between
     the overlay cube footprint and the available data coverage.
   * A positive float — restricts the grid to this radius, useful when the
     overlay is large but only the central region is of interest.

   *Default:* ``auto``

----

Shuffled Spectrum
~~~~~~~~~~~~~~~~~

.. code-block:: ini

   [resolution]
   NAXIS_shuff = 200
   CDELT_SHUFF = 4000.0

``NAXIS_shuff``
   Number of spectral channels in the shuffled spectrum output.  Each
   sightline's spectrum is shifted so that the systemic velocity sits at
   channel zero; the shuffled spectrum then runs from
   ``-NAXIS_shuff/2 × CDELT_SHUFF`` to ``+NAXIS_shuff/2 × CDELT_SHUFF``
   metres per second.

   Choose a value wide enough to cover the full line width of the fastest-
   rotating targets in your sample.

   *Default:* ``200``

``CDELT_SHUFF``
   Channel width of the shuffled spectrum in **m/s**.  This sets the
   velocity resolution of all shuffled spectra regardless of the native
   channel width of the input cubes.  A coarser value reduces the output
   file size and is sufficient for stacking studies where the native
   resolution is not needed.

   *Default:* ``4000.0``

----


[masking]
---------

Mask table
~~~~~~~~~~

External masks are defined as rows after the ``# ---- mask ----`` comment
marker. Three types of entry are supported, each distinguished by its key:

.. code-block:: ini

   # ---- mask ----

   # File mask — a pre-computed binary FITS mask sampled onto the hex grid:
   # input_mask = name, description, file_ext, directory
   # input_mask = co_mask, CO signal mask, _co_mask.fits, data/

   # Velocity-window mask — channels within a fixed velocity range:
   # window_mask = name, description, v_start, v_end, unit
   # window_mask = win, Fixed velocity window, 400, 600, km/s

   # Noise velocity windows — line-free channels for RMS estimation:
   # noise_mask = name, description, v_start, v_end, unit
   # noise_mask = noise_b, Noise blue, -300, -150, km/s
   # noise_mask = noise_r, Noise red,   150,  300, km/s

``input_mask``
   An external binary FITS file that is sampled onto the hexagonal grid at
   the regrid stage and stored as a ``SPEC_<name>`` column. To use it as
   part of the signal mask, add the ``input`` token to ``ref_line``.

   *Columns:* ``name``, ``description``, ``file_ext``, ``directory``

``window_mask``
   A fixed velocity window defined by start and end velocities. All channels
   within the window are set to 1; all others to 0. To use it as part of the
   signal mask, add the ``window`` token to ``ref_line``.

   *Columns:* ``name``, ``description``, ``v_start``, ``v_end``, ``unit``
   (``unit`` must be an ``astropy.units``-readable velocity unit, e.g.
   ``km/s``).

``noise_mask``
   One or more line-free velocity windows used for per-sightline RMS
   estimation. Multiple rows are OR-combined into a single noise channel mask.
   Enabled by ``use_fixed_noise_mask = true`` in ``[masking]``.

   *Columns:* ``name``, ``description``, ``v_start``, ``v_end``, ``unit``

   The pipeline automatically excludes any noise channels that overlap with
   the signal integration mask, preventing signal contamination of the RMS
   estimate.

----


Reference Line and Mask Combinations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: ini

   [masking]
   ref_line = first

``ref_line``
   The central masking control.  A comma-separated list of tokens that
   specifies which masks are built, what external data are included, and
   how everything is combined.  At least one token must be present.

   **Line-selection tokens** — choose *exactly one*:

   ``first``
      Build the S/N mask from the first cube in the ``# ---- cubes ----``
      table.  This is the default when no token is given.  Put your
      brightest, highest-SNR line first in the cube list.

   ``<LINE_NAME>``
      A specific named line (e.g. ``12co21``), matched case-insensitively
      against the cube names.  For two or more named lines, give them as a
      comma-separated sub-list (e.g. ``12co21, 12co10``); their S/N masks
      are OR-combined into the master mask.

   ``all``
      OR-combine the S/N masks from every cube in the list into a single
      master mask.  Useful when you want any detected line to define the
      integration window.

   ``<n>`` (positive integer)
      OR-combine the S/N masks from the first *n* cubes in the list
      (``1`` is equivalent to ``first``).

   ``individual``
      Build one independent S/N mask per cube and apply each mask only to
      that line's moment computation.  No master mask is created.  External
      masks (``input``, ``window``) are combined with each per-line mask
      individually; when HFS correction is active, the external mask is also
      extended to the satellite frequencies for that line.

   **External-mask tokens** — zero or more, appended after the line-selection
   token:

   ``input``
      Include the pre-sampled binary mask defined by an ``input_mask`` row
      in the mask table.  The mask is OR- or AND-combined with the S/N mask
      according to the combinator.

   ``window``
      Include the fixed velocity-window mask defined by a ``window_mask``
      row in the mask table.  Channels inside the window are 1, all others 0.

   **Combinator token** — optional, default ``OR``:

   ``OR``
      A sightline/channel is included in the mask if it passes *any* of
      the selected masks (logical union).  This is the most inclusive option.

   ``AND``
      A sightline/channel is included only if it passes *all* selected
      masks (logical intersection).  Use this to restrict integration to
      the overlap of a bright-line S/N mask and an external spatial mask.

   **Examples:**

   .. code-block:: ini

      ref_line = first                   # default: S/N mask from first cube
      ref_line = 12co21                  # S/N mask from 12co21
      ref_line = all                     # OR of all cube S/N masks
      ref_line = 2                       # OR of first two cube S/N masks
      ref_line = individual              # one independent mask per cube
      ref_line = first, input            # OR of first-cube S/N + input mask
      ref_line = 12co21, input, AND      # 12co21 S/N AND input mask
      ref_line = first, window, AND      # first-cube S/N AND velocity window
      ref_line = all, input, window      # OR of all-cube + input + window
      ref_line = individual, input, AND  # per-line S/N AND input mask

   *Default:* ``first``

----

Mask Tuning
~~~~~~~~~~~

.. code-block:: ini

   [masking]
   SN_processing = 2, 4
   use_fixed_noise_mask = false
   fov_erosion_beams = 0.5
   conseq_channels = 3
   strict_mask = false

``SN_processing``
   Two S/N thresholds ``low_SN, high_SN`` for the two-level mask
   construction.

   1. A **core mask** is built from all channels where the local S/N exceeds
      ``high_SN`` in at least ``conseq_channels`` consecutive channels.
   2. The core mask is **grown** into adjacent channels that exceed ``low_SN``
      in at least ``conseq_channels`` consecutive channels, capturing the
      line wings connected to bright cores.

   A higher ``high_SN`` reduces false positives; a lower ``low_SN`` captures
   more of the faint line wings. Typical values for extragalactic surveys are
   ``2, 4`` or ``3, 5``.

   *Default:* ``2, 4``

``use_fixed_noise_mask``
   When ``true``, use the velocity windows defined by ``noise_mask`` rows
   in the mask table for per-sightline RMS estimation, instead of the
   channels outside the integration mask.  Useful when the baseline
   contains emission from other lines or instrumental artefacts that would
   bias the noise estimate upward.

   The pipeline automatically removes any noise channels that overlap with
   the signal integration mask before computing RMS, preventing signal
   contamination even when noise windows are defined broadly.

   Requires at least one ``noise_mask`` row in the mask table.

   *Default:* ``false``

``fov_erosion_beams``
   Trim the effective field of view by this multiple of the beam FWHM.
   Pixels near the map edge where the convolution kernel extends beyond the
   observed area are biased; erosion removes them.

   * ``0`` — disable; keep the full overlay footprint
   * ``0.5`` *(default)* — recommended minimum; trims half a beam
   * ``1.0`` — conservative; trims one full beam

   The same erosion mask is applied consistently to the hex-grid footprint,
   all moment map FITS files, and all PPV mask cubes, ensuring a uniform
   effective FOV across all output products.

   *Default:* ``0.5``

``conseq_channels``
   Minimum number of consecutive channels that must exceed the S/N threshold
   (at both the core and wing levels) for a detection to be considered
   genuine.  Isolated single-channel spikes are rejected even when they
   exceed ``SN_processing[1]``.

   Increasing this value makes the mask more conservative; decreasing it
   to ``1`` effectively removes the spectral coherence requirement.

   *Default:* ``3``

``strict_mask``
   Optional post-processing coherence filter applied to the signal mask
   after it is fully assembled.  See the dedicated section below for
   details and caveats.

   * ``false`` *(default)* — no additional filtering
   * ``strict`` — remove spatially isolated detections smaller than ~1 beam
     area per channel
   * ``broad`` — re-derive the mask from a spatially smoothed cube with
     two-level S/N dilation

   *Default:* ``false``

.. warning::

   The ``strict`` and ``broad`` options are experimental and have not been
   thoroughly tested across a wide range of datasets.  Use the default
   value (``false``) unless you have a specific reason to apply a coherence
   filter.  If you do use ``strict`` or ``broad``, treat the results with
   care and validate them against your expectations for the source structure
   and noise level.

The motivation for both modes is the same: real emission is spatially
coherent across the beam and persists over multiple spectral channels, while
noise peaks that happen to exceed the S/N threshold are typically isolated
in space or in velocity.

``strict``
   After assembling the mask, connected components are identified in each
   channel independently using 2D spatial connectivity (hex-grid path) or
   4-connected pixels (PPV path).  Any component smaller than approximately
   one beam area is removed.  This is analogous to the ACES beam-area
   pruning step and suppresses isolated noise spikes while preserving
   spatially extended emission.

   Best suited for high-angular-resolution interferometric data where the
   beam covers many independent pixels and isolated spikes are common.

``broad``
   The mask is re-derived from scratch using a spatially smoothed version
   of the cube, following the PHANGS-ALMA broad-mask strategy:

   1. Each channel is convolved with a Gaussian of σ ≈ 1 beam to enhance
      the S/N of faint, spatially coherent emission.
   2. A core mask is built from the smoothed cube at the ``high_SN``
      threshold with the consecutive-channel requirement.
   3. The core mask is dilated into adjacent voxels that exceed ``low_SN``
      in the smoothed cube (5 dilation passes + 2 channel-grow passes).
   4. The resulting mask is applied to the **original** (unsmoothed) spectra
      for moment computation.

   The broad mask captures faint line wings that the unsmoothed S/N mask
   would miss, at the cost of slightly reduced spatial resolution in the
   mask boundary.  It is better suited for faint, extended emission in
   single-dish data.

----

Hyperfine Structure Lines
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: ini

   [masking]
   use_hfs_lines = false

``use_hfs_lines``
   When ``true``, extend the signal mask to the hyperfine satellite
   frequencies of lines defined in ``hfs_file``.  Requires ``hfs_file`` to
   be set in ``[paths]``.  See the :doc:`hfs_lines` page for the file
   format and a description of how the extension works.
   
   *Default:* ``false``

The behaviour depends on the masking mode:

**Combined mode** (any ``ref_line`` other than ``individual``):
   After the master mask is assembled, HexMaps loops over every HFS-capable
   line listed in ``hfs_lines.txt``.  For each satellite component, the
   master mask is shifted by the velocity offset between the satellite and
   main frequency and OR-combined with the original mask.  The result is
   stored as ``SPEC_MASK_<LINE>`` in the database; moments for that line
   are computed using its own extended mask rather than the global master
   mask.  If external masks (``input``, ``window``) are active, they are
   extended to the satellite frequencies as well before being combined with
   the S/N mask.

**Individual mode** (``ref_line = individual``):
   Each line's S/N mask is built independently from that line's own cube,
   so the satellite emission is naturally detected if it exceeds the S/N
   threshold.  The S/N mask itself is therefore *not* extended further.
   External masks (``input``, ``window``) *are* still extended to the
   satellite frequencies and combined with the per-line S/N mask.

----

Moment Computation
~~~~~~~~~~~~~~~~~~

.. code-block:: ini

   [masking]
   mom_thresh = 5
   mom2_method = fwhm

``mom_thresh``
   S/N threshold for moment-1 (mean velocity), moment-2 (line width), and
   equivalent-width computation.  Sightlines whose peak S/N falls below
   this value are excluded from those quantities and receive ``NaN`` in
   the output.  Moment-0 (integrated intensity) is computed for all
   sightlines that have any masked channels.

   *Default:* ``5``

``mom2_method``
   Definition used for the line-width (moment-2) output column:

   * ``fwhm`` *(default)* — converts the intensity-weighted second moment
     to a full-width at half-maximum by multiplying by
     :math:`2\sqrt{2\ln 2} \approx 2.355`.  Comparable to the FWHM of a
     Gaussian fit.
   * ``sqrt`` — returns :math:`\sqrt{\mu_2}`, the intensity-weighted
     velocity dispersion (σ).  Related to FWHM by a factor of 2.355.
   * ``math`` — returns the raw mathematical second moment :math:`\mu_2`
     in (km/s)².

   *Default:* ``fwhm``

----


[spectral]
----------

.. code-block:: ini

   [spectral]
   spec_smooth = default
   spec_smooth_method = binned

``spec_smooth``
   Spectral smoothing applied to each cube before the spectra are sampled
   onto the hex grid.

   * ``default`` *(default)* — no smoothing; the native spectral resolution
     of each cube is preserved.
   * ``overlay`` — smooth to the spectral resolution of the overlay cube.
     Useful when different cubes have different native channel widths and you
     want a uniform spectral resolution across all lines.
   * A positive float (e.g. ``5.0``) — convolve to the specified velocity
     resolution in km/s.  All cubes are smoothed to the same resolution
     regardless of their native channel width.

   *Default:* ``default``

``spec_smooth_method``
   Algorithm used when ``spec_smooth`` is not ``default``:

   ``binned``
      Bin adjacent channels by the nearest integer ratio.  Computationally
      fast; produces independent output channels.  This is the recommended
      method for science-quality results. *(default)*

   ``gauss``
      Convolve with a Gaussian kernel.  Note that this can underestimate
      the per-channel RMS by 10–15 % in low-S/N regions because adjacent
      output channels become correlated.  Avoid for precision noise
      characterisation.

   ``combined``
      Bin first (integer ratio), then apply a small Gaussian to handle the
      fractional remainder.  A compromise between ``binned`` and ``gauss``.

   *Default:* ``binned``

----


[output]
--------

.. code-block:: ini

   [output]
   save_cubes = false
   save_mom_maps = true
   save_maps = true
   save_mask = false

``save_cubes``
   When ``true``, write each convolved PPV cube to a FITS file in
   ``folder_savefits``. Only applies to the *fits* stage
   (``--stages all``).

   *Default:* ``false``

``save_mom_maps``
   When ``true``, write moment maps (mom0, mom1, mom2, rms, Tpeak, EW, and
   their error maps) to FITS files in ``folder_savefits``.

   *Default:* ``true``

``save_maps``
   When ``true``, write the convolved and reprojected 2D band maps to FITS
   files in ``folder_savefits``.

   *Default:* ``true``

``save_mask``
   When ``true``, write the velocity-integration mask(s) to FITS cubes in
   ``folder_savefits``. One file is written for the combined mask; additional
   files are written per line when HFS masks or per-line masks are active.

   *Default:* ``false``

----


[structure]
-----------

.. code-block:: ini

   [structure]
   structure_creation = default
   # fname_fill = ngc5194_hexmaps_27p0as_2025_01_01.ecsv

``structure_creation``
   Controls how the pipeline handles existing ``.ecsv`` files in ``out_dir``:

   * ``default`` *(default)* — create or overwrite the output file each run.
   * ``fill`` — open an existing file and add only maps/cubes not yet present.
     The file to open is identified by ``fname_fill`` (or by searching
     ``out_dir`` for the most recent matching file if ``fname_fill`` is not
     set).  Useful for incrementally building up a database as new
     observations are reduced.
   * ``archive`` — append a timestamp to the filename and write a new
     versioned copy each run.  The original file is never overwritten.

   *Default:* ``default``

``fname_fill``
   Explicit path or filename of the ``.ecsv`` to open when
   ``structure_creation = fill``.  If left blank, HexMaps searches ``out_dir``
   for the most recent ``.ecsv`` matching the target name and resolution
   suffix.

   *Only used when* ``structure_creation = fill``.
