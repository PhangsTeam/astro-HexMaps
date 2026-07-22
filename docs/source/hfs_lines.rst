Hyperfine Structure Lines File
==============================

The ``keys/`` directory contains the two reference files that describe your
targets and spectral lines.  This page covers ``hfs_lines.txt``; see
:doc:`target_definitions` for the companion target geometry file.

The ``keys/hfs_lines.txt`` file defines the hyperfine structure (HFS)
satellite frequencies for spectral lines that have such transitions.  When
``use_hfs_lines = true`` is set in ``[masking]``, HexMaps reads this file
and extends the signal mask for each listed line to cover its satellite
components, ensuring that the integration window captures all hyperfine
emission.

This file is **optional**.  It is only read by the pipeline when
``use_hfs_lines = true``.  If you are not working with HFS lines (e.g.
pure CO surveys), you can leave this file empty or omit it entirely.


File Format
-----------

The file is a comma-separated text file with no header row.  Lines beginning
with ``#`` are comments and are ignored.  Whitespace around commas is ignored,
so columns can be aligned freely for readability.

Each row defines one satellite transition of a line.  A line with *n*
satellite components requires *n* rows, all sharing the same ``line`` name
and ``ref_freq``.

.. code-block:: text

   # line,     ref_freq,       hfs_freq,       unit
   hcn10,      88.6316023,     88.6304156,     GHz
   hcn10,      88.6316023,     88.6318475,     GHz
   hcn10,      88.6316023,     88.6339357,     GHz
   n2hp10,     93.17370000,    93.17188000,    GHz
   n2hp10,     93.17370000,    93.17613000,    GHz


Column Descriptions
-------------------

``line``
   Name of the spectral line.  Must match exactly (case-insensitive) the
   line name used in the ``# ---- cubes ----`` table of ``config.txt`` and
   in the ``.ecsv`` database column names.  For example, if the cube is
   listed as ``hcn10`` in ``config.txt``, use ``hcn10`` here.

``ref_freq``
   Rest frequency of the **main (brightest) hyperfine component**, in the
   units given by ``unit``.  This is the frequency at which the cube data
   are calibrated and to which the spectral axis refers.  It corresponds to
   the ``RESTFRQ`` keyword in the input FITS header.

``hfs_freq``
   Rest frequency of **one satellite hyperfine component**, in the same units
   as ``ref_freq``.  The pipeline computes the velocity offset between
   ``hfs_freq`` and ``ref_freq`` and shifts the mask by the corresponding
   number of channels.

   Add one row per satellite component.  The pipeline accumulates all shifts
   for a given line and OR-combines the shifted masks into the final per-line
   integration mask.

``unit``
   Frequency unit, readable by ``astropy.units`` (e.g. ``GHz``, ``MHz``).
   Both ``ref_freq`` and ``hfs_freq`` must be given in the same unit.


How the Pipeline Uses This File
--------------------------------

1. After the main signal mask is built (from the S/N threshold and any
   external masks), the pipeline looks up each HFS-capable line in this
   file.
2. For each satellite component of that line, it computes the velocity
   offset :math:`\Delta v = c \cdot (f_{\rm ref} - f_{\rm hfs}) / f_{\rm ref}`
   and shifts the mask by the corresponding number of channels.
3. The shifted copies are OR-combined with the original mask, producing a
   per-line mask (``SPEC_MASK_<LINE>``) that covers both the main line and
   all satellite components.
4. Moments for that line are computed using the extended per-line mask
   rather than the global master mask.

In **individual mode** (``ref_line = individual``), the per-line S/N mask
already detects all emission including satellite lines if their S/N exceeds
the threshold.  HFS extension is therefore not applied to the S/N mask
itself in individual mode.  External masks (``input``, ``window``) are still
extended to the satellite frequencies when ``use_hfs_lines = true``.


Supported Lines (Built-in Examples)
------------------------------------

The following lines have known hyperfine structure and are commonly included
in ``hfs_lines.txt``.  Frequencies are taken from the Cologne Database for
Molecular Spectroscopy (CDMS) and JPL Molecular Spectroscopy Catalog.

+------------+----------------------+-------------------------------------------+
| Line       | Main freq. (GHz)     | Notes                                     |
+============+======================+===========================================+
| HCN(1–0)   | 88.6316023           | 3 hyperfine components (F = 0–1, 1–1,    |
|            |                      | 2–1); strongest is F = 2–1               |
+------------+----------------------+-------------------------------------------+
| N₂H⁺(1–0) | 93.1737000           | 7 hyperfine components grouped in three  |
|            |                      | clusters; include at least 2 clusters    |
+------------+----------------------+-------------------------------------------+
| CN(1–0)    | 113.4909             | Multiple hyperfine components spread     |
|            |                      | over ~40 km/s                            |
+------------+----------------------+-------------------------------------------+
| CCH(1–0)   | 87.3169              | Two groups separated by ~50 km/s         |
+------------+----------------------+-------------------------------------------+

.. note::

   HexMaps does not have built-in frequency tables.  All HFS frequencies must
   be provided explicitly in ``hfs_lines.txt``.  We recommend cross-checking
   frequencies against the CDMS (https://cdms.astro.uni-koeln.de) or JPL
   (https://spec.jpl.nasa.gov) catalogues for your target transitions.


Adding a New Line
-----------------

To add HFS support for a new line:

1. Look up the rest frequency of the main component (used as
   ``ref_freq``) and the frequencies of all relevant satellite components
   from CDMS or JPL.

2. Verify that the ``line`` name matches the cube name in ``config.txt``
   exactly (case-insensitive).

3. Add one row per satellite component:

   .. code-block:: text

      # line,     ref_freq,       hfs_freq,       unit
      myline,     100.0000000,    99.9985000,     GHz
      myline,     100.0000000,    100.0015000,    GHz

4. Enable HFS correction in ``config.txt``:

   .. code-block:: ini

      [masking]
      use_hfs_lines = true

5. Ensure ``hfs_file`` in ``[paths]`` points to the correct file (or leave
   it at the default ``keys/hfs_lines.txt``).


Example File
------------

The following is a minimal example ``hfs_lines.txt`` covering HCN(1–0) and
N₂H⁺(1–0):

.. code-block:: text

   # =============================================================================
   # HexMaps hfs_lines.txt  —  Hyperfine structure line definitions
   # =============================================================================
   # Comma-separated; spaces and tabs around each comma are ignored.
   # Each row defines ONE satellite component.
   # Lines with n satellite components need n rows (same line + ref_freq).
   #
   # Columns:
   #   line      — line name (must match config.txt cube name, case-insensitive)
   #   ref_freq  — rest frequency of the main component
   #   hfs_freq  — rest frequency of the satellite component
   #   unit      — frequency unit (astropy.units readable, e.g. GHz)
   # =============================================================================
   # line,     ref_freq,       hfs_freq,       unit
   hcn10,      88.6316023,     88.6304156,     GHz
   hcn10,      88.6316023,     88.6318475,     GHz
   hcn10,      88.6316023,     88.6339357,     GHz
   n2hp10,     93.17370000,    93.17188000,    GHz
   n2hp10,     93.17370000,    93.17613000,    GHz
