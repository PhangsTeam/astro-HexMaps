.. _Analysis:

Working with HexMaps Output
============================

The Output File
---------------

Running HexMaps produces an Astropy Enhanced CSV (``.ecsv``) file in the
``output/`` directory. The filename follows the pattern::

   <source>_hexmaps_<res_suffix>_<date>.ecsv

For example: ``ngc5194_hexmaps_27p0as_2025_01_01.ecsv``

The ``.ecsv`` format is human-readable plain text; units and metadata are
stored in a comment header. It can be opened with any tool that reads CSV,
but the recommended approach is via Astropy:

.. code-block:: python

   from astropy.table import Table

   table = Table.read("output/ngc5194_hexmaps_27p0as_2025_01_01.ecsv")
   print(table.colnames)

Or using the HexMaps convenience loader:

.. code-block:: python

   from hexmaps.utils_table import load_hexmaps

   table = load_hexmaps("output/ngc5194_hexmaps_27p0as_2025_01_01.ecsv")


Column Naming Conventions
--------------------------

Columns follow a consistent naming pattern:

+------------------------+-------------------------------------------------------+
| Prefix                 | Content                                               |
+========================+=======================================================+
| ``ra_deg``, ``dec_deg``| Hexagonal sightline coordinates                      |
+------------------------+-------------------------------------------------------+
| ``rgal_kpc``           | Deprojected galactocentric radius                     |
+------------------------+-------------------------------------------------------+
| ``SPEC_<LINE>``        | Full spectrum at each sightline (n_pts × n_chan)      |
+------------------------+-------------------------------------------------------+
| ``MOM0_<LINE>``        | Integrated intensity (moment 0)                       |
+------------------------+-------------------------------------------------------+
| ``MOM1_<LINE>``        | Intensity-weighted mean velocity (moment 1)           |
+------------------------+-------------------------------------------------------+
| ``MOM2_<LINE>``        | Intensity-weighted line width (moment 2)              |
+------------------------+-------------------------------------------------------+
| ``RMS_<LINE>``         | Per-sightline RMS noise                               |
+------------------------+-------------------------------------------------------+
| ``TPEAK_<LINE>``       | Peak brightness temperature                           |
+------------------------+-------------------------------------------------------+
| ``EW_<LINE>``          | Equivalent width                                      |
+------------------------+-------------------------------------------------------+
| ``MAP_<NAME>``         | 2D band map values                                    |
+------------------------+-------------------------------------------------------+
| ``SPEC_VAXIS``         | Velocity axis array (km/s)                            |
+------------------------+-------------------------------------------------------+

For example: ``MOM0_12CO21``, ``RMS_12CO10``, ``MAP_SPIRE250``.


The HexMapsAnalysis Class
--------------------------

``analysis/hexmaps_analysis.py`` wraps the table in a convenience class with
quicklook plotting and data extraction methods:

.. code-block:: python

   import sys
   sys.path.append("analysis/")
   from hexmaps_analysis import HexMapsAnalysis

   db = HexMapsAnalysis("output/ngc5194_hexmaps_27p0as_2025_01_01.ecsv")
   print(db)
   # HexMapsAnalysis(source='ngc5194', n_pts=939, lines=['12CO21', '12CO10'])


Quick Examples
--------------

**List spectral lines in the database:**

.. code-block:: python

   print(db.lines)
   # ['12CO21', '12CO10']

**Plot a 2D moment map:**

.. code-block:: python

   db.quickplot_map("12CO21")

.. image:: quicklook2.png
   :width: 400

**Plot a spectrum at the brightest sightline:**

.. code-block:: python

   db.quickplot_spectrum("12CO21")

.. image:: spec.png
   :width: 600

**Custom 2D scatter map:**

.. code-block:: python

   import matplotlib.pyplot as plt

   ra, dec = db.get_coordinates("13:29:52.7 47:11:43")
   mom0 = db.struct["MOM0_12CO21"]

   fig, ax = plt.subplots(figsize=(5, 5))
   sc = ax.scatter(ra, dec, c=mom0, s=90, marker="h", cmap="inferno")
   ax.invert_xaxis()
   ax.set_xlabel(r"$\Delta$R.A. [arcsec]")
   ax.set_ylabel(r"$\Delta$Decl. [arcsec]")
   plt.colorbar(sc, label="MOM0 [K km/s]")
   plt.show()

.. image:: map_2D.png
   :width: 400

**Compute a line ratio:**

.. code-block:: python

   ratio = db.get_ratio("12CO21", "12CO10", sn=5.0)
   print(ratio["ratio"])   # array of CO(2-1)/CO(1-0) ratios

**Plot a radial profile:**

.. code-block:: python

   db.quickplot_radial_profile("12CO21")

**Recover provenance information:**

.. code-block:: python

   # Print the config.txt that was used to produce this database
   print(db.get_config())

   # Save it to a file
   db.get_config(save_to="recovered_config.txt")

   # List all embedded raw FITS headers
   print(db.list_input_headers())
   # ['12CO10', '12CO21', 'OVERLAY', 'SPIRE250']

   # Recover a specific header (returns an astropy.io.fits.Header object)
   hdr = db.get_input_header("12CO21")
   print(f"Native beam: {hdr['BMAJ'] * 3600:.1f} arcsec")


Accessing the Raw Table
------------------------

The underlying Astropy table is always accessible via ``db.struct``:

.. code-block:: python

   # All column names
   print(db.struct.colnames)

   # Galactocentric radii in kpc
   print(db.struct["rgal_kpc"])

   # Spectrum of the brightest CO(2-1) sightline
   import numpy as np
   idx = np.argmax(db.struct["MOM0_12CO21"])
   spec = db.struct["SPEC_12CO21"][idx]
   vaxis = db.struct["SPEC_VAXIS"][idx]
