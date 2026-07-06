Advanced Setup
==============

The configure file comes with a set of advanced specification. These can be adjusted in :ref:`step2` of the file. Here we provide details on the individual advanced options.

Guide for advanced settings
---------------------------

Overall Parameters
^^^^^^^^^^^^^^^^^^


..  code-block::

	NAXIS_shuff = 200
	CDELT_SHUFF = 4000.  #m/s
	spacing_per_beam = 2 #default, use half beam spacing

	# give number (in units deg) or set to "auto"
	max_rad = "auto" #default extension of the map in deg (increase, if you map is larger)



Resolution
^^^^^^^^^^

..  code-block::

	#angular: use target_res in as
	#physical: convert target_res (in pc) to as
	#native: use the angular resolution of the overlay image
	resolution = 'angular'

Save Products
^^^^^^^^^^^^^

..  code-block::

	# Save the convolved cubes & bands
	save_fits = False

	# Save the moment maps
	save_mom_maps = True

	#folder to save fits files in
	folder_savefits="./saved_FITS_files/"


Line Masking
^^^^^^^^^^^^

.. code-block::

	#Define which line to use as reference line for the spectral processing
	#"first": use first line in cube_list as reference line
	#"<LINE_NAME>": Use line name as reference line
	#"all": Use all lines in cube for mask
	#n: (integer) use first n lines as reference. n=0 is same result as "first".
	ref_line = "first"

	#define upper and lower mask threshold (S/N)
	SN_processing = [2,4]
	strict_mask= False

Moment Calculations
^^^^^^^^^^^^^^^^^^^

.. code-block::

	#define SN threshold for Mom1, Mom2 and EW calculation (for individual lines)
	mom_thresh = 5

	#differentiate between "fwhm", "sqrt", or "math"
	# math: use mathematical definition
	# sqrt: take square-root of mom2
	# fwhm: convert sqrt(mom2) to fwhm
	mom2_method = "fwhm"

Spectral Smoothing
^^^^^^^^^^^^^^^^^^

.. code-block::

	"""
	"default": Do not perform any spectral smoothing
	"overlay": Perform spectral smoothing to spectral resolution of overlay cube
	n: float – convolve to spectral resolution n [km/s] !!!Not yet correctly implemented -> 	highly oversampled
	"""
	spec_smooth = "default"

	"""
	define the way the spectral smoothing should be performed:
	"binned": binn channels together (to nearest integer of ratio theta_target/theta_nat)
	"gauss": perform convolution with gaussian kernel (theta_target^2-theta_nat^2)**0.5
	!!!! Warning, gaussian smoothing seems to systematicaly underestimate the rms by 10-15%
	"combined": do the binned smoothing first (to nearest integer ratio) and then the rest via 	Gauss
	"""
	spec_smooth_method = "binned"
