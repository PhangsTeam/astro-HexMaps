PyStructure Functions
---------------------
The PyStructure class contains a set of functions that help handle and extract the data stored in the dictionary.
You can run the functions using ``database.<function>()``.

.. function:: database.get_coordinates(center = None)

   A function extracting the rightasencion and declination coordinates. If center coordinate is provided, the ra and dec coordimnates are returned as offset in arcsec.

   :param center: optional reference coordinate (e.g., ``"13:29:52.7 47:11:43"``). If provided, the returned values will represent the offset in arccsec with respect to this coordinate.
   :type center: str
   :return: ``ra``, ``dec`` ; two 1D arrays, one for the rightasencion and on for the declination.
   :rtype: np.array


.. function:: database.quickplot_2Dmap(line, s = 50, cmap = None)

    A function that generates a plot showing the integarted intensities of a user defined line.

    :param line: Name of the line in the PyStructure, e.g. "12CO21".
    :type line: str
    :param s: Marker size. User needs to vary if points do overlap or if there is too much space between the scatter points.
    :type s: int
    :param cmap: Colormap, default is "RdYlBu_r"
    :type cmap: int
    :return: 2D scatter plot illustrating the integrated intensities.

.. function:: database.get_vaxis(get_shuff = False)

    A function extracting the spectral velocities

    :param get_shuff: If ``True``, return the shuffled spectral axis.
    :type get_shuff: bool
    :return: ``vaxis`` ; 1D arrays with the spectral axis values.
    :rtype: array

.. function:: database.get_ratio(line,sn = 5)

    A function computing the line ratio between two lines

    :param line: A list of two string of line1 and line2, with the ratio being line1/line2 (e.g. ``["12CO21","12CO10"]``)
    :type line: list
    :param sn: signal-to-noise ratio used for sigma clipping (can be ``float`` as well).
    :type sn: int
    :return: ``ratio`` ; Dictionary with the ratio (extract using ``ratio['ratio']``).
    :rtype: dict

.. function:: database.export_fits(data_array,fname,adjust_header=None, verbose=False)

    A function that exports a 2D PyStructure map back to a FITS file

    :param data_array: The data array that will be exported back to a FITS file.
    :type data_array: array
    :param fname: Name of the FITS file to be saved.
    :type fname: str
    :param adjust_header: Dictionary with header keys and the corresponding value.
    :type adjust_header: dict
    :param verbose: If ``True``, print progress outpur.
    :type verbose: bool
    :return: Saves the FITS file on the disk
