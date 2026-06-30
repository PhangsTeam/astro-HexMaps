"""
HexMapsAnalysis: helper class for loading and analysing HexMaps .ecsv files.

Usage
-----
    from hexmaps_analysis import HexMapsAnalysis

    hm = HexMapsAnalysis("Output/ngc5194_hexmaps_27as_2025_01_01.ecsv")
    hm.quickplot_map("12CO21")
    hm.quickplot_spectrum("12CO21")
    hm.quickplot_shuffled_spectrum("12CO21")
"""

__author__ = "J. den Brok & L. Neumann"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import astropy.units as au
from astropy.coordinates import SkyCoord, FK5
from astropy.table import Table


class HexMapsAnalysis:
    """
    Load and analyse a HexMaps .ecsv database.

    Parameters
    ----------
    path : str
        Path to the .ecsv file produced by the pipeline.

    Attributes
    ----------
    struct : astropy.table.Table
    lines  : list of str — spectral line names found in the database
    rgal   : np.ndarray  — deprojected galactocentric radius [kpc]
    theta  : np.ndarray  — polar angle [rad]
    """

    def __init__(self, path: str):
        self.struct = Table.read(path)
        self.lines = self._find_lines()

        # rgal and theta are only available for galaxy targets
        self.rgal  = (np.array(self.struct["RGAL_KPC"])
                      if "RGAL_KPC" in self.struct.colnames else None)
        self.theta = (np.array(self.struct["THETA_RAD"]) + np.pi
                      if "THETA_RAD" in self.struct.colnames else None)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _find_lines(self) -> list:
        """Return line names from 3D or 2D columns."""
        lines = []
        for key in self.struct.keys():
            if key.startswith("SHUFF"):
                lines.append(key[len("SHUFF"):])
            elif key.startswith("MOM0"):
                lines.append(key.split("_")[1])
        lines = list(dict.fromkeys(lines))  # only keep unique elements
        return lines

    def _coord_cols(self):
        """
        Return the names of the two spatial coordinate columns.

        Column names depend on the WCS of the overlay cube used to produce this
        database (e.g. ``RA``/``DEC`` for equatorial,
        ``GLON``/``GLAT`` for galactic).  Falls back to
        ``RA``/``DEC`` if no matching columns are found.
        """
        _skip = {"incl_deg", "posang_deg"}
        cols = [c for c in self.struct.colnames
                if c.endswith("_deg") and c not in _skip]
        if len(cols) >= 2:
            return cols[0], cols[1]
        return "RA", "DEC"

    def _get_vaxis(self, shuffled: bool = False) -> au.Quantity:
        """Return the velocity axis for the first line."""
        if shuffled:
            return self.struct["SPEC_VAXIS_SHUFF"][0]
        else:
            return self.struct["SPEC_VAXIS"][0]

    def _centre_pixel(self) -> int:
        """Return the index of the sampling point closest to the source centre."""
        if self.rgal is not None:
            return int(np.argmin(self.rgal))
        # Fallback: centre of the bounding box
        col1, col2 = self._coord_cols()
        x = np.array(self.struct[col1])
        y = np.array(self.struct[col2])
        cx, cy = np.nanmean(x), np.nanmean(y)
        return int(np.argmin((x - cx)**2 + (y - cy)**2))

    def _get_vaxis(self, shuffled: bool = False) -> au.Quantity:
        """Return the velocity axis for the first line."""
        if shuffled:
            return self.struct["SPEC_VAXIS_SHUFF"][0]
        else:
            return self.struct["SPEC_VAXIS"][0]

    def _centre_pixel(self) -> int:
        """Return the index of the sampling point closest to the source centre."""
        return int(np.argmin(self.rgal))

    # ------------------------------------------------------------------
    # Coordinate helpers
    # ------------------------------------------------------------------

    def get_coordinates(self, center: str = None):
        """
        Return sky coordinate arrays, or offset coordinates relative to *center*.

        Parameters
        ----------
        center : str, optional
            Sky coordinate string, e.g. ``"13:29:52.7 47:11:43"``.
            If given, returns (delta_axis1, delta_axis2) in arcseconds.
            If None, returns absolute coordinates in degrees.
        """
        col1, col2 = self._coord_cols()
        ra = np.array(self.struct[col1])
        dec = np.array(self.struct[col2])

        if center is None:
            return ra, dec

        ref = SkyCoord(center, frame=FK5, unit=(au.hourangle, au.deg))
        pts = SkyCoord(ra=ra * au.deg, dec=dec * au.deg, frame=FK5)
        aframe = ref.skyoffset_frame()
        delta_ra = pts.transform_to(aframe).lon.arcsec
        delta_dec = pts.transform_to(aframe).lat.arcsec
        return delta_ra, delta_dec

    def sky_aspect_factor(self):
        """
        Aspect ratio correction for plotting axis 1 vs axis 2 in degrees.

        Returns a value suitable for ``ax.set_aspect(factor)``.
        For RA/Dec this is 1/cos(Dec); for galactic or ecliptic coordinates
        the same formula applies to the latitude axis.
        """
        _, col2 = self._coord_cols()
        lat = np.array(self.struct[col2])
        lat0 = np.median(lat)
        return 1.0 / np.cos(np.deg2rad(lat0))

    # ------------------------------------------------------------------
    # Quick-look plots
    # ------------------------------------------------------------------

    def quickplot_map(
        self,
        line: str,
        quantity: str = "MOM0",
        s: int = 100,
        cmap: str = "viridis",
        stretch: str = "lin",
        center: str = None,
        ax=None,
    ):
        """
        Scatter-plot a 2D moment map on the hexagonal grid.

        Parameters
        ----------
        line     : str   — line name, e.g. ``"12CO21"``
        quantity : str   — column prefix: ``"MOM0"``, ``"MOM1"``, ``"MOM2"``,
                           ``"TPEAK"``, ``"RMS"``, or ``"MAP"`` for a 2D map
        s        : int   — marker size
        cmap     : str   — matplotlib colormap name
        stretch  : str   — ``"lin"``, ``"log"``, or ``"symlog"``
        center   : str   — if given, plot offset coords (arcsec)
        ax       : Axes  — existing axes to draw into (creates new figure if None)
        """
        col = f"{quantity}_{line.upper()}"
        if col not in self.struct.colnames:
            raise KeyError(
                f"Column '{col}' not found. Available: {self.struct.colnames}"
            )

        values = np.array(self.struct[col])
        unit = str(self.struct[col].unit) if hasattr(self.struct[col], "unit") else ""

        if center is not None:
            x, y = self.get_coordinates(center)
            xlabel, ylabel = r"$\Delta$Axis1 [arcsec]", r"$\Delta$Axis2 [arcsec]"
        else:
            col1, col2 = self._coord_cols()
            x = np.array(self.struct[col1])
            y = np.array(self.struct[col2])
            xlabel = f"{col1.replace('_deg', '')} [deg]"
            ylabel = f"{col2.replace('_deg', '')} [deg]"

        finite = values[np.isfinite(values)]
        if len(finite) == 0:
            print(f"[WARNING] All values are NaN for {col}.")
            return

        if stretch == "lin":
            norm = mcolors.Normalize(vmin=np.nanmin(values), vmax=np.nanmax(values))
        elif stretch == "log":
            norm = mcolors.LogNorm(
                vmin=finite[finite > 0].min(), vmax=np.nanmax(values)
            )
        elif stretch == "symlog":
            norm = mcolors.SymLogNorm(
                linthresh=np.nanmax(np.abs(values)) * 0.1,
                vmin=np.nanmin(values),
                vmax=np.nanmax(values),
            )
        else:
            raise ValueError(
                f"Unknown stretch '{stretch}'. Use 'lin', 'log', or 'symlog'."
            )

        own_fig = ax is None
        if own_fig:
            fig, ax = plt.subplots(figsize=(6, 6))
        else:
            fig = ax.get_figure()

        im = ax.scatter(x, y, c=values, s=s, marker="h", cmap=cmap, norm=norm)
        ax.invert_xaxis()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{line}  {quantity}")

        if center is not None:
            ax.set_aspect("equal")
        else:
            ax.set_aspect(self.sky_aspect_factor())

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(f"{quantity} [{unit}]" if unit else quantity)

        if own_fig:
            plt.tight_layout()
            plt.show()

    def quickplot_spectrum(
        self,
        line: str,
        idx: int = None,
        show_mask: bool = True,
        show_rms: bool = True,
        ax=None,
    ):
        """
        Plot a single spectrum from the native velocity grid.

        Parameters
        ----------
        line      : str  — line name, e.g. ``"12CO21"``
        idx       : int  — sampling-point index (default: point closest to centre)
        show_mask : bool — shade the integration mask region
        show_rms  : bool — show the RMS level
        ax        : Axes — existing axes to draw into
        """
        col = f"SPEC_{line.upper()}"
        if col not in self.struct.colnames:
            raise KeyError(f"Column '{col}' not found.")

        if idx is None:
            idx = self._centre_pixel()

        vaxis = self._get_vaxis(shuffled=False)
        spec = np.array(self.struct[col][idx])
        unit = str(self.struct[col].unit) if hasattr(self.struct[col], "unit") else ""
        vunit = str(vaxis.unit) if hasattr(vaxis, "unit") else "km/s"
        vvals = vaxis.value if hasattr(vaxis, "value") else np.array(vaxis)

        own_fig = ax is None
        if own_fig:
            fig, ax = plt.subplots(figsize=(7, 4))

        ax.step(vvals, spec, where="mid", color="steelblue", linewidth=1.5, label=line)

        if show_mask and "SPEC_MASK" in self.struct.colnames:
            mask = np.array(self.struct["SPEC_MASK"][idx])
            ylo, yhi = ax.get_ylim()
            ax.fill_between(
                vvals,
                ylo,
                yhi,
                where=(mask == 1),
                color="lightblue",
                alpha=0.4,
                label="mask",
            )
            ax.set_ylim(ylo, yhi)

        if show_rms and f"RMS_{line.upper()}" in self.struct.colnames:
            rms = self.struct[f"RMS_{line.upper()}"][idx]
            ax.axhline(rms, color="r", linewidth=0.8, linestyle=":", label="RMS")

        ax.axhline(0, color="k", linewidth=0.6, linestyle="--")
        ax.set_xlabel(f"Velocity [{vunit}]")
        ax.set_ylabel(f"T$_{{\\rm b}}$ [{unit}]" if unit else "Brightness temperature")
        ax.set_title(
            f"{line}  —  pixel {idx}"
            + (f"  (r$_{{\\rm gal}}$ = {self.rgal[idx]:.2f} kpc)"
               if self.rgal is not None else "")
        )
        ax.legend(fontsize=9)

        if own_fig:
            plt.tight_layout()
            plt.show()

    def quickplot_shuffled_spectrum(self, line: str, idx: int = None, ax=None):
        """
        Plot the velocity-shifted (shuffled) spectrum for a single pixel.

        Parameters
        ----------
        line : str  — line name, e.g. ``"12CO21"``
        idx  : int  — sampling-point index (default: point closest to centre)
        ax   : Axes — existing axes to draw into
        """
        col = f"SPEC_SHUFF_{line.upper()}"
        if col not in self.struct.colnames:
            raise KeyError(f"Column '{col}' not found.")

        if idx is None:
            idx = self._centre_pixel()

        vaxis = self._get_vaxis(shuffled=True)
        spec = np.array(self.struct[col][idx])
        unit = str(self.struct[col].unit) if hasattr(self.struct[col], "unit") else ""
        vvals = vaxis.value if hasattr(vaxis, "value") else np.array(vaxis)

        own_fig = ax is None
        if own_fig:
            fig, ax = plt.subplots(figsize=(7, 4))

        ax.step(vvals, spec, where="mid", color="darkorange", linewidth=1.5)
        ax.axhline(0, color="k", linewidth=0.6, linestyle="--")
        ax.axvline(0, color="grey", linewidth=0.8, linestyle=":")
        ax.set_xlabel("Velocity offset [km/s]")
        ax.set_ylabel(f"T$_{{\\rm b}}$ [{unit}]" if unit else "Brightness temperature")
        ax.set_title(f"{line} (shuffled)  —  pixel {idx}")

        if own_fig:
            plt.tight_layout()
            plt.show()

    def quickplot_radial_profile(
        self, line: str, quantity: str = "MOM0", nbins: int = 10, ax=None
    ):
        """
        Plot a radial profile (binned median) of a moment map.

        Parameters
        ----------
        line     : str — line name
        quantity : str — column prefix (``"MOM0"``, ``"TPEAK"``, etc.)
        nbins    : int — number of radial bins
        ax       : Axes
        """
        col = f"{quantity}_{line.upper()}"
        if col not in self.struct.colnames:
            raise KeyError(f"Column '{col}' not found.")

        if self.rgal is None:
            raise RuntimeError(
                "quickplot_radial_profile requires galaxy geometry "
                "(RGAL_KPC), which is not available for this source. "
                "Add incl_deg, posang_deg, and r25 to target_definitions.txt."
            )

        values = np.array(self.struct[col])
        unit = str(self.struct[col].unit) if hasattr(self.struct[col], "unit") else ""

        r_max = np.nanmax(self.rgal)
        edges = np.linspace(0, r_max, nbins + 1)
        r_cen = 0.5 * (edges[:-1] + edges[1:])
        medians = np.array(
            [
                np.nanmedian(
                    values[(self.rgal >= edges[i]) & (self.rgal < edges[i + 1])]
                )
                for i in range(nbins)
            ]
        )

        own_fig = ax is None
        if own_fig:
            fig, ax = plt.subplots(figsize=(6, 4))

        ax.plot(r_cen, medians, marker="o", color="steelblue", linewidth=1.5)
        ax.set_xlabel("Galactocentric radius [kpc]")
        ax.set_ylabel(f"Median {quantity} [{unit}]" if unit else f"Median {quantity}")
        ax.set_title(f"{line}  —  radial profile")

        if own_fig:
            plt.tight_layout()
            plt.show()

    # ------------------------------------------------------------------
    # Data access helpers
    # ------------------------------------------------------------------

    def get_mom0(self, line: str) -> np.ndarray:
        """Return the moment-0 array for *line*."""
        return np.array(self.struct[f"MOM0_{line.upper()}"])

    def get_ratio(self, line1: str, line2: str, sn: float = 5.0) -> dict:
        """
        Compute the line ratio line1 / line2 with upper/lower limits.

        Returns a dict with keys: ``ratio``, ``uc``, ``ulimit``, ``llimit``.
        """
        i1 = np.array(self.struct[f"MOM0_{line1.upper()}"])
        e1 = np.array(self.struct[f"EMOM0_{line1.upper()}"])
        i2 = np.array(self.struct[f"MOM0_{line2.upper()}"])
        e2 = np.array(self.struct[f"EMOM0_{line2.upper()}"])

        ratio = np.full_like(i1, np.nan)
        uc = np.full_like(i1, np.nan)
        ulim = np.full_like(i1, np.nan)
        llim = np.full_like(i1, np.nan)

        det = (i1 / e1 > sn) & (i2 / e2 > sn)
        ratio[det] = i1[det] / i2[det]
        uc[det] = ratio[det] * np.sqrt(
            (e1[det] / i1[det]) ** 2 + (e2[det] / i2[det]) ** 2
        )

        ul = (~det) & (i2 / e2 > sn)
        ulim[ul] = (2 / 3) * sn * e1[ul] / i2[ul]

        ll = (i1 / e1 > sn) & (~det)
        llim[ll] = i1[ll] / e2[ll] / ((2 / 3) * sn)

        return {"ratio": ratio, "uc": uc, "ulimit": ulim, "llimit": llim}

    def get_2D_database(self, fname: str = None, save: bool = False):
        """
        Return (and optionally save) a version of the table with all SPEC_ columns removed.
        """
        tbl = self.struct.copy()
        tbl.remove_columns([c for c in tbl.colnames if c.startswith("SPEC_")])
        for key in [k for k in tbl.meta if k.startswith("SPEC_")]:
            del tbl.meta[key]
        if save:
            if fname is None:
                fname = f'{tbl.meta.get("Source", "source")}_hexmaps_2D.ecsv'
            tbl.write(fname, format="ascii.ecsv", overwrite=True)
            print(f"[INFO] 2D table saved to {fname}")
        return tbl

    def get_config(self, save_to: str = None) -> str:
        """
        Return the config.txt content that was embedded in the .ecsv at
        pipeline run time.

        The content is stored in ``table.meta["config_file"]`` as a single
        line with newlines encoded as the two-character sequence ``\\n``.
        This method decodes that back to the original multi-line string.

        Parameters
        ----------
        save_to : str, optional
            If given, write the config content to this file path.  Useful
            for inspecting or reproducing a previous pipeline run.

        Returns
        -------
        str
            The full content of the config.txt that was used to produce
            this database, or an empty string if the metadata key is absent
            (e.g. files produced by an older pipeline version).

        Examples
        --------
        >>> hm = HexMapsAnalysis("ngc5194_hexmaps_27p0as_2025_01_01.ecsv")
        >>> print(hm.get_config())
        # HexMaps configuration file
        ...
        >>> hm.get_config(save_to="recovered_config.txt")
        """
        raw = self.struct.meta.get("config_file", "")
        if not raw:
            print(
                "[WARNING] No config_file entry found in the database metadata. "
                "This file may have been produced by an older pipeline version."
            )
            return ""

        # Decode the newline escape used when storing in the ECSV header
        content = raw.replace("\\n", "\n")

        if save_to is not None:
            from pathlib import Path

            Path(save_to).write_text(content, encoding="utf-8")
            print(f"[INFO] Config file written to {save_to}")

        return content

    def get_log(self, save_to: str = None) -> str:
        """
        Return the pipeline log that was embedded in the .ecsv at run time.

        The log is stored in ``table.meta["pipeline_log"]`` as a single line
        with newlines encoded as the two-character sequence ``\\n`` (same
        encoding as ``config_file``).  This method decodes it back to the
        original multi-line string.

        Parameters
        ----------
        save_to : str, optional
            If given, write the log content to this file path.

        Returns
        -------
        str
            The full pipeline log produced during the run that created this
            database, or an empty string if the metadata key is absent
            (e.g. files produced by an older pipeline version).

        Examples
        --------
        >>> hm = HexMapsAnalysis("ngc5194_hexmaps_27p0as_2025_01_01.ecsv")
        >>> print(hm.get_log())
        2025-01-01 12:00:00 [Loading]  [INFO]    Loading config file ...
        ...
        >>> hm.get_log(save_to="run.log")
        """
        raw = self.struct.meta.get("pipeline_log", "")
        if not raw:
            print(
                "[WARNING] No pipeline_log entry found in the database metadata. "
                "This file may have been produced by an older pipeline version."
            )
            return ""

        # Decode the newline escape used when storing in the ECSV header
        content = raw.replace("\\n", "\n")

        if save_to is not None:
            from pathlib import Path
            Path(save_to).write_text(content, encoding="utf-8")
            print(f"[INFO] Pipeline log written to {save_to}")

        return content

    def list_input_headers(self) -> list:
        """
        Return the labels of all input FITS headers embedded in this database.

        Labels match the table column keys used for each input, uppercased:

        - ``"OVERLAY"`` — the overlay cube
        - ``"<LINE_NAME>"`` — spectral cube (e.g. ``"12CO21"``)
        - ``"MAP_<LINE_NAME>"`` — 2D companion map for a cube (e.g. ``"MAP_12CO21"``)
        - ``"EMAP_<LINE_NAME>"`` — uncertainty map (e.g. ``"EMAP_12CO21"``)
        - ``"<MAP_NAME>"`` — standalone 2D map (e.g. ``"SPIRE250"``)
        - ``"SPEC_<MASK_NAME>"`` — external FITS mask (e.g. ``"SPEC_HEXMASK"``)

        Returns
        -------
        list of str — sorted list of label strings, or an empty list if
        no input headers were embedded (e.g. older pipeline version).

        Examples
        --------
        >>> hm = HexMapsAnalysis("ngc5194_hexmaps_27p0as_2025_01_01.ecsv")
        >>> hm.list_input_headers()
        ['12CO21', '12CO10', 'OVERLAY', 'SPIRE250']
        """
        prefix = "input_header_"
        return sorted(
            k[len(prefix) :] for k in self.struct.meta if k.startswith(prefix)
        )

    def get_input_header(self, label: str):
        """
        Return the raw FITS header for the input file identified by *label*.

        The header is stored in the .ecsv metadata as a compact 80-char-per-card
        string (FITS standard ``header.tostring()`` format). This method
        deserialises it back to an ``astropy.io.fits.Header`` object.

        Parameters
        ----------
        label : str
            Label matching the table column key, uppercased — as returned by
            ``list_input_headers()``.  Examples: ``"12CO21"``, ``"SPIRE250"``,
            ``"OVERLAY"``, ``"MAP_12CO21"``, ``"EMAP_SPIRE250"``.
            You can also pass the full metadata key
            (``"input_header_12CO21"``); the prefix will be stripped.

        Returns
        -------
        astropy.io.fits.Header
            The original FITS header before any pipeline processing.

        Raises
        ------
        KeyError if *label* is not found in the embedded headers.

        Examples
        --------
        >>> hm = HexMapsAnalysis("ngc5194_hexmaps_27p0as_2025_01_01.ecsv")
        >>> hdr = hm.get_input_header("12CO21")
        >>> print(hdr["BMAJ"] * 3600, "arcsec")
        12.82 arcsec
        >>> hdr_ov = hm.get_input_header("OVERLAY")
        >>> print(repr(hdr_ov))    # prints all header cards
        """
        from astropy.io import fits as _fits

        # Accept both the bare label and the full meta key
        if label.startswith("input_header_"):
            key = label
            label = label[len("input_header_") :]
        else:
            key = f"input_header_{label}"

        if key not in self.struct.meta:
            available = self.list_input_headers()
            raise KeyError(
                f"No embedded header found for label {label!r}. "
                f"Available: {available}"
            )

        return _fits.Header.fromstring(self.struct.meta[key])

    def __repr__(self):
        return (
            f"HexMapsAnalysis(source='{self.struct.meta.get('Source', '?')}', "
            f"n_pts={len(self.struct)}, lines={self.lines})"
        )
