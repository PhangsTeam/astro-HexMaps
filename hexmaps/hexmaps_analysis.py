"""
HexMapsAnalysis: helper class for loading and analysing HexMaps .ecsv files.

Usage
-----
    from hexmaps import HexMapsAnalysis

    hm = HexMapsAnalysis("output/ngc5194_hexmaps_27p0as_2025_01_01.ecsv")
    hm.quickplot_map("12CO21")
    hm.quickplot_spectrum("12CO21")
    hm.quickplot_shuffled_spectrum("12CO21")
"""

__author__ = "J. den Brok & L. Neumann"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import astropy.units as au
from astropy.coordinates import SkyCoord, FK5, Galactic
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
    rgal   : np.ndarray or None — deprojected galactocentric radius [kpc]
    theta  : np.ndarray or None — polar angle [rad]
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
            if key.startswith("SHUFF_"):
                lines.append(key[len("SHUFF_"):])
            elif key.startswith("MOM0_"):
                lines.append(key.split("_", 1)[1])
        lines = list(dict.fromkeys(lines))  # preserve order, drop duplicates
        return lines

    def _coord_cols(self):
        """
        Return the names of the two spatial coordinate columns.

        Works for any coordinate system the overlay cube used:
        equatorial (RA/DEC), galactic (GLON/GLAT), ecliptic, etc.
        Falls back to ('RA', 'DEC') if nothing else is found.
        """
        # Prefer columns explicitly named after known coordinate systems
        _candidates = [
            ("RA",   "DEC"),
            ("GLON", "GLAT"),
            ("ELON", "ELAT"),
        ]
        for c1, c2 in _candidates:
            if c1 in self.struct.colnames and c2 in self.struct.colnames:
                return c1, c2

        # Fall back to any pair of *_deg columns (excluding incl/posang)
        _skip = {"incl_deg", "posang_deg", "INCL_DEG", "POSANG_DEG"}
        deg_cols = [c for c in self.struct.colnames
                    if c.lower().endswith("_deg") and c not in _skip]
        if len(deg_cols) >= 2:
            return deg_cols[0], deg_cols[1]

        return "RA", "DEC"

    def _is_ra_dec(self) -> bool:
        """Return True if the spatial axes are equatorial (RA/DEC)."""
        c1, _ = self._coord_cols()
        return c1.upper() in ("RA",)

    def _is_galactic(self) -> bool:
        """Return True if the spatial axes are galactic (GLON/GLAT)."""
        c1, _ = self._coord_cols()
        return c1.upper() in ("GLON", "L")

    def _get_vaxis(self, shuffled: bool = False) -> au.Quantity:
        """Return the velocity axis for the first line."""
        key = "SPEC_VAXIS_SHUFF" if shuffled else "SPEC_VAXIS"
        return self.struct[key][0]

    def _centre_pixel(self) -> int:
        """Return the index of the sampling point closest to the target centre."""
        if self.rgal is not None:
            return int(np.argmin(self.rgal))
        # Fallback: sightline closest to the bounding-box centroid
        col1, col2 = self._coord_cols()
        x = np.array(self.struct[col1])
        y = np.array(self.struct[col2])
        cx, cy = np.nanmean(x), np.nanmean(y)
        return int(np.argmin((x - cx)**2 + (y - cy)**2))

    # ------------------------------------------------------------------
    # Coordinate helpers
    # ------------------------------------------------------------------

    def get_coordinates(self, center: str = None):
        """
        Return sky coordinate arrays, or offset coordinates relative to *center*.

        Works for equatorial (RA/DEC) and galactic (GLON/GLAT) coordinate
        systems.

        Parameters
        ----------
        center : str, optional
            Sky coordinate string.  For equatorial coordinates use the
            sexagesimal form ``"13:29:52.7 +47:11:43"``; for galactic
            coordinates use ``"202.47 +47.19"`` (decimal degrees).
            If None, returns absolute coordinates in degrees.

        Returns
        -------
        (axis1, axis2) : pair of np.ndarray
            Absolute coordinates in degrees (if center is None), or offsets
            in arcseconds relative to *center*.
        """
        col1, col2 = self._coord_cols()
        x = np.array(self.struct[col1])
        y = np.array(self.struct[col2])

        if center is None:
            return x, y

        if self._is_galactic():
            # Galactic offsets
            parts = center.split()
            cen_l = float(parts[0])
            cen_b = float(parts[1])
            ref = SkyCoord(l=cen_l * au.deg, b=cen_b * au.deg, frame=Galactic)
            pts = SkyCoord(l=x * au.deg, b=y * au.deg, frame=Galactic)
            aframe = ref.skyoffset_frame()
            off = pts.transform_to(aframe)
            return off.lon.arcsec, off.lat.arcsec
        else:
            # Equatorial offsets
            ref = SkyCoord(center, frame=FK5, unit=(au.hourangle, au.deg))
            pts = SkyCoord(ra=x * au.deg, dec=y * au.deg, frame=FK5)
            aframe = ref.skyoffset_frame()
            off = pts.transform_to(aframe)
            return off.lon.arcsec, off.lat.arcsec

    def sky_aspect_factor(self) -> float:
        """
        Aspect-ratio correction for plotting axis1 vs axis2 in degrees.

        For RA/Dec: ``1 / cos(Dec_median)`` so that one degree of RA and
        one degree of Dec have the same physical size on screen.
        For Galactic or other latitude/longitude systems the same formula
        applies to the latitude axis.

        Returns a value suitable for ``ax.set_aspect(factor)``.
        """
        _, col2 = self._coord_cols()
        lat = np.array(self.struct[col2])
        lat0 = np.nanmedian(lat)
        return 1.0 / np.cos(np.deg2rad(lat0))

    def _axis_labels(self, center: str = None):
        """
        Return (xlabel, ylabel) strings for the coordinate axes.

        When *center* is given the labels use arcsecond offsets.
        Otherwise they reflect the native coordinate system.
        """
        if center is not None:
            col1, col2 = self._coord_cols()
            name1 = col1.replace("_deg", "").replace("_DEG", "")
            name2 = col2.replace("_deg", "").replace("_DEG", "")
            return (rf"$\Delta${name1} [arcsec]", rf"$\Delta${name2} [arcsec]")

        col1, col2 = self._coord_cols()
        _label_map = {
            "RA":   "R.A. [deg]",
            "DEC":  "Dec. [deg]",
            "GLON": "Galactic longitude [deg]",
            "GLAT": "Galactic latitude [deg]",
            "ELON": "Ecliptic longitude [deg]",
            "ELAT": "Ecliptic latitude [deg]",
        }
        xlabel = _label_map.get(col1.upper(), f"{col1} [deg]")
        ylabel = _label_map.get(col2.upper(), f"{col2} [deg]")
        return xlabel, ylabel

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

        Works for any coordinate system (equatorial or galactic).

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
        else:
            col1, col2 = self._coord_cols()
            x = np.array(self.struct[col1])
            y = np.array(self.struct[col2])

        xlabel, ylabel = self._axis_labels(center)

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

        # Invert x-axis only for RA (east is left); galactic longitude
        # increases to the left too, so invert there as well.
        # For generic coordinate systems without a clear convention, invert
        # only when the first axis name suggests a "right ascension"-like axis.
        col1, _ = self._coord_cols()
        if col1.upper() in ("RA", "GLON", "L", "ELON"):
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
                "(RGAL_KPC), which is not available for this target. "
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
        uc    = np.full_like(i1, np.nan)
        ulim  = np.full_like(i1, np.nan)
        llim  = np.full_like(i1, np.nan)

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
        Return (and optionally save) a version of the table with all
        SPEC_ columns removed.
        """
        tbl = self.struct.copy()
        tbl.remove_columns([c for c in tbl.colnames if c.startswith("SPEC_")])
        for key in [k for k in tbl.meta if k.startswith("SPEC_")]:
            del tbl.meta[key]
        if save:
            if fname is None:
                fname = f'{tbl.meta.get("Target", "target")}_hexmaps_2D.ecsv'
            tbl.write(fname, format="ascii.ecsv", overwrite=True)
            print(f"[INFO] 2D table saved to {fname}")
        return tbl

    def get_config(self, save_to: str = None) -> str:
        """
        Return the config.txt content embedded in the .ecsv at run time.

        Parameters
        ----------
        save_to : str, optional
            If given, write the config content to this file path.
        """
        raw = self.struct.meta.get("config_file", "")
        if not raw:
            print(
                "[WARNING] No config_file entry found in the database metadata. "
                "This file may have been produced by an older pipeline version."
            )
            return ""
        content = raw.replace("\\n", "\n")
        if save_to is not None:
            from pathlib import Path
            Path(save_to).write_text(content, encoding="utf-8")
            print(f"[INFO] Config file written to {save_to}")
        return content

    def get_log(self, save_to: str = None) -> str:
        """
        Return the pipeline log embedded in the .ecsv at run time.

        Parameters
        ----------
        save_to : str, optional
            If given, write the log content to this file path.
        """
        raw = self.struct.meta.get("pipeline_log", "")
        if not raw:
            print(
                "[WARNING] No pipeline_log entry found in the database metadata. "
                "This file may have been produced by an older pipeline version."
            )
            return ""
        content = raw.replace("\\n", "\n")
        if save_to is not None:
            from pathlib import Path
            Path(save_to).write_text(content, encoding="utf-8")
            print(f"[INFO] Pipeline log written to {save_to}")
        return content

    def list_input_headers(self) -> list:
        """
        Return the labels of all input FITS headers embedded in this database.

        Returns
        -------
        list of str — sorted list of label strings.
        """
        prefix = "input_header_"
        return sorted(
            k[len(prefix):] for k in self.struct.meta if k.startswith(prefix)
        )

    def get_input_header(self, label: str):
        """
        Return the raw FITS header for the input file identified by *label*.

        Parameters
        ----------
        label : str — as returned by ``list_input_headers()``.

        Returns
        -------
        astropy.io.fits.Header
        """
        from astropy.io import fits as _fits

        if label.startswith("input_header_"):
            key = label
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
            f"HexMapsAnalysis(target='{self.struct.meta.get('Target', '?')}', "
            f"n_pts={len(self.struct)}, lines={self.lines})"
        )
