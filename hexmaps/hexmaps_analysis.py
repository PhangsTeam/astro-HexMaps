"""
HexMapsAnalysis: helper class for loading and analysing HexMaps .ecsv files.

Usage
-----
    from hexmaps import HexMapsAnalysis

    hm = HexMapsAnalysis("output/ngc5194_hexmaps_27p0as_2025_01_01.ecsv")
    hm.quickplot_map("12CO21")
    hm.quickplot_spectrum("12CO21")
    hm.quickplot_shuffled_spectrum("12CO21")

    # Binned trends with uncertainties (see hexmaps.utils_binning)
    profile = hm.bin_trend("RGAL_KPC", "12CO21", nbins=8)
    ratio = hm.bin_ratio("RGAL_KPC", "12CO21", "12CO10", nbins=8)
    hm.quickplot_binned(ratio)
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
        self.rgal = (
            np.array(self.struct["RGAL_KPC"])
            if "RGAL_KPC" in self.struct.colnames
            else None
        )
        self.theta = (
            np.array(self.struct["THETA_RAD"]) + np.pi
            if "THETA_RAD" in self.struct.colnames
            else None
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _find_lines(self) -> list:
        """Return line names from 3D or 2D columns."""
        lines = []
        for key in self.struct.keys():
            if key.startswith("SHUFF_"):
                lines.append(key[len("SHUFF_") :])
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
            ("RA", "DEC"),
            ("GLON", "GLAT"),
            ("ELON", "ELAT"),
        ]
        for c1, c2 in _candidates:
            if c1 in self.struct.colnames and c2 in self.struct.colnames:
                return c1, c2

        # Fall back to any pair of *_deg columns (excluding incl/posang)
        _skip = {"incl_deg", "posang_deg", "INCL_DEG", "POSANG_DEG"}
        deg_cols = [
            c
            for c in self.struct.colnames
            if c.lower().endswith("_deg") and c not in _skip
        ]
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
        return int(np.argmin((x - cx) ** 2 + (y - cy) ** 2))

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
            "RA": "R.A. [deg]",
            "DEC": "Dec. [deg]",
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
            + (
                f"  (r$_{{\\rm gal}}$ = {self.rgal[idx]:.2f} kpc)"
                if self.rgal is not None
                else ""
            )
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

    # ------------------------------------------------------------------
    # Binned trends
    # ------------------------------------------------------------------

    def _error_column(self, col: str):
        """Name of the uncertainty column belonging to *col*, or None."""
        candidates = [f"E{col}"]  # MOM0_X -> EMOM0_X, MAP_X -> EMAP_X, EW_X -> EEW_X
        if col.startswith("TPEAK_"):
            candidates.append("RMS_" + col[len("TPEAK_") :])  # peak error = channel rms
        return next((c for c in candidates if c in self.struct.colnames), None)

    def _resolve_quantity(self, quantity, errors="auto"):
        """
        Resolve a quantity for binning.

        *quantity* is a column name (e.g. ``"RGAL_KPC"``, ``"TPEAK_12CO21"``),
        a line name meaning its integrated intensity (``"12co21"`` ->
        ``MOM0_12CO21``), or an array. *errors* is ``"auto"`` (look up the
        uncertainty column of a named quantity), ``None``, a column name or
        an array.

        Returns (values, errors or None, column name or None).
        """
        col = None
        if isinstance(quantity, str):
            names = {c.upper(): c for c in self.struct.colnames}
            key = quantity.upper()
            if key in names:
                col = names[key]
            elif key in (ln.upper() for ln in self.lines) and f"MOM0_{key}" in names:
                col = names[f"MOM0_{key}"]
            else:
                raise KeyError(
                    f"{quantity!r} is neither a column nor a line name. "
                    f"Lines: {self.lines}"
                )
            values = self.struct[col]
        else:
            values = quantity

        if isinstance(errors, str) and errors == "auto":
            err_col = self._error_column(col) if col is not None else None
            errors = self.struct[err_col] if err_col else None
        elif isinstance(errors, str):
            errors = self._resolve_quantity(errors, errors=None)[0]
        return values, errors, col

    def _sampling_spacing_as(self) -> float:
        """Median distance [arcsec] between neighbouring sampling points."""
        from scipy.spatial import cKDTree

        col1, col2 = self._coord_cols()
        lon = np.asarray(self.struct[col1], dtype=float)
        lat = np.asarray(self.struct[col2], dtype=float)
        lat0 = np.nanmedian(lat)
        dlon = (lon - np.nanmedian(lon) + 180.0) % 360.0 - 180.0  # safe across 0/360
        xy = np.column_stack(
            [dlon * np.cos(np.deg2rad(lat0)) * 3600.0, (lat - lat0) * 3600.0]
        )
        xy = xy[np.all(np.isfinite(xy), axis=1)]
        if len(xy) < 2:
            return np.nan
        dist, _ = cKDTree(xy).query(xy, k=2)
        return float(np.median(dist[:, 1]))

    def oversampling_factor(self) -> float:
        """
        Number of sampling points per independent resolution element.

        Beam area divided by the area per hexagonal sampling point,
        ``(pi / (4 ln 2)) * beam**2 / ((sqrt(3) / 2) * spacing**2)``, with the
        beam FWHM from the table metadata (``beam_as``) and the spacing
        measured as the median nearest-neighbour distance of the points.
        For the default half-beam sampling (``pixels_per_beam = 2``) this is
        about 5.2. Neighbouring points are correlated, so a bin of *n*
        points contains only ``n / oversampling_factor`` independent
        measurements; the binning methods use this by default.

        Returns 1.0 (with a warning) if the beam size is not in the metadata.
        """
        beam = self.struct.meta.get("beam_as")
        spacing = self._sampling_spacing_as()
        if beam is None or not np.isfinite(spacing) or spacing <= 0:
            print(
                "[WARNING] Cannot determine the oversampling (no beam_as in the "
                "metadata or no sampling points); using 1."
            )
            return 1.0
        beam = beam.to_value(au.arcsec) if hasattr(beam, "unit") else float(beam)
        beam_area = np.pi / (4 * np.log(2)) * beam**2
        point_area = np.sqrt(3) / 2 * spacing**2
        return float(beam_area / point_area)

    def _resolve_extra(self, extra):
        """Resolve ``extra`` (names, or a dict name -> array/(array, errors))."""
        if extra is None:
            return None
        if isinstance(extra, dict):
            return extra
        out = {}
        for name in extra:
            values, errors, col = self._resolve_quantity(name)
            out[col or name] = (values, errors)
        return out

    def bin_trend(
        self,
        x,
        y,
        bins=None,
        nbins: int = 10,
        logbins: bool = False,
        method: str = "median",
        oversample="auto",
        x_err="auto",
        y_err="auto",
        extra=None,
        keep_empty: bool = False,
    ):
        """
        Bin a quantity *y* in bins of *x*, with uncertainties.

        For example the CO intensity as a function of galactocentric radius::

            tab = hm.bin_trend("RGAL_KPC", "12co21", nbins=8)
            hm.quickplot_binned(tab)

        Parameters
        ----------
        x, y : str or array-like
            Column name (``"RGAL_KPC"``, ``"TPEAK_12CO21"``, ``"MAP_SPIRE250"``,
            ...), line name for its integrated intensity (``"12co21"`` ->
            ``MOM0_12CO21``), or an array with one value per sampling point.
        bins : array-like, optional
            Explicit bin edges. Otherwise *nbins* bins span the range of *x*.
        nbins : int
            Number of bins if *bins* is not given (default 10).
        logbins : bool
            Logarithmically spaced bins (over the positive *x* values).
        method : ``"median"`` (default) or ``"mean"``
            Statistic used for the bin value.
        oversample : float or ``"auto"``
            Sampling points per independent beam, used for the bin
            uncertainties. ``"auto"`` (default) uses
            :meth:`oversampling_factor`; pass 1 for independent points.
        x_err, y_err : ``"auto"``, None, str or array-like
            Per-point uncertainties. ``"auto"`` uses the matching error
            column (``EMOM0_*`` for ``MOM0_*``, ``EMAP_*`` for ``MAP_*``,
            ``RMS_*`` for ``TPEAK_*``) if there is one.
        extra : list of str, or dict, optional
            More quantities to average in the same bins, e.g.
            ``["MAP_SPIRE250", "MOM2_12CO21"]`` (errors looked up
            automatically) or ``{name: values}`` / ``{name: (values, errors)}``.
        keep_empty : bool
            Keep rows (filled with NaN) for bins without data points.

        Returns
        -------
        astropy.table.Table
            One row per bin: ``bin_lo``, ``bin_hi``, ``n``, ``x``, ``x_err``,
            ``y``, ``y_err``, ``snr``, ``y_p16``, ``y_p84`` (16th/84th
            percentile of the values in the bin), ``n_pos``, ``n_neg``, and
            ``<name>``, ``<name>_err`` for *extra* quantities. See
            :mod:`hexmaps.utils_binning` for the definitions.
        """
        from hexmaps.utils_binning import bin_trend

        xv, xe, xcol = self._resolve_quantity(x, x_err)
        yv, ye, ycol = self._resolve_quantity(y, y_err)
        if isinstance(oversample, str) and oversample == "auto":
            oversample = self.oversampling_factor()
        tab = bin_trend(
            xv, yv, y_err=ye, x_err=xe, bins=bins, nbins=nbins, logbins=logbins,
            method=method, oversample=oversample, extra=self._resolve_extra(extra),
            keep_empty=keep_empty,
        )
        tab.meta.update(x_name=xcol or "x", y_name=ycol or "y")
        return tab

    def bin_ratio(
        self,
        x,
        numerator,
        denominator,
        bins=None,
        nbins: int = 10,
        logbins: bool = False,
        method: str = "median",
        oversample="auto",
        snr_limit: float = 3.0,
        x_err="auto",
        numerator_err="auto",
        denominator_err="auto",
        extra=None,
        keep_empty: bool = False,
    ):
        """
        Bin the ratio ``numerator / denominator`` in bins of *x*.

        Numerator and denominator are binned separately and divided
        ("stacked" ratio), so the ratio stays defined in bins where
        individual sightlines are not detected. Where the numerator
        (denominator) is not detected in a bin, an upper (lower) limit is
        given instead. For example a line ratio versus radius::

            tab = hm.bin_ratio("RGAL_KPC", "12co21", "12co10", nbins=8)
            hm.quickplot_binned(tab)

        Parameters
        ----------
        x, numerator, denominator : str or array-like
            As *x* / *y* in :meth:`bin_trend`; line names select ``MOM0``.
        snr_limit : float
            S/N threshold for the limits (default 3).
        numerator_err, denominator_err : ``"auto"``, str or array-like
            Per-point uncertainties (required; ``"auto"`` uses ``EMOM0_*``
            etc.).
        bins, nbins, logbins, method, oversample, x_err, extra, keep_empty
            As in :meth:`bin_trend`.

        Returns
        -------
        astropy.table.Table
            The columns of :meth:`bin_trend` (``y`` is the binned ratio,
            ``y_p16``/``y_p84`` are percentiles of the per-point ratios),
            plus ``y1``, ``y1_err``, ``y2``, ``y2_err`` (binned numerator and
            denominator), ``upper_limit`` and ``lower_limit`` (NaN where the
            ratio is significant). See :mod:`hexmaps.utils_binning`.
        """
        from hexmaps.utils_binning import bin_ratio_trend

        xv, xe, xcol = self._resolve_quantity(x, x_err)
        y1, e1, c1 = self._resolve_quantity(numerator, numerator_err)
        y2, e2, c2 = self._resolve_quantity(denominator, denominator_err)
        for name, err in (("numerator", e1), ("denominator", e2)):
            if err is None:
                raise ValueError(
                    f"No uncertainties found for the {name}; pass {name}_err."
                )
        if isinstance(oversample, str) and oversample == "auto":
            oversample = self.oversampling_factor()
        tab = bin_ratio_trend(
            xv, y1, y2, e1, e2, x_err=xe, bins=bins, nbins=nbins, logbins=logbins,
            method=method, oversample=oversample, snr_limit=snr_limit,
            extra=self._resolve_extra(extra), keep_empty=keep_empty,
        )
        tab.meta.update(
            x_name=xcol or "x",
            y_name=f"{c1 or 'y1'} / {c2 or 'y2'}",
            y1_name=c1,
            y2_name=c2,
        )
        return tab

    def quickplot_binned(
        self,
        binned,
        ax=None,
        show_points: bool = True,
        show_scatter: bool = True,
        logx: bool = False,
        logy: bool = False,
        color: str = "C0",
        label: str = None,
    ):
        """
        Plot a binned trend from :meth:`bin_trend` or :meth:`bin_ratio`.

        Draws the bin values with their uncertainties, the 16th-84th
        percentile range of the individual values as a shaded band, and,
        for ratios, upper (``v``) and lower (``^``) limits in place of
        non-significant bins.

        Parameters
        ----------
        binned       : astropy.table.Table — result of a binning method
        ax           : Axes — existing axes to draw into (new figure if None)
        show_points  : bool — show the individual sampling points in grey
                       (when *x* and *y* were given as column/line names; for
                       ratios only points where both lines exceed
                       ``snr_limit``)
        show_scatter : bool — shade the 16th-84th percentile range
        logx, logy   : bool — logarithmic axes
        color, label : matplotlib colour and legend label for the bins

        Returns
        -------
        matplotlib.axes.Axes
        """
        meta = binned.meta
        own_fig = ax is None
        if own_fig:
            fig, ax = plt.subplots(figsize=(6, 4.5))

        cols = self.struct.colnames
        if show_points and meta.get("x_name") in cols:
            px = np.asarray(self.struct[meta["x_name"]], dtype=float)
            py = None
            if meta.get("y_name") in cols:
                py = np.asarray(self.struct[meta["y_name"]], dtype=float)
            elif meta.get("y1_name") in cols and meta.get("y2_name") in cols:
                y1, e1, _ = self._resolve_quantity(meta["y1_name"])
                y2, e2, _ = self._resolve_quantity(meta["y2_name"])
                y1, y2 = np.asarray(y1, dtype=float), np.asarray(y2, dtype=float)
                lim = meta.get("snr_limit", 3.0)
                with np.errstate(divide="ignore", invalid="ignore"):
                    det = (y1 / np.asarray(e1, dtype=float) > lim) & (
                        y2 / np.asarray(e2, dtype=float) > lim
                    )
                    py = np.where(det, y1 / y2, np.nan)
            if py is not None:
                ax.scatter(px, py, s=6, color="0.7", alpha=0.6, lw=0, zorder=1,
                           label="sampling points")

        def _col(name):
            if name not in binned.colnames:
                return np.full(len(binned), np.nan)
            return np.asarray(binned[name], dtype=float)

        x, y = _col("x"), _col("y")
        p16, p84 = _col("y_p16"), _col("y_p84")
        upper, lower = _col("upper_limit"), _col("lower_limit")
        # Values that cannot be drawn on a logarithmic axis are left out
        # (e.g. bins whose median is 0 because the sightlines are undetected).
        drawable = np.isfinite(x) & np.isfinite(y)
        if logx:
            drawable &= x > 0
        if logy:
            drawable &= y > 0
            p16 = np.where(p16 > 0, p16, np.nan)

        if show_scatter:
            band = np.isfinite(x) & np.isfinite(p16) & np.isfinite(p84)
            ax.fill_between(x[band], p16[band], p84[band], color=color,
                            alpha=0.2, lw=0, zorder=2,
                            label="16-84th percentile" if label is None else None)

        det = drawable & ~(np.isfinite(upper) | np.isfinite(lower))
        ax.errorbar(
            x[det], y[det], xerr=_col("x_err")[det], yerr=_col("y_err")[det],
            fmt="o", color=color, mec="k", mew=0.6, ms=6, capsize=2, zorder=3,
            label=label or f"binned {meta.get('method', '')}".strip(),
        )
        if np.isfinite(upper).any():
            ax.scatter(x, upper, marker="v", color=color, edgecolor="k", lw=0.6,
                       s=45, zorder=3, label="upper limit")
        if np.isfinite(lower).any():
            ax.scatter(x, lower, marker="^", color=color, edgecolor="k", lw=0.6,
                       s=45, zorder=3, label="lower limit")

        def _label(name, unit):
            return f"{name} [{unit}]" if unit is not None else name

        ax.set_xlabel(_label(meta.get("x_name", "x"), binned["x"].unit))
        ax.set_ylabel(_label(meta.get("y_name", "y"), binned["y"].unit))
        from matplotlib.ticker import FuncFormatter, LogLocator

        for axis, use_log, get_lim in (
            (ax.xaxis, logx, ax.get_xlim), (ax.yaxis, logy, ax.get_ylim)
        ):
            if not use_log:
                continue
            (ax.set_xscale if axis is ax.xaxis else ax.set_yscale)("log")
            lo, hi = get_lim()
            # Short log ranges: label 1, 2, 5 x 10^n as plain numbers instead of
            # matplotlib's crowded minor-tick labels.
            if lo > 0 and np.log10(hi / lo) <= 2.5:
                plain = FuncFormatter(lambda v, _: f"{v:g}")
                axis.set_major_formatter(plain)
                axis.set_minor_locator(LogLocator(subs=(2, 5)))
                axis.set_minor_formatter(plain)
        ax.legend(fontsize=8)

        if own_fig:
            plt.tight_layout()
            plt.show()
        return ax

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
            k[len(prefix) :] for k in self.struct.meta if k.startswith(prefix)
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
