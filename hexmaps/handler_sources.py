"""
handler_sources.py — SourceHandler: manages the source list and geometry.

Wraps the source geometry table loaded by KeyHandler and provides
named-parameter lookups used by every pipeline stage that needs
per-source geometry (position angle, inclination, distance, r25).
"""

import pandas as pd

from hexmaps.logger import get_logger

LOG = get_logger("Loading")


class SourceHandler:
    """
    Manages the source list and their geometric parameters.

    Parameters
    ----------
    source_table : pd.DataFrame
        Full geometry table from target_definitions.txt, as loaded by
        KeyHandler.  Must contain a 'source' column plus the geometry columns
        defined in handler_keys.TARGET_COLUMNS.
    sources : list of str
        Ordered list of source names to process.  Must be a subset of the
        names in source_table['source'].

    Raises
    ------
    ValueError
        If any name in *sources* is not present in *source_table*.

    Example
    -------
    >>> kh = KeyHandler("./keys/")
    >>> sh = SourceHandler(kh.get_source_table(), kh.get_sources())
    >>> params = sh.get_source_params("ngc5194")
    >>> print(params["dist_mpc"])
    """

    def __init__(self, source_table: pd.DataFrame, sources: list):
        self.source_table = source_table.copy()
        self.sources = list(sources)
        # Build a name → row-index lookup for O(1) access
        self._index = {row["source"]: idx for idx, row in source_table.iterrows()}
        self._validate()

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate(self):
        """
        Check that every requested source is present in the geometry table.

        Raises ValueError listing all missing names so the user can fix
        target_definitions.txt or config.txt in one go.
        """
        missing = [s for s in self.sources if s not in self._index]
        if missing:
            LOG.error(
                f"The following sources are not in "
                f"target_definitions.txt: {missing}"
            )
            raise ValueError(
                f"The following sources are not in "
                f"target_definitions.txt: {missing}"
            )

    # ------------------------------------------------------------------
    # Parameter lookups
    # ------------------------------------------------------------------

    def get_source_params(self, source: str) -> dict:
        """
        Return a dict of all geometric parameters for *source*.

        Keys
        ----
        ra_ctr, dec_ctr   — centre coordinates (degrees, in the WCS frame of
                            the overlay cube — may be RA/Dec, galactic l/b, etc.)
        dist_mpc          — adopted distance (Mpc)
        e_dist_mpc        — distance uncertainty (Mpc)
        incl_deg          — inclination (degrees); NaN if not provided
        e_incl_deg        — inclination uncertainty (degrees); NaN if not provided
        posang_deg        — position angle E of N (degrees); NaN if not provided
        e_posang_deg      — PA uncertainty (degrees); NaN if not provided
        r25               — optical radius at 25 mag/arcsec² (arcmin); NaN if not provided
        e_r25             — r25 uncertainty (arcmin); NaN if not provided

        The galaxy-geometry columns (incl_deg, posang_deg, r25 and their
        uncertainties) are optional — they are NaN when absent from
        target_definitions.txt.  Use has_galaxy_geometry() to check whether
        deprojected radii and polar angles can be computed.

        Raises
        ------
        KeyError if *source* is not in the geometry table.
        """
        if source not in self._index:
            LOG.error(f"Source '{source}' not found in geometry table.")
            raise KeyError(f"Source '{source}' not found in geometry table.")
        row = self.source_table.loc[self._index[source]].to_dict()
        # Fill any missing optional keys with NaN so callers always get a
        # complete dict regardless of how many columns the file had
        for key in ["incl_deg", "e_incl_deg", "posang_deg", "e_posang_deg",
                    "r25", "e_r25", "dist_mpc", "e_dist_mpc"]:
            if key not in row or row[key] is None:
                row.setdefault(key, float("nan"))
        return row

    def has_galaxy_geometry(self, source: str) -> bool:
        """
        Return True if *source* has all three galaxy-geometry values needed for
        deprojection: ``incl_deg``, ``posang_deg``, and ``r25``.

        When any of these is NaN, rgal/theta columns cannot be computed and the
        corresponding pipeline steps are skipped with a warning.
        """
        import math
        p = self.get_source_params(source)
        return not any(
            math.isnan(float(p.get(k, float("nan"))))
            for k in ("incl_deg", "posang_deg", "r25")
        )

    # Convenience accessors for the most commonly needed parameters

    def get_ra_ctr(self, source: str) -> float:
        """RA of source centre (degrees, J2000)."""
        return self.get_source_params(source)["ra_ctr"]

    def get_dec_ctr(self, source: str) -> float:
        """Dec of source centre (degrees, J2000)."""
        return self.get_source_params(source)["dec_ctr"]

    def get_dist_mpc(self, source: str) -> float:
        """Adopted distance in Mpc."""
        return self.get_source_params(source)["dist_mpc"]

    def get_incl_deg(self, source: str) -> float:
        """Inclination in degrees."""
        return self.get_source_params(source)["incl_deg"]

    def get_posang_deg(self, source: str) -> float:
        """Position angle (E of N) in degrees."""
        return self.get_source_params(source)["posang_deg"]

    def get_r25(self, source: str) -> float:
        """Optical radius r25 in arcmin."""
        return self.get_source_params(source)["r25"]

    def n_sources(self) -> int:
        """Number of sources to process."""
        return len(self.sources)

    def all_sources(self) -> list:
        """Return a copy of the ordered source list."""
        return list(self.sources)

    def __repr__(self):
        return (
            f"SourceHandler(n_sources={self.n_sources()}, " f"sources={self.sources})"
        )
