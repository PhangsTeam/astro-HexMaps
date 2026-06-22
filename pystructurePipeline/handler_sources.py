"""
handler_sources.py — SourceHandler: manages the source list and geometry.

Wraps the source geometry table loaded by KeyHandler and provides
named-parameter lookups used by every pipeline stage that needs
per-source geometry (position angle, inclination, distance, r25).
"""

import pandas as pd

from pystructurePipeline.pystructureLogger import get_logger

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
        ra_ctr, dec_ctr   — J2000 centre coordinates (degrees)
        dist_mpc          — adopted distance (Mpc)
        e_dist_mpc        — distance uncertainty (Mpc)
        incl_deg          — inclination (degrees)
        e_incl_deg        — inclination uncertainty (degrees)
        posang_deg        — position angle E of N (degrees)
        e_posang_deg      — PA uncertainty (degrees)
        r25               — optical radius at 25 mag/arcsec² (arcmin)
        e_r25             — r25 uncertainty (arcmin)

        Raises
        ------
        KeyError if *source* is not in the geometry table.
        """
        if source not in self._index:
            LOG.error(f"Source '{source}' not found in geometry table.")
            raise KeyError(f"Source '{source}' not found in geometry table.")
        return self.source_table.loc[self._index[source]].to_dict()

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
