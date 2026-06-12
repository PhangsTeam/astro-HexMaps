"""
SourceHandler: manages the list of targets and their geometric parameters.

Wraps the source table loaded by KeyHandler and provides convenient
per-source lookups used throughout the pipeline stages.
"""

import numpy as np
import pandas as pd
from pathlib import Path


class SourceHandler:
    """
    Manages the source list and their geometry parameters.

    Parameters
    ----------
    source_table : pd.DataFrame
        DataFrame with columns as defined in TARGET_COLUMNS (from KeyHandler).
    sources : list of str
        Ordered list of source names to process.  Must be a subset of the
        source_table 'source' column.

    Example
    -------
    >>> from pystructurePipeline.handlerKeys import KeyHandler, SourceHandler
    >>> kh = KeyHandler("./keys/")
    >>> th = SourceHandler(kh.get_source_table(), kh.get_sources())
    >>> params = th.get_source_params("ngc5194")
    """

    def __init__(self, source_table: pd.DataFrame, sources: list):
        self.source_table = source_table.copy()
        self.sources = list(sources)
        self._index = {
            row["source"]: idx
            for idx, row in source_table.iterrows()
        }
        self._validate()

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate(self):
        missing = [s for s in self.sources if s not in self._index]
        if missing:
            raise ValueError(
                f"The following sources are not in the geometry table: {missing}"
            )

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def get_source_params(self, source: str) -> dict:
        """
        Return a dict of geometric parameters for *source*.

        Keys
        ----
        ra_ctr, dec_ctr, dist_mpc, e_dist_mpc,
        incl_deg, e_incl_deg, posang_deg, e_posang_deg, r25, e_r25
        """
        if source not in self._index:
            raise KeyError(f"Source '{source}' not found in geometry table.")
        row = self.source_table.loc[self._index[source]]
        return row.to_dict()

    def get_ra_ctr(self, source: str) -> float:
        return self.get_source_params(source)["ra_ctr"]

    def get_dec_ctr(self, source: str) -> float:
        return self.get_source_params(source)["dec_ctr"]

    def get_dist_mpc(self, source: str) -> float:
        return self.get_source_params(source)["dist_mpc"]

    def get_incl_deg(self, source: str) -> float:
        return self.get_source_params(source)["incl_deg"]

    def get_posang_deg(self, source: str) -> float:
        return self.get_source_params(source)["posang_deg"]

    def get_r25(self, source: str) -> float:
        return self.get_source_params(source)["r25"]

    def n_sources(self) -> int:
        return len(self.sources)

    def all_sources(self) -> list:
        return list(self.sources)

    def __repr__(self):
        return (
            f"SourceHandler(n_sources={self.n_sources()}, "
            f"sources={self.sources})"
        )
