"""
handler_sources.py — TargetHandler: manages the target list and geometry.

Wraps the target geometry table loaded by KeyHandler and provides
named-parameter lookups used by every pipeline stage that needs
per-target geometry (position angle, inclination, distance, r25).
"""

import pandas as pd

from hexmaps.logger import get_logger

LOG = get_logger("Loading")


class TargetHandler:
    """
    Manages the target list and their geometric parameters.

    Parameters
    ----------
    target_table : pd.DataFrame
        Full geometry table from target_definitions.txt, as loaded by
        KeyHandler.  Must contain a 'target' column plus the geometry columns
        defined in handler_keys.TARGET_COLUMNS.
    targets : list of str
        Ordered list of target names to process.  Must be a subset of the
        names in target_table['target'].

    Raises
    ------
    ValueError
        If any name in *targets* is not present in *target_table*.

    Example
    -------
    >>> kh = KeyHandler("./keys/")
    >>> sh = TargetHandler(kh.get_source_table(), kh.get_sources())
    >>> params = sh.get_target_params("ngc5194")
    >>> print(params["dist_mpc"])
    """

    def __init__(self, target_table: pd.DataFrame, targets: list):
        self.target_table = target_table.copy()
        self.targets = list(targets)
        # Build a name → row-index lookup for O(1) access
        self._index = {row["target"]: idx for idx, row in target_table.iterrows()}
        self._validate()

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate(self):
        """
        Check that every requested target is present in the geometry table.

        Raises ValueError listing all missing names so the user can fix
        target_definitions.txt or config.txt in one go.
        """
        missing = [s for s in self.targets if s not in self._index]
        if missing:
            LOG.error(
                f"The following targets are not in "
                f"target_definitions.txt: {missing}"
            )
            raise ValueError(
                f"The following targets are not in "
                f"target_definitions.txt: {missing}"
            )

    # ------------------------------------------------------------------
    # Parameter lookups
    # ------------------------------------------------------------------

    def get_target_params(self, target: str) -> dict:
        """
        Return a dict of all geometric parameters for *target*.

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
        KeyError if *target* is not in the geometry table.
        """
        if target not in self._index:
            LOG.error(f"Target '{target}' not found in geometry table.")
            raise KeyError(f"Target '{target}' not found in geometry table.")
        row = self.target_table.loc[self._index[target]].to_dict()
        # Fill any missing optional keys with NaN so callers always get a
        # complete dict regardless of how many columns the file had
        for key in ["incl_deg", "e_incl_deg", "posang_deg", "e_posang_deg",
                    "r25", "e_r25", "dist_mpc", "e_dist_mpc"]:
            if key not in row or row[key] is None:
                row.setdefault(key, float("nan"))
        return row

    def has_galaxy_geometry(self, target: str) -> bool:
        """
        Return True if *target* has all three galaxy-geometry values needed for
        deprojection: ``incl_deg``, ``posang_deg``, and ``r25``.

        When any of these is NaN, rgal/theta columns cannot be computed and the
        corresponding pipeline steps are skipped with a warning.
        """
        import math
        p = self.get_target_params(target)
        return not any(
            math.isnan(float(p.get(k, float("nan"))))
            for k in ("incl_deg", "posang_deg", "r25")
        )

    # Convenience accessors for the most commonly needed parameters

    def get_ra_ctr(self, target: str) -> float:
        """RA of target centre (degrees, J2000)."""
        return self.get_target_params(target)["ra_ctr"]

    def get_dec_ctr(self, target: str) -> float:
        """Dec of target centre (degrees, J2000)."""
        return self.get_target_params(target)["dec_ctr"]

    def get_dist_mpc(self, target: str) -> float:
        """Adopted distance in Mpc."""
        return self.get_target_params(target)["dist_mpc"]

    def get_incl_deg(self, target: str) -> float:
        """Inclination in degrees."""
        return self.get_target_params(target)["incl_deg"]

    def get_posang_deg(self, target: str) -> float:
        """Position angle (E of N) in degrees."""
        return self.get_target_params(target)["posang_deg"]

    def get_r25(self, target: str) -> float:
        """Optical radius r25 in arcmin."""
        return self.get_target_params(target)["r25"]

    def n_sources(self) -> int:
        """Number of targets to process."""
        return len(self.targets)

    def all_sources(self) -> list:
        """Return a copy of the ordered target list."""
        return list(self.targets)

    def __repr__(self):
        return (
            f"TargetHandler(n_sources={self.n_sources()}, " f"targets={self.targets})"
        )
