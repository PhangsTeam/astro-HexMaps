"""
handler_pipeline.py — PipelineHandler: orchestrates the HexMaps pipeline.

This is the main entry point for programmatic use of HexMapsPipeline.
It loads all key files, validates the configuration, creates the output
directory, and then dispatches the requested pipeline stages for each target.

Pipeline stages (in execution order)
--------------------------------------
regrid
    Generate the hexagonal sampling grid from the overlay cube, then convolve
    each input map and cube to the target resolution, reproject onto the
    overlay WCS, and sample at the hex-grid points.
    Output: .ecsv file per target written to out_dir.

products
    Read the .ecsv file, construct the S/N mask from the reference line,
    compute moment maps (mom0/1/2, Tpeak, rms, EW) and shuffled spectra
    for every line, then overwrite the .ecsv.

fits
    Read the final .ecsv and regrid the moment maps and 2D maps back
    onto a rectangular pixel grid, writing one FITS file per quantity.

Usage
-----
Programmatic (from Python)::

    from hexmaps import PipelineHandler
    handler = PipelineHandler(conf_path="config.txt")
    handler.run_all()                              # regrid + products (default), all targets
    handler.run_stages(["regrid", "products"])     # subset of stages
    handler.run_stages(["regrid"], targets=["ngc5194"])  # subset of targets

Command-line (after pip install)::

    hexmaps --conf config.txt
    hexmaps --conf config.txt --stages regrid products fits
    hexmaps --conf config.txt --targets ngc5194 ngc5457
"""

import os
import numpy as np
from pathlib import Path
from datetime import date

from hexmaps.handler_keys import KeyHandler
from hexmaps.handler_targets import TargetHandler
from hexmaps.logger import get_logger, logger

ALL_STAGES = ["regrid", "products", "fits"]
DATABASE_STAGES = ["regrid", "products"]  # default: omits the optional fits stage

# Two loggers bookend a run: "Loading" covers configuration/setup before the
# regrid stage begins, and "Return" covers the run summary. Per-stage progress
# messages are logged by the stage modules themselves (stage_regrid -> "Regrid",
# stage_products -> "Products", stage_fits -> "FITS").
LOG_LOADING = get_logger("Loading")
LOG_RETURN = get_logger("Return")


class PipelineHandler:
    """
    Orchestrates the HexMaps pipeline stages.

    Parameters
    ----------
    conf_path : str or Path
        Path to config.txt — the single configuration file containing paths,
        metadata, the target/overlay/maps/cubes/mask tables, and all pipeline
        settings. The target geometry table (keys/target_definitions.txt) and
        optional hyperfine-structure file (keys/hfs_lines.txt) are looked up
        in a `keys/` subfolder next to config.txt.
    verbose : bool, optional
        If True (default), print progress messages to stdout.
    log_file : str, optional
        If given, write every log message to this file as it is produced
        (in addition to printing, if verbose=True).  The file is created
        (overwriting any existing content) when the handler is constructed.

    Attributes
    ----------
    key_handler    : KeyHandler    — loaded configuration
    target_handler : TargetHandler — target geometry lookups
    run_success    : dict          — maps target name → bool after a run
    """

    def __init__(self, conf_path: str, verbose: bool = True, log_file: str = None):
        self.conf_path = Path(conf_path)
        self.verbose = verbose
        self.log_file = log_file
        self.run_success = {}

        # Configure the shared logger: controls whether messages are printed
        # to stdout and, if log_file is given, streams every message to that
        # file as it is logged (in addition to printing).
        logger.configure(verbose=verbose, log_file=log_file)

        LOG_LOADING.info("Loading config file ...")
        self.key_handler = KeyHandler(conf_path)
        self.key_handler.validate()

        self.target_handler = TargetHandler(
            self.key_handler.get_target_table(),
            self.key_handler.get_targets(),
        )
        LOG_LOADING.info(
            f"Loaded {self.target_handler.n_targets()} target(s): "
            f"{self.target_handler.all_targets()}"
        )

        # Ensure the output directory exists before any stage tries to write
        out_dir = self.key_handler.meta.get("out_dir", "output/")
        os.makedirs(out_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Public run interface
    # ------------------------------------------------------------------

    def run_all(self, targets: list = None):
        """
        Run the default pipeline stages (regrid + products) for the given targets.

        The ``fits`` stage is intentionally excluded from the default run
        because it is optional — it produces FITS moment maps and band images
        as a convenience output, but the primary pipeline deliverable is the
        .ecsv database written by the products stage.  To include FITS output,
        call ``run_stages(["regrid", "products", "fits"])`` explicitly, or pass
        ``--stages regrid products fits`` on the command line.

        Parameters
        ----------
        targets : list of str, optional
            Restrict to these target names.  Defaults to all targets in
            config.txt.
        """
        self.run_stages(DATABASE_STAGES, targets=targets)

    def run_stages(self, stages: list, targets: list = None):
        """
        Run a specified subset of pipeline stages.

        Stages are always executed in the canonical order (regrid
        → products → fits) regardless of the order they appear in *stages*.
        This means you can safely pass ["products", "regrid"] and the regrid
        stage will still run before products.

        Passing ``"all"`` (or ``["all"]``) as *stages* is a shorthand for
        all stages in ``ALL_STAGES``, equivalent to
        ``["regrid", "products", "fits"]``.

        Parameters
        ----------
        stages : list of str or str
            Stage names to execute.  Must be a subset of:
            "regrid", "products", "fits".
            The special value ``"all"`` expands to all available stages.
        targets : list of str, optional
            Restrict processing to these target names.  Defaults to all
            targets defined in config.txt.

        Raises
        ------
        ValueError if any element of *stages* is not a valid stage name.
        """
        # Accept a bare string as well as a list
        if isinstance(stages, str):
            stages = [stages]

        # Expand "all" to every available stage
        if "all" in stages:
            stages = ALL_STAGES

        unknown = [s for s in stages if s not in ALL_STAGES]
        if unknown:
            LOG_LOADING.error(
                f"Unknown stage(s): {unknown}. Valid stages: {ALL_STAGES} (or 'all')"
            )
            raise ValueError(
                f"Unknown stage(s): {unknown}. Valid stages: {ALL_STAGES} (or 'all')"
            )

        # Preserve canonical stage order
        ordered = [s for s in ALL_STAGES if s in stages]

        target_list = targets if targets else self.target_handler.all_targets()
        self.run_success = {s: True for s in target_list}

        LOG_LOADING.info(f"Running stages: {ordered}")

        if ordered == DATABASE_STAGES and "fits" not in ordered:
            LOG_LOADING.info(
                "The 'fits' stage is not included in the current run. "
                "To enable FITS output, add 'fits' to your stage list."
            )

        for target in target_list:
            LOG_LOADING.info(f"--- Processing target: {target} ---")
            try:
                if "regrid" in ordered:
                    self._run_regrid(target)
                if "products" in ordered:
                    self._run_products(target)
                if "fits" in ordered:
                    self._run_fits(target)
            except Exception as exc:
                self.run_success[target] = False
                LOG_RETURN.error(f"Stage failed for {target}: {exc}")
                import traceback

                traceback.print_exc()

        self._print_summary()

        # If a log file was configured, ensure the full record is flushed
        # (the per-message streaming already wrote each line, but save()
        # rewrites the complete, ordered log in one go).
        if self.log_file:
            logger.save(self.log_file)

    # ------------------------------------------------------------------
    # Stage dispatch
    #
    # Each method imports its stage module lazily so that importing
    # PipelineHandler does not pull in all stage dependencies (astropy,
    # reproject, scipy, …) until they are actually needed.
    # ------------------------------------------------------------------

    def _run_regrid(self, target: str):
        """
        Dispatch the regrid stage for *target*.

        Generates the hexagonal sampling grid from the overlay cube, then
        convolves every input map and cube to the target resolution, reprojects
        onto the overlay WCS, samples at the hex-grid points, computes
        deprojected coordinates, and writes the result to a .ecsv file.
        """
        from hexmaps.stage_regrid import run_regrid, LOG as REGRID_LOG

        REGRID_LOG.info(f"Convolving and sampling data for {target}.")
        run_regrid(
            target=target,
            params=self.target_handler.get_target_params(target),
            meta=self.key_handler.meta,
            maps=self.key_handler.get_maps(),
            cubes=self.key_handler.get_cubes(),
            input_mask=self.key_handler.get_input_mask(),
            window_mask=self.key_handler.get_window_mask(),
        )

    def _run_products(self, target: str):
        """
        Dispatch the products stage for *target*.

        Reads the .ecsv written by regrid (discovers the most recent existing
        file via _find_output_fname so products can run independently without
        regrid in the same session), builds the S/N mask, computes moment
        maps and shuffled spectra for every line, and overwrites the .ecsv
        with the enriched table.
        """
        from hexmaps.stage_products import run_products, LOG as PRODUCTS_LOG

        PRODUCTS_LOG.info(f"Create products for target: {target} ...")
        run_products(
            target=target,
            fname=self._find_output_fname(target),
            meta=self.key_handler.meta,
            cubes=self.key_handler.get_cubes(),
            input_mask=self.key_handler.get_input_mask(),
            window_mask=self.key_handler.get_window_mask(),
            hfs_data=self.key_handler.get_hfs_data(),
            noise_mask_df=self.key_handler.get_noise_mask(),
        )

    def _run_fits(self, target: str):
        """
        Dispatch the fits stage for *target*.

        Computes moment maps directly on the convolved PPV cubes (not the
        hex-grid .ecsv table), regrid the 2D map columns onto a rectangular
        pixel grid, and optionally writes the velocity-integration mask(s)
        as FITS cubes — all into folder_savefits. Uses _find_output_fname
        so the fits stage can run independently without regrid or products
        in the same session.
        """
        from hexmaps.stage_fits import run_fits, LOG as FITS_LOG

        FITS_LOG.info(f"Creating FITS files for target: {target} ...")
        run_fits(
            target=target,
            fname=self._find_output_fname(target),
            meta=self.key_handler.meta,
            maps=self.key_handler.get_maps(),
            cubes=self.key_handler.get_cubes(),
            params=self.target_handler.get_target_params(target),
            input_mask=self.key_handler.get_input_mask(),
            window_mask=self.key_handler.get_window_mask(),
            hfs_data=self.key_handler.get_hfs_data(),
            noise_mask_df=self.key_handler.get_noise_mask(),
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _get_output_fname(self, target: str) -> str:
        """
        Build the .ecsv output filename for *target*.

        Reads ``meta["res_suffix"]`` — set by handler_keys._resolve_resolution
        at config-load time and overwritten per-target by run_sampling —
        so the suffix is always consistent with the resolution mode.

        Examples
        --------
        angular 27 arcsec  → ngc5194_hexmaps_27p0as_2025_06_01.ecsv
        physical 100 pc    → ngc5194_hexmaps_100pc_2025_06_01.ecsv
        native 12.8 arcsec → ngc5194_hexmaps_12p8as_2025_06_01.ecsv
        """
        meta = self.key_handler.meta
        out_dir = meta.get("out_dir", "output/")
        res_suffix = meta.get("res_suffix", "27p0as")
        date_str = date.today().strftime("%Y_%m_%d")
        fname = os.path.join(out_dir, f"{target}_hexmaps_{res_suffix}_{date_str}.ecsv")

        # In archive mode, bump the version number if the file already exists
        if "archive" in meta.get("structure_creation", "") and os.path.exists(fname):
            version = 1
            base = fname[:-5]
            while os.path.exists(f"{base}_v{version}.ecsv"):
                version += 1
            fname = f"{base}_v{version}.ecsv"

        return fname

    def _find_output_fname(self, target: str) -> str:
        """
        Find the most recent existing .ecsv for *target*, or fall back to
        ``_get_output_fname`` (today's date) if none exists yet.

        Used by products and fits stages when running independently (without
        regrid in the same session). Globs for
        ``{target}_hexmaps_{res_suffix}_*.ecsv`` in out_dir and returns
        the most recently modified match, so re-running products or fits
        after a prior regrid just works without needing to supply the date.
        """
        import glob

        meta = self.key_handler.meta
        out_dir = meta.get("out_dir", "output/")
        res_suffix = meta.get("res_suffix", "27p0as")
        pattern = os.path.join(out_dir, f"{target}_hexmaps_{res_suffix}_*.ecsv")
        matches = sorted(glob.glob(pattern), key=os.path.getmtime, reverse=True)
        if matches:
            LOG_LOADING.info(f"Found existing database for {target}: {matches[0]}")
            return matches[0]
        return self._get_output_fname(target)

    def _print_summary(self):
        """Print a per-target pass/fail summary after all stages complete."""
        LOG_RETURN.info("--- Run summary ---")
        all_ok = True
        for target, ok in self.run_success.items():
            status = "OK" if ok else "FAILED"
            LOG_RETURN.info(f"  {target}: {status}")
            if not ok:
                all_ok = False
        if all_ok:
            LOG_RETURN.info("All targets completed successfully.")
        else:
            LOG_RETURN.warning("Some targets failed — check errors above.")

    def save_log(self, path: str):
        """
        Write the full log history (all messages from this run) to *path*.

        This can be called at any time, not just at the end of a run, and is
        useful if PipelineHandler was created without log_file but you decide
        afterwards that you want to keep a record.

        Parameters
        ----------
        path : str — output file path
        """
        logger.save(path)

    def __repr__(self):
        return (
            f"PipelineHandler(conf_path='{self.conf_path}', "
            f"targets={self.target_handler.all_targets()})"
        )
