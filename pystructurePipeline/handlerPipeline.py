"""
handlerPipeline.py — PipelineHandler: orchestrates the PyStructure pipeline.

This is the main entry point for programmatic use of PyStructurePipeline.
It loads all key files, validates the configuration, creates the output
directory, and then dispatches the requested pipeline stages for each source.

Pipeline stages (in execution order)
--------------------------------------
sampling
    Generate the hexagonal sampling grid from the overlay cube.
    Output: stored in memory and passed directly to regrid.

regrid
    Convolve each input map and cube to the target resolution, reproject
    onto the overlay WCS, and sample at the hex-grid points.
    Output: .ecsv file per source written to out_dir.

spectra
    Read the .ecsv file, construct the S/N mask from the reference line,
    compute moment maps (mom0/1/2, Tpeak, rms, EW) and shuffled spectra
    for every line, then overwrite the .ecsv.

output
    Read the final .ecsv and regrid the moment maps and 2D maps back
    onto a rectangular pixel grid, writing one FITS file per quantity.

Usage
-----
Programmatic (from Python)::

    from pystructurePipeline import PipelineHandler
    handler = PipelineHandler(key_dir="keys/")
    handler.run_all()                              # all stages, all sources
    handler.run_stages(["regrid", "spectra"])      # subset of stages
    handler.run_stages(["regrid"], targets=["ngc5194"])  # subset of sources

Command-line (after pip install)::

    pystructure --key_dir keys/
    pystructure --key_dir keys/ --stages sampling regrid
    pystructure --key_dir keys/ --targets ngc5194 ngc5457
"""

import os
from pathlib import Path
from datetime import date

from pystructurePipeline.handlerKeys    import KeyHandler
from pystructurePipeline.handlerSources import SourceHandler
from pystructurePipeline.pystructureLogger import get_logger, logger

ALL_STAGES = ["sampling", "regrid", "spectra", "output"]
_log = get_logger("Pipeline")


class PipelineHandler:
    """
    Orchestrates the PyStructure pipeline stages.

    Parameters
    ----------
    key_dir : str or Path
        Directory containing master_key.txt (and typically the other key
        files, unless they are located elsewhere as configured in master_key).
    verbose : bool, optional
        If True (default), print progress messages to stdout.
    log_file : str, optional
        If given, write every log message to this file as it is produced
        (in addition to printing, if verbose=True).  The file is created
        (overwriting any existing content) when the handler is constructed.

    Attributes
    ----------
    key_handler    : KeyHandler    — loaded configuration
    source_handler : SourceHandler — source geometry lookups
    run_success    : dict          — maps source name → bool after a run
    """

    def __init__(self, key_dir: str, verbose: bool = True, log_file: str = None):
        self.key_dir = Path(key_dir)
        self.verbose = verbose
        self.log_file = log_file
        self.run_success = {}

        # Configure the shared logger: controls whether messages are printed
        # to stdout and, if log_file is given, streams every message to that
        # file as it is logged (in addition to printing).
        logger.configure(verbose=verbose, log_file=log_file)

        _log.info("Loading key files...")
        self.key_handler = KeyHandler(key_dir)
        self.key_handler.validate()

        self.source_handler = SourceHandler(
            self.key_handler.get_source_table(),
            self.key_handler.get_sources(),
        )
        _log.info(f"Loaded {self.source_handler.n_sources()} source(s): "
                  f"{self.source_handler.all_sources()}")

        # Ensure the output directory exists before any stage tries to write
        out_dir = self.key_handler.meta.get("out_dir", "output/")
        os.makedirs(out_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Public run interface
    # ------------------------------------------------------------------

    def run_all(self, targets: list = None):
        """
        Run all four pipeline stages in order for the given sources.

        Parameters
        ----------
        targets : list of str, optional
            Restrict to these source names.  Defaults to all sources in
            the key files.
        """
        self.run_stages(ALL_STAGES, targets=targets)

    def run_stages(self, stages: list, targets: list = None):
        """
        Run a specified subset of pipeline stages.

        Stages are always executed in the canonical order (sampling → regrid
        → spectra → output) regardless of the order they appear in *stages*.
        This means you can safely pass ["spectra", "regrid"] and the regrid
        stage will still run before spectra.

        Parameters
        ----------
        stages : list of str
            Stage names to execute.  Must be a subset of:
            "sampling", "regrid", "spectra", "output".
        targets : list of str, optional
            Restrict processing to these source names.  Defaults to all
            sources defined in the key files.

        Raises
        ------
        ValueError if any element of *stages* is not a valid stage name.
        """
        unknown = [s for s in stages if s not in ALL_STAGES]
        if unknown:
            _log.error(f"Unknown stage(s): {unknown}. Valid stages: {ALL_STAGES}")
            raise ValueError(f"Unknown stage(s): {unknown}. Valid stages: {ALL_STAGES}")

        # Preserve canonical stage order
        ordered = [s for s in ALL_STAGES if s in stages]

        source_list = targets if targets else self.source_handler.all_sources()
        self.run_success = {s: True for s in source_list}

        _log.info(f"Running stages {ordered} for source(s): {source_list}")

        for source in source_list:
            _log.info(f"--- Processing source: {source} ---")
            try:
                if "sampling" in ordered:
                    self._run_sampling(source)
                if "regrid" in ordered:
                    self._run_regrid(source)
                if "spectra" in ordered:
                    self._run_spectra(source)
                if "output" in ordered:
                    self._run_output(source)
            except Exception as exc:
                self.run_success[source] = False
                _log.error(f"Stage failed for {source}: {exc}")
                import traceback; traceback.print_exc()

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

    def _run_sampling(self, source: str):
        """
        Dispatch the sampling stage for *source*.

        Generates the hexagonal sampling grid from the overlay cube and stores
        the result for use by the regrid stage.  The grid spacing is
        target_res / spacing_per_beam degrees.
        """
        from pystructurePipeline.pipelineSampling import run_sampling
        _log.info(f"Generating hexagonal grid for {source}.")
        run_sampling(
            source = source,
            params = self.source_handler.get_source_params(source),
            meta   = self.key_handler.meta,
        )

    def _run_regrid(self, source: str):
        """
        Dispatch the regrid stage for *source*.

        Convolves every input map and cube to the target resolution, reprojects
        onto the overlay WCS, samples at the hex-grid points, computes
        deprojected coordinates, and writes the result to a .ecsv file.
        """
        from pystructurePipeline.pipelineRegrid import run_regrid
        _log.info(f"Convolving and sampling data for {source}.")
        run_regrid(
            source     = source,
            params     = self.source_handler.get_source_params(source),
            meta       = self.key_handler.meta,
            maps       = self.key_handler.get_maps(),
            cubes      = self.key_handler.get_cubes(),
            input_mask = self.key_handler.get_input_mask(),
        )

    def _run_spectra(self, source: str):
        """
        Dispatch the spectra stage for *source*.

        Reads the .ecsv written by regrid, builds the S/N mask, computes
        moment maps and shuffled spectra for every line, and overwrites the
        .ecsv with the enriched table.
        """
        from pystructurePipeline.pipelineProducts import run_spectra
        _log.info(f"Processing spectra for {source}.")
        run_spectra(
            source     = source,
            fname      = self._get_output_fname(source),
            meta       = self.key_handler.meta,
            cubes      = self.key_handler.get_cubes(),
            input_mask = self.key_handler.get_input_mask(),
            hfs_data   = self.key_handler.get_hfs_data(),
        )

    def _run_output(self, source: str):
        """
        Dispatch the output stage for *source*.

        Reads the final .ecsv and regrid the moment maps and 2D maps onto a
        rectangular pixel grid, writing one FITS file per quantity into
        folder_savefits.
        """
        from pystructurePipeline.pipelineFITS import run_output
        _log.info(f"Writing FITS maps for {source}.")
        run_output(
            source = source,
            fname  = self._get_output_fname(source),
            meta   = self.key_handler.meta,
            maps   = self.key_handler.get_maps(),
            cubes  = self.key_handler.get_cubes(),
            params = self.source_handler.get_source_params(source),
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _get_output_fname(self, source: str) -> str:
        """
        Build the .ecsv output filename for *source*.

        The filename encodes the source name, resolution and units, and the
        current date so that successive runs do not silently overwrite each
        other (unless structure_creation = "default").

        Examples
        --------
        angular mode, 27 arcsec → ngc5194_data_struct_27as_2025_06_01.ecsv
        physical mode, 100 pc   → ngc5194_data_struct_100pc_2025_06_01.ecsv
        """
        meta       = self.key_handler.meta
        out_dir    = meta.get("out_dir", "output/")
        resolution = meta.get("resolution", "angular")
        target_res = meta.get("target_res", 27.0)

        suffix = (str(int(target_res)) + "as" if resolution == "angular"
                  else str(int(target_res)) + "pc" if resolution == "physical"
                  else "native")

        date_str = date.today().strftime("%Y_%m_%d")
        fname    = os.path.join(out_dir, f"{source}_data_struct_{suffix}_{date_str}.ecsv")

        # In archive mode, bump the version number if the file already exists
        if "archive" in meta.get("structure_creation", "") and os.path.exists(fname):
            version = 1
            base = fname[:-5]
            while os.path.exists(f"{base}_v{version}.ecsv"):
                version += 1
            fname = f"{base}_v{version}.ecsv"

        return fname

    def _print_summary(self):
        """Print a per-source pass/fail summary after all stages complete."""
        _log.info("--- Run summary ---")
        all_ok = True
        for source, ok in self.run_success.items():
            status = "OK" if ok else "FAILED"
            _log.info(f"  {source}: {status}")
            if not ok:
                all_ok = False
        if all_ok:
            _log.info("All sources completed successfully.")
        else:
            _log.warning("Some sources failed — check errors above.")

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
            f"PipelineHandler(key_dir='{self.key_dir}', "
            f"sources={self.source_handler.all_sources()})"
        )
