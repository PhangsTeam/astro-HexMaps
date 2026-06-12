"""
PipelineHandler: orchestrates the PyStructure pipeline stages.

This is the main entry point for programmatic use of PyStructure.
Users instantiate a PipelineHandler with a key directory, then call
run_stages() to execute any combination of pipeline stages.

Available stages (in order)
----------------------------
  sampling      - Generate hexagonal sampling grid for each source
  regrid        - Convolve and sample bands and cubes onto the grid
  spectra       - Process spectra: masking, moments, shuffling
  output        - Write FITS moment / 2D maps

Example
-------
>>> from pystructurePipeline import PipelineHandler
>>> handler = PipelineHandler(key_dir="./keys/")
>>> handler.run_stages(["sampling", "regrid", "spectra", "output"])

Or run everything:
>>> handler.run_all()

Or a single stage for one source:
>>> handler.run_stages(["regrid"], targets=["ngc5194"])
"""

import os
import shutil
from pathlib import Path
from datetime import date

from pystructurePipeline.handlerKeys import KeyHandler
from pystructurePipeline.handlerSources import SourceHandler


# Ordered list of valid stage names
ALL_STAGES = ["sampling", "regrid", "spectra", "output"]


class PipelineHandler:
    """
    Orchestrates the PyStructure pipeline.

    Parameters
    ----------
    key_dir : str or Path
        Directory containing all key files.
    verbose : bool
        If True, print progress messages.  Default True.

    Attributes
    ----------
    key_handler : KeyHandler
    source_handler : SourceHandler
    run_success : dict
        Maps source name → bool after a run.
    """

    def __init__(self, key_dir: str, verbose: bool = True):
        self.key_dir = Path(key_dir)
        self.verbose = verbose
        self.run_success = {}

        self._log("[INFO]", "Loading key files...")
        self.key_handler = KeyHandler(key_dir)
        self.key_handler.validate()

        self.source_handler = SourceHandler(
            self.key_handler.get_source_table(),
            self.key_handler.get_sources(),
        )
        self._log("[INFO]", f"Loaded {self.source_handler.n_sources()} source(s): "
                            f"{self.source_handler.all_sources()}")

        # Ensure output directory exists
        out_dir = self.key_handler.meta.get("out_dir", "Output/")
        os.makedirs(out_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Public run interface
    # ------------------------------------------------------------------

    def run_all(self, targets: list = None):
        """Run all pipeline stages in order for the given targets."""
        self.run_stages(ALL_STAGES, targets=targets)

    def run_stages(self, stages: list, targets: list = None):
        """
        Run the specified subset of pipeline stages.

        Parameters
        ----------
        stages : list of str
            Stage names to execute, e.g. ["sampling", "regrid"].
            Must be a subset of: sampling, regrid, spectra, output.
        targets : list of str, optional
            Restrict processing to these source names.
            Defaults to all sources defined in the key files.
        """
        unknown = [s for s in stages if s not in ALL_STAGES]
        if unknown:
            raise ValueError(f"Unknown stage(s): {unknown}. Valid: {ALL_STAGES}")

        source_list = targets if targets else self.source_handler.all_sources()
        self.run_success = {s: True for s in source_list}

        self._log("[INFO]", f"Running stages {stages} for source(s): {source_list}")

        for source in source_list:
            self._log("[INFO]", f"--- Processing source: {source} ---")
            try:
                if "sampling" in stages:
                    self._run_sampling(source)
                if "regrid" in stages:
                    self._run_regrid(source)
                if "spectra" in stages:
                    self._run_spectra(source)
                if "output" in stages:
                    self._run_output(source)
            except Exception as exc:
                self.run_success[source] = False
                self._log("[ERROR]", f"Stage failed for {source}: {exc}")

        self._print_summary()

    # ------------------------------------------------------------------
    # Stage dispatch methods
    # These are thin wrappers that import the stage modules lazily so that
    # individual stages can be used without importing the whole pipeline.
    # ------------------------------------------------------------------

    def _run_sampling(self, source: str):
        """Generate hexagonal sampling grid for *source*."""
        from pystructurePipeline.pipelineSampling import run_sampling
        self._log("[INFO]", f"[sampling] Generating hexagonal grid for {source}.")
        meta = self.key_handler.meta
        params = self.source_handler.get_source_params(source)
        run_sampling(source=source, params=params, meta=meta)

    def _run_regrid(self, source: str):
        """Convolve and sample bands/cubes for *source*."""
        from pystructurePipeline.pipelineRegrid import run_regrid
        self._log("[INFO]", f"[regrid] Convolving and sampling data for {source}.")
        meta  = self.key_handler.meta
        maps = self.key_handler.get_maps()
        cubes = self.key_handler.get_cubes()
        params = self.source_handler.get_source_params(source)
        run_regrid(
            source=source,
            params=params,
            meta=meta,
            maps=maps,
            cubes=cubes,
            input_mask=self.key_handler.get_input_mask(),
        )

    def _run_spectra(self, source: str):
        """Process spectra (masking, moments) for *source*."""
        from pystructurePipeline.pipelineSpectra import run_spectra
        self._log("[INFO]", f"[spectra] Processing spectra for {source}.")
        meta  = self.key_handler.meta
        cubes = self.key_handler.get_cubes()
        fname = self._get_output_fname(source)
        run_spectra(
            source=source,
            fname=fname,
            meta=meta,
            cubes=cubes,
            input_mask=self.key_handler.get_input_mask(),
            hfs_data=self.key_handler.get_hfs_data(),
        )

    def _run_output(self, source: str):
        """Write FITS moment/2D maps for *source*."""
        from pystructurePipeline.pipelineOutput import run_output
        self._log("[INFO]", f"[output] Writing FITS maps for {source}.")
        meta   = self.key_handler.meta
        maps   = self.key_handler.get_maps()
        cubes  = self.key_handler.get_cubes()
        fname  = self._get_output_fname(source)
        params = self.source_handler.get_source_params(source)
        run_output(
            source=source,
            fname=fname,
            meta=meta,
            maps=maps,
            cubes=cubes,
            params=params,
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _get_output_fname(self, source: str) -> str:
        """Build the .ecsv output filename for *source*."""
        meta = self.key_handler.meta
        out_dir = meta.get("out_dir", "Output/")
        resolution = meta.get("resolution", "angular")
        target_res = meta.get("target_res", 27.0)

        if resolution == "angular":
            res_suffix = str(int(target_res)) + "as"
        elif resolution == "physical":
            res_suffix = str(int(target_res)) + "pc"
        else:
            res_suffix = "native"

        date_str = date.today().strftime("%Y_%m_%d")
        fname = os.path.join(out_dir, f"{source}_data_struct_{res_suffix}_{date_str}.ecsv")

        structure_creation = meta.get("structure_creation", "default")
        if "archive" in structure_creation and os.path.exists(fname):
            version = 1
            base = fname[:-5]
            while os.path.exists(f"{base}_v{version}.ecsv"):
                version += 1
            fname = f"{base}_v{version}.ecsv"

        return fname

    def _log(self, level: str, msg: str):
        if self.verbose:
            print(f"{level:<10} {msg}")

    def _print_summary(self):
        self._log("[INFO]", "--- Run summary ---")
        all_ok = True
        for source, ok in self.run_success.items():
            status = "OK" if ok else "FAILED"
            self._log("[INFO]", f"  {source}: {status}")
            if not ok:
                all_ok = False
        if all_ok:
            self._log("[INFO]", "All sources completed successfully.")
        else:
            self._log("[WARNING]", "Some sources failed — check errors above.")

    def __repr__(self):
        return (
            f"PipelineHandler(key_dir='{self.key_dir}', "
            f"sources={self.source_handler.all_sources()})"
        )
