"""
Tests for the PyStructure Pipeline package.
Run with:  pytest tests/ -v
"""

import sys
import pytest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


# ---------------------------------------------------------------------------
# KeyHandler
# ---------------------------------------------------------------------------

class TestKeyHandler:

    def _write_minimal_keys(self, tmpdir: Path):
        (tmpdir / "master_key.txt").write_text(
            "[paths]\n"
            f"data_dir    = {tmpdir}/data/\n"
            f"out_dir     = {tmpdir}/Output/\n"
            f"geom_file   = {tmpdir}/target_definitions.txt\n"
            f"data_key    = {tmpdir}/data_key.txt\n"
            f"config_key  = {tmpdir}/config_key.txt\n"
            "[meta]\nuser = Test\ncomments = test\n"
        )
        (tmpdir / "target_definitions.txt").write_text(
            "ngc5194\t202.4696\t47.1952\t8.58\t0.10\t22.0\t3.0\t173.0\t3.0\t3.54\t0.05\n"
        )
        (tmpdir / "data_key.txt").write_text(
            "[sources]\nsources = ngc5194\n"
            "[overlay]\noverlay_file = _12co21.fits\n"
            "# ---- maps ----\n"
            "spire250, SPIRE250, MJy/sr, _spire250.fits, data/\n"
            "# ---- cubes ----\n"
            "12co21, 12CO(2-1), K, _12co21.fits, data/\n"
            "# ---- mask ----\n"
        )
        (tmpdir / "config_key.txt").write_text(
            "[resolution]\ntarget_res = 27.0\nresolution = angular\n"
            "spacing_per_beam = 2\nmax_rad = auto\n"
            "NAXIS_shuff = 200\nCDELT_SHUFF = 4000\n"
            "[masking]\nref_line = first\nSN_processing = 2,4\n"
            "strict_mask = false\nuse_input_mask = false\n"
            "use_fixed_vel_mask = false\nuse_hfs_lines = false\n"
            "mom_thresh = 5\nconseq_channels = 3\nmom2_method = fwhm\n"
            "[spectral]\nspec_smooth = default\nspec_smooth_method = binned\n"
            "[output]\nsave_fits = false\nsave_mom_maps = true\n"
            "save_maps = true\nfolder_savefits = ./saved_FITS_files/\n"
            "[structure]\nstructure_creation = default\n"
        )

    def test_load_basic(self, tmp_path):
        self._write_minimal_keys(tmp_path)
        from pystructurePipeline.handlerKeys import KeyHandler
        kh = KeyHandler(str(tmp_path))
        assert kh.sources == ["ngc5194"]
        assert len(kh.maps) == 1
        assert len(kh.cubes) == 1
        assert kh.meta["target_res"] == 27.0

    def test_validate_passes(self, tmp_path):
        self._write_minimal_keys(tmp_path)
        from pystructurePipeline.handlerKeys import KeyHandler
        assert KeyHandler(str(tmp_path)).validate() is True

    def test_missing_key_dir_raises(self):
        from pystructurePipeline.handlerKeys import KeyHandler
        with pytest.raises(FileNotFoundError):
            KeyHandler("/nonexistent/path/")

    def test_repr(self, tmp_path):
        self._write_minimal_keys(tmp_path)
        from pystructurePipeline.handlerKeys import KeyHandler
        kh = KeyHandler(str(tmp_path))
        assert "KeyHandler" in repr(kh)
        assert "ngc5194" in repr(kh)

    def test_multi_source_list(self, tmp_path):
        self._write_minimal_keys(tmp_path)
        (tmp_path / "data_key.txt").write_text(
            "[sources]\nsources = ngc5194, ngc5457\n"
            "[overlay]\noverlay_file = _12co21.fits\n"
            "# ---- maps ----\nspire250, SPIRE250, MJy/sr, _s.fits, data/\n"
            "# ---- cubes ----\n12co21, 12CO(2-1), K, _12co21.fits, data/\n"
            "# ---- mask ----\n"
        )
        (tmp_path / "target_definitions.txt").write_text(
            "ngc5194\t202.4696\t47.1952\t8.58\t0.10\t22.0\t3.0\t173.0\t3.0\t3.54\t0.05\n"
            "ngc5457\t210.8025\t54.3492\t6.70\t0.32\t18.0\t5.0\t39.0\t5.0\t13.46\t0.50\n"
        )
        from pystructurePipeline.handlerKeys import KeyHandler
        kh = KeyHandler(str(tmp_path))
        assert kh.sources == ["ngc5194", "ngc5457"]


# ---------------------------------------------------------------------------
# SourceHandler
# ---------------------------------------------------------------------------

class TestSourceHandler:

    def _make_table(self):
        import pandas as pd
        return pd.DataFrame([{
            "source": "ngc5194", "ra_ctr": 202.47, "dec_ctr": 47.20,
            "dist_mpc": 8.58, "e_dist_mpc": 0.1, "incl_deg": 22.0,
            "e_incl_deg": 3.0, "posang_deg": 173.0, "e_posang_deg": 3.0,
            "r25": 3.54, "e_r25": 0.05,
        }])

    def test_get_source_params(self):
        from pystructurePipeline.handlerSources import SourceHandler
        th = SourceHandler(self._make_table(), ["ngc5194"])
        assert abs(th.get_source_params("ngc5194")["ra_ctr"] - 202.47) < 1e-6

    def test_unknown_source_raises(self):
        from pystructurePipeline.handlerSources import SourceHandler
        th = SourceHandler(self._make_table(), ["ngc5194"])
        with pytest.raises(KeyError):
            th.get_source_params("ngc9999")

    def test_source_not_in_table_raises(self):
        from pystructurePipeline.handlerSources import SourceHandler
        with pytest.raises(ValueError):
            SourceHandler(self._make_table(), ["ngc9999"])


# ---------------------------------------------------------------------------
# utils.fits_utils
# ---------------------------------------------------------------------------

class TestFitsUtils:

    def test_get_beam_arcsec_missing_file(self):
        from pystructurePipeline.utilsFits import get_beam_arcsec
        with pytest.raises(FileNotFoundError):
            get_beam_arcsec("/nonexistent/file.fits")

    def test_read_fits_cube_missing_file(self):
        from pystructurePipeline.utilsFits import read_fits_cube
        with pytest.raises(FileNotFoundError):
            read_fits_cube("/nonexistent/file.fits")

    def test_hex_grid_basic(self):
        from pystructurePipeline.utilsFits import hex_grid
        x, y = hex_grid(0.0, 0.0, 0.01, radec=False, r_limit=0.05)
        assert len(x) > 0

    def test_deproject_shape(self):
        import numpy as np
        from pystructurePipeline.utilsFits import deproject
        ra  = np.linspace(202.0, 203.0, 10)
        dec = np.linspace(47.0, 48.0, 10)
        r, t = deproject(ra, dec, [173.0, 22.0, 202.47, 47.20], vector=True)
        assert r.shape == ra.shape

    def test_gaussian_PSF_2D_shape(self):
        import numpy as np
        from pystructurePipeline.utilsFits import gaussian_PSF_2D
        psf = gaussian_PSF_2D(11, [0., 1., 3., 3., 0., 0., 0.], center=True, normalize=True)
        assert psf.shape == (11, 11)
        assert abs(np.sum(psf) - 1.0) < 1e-6

    def test_deconvolve_gauss_basic(self):
        from pystructurePipeline.utilsFits import deconvolve_gauss
        maj, minn, pa, info = deconvolve_gauss(30.0, 20.0, 30.0, 0.0, 20.0, 0.0)
        assert info[0]   # worked
        assert maj > 0


# ---------------------------------------------------------------------------
# utils.table_utils
# ---------------------------------------------------------------------------

class TestTableUtils:

    def test_load_missing_file(self):
        from pystructurePipeline.utilsTable import load_pystructure
        with pytest.raises(FileNotFoundError):
            load_pystructure("/nonexistent/file.ecsv")

    def test_find_latest_missing(self, tmp_path):
        from pystructurePipeline.utilsTable import find_latest_pystructure
        with pytest.raises(FileNotFoundError):
            find_latest_pystructure(str(tmp_path), "ngc5194")

    def test_shuffle_roundtrip(self):
        import numpy as np
        from pystructurePipeline.utilsTable import shuffle
        vaxis    = np.arange(-100, 101, 1.0)
        spec     = np.exp(-0.5 * (vaxis / 20.0) ** 2)
        shuffled = shuffle(spec, vaxis, zero=0.0, new_vaxis=vaxis)
        # Should be identical (same axis, zero shift)
        assert np.allclose(shuffled, spec, equal_nan=True)

    def test_get_mom_maps_runs(self):
        import numpy as np
        from astropy import units as au
        from pystructurePipeline.utilsTable import get_mom_maps

        n_pts, n_chan = 5, 50
        vaxis = np.linspace(-100, 100, n_chan) * au.km / au.s
        spec  = np.zeros((n_pts, n_chan)) * au.K
        # Put a Gaussian signal in one spectrum
        spec[2, :] = np.exp(-0.5 * (np.linspace(-100, 100, n_chan) / 15.0) ** 2) * au.K
        mask = (spec.value > 0.1).astype(float)

        moms = get_mom_maps(spec, mask, vaxis, mom_calc=[2, 3, "fwhm"])
        assert moms["mom0"].shape == (n_pts,)
        assert np.isfinite(moms["mom0"][2].value)


# ---------------------------------------------------------------------------
# init_workdir
# ---------------------------------------------------------------------------

class TestInitWorkdir:

    def test_creates_expected_files(self, tmp_path):
        from pystructurePipeline.init_workdir import init_workdir
        init_workdir(str(tmp_path))
        assert (tmp_path / "keys" / "master_key.txt").exists()
        assert (tmp_path / "keys" / "config_key.txt").exists()
        assert (tmp_path / "run_pystructure.py").exists()

    def test_overwrite_false_raises(self, tmp_path):
        from pystructurePipeline.init_workdir import init_workdir
        init_workdir(str(tmp_path))
        with pytest.raises(FileExistsError):
            init_workdir(str(tmp_path), overwrite=False)

    def test_overwrite_true_replaces(self, tmp_path):
        from pystructurePipeline.init_workdir import init_workdir
        init_workdir(str(tmp_path))
        (tmp_path / "run_pystructure.py").write_text("# corrupted")
        init_workdir(str(tmp_path), overwrite=True)
        assert "PipelineHandler" in (tmp_path / "run_pystructure.py").read_text()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

class TestCLI:

    def test_init_creates_files(self, tmp_path):
        from pystructurePipeline.cli import main
        main(["--init", "--workdir", str(tmp_path)])
        assert (tmp_path / "run_pystructure.py").exists()

    def test_init_overwrite_conflict(self, tmp_path):
        from pystructurePipeline.cli import main
        main(["--init", "--workdir", str(tmp_path)])
        with pytest.raises(SystemExit):
            main(["--init", "--workdir", str(tmp_path)])

    def test_missing_key_dir_exits(self):
        from pystructurePipeline.cli import main
        with pytest.raises((SystemExit, FileNotFoundError)):
            main(["--key_dir", "/nonexistent/"])

    def test_no_args_exits(self):
        from pystructurePipeline.cli import main
        with pytest.raises(SystemExit):
            main([])

    def test_invalid_stage_exits(self):
        from pystructurePipeline.cli import main
        with pytest.raises(SystemExit):
            main(["--key_dir", ".", "--stages", "invalid_stage"])


# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------

class TestLogger:

    def test_get_logger_prints_formatted_message(self, capsys):
        from pystructurePipeline.pystructureLogger import logger, get_logger
        logger.configure(verbose=True, log_file=None)
        log = get_logger("Regrid")
        log.info("hello world")
        captured = capsys.readouterr()
        assert "[pyStructure] [Regrid]  [INFO]  hello world" in captured.out

    def test_verbose_false_suppresses_print(self, capsys):
        from pystructurePipeline.pystructureLogger import logger, get_logger
        logger.configure(verbose=False, log_file=None)
        log = get_logger("Products")
        log.warning("should not print")
        captured = capsys.readouterr()
        assert captured.out == ""
        # but still recorded
        assert any(r["message"] == "should not print" for r in logger.get_records())
        logger.configure(verbose=True, log_file=None)  # restore default

    def test_log_file_written(self, tmp_path):
        from pystructurePipeline.pystructureLogger import logger, get_logger
        log_path = tmp_path / "run.log"
        logger.configure(verbose=False, log_file=str(log_path))
        log = get_logger("FITS")
        log.error("file not found")
        content = log_path.read_text()
        assert "[pyStructure] [FITS]  [ERROR]  file not found" in content
        logger.configure(verbose=True, log_file=None)  # restore default

    def test_save_writes_all_records(self, tmp_path):
        from pystructurePipeline.pystructureLogger import logger, get_logger
        logger.configure(verbose=False, log_file=None)
        log = get_logger("Sampling")
        log.info("a")
        log.warning("b")
        save_path = tmp_path / "saved.log"
        logger.save(str(save_path))
        content = save_path.read_text()
        assert "[Sampling]  [INFO]  a" in content
        assert "[Sampling]  [WARNING]  b" in content
        logger.configure(verbose=True, log_file=None)  # restore default

    def test_get_records_filtering(self):
        from pystructurePipeline.pystructureLogger import logger, get_logger
        logger.configure(verbose=False, log_file=None)
        log = get_logger("Keys")
        log.info("info msg")
        log.error("error msg")
        errors = logger.get_records(stage="Keys", level="ERROR")
        assert len(errors) >= 1
        assert all(r["level"] == "ERROR" for r in errors)
        logger.configure(verbose=True, log_file=None)  # restore default


class TestPipelineHandlerLogging:

    def test_log_file_created_on_init(self, tmp_path):
        """PipelineHandler(log_file=...) should create the log file immediately."""
        # Re-use the minimal key set from TestKeyHandler
        kh = TestKeyHandler()
        kh._write_minimal_keys(tmp_path)

        from pystructurePipeline.handlerPipeline import PipelineHandler
        log_path = tmp_path / "run.log"
        handler = PipelineHandler(key_dir=str(tmp_path), verbose=False,
                                  log_file=str(log_path))
        assert log_path.exists()
        content = log_path.read_text()
        assert "[pyStructure] [Pipeline]  [INFO]  Loading key files..." in content
