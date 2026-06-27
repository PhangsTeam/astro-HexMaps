"""
Tests for the HexMaps Pipeline package.
Run with:  pytest tests/ -v
"""

import sys
import pytest
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


# ---------------------------------------------------------------------------
# KeyHandler
# ---------------------------------------------------------------------------


class TestKeyHandler:

    def _write_minimal_config(self, tmpdir: Path) -> Path:
        """
        Write a minimal config.txt + keys/target_definitions.txt into *tmpdir*
        and return the path to config.txt.
        """
        keys_dir = tmpdir / "keys"
        keys_dir.mkdir(exist_ok=True)

        (keys_dir / "target_definitions.txt").write_text(
            "ngc5194,\t202.4696,  47.1952,\t8.58, 0.10,\t22.0, 3.0, 173.0, 3.0, 3.54, 0.05\n"
        )

        conf_path = tmpdir / "config.txt"
        conf_path.write_text(
            "[paths]\n"
            f"data_dir    = {tmpdir}/data/\n"
            f"out_dir     = {tmpdir}/Output/\n"
            "[meta]\nuser = Test\ncomments = test\n"
            "[sources]\nsources = ngc5194\n"
            "[overlay]\noverlay_file = _12co21.fits\n"
            "# ---- maps ----\n"
            "spire250, SPIRE250, MJy/sr, _spire250.fits, data/\n"
            "# ---- cubes ----\n"
            "12co21, 12CO(2-1), K, _12co21.fits, data/\n"
            "# ---- mask ----\n"
            # NOTE: [resolution]/[masking]/[spectral]/[output]/[structure]
            # are placed AFTER the maps/cubes/mask tables here, matching the
            # real config.txt template. This ordering is what previously
            # triggered the bug where every setting in these sections
            # silently fell back to its default — see
            # test_settings_after_tables_are_not_silently_dropped below.
            "[resolution]\ntarget_res = 27.0\nresolution = angular\n"
            "pixels_per_beam = 2\nmax_rad = auto\n"
            "NAXIS_shuff = 200\nCDELT_SHUFF = 4000\n"
            "[masking]\nref_line = first\nSN_processing = 2,4\n"
            "strict_mask = false\nuse_input_mask = false\n"
            "use_fixed_vel_mask = false\nuse_hfs_lines = false\n"
            "mom_thresh = 5\nconseq_channels = 3\n"
            "[spectral]\nspec_smooth = default\nspec_smooth_method = binned\n"
            "[output]\nsave_cubes = false\nsave_mom_maps = true\n"
            "save_maps = true\nfolder_savefits = ./saved_fits_files/\n"
            "[structure]\nstructure_creation = default\n"
        )
        return conf_path

    def test_load_basic(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        assert kh.sources == ["ngc5194"]
        assert len(kh.maps) == 1
        assert len(kh.cubes) == 1
        assert kh.meta["target_res"] == 27.0

    def test_resolve_resolution_sets_all_three_keys_angular(self, tmp_path):
        """
        After KeyHandler.load(), meta must contain target_res (arcsec),
        target_res_pc (parsec), and res_suffix for all three modes.
        Angular mode is the simplest case: no FITS file needed.
        """
        conf_path = self._write_minimal_config(tmp_path)
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        assert kh.meta["target_res"] == 27.0  # arcsec
        assert kh.meta["target_res_pc"] > 0  # parsec — positive
        assert kh.meta["res_suffix"] == "27p0as"

    def test_resolve_resolution_suffix_physical(self, tmp_path):
        """Physical mode: res_suffix ends in 'pc'."""
        conf_path = self._write_minimal_config(tmp_path)
        conf_path.write_text(
            conf_path.read_text()
            .replace("resolution = angular", "resolution = physical")
            .replace("target_res = 27.0", "target_res = 100.0")
        )
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        assert kh.meta["res_suffix"].endswith("pc")
        assert kh.meta["target_res"] > 0  # converted to arcsec
        assert kh.meta["target_res_pc"] == pytest.approx(100.0, rel=1e-3)

    def test_find_output_fname_discovers_existing_file(self, tmp_path):
        """
        _find_output_fname must return the most recent existing .ecsv
        rather than today's date when a matching file already exists.
        """
        import os
        from hexmaps.handler_pipeline import PipelineHandler

        handler = (
            PipelineHandler(conf_path=str(tmp_path / "config.txt"), verbose=False)
            if False
            else None
        )

        conf_path = self._write_minimal_config(tmp_path)
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        out_dir = kh.meta["out_dir"]
        res_suffix = kh.meta["res_suffix"]
        os.makedirs(out_dir, exist_ok=True)

        # Write a fake .ecsv with an older date
        existing = os.path.join(
            out_dir, f"ngc5194_hexmaps_{res_suffix}_2024_01_01.ecsv"
        )
        with open(existing, "w") as f:
            f.write("# fake ecsv\n")

        # _find_output_fname should discover this file without needing regrid
        from hexmaps.handler_pipeline import PipelineHandler

        handler = PipelineHandler(conf_path=str(conf_path), verbose=False)
        found = handler._find_output_fname("ngc5194")
        assert found == existing

    # ------------------------------------------------------------------
    # Mandatory-key error tests
    # ------------------------------------------------------------------

    @pytest.mark.parametrize(
        "missing_key,section,replace_line",
        [
            ("target_res", "resolution", "target_res = 27.0"),
            ("resolution", "resolution", "resolution = angular"),
            ("ref_line", "masking", "ref_line = first"),
        ],
    )
    def test_mandatory_setting_raises_on_missing_key(
        self, tmp_path, missing_key, section, replace_line
    ):
        """
        target_res, resolution, and ref_line have no fallback defaults.
        Omitting any of them must raise a KeyError with a clear message.
        """
        conf_path = self._write_minimal_config(tmp_path)
        text = conf_path.read_text().replace(replace_line + "\n", "")
        conf_path.write_text(text)

        from hexmaps.handler_keys import KeyHandler

        with pytest.raises(KeyError, match="Mandatory key"):
            KeyHandler(str(conf_path))

    def test_mandatory_sources_raises_on_missing(self, tmp_path):
        """
        [sources] / sources is mandatory; omitting it must raise a KeyError.
        """
        conf_path = self._write_minimal_config(tmp_path)
        text = conf_path.read_text().replace("[sources]\nsources = ngc5194\n", "")
        conf_path.write_text(text)

        from hexmaps.handler_keys import KeyHandler

        with pytest.raises(KeyError, match="Mandatory key"):
            KeyHandler(str(conf_path))

    def test_mandatory_overlay_file_raises_on_missing(self, tmp_path):
        """
        [overlay] / overlay_file is mandatory; omitting it must raise a KeyError.
        """
        conf_path = self._write_minimal_config(tmp_path)
        text = conf_path.read_text().replace(
            "[overlay]\noverlay_file = _12co21.fits\n", ""
        )
        conf_path.write_text(text)

        from hexmaps.handler_keys import KeyHandler

        with pytest.raises(KeyError, match="Mandatory key"):
            KeyHandler(str(conf_path))

    def test_optional_galaxy_columns_load_as_nan(self, tmp_path):
        """
        Galaxy geometry columns (incl_deg, posang_deg, r25 etc.) are optional.
        Rows that omit them should load with NaN rather than raising an error.
        """
        conf_path = self._write_minimal_config(tmp_path)

        # Write a target_definitions.txt with only the 4 mandatory columns
        geom = tmp_path / "keys" / "target_definitions.txt"
        geom.write_text("ngc5194, 202.4696, 47.1952, 8.58\n")

        from hexmaps.handler_keys import KeyHandler
        kh = KeyHandler(str(conf_path))
        row = kh.source_table.iloc[0]

        import math
        assert math.isnan(float(row["incl_deg"]))
        assert math.isnan(float(row["posang_deg"]))
        assert math.isnan(float(row["r25"]))

    def test_has_galaxy_geometry_false_when_nan(self, tmp_path):
        """has_galaxy_geometry() must return False when any column is NaN."""
        conf_path = self._write_minimal_config(tmp_path)
        geom = tmp_path / "keys" / "target_definitions.txt"
        geom.write_text("ngc5194, 202.4696, 47.1952, 8.58\n")

        from hexmaps.handler_keys import KeyHandler
        from hexmaps.handler_sources import SourceHandler
        kh = KeyHandler(str(conf_path))
        sh = SourceHandler(kh.source_table, kh.sources)
        assert sh.has_galaxy_geometry("ngc5194") is False

    def test_has_galaxy_geometry_true_when_full(self, tmp_path):
        """has_galaxy_geometry() must return True when all three columns are present."""
        conf_path = self._write_minimal_config(tmp_path)

        from hexmaps.handler_keys import KeyHandler
        from hexmaps.handler_sources import SourceHandler
        kh = KeyHandler(str(conf_path))
        sh = SourceHandler(kh.source_table, kh.sources)
        assert sh.has_galaxy_geometry("ngc5194") is True

    def test_save_mask_defaults_false(self, tmp_path):
        """save_mask is not set in the minimal fixture, so it must default to False."""
        conf_path = self._write_minimal_config(tmp_path)
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        assert kh.meta["save_mask"] is False

    def test_save_mask_explicit_true(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        conf_path.write_text(
            conf_path.read_text().replace(
                "[output]\nsave_cubes = false",
                "[output]\nsave_mask = true\nsave_cubes = false",
            )
        )
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        assert kh.meta["save_mask"] is True

    def test_settings_after_tables_are_not_silently_dropped(self, tmp_path):
        """
        Regression test: in the real config.txt template, [resolution] /
        [masking] / [spectral] / [output] / [structure] all come AFTER the
        maps/cubes/mask tables. A previous bug stopped feeding configparser
        at the first "# ---- maps ----" divider and never resumed, so every
        setting in those later sections silently fell back to its default
        no matter what the file said. Using a non-default value here makes
        sure that bug can't return unnoticed.
        """
        conf_path = self._write_minimal_config(tmp_path)
        conf_path.write_text(
            conf_path.read_text()
            .replace("target_res = 27.0", "target_res = 45.0")
            .replace("save_cubes = false", "save_cubes = true")
        )
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        assert kh.meta["target_res"] == 45.0
        assert kh.meta["save_cubes"] is True

    def test_target_res_written_back_after_resolution_block(self, tmp_path):
        """
        run_sampling must populate two keys in meta for all three modes:

          meta["target_res"]    → always arcseconds (single source of truth for math)
          meta["target_res_pc"] → always parsecs (for display/filenames in physical mode)

        For angular and native modes meta["target_res"] stays in arcseconds;
        for physical mode it is converted from the config parsec value.
        """
        import numpy as np
        from astropy.io import fits
        from hexmaps.stage_regrid import run_sampling

        # Minimal 3-D overlay cube with known BMAJ = 30 arcsec
        ny, nx, nv = 10, 10, 5
        cube = np.ones((nv, ny, nx))
        hdr = fits.Header()
        hdr["NAXIS"] = 3
        hdr["NAXIS1"], hdr["NAXIS2"], hdr["NAXIS3"] = nx, ny, nv
        hdr["CTYPE1"], hdr["CTYPE2"], hdr["CTYPE3"] = "RA---TAN", "DEC--TAN", "VELO"
        hdr["CRVAL1"], hdr["CRVAL2"], hdr["CRVAL3"] = 10.0, 20.0, 0.0
        hdr["CDELT1"], hdr["CDELT2"], hdr["CDELT3"] = -0.002, 0.002, 1000.0
        hdr["CRPIX1"], hdr["CRPIX2"], hdr["CRPIX3"] = 5, 5, 1
        hdr["CUNIT3"] = "m/s"
        hdr["BMAJ"] = 30.0 / 3600.0
        hdr["BMIN"] = 30.0 / 3600.0
        fits.writeto(str(tmp_path / "testsrc_12co21.fits"), cube, hdr)

        base_meta = {
            "data_dir": str(tmp_path),
            "out_dir": str(tmp_path),
            "overlay_file": "_12co21.fits",
            "pixels_per_beam": 2.0,
            "max_rad": "auto",
        }
        dist_mpc = 10.0
        params = {
            "ra_ctr": 10.0,
            "dec_ctr": 20.0,
            "dist_mpc": dist_mpc,
            "incl_deg": 0.0,
            "posang_deg": 0.0,
            "r25": 0.05,
        }

        # angular: target_res in arcsec; target_res_pc computed from distance
        meta = {**base_meta, "resolution": "angular", "target_res": 27.0}
        run_sampling("testsrc", params, meta)
        assert meta["target_res"] == 27.0
        expected_pc = 27.0 / 3600.0 * np.pi / 180.0 * dist_mpc * 1e6
        assert abs(meta["target_res_pc"] - expected_pc) < 0.01
        assert meta["res_suffix"] == "27p0as"

        # physical: target_res converted to arcsec; target_res_pc from distance
        meta = {**base_meta, "resolution": "physical", "target_res": 100.0}
        run_sampling("testsrc", params, meta)
        expected_as = 3600.0 * 180.0 / np.pi * 1e-6 * 100.0 / dist_mpc
        assert abs(meta["target_res"] - expected_as) < 0.01  # now in arcsec
        assert (
            abs(
                meta["target_res_pc"]
                - expected_as / 3600.0 * np.pi / 180.0 * dist_mpc * 1e6
            )
            < 0.01
        )
        assert meta["res_suffix"] == "100pc"

        # native: target_res = BMAJ = 30 arcsec; target_res_pc from distance
        meta = {**base_meta, "resolution": "native", "target_res": 99.0}
        run_sampling("testsrc", params, meta)
        assert meta["target_res"] == 30.0
        expected_pc_nat = 30.0 / 3600.0 * np.pi / 180.0 * dist_mpc * 1e6
        assert abs(meta["target_res_pc"] - expected_pc_nat) < 0.01
        assert meta["res_suffix"] == "30p0as"

    def test_fallback_use_logs_warning(self, tmp_path, capsys):
        """
        Whenever a [resolution]/[masking]/[spectral]/[output]/[structure]
        setting is absent from config.txt and its hardcoded default is used
        instead, a [WARNING] should be logged so this doesn't go unnoticed
        (the original motivation: a typo'd or misplaced setting should never
        silently and quietly fall back without a trace).
        """
        from hexmaps.logger import logger

        logger.configure(verbose=True, log_file=None)

        conf_path = self._write_minimal_config(tmp_path)
        from hexmaps.handler_keys import KeyHandler

        KeyHandler(str(conf_path))
        captured = capsys.readouterr()
        assert "[WARNING]" in captured.out
        # "mom2_method" is absent from the minimal config fixture's
        # [masking] section, so it must fall back and warn.
        assert "mom2_method" in captured.out
        assert "using default" in captured.out

    def test_fname_fill_fallback_does_not_warn(self, tmp_path, capsys):
        """
        fname_fill is an optional, rarely-used parameter (only relevant when
        structure_creation = "fill"), so its fallback to "" must NOT log a
        warning, unlike every other [resolution]/[masking]/[spectral]/
        [output]/[structure] setting.
        """
        from hexmaps.logger import logger

        logger.configure(verbose=True, log_file=None)

        conf_path = self._write_minimal_config(tmp_path)
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        captured = capsys.readouterr()
        assert kh.meta["fname_fill"] == ""
        assert "fname_fill" not in captured.out

    def test_explicit_setting_does_not_log_warning(self, tmp_path, capsys):
        """An explicitly-set value should not trigger a fallback warning."""
        from hexmaps.logger import logger

        logger.configure(verbose=True, log_file=None)

        conf_path = self._write_minimal_config(tmp_path)
        from hexmaps.handler_keys import KeyHandler

        KeyHandler(str(conf_path))
        captured = capsys.readouterr()
        # target_res IS set explicitly in the minimal config fixture, so it
        # must not appear in any fallback warning.
        assert "target_res not set" not in captured.out

    def test_target_definitions_ignores_whitespace_around_commas(self, tmp_path):
        """
        target_definitions.txt is comma-separated, but mixed tabs/spaces
        around each comma (for column alignment) must be ignored, and
        numeric columns must come back as floats, not whitespace-padded
        strings.
        """
        conf_path = self._write_minimal_config(tmp_path)
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        row = kh.source_table.iloc[0]
        assert row["source"] == "ngc5194"
        assert row["ra_ctr"] == 202.4696
        assert row["dec_ctr"] == 47.1952
        assert isinstance(row["dist_mpc"], float)

    def test_validate_passes(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        from hexmaps.handler_keys import KeyHandler

        assert KeyHandler(str(conf_path)).validate() is True

    def test_missing_conf_path_raises(self):
        from hexmaps.handler_keys import KeyHandler

        with pytest.raises(FileNotFoundError):
            KeyHandler("/nonexistent/path/config.txt")

    def test_missing_target_definitions_raises(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        (tmp_path / "keys" / "target_definitions.txt").unlink()
        from hexmaps.handler_keys import KeyHandler

        with pytest.raises(FileNotFoundError):
            KeyHandler(str(conf_path))

    def test_repr(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        assert "KeyHandler" in repr(kh)
        assert "ngc5194" in repr(kh)

    def test_multi_source_list(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        conf_path.write_text(
            conf_path.read_text().replace(
                "[sources]\nsources = ngc5194\n",
                "[sources]\nsources = ngc5194, ngc5457\n",
            )
        )
        (tmp_path / "keys" / "target_definitions.txt").write_text(
            "ngc5194, 202.4696, 47.1952, 8.58, 0.10, 22.0, 3.0, 173.0, 3.0, 3.54, 0.05\n"
            "ngc5457, 210.8025, 54.3492, 6.70, 0.32, 18.0, 5.0, 39.0, 5.0, 13.46, 0.50\n"
        )
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        assert kh.sources == ["ngc5194", "ngc5457"]

    def test_hfs_file_loaded_when_present(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        (tmp_path / "keys" / "hfs_lines.txt").write_text(
            "hcn10,\t88.6316023,  88.6304156,\tGHz\n"
        )
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        assert kh.hfs_data is not None
        assert len(kh.hfs_data) == 1
        row = kh.hfs_data.iloc[0]
        assert row["hfs_name"] == "hcn10"
        assert row["hfs_ref_freq"] == 88.6316023
        assert row["unit"] == "GHz"

    def test_hfs_file_none_when_absent(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        assert kh.hfs_data is None

    def test_geom_file_custom_path(self, tmp_path):
        """geom_file should be configurable via [paths], just like hfs_file."""
        conf_path = self._write_minimal_config(tmp_path)
        custom_geom = tmp_path / "shared" / "my_targets.txt"
        custom_geom.parent.mkdir(parents=True)
        custom_geom.write_text(
            "ngc1234, 10.0, 20.0, 5.0, 0.1, 30.0, 2.0, 90.0, 2.0, 2.0, 0.1\n"
        )
        conf_path.write_text(
            conf_path.read_text()
            .replace(
                "[paths]\n",
                f"[paths]\ngeom_file = {custom_geom}\n",
                1,
            )
            .replace(
                "[sources]\nsources = ngc5194\n",
                "[sources]\nsources = ngc1234\n",
            )
        )
        from hexmaps.handler_keys import KeyHandler

        kh = KeyHandler(str(conf_path))
        assert kh.sources == ["ngc1234"]
        assert "ngc1234" in list(kh.source_table["source"])

    def test_geom_file_missing_raises(self, tmp_path):
        """Unlike hfs_file, geom_file is required: a missing file must raise."""
        conf_path = self._write_minimal_config(tmp_path)
        (tmp_path / "keys" / "target_definitions.txt").unlink()
        from hexmaps.handler_keys import KeyHandler

        with pytest.raises(FileNotFoundError):
            KeyHandler(str(conf_path))

    def test_geom_file_missing_with_custom_path_raises(self, tmp_path):
        """A configured-but-nonexistent geom_file path must also raise."""
        conf_path = self._write_minimal_config(tmp_path)
        conf_path.write_text(
            conf_path.read_text().replace(
                "[paths]\n",
                "[paths]\ngeom_file = does_not_exist.txt\n",
                1,
            )
        )
        from hexmaps.handler_keys import KeyHandler

        with pytest.raises(FileNotFoundError):
            KeyHandler(str(conf_path))


# ---------------------------------------------------------------------------
# SourceHandler
# ---------------------------------------------------------------------------


class TestSourceHandler:

    def _make_table(self):
        import pandas as pd

        return pd.DataFrame(
            [
                {
                    "source": "ngc5194",
                    "ra_ctr": 202.47,
                    "dec_ctr": 47.20,
                    "dist_mpc": 8.58,
                    "e_dist_mpc": 0.1,
                    "incl_deg": 22.0,
                    "e_incl_deg": 3.0,
                    "posang_deg": 173.0,
                    "e_posang_deg": 3.0,
                    "r25": 3.54,
                    "e_r25": 0.05,
                }
            ]
        )

    def test_get_source_params(self):
        from hexmaps.handler_sources import SourceHandler

        th = SourceHandler(self._make_table(), ["ngc5194"])
        assert abs(th.get_source_params("ngc5194")["ra_ctr"] - 202.47) < 1e-6

    def test_unknown_source_raises(self):
        from hexmaps.handler_sources import SourceHandler

        th = SourceHandler(self._make_table(), ["ngc5194"])
        with pytest.raises(KeyError):
            th.get_source_params("ngc9999")

    def test_source_not_in_table_raises(self):
        from hexmaps.handler_sources import SourceHandler

        with pytest.raises(ValueError):
            SourceHandler(self._make_table(), ["ngc9999"])


# ---------------------------------------------------------------------------
# utils.fits_utils
# ---------------------------------------------------------------------------


class TestFitsUtils:

    def test_get_beam_arcsec_missing_file(self):
        from hexmaps.utils_fits import get_beam_arcsec

        with pytest.raises(FileNotFoundError):
            get_beam_arcsec("/nonexistent/file.fits")

    def test_read_fits_cube_missing_file(self):
        from hexmaps.utils_fits import read_fits_cube

        with pytest.raises(FileNotFoundError):
            read_fits_cube("/nonexistent/file.fits")

    def test_hex_grid_basic(self):
        from hexmaps.utils_fits import hex_grid

        x, y = hex_grid(0.0, 0.0, 0.01, radec=False, r_limit=0.05)
        assert len(x) > 0

    def test_deproject_shape(self):
        import numpy as np
        from hexmaps.utils_fits import deproject

        ra = np.linspace(202.0, 203.0, 10)
        dec = np.linspace(47.0, 48.0, 10)
        r, t = deproject(ra, dec, [173.0, 22.0, 202.47, 47.20], vector=True)
        assert r.shape == ra.shape

    def test_gaussian_PSF_2D_shape(self):
        import numpy as np
        from hexmaps.utils_fits import gaussian_PSF_2D

        psf = gaussian_PSF_2D(
            11, [0.0, 1.0, 3.0, 3.0, 0.0, 0.0, 0.0], center=True, normalize=True
        )
        assert psf.shape == (11, 11)
        assert abs(np.sum(psf) - 1.0) < 1e-6

    def test_deconvolve_gauss_basic(self):
        from hexmaps.utils_fits import deconvolve_gauss

        maj, minn, pa, info = deconvolve_gauss(30.0, 20.0, 30.0, 0.0, 20.0, 0.0)
        assert info[0]  # worked
        assert maj > 0


# ---------------------------------------------------------------------------
# stage_fits — mask cube output
# ---------------------------------------------------------------------------


class TestStageFits:

    def test_get_coord_names_radec(self):
        """RA/Dec CTYPE values produce RA / DEC column names."""
        from astropy.io import fits
        from hexmaps.stage_regrid import _get_coord_names
        hdr = fits.Header()
        hdr["CTYPE1"] = "RA---TAN"
        hdr["CTYPE2"] = "DEC--TAN"
        c1, c2, d1, d2 = _get_coord_names(hdr)
        assert c1 == "RA"
        assert c2 == "DEC"
        assert "Right ascension" in d1
        assert "Declination" in d2

    def test_get_coord_names_galactic(self):
        """Galactic CTYPE values produce GLON / GLAT column names."""
        from astropy.io import fits
        from hexmaps.stage_regrid import _get_coord_names
        hdr = fits.Header()
        hdr["CTYPE1"] = "GLON-CAR"
        hdr["CTYPE2"] = "GLAT-CAR"
        c1, c2, d1, d2 = _get_coord_names(hdr)
        assert c1 == "GLON"
        assert c2 == "GLAT"
        assert "Galactic longitude" in d1
        assert "Galactic latitude" in d2

    def test_build_edge_mask_erodes_footprint(self):
        """
        build_edge_mask must return a 2-D float array smaller than the input
        footprint, eroded by the expected number of pixels on all sides.
        """
        from hexmaps.stage_fits import build_edge_mask
        from astropy.io import fits

        # Full 30x30 observed footprint
        ov_footprint = np.ones((30, 30), dtype=bool)
        hdr = fits.Header()
        hdr["CDELT1"] = -1.0 / 3600.0  # 1 arcsec/pixel

        # target_res = 10 arcsec -> trim radius = floor(10/2 / 1) = 5 pixels
        edge_mask = build_edge_mask(ov_footprint, hdr, target_res_as=10.0)
        assert edge_mask.shape == (30, 30)
        # Centre must be kept, corners must be removed
        assert edge_mask[15, 15] == 1.0
        assert edge_mask[0, 0] == 0.0
        assert edge_mask[29, 29] == 0.0
        # Total kept pixels must be strictly fewer than full 30x30
        assert edge_mask.sum() < 30 * 30

    def test_build_edge_mask_zero_radius_returns_footprint(self):
        """
        When the pixel scale is coarser than half the beam (trim radius < 1),
        no trimming is applied and the full footprint is returned.
        """
        from hexmaps.stage_fits import build_edge_mask
        from astropy.io import fits

        ov_footprint = np.ones((10, 10), dtype=bool)
        hdr = fits.Header()
        hdr["CDELT1"] = -10.0 / 3600.0  # 10 arcsec/pixel — half beam is < 1 px

        edge_mask = build_edge_mask(ov_footprint, hdr, target_res_as=10.0)
        assert edge_mask.sum() == 100  # full footprint returned

    def test_build_edge_mask_erodes_non_rectangular_blob(self):
        """
        build_edge_mask must erode the irregular non-NaN blob defined by
        ov_footprint, not a rectangular grid extent. A circular island
        of observed pixels should be eroded at its own boundary.
        """
        from hexmaps.stage_fits import build_edge_mask
        from astropy.io import fits

        ny, nx = 20, 20
        # Circular island: observed area is a disc of radius 8 at centre
        y, x = np.ogrid[:ny, :nx]
        ov_footprint = ((y - 10) ** 2 + (x - 10) ** 2) < 8**2

        hdr = fits.Header()
        hdr["CDELT1"] = -1.0 / 3600.0  # 1 arcsec/pixel

        # Erode by 2 pixels
        edge_mask = build_edge_mask(ov_footprint, hdr, target_res_as=4.0)

        # Pixels outside the overlay footprint must never be in the edge mask
        assert np.all(edge_mask[~ov_footprint] == 0.0)
        # Centre of the disc (far from boundary) must be kept
        assert edge_mask[10, 10] == 1.0
        # Pixels on the disc boundary must be removed by erosion
        assert edge_mask[2, 10] == 0.0  # top of disc

    def test_overlay_footprint_constrains_erosion(self):
        """
        Passing an ov_footprint that is smaller than the full grid (as in
        the real pipeline, where the overlay has NaN outside the observed
        area) must constrain the edge mask to that footprint. Pixels outside
        the overlay footprint must be 0, and the eroded boundary must follow
        the overlay footprint boundary, not the full grid boundary.
        """
        from hexmaps.stage_fits import build_edge_mask
        from astropy.io import fits

        ny, nx = 20, 20
        # Only the centre 10x10 block is observed
        ov_footprint = np.zeros((ny, nx), dtype=bool)
        ov_footprint[5:15, 5:15] = True

        hdr = fits.Header()
        hdr["CDELT1"] = -1.0 / 3600.0
        # Erode by 2 pixels (target_res=4 arcsec at 1 arcsec/px)
        edge_mask = build_edge_mask(ov_footprint, hdr, target_res_as=4.0)

        # Outside overlay footprint: must be zero
        assert np.all(edge_mask[~ov_footprint] == 0.0)
        # Well inside the centre block: must be kept
        assert edge_mask[10, 10] == 1.0
        # Pixels at the overlay footprint boundary must be removed by erosion
        assert edge_mask[5, 10] == 0.0
        assert edge_mask[14, 10] == 0.0

    def test_build_edge_mask_fov_erosion_beams_scales_radius(self):
        """
        fov_erosion_beams scales the trim radius linearly: 1.0 beam gives
        twice the trim radius of 0.5 beam.
        """
        from hexmaps.stage_fits import build_edge_mask
        from astropy.io import fits

        ov_footprint = np.ones((40, 40), dtype=bool)
        hdr = fits.Header()
        hdr["CDELT1"] = -1.0 / 3600.0  # 1 arcsec/pixel, beam = 10 arcsec

        # 0.5 beam → trim 5 px; 1.0 beam → trim 10 px
        mask_half = build_edge_mask(ov_footprint, hdr, 10.0, fov_erosion_beams=0.5)
        mask_full = build_edge_mask(ov_footprint, hdr, 10.0, fov_erosion_beams=1.0)

        # full-beam erosion must keep fewer pixels than half-beam
        assert mask_full.sum() < mask_half.sum()
        # centre must be kept in both cases
        assert mask_half[20, 20] == 1.0
        assert mask_full[20, 20] == 1.0

    def test_build_edge_mask_zero_erosion_returns_full_footprint(self):
        """
        fov_erosion_beams = 0 must disable erosion and return the full
        footprint unchanged, without a warning about <= 0 pixels.
        """
        from hexmaps.stage_fits import build_edge_mask
        from astropy.io import fits

        ov_footprint = np.ones((15, 15), dtype=bool)
        hdr = fits.Header()
        hdr["CDELT1"] = -1.0 / 3600.0

        mask = build_edge_mask(ov_footprint, hdr, 10.0, fov_erosion_beams=0.0)
        assert mask.sum() == 15 * 15  # full footprint, no pixels removed

    def test_save_ppv_mask_to_fits_writes_cube_as_is(self, tmp_path):
        """
        save_ppv_mask_to_fits writes a plain numpy mask array directly,
        with no regridding/resampling/re-binarization step (the mask is
        already on the overlay's native PPV grid, unlike the old hex-grid
        save_to_fits_cube path).
        """
        from hexmaps.stage_fits import save_ppv_mask_to_fits
        from astropy.io import fits

        mask = np.zeros((4, 5, 5))
        mask[1:3, 2, 2] = 1

        ov_hdr = fits.Header()
        ov_hdr["NAXIS"] = 3
        ov_hdr["NAXIS1"], ov_hdr["NAXIS2"], ov_hdr["NAXIS3"] = 5, 5, 4
        ov_hdr["CTYPE1"], ov_hdr["CTYPE2"], ov_hdr["CTYPE3"] = (
            "RA---TAN",
            "DEC--TAN",
            "VELO",
        )
        ov_hdr["CRVAL3"], ov_hdr["CDELT3"], ov_hdr["CRPIX3"] = 0.0, 1000.0, 1

        save_ppv_mask_to_fits(mask, ov_hdr, "testsrc", "mask", str(tmp_path))

        out_path = tmp_path / "testsrc_mask.fits"
        assert out_path.exists()
        data, out_hdr = fits.getdata(str(out_path), header=True)
        assert np.array_equal(data, mask)
        assert out_hdr["NAXIS3"] == 4
        assert out_hdr["CRVAL3"] == 0.0

    # -----------------------------------------------------------------
    # PPV-native moment pipeline
    # -----------------------------------------------------------------

    def _make_synthetic_ppv_cube(self, n_chan=40, ny=4, nx=4, seed=0):
        """
        Build a small synthetic (n_chan, ny, nx) cube with a clean Gaussian
        emission line (high S/N, centred at channel 20) plus unit-variance
        noise, identical at every spatial pixel. Returns (cube, vaxis_kms).
        """
        import numpy as np

        rng = np.random.RandomState(seed)
        vaxis = np.arange(n_chan, dtype=float)  # channel index stands in for km/s
        cube = rng.normal(0, 1.0, size=(n_chan, ny, nx))
        line = 20 * np.exp(-0.5 * ((vaxis - 20) / 1.7) ** 2)
        cube += line[:, None, None]
        return cube, vaxis

    def test_construct_mask_ppv_recovers_known_line(self):
        """
        construct_mask_ppv should mask channels around the injected line
        centre and leave most far-from-line channels unmasked, for every
        spatial pixel (since the synthetic cube is spatially uniform).
        """
        from hexmaps.stage_fits import construct_mask_ppv

        cube, vaxis = self._make_synthetic_ppv_cube()
        mask = construct_mask_ppv(cube, SN_processing=[2, 4])

        assert mask.shape == cube.shape
        assert set(np.unique(mask)).issubset({0, 1})

        # The line centre (channel 20) should be masked everywhere
        assert np.all(mask[20] == 1)
        # Far from the line (channel 0) should be unmasked everywhere
        assert np.all(mask[0] == 0)

    def test_construct_mask_ppv_matches_hex_grid_construct_mask(self):
        """
        construct_mask_ppv must reproduce stage_products.construct_mask's
        mask exactly when run on the same underlying data, just reshaped:
        one hex-grid "point" per spatial pixel of the PPV cube.
        """
        import astropy.units as au
        from astropy.table import Table, Column
        from hexmaps.stage_fits import construct_mask_ppv
        from hexmaps.stage_products import construct_mask

        cube, vaxis = self._make_synthetic_ppv_cube(ny=3, nx=3)
        n_chan, ny, nx = cube.shape

        # Build the equivalent hex-grid table: one row per spatial pixel
        spec = np.moveaxis(cube, 0, -1).reshape(ny * nx, n_chan)
        t = Table()
        t["SPEC_LINE"] = Column(spec)
        t.meta["SPEC_VCHAN0"] = 0.0 * au.km / au.s
        t.meta["SPEC_DELTAV"] = 1.0 * au.km / au.s
        t.meta["SPEC_CRPIX"] = 1

        mask_ppv = construct_mask_ppv(cube, SN_processing=[2, 4])
        mask_hex, _, _ = construct_mask("LINE", t, SN_processing=[2, 4])

        mask_hex_reshaped = np.moveaxis(mask_hex.value.reshape(ny, nx, n_chan), -1, 0)
        assert np.array_equal(mask_ppv, mask_hex_reshaped)

    def test_apply_strict_mask_ppv_removes_small_components(self):
        from hexmaps.stage_fits import apply_strict_mask_ppv

        mask = np.zeros((1, 10, 10), dtype=int)
        mask[0, 5, 5] = 1  # isolated single pixel: too small, removed
        mask[0, 0:3, 0:3] = 1  # 3x3 block = 9 pixels: kept

        filtered = apply_strict_mask_ppv(mask, min_pixels=5)
        assert filtered[0, 5, 5] == 0
        assert np.all(filtered[0, 0:3, 0:3] == 1)

    def test_get_mom_maps_ppv_matches_get_mom_maps(self):
        """
        get_mom_maps_ppv must return the literal same values as calling
        utils_table.get_mom_maps directly on the reshaped (n_pix, n_chan)
        array -- it is a pure reshape wrapper, not a re-implementation.
        """
        import astropy.units as au
        from hexmaps.stage_fits import get_mom_maps_ppv
        from hexmaps.utils_table import get_mom_maps

        cube, vaxis_arr = self._make_synthetic_ppv_cube(ny=2, nx=2)
        n_chan, ny, nx = cube.shape

        cube_q = cube * au.K
        vaxis_q = vaxis_arr * au.km / au.s
        mask = (cube > 3).astype(int)
        # widen the mask a bit so get_mom_maps' high-S/N submask has enough
        # consecutive channels to compute mom1/mom2 (mirrors construct_mask's
        # dilation in spirit, simplified for this unit test)
        for shift in (1, 2, -1, -2):
            mask = np.maximum(mask, np.roll(cube, shift, axis=0) > 3)
        mom_calc = (3, 3, "fwhm")

        ppv_maps = get_mom_maps_ppv(cube_q, mask, vaxis_q, mom_calc)

        cube_pts = np.moveaxis(cube, 0, -1).reshape(ny * nx, n_chan) * au.K
        mask_pts = np.moveaxis(mask, 0, -1).reshape(ny * nx, n_chan)
        flat_maps = get_mom_maps(cube_pts, mask_pts, vaxis_q, mom_calc)

        for key in ppv_maps:
            assert ppv_maps[key].shape == (ny, nx)
            np.testing.assert_allclose(
                ppv_maps[key].value.ravel(),
                flat_maps[key].value,
                equal_nan=True,
            )

    def test_convolve_cube_to_target_skips_when_already_at_resolution(self):
        from hexmaps.stage_fits import convolve_cube_to_target
        from astropy.io import fits

        cube = np.ones((5, 6, 6))
        hdr = fits.Header()
        hdr["BMAJ"] = 30.0 / 3600.0  # already coarser than the 27" target
        hdr["BMIN"] = 30.0 / 3600.0

        out_data, out_hdr = convolve_cube_to_target(cube, hdr, target_res_as=27.0)
        assert np.array_equal(out_data, cube)

    def test_get_convolved_ppv_cube_convolves_from_raw_input(self, tmp_path):
        """
        get_convolved_ppv_cube must always read the raw input file and
        convolve from scratch — there is no cache lookup.
        """
        from hexmaps.stage_fits import get_convolved_ppv_cube
        from astropy.io import fits

        ny, nx, nv = 3, 3, 2
        cube = np.ones((nv, ny, nx))
        hdr = fits.Header()
        hdr["NAXIS"] = 3
        hdr["NAXIS1"], hdr["NAXIS2"], hdr["NAXIS3"] = nx, ny, nv
        hdr["CTYPE1"], hdr["CTYPE2"], hdr["CTYPE3"] = "RA---TAN", "DEC--TAN", "VELO"
        hdr["CRVAL1"], hdr["CRVAL2"], hdr["CRVAL3"] = 10.0, 20.0, 0.0
        hdr["CDELT1"], hdr["CDELT2"], hdr["CDELT3"] = -0.01, 0.01, 1000.0
        hdr["CRPIX1"], hdr["CRPIX2"], hdr["CRPIX3"] = 2, 2, 1
        hdr["CUNIT3"] = "m/s"
        hdr["BMAJ"] = 30.0 / 3600.0
        hdr["BMIN"] = 30.0 / 3600.0
        raw_path = tmp_path / "testsrc_co.fits"
        fits.writeto(str(raw_path), cube, hdr)

        meta = {
            "target_res": 27.0,
            "target_res_pc": 1000.0,
            "resolution": "angular",
            "res_suffix": "27p0as",
        }
        data, _ = get_convolved_ppv_cube(
            "testsrc",
            "co",
            str(tmp_path),
            "_co.fits",
            meta,
            hdr,
        )
        assert data.shape == (nv, ny, nx)

    def test_get_convolved_ppv_cube_raises_if_nothing_available(self, tmp_path):
        from hexmaps.stage_fits import get_convolved_ppv_cube
        from astropy.io import fits

        hdr = fits.Header()
        hdr["NAXIS"] = 3
        hdr["NAXIS1"], hdr["NAXIS2"], hdr["NAXIS3"] = 3, 3, 2

        meta = {
            "target_res": 27.0,
            "target_res_pc": 1000.0,
            "resolution": "angular",
            "res_suffix": "27p0as",
        }
        with pytest.raises(FileNotFoundError):
            get_convolved_ppv_cube(
                "testsrc",
                "co",
                str(tmp_path),
                ".fits",
                meta,
                hdr,
            )


# ---------------------------------------------------------------------------
# utils.table_utils
# ---------------------------------------------------------------------------


class TestTableUtils:

    def test_load_missing_file(self):
        from hexmaps.utils_table import load_hexmaps

        with pytest.raises(FileNotFoundError):
            load_hexmaps("/nonexistent/file.ecsv")

    def test_find_latest_missing(self, tmp_path):
        from hexmaps.utils_table import find_latest_hexmaps

        with pytest.raises(FileNotFoundError):
            find_latest_hexmaps(str(tmp_path), "ngc5194")

    def test_shuffle_roundtrip(self):
        import numpy as np
        from hexmaps.utils_table import shuffle

        vaxis = np.arange(-100, 101, 1.0)
        spec = np.exp(-0.5 * (vaxis / 20.0) ** 2)
        shuffled = shuffle(spec, vaxis, zero=0.0, new_vaxis=vaxis)
        # Should be identical (same axis, zero shift)
        assert np.allclose(shuffled, spec, equal_nan=True)

    def test_get_mom_maps_runs(self):
        import numpy as np
        from astropy import units as u
        from hexmaps.utils_table import get_mom_maps

        n_pts, n_chan = 5, 50
        vaxis = np.linspace(-100, 100, n_chan) * u.km / u.s
        spec = np.zeros((n_pts, n_chan)) * u.K
        # Put a Gaussian signal in one spectrum
        spec[2, :] = np.exp(-0.5 * (np.linspace(-100, 100, n_chan) / 15.0) ** 2) * u.K
        mask = (spec.value > 0.1).astype(float)

        moms = get_mom_maps(spec, mask, vaxis, mom_calc=[2, 3, "fwhm"])
        assert moms["mom0"].shape == (n_pts,)
        assert np.isfinite(moms["mom0"][2].value)


# ---------------------------------------------------------------------------
# init_workdir
# ---------------------------------------------------------------------------


class TestInitWorkdir:

    def test_creates_expected_files(self, tmp_path):
        from hexmaps.init_workdir import init_workdir

        init_workdir(str(tmp_path))
        assert (tmp_path / "config.txt").exists()
        assert (tmp_path / "keys" / "target_definitions.txt").exists()
        assert (tmp_path / "run_hexmaps.py").exists()

    def test_overwrite_false_raises(self, tmp_path):
        from hexmaps.init_workdir import init_workdir

        init_workdir(str(tmp_path))
        with pytest.raises(FileExistsError):
            init_workdir(str(tmp_path), overwrite=False)

    def test_overwrite_true_replaces(self, tmp_path):
        from hexmaps.init_workdir import init_workdir

        init_workdir(str(tmp_path))
        (tmp_path / "run_hexmaps.py").write_text("# corrupted")
        init_workdir(str(tmp_path), overwrite=True)
        assert "PipelineHandler" in (tmp_path / "run_hexmaps.py").read_text()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


class TestCLI:

    def test_init_creates_files(self, tmp_path):
        from hexmaps.cli import main

        main(["--init", "--workdir", str(tmp_path)])
        assert (tmp_path / "run_hexmaps.py").exists()

    def test_init_overwrite_conflict(self, tmp_path):
        from hexmaps.cli import main

        main(["--init", "--workdir", str(tmp_path)])
        with pytest.raises(SystemExit):
            main(["--init", "--workdir", str(tmp_path)])

    def test_missing_conf_exits(self):
        from hexmaps.cli import main

        with pytest.raises((SystemExit, FileNotFoundError)):
            main(["--conf", "/nonexistent/config.txt"])

    def test_no_args_exits(self):
        from hexmaps.cli import main

        with pytest.raises(SystemExit):
            main([])

    def test_invalid_stage_exits(self):
        from hexmaps.cli import main

        with pytest.raises(SystemExit):
            main(["--conf", "config.txt", "--stages", "invalid_stage"])


# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------


class TestLogger:

    def test_get_logger_prints_formatted_message(self, capsys):
        from hexmaps.logger import logger, get_logger

        logger.configure(verbose=True, log_file=None)
        log = get_logger("Regrid")
        log.info("hello world")
        captured = capsys.readouterr()
        assert "[HexMaps] [Regrid]    [INFO]     hello world" in captured.out

    def test_verbose_false_suppresses_print(self, capsys):
        from hexmaps.logger import logger, get_logger

        logger.configure(verbose=False, log_file=None)
        log = get_logger("Products")
        log.warning("should not print")
        captured = capsys.readouterr()
        assert captured.out == ""
        # but still recorded
        assert any(r["message"] == "should not print" for r in logger.get_records())
        logger.configure(verbose=True, log_file=None)  # restore default

    def test_log_file_written(self, tmp_path):
        from hexmaps.logger import logger, get_logger

        log_path = tmp_path / "run.log"
        logger.configure(verbose=False, log_file=str(log_path))
        log = get_logger("FITS")
        log.error("file not found")
        content = log_path.read_text()
        assert "[HexMaps] [FITS]      [ERROR]    file not found" in content
        logger.configure(verbose=True, log_file=None)  # restore default

    def test_save_writes_all_records(self, tmp_path):
        from hexmaps.logger import logger, get_logger

        logger.configure(verbose=False, log_file=None)
        log = get_logger("Sampling")
        log.info("a")
        log.warning("b")
        save_path = tmp_path / "saved.log"
        logger.save(str(save_path))
        content = save_path.read_text()
        assert "[Sampling]  [INFO]     a" in content
        assert "[Sampling]  [WARNING]  b" in content
        logger.configure(verbose=True, log_file=None)  # restore default

    def test_get_records_filtering(self):
        from hexmaps.logger import logger, get_logger

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
        # Re-use the minimal config from TestKeyHandler
        kh = TestKeyHandler()
        conf_path = kh._write_minimal_config(tmp_path)

        from hexmaps.handler_pipeline import PipelineHandler

        log_path = tmp_path / "run.log"
        handler = PipelineHandler(
            conf_path=str(conf_path), verbose=False, log_file=str(log_path)
        )
        assert log_path.exists()
        content = log_path.read_text()
        assert "[Loading]" in content
        assert "Loading configuration..." in content
