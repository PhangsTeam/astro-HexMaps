"""
Tests for the PyStructure Pipeline package.
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
            "[output]\nsave_fits = false\nsave_mom_maps = true\n"
            "save_maps = true\nfolder_savefits = ./saved_fits_files/\n"
            "[structure]\nstructure_creation = default\n"
        )
        return conf_path

    def test_load_basic(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        from pystructurePipeline.handler_keys import KeyHandler
        kh = KeyHandler(str(conf_path))
        assert kh.sources == ["ngc5194"]
        assert len(kh.maps) == 1
        assert len(kh.cubes) == 1
        assert kh.meta["target_res"] == 27.0

    def test_save_mask_defaults_false(self, tmp_path):
        """save_mask is not set in the minimal fixture, so it must default to False."""
        conf_path = self._write_minimal_config(tmp_path)
        from pystructurePipeline.handler_keys import KeyHandler
        kh = KeyHandler(str(conf_path))
        assert kh.meta["save_mask"] is False

    def test_save_mask_explicit_true(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        conf_path.write_text(
            conf_path.read_text().replace(
                "[output]\nsave_fits = false",
                "[output]\nsave_mask = true\nsave_fits = false",
            )
        )
        from pystructurePipeline.handler_keys import KeyHandler
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
            conf_path.read_text().replace(
                "target_res = 27.0", "target_res = 45.0"
            ).replace(
                "save_fits = false", "save_fits = true"
            )
        )
        from pystructurePipeline.handler_keys import KeyHandler
        kh = KeyHandler(str(conf_path))
        assert kh.meta["target_res"] == 45.0
        assert kh.meta["save_fits"] is True

    def test_fallback_use_logs_warning(self, tmp_path, capsys):
        """
        Whenever a [resolution]/[masking]/[spectral]/[output]/[structure]
        setting is absent from config.txt and its hardcoded default is used
        instead, a [WARNING] should be logged so this doesn't go unnoticed
        (the original motivation: a typo'd or misplaced setting should never
        silently and quietly fall back without a trace).
        """
        from pystructurePipeline.pystructureLogger import logger
        logger.configure(verbose=True, log_file=None)

        conf_path = self._write_minimal_config(tmp_path)
        from pystructurePipeline.handler_keys import KeyHandler
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
        from pystructurePipeline.pystructureLogger import logger
        logger.configure(verbose=True, log_file=None)

        conf_path = self._write_minimal_config(tmp_path)
        from pystructurePipeline.handler_keys import KeyHandler
        kh = KeyHandler(str(conf_path))
        captured = capsys.readouterr()
        assert kh.meta["fname_fill"] == ""
        assert "fname_fill" not in captured.out

    def test_explicit_setting_does_not_log_warning(self, tmp_path, capsys):
        """An explicitly-set value should not trigger a fallback warning."""
        from pystructurePipeline.pystructureLogger import logger
        logger.configure(verbose=True, log_file=None)

        conf_path = self._write_minimal_config(tmp_path)
        from pystructurePipeline.handler_keys import KeyHandler
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
        from pystructurePipeline.handler_keys import KeyHandler
        kh = KeyHandler(str(conf_path))
        row = kh.source_table.iloc[0]
        assert row["source"] == "ngc5194"
        assert row["ra_ctr"] == 202.4696
        assert row["dec_ctr"] == 47.1952
        assert isinstance(row["dist_mpc"], float)

    def test_validate_passes(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        from pystructurePipeline.handler_keys import KeyHandler
        assert KeyHandler(str(conf_path)).validate() is True

    def test_missing_conf_path_raises(self):
        from pystructurePipeline.handler_keys import KeyHandler
        with pytest.raises(FileNotFoundError):
            KeyHandler("/nonexistent/path/config.txt")

    def test_missing_target_definitions_raises(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        (tmp_path / "keys" / "target_definitions.txt").unlink()
        from pystructurePipeline.handler_keys import KeyHandler
        with pytest.raises(FileNotFoundError):
            KeyHandler(str(conf_path))

    def test_repr(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        from pystructurePipeline.handler_keys import KeyHandler
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
        from pystructurePipeline.handler_keys import KeyHandler
        kh = KeyHandler(str(conf_path))
        assert kh.sources == ["ngc5194", "ngc5457"]

    def test_hfs_file_loaded_when_present(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        (tmp_path / "keys" / "hfs_lines.txt").write_text(
            "hcn10,\t88.6316023,  88.6304156,\tGHz\n"
        )
        from pystructurePipeline.handler_keys import KeyHandler
        kh = KeyHandler(str(conf_path))
        assert kh.hfs_data is not None
        assert len(kh.hfs_data) == 1
        row = kh.hfs_data.iloc[0]
        assert row["hfs_name"] == "hcn10"
        assert row["hfs_ref_freq"] == 88.6316023
        assert row["unit"] == "GHz"

    def test_hfs_file_none_when_absent(self, tmp_path):
        conf_path = self._write_minimal_config(tmp_path)
        from pystructurePipeline.handler_keys import KeyHandler
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
            conf_path.read_text().replace(
                "[paths]\n",
                f"[paths]\ngeom_file = {custom_geom}\n",
                1,
            ).replace(
                "[sources]\nsources = ngc5194\n",
                "[sources]\nsources = ngc1234\n",
            )
        )
        from pystructurePipeline.handler_keys import KeyHandler
        kh = KeyHandler(str(conf_path))
        assert kh.sources == ["ngc1234"]
        assert "ngc1234" in list(kh.source_table["source"])

    def test_geom_file_missing_raises(self, tmp_path):
        """Unlike hfs_file, geom_file is required: a missing file must raise."""
        conf_path = self._write_minimal_config(tmp_path)
        (tmp_path / "keys" / "target_definitions.txt").unlink()
        from pystructurePipeline.handler_keys import KeyHandler
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
        from pystructurePipeline.handler_keys import KeyHandler
        with pytest.raises(FileNotFoundError):
            KeyHandler(str(conf_path))


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
        from pystructurePipeline.handler_sources import SourceHandler
        th = SourceHandler(self._make_table(), ["ngc5194"])
        assert abs(th.get_source_params("ngc5194")["ra_ctr"] - 202.47) < 1e-6

    def test_unknown_source_raises(self):
        from pystructurePipeline.handler_sources import SourceHandler
        th = SourceHandler(self._make_table(), ["ngc5194"])
        with pytest.raises(KeyError):
            th.get_source_params("ngc9999")

    def test_source_not_in_table_raises(self):
        from pystructurePipeline.handler_sources import SourceHandler
        with pytest.raises(ValueError):
            SourceHandler(self._make_table(), ["ngc9999"])


# ---------------------------------------------------------------------------
# utils.fits_utils
# ---------------------------------------------------------------------------

class TestFitsUtils:

    def test_get_beam_arcsec_missing_file(self):
        from pystructurePipeline.utils_fits import get_beam_arcsec
        with pytest.raises(FileNotFoundError):
            get_beam_arcsec("/nonexistent/file.fits")

    def test_read_fits_cube_missing_file(self):
        from pystructurePipeline.utils_fits import read_fits_cube
        with pytest.raises(FileNotFoundError):
            read_fits_cube("/nonexistent/file.fits")

    def test_hex_grid_basic(self):
        from pystructurePipeline.utils_fits import hex_grid
        x, y = hex_grid(0.0, 0.0, 0.01, radec=False, r_limit=0.05)
        assert len(x) > 0

    def test_deproject_shape(self):
        import numpy as np
        from pystructurePipeline.utils_fits import deproject
        ra  = np.linspace(202.0, 203.0, 10)
        dec = np.linspace(47.0, 48.0, 10)
        r, t = deproject(ra, dec, [173.0, 22.0, 202.47, 47.20], vector=True)
        assert r.shape == ra.shape

    def test_gaussian_PSF_2D_shape(self):
        import numpy as np
        from pystructurePipeline.utils_fits import gaussian_PSF_2D
        psf = gaussian_PSF_2D(11, [0., 1., 3., 3., 0., 0., 0.], center=True, normalize=True)
        assert psf.shape == (11, 11)
        assert abs(np.sum(psf) - 1.0) < 1e-6

    def test_deconvolve_gauss_basic(self):
        from pystructurePipeline.utils_fits import deconvolve_gauss
        maj, minn, pa, info = deconvolve_gauss(30.0, 20.0, 30.0, 0.0, 20.0, 0.0)
        assert info[0]   # worked
        assert maj > 0


# ---------------------------------------------------------------------------
# stage_fits — mask cube output
# ---------------------------------------------------------------------------

class TestStageFits:

    def _make_mask_table_and_header(self):
        """
        Build a small synthetic table with a SPEC_MASK column plus a matching
        2-D overlay header, small enough to regrid quickly in tests.
        """
        import numpy as np
        import astropy.units as au
        from astropy.table import Table, Column
        from astropy.io import fits

        n_pts, n_chan = 9, 4
        ra  = np.array([10.0, 10.001, 10.002] * 3)
        dec = np.array([20.0] * 3 + [20.001] * 3 + [20.002] * 3)
        mask = np.zeros((n_pts, n_chan))
        mask[:, 1:3] = 1  # channels 1-2 "in mask" for every point

        t = Table()
        t["ra_deg"]    = Column(ra,  unit=au.deg)
        t["dec_deg"]   = Column(dec, unit=au.deg)
        t["SPEC_MASK"] = Column(mask)
        t.meta["SPEC_VCHAN0"] = 100.0 * au.km / au.s
        t.meta["SPEC_DELTAV"] = 10.0  * au.km / au.s
        t.meta["SPEC_CRPIX"]  = 1

        hdr = fits.Header()
        hdr["NAXIS"]  = 2
        hdr["NAXIS1"] = 5
        hdr["NAXIS2"] = 5
        hdr["CTYPE1"] = "RA---TAN"
        hdr["CTYPE2"] = "DEC--TAN"
        hdr["CRVAL1"] = 10.001
        hdr["CRVAL2"] = 20.001
        hdr["CRPIX1"] = 3
        hdr["CRPIX2"] = 3
        hdr["CDELT1"] = -0.001
        hdr["CDELT2"] = 0.001
        hdr["BMAJ"]   = 0.001
        hdr["BMIN"]   = 0.001

        ov_slice = np.ones((5, 5))
        return ra, dec, hdr, ov_slice, t

    def test_save_to_fits_cube_writes_3d_cube(self, tmp_path):
        from pystructurePipeline.stage_fits import save_to_fits_cube
        from astropy.io import fits

        ra, dec, hdr, ov_slice, t = self._make_mask_table_and_header()
        save_to_fits_cube(ra, dec, hdr, ov_slice, "SPEC_MASK", "mask",
                          "testsrc", t, str(tmp_path), target_res=3.6)

        out_path = tmp_path / "testsrc_mask.fits"
        assert out_path.exists()

        data, out_hdr = fits.getdata(str(out_path), header=True)
        assert data.shape == (4, 5, 5)
        assert out_hdr["NAXIS3"] == 4
        assert out_hdr["CRVAL3"] == 100.0
        assert out_hdr["CDELT3"] == 10.0

    def test_save_to_fits_cube_output_is_binary(self, tmp_path):
        """
        Resampling onto a coarser pixel grid (when target_res is larger than
        the native beam) uses bilinear-style interpolation internally, which
        can introduce fractional values between 0 and 1. The final output
        must always be re-thresholded back to strictly 0/1 (NaN allowed
        outside the footprint), regardless of whether that resampling step
        ran.
        """
        from pystructurePipeline.stage_fits import save_to_fits_cube
        from astropy.io import fits
        import numpy as np

        ra, dec, hdr, ov_slice, t = self._make_mask_table_and_header()
        # hdr's native beam is 0.001 deg = 3.6 arcsec; request a much coarser
        # target_res so resample_hdr's reprojection branch is exercised.
        save_to_fits_cube(ra, dec, hdr, ov_slice, "SPEC_MASK", "mask",
                          "testsrc", t, str(tmp_path), target_res=20.0)

        data = fits.getdata(str(tmp_path / "testsrc_mask.fits"))
        finite = data[np.isfinite(data)]
        assert len(finite) > 0
        assert set(np.unique(finite)).issubset({0.0, 1.0})

    def test_save_to_fits_cube_skips_missing_column(self, tmp_path):
        from pystructurePipeline.stage_fits import save_to_fits_cube

        ra, dec, hdr, ov_slice, t = self._make_mask_table_and_header()
        # No "SPEC_MASK_HCN10" column exists in this table
        save_to_fits_cube(ra, dec, hdr, ov_slice, "SPEC_MASK_HCN10", "mask_hcn10",
                          "testsrc", t, str(tmp_path), target_res=3.6)
        assert not (tmp_path / "testsrc_mask_hcn10.fits").exists()

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
        from pystructurePipeline.stage_fits import construct_mask_ppv

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
        from pystructurePipeline.stage_fits import construct_mask_ppv
        from pystructurePipeline.stage_products import construct_mask

        cube, vaxis = self._make_synthetic_ppv_cube(ny=3, nx=3)
        n_chan, ny, nx = cube.shape

        # Build the equivalent hex-grid table: one row per spatial pixel
        spec = np.moveaxis(cube, 0, -1).reshape(ny * nx, n_chan)
        t = Table()
        t["SPEC_LINE"] = Column(spec)
        t.meta["SPEC_VCHAN0"] = 0.0 * au.km / au.s
        t.meta["SPEC_DELTAV"] = 1.0 * au.km / au.s
        t.meta["SPEC_CRPIX"]  = 1

        mask_ppv = construct_mask_ppv(cube, SN_processing=[2, 4])
        mask_hex, _, _ = construct_mask("LINE", t, SN_processing=[2, 4])

        mask_hex_reshaped = np.moveaxis(
            mask_hex.value.reshape(ny, nx, n_chan), -1, 0
        )
        assert np.array_equal(mask_ppv, mask_hex_reshaped)

    def test_apply_strict_mask_ppv_removes_small_components(self):
        from pystructurePipeline.stage_fits import apply_strict_mask_ppv

        mask = np.zeros((1, 10, 10), dtype=int)
        mask[0, 5, 5] = 1            # isolated single pixel: too small, removed
        mask[0, 0:3, 0:3] = 1        # 3x3 block = 9 pixels: kept

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
        from pystructurePipeline.stage_fits import get_mom_maps_ppv
        from pystructurePipeline.utils_table import get_mom_maps

        cube, vaxis_arr = self._make_synthetic_ppv_cube(ny=2, nx=2)
        n_chan, ny, nx = cube.shape

        cube_q  = cube * au.K
        vaxis_q = vaxis_arr * au.km / au.s
        mask    = (cube > 3).astype(int)
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
                ppv_maps[key].value.ravel(), flat_maps[key].value,
                equal_nan=True,
            )

    def test_convolve_cube_to_target_skips_when_already_at_resolution(self):
        from pystructurePipeline.stage_fits import convolve_cube_to_target
        from astropy.io import fits

        cube = np.ones((5, 6, 6))
        hdr = fits.Header()
        hdr["BMAJ"] = 30.0 / 3600.0  # already coarser than the 27" target
        hdr["BMIN"] = 30.0 / 3600.0

        out_data, out_hdr = convolve_cube_to_target(cube, hdr, target_res_as=27.0)
        assert np.array_equal(out_data, cube)

    def test_get_convolved_ppv_cube_uses_cached_file(self, tmp_path):
        from pystructurePipeline.stage_fits import get_convolved_ppv_cube
        from astropy.io import fits

        cube = np.arange(2 * 3 * 3, dtype=float).reshape(2, 3, 3)
        hdr = fits.Header()
        hdr["NAXIS"] = 3
        hdr["NAXIS1"], hdr["NAXIS2"], hdr["NAXIS3"] = 3, 3, 2
        cached_path = tmp_path / "testsrc_co_27.0as.fits"
        fits.writeto(str(cached_path), cube, hdr)

        data, _ = get_convolved_ppv_cube(
            "testsrc", "co", "/nonexistent_dir", ".fits",
            27.0, hdr, str(tmp_path),
        )
        assert np.array_equal(data, cube)

    def test_get_convolved_ppv_cube_raises_if_nothing_available(self, tmp_path):
        from pystructurePipeline.stage_fits import get_convolved_ppv_cube
        from astropy.io import fits

        hdr = fits.Header()
        hdr["NAXIS"] = 3
        hdr["NAXIS1"], hdr["NAXIS2"], hdr["NAXIS3"] = 3, 3, 2

        with pytest.raises(FileNotFoundError):
            get_convolved_ppv_cube(
                "testsrc", "co", str(tmp_path), ".fits",
                27.0, hdr, str(tmp_path),
            )


# ---------------------------------------------------------------------------
# utils.table_utils
# ---------------------------------------------------------------------------

class TestTableUtils:

    def test_load_missing_file(self):
        from pystructurePipeline.utils_table import load_pystructure
        with pytest.raises(FileNotFoundError):
            load_pystructure("/nonexistent/file.ecsv")

    def test_find_latest_missing(self, tmp_path):
        from pystructurePipeline.utils_table import find_latest_pystructure
        with pytest.raises(FileNotFoundError):
            find_latest_pystructure(str(tmp_path), "ngc5194")

    def test_shuffle_roundtrip(self):
        import numpy as np
        from pystructurePipeline.utils_table import shuffle
        vaxis    = np.arange(-100, 101, 1.0)
        spec     = np.exp(-0.5 * (vaxis / 20.0) ** 2)
        shuffled = shuffle(spec, vaxis, zero=0.0, new_vaxis=vaxis)
        # Should be identical (same axis, zero shift)
        assert np.allclose(shuffled, spec, equal_nan=True)

    def test_get_mom_maps_runs(self):
        import numpy as np
        from astropy import units as u
        from pystructurePipeline.utils_table import get_mom_maps

        n_pts, n_chan = 5, 50
        vaxis = np.linspace(-100, 100, n_chan) * u.km / u.s
        spec  = np.zeros((n_pts, n_chan)) * u.K
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
        from pystructurePipeline.init_workdir import init_workdir
        init_workdir(str(tmp_path))
        assert (tmp_path / "config.txt").exists()
        assert (tmp_path / "keys" / "target_definitions.txt").exists()
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

    def test_missing_conf_exits(self):
        from pystructurePipeline.cli import main
        with pytest.raises((SystemExit, FileNotFoundError)):
            main(["--conf", "/nonexistent/config.txt"])

    def test_no_args_exits(self):
        from pystructurePipeline.cli import main
        with pytest.raises(SystemExit):
            main([])

    def test_invalid_stage_exits(self):
        from pystructurePipeline.cli import main
        with pytest.raises(SystemExit):
            main(["--conf", "config.txt", "--stages", "invalid_stage"])


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
        assert "[pyStructure] [Regrid]    [INFO]     hello world" in captured.out

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
        assert "[pyStructure] [FITS]      [ERROR]    file not found" in content
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
        assert "[Sampling]  [INFO]     a" in content
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
        # Re-use the minimal config from TestKeyHandler
        kh = TestKeyHandler()
        conf_path = kh._write_minimal_config(tmp_path)

        from pystructurePipeline.handler_pipeline import PipelineHandler
        log_path = tmp_path / "run.log"
        handler = PipelineHandler(conf_path=str(conf_path), verbose=False,
                                  log_file=str(log_path))
        assert log_path.exists()
        content = log_path.read_text()
        assert "[Loading]" in content
        assert "Loading configuration..." in content
