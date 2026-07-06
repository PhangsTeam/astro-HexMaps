"""
hexmaps.download_example
========================
Download the NGC 5194 example dataset into a working directory so users
can run the pipeline immediately after ``hexmaps --init``.

CLI usage::

    hexmaps --download-example
    hexmaps --download-example --workdir ./my_project
    hexmaps --download-example --workdir ./my_project --force

Python usage::

    from hexmaps.download_example import download_example_data
    download_example_data("./my_project")

The files are fetched from the public GitHub repository of HexMaps at
https://github.com/lukas-neumann-astro/astro-HexMaps and placed in the
``data/`` sub-directory of *workdir*, which is created if it does not
exist.

Only the raw input files required to run the example pipeline are
downloaded (the pre-convolved ``*_27.0as.fits`` intermediate files are
not included — they are produced by the pipeline itself).
"""

import os
import sys
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Files to download
# ---------------------------------------------------------------------------

#: Base URL of the raw data files on GitHub.
_BASE_URL = (
    "https://raw.githubusercontent.com/"
    "lukas-neumann-astro/astro-HexMaps/main/data/"
)

#: Files that are actually needed to run the example pipeline.
#: Keys are destination filenames; values are source filenames on GitHub
#: (same in this case, but kept explicit for future flexibility).
_EXAMPLE_FILES = {
    "ngc5194_12co21.fits":              "ngc5194_12co21.fits",
    "ngc5194_12co10.fits":              "ngc5194_12co10.fits",
    "ngc5194_spire250_gauss21.fits":    "ngc5194_spire250_gauss21.fits",
    "ngc5194_spire250_gauss21_unc.fits":"ngc5194_spire250_gauss21_unc.fits",
    "ngc5194_12co21_ii.fits":           "ngc5194_12co21_ii.fits",
    "ngc5194_12co21_ii_uc.fits":        "ngc5194_12co21_ii_uc.fits",
}

# Total approximate download size (shown to user before starting)
_TOTAL_MB = 46


# ---------------------------------------------------------------------------
# Progress hook
# ---------------------------------------------------------------------------

def _progress_hook(filename: str):
    """Return a reporthook function that prints a simple progress bar."""

    def _hook(block_count, block_size, total_size):
        if total_size <= 0:
            downloaded = block_count * block_size
            sys.stdout.write(f"\r  Downloading {filename} … {downloaded // 1024:>6} KB")
        else:
            downloaded = min(block_count * block_size, total_size)
            pct = int(100 * downloaded / total_size)
            bar = "█" * (pct // 5) + "░" * (20 - pct // 5)
            sys.stdout.write(
                f"\r  Downloading {filename} … [{bar}] {pct:3d}%"
                f"  ({downloaded // 1024:>6} / {total_size // 1024} KB)"
            )
        sys.stdout.flush()

    return _hook


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def download_example_data(workdir: str = ".", force: bool = False) -> None:
    """
    Download the NGC 5194 example dataset into *workdir*/data/.

    Parameters
    ----------
    workdir : str or Path
        Root of the HexMaps working directory (the same path you passed to
        ``hexmaps --init``).  A ``data/`` sub-directory is created inside it.
    force   : bool
        If False (default), skip files that already exist.
        If True, overwrite existing files.

    Raises
    ------
    urllib.error.URLError
        If a file cannot be fetched (no network, URL changed, etc.).
    """
    workdir = Path(workdir).resolve()
    data_dir = workdir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"[INFO]     Downloading HexMaps example data (NGC 5194) into:\n"
        f"[INFO]       {data_dir}\n"
        f"[INFO]     Approximate download size: ~{_TOTAL_MB} MB\n"
    )

    downloaded, skipped, failed = [], [], []

    for dst_name, src_name in _EXAMPLE_FILES.items():
        dst = data_dir / dst_name
        if dst.exists() and not force:
            print(f"  [skip]     {dst_name}  (already exists; use --force to overwrite)")
            skipped.append(dst_name)
            continue

        url = _BASE_URL + src_name
        try:
            urllib.request.urlretrieve(url, dst, reporthook=_progress_hook(dst_name))
            sys.stdout.write("\n")   # newline after progress bar
            downloaded.append(dst_name)
        except Exception as exc:
            sys.stdout.write("\n")
            print(f"  [ERROR]    Failed to download {dst_name}: {exc}")
            if dst.exists():
                dst.unlink()         # remove partial file
            failed.append(dst_name)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print()
    if downloaded:
        print(f"[INFO]     Downloaded {len(downloaded)} file(s):")
        for f in downloaded:
            size_kb = (data_dir / f).stat().st_size // 1024
            print(f"[INFO]       {f}  ({size_kb:,} KB)")
    if skipped:
        print(f"[INFO]     Skipped {len(skipped)} existing file(s) (use --force to overwrite).")
    if failed:
        print(
            f"[ERROR]    {len(failed)} file(s) could not be downloaded:\n"
            + "".join(f"[ERROR]      {f}\n" for f in failed)
            + "[ERROR]    Check your internet connection and try again."
        )
        sys.exit(1)

    if not failed:
        print(
            f"\n[INFO]     Example data ready.  You can now run:\n"
            f"[INFO]       hexmaps --conf {workdir / 'config.txt'}\n"
            f"[INFO]     or:\n"
            f"[INFO]       python {workdir / 'run_hexmaps.py'}"
        )
