"""
hexmaps.init_workdir
========================
Copies the bundled config/key templates and a run script into a user-chosen
working directory so they can get started without hunting for example files.

Called via the CLI:
    hexmaps --init [--workdir ./my_project]

Or from Python:
    from hexmaps import init_workdir
    init_workdir("./my_project")
"""

import shutil
from pathlib import Path

# Templates are bundled inside the installed package at templates/.
# This directory ships as package-data (see pyproject.toml), so it is
# present regardless of how hexmaps was installed (PyPI wheel, git URL,
# or editable local clone).
_PACKAGE_DIR = Path(__file__).parent
_TEMPLATES_DIR = _PACKAGE_DIR / "templates"


def init_workdir(workdir: str = ".", overwrite: bool = False) -> None:
    """
    Initialise a HexMaps working directory.

    Copies the following into *workdir*:
      config.txt                    ← the file you edit on every run
      keys/target_definitions.txt   ← target geometry (edit once, reuse)
      keys/hfs_lines.txt            ← hyperfine structure lines (optional, edit once)
      run_hexmaps.py                ← ready-to-edit run script

    Parameters
    ----------
    workdir   : str or Path — destination directory (created if absent)
    overwrite : bool — if False, raise if any existing file already exists
    """
    if not _TEMPLATES_DIR.exists():
        raise FileNotFoundError(
            f"Bundled templates not found at {_TEMPLATES_DIR}. "
            "This usually means hexmaps was installed without its package "
            "data (e.g. a broken build). Try reinstalling with: "
            "pip install --force-reinstall astro-hexmaps"
        )

    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    keys_dst = workdir / "keys"
    keys_dst.mkdir(exist_ok=True)
    keys_src = _TEMPLATES_DIR / "keys"

    copied = []

    # --- Unified config file ---
    conf_src = _TEMPLATES_DIR / "config.txt"
    conf_dst = workdir / "config.txt"
    if conf_dst.exists() and not overwrite:
        raise FileExistsError(
            f"{conf_dst} already exists. Use overwrite=True to replace it."
        )
    shutil.copy2(conf_src, conf_dst)
    copied.append("config.txt")

    # --- keys/ subfolder: target_definitions.txt + hfs_lines.txt ---
    for key_file in keys_src.iterdir():
        dst = keys_dst / key_file.name
        if dst.exists() and not overwrite:
            raise FileExistsError(
                f"{dst} already exists. Use overwrite=True to replace it."
            )
        shutil.copy2(key_file, dst)
        copied.append(str(dst.relative_to(workdir)))

    # --- Run script ---
    run_script_src = _TEMPLATES_DIR / "run_hexmaps.py"
    run_script_dst = workdir / "run_hexmaps.py"
    if run_script_dst.exists() and not overwrite:
        raise FileExistsError(
            f"{run_script_dst} already exists. Use overwrite=True to replace it."
        )
    shutil.copy2(run_script_src, run_script_dst)
    copied.append("run_hexmaps.py")

    print(f"[INFO]     HexMaps working directory initialised at: {workdir}")
    print(f"[INFO]     Files created:")
    for f in copied:
        print(f"[INFO]       {f}")
    print(f"[INFO]     Next steps:")
    print(
        f"[INFO]       1. (optional) Download the NGC 5194 example dataset (~46 MB):\n"
        f"[INFO]              hexmaps --download-example --workdir {workdir}"
    )
    print(
        f"[INFO]       2. Edit config.txt  — paths, targets, maps/cubes, resolution, masking"
    )
    print(
        f"[INFO]       3. Edit keys/target_definitions.txt  — add your targets"
    )
    print(
        f"[INFO]       4. (optional) Edit keys/hfs_lines.txt  — hyperfine structure lines"
    )
    print(
        f"[INFO]       5. Run:  hexmaps --conf config.txt  (or:  python run_hexmaps.py)"
    )
