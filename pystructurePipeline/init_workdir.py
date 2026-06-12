"""
pystructure.init_workdir
========================
Copies the bundled key-file templates and a run script into a user-chosen
working directory so they can get started without hunting for example files.

Called via the CLI:
    pystructure --init [--workdir ./my_project]

Or from Python:
    from pystructurePipeline import init_workdir
    init_workdir("./my_project")
"""

import shutil
import os
from pathlib import Path


# The templates are bundled inside the installed package
_TEMPLATES_DIR = Path(__file__).parent / "templates"


def init_workdir(workdir: str = ".", overwrite: bool = False) -> None:
    """
    Initialise a PyStructure working directory.

    Copies the following into *workdir*:
      keys/master_key.txt
      keys/target_definitions.txt
      keys/imaging_key.txt
      keys/config_key.txt
      run_pystructure.py          ← ready-to-edit run script

    Parameters
    ----------
    workdir   : str or Path — destination directory (created if absent)
    overwrite : bool — if False, raise if any existing file already exists
    """
    workdir = Path(workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    keys_dst = workdir / "keys"
    keys_dst.mkdir(exist_ok=True)

    keys_src = _TEMPLATES_DIR / "keys"

    copied = []

    # --- Key files ---
    for key_file in keys_src.iterdir():
        dst = keys_dst / key_file.name
        if dst.exists() and not overwrite:
            raise FileExistsError(
                f"{dst} already exists. Use overwrite=True to replace it."
            )
        shutil.copy2(key_file, dst)
        copied.append(str(dst.relative_to(workdir)))

    # --- Run script ---
    run_script_src = _TEMPLATES_DIR / "run_pystructure.py"
    run_script_dst = workdir / "run_pystructure.py"
    if run_script_dst.exists() and not overwrite:
        raise FileExistsError(
            f"{run_script_dst} already exists. Use overwrite=True to replace it."
        )
    shutil.copy2(run_script_src, run_script_dst)
    copied.append("run_pystructure.py")

    print(f"[INFO]     PyStructure working directory initialised at: {workdir}")
    print(f"[INFO]     Files created:")
    for f in copied:
        print(f"[INFO]       {f}")
    print(f"[INFO]     Next steps:")
    print(f"[INFO]       1. Edit keys/master_key.txt  — set your data_dir and out_dir")
    print(f"[INFO]       2. Edit keys/target_definitions.txt  — add your sources")
    print(f"[INFO]       3. Edit keys/imaging_key.txt  — list your bands and cubes")
    print(f"[INFO]       4. Edit keys/config_key.txt  — adjust resolution and masking")
    print(f"[INFO]       5. Run:  python run_pystructure.py  (or:  pystructure --key_dir keys/)")
