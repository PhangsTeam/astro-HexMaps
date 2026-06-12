"""
pystructurePipeline.cli — entry point for the `pystructure` console script.

Installed by pip via pyproject.toml [project.scripts]:
    pystructure = "pystructurePipeline.cli:main"

Usage
-----
Initialise a working directory (copies key templates + run script):
    pystructure --init
    pystructure --init --workdir ./my_project

Run the pipeline (from inside a working directory):
    pystructure --key_dir keys/
    pystructure --key_dir keys/ --stages sampling regrid
    pystructure --key_dir keys/ --targets ngc5194
"""

import argparse
import sys

ALL_STAGES = ["sampling", "regrid", "spectra", "output"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="PyStructure: homogenize and analyze multi-wavelength astronomical datasets.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # --- Init mode ---
    parser.add_argument(
        "--init",
        action="store_true",
        help=(
            "Initialise a new working directory: copies key-file templates "
            "and a run script. Use --workdir to set the destination."
        ),
    )
    parser.add_argument(
        "--workdir",
        default=".",
        help="Working directory to initialise (used with --init). Default: current directory.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing files when using --init.",
    )

    # --- Run mode ---
    parser.add_argument(
        "--key_dir",
        default=None,
        help="Path to the directory containing your key files.",
    )
    parser.add_argument(
        "--stages",
        nargs="+",
        choices=ALL_STAGES,
        default=None,
        help=(
            f"Pipeline stage(s) to run: {', '.join(ALL_STAGES)}. Default: all."
        ),
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        default=None,
        help="Source name(s) to process. Default: all sources in imaging_key.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress informational output.",
    )

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # ------------------------------------------------------------------
    # --init mode: scaffold a working directory and exit
    # ------------------------------------------------------------------
    if args.init:
        from pystructurePipeline.init_workdir import init_workdir
        try:
            init_workdir(workdir=args.workdir, overwrite=args.overwrite)
        except FileExistsError as exc:
            print(f"[ERROR]    {exc}")
            print("[ERROR]    Use --overwrite to replace existing files.")
            sys.exit(1)
        return

    # ------------------------------------------------------------------
    # Run mode
    # ------------------------------------------------------------------
    if args.key_dir is None:
        print(
            "[ERROR]    --key_dir is required when not using --init.\n"
            "           Example: pystructure --key_dir keys/\n"
            "           To set up a new project: pystructure --init"
        )
        sys.exit(1)

    try:
        from pystructurePipeline import PipelineHandler
    except ImportError as exc:
        print(f"[ERROR]    Could not import pystructurePipeline: {exc}")
        sys.exit(1)

    stages = args.stages if args.stages else ALL_STAGES
    handler = PipelineHandler(key_dir=args.key_dir, verbose=not args.quiet)

    try:
        handler.run_stages(stages=stages, targets=args.targets)
    except Exception as exc:
        print(f"[ERROR]    Pipeline terminated unexpectedly: {exc}")
        sys.exit(1)

    if not all(handler.run_success.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
