"""
hexmaps.cli — entry point for the `hexmaps` console script.

Installed by pip via pyproject.toml [project.scripts]:
    hexmaps = "hexmaps.cli:main"

Usage
-----
Initialise a working directory (copies key templates + run script):
    hexmaps --init
    hexmaps --init --workdir ./my_project

Run the pipeline (from inside a working directory):
    hexmaps --conf config.txt
    hexmaps --conf config.txt --stages regrid products
    hexmaps --conf config.txt --targets ngc5194
    hexmaps --conf config.txt --log_file hexmaps_run.log
"""

import argparse
import sys

ALL_STAGES = ["regrid", "products", "fits"]
DATABASE_STAGES = ["regrid", "products"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="HexMaps: homogenize and analyze multi-wavelength astronomical datasets.",
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
        "--conf",
        default=None,
        help="Path to your config.txt configuration file.",
    )
    parser.add_argument(
        "--stages",
        nargs="+",
        choices=ALL_STAGES + ["all"],
        default=None,
        help=(
            f"Pipeline stage(s) to run: {', '.join(ALL_STAGES)}, or 'all' to run every stage. "
            f"Default: {' '.join(DATABASE_STAGES)} (fits is optional and not run by default)."
        ),
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        default=None,
        help="Source name(s) to process. Default: all sources in config.txt.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress informational output.",
    )
    parser.add_argument(
        "--log_file",
        default=None,
        help="Write all log messages to this file in addition to stdout.",
    )

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # ------------------------------------------------------------------
    # --init mode: scaffold a working directory and exit
    # ------------------------------------------------------------------
    if args.init:
        from hexmaps.init_workdir import init_workdir

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
    if args.conf is None:
        print(
            "[ERROR]    --conf is required when not using --init.\n"
            "           Example: hexmaps --conf config.txt\n"
            "           To set up a new project: hexmaps --init"
        )
        sys.exit(1)

    try:
        from hexmaps import PipelineHandler
    except ImportError as exc:
        print(f"[ERROR]    Could not import hexmaps: {exc}")
        sys.exit(1)

    raw_stages = args.stages if args.stages else DATABASE_STAGES
    # "all" is a convenience alias for ALL_STAGES; run_stages also handles
    # it, but expanding here gives a cleaner log line.
    stages = ALL_STAGES if raw_stages == ["all"] else raw_stages
    handler = PipelineHandler(
        conf_path=args.conf, verbose=not args.quiet, log_file=args.log_file
    )

    try:
        handler.run_stages(stages=stages, targets=args.targets)
    except Exception as exc:
        print(f"[ERROR]    Pipeline terminated unexpectedly: {exc}")
        sys.exit(1)

    if not all(handler.run_success.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
