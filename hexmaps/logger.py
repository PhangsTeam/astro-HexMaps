"""
logger.py — centralized logging for the HexMaps pipeline.

All pipeline modules log through a single shared PipelineLogger instance,
obtained via get_logger(stage). This gives every message a consistent format:

    YYYY-MM-DD HH:MM:SS [<Stage>] [<LEVEL>] <message>

The closing bracket of [<Stage>] and [<LEVEL>] always follows directly after
the stage/level name (no padding inside the brackets); instead, the space
*after* each bracket is padded out to a fixed column width (based on the
length of the longest known stage name / level name) so that the message
column lines up across all log lines, regardless of which stage or level
produced them:

e.g.
    YYYY-MM-DD HH:MM:SS [Regrid]   [INFO]    Map SPIRE250 sampled successfully.
    YYYY-MM-DD HH:MM:SS [Products] [ERROR]   12CO21 spectrum is all zeros; skipping.
    YYYY-MM-DD HH:MM:SS [Loading]  [INFO]    Loading key files...
    YYYY-MM-DD HH:MM:SS [FITS]     [WARNING] pixels_per_beam < 4; expect artefacts.

In addition to printing, every message is stored as a structured record
(timestamp, stage, level, message). The full log can optionally be written
to a file at the end of a run via PipelineHandler(..., log_file=...) or by
calling logger.save(path) directly.

Usage
-----
Module-level logger bound to a fixed stage name::

    from hexmaps.logger import get_logger
    LOG = get_logger("Regrid")

    LOG.info("Map SPIRE250 sampled successfully.")
    LOG.warning("pixels_per_beam < 4; expect artefacts.")
    LOG.error(f"Cube {name} not found: {path}")

Configuring verbosity / file output (done once, typically by PipelineHandler)::

    from hexmaps.logger import logger
    logger.configure(verbose=True, log_file="hexmaps_run.log")
    ...
    logger.save("hexmaps_run.log")   # also done automatically if log_file is set
"""

import os
import datetime

# Column widths for the [<Stage>] and [<LEVEL>] fields, INCLUDING the
# brackets and a 2-space separator. The closing bracket always follows
# directly after the stage/level name; any extra space needed for alignment
# is added *after* the bracket. Update these if a longer stage name or level
# name is introduced.
_STAGE_COL_WIDTH = len("Products") + 2 + 1  # "[Products]" + 2 spaces
_LEVEL_COL_WIDTH = len("WARNING") + 2 + 1  # "[WARNING]"  + 2 spaces


class PipelineLogger:
    """
    Shared logger for all HexMaps pipeline modules.

    A single instance of this class (``logger``, defined at the bottom of this
    module) is imported by every pipeline module.  Each module obtains a
    stage-bound view of it via get_logger(stage_name), which is a thin wrapper
    that prefixes every call with the module's stage name.

    Attributes
    ----------
    verbose  : bool — if True (default), print every message to stdout
    log_file : str or None — if set, every message is appended to this file
               as it is logged (in addition to being printed)
    records  : list of dict — all messages logged so far, each with keys
               'time', 'stage', 'level', 'message'
    """

    def __init__(self):
        self.verbose = True
        self.log_file = None
        self.records = []

    # ------------------------------------------------------------------
    # Configuration
    # ------------------------------------------------------------------

    def configure(self, verbose: bool = True, log_file: str = None, reset: bool = True):
        """
        Configure global logging behaviour.

        Parameters
        ----------
        verbose  : bool — print messages to stdout (default True)
        log_file : str or None — if given, append every message to this file
                   as it is logged. The file is created (and any existing
                   content overwritten) when configure() is called.
        reset    : bool — if True (default), clear any previously recorded
                   messages from a prior run.
        """
        self.verbose = verbose
        self.log_file = log_file
        if reset:
            self.records = []
        if self.log_file:
            os.makedirs(
                os.path.dirname(os.path.abspath(self.log_file)) or ".", exist_ok=True
            )
            with open(self.log_file, "w") as f:
                f.write(
                    f"HexMaps log started at "
                    f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                )

    # ------------------------------------------------------------------
    # Core logging
    # ------------------------------------------------------------------

    def log(self, stage: str, level: str, message: str):
        """
        Record and (optionally) print a log message.

        Parameters
        ----------
        stage   : str — pipeline stage name, e.g. "Regrid", "Spectra"
        level   : str — "INFO", "WARNING", or "ERROR"
        message : str — the message text
        """
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.records.append(
            {
                "time": timestamp,
                "stage": stage,
                "level": level,
                "message": message,
            }
        )

        stage_field = f"[{stage}]".ljust(_STAGE_COL_WIDTH)
        level_field = f"[{level}]".ljust(_LEVEL_COL_WIDTH)
        time_field = f"{timestamp} "
        formatted = f"{time_field}{stage_field}{level_field}{message}"

        if self.verbose:
            print(formatted)

        if self.log_file:
            with open(self.log_file, "a") as f:
                f.write(f"{timestamp}  {formatted}\n")

    def info(self, stage: str, message: str):
        """Log an INFO-level message for *stage*."""
        self.log(stage, "INFO", message)

    def warning(self, stage: str, message: str):
        """Log a WARNING-level message for *stage*."""
        self.log(stage, "WARNING", message)

    def error(self, stage: str, message: str):
        """Log an ERROR-level message for *stage*."""
        self.log(stage, "ERROR", message)

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def as_text(self) -> str:
        """
        Return the full log history as a plain-text string.

        Same format as ``save()`` — one line per record:
            YYYY-MM-DD HH:MM:SS [Stage]   [LEVEL]   <message>

        Returns
        -------
        str — the complete log output, with real newlines between lines.
        """
        lines = []
        for r in self.records:
            stage_field = f"[{r['stage']}]".ljust(_STAGE_COL_WIDTH)
            level_field = f"[{r['level']}]".ljust(_LEVEL_COL_WIDTH)
            lines.append(f"{r['time']} {stage_field}{level_field}{r['message']}")
        return "\n".join(lines)

    def save(self, path: str):
        """
        Write the full log history to *path*, overwriting any existing file.

        Each line has the format::

            YYYY-MM-DD HH:MM:SS [Stage] [LEVEL] <message>

        Parameters
        ----------
        path : str — output file path; parent directories are created if needed
        """
        path = str(path)
        os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
        with open(path, "w") as f:
            for r in self.records:
                stage_field = f"[{r['stage']}]".ljust(_STAGE_COL_WIDTH)
                level_field = f"[{r['level']}]".ljust(_LEVEL_COL_WIDTH)
                time_field = f"{r['time']} "
                f.write(f"{time_field}{stage_field}{level_field}{r['message']}\n")

    def get_records(self, stage: str = None, level: str = None) -> list:
        """
        Return logged records, optionally filtered by stage and/or level.

        Parameters
        ----------
        stage : str, optional — only return records for this stage
        level : str, optional — only return records at this level
                ("INFO", "WARNING", or "ERROR")

        Returns
        -------
        list of dict — matching records
        """
        records = self.records
        if stage is not None:
            records = [r for r in records if r["stage"] == stage]
        if level is not None:
            records = [r for r in records if r["level"] == level]
        return records


class StageLogger:
    """
    Stage-bound view of a PipelineLogger.

    Returned by get_logger(stage); every call is forwarded to the shared
    PipelineLogger with the stage name pre-filled, so individual modules
    never need to repeat their stage name.

    Example
    -------
    >>> LOG = get_logger("Regrid")
    >>> LOG.info("Map SPIRE250 sampled successfully.")
    YYYY-MM-DD HH:MM:SS [Regrid] [INFO] Map SPIRE250 sampled successfully.
    """

    def __init__(self, parent: PipelineLogger, stage: str):
        self._parent = parent
        self._stage = stage

    def info(self, message: str):
        """Log an INFO-level message for this stage."""
        self._parent.info(self._stage, message)

    def warning(self, message: str):
        """Log a WARNING-level message for this stage."""
        self._parent.warning(self._stage, message)

    def error(self, message: str):
        """Log an ERROR-level message for this stage."""
        self._parent.error(self._stage, message)


# ---------------------------------------------------------------------------
# Module-level singleton
# ---------------------------------------------------------------------------

#: The single shared logger instance used by the whole pipeline.
logger = PipelineLogger()


def get_logger(stage: str) -> StageLogger:
    """
    Return a StageLogger bound to *stage*, backed by the shared PipelineLogger.

    Parameters
    ----------
    stage : str — pipeline stage name, e.g. "Regrid", "Products", "FITS"

    Returns
    -------
    StageLogger — call .info(msg) / .warning(msg) / .error(msg) on this object
    """
    return StageLogger(logger, stage)
