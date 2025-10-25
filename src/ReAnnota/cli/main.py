"""Provide a command-line interface (CLI) for ReAnnota functionality."""

import logging
from enum import Enum, unique
from pathlib import Path
from typing import Optional

import typer

logger = logging.getLogger("ReAnnota")


@unique
class LogLevel(str, Enum):
    """Define the choices for the log level option."""

    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"


# Create the main Typer app
app = typer.Typer(
    help="Cyanobacteria genome annotation enhancement pipeline",
    context_settings={"help_option_names": ["-h", "--help"]},
    add_completion=False,  # Disable shell completion options
)


def version_callback(is_set: bool) -> None:
    """
    Print the tool version if desired.

    Args:
        is_set: Whether the version was requested as a command line option.

    Raises:
        typer.Exit: With default code 0 to signal normal program end.
    """
    if is_set:
        from ReAnnota import __version__

        typer.echo(f"ReAnnota v{__version__}")
        raise typer.Exit(code=0)


def setup_logger(log_file: Path, log_level: str = "INFO") -> None:
    """
    Configure logging system.

    Args:
        log_file: Path to the log file.
        log_level: The desired log level.
    """
    logging.basicConfig(
        filename=log_file,
        filemode="w",  # overwrite each run
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=getattr(logging, log_level),
    )
    # Log to console
    console = logging.StreamHandler()
    console.setLevel(getattr(logging, log_level))
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)


@app.callback(invoke_without_command=True)
def initialize(
    context: typer.Context,
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        "-v",
        callback=version_callback,
        is_eager=True,
        help="Print the current tool version and exit.",
    ),
) -> None:
    """Initialize the ReAnnota CLI."""
    # If no command was provided, show help
    if context.invoked_subcommand is None:
        typer.echo(context.get_help())
        raise typer.Exit()

