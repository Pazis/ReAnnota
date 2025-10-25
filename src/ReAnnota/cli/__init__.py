"""Provide the ReAnnota command-line interface (CLI)."""

# Import commands to register them with the app
from . import annotate  # noqa: F401
from .main import app


def main() -> None:
    """Entry point for the CLI."""
    app()


__all__ = ["app", "main"]

