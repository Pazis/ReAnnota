"""Provide the cyano-annotation command-line interface (CLI)."""

from .main import app

# Import commands to register them with the app
from . import annotate  # noqa: F401


def main() -> None:
    """Entry point for the CLI."""
    app()


__all__ = ["app", "main"]

