"""ReAnnota: Genome annotation enhancement pipeline."""

__version__ = "0.1.0"
__author__ = "Aggelos Pazis, Efthymios Parisis"

# Re-export main CLI and app for convenience
from reannota.cli import app, main

__all__ = ["app", "main", "__version__"]
