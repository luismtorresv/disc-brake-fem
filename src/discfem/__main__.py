"""Punto de entrada ``python -m discfem``."""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
