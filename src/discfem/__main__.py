"""Punto de entrada `python -m discfem`. La CLI real se implementa en la Fase 5."""

import sys


def main() -> int:
    """Informa que la CLI aún no existe y devuelve un código de error."""
    print("CLI pendiente (Fase 5 de PLAN.md).", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
