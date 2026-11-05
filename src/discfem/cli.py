"""Línea de comandos: ``python -m discfem run | mesh | study <config.yaml>``.

En la Fase 1 solo está implementado ``mesh``; ``run`` llega en la Fase 5 y ``study`` en la 6.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import yaml

from .geometry import BOUNDARY_GROUPS, RectangleGeometry
from .mesh import MeshOptions, build_mesh


def _load_yaml(path: Path) -> dict:
    """Lee el YAML de configuración; error claro si no existe o no es un diccionario."""
    if not path.is_file():
        raise SystemExit(f"Error: no existe el archivo de configuración '{path}'.")
    with path.open(encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise SystemExit(f"Error: '{path}' no contiene un diccionario YAML.")
    return cfg


def cmd_mesh(config_path: Path) -> int:
    """Genera la malla del caso, imprime un resumen y guarda ``mesh.png`` y una copia de la
    configuración en ``output.dir``."""
    cfg = _load_yaml(config_path)
    try:
        geo = RectangleGeometry.from_dict(cfg.get("geometry", {}))
        opts = MeshOptions.from_dict(cfg.get("mesh", {}))
    except ValueError as exc:
        raise SystemExit(f"Error en '{config_path}': {exc}") from None
    mesh = build_mesh(geo, opts)

    out = Path(cfg.get("output", {}).get("dir", f"results/{cfg.get('case', 'caso')}"))
    out.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(config_path, out / "config.yaml")
    plot_path = out / "mesh.png"
    from .postprocess import plot_mesh

    plot_mesh(mesh, plot_path, title=f"{cfg.get('case', '')} [{opts.backend}]")

    A = mesh.areas()
    print(f"Malla '{opts.backend}' (refine_factor = {opts.refine_factor:g})")
    print(f"  nodos: {mesh.n_nodes}   triángulos: {mesh.n_elements}")
    print(f"  área mín/máx: {A.min():.3e} / {A.max():.3e} m²")
    print(f"  Σ A·r̄ = {np.sum(A * mesh.r_centroid()):.12e} m³ "
          f"(exacto e(r_o²−r_i²)/2 = "
          f"{geo.half_thickness * (geo.r_outer**2 - geo.r_inner**2) / 2:.12e})")
    for g in BOUNDARY_GROUPS:
        print(f"  {g:10s}: {len(mesh.edges[g])} aristas")
    print(f"  figura: {plot_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Punto de entrada de la CLI."""
    # La consola de Windows usa cp1252 por defecto y no puede imprimir Σ, r̄, ²...
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(prog="python -m discfem", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, help_ in [
        ("run", "corrida completa (Fase 5)"),
        ("mesh", "solo malla + figura"),
        ("study", "estudio de convergencia (Fase 6)"),
    ]:
        p = sub.add_parser(name, help=help_)
        p.add_argument("config", type=Path, help="archivo YAML de configuración")
    args = parser.parse_args(argv)

    if args.command == "mesh":
        return cmd_mesh(args.config)
    print(f"El comando '{args.command}' aún no está implementado (ver PLAN.md).", file=sys.stderr)
    return 1
