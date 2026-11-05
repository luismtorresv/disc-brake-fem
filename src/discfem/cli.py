"""Línea de comandos: ``python -m discfem run | mesh | study <config.yaml>``.

- ``run``: corrida completa; ``--set sección.clave=valor`` (repetible) aplica overrides.
- ``mesh``: solo malla + figura.
- ``study``: estudios de convergencia (Fase 6, aún no implementado).
"""

from __future__ import annotations

import argparse
import shutil
import sys
import warnings
from pathlib import Path

import numpy as np
import yaml

from .config import set_dotted
from .geometry import BOUNDARY_GROUPS, RectangleGeometry
from .materials import PropertyRangeWarning, UnverifiedMaterialWarning
from .mesh import MeshOptions, build_mesh


def _parse_set(items: list[str]) -> dict:
    """``["mesh.refine_factor=2", "time.theta=0.5"]`` → ``{"mesh.refine_factor": 2, ...}``
    (el valor se interpreta como YAML: números, listas, true/false…)."""
    out = {}
    for item in items or []:
        if "=" not in item:
            raise SystemExit(f"Error: --set espera 'clave=valor'; se recibió '{item}'.")
        key, value = item.split("=", 1)
        out[key.strip()] = yaml.safe_load(value)
    return out


def _format_warning(message, category, filename, lineno, line=None) -> str:
    return f"ADVERTENCIA ({category.__name__}): {message}\n"


def cmd_mesh(config_path: Path, overrides: dict) -> int:
    """Genera la malla del caso, imprime un resumen y guarda ``mesh.png`` y una copia de la
    configuración en ``output.dir``. Solo valida ``geometry``, ``mesh`` y ``output``: basta un
    YAML con esas secciones."""
    if not config_path.is_file():
        raise ValueError(f"No existe el archivo de configuración '{config_path}'.")
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError(f"'{config_path}' no contiene un diccionario YAML.")
    for key, value in overrides.items():
        set_dotted(raw, key, value)
    geo = RectangleGeometry.from_dict(raw.get("geometry") or {})
    opts = MeshOptions.from_dict(raw.get("mesh") or {})
    case = raw.get("case", "caso")
    mesh = build_mesh(geo, opts)
    out = Path((raw.get("output") or {}).get("dir", f"results/{case}"))
    out.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(config_path, out / "config.yaml")
    plot_path = out / "mesh.png"
    from .postprocess import plot_mesh

    plot_mesh(mesh, plot_path, title=f"{case} [{opts.backend}]")

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


def cmd_run(config_path: Path, overrides: dict) -> int:
    """Corrida completa: escribe las salidas y muestra un resumen."""
    from .solver import run

    res = run(config_path, overrides)
    s = res.summary
    print(f"Caso '{s['case']}': {s['mesh']['n_nodes']} nodos, {s['mesh']['n_elements']} "
          f"triángulos, {s['time']['n_steps']} pasos (θ = {s['time']['theta']:g}, "
          f"Rannacher = {s['time']['rannacher_steps']}, Δt = {s['time']['dt']:.4g} s)")
    if s["load"]["gamma"] is not None:
        print(f"  γ = {s['load']['gamma']:.4f}   t_f = {s['time']['t_f']:.4g} s")
    tm = s["T_max"]
    print(f"  T_max = {tm['value']:.3f} en t = {tm['time']:.3f} s (r = {tm['r']:.4f}, "
          f"z = {tm['z']:.4f})")
    print(f"  T_max(t_end = {s['T_final']['time']:g} s) = {s['T_final']['max']:.3f}")
    e = s["energy"]
    print(f"  energía [J, con 2π]: entra {e['brake_in_J']:.1f}, sale por convección "
          f"{e['convection_out_J']:.1f}, almacenada {e['stored_J']:.1f}; "
          f"cierre del balance {e['balance_relative_residual']:.1e}")
    for p in s["probes"]:
        print(f"  sonda (r={p['r']:g}, z={p['z']:g}): T_max = {p['T_max']:.3f} en "
              f"t = {p['time_of_max']:.3f} s; T_final = {p['T_final']:.3f}")
    if s["unverified_materials"]:
        print(f"  Materiales NO verificados: {', '.join(s['unverified_materials'])}")
    print(f"  tiempo de cómputo: {s['compute_time_s']:.2f} s   salidas: {res.output_dir}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Punto de entrada de la CLI."""
    # La consola de Windows usa cp1252 por defecto y no puede imprimir Σ, r̄, ²...
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(prog="python -m discfem", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, help_ in [
        ("run", "corrida completa"),
        ("mesh", "solo malla + figura"),
        ("study", "estudio de convergencia (Fase 6)"),
    ]:
        p = sub.add_parser(name, help=help_)
        p.add_argument("config", type=Path, help="archivo YAML de configuración")
        p.add_argument("--set", action="append", default=[], metavar="CLAVE=VALOR",
                       help="override con clave con puntos, p. ej. mesh.refine_factor=2")
    args = parser.parse_args(argv)
    overrides = _parse_set(args.set)

    warnings.formatwarning = _format_warning
    for cat in (UnverifiedMaterialWarning, PropertyRangeWarning):
        warnings.simplefilter("always", cat)
    try:
        if args.command == "mesh":
            return cmd_mesh(args.config, overrides)
        if args.command == "run":
            return cmd_run(args.config, overrides)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    print(f"El comando '{args.command}' aún no está implementado (ver PLAN.md).", file=sys.stderr)
    return 1
