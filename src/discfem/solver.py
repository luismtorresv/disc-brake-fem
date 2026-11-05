"""Orquestación de una corrida: ``run(config, overrides=None) -> Result``.

Pasos: configuración → librería de materiales (con advertencia si hay no verificados) → malla →
γ y carga → ensamblaje (M, K, H, F) → malla temporal con t_f → esquema θ → sondas → resumen →
salidas (``postprocess.write_outputs``).

Energías: el ensamblaje omite 2π (CLAUDE.md). En ``summary`` las energías físicas [J] se
reportan con el 2π explícito; las de ``TransientResult`` están en J/rad.
"""

from __future__ import annotations

import math
import time as _time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from . import __version__
from .config import Config, load_config
from .fem import assemble_robin, assemble_stiffness_mass
from .loads import BrakingLoad, Load, brake_load_vector, make_load
from .materials import (
    Material,
    add_custom_materials,
    check_range,
    load_library,
    load_presets,
    resolve_gamma,
    resolve_library_path,
    resolve_region_materials,
    unverified,
    warn_unverified,
)
from .mesh import Mesh, build_mesh
from .timestepping import TransientResult, integrate, time_grid


@dataclass
class Result:
    """Resultado completo de una corrida.

    Atributos:
        config: configuración efectiva (con overrides).
        mesh: malla usada.
        transient: resultado del esquema θ (campos, T_max(t), energías en J/rad).
        materials: material por región (con overrides aplicados); pad: material de la pastilla.
        gamma: γ usado (None si la carga no es ``braking``); load: modelo de carga.
        probe_T: (n_pasos + 1, n_sondas) temperatura en las sondas en TODOS los pasos.
        summary: diccionario que se escribe como ``summary.json``.
        output_dir: carpeta de resultados (None si no se escribió nada).
    """

    config: Config
    mesh: Mesh
    transient: TransientResult
    materials: dict[str, Material]
    pad: Material | None
    gamma: float | None
    load: Load | None
    probe_T: np.ndarray
    summary: dict
    output_dir: Path | None = None


def run(
    config: Config | str | Path | dict,
    overrides: dict[str, Any] | None = None,
    write: bool = True,
) -> Result:
    """Ejecuta una corrida completa.

    Args:
        config: :class:`Config`, ruta al YAML o diccionario.
        overrides: claves con puntos, p. ej. ``{"mesh.refine_factor": 2}``.
        write: si es True, escribe las salidas en ``output.dir`` (ver ``output.formats``).
    """
    t_start = _time.perf_counter()
    if isinstance(config, Config):
        cfg = load_config(config.raw, overrides) if overrides else config
        if overrides:
            cfg = _with_dir(cfg, config.config_dir)
    else:
        cfg = load_config(config, overrides)

    # --- materiales (advertencia al inicio de cada corrida) ---------------------------
    lib_path = resolve_library_path(cfg.materials.library, cfg.config_dir)
    library = add_custom_materials(load_library(lib_path), cfg.materials.custom)
    regions = dict(cfg.materials.regions)
    pad_name = cfg.materials.pad
    if cfg.materials.preset is not None:
        presets = load_presets(lib_path, library)
        if cfg.materials.preset not in presets:
            raise ValueError(f"materials.preset = '{cfg.materials.preset}' no existe; presets "
                             f"disponibles en la librería: {sorted(presets)}.")
        pr = presets[cfg.materials.preset]
        regions.setdefault("disc", pr["disc"])  # lo explícito en la config tiene prioridad
        pad_name = pad_name or pr.get("pad")
    materials = resolve_region_materials(regions, library, cfg.materials.overrides)
    pad = None
    if pad_name is not None:
        if pad_name not in library:
            raise ValueError(f"materials.pad = '{pad_name}' no está en la librería.")
        pad = library[pad_name]
    if pad is None and cfg.load.get("type") == "braking" \
            and cfg.load.get("gamma", "auto") == "auto" and cfg.groups_with("flux"):
        origin = f"el preset '{cfg.materials.preset}'" if cfg.materials.preset else "la config"
        raise ValueError(f"gamma: auto necesita una pastilla y {origin} no la define; indique "
                         "materials.pad (o un gamma numérico).")
    used = [*materials.values(), *([pad] if pad else [])]
    warn_unverified(used)
    for g, (h_g, _) in cfg.convection.items():
        check_range("h", h_g, f"convección en '{g}'")

    # --- malla, γ y carga ----------------------------------------------------------------
    mesh = build_mesh(cfg.geometry, cfg.mesh)
    flux_groups = cfg.groups_with("flux")
    gamma, load = None, None
    if flux_groups:
        if cfg.load.get("type") == "braking":
            gamma = resolve_gamma(cfg.load.get("gamma", "auto"), mesh, materials, pad)
        load = make_load(cfg.load, cfg.geometry.r_a, cfg.geometry.r_b, gamma)

    # --- ensamblaje ----------------------------------------------------------------------
    axi = cfg.geometry.axisymmetric
    k = {r: m.k for r, m in materials.items()}
    rho_c = {r: m.rho_cp for r, m in materials.items()}
    K, M = assemble_stiffness_mass(mesh, k, rho_c, axi)
    H, F_conv = assemble_convection(mesh, cfg.convection, axi)
    F_brake = brake_load_vector(mesh, load, axi, flux_groups) if load is not None else None

    # --- tiempo ----------------------------------------------------------------------------
    t_f = load.t_f if load is not None else None
    from .postprocess import interpolation_matrix, write_outputs

    P = interpolation_matrix(mesh, cfg.output.probes)
    times = time_grid(cfg.time.t_end, t_f=t_f, dt=cfg.time.dt,
                      steps_per_tf=cfg.time.steps_per_tf)
    tr = integrate(
        M, K, H, cfg.T0, times, F_brake=F_brake, t_f=t_f, F_conv=F_conv,
        theta=cfg.time.theta, rannacher_steps=cfg.time.rannacher_steps,
        save_every=cfg.output.save_every, probe_matrix=P,
    )

    # --- resumen ---------------------------------------------------------------------------
    probe_T = tr.probe_history
    elapsed = _time.perf_counter() - t_start

    result = Result(
        config=cfg, mesh=mesh, transient=tr, materials=materials, pad=pad, gamma=gamma,
        load=load, probe_T=probe_T, summary={},
    )
    result.summary = make_summary(result, used, elapsed)
    if write:
        result.output_dir = write_outputs(result, lib_path)
    return result


def assemble_convection(mesh: Mesh, convection: dict[str, tuple[float, float]], axi: bool):
    """H = Σ_g H_g(h_g) y F_conv = Σ_g H_g·(T∞_g 1), con (h_g, T∞_g) por grupo.

    Si todos los grupos comparten h y T∞ se ensambla en una sola llamada, exactamente como en
    las fases 1–6 (mismo orden de sumas en punto flotante ⇒ resultados idénticos bit a bit).
    """
    if not convection:
        return None, None
    pairs = set(convection.values())
    if len(pairs) == 1:
        h, T_inf = pairs.pop()
        H = assemble_robin(mesh, dict.fromkeys(convection, h), axi)
        return H, H @ np.full(mesh.n_nodes, T_inf)
    H, F = None, np.zeros(mesh.n_nodes)
    for g, (h, T_inf) in convection.items():
        Hg = assemble_robin(mesh, {g: h}, axi)
        H = Hg if H is None else H + Hg
        F += Hg @ np.full(mesh.n_nodes, T_inf)
    return H, F


def _with_dir(cfg: Config, config_dir: Path) -> Config:
    from dataclasses import replace

    return replace(cfg, config_dir=config_dir)


def make_summary(res: Result, used: list[Material], elapsed: float) -> dict:
    """Diccionario de ``summary.json``. Todo es determinista salvo ``compute_time_s``."""
    cfg, tr, mesh = res.config, res.transient, res.mesh
    i_max = int(np.argmax(tr.T_max))
    node = int(tr.T_max_node[i_max])
    two_pi = 2 * math.pi if cfg.geometry.axisymmetric else 1.0
    times = tr.times
    summary: dict[str, Any] = {
        "case": cfg.case,
        "discfem_version": __version__,
        "T_max": {
            "value": float(tr.T_max[i_max]),
            "time": float(times[i_max]),
            "node": node,
            "r": float(mesh.nodes[node, 0]),
            "z": float(mesh.nodes[node, 1]),
        },
        "T_final": {
            "time": float(times[-1]),
            "max": float(tr.T_max[-1]),
            "min": float(tr.T_final.min()),
        },
        "T0": cfg.T0,
        "energy": {
            "note": "Energías físicas en J (con 2π en modo axisimétrico)",
            "balance_relative_residual": tr.balance_relative(),
            "brake_in_J": two_pi * float(tr.energy_brake[-1]),
            "convection_out_J": two_pi * float(tr.energy_conv[-1]),
            "stored_J": two_pi * float(tr.energy_stored[-1]),
        },
        "mesh": {
            "backend": cfg.mesh.backend,
            "refine_factor": cfg.mesh.refine_factor,
            "n_nodes": mesh.n_nodes,
            "n_elements": mesh.n_elements,
        },
        "time": {
            "theta": cfg.time.theta,
            "rannacher_steps": cfg.time.rannacher_steps,
            "dt": float(np.diff(times).max()),
            "n_steps": int(times.size - 1),
            "t_end": float(times[-1]),
            "t_f": None if res.load is None or res.load.t_f is None else float(res.load.t_f),
            "n_factorizations": tr.n_factorizations,
        },
        "load": {
            "type": cfg.load.get("type") if res.load is not None else None,
            "gamma": res.gamma,
        },
        "materials": {
            region: {"name": m.name, "k": m.k, "rho": m.rho, "cp": m.cp,
                     "source": m.source, "verified": m.verified}
            for region, m in res.materials.items()
        },
        "pad": None if res.pad is None else {
            "name": res.pad.name, "k": res.pad.k, "rho": res.pad.rho, "cp": res.pad.cp,
            "source": res.pad.source, "verified": res.pad.verified,
        },
        "unverified_materials": [m.name for m in unverified(used)],
        "probes": [
            {
                "r": r, "z": z,
                "T_max": float(res.probe_T[:, j].max()),
                "time_of_max": float(times[int(np.argmax(res.probe_T[:, j]))]),
                "T_final": float(res.probe_T[-1, j]),
            }
            for j, (r, z) in enumerate(cfg.output.probes)
        ],
        "compute_time_s": elapsed,
        "geometry": {
            "type": cfg.geometry.type,
            "r_a": cfg.geometry.r_a,
            "r_b": cfg.geometry.r_b,
            "regions": mesh.region_names,
            "area_r_moment": cfg.geometry.r_moment(),
        },
        "convection": {g: {"h": h, "T_inf": T} for g, (h, T) in cfg.convection.items()},
    }
    if isinstance(res.load, BrakingLoad):
        summary["load"].update({
            "pressure_profile": res.load.profile,
            "power_per_face_t0_W": float(res.load.power_per_face(0.0)),
            "energy_per_face_J": res.load.energy_per_face,
            "gamma_energy_per_face_J": res.gamma * res.load.energy_per_face,
            "q_max_t0_W_m2": float(np.max(res.load.shape(
                np.linspace(cfg.geometry.r_a, cfg.geometry.r_b, 201)))),
        })
    return summary
