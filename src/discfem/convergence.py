"""Estudios de convergencia (malla y Δt) y barridos paramétricos (PLAN.md Fase 6).

Definiciones (PLAN.md §5, Fase 6). Para una sucesión de niveles i = 0, 1, 2, … con una
magnitud Q_i (T_max, o T_max en una sonda):

    ε_i = |Q_i − Q_{i−1}| / |Q_i| × 100 %                          (i ≥ 1)
    p_i = log(|Q_{i−2} − Q_{i−1}| / |Q_{i−1} − Q_i|) / log(ρ)       (i ≥ 2)

con ρ la razón de refinamiento entre niveles consecutivos (ρ = 2 para refine_factor 1, 2, 4,
8 o Δt, Δt/2, Δt/4), que debe ser constante en el trío. Con ρ = 2 es exactamente la fórmula
del PLAN, p = log2(|T_h − T_h/2| / |T_h/2 − T_h/4|).

Nivel seleccionado: el primer nivel i ≥ 1 cuyo ε_i < tolerancia en T_max Y en todas las sondas.

Nota: ε se mide sobre la temperatura (en las unidades del caso), como define el PLAN. Como T
incluye T0, ε sobre el incremento T − T0 es mayor; se reporta también (columna ``eps_rise``).
"""

from __future__ import annotations

import csv
import json
import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from .config import load_config
from .materials import resolve_library_path
from .solver import run

_STUDY_KEYS = {"name", "base", "base_overrides", "output_dir", "mesh", "timestep", "sweeps",
               "gamma_map",
               "use_selected_mesh"}


# ------------------------------------------------------------------------------ fórmulas

def relative_change_percent(prev: float, cur: float) -> float:
    """ε = |cur − prev| / |cur| × 100."""
    return abs(cur - prev) / abs(cur) * 100.0


def observed_order(q0: float, q1: float, q2: float, ratio: float = 2.0) -> float | None:
    """p = log(|q0 − q1| / |q1 − q2|) / log(ratio); None si alguna diferencia es 0."""
    d01, d12 = abs(q0 - q1), abs(q1 - q2)
    if d01 == 0 or d12 == 0 or ratio <= 1:
        return None
    return math.log(d01 / d12) / math.log(ratio)


@dataclass
class Level:
    """Un nivel de un estudio o un punto de un barrido."""

    label: str
    value: Any
    summary: dict
    times: np.ndarray
    T_max_t: np.ndarray
    fields: dict = field(default_factory=dict)

    @property
    def T_max(self) -> float:
        return self.summary["T_max"]["value"]

    @property
    def probe_max(self) -> list[float]:
        return [p["T_max"] for p in self.summary["probes"]]


def _run_level(base: Path, overrides: dict, label: str, value: Any) -> Level:
    res = run(base, {**overrides, "output.save_every": 10**9}, write=False)
    return Level(label=label, value=value, summary=res.summary, times=res.transient.times,
                 T_max_t=res.transient.T_max)


def _ratios(values: list[float]) -> list[float]:
    return [values[i] / values[i - 1] for i in range(1, len(values))]


def convergence_table(levels: list[Level], refine: list[float], T0: float,
                      tol: float | None = None) -> tuple[list[dict], int | None]:
    """Filas de la tabla de convergencia y el índice del nivel seleccionado (o None).

    ``refine`` es el factor de refinamiento de cada nivel (refine_factor o Δt₀/Δt).
    """
    rows, selected = [], None
    ratios = _ratios(refine)
    for i, lv in enumerate(levels):
        s = lv.summary
        row = {
            "level": i, "label": lv.label, "value": lv.value,
            "n_nodes": s["mesh"]["n_nodes"], "n_elements": s["mesh"]["n_elements"],
            "dt": s["time"]["dt"], "n_steps": s["time"]["n_steps"],
            "T_max": lv.T_max, "t_of_T_max": s["T_max"]["time"],
            "T_final_max": s["T_final"]["max"],
            "eps_percent": None, "eps_rise_percent": None, "eps_probes_max_percent": None,
            "order_p": None, "compute_time_s": s["compute_time_s"],
        }
        for j, v in enumerate(lv.probe_max):
            row[f"probe{j}_T_max"] = v
        if i >= 1:
            prev = levels[i - 1]
            row["eps_percent"] = relative_change_percent(prev.T_max, lv.T_max)
            row["eps_rise_percent"] = relative_change_percent(prev.T_max - T0, lv.T_max - T0)
            eps_p = [relative_change_percent(a, b) for a, b in zip(prev.probe_max, lv.probe_max,
                                                                    strict=True)]
            row["eps_probes_max_percent"] = max(eps_p) if eps_p else 0.0
            if tol is not None and selected is None and row["eps_percent"] < tol \
                    and row["eps_probes_max_percent"] < tol:
                selected = i
        if i >= 2 and math.isclose(ratios[i - 2], ratios[i - 1], rel_tol=1e-9):
            row["order_p"] = observed_order(levels[i - 2].T_max, levels[i - 1].T_max, lv.T_max,
                                            ratios[i - 1])
        rows.append(row)
    return rows, selected


# ------------------------------------------------------------------------------ estudios

def mesh_study(base: Path, refine_factors: list[float], overrides: dict, tol: float):
    """Niveles con ``mesh.refine_factor`` = cada valor, mismo Δt (el de la config base)."""
    levels = [_run_level(base, {**overrides, "mesh.refine_factor": rf}, f"rf={rf:g}", rf)
              for rf in refine_factors]
    T0 = levels[0].summary["T0"]
    rows, sel = convergence_table(levels, list(refine_factors), T0, tol)
    return levels, rows, sel


def timestep_study(base: Path, factors: list[int], overrides: dict):
    """Niveles con Δt/f para cada factor f (Δt base = el de la config)."""
    cfg = load_config(base, overrides)
    levels = []
    for f in factors:
        if cfg.time.dt == "auto":
            ov = {"time.steps_per_tf": int(cfg.time.steps_per_tf * f)}
        else:
            ov = {"time.dt": float(cfg.time.dt) / f}
        levels.append(_run_level(base, {**overrides, **ov}, f"Δt/{f:g}", f))
    rows, _ = convergence_table(levels, [float(f) for f in factors], cfg.T0, None)
    return levels, rows


def parameter_sweep(base: Path, key: str, values: list, overrides: dict):
    """Una corrida por valor de ``key`` (clave con puntos)."""
    levels = [_run_level(base, {**overrides, key: v}, f"{key}={v}", v) for v in values]
    rows = []
    for lv in levels:
        s = lv.summary
        row = {"value": lv.value, "T_max": lv.T_max, "t_of_T_max": s["T_max"]["time"],
               "T_final_max": s["T_final"]["max"], "t_f": s["time"]["t_f"],
               "brake_in_J": s["energy"]["brake_in_J"], "gamma": s["load"]["gamma"],
               "n_nodes": s["mesh"]["n_nodes"], "compute_time_s": s["compute_time_s"]}
        for j, v in enumerate(lv.probe_max):
            row[f"probe{j}_T_max"] = v
        rows.append(row)
    return levels, rows


# ------------------------------------------------------------------------------ salidas

def _fmt(v) -> str:
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:.6g}" if abs(v) >= 1e-3 or v == 0 else f"{v:.3e}"
    return str(v)


def write_table(rows: list[dict], path_stem: Path, columns: list[str] | None = None) -> None:
    """Escribe la tabla como CSV (todas las columnas) y Markdown (``columns``)."""
    cols_all = list(rows[0].keys()) if rows else []
    with path_stem.with_suffix(".csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols_all)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r[k] is None else r[k]) for k in cols_all})
    cols = columns or cols_all
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(_fmt(r.get(c)) for c in cols) + " |" for r in rows]
    path_stem.with_suffix(".md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _save(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=150)


def plot_convergence(levels: list[Level], rows: list[dict], stem: Path, xlabel: str,
                     x_key: str, t_f: float | None, title: str) -> None:
    """Panel izquierdo: T_max(t) de cada nivel, ampliado alrededor del pico. Panel derecho:
    ε (%) frente a la resolución, en escala log-log."""
    plt = _plt()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for lv in levels:
        a1.plot(lv.times, lv.T_max_t, lw=1.4, label=lv.label)
    t_pk = levels[-1].summary["T_max"]["time"]
    a1.set_xlim(0, (t_f or t_pk) * 1.3 if (t_f or t_pk) else None)
    hi = max(lv.T_max for lv in levels)
    lo = min(lv.T_max for lv in levels)
    span = max(hi - lo, 1e-6)
    a1.set_ylim(lo - 6 * span, hi + 2 * span)
    a1.set_xlabel("t [s]")
    a1.set_ylabel("T máx")
    a1.set_title("T máx(t), ampliado cerca del pico")
    a1.grid(alpha=0.3)
    a1.legend(fontsize=8)
    xs = [r[x_key] for r in rows[1:]]
    eps = [r["eps_percent"] for r in rows[1:]]
    eps_p = [r["eps_probes_max_percent"] for r in rows[1:]]
    a2.loglog(xs, eps, "o-", label="ε T máx")
    if any(v for v in eps_p):
        a2.loglog(xs, eps_p, "s--", label="ε máx. en sondas")
    a2.axhline(1.0, color="0.4", ls=":", lw=1, label="1 %")
    from matplotlib.ticker import NullFormatter

    a2.xaxis.set_minor_formatter(NullFormatter())
    a2.set_xticks(xs, [f"{x:g}" for x in xs])
    a2.set_xlabel(xlabel)
    a2.set_ylabel("ε [%]")
    a2.set_title("Cambio relativo entre niveles consecutivos")
    a2.grid(alpha=0.3, which="both")
    a2.legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    _save(fig, stem)
    plt.close(fig)


def plot_sweep(levels: list[Level], rows: list[dict], stem: Path, label: str, title: str) -> None:
    """Panel izquierdo: T_max(t) para cada valor. Panel derecho: T_max frente al parámetro."""
    plt = _plt()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for lv in levels:
        a1.plot(lv.times, lv.T_max_t, lw=1.4, label=f"{label} = {lv.value}")
    a1.set_xlabel("t [s]")
    a1.set_ylabel("T máx")
    a1.grid(alpha=0.3)
    a1.legend(fontsize=8)
    vals = [r["value"] for r in rows]
    numeric = all(isinstance(v, int | float) for v in vals)
    x = vals if numeric else list(range(len(vals)))
    a2.plot(x, [r["T_max"] for r in rows], "o-", label="T máx (pico)")
    a2.plot(x, [r["T_final_max"] for r in rows], "s--", label="T máx al final")
    if not numeric:
        a2.set_xticks(x, [str(v) for v in vals], rotation=20)
    a2.set_xlabel(label)
    a2.set_ylabel("T")
    a2.grid(alpha=0.3)
    a2.legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    _save(fig, stem)
    plt.close(fig)


def plot_gamma_map(sweeps: dict, pairs: list[dict], stem: Path, title: str) -> None:
    """T_max y T(t_end) frente a γ. Para cada par: una curva (barrido continuo de ``load.gamma``)
    y, encima, los puntos de un barrido de materiales (cada uno en su γ calculado), con el
    nombre del material. Así se ve cuánto mueve la temperatura la elección de la pastilla."""
    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i, pr in enumerate(pairs):
        c = colors[i % len(colors)]
        curve = sorted(sweeps[pr["curve"]]["rows"], key=lambda r: r["gamma"])
        pts = sweeps[pr["points"]]["rows"]
        for ax, key in zip(axes, ("T_max", "T_final_max"), strict=True):
            ax.plot([r["gamma"] for r in curve], [r[key] for r in curve], "-", color=c,
                    lw=1.6, label=pr.get("label", pr["curve"]))
            ax.plot([r["gamma"] for r in pts], [r[key] for r in pts], "o", color=c, ms=7,
                    mec="k")
            names = pr.get("point_labels") or {}
            for r in pts:
                # Etiquetas del primer par abajo a la derecha y del segundo arriba a la
                # izquierda, para que no se encimen cuando dos puntos tienen γ parecido.
                right = i % 2 == 0
                ax.annotate(names.get(str(r["value"]), str(r["value"])), (r["gamma"], r[key]),
                            textcoords="offset points", xytext=(7, -12) if right else (-7, 7),
                            ha="left" if right else "right", fontsize=7.5, color=c)
    axes[0].set_ylabel("T máx (pico)")
    axes[1].set_ylabel("T máx al final de la simulación")
    for ax in axes:
        ax.set_xlabel("γ (fracción del calor de fricción que entra al disco)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    _save(fig, stem)
    plt.close(fig)


# ------------------------------------------------------------------------------ orquestación

MESH_COLS = ["label", "n_nodes", "n_elements", "T_max", "t_of_T_max", "eps_percent",
             "eps_rise_percent", "eps_probes_max_percent", "order_p", "compute_time_s"]
DT_COLS = ["label", "dt", "n_steps", "T_max", "t_of_T_max", "T_final_max", "eps_percent",
           "order_p", "compute_time_s"]
SWEEP_COLS = ["value", "gamma", "T_max", "t_of_T_max", "T_final_max", "t_f", "brake_in_J",
              "compute_time_s"]


def load_study(path: str | Path) -> dict:
    """Lee y valida el YAML de estudio (sección ``study``)."""
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"No existe el archivo de estudio '{path}'.")
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or not isinstance(raw.get("study"), dict):
        raise ValueError(f"'{path}' debe tener una sección 'study'.")
    st = raw["study"]
    unknown = set(st) - _STUDY_KEYS
    if unknown:
        raise ValueError(f"Claves desconocidas en 'study': {sorted(unknown)}.")
    if "base" not in st:
        raise ValueError("'study.base' (config base del caso) es obligatoria.")
    base = Path(st["base"])
    if not base.is_file() and (path.parent / base).is_file():
        base = path.parent / base
    if not base.is_file():
        raise ValueError(f"No existe la config base '{st['base']}'.")
    st = dict(st, base=base, path=path)
    if "mesh" in st:
        rfs = st["mesh"].get("refine_factors")
        if not rfs or any(float(r) <= 0 for r in rfs) or list(rfs) != sorted(rfs):
            raise ValueError("study.mesh.refine_factors debe ser una lista creciente de números > 0.")
    if "timestep" in st:
        fs = st["timestep"].get("factors")
        if not fs or any(int(f) != f or f < 1 for f in fs) or list(fs) != sorted(fs):
            raise ValueError("study.timestep.factors debe ser una lista creciente de enteros ≥ 1.")
    names = []
    for sw in st.get("sweeps", []) or []:
        if not {"name", "key", "values"} <= set(sw):
            raise ValueError("Cada barrido necesita 'name', 'key' y 'values'.")
        unknown = set(sw) - {"name", "key", "values", "label", "title", "overrides"}
        if unknown:
            raise ValueError(f"Claves desconocidas en el barrido '{sw['name']}': {sorted(unknown)}.")
        if not isinstance(sw.get("overrides") or {}, dict):
            raise ValueError(f"'overrides' del barrido '{sw['name']}' debe ser un diccionario.")
        names.append(sw["name"])
    gm = st.get("gamma_map")
    if gm is not None:
        pairs = gm.get("pairs") if isinstance(gm, dict) else None
        if not pairs:
            raise ValueError("study.gamma_map necesita 'pairs': [{curve, points, label}, ...].")
        for pr in pairs:
            if not isinstance(pr.get("point_labels") or {}, dict):
                raise ValueError("gamma_map: 'point_labels' debe ser un diccionario {valor: texto}.")
            for role in ("curve", "points"):
                if pr.get(role) not in names:
                    raise ValueError(f"gamma_map: '{pr.get(role)}' ({role}) no es un barrido "
                                     f"definido; barridos: {names}.")
    return st


def run_study(path: str | Path, log=print) -> dict:
    """Ejecuta el estudio descrito en ``path`` y escribe tablas (CSV, Markdown), figuras
    (PDF vectorial + PNG) y ``study_summary.json`` en ``study.output_dir``."""
    st = load_study(path)
    base: Path = st["base"]
    name = st.get("name", Path(path).stem)
    out = Path(st.get("output_dir", f"results/{name}"))
    out.mkdir(parents=True, exist_ok=True)
    ov = dict(st.get("base_overrides") or {})
    base_cfg = load_config(base, ov)

    # Reproducibilidad: copia del estudio, de la config base y de la librería.
    shutil.copyfile(st["path"], out / "study.yaml")
    shutil.copyfile(base, out / "base_config.yaml")
    shutil.copyfile(resolve_library_path(base_cfg.materials.library, base_cfg.config_dir),
                    out / "materiales.yaml")

    summary: dict[str, Any] = {"name": name, "base": str(base), "base_overrides": ov}
    t_f_ref = None

    if "mesh" in st:
        rfs = [float(r) for r in st["mesh"]["refine_factors"]]
        tol = float(st["mesh"].get("tolerance_percent", 1.0))
        log(f"Estudio de malla: refine_factor = {rfs}")
        levels, rows, sel = mesh_study(base, rfs, ov, tol)
        t_f_ref = levels[0].summary["time"]["t_f"]
        write_table(rows, out / "mesh_study", MESH_COLS)
        plot_convergence(levels, rows, out / "mesh_study", "nodos", "n_nodes", t_f_ref,
                         f"Estudio de malla ({base_cfg.mesh.backend}, Δt fijo)")
        summary["mesh"] = {"tolerance_percent": tol, "rows": rows,
                           "selected_level": sel,
                           "selected_refine_factor": None if sel is None else rfs[sel]}
        log((out / "mesh_study.md").read_text(encoding="utf-8"))
        log("Nivel seleccionado: " + ("ninguno cumple la tolerancia" if sel is None
                                      else f"{rows[sel]['label']} (ε < {tol:g} %)"))
        if st.get("use_selected_mesh", True) and sel is not None:
            ov["mesh.refine_factor"] = rfs[sel]

    if "timestep" in st:
        fs = [int(f) for f in st["timestep"]["factors"]]
        log(f"Estudio de Δt: factores {fs}")
        levels, rows = timestep_study(base, fs, ov)
        t_f_ref = t_f_ref or levels[0].summary["time"]["t_f"]
        write_table(rows, out / "timestep_study", DT_COLS)
        plot_convergence(levels, rows, out / "timestep_study", "Δt₀/Δt", "value", t_f_ref,
                         f"Estudio de Δt (θ = {base_cfg.time.theta:g})")
        summary["timestep"] = {"rows": rows}
        log((out / "timestep_study.md").read_text(encoding="utf-8"))

    summary["sweeps"] = {}
    for sw in st.get("sweeps", []) or []:
        log(f"Barrido '{sw['name']}': {sw['key']} = {sw['values']}")
        levels, rows = parameter_sweep(base, sw["key"], list(sw["values"]),
                                       {**ov, **(sw.get("overrides") or {})})
        stem = out / f"sweep_{sw['name']}"
        write_table(rows, stem, SWEEP_COLS)
        plot_sweep(levels, rows, stem, sw.get("label", sw["key"]), sw.get("title", sw["name"]))
        summary["sweeps"][sw["name"]] = {"key": sw["key"], "rows": rows,
                                         "overrides": sw.get("overrides") or {}}
        log(stem.with_suffix(".md").read_text(encoding="utf-8"))

    if st.get("gamma_map"):
        gm = st["gamma_map"]
        plot_gamma_map(summary["sweeps"], gm["pairs"], out / "gamma_map",
                       gm.get("title", "Temperatura frente a γ"))
        log(f"Figura γ: {out / 'gamma_map.pdf'}")

    summary["effective_overrides"] = ov
    (out / "study_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False),
                                            encoding="utf-8")
    summary["output_dir"] = str(out)
    return summary
