"""Posproceso: sondas, CSV, figuras, VTK y ``summary.json`` de una corrida (PLAN.md Fase 5).

Salidas en ``output.dir``:
- ``config.yaml``: configuración EFECTIVA (YAML + overrides) con ``materials.library`` apuntando
  a la copia ``materiales.yaml`` de la misma carpeta: basta con ella para repetir la corrida.
- ``tmax.csv`` (t, T_max, nodo, r, z) y ``probes.csv`` (t y T en cada sonda).
- ``summary.json`` (ver ``solver.make_summary``).
- PNG: ``mesh.png``, ``field.png`` (T(r, z) en 4 instantes, reflejado en z para ver el espesor
  completo), ``tmax.png`` y ``depth_profiles.png``.
- VTK: ``field_XXXX.vtu`` en los mismos 4 instantes y ``field.pvd`` para abrirlos en ParaView.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import yaml

from .geometry import BOUNDARY_GROUPS
from .mesh import Mesh

_GROUP_COLORS = {
    "brake": "#d62728",
    "top_inner": "#1f77b4",
    "top_outer": "#17becf",
    "r_inner": "#2ca02c",
    "r_outer": "#9467bd",
    "symmetry": "#ff7f0e",
}


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def plot_mesh(mesh: Mesh, path: str | Path, title: str = "", zoom_r: float | None = None) -> Path:
    """Dibuja la malla con cada grupo de frontera en un color y guarda la figura en ``path``.

    Panel superior: sección completa (escala real r:z). Panel inferior: ampliación alrededor
    de ``zoom_r`` (por defecto el extremo interior de la pista), para ver el refinamiento bajo
    ``brake`` y su transición.
    """
    plt = _plt()
    from matplotlib.collections import LineCollection

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    r, z = mesh.nodes[:, 0], mesh.nodes[:, 1]
    e = z.max()
    if zoom_r is None:
        zoom_r = mesh.nodes[np.unique(mesh.edges["brake"]), 0].min()

    fig, axes = plt.subplots(2, 1, figsize=(12, 6.5), gridspec_kw={"height_ratios": [1, 2]})
    multi = len(mesh.region_names) > 1
    for ax in axes:
        if multi:  # sombrea cada región de material
            ax.tripcolor(r, z, mesh.tris, facecolors=mesh.region.astype(float), cmap="Pastel1",
                         vmin=0, vmax=8, alpha=0.8)
        ax.triplot(r, z, mesh.tris, color="0.55", lw=0.3)
        for g in BOUNDARY_GROUPS:
            seg = mesh.nodes[mesh.edges[g]]
            ax.add_collection(LineCollection(seg, colors=_GROUP_COLORS[g], lw=2.2, label=g))
        ax.set_aspect("equal")
        ax.set_xlabel("r [m]")
        ax.set_ylabel("z [m]")
    axes[0].set_xlim(r.min() - 0.01 * e, r.max() + 0.01 * e)
    axes[0].set_title(
        f"{title}  —  {mesh.n_nodes} nodos, {mesh.n_elements} triángulos".strip(" —")
    )
    zmin = z.min()
    axes[1].set_xlim(zoom_r - 1.5 * (e - zmin), zoom_r + 1.5 * (e - zmin))
    axes[1].set_ylim(zmin - 0.05 * (e - zmin), e + 0.05 * (e - zmin))
    axes[1].set_title(f"Ampliación alrededor de r = {zoom_r:g} m")
    handles, labels = axes[0].get_legend_handles_labels()
    if multi:
        from matplotlib.patches import Patch

        cmap = plt.get_cmap("Pastel1")
        handles += [Patch(color=cmap(j / 8), label=f"región {n}")
                    for j, n in enumerate(mesh.region_names)]
        labels += [f"región {n}" for n in mesh.region_names]
    fig.legend(handles, labels, loc="lower center", ncol=min(len(labels), 8), fontsize=9,
               frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


# ------------------------------------------------------------------------------ sondas

def interpolation_matrix(mesh: Mesh, points, allow_outside: bool = False):
    """Matriz P (n_puntos × N) tal que P @ T es la interpolación P1 de T en ``points``.

    Para cada punto se calculan sus coordenadas baricéntricas en TODOS los triángulos
    (vectorizado) y se usa el primero que lo contiene (todas ≥ −1e-10). Un punto fuera del
    dominio produce ``ValueError``; con ``allow_outside=True`` su fila queda vacía y se
    devuelve ``(P, dentro)`` con la máscara de puntos interiores.
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 2)
    inside_mask = np.ones(pts.shape[0], dtype=bool)
    if pts.shape[0] == 0:
        P0 = sp.csr_matrix((0, mesh.n_nodes))
        return (P0, inside_mask) if allow_outside else P0
    p = mesh.nodes[mesh.tris]
    v0, v1, v2 = p[:, 0], p[:, 1], p[:, 2]
    d1, d2 = v1 - v0, v2 - v0
    det = d1[:, 0] * d2[:, 1] - d2[:, 0] * d1[:, 1]
    rows, cols, vals = [], [], []
    for i, (r, z) in enumerate(pts):
        dr, dz = r - v0[:, 0], z - v0[:, 1]
        l1 = (dr * d2[:, 1] - d2[:, 0] * dz) / det
        l2 = (d1[:, 0] * dz - dr * d1[:, 1]) / det
        lam = np.stack([1.0 - l1 - l2, l1, l2], axis=1)
        inside = np.flatnonzero(np.all(lam >= -1e-10, axis=1))
        if inside.size == 0:
            if allow_outside:
                inside_mask[i] = False
                continue
            raise ValueError(f"La sonda (r={r:g}, z={z:g}) está fuera del dominio.")
        e = inside[0]
        w = np.clip(lam[e], 0.0, None)
        rows += [i, i, i]
        cols += list(mesh.tris[e])
        vals += list(w / w.sum())
    P = sp.csr_matrix((vals, (rows, cols)), shape=(pts.shape[0], mesh.n_nodes))
    return (P, inside_mask) if allow_outside else P


def snapshot_indices(saved_times: np.ndarray, t_peak: float, t_f: float | None,
                     t_end: float) -> list[int]:
    """Índices (en ``saved_times``) de los instantes para figuras y VTK: t_f/4, instante de
    T_max, t_f y t_end (sin t_f: t_end/4, T_max, t_end/2, t_end). Se toma el más cercano."""
    targets = [t_f / 4, t_peak, t_f, t_end] if t_f else [t_end / 4, t_peak, t_end / 2, t_end]
    return sorted({int(np.argmin(np.abs(saved_times - t))) for t in targets})


# ------------------------------------------------------------------------------ figuras

def plot_field_snapshots(mesh: Mesh, times, fields, path, T_ref: float | None = None) -> Path:
    """T(r, z) en varios instantes, reflejado respecto a z = 0 para ver el espesor completo
    del disco; misma escala de color en todos los paneles."""
    plt = _plt()
    from matplotlib.tri import Triangulation

    r, z = mesh.nodes[:, 0], mesh.nodes[:, 1]
    tri_up = Triangulation(r, z, mesh.tris)
    tri_dn = Triangulation(r, -z, mesh.tris[:, [0, 2, 1]])
    # Reflejar solo si el plano de simetría es z = 0 (rectángulo y perfiles con 'symmetry'
    # sobre z = 0); si no, se dibuja solo la sección.
    sym = np.unique(mesh.edges["symmetry"])
    mirror = sym.size > 0 and np.all(mesh.nodes[sym, 1] == 0.0)
    vmin = min(f.min() for f in fields) if T_ref is None else min(T_ref, *(f.min() for f in fields))
    vmax = max(f.max() for f in fields)
    levels = np.linspace(vmin, vmax if vmax > vmin else vmin + 1e-9, 41)
    fig, axes = plt.subplots(len(fields), 1, figsize=(10, 1.9 * len(fields) + 0.8),
                             squeeze=False, sharex=True, layout="constrained")
    cs = None
    for ax, t, f in zip(axes[:, 0], times, fields, strict=True):
        cs = ax.tricontourf(tri_up, f, levels=levels, cmap="inferno")
        if mirror:
            ax.tricontourf(tri_dn, f, levels=levels, cmap="inferno")
            ax.axhline(0.0, color="w", lw=0.5, ls=":")
        ax.set_aspect("equal")
        ax.set_ylabel("z [m]")
        ax.set_title(f"t = {t:.3f} s   (T máx = {f.max():.1f})", fontsize=10)
    axes[-1, 0].set_xlabel("r [m]")
    fig.colorbar(cs, ax=axes[:, 0].tolist(), label="T", shrink=0.9)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return Path(path)


def plot_tmax(times, T_max, t_f, path, probes=(), probe_times=None, probe_T=None) -> Path:
    """T_max(t) y la temperatura en las sondas; línea vertical en t_f."""
    plt = _plt()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(times, T_max, lw=2, label="T máx (nodal)")
    for j, (r, z) in enumerate(probes):
        ax.plot(probe_times, probe_T[:, j], lw=1.2, label=f"sonda r={r:g}, z={z:g}")
    if t_f:
        ax.axvline(t_f, color="0.4", ls="--", lw=1, label=f"t_f = {t_f:.3g} s")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("T")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return Path(path)


def plot_depth_profiles(mesh: Mesh, r_line: float, times, fields, path) -> Path:
    """Perfiles T(z) en r = ``r_line`` (el centro de la pista) en varios instantes."""
    plt = _plt()
    zs = np.linspace(mesh.nodes[:, 1].min(), mesh.nodes[:, 1].max(), 101)
    P, inside = interpolation_matrix(mesh, np.column_stack([np.full_like(zs, r_line), zs]),
                                     allow_outside=True)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for t, f in zip(times, fields, strict=True):
        ax.plot((P @ f)[inside], zs[inside] * 1e3, label=f"t = {t:.2f} s")
    ax.set_xlabel("T")
    ax.set_ylabel("z [mm]  (z = e: pista; z = 0: plano medio)")
    ax.set_title(f"Perfiles en profundidad, r = {r_line:g} m")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return Path(path)


# ------------------------------------------------------------------------------ escritura

def write_vtk(mesh: Mesh, times, fields, out_dir: Path) -> list[Path]:
    """Un ``.vtu`` por instante (puntos (r, z, 0), campo nodal T, región por celda) y un
    ``field.pvd`` que los agrupa como serie temporal para ParaView."""
    import meshio

    pts = np.column_stack([mesh.nodes, np.zeros(mesh.n_nodes)])
    files = []
    for k, (t, f) in enumerate(zip(times, fields, strict=True)):
        name = f"field_{k:04d}.vtu"
        meshio.write(
            out_dir / name,
            meshio.Mesh(pts, [("triangle", mesh.tris)], point_data={"T": np.asarray(f)},
                        cell_data={"region": [mesh.region]}),
        )
        files.append((t, name))
    entries = "\n".join(f'    <DataSet timestep="{t:.12g}" part="0" file="{n}"/>' for t, n in files)
    pvd = ('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n  <Collection>\n'
           f"{entries}\n  </Collection>\n</VTKFile>\n")
    (out_dir / "field.pvd").write_text(pvd, encoding="utf-8")
    return [out_dir / n for _, n in files]


def write_outputs(result, library_path: str | Path) -> Path:
    """Escribe todas las salidas de la corrida en ``config.output.dir`` y devuelve la carpeta."""
    cfg, tr, mesh = result.config, result.transient, result.mesh
    out = Path(cfg.output.dir)
    out.mkdir(parents=True, exist_ok=True)

    # Configuración efectiva + copia de la librería: la corrida se repite solo con esto.
    raw = json.loads(json.dumps(cfg.raw))  # copia profunda
    raw.setdefault("materials", {})["library"] = "materiales.yaml"
    if Path(library_path).resolve() != (out / "materiales.yaml").resolve():
        shutil.copyfile(library_path, out / "materiales.yaml")
    with (out / "config.yaml").open("w", encoding="utf-8") as fh:
        fh.write("# Configuración efectiva de la corrida (YAML original + overrides).\n")
        yaml.safe_dump(raw, fh, allow_unicode=True, sort_keys=False)

    with (out / "summary.json").open("w", encoding="utf-8") as fh:
        json.dump(result.summary, fh, indent=2, ensure_ascii=False)

    formats = set(cfg.output.formats)
    if "csv" in formats:
        with (out / "tmax.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["t", "T_max", "node", "r", "z"])
            for t, T, n in zip(tr.times, tr.T_max, tr.T_max_node, strict=True):
                w.writerow([repr(float(t)), repr(float(T)), int(n),
                            repr(float(mesh.nodes[n, 0])), repr(float(mesh.nodes[n, 1]))])
        with (out / "probes.csv").open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["t", *[f"T(r={r:g};z={z:g})" for r, z in cfg.output.probes]])
            for i, t in enumerate(tr.times):
                w.writerow([repr(float(t)), *[repr(float(v)) for v in result.probe_T[i]]])

    t_peak = result.summary["T_max"]["time"]
    t_f = result.summary["time"]["t_f"]
    idx = snapshot_indices(tr.saved_times, t_peak, t_f, float(tr.times[-1]))
    snap_t = [float(tr.saved_times[i]) for i in idx]
    snap_T = [tr.saved_T[i] for i in idx]
    if "png" in formats:
        plot_mesh(mesh, out / "mesh.png", title=f"{cfg.case} [{cfg.mesh.backend}]")
        plot_field_snapshots(mesh, snap_t, snap_T, out / "field.png", T_ref=cfg.T0)
        plot_tmax(tr.times, tr.T_max, t_f, out / "tmax.png", cfg.output.probes,
                  tr.times, result.probe_T)
        r_mid = 0.5 * (cfg.geometry.r_a + cfg.geometry.r_b)
        plot_depth_profiles(mesh, r_mid, snap_t, snap_T, out / "depth_profiles.png")
    if "vtk" in formats:
        write_vtk(mesh, snap_t, snap_T, out)
    return out
