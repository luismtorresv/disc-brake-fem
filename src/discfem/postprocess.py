"""Posproceso y figuras. Por ahora solo el dibujo de la malla (Fase 1); el resto llega en la
Fase 5."""

from __future__ import annotations

from pathlib import Path

import numpy as np

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


def plot_mesh(mesh: Mesh, path: str | Path, title: str = "", zoom_r: float | None = None) -> Path:
    """Dibuja la malla con cada grupo de frontera en un color y guarda la figura en ``path``.

    Panel superior: sección completa (escala real r:z). Panel inferior: ampliación alrededor
    de ``zoom_r`` (por defecto el extremo interior de la pista), para ver el refinamiento bajo
    ``brake`` y su transición.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    r, z = mesh.nodes[:, 0], mesh.nodes[:, 1]
    e = z.max()
    if zoom_r is None:
        zoom_r = mesh.nodes[np.unique(mesh.edges["brake"]), 0].min()

    fig, axes = plt.subplots(2, 1, figsize=(12, 6.5), gridspec_kw={"height_ratios": [1, 2]})
    for ax in axes:
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
    axes[1].set_xlim(zoom_r - 1.5 * e, zoom_r + 1.5 * e)
    axes[1].set_ylim(-0.05 * e, 1.05 * e)
    axes[1].set_title(f"Ampliación alrededor de r = {zoom_r:g} m")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=9, frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path
