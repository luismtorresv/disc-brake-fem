"""Generador de malla con la API de Python de Gmsh (sin pasar por archivo).

Portado de ``legacy/mesh/create_mesh.py``: misma sección y mismo par de campos
``Distance`` + ``Threshold`` sobre la curva ``brake``, pero con los seis grupos físicos de
CLAUDE.md, tamaños leídos del YAML y divididos por ``refine_factor``.
"""

from __future__ import annotations

import numpy as np

from ..geometry import BOUNDARY_GROUPS, RectangleGeometry
from .base import Mesh, MeshOptions, orient_edges


def gmsh_mesh(geo: RectangleGeometry, opts: MeshOptions) -> Mesh:
    """Malla triangular no estructurada con refinamiento junto a ``brake``.

    El tamaño de elemento es h(d) = h_min + (h_max − h_min)·min(d/grading_distance, 1), con d
    la distancia a la curva ``brake`` (campo ``Threshold`` con DistMin = 0). Los vértices
    (r_a, e) y (r_b, e) son puntos de la geometría, así que hay nodos exactamente ahí.
    """
    import gmsh

    h_min, h_max = opts.h_min, opts.h_max
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)  # malla reproducible
        gmsh.model.add("disc_section")
        geo_k = gmsh.model.geo

        pts = [geo_k.addPoint(r, z, 0.0, h_max) for r, z in geo.corners]
        lines = [geo_k.addLine(pts[k], pts[(k + 1) % 6]) for k in range(6)]
        surf = geo_k.addPlaneSurface([geo_k.addCurveLoop(lines)])
        geo_k.synchronize()

        for k, name in enumerate(BOUNDARY_GROUPS):
            gmsh.model.addPhysicalGroup(1, [lines[k]], k + 1, name)
        gmsh.model.addPhysicalGroup(2, [surf], 100, "disc")

        brake_line = lines[BOUNDARY_GROUPS.index("brake")]
        fld = gmsh.model.mesh.field
        fld.add("Distance", 1)
        fld.setNumbers(1, "CurvesList", [brake_line])
        fld.setNumber(1, "Sampling", max(200, int(4 * (geo.r_b - geo.r_a) / h_min)))
        fld.add("Threshold", 2)
        fld.setNumber(2, "InField", 1)
        fld.setNumber(2, "SizeMin", h_min)
        fld.setNumber(2, "SizeMax", h_max)
        fld.setNumber(2, "DistMin", 0.0)
        fld.setNumber(2, "DistMax", opts.grading_distance)
        fld.setAsBackgroundMesh(2)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.ElementOrder", 1)

        gmsh.model.mesh.generate(2)

        tags, coords, _ = gmsh.model.mesh.getNodes()
        xyz = coords.reshape(-1, 3)
        tag_to_idx = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
        tag_to_idx[tags.astype(np.int64)] = np.arange(tags.size)

        _, tri_nodes = gmsh.model.mesh.getElementsByType(2, tag=surf)
        tris = tag_to_idx[tri_nodes.astype(np.int64)].reshape(-1, 3)

        raw_edges = {}
        for k, name in enumerate(BOUNDARY_GROUPS):
            _, line_nodes = gmsh.model.mesh.getElementsByType(1, tag=lines[k])
            raw_edges[name] = tag_to_idx[line_nodes.astype(np.int64)].reshape(-1, 2)
    finally:
        gmsh.finalize()

    # Compactar: solo los nodos usados por triángulos, en orden de aparición estable.
    used = np.unique(tris)
    remap = np.full(xyz.shape[0], -1, dtype=np.int64)
    remap[used] = np.arange(used.size)
    nodes = xyz[used, :2].copy()
    tris = remap[tris]
    raw_edges = {k: remap[v] for k, v in raw_edges.items()}

    # Orientación antihoraria (área > 0).
    p = nodes[tris]
    area2 = (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - (
        p[:, 2, 0] - p[:, 0, 0]
    ) * (p[:, 1, 1] - p[:, 0, 1])
    flip = area2 < 0
    tris[flip] = tris[flip][:, [0, 2, 1]]

    edges = {k: orient_edges(tris, v, nodes.shape[0]) for k, v in raw_edges.items()}
    mesh = Mesh(
        nodes=nodes,
        tris=tris,
        region=np.zeros(tris.shape[0], dtype=np.int64),
        region_names=["disc"],
        edges=edges,
    )
    mesh.validate(BOUNDARY_GROUPS)
    return mesh
