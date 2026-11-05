"""Generador de malla con la API de Python de Gmsh (sin pasar por archivo).

Dos caminos:

1. **Sin regiones** (rectángulo o perfil homogéneo): kernel ``geo``. Se crean los n vértices
   del borde, n líneas, un grupo físico por grupo de frontera y los campos ``Distance`` +
   ``Threshold`` sobre las líneas ``brake``. Para el rectángulo, la secuencia de llamadas es la
   misma de las fases 1–6, así que la malla es idéntica bit a bit (prueba
   ``test_rectangle_unchanged``).
2. **Con regiones** (Fase 7): kernel OCC. Se crean la sección y cada región como superficies
   y se aplica ``fragment``, que parte la sección en superficies con aristas compartidas: la
   malla resulta conforme en las interfaces entre materiales. Tras mallar, la clasificación es
   geométrica: cada triángulo va a la región que contiene su centroide y cada arista de borde
   al lado del perfil que contiene su punto medio. Como control, el área de cada región en la
   malla debe coincidir con la del polígono (detecta regiones solapadas o fuera del dominio).
"""

from __future__ import annotations

import numpy as np

from ..geometry import (
    BOUNDARY_GROUPS,
    Geometry,
    point_segment_distance,
    points_in_polygon,
    polygon_area,
)
from .base import Mesh, MeshOptions, boundary_edges, orient_edges


def _set_size_field(gmsh, brake_curves: list[int], geo: Geometry, opts: MeshOptions) -> None:
    """h(d) = h_min + (h_max − h_min)·min(d/grading_distance, 1), d = distancia a ``brake``."""
    fld = gmsh.model.mesh.field
    fld.add("Distance", 1)
    fld.setNumbers(1, "CurvesList", brake_curves)
    fld.setNumber(1, "Sampling", max(200, int(4 * (geo.r_b - geo.r_a) / opts.h_min)))
    fld.add("Threshold", 2)
    fld.setNumber(2, "InField", 1)
    fld.setNumber(2, "SizeMin", opts.h_min)
    fld.setNumber(2, "SizeMax", opts.h_max)
    fld.setNumber(2, "DistMin", 0.0)
    fld.setNumber(2, "DistMax", opts.grading_distance)
    fld.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)


def _read_mesh(gmsh, surf: int | None):
    """Nodos (todos) y triángulos (de ``surf`` o de todas las superficies) como índices."""
    tags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = coords.reshape(-1, 3)
    tag_to_idx = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
    tag_to_idx[tags.astype(np.int64)] = np.arange(tags.size)
    if surf is None:
        _, tri_nodes = gmsh.model.mesh.getElementsByType(2)
    else:
        _, tri_nodes = gmsh.model.mesh.getElementsByType(2, tag=surf)
    tris = tag_to_idx[tri_nodes.astype(np.int64)].reshape(-1, 3)
    return xyz, tris, tag_to_idx


def _finish(xyz, tris, raw_edges: dict | None, geo: Geometry, region=None) -> Mesh:
    """Compacta nodos, orienta triángulos (área > 0) y aristas, y valida."""
    used = np.unique(tris)
    remap = np.full(xyz.shape[0], -1, dtype=np.int64)
    remap[used] = np.arange(used.size)
    nodes = xyz[used, :2].copy()
    tris = remap[tris]

    p = nodes[tris]
    area2 = (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - (
        p[:, 2, 0] - p[:, 0, 0]
    ) * (p[:, 1, 1] - p[:, 0, 1])
    flip = area2 < 0
    tris[flip] = tris[flip][:, [0, 2, 1]]

    if raw_edges is None:
        edges = _classify_boundary(nodes, tris, geo)
    else:
        edges = {k: orient_edges(tris, remap[v], nodes.shape[0]) for k, v in raw_edges.items()}
    mesh = Mesh(
        nodes=nodes,
        tris=tris,
        region=np.zeros(tris.shape[0], dtype=np.int64) if region is None else region,
        region_names=geo.region_names,
        edges=edges,
    )
    if region is None and geo.regions:
        mesh.region = _classify_regions(mesh, geo)
    mesh.validate(BOUNDARY_GROUPS)
    return mesh


def _classify_boundary(nodes, tris, geo: Geometry) -> dict[str, np.ndarray]:
    """Asigna cada arista de borde al lado del perfil que contiene sus dos extremos."""
    bnd = boundary_edges(tris)
    a, b = nodes[bnd[:, 0]], nodes[bnd[:, 1]]
    scale = float(np.ptp(nodes[:, 0])) or 1.0
    tol = 1e-9 * scale
    owner = np.full(bnd.shape[0], -1)
    segs = geo.boundary_segments()
    for k, (p, q, _) in enumerate(segs):
        on = (point_segment_distance(a, p, q) <= tol) & (point_segment_distance(b, p, q) <= tol)
        owner[(owner < 0) & on] = k
    if np.any(owner < 0):
        raise ValueError(f"{int((owner < 0).sum())} aristas de borde no caen sobre ningún lado "
                         "del perfil (¿una región sale del dominio?).")
    edges = {}
    for g in BOUNDARY_GROUPS:
        ks = [k for k, (_, _, gg) in enumerate(segs) if gg == g]
        edges[g] = bnd[np.isin(owner, ks)].reshape(-1, 2)
    return edges


def _classify_regions(mesh: Mesh, geo: Geometry) -> np.ndarray:
    """Región de cada triángulo por su centroide; verifica áreas por región."""
    c = mesh.nodes[mesh.tris].mean(axis=1)
    region = np.zeros(mesh.n_elements, dtype=np.int64)
    A = mesh.areas()
    for j, (name, poly) in enumerate(geo.regions.items(), start=1):
        inside = points_in_polygon(c, poly)
        if np.any(inside & (region != 0)):
            raise ValueError(f"La región '{name}' se solapa con otra región.")
        region[inside] = j
        exact = polygon_area(poly)
        if abs(A[inside].sum() - exact) > 1e-9 * exact:
            raise ValueError(
                f"El área mallada de la región '{name}' ({A[inside].sum():.6e}) no coincide con "
                f"la del polígono ({exact:.6e}): ¿se solapa con otra región o sale del dominio?"
            )
    return region


def _mesh_geo_kernel(gmsh, geo: Geometry, opts: MeshOptions) -> Mesh:
    """Camino sin regiones (kernel ``geo``). Idéntico a las fases 1–6 para el rectángulo."""
    geo_k = gmsh.model.geo
    pts_rz, segs = geo.points, geo.segments
    n = len(pts_rz)
    pts = [geo_k.addPoint(r, z, 0.0, opts.h_max) for r, z in pts_rz]
    lines = [geo_k.addLine(pts[k], pts[(k + 1) % n]) for k in range(n)]
    surf = geo_k.addPlaneSurface([geo_k.addCurveLoop(lines)])
    geo_k.synchronize()

    for tag, name in enumerate(BOUNDARY_GROUPS, start=1):
        group_lines = [lines[k] for k in range(n) if segs[k] == name]
        if group_lines:
            gmsh.model.addPhysicalGroup(1, group_lines, tag, name)
    gmsh.model.addPhysicalGroup(2, [surf], 100, "disc")
    _set_size_field(gmsh, [lines[k] for k in range(n) if segs[k] == "brake"], geo, opts)

    gmsh.model.mesh.generate(2)
    xyz, tris, tag_to_idx = _read_mesh(gmsh, surf)
    raw_edges = {}
    for name in BOUNDARY_GROUPS:
        chunks = [gmsh.model.mesh.getElementsByType(1, tag=lines[k])[1]
                  for k in range(n) if segs[k] == name]
        nodes_ = np.concatenate(chunks) if chunks else np.zeros(0)
        raw_edges[name] = tag_to_idx[nodes_.astype(np.int64)].reshape(-1, 2)
    return _finish(xyz, tris, raw_edges, geo)


def _occ_polygon(occ, pts) -> int:
    p = [occ.addPoint(r, z, 0.0) for r, z in pts]
    n = len(p)
    lines = [occ.addLine(p[k], p[(k + 1) % n]) for k in range(n)]
    return occ.addPlaneSurface([occ.addCurveLoop(lines)])


def _mesh_occ_kernel(gmsh, geo: Geometry, opts: MeshOptions) -> Mesh:
    """Camino con regiones (kernel OCC + fragment; clasificación geométrica)."""
    occ = gmsh.model.occ
    outer = _occ_polygon(occ, geo.points)
    tools = [(2, _occ_polygon(occ, poly)) for poly in geo.regions.values()]
    occ.fragment([(2, outer)], tools)
    occ.synchronize()

    # Curvas de la pista: ambos extremos sobre algún lado 'brake'.
    scale = float(np.ptp(np.asarray(geo.points)[:, 0])) or 1.0
    brake_segs = [(a, b) for a, b, g in geo.boundary_segments() if g == "brake"]
    brake_curves = []
    for _, c in gmsh.model.getEntities(1):
        ends = gmsh.model.getBoundary([(1, c)], oriented=False)
        xyz = np.array([gmsh.model.getValue(0, abs(t), []) for _, t in ends])[:, :2]
        if any(np.all(point_segment_distance(xyz, a, b) <= 1e-9 * scale) for a, b in brake_segs):
            brake_curves.append(c)
    if not brake_curves:
        raise ValueError("No se encontró la curva 'brake' después de fragmentar las regiones.")
    _set_size_field(gmsh, brake_curves, geo, opts)

    gmsh.model.mesh.generate(2)
    xyz, tris, _ = _read_mesh(gmsh, None)
    return _finish(xyz, tris, None, geo)


def gmsh_mesh(geo: Geometry, opts: MeshOptions) -> Mesh:
    """Malla triangular no estructurada con refinamiento junto a ``brake``.

    El tamaño de elemento es h(d) = h_min + (h_max − h_min)·min(d/grading_distance, 1), con d
    la distancia a las curvas ``brake`` (campo ``Threshold`` con DistMin = 0). Los extremos de
    la pista son vértices de la geometría, así que hay nodos exactamente en r_a y r_b.
    """
    import gmsh

    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)  # malla reproducible
        gmsh.model.add("disc_section")
        if geo.regions:
            return _mesh_occ_kernel(gmsh, geo, opts)
        return _mesh_geo_kernel(gmsh, geo, opts)
    finally:
        gmsh.finalize()
