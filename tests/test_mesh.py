"""Pruebas T1 (PLAN.md §6): malla, con ambos backends.

- Σ A r̄ = e (r_o² − r_i²)/2 y Σ A = e (r_o − r_i), a 1e-12 relativo.
- Cada arista de borde en exactamente un grupo; `brake` abarca exactamente [r_a, r_b].
- Triángulos antihorarios; aristas de borde con normal exterior correcta.
- La celda más pequeña toca `brake` (graduación no invertida, trampa 1 de la §7).
"""

from __future__ import annotations

import numpy as np
import pytest

from discfem.geometry import BOUNDARY_GROUPS, RectangleGeometry
from discfem.mesh import MeshOptions, boundary_edges, build_mesh
from discfem.mesh.structured import graded_spacing, structured_mesh

GEO = RectangleGeometry(r_inner=0.07, r_outer=0.15, half_thickness=0.01, r_a=0.08, r_b=0.14)
BACKENDS = ["structured", "gmsh"]

_cache: dict = {}


def get_mesh(backend: str, refine: float = 1.0):
    """Malla del caso de referencia (cacheada: Gmsh tarda ~0,1 s por malla)."""
    key = (backend, refine)
    if key not in _cache:
        _cache[key] = build_mesh(GEO, MeshOptions(backend=backend, refine_factor=refine))
    return _cache[key]


@pytest.fixture(params=BACKENDS)
def mesh(request):
    return get_mesh(request.param)


def test_weighted_area_exact(mesh):
    e, ri, ro = GEO.half_thickness, GEO.r_inner, GEO.r_outer
    A = mesh.areas()
    exact_r = e * (ro**2 - ri**2) / 2
    assert abs(np.sum(A * mesh.r_centroid()) - exact_r) / exact_r < 1e-12
    exact = e * (ro - ri)
    assert abs(A.sum() - exact) / exact < 1e-12


def test_triangles_counterclockwise(mesh):
    assert np.all(mesh.areas() > 0)


def test_groups_partition_boundary(mesh):
    n = mesh.n_nodes
    keys = {g: set(map(tuple, np.sort(mesh.edges[g], axis=1))) for g in BOUNDARY_GROUPS}
    assert set(mesh.edges) == set(BOUNDARY_GROUPS)
    assert all(len(k) == len(mesh.edges[g]) > 0 for g, k in keys.items())
    union = set().union(*keys.values())
    assert sum(len(k) for k in keys.values()) == len(union)  # disjuntos
    bnd = set(map(tuple, np.sort(boundary_edges(mesh.tris), axis=1)))
    assert union == bnd
    assert n == np.unique(mesh.tris).size


EXPECTED_SPAN = {
    # grupo: (coordenada fija, valor fijo, coordenada libre, mínimo, máximo)
    "symmetry": (1, 0.0, 0, GEO.r_inner, GEO.r_outer),
    "r_outer": (0, GEO.r_outer, 1, 0.0, GEO.half_thickness),
    "top_outer": (1, GEO.half_thickness, 0, GEO.r_b, GEO.r_outer),
    "brake": (1, GEO.half_thickness, 0, GEO.r_a, GEO.r_b),
    "top_inner": (1, GEO.half_thickness, 0, GEO.r_inner, GEO.r_a),
    "r_inner": (0, GEO.r_inner, 1, 0.0, GEO.half_thickness),
}


@pytest.mark.parametrize("group", BOUNDARY_GROUPS)
def test_group_spans_exactly(mesh, group):
    fixed, value, free, lo, hi = EXPECTED_SPAN[group]
    pts = mesh.nodes[np.unique(mesh.edges[group])]
    assert np.all(pts[:, fixed] == value)
    assert pts[:, free].min() == lo
    assert pts[:, free].max() == hi
    seg = mesh.nodes[mesh.edges[group]]
    length = np.linalg.norm(seg[:, 1] - seg[:, 0], axis=1).sum()
    assert abs(length - (hi - lo)) / (hi - lo) < 1e-12


def test_nodes_exactly_at_track_ends(mesh):
    e = GEO.half_thickness
    for r in (GEO.r_a, GEO.r_b):
        assert np.any((mesh.nodes[:, 0] == r) & (mesh.nodes[:, 1] == e))


OUTWARD = {
    "symmetry": (0, -1), "r_outer": (1, 0), "top_outer": (0, 1),
    "brake": (0, 1), "top_inner": (0, 1), "r_inner": (-1, 0),
}


@pytest.mark.parametrize("group", BOUNDARY_GROUPS)
def test_edges_oriented_with_outward_normal(mesh, group):
    seg = mesh.nodes[mesh.edges[group]]
    t = seg[:, 1] - seg[:, 0]
    normal = np.column_stack([t[:, 1], -t[:, 0]]) / np.linalg.norm(t, axis=1)[:, None]
    np.testing.assert_allclose(normal, np.broadcast_to(OUTWARD[group], normal.shape), atol=1e-12)


def _brake_distance(mesh):
    c = mesh.nodes[mesh.tris].mean(axis=1)
    dr = np.maximum(0.0, np.maximum(GEO.r_a - c[:, 0], c[:, 0] - GEO.r_b))
    return np.hypot(dr, GEO.half_thickness - c[:, 1])


def test_smallest_cell_touches_brake_structured():
    """Estructurada: el triángulo de menor área tiene un nodo en `brake`, y los espesores de
    capa en z crecen de forma monótona al alejarse de z = e."""
    m = get_mesh("structured")
    A = m.areas()
    brake_nodes = np.unique(m.edges["brake"])
    assert np.isin(m.tris[A.argmin()], brake_nodes).any()
    z = np.unique(m.nodes[:, 1])[::-1]  # desde z = e hacia z = 0
    dz = -np.diff(z)
    assert np.all(np.diff(dz) >= -1e-15)
    assert dz[0] < dz[-1]


def test_smallest_cell_near_brake_gmsh():
    """Gmsh (no estructurada): el triángulo de menor área está a menos de 2·h_min de `brake`
    (primera o segunda fila), y el área media junto a `brake` es < 1/10 de la del fondo."""
    m = get_mesh("gmsh")
    h_min = MeshOptions().h_min
    A = m.areas()
    assert _brake_distance(m)[A.argmin()] < 2 * h_min
    touch_brake = np.isin(m.tris, np.unique(m.edges["brake"])).any(axis=1)
    touch_sym = np.isin(m.tris, np.unique(m.edges["symmetry"])).any(axis=1)
    assert A[touch_brake].mean() < 0.1 * A[touch_sym].mean()


@pytest.mark.parametrize("backend", BACKENDS)
def test_refine_factor_quadruples_cells(backend):
    ratio = get_mesh(backend, 2.0).n_elements / get_mesh(backend, 1.0).n_elements
    assert 3.0 < ratio < 5.0


def test_gmsh_is_reproducible():
    a = get_mesh("gmsh")
    b = build_mesh(GEO, MeshOptions(backend="gmsh"))
    assert np.array_equal(a.nodes, b.nodes) and np.array_equal(a.tris, b.tris)


def test_graded_spacing_monotone_and_bounded():
    d = graded_spacing(0.01, 3e-4, 4e-3, 0.01)
    assert d[0] == 0.0 and d[-1] == 0.01
    h = np.diff(d)
    assert np.all(np.diff(h) >= -1e-15)
    assert h[0] <= 3e-4 * (1 + 1e-12)


def test_structured_uniform_when_sizes_equal():
    m = structured_mesh(GEO, MeshOptions(backend="structured", size_track=2e-3, size_far=2e-3))
    dz = np.diff(np.unique(m.nodes[:, 1]))
    np.testing.assert_allclose(dz, dz[0], rtol=1e-10)


@pytest.mark.parametrize(
    "bad",
    [
        # (r_inner, r_outer, half_thickness, r_a, r_b)
        (0.07, 0.15, 0.01, 0.06, 0.14),  # pista empieza antes de r_i
        (0.07, 0.15, 0.01, 0.08, 0.16),  # pista termina después de r_o
        (0.07, 0.15, 0.0, 0.08, 0.14),   # espesor nulo
    ],
)
def test_invalid_geometry_raises(bad):
    with pytest.raises(ValueError):
        RectangleGeometry(*bad)


def test_invalid_mesh_options_raise():
    with pytest.raises(ValueError, match="backend"):
        MeshOptions(backend="triangle")
    with pytest.raises(ValueError, match="size_track"):
        MeshOptions(size_track=1e-2, size_far=1e-3)
    with pytest.raises(ValueError, match="desconocidas"):
        MeshOptions.from_dict({"refine": 2})


def test_cli_mesh(tmp_path):
    import yaml

    from discfem.cli import main

    cfg = {
        "case": "prueba_cli",
        "geometry": {"r_inner": 0.07, "r_outer": 0.15, "half_thickness": 0.01,
                     "track": {"r_a": 0.08, "r_b": 0.14}},
        "mesh": {"backend": "structured", "size_track": 1e-3, "size_far": 4e-3},
        "output": {"dir": str(tmp_path / "out")},
    }
    path = tmp_path / "c.yaml"
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    assert main(["mesh", str(path)]) == 0
    assert (tmp_path / "out" / "mesh.png").stat().st_size > 0
    assert (tmp_path / "out" / "config.yaml").read_text(encoding="utf-8") == path.read_text(
        encoding="utf-8"
    )

    # Regresión: en una consola cp1252 (Windows) imprimir Σ o r̄ fallaba con
    # UnicodeEncodeError. pytest captura en UTF-8, así que se prueba en un subproceso.
    import os
    import subprocess
    import sys

    env = {**os.environ, "PYTHONIOENCODING": "cp1252"}
    proc = subprocess.run(
        [sys.executable, "-m", "discfem", "mesh", str(path)],
        check=False, env=env, capture_output=True, cwd=tmp_path,
    )
    assert proc.returncode == 0, proc.stderr.decode("utf-8", "replace")
