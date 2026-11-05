"""Generador estructurado de respaldo (solo numpy).

Rejilla tensorial r × z partida en triángulos. En r se construye por tramos [r_i, r_a],
[r_a, r_b], [r_b, r_o] para que haya nodos exactamente en r_a y r_b; en z la graduación es
geométrica-lineal con la celda MÁS pequeña tocando ``brake`` (z = e). Ver PLAN.md §5, Fase 1, y
la trampa 1 de la §7 (graduación al revés).
"""

from __future__ import annotations

import numpy as np

from ..geometry import BOUNDARY_GROUPS, RectangleGeometry
from .base import Mesh, MeshOptions


def graded_spacing(length: float, h_min: float, h_max: float, distance: float) -> np.ndarray:
    """Distancias d_0 = 0 < d_1 < ... < d_n = ``length`` medidas desde la pista.

    El tamaño objetivo crece linealmente con la distancia d a la pista,
    h(d) = h_min + (h_max − h_min)·min(d/distance, 1), igual que el campo ``Threshold`` de
    Gmsh. Se avanza d_{k+1} = d_k + h(d_k) hasta pasar ``length`` y se escala todo por
    length/d_n ≤ 1, de modo que ninguna celda supera su tamaño objetivo y las celdas son no
    decrecientes al alejarse de la pista (la primera es la más pequeña).
    """
    d = [0.0]
    while d[-1] < length:
        x = d[-1]
        d.append(x + h_min + (h_max - h_min) * min(x / distance, 1.0))
    d = np.asarray(d)
    d = d * (length / d[-1])  # escalado uniforme: conserva la monotonía de los tamaños
    d[-1] = length
    return d


def _r_grid(geo: RectangleGeometry, opts: MeshOptions) -> tuple[np.ndarray, int, int]:
    """Coordenadas r por tramos; devuelve (r, índice de r_a, índice de r_b)."""
    h_min, h_max, dist = opts.h_min, opts.h_max, opts.grading_distance
    left = geo.r_a - graded_spacing(geo.r_a - geo.r_inner, h_min, h_max, dist)[::-1]
    n_track = max(1, int(np.ceil((geo.r_b - geo.r_a) / h_min - 1e-9)))
    track = np.linspace(geo.r_a, geo.r_b, n_track + 1)
    right = geo.r_b + graded_spacing(geo.r_outer - geo.r_b, h_min, h_max, dist)
    left[0], left[-1] = geo.r_inner, geo.r_a
    right[0], right[-1] = geo.r_b, geo.r_outer
    r = np.concatenate([left[:-1], track, right[1:]])
    ia = left.size - 1
    ib = ia + n_track
    return r, ia, ib


def _z_grid(geo: RectangleGeometry, opts: MeshOptions) -> np.ndarray:
    """Coordenadas z crecientes de 0 a e, con la celda más fina en z = e."""
    e = geo.half_thickness
    z = e - graded_spacing(e, opts.h_min, opts.h_max, opts.grading_distance)[::-1]
    z[0], z[-1] = 0.0, e
    return z


def structured_mesh(geo: RectangleGeometry, opts: MeshOptions) -> Mesh:
    """Malla estructurada triangular de la sección rectangular, graduada según ``opts``."""
    r, ia, ib = _r_grid(geo, opts)
    z = _z_grid(geo, opts)
    return tensor_mesh(r, z, ia, ib)


def tensor_mesh(r: np.ndarray, z: np.ndarray, ia: int, ib: int) -> Mesh:
    """Malla triangular de la rejilla r × z dada por coordenadas explícitas.

    Args:
        r: coordenadas radiales estrictamente crecientes (r[0] = r_i, r[-1] = r_o).
        z: coordenadas axiales estrictamente crecientes (z[0] = 0, z[-1] = e).
        ia, ib: índices en ``r`` de r_a y r_b (0 ≤ ia < ib ≤ len(r) − 1); las aristas
            superiores entre ellos forman ``brake``. Con ia = 0 queda ``top_inner`` vacío y con
            ib = len(r) − 1 queda ``top_outer`` vacío (útil en T7).

    Cada celda (i, j) con nodos n00, n10, n01, n11 se parte en los triángulos (n00, n10, n11)
    y (n00, n11, n01), ambos antihorarios. El nodo (i, j) tiene índice j·n_r + i. Útil también
    para las pruebas de validación (T6, T7) que piden mallas uniformes o con capas dadas.
    """
    r = np.asarray(r, dtype=float)
    z = np.asarray(z, dtype=float)
    if np.any(np.diff(r) <= 0) or np.any(np.diff(z) <= 0):
        raise ValueError("Las coordenadas r y z de la rejilla deben ser estrictamente crecientes.")
    if not (0 <= ia < ib <= r.size - 1):
        raise ValueError(
            f"Se requiere 0 ≤ ia < ib ≤ len(r) − 1; se recibió ia={ia}, ib={ib}, len(r)={r.size}."
        )
    nr, nz = r.size, z.size
    R, Z = np.meshgrid(r, z)
    nodes = np.column_stack([R.ravel(), Z.ravel()])

    ii, jj = np.meshgrid(np.arange(nr - 1), np.arange(nz - 1))
    n00 = (jj * nr + ii).ravel()
    n10, n01 = n00 + 1, n00 + nr
    n11 = n01 + 1
    tris = np.concatenate(
        [np.column_stack([n00, n10, n11]), np.column_stack([n00, n11, n01])]
    ).astype(np.int64)

    i = np.arange(nr - 1)
    j = np.arange(nz - 1)
    top = (nz - 1) * nr
    bottom_e = np.column_stack([i, i + 1])                       # +r en z = 0
    outer_e = np.column_stack([j * nr + nr - 1, (j + 1) * nr + nr - 1])  # +z en r = r_o
    top_e = np.column_stack([top + i + 1, top + i])              # −r en z = e
    inner_e = np.column_stack([(j + 1) * nr, j * nr])            # −z en r = r_i
    edges = {
        "symmetry": bottom_e,
        "r_outer": outer_e,
        "top_outer": top_e[ib:],
        "brake": top_e[ia:ib],
        "top_inner": top_e[:ia],
        "r_inner": inner_e,
    }
    edges = {k: np.reshape(v, (-1, 2)).astype(np.int64) for k, v in edges.items()}

    mesh = Mesh(
        nodes=nodes,
        tris=tris,
        region=np.zeros(tris.shape[0], dtype=np.int64),
        region_names=["disc"],
        edges=edges,
    )
    mesh.validate(BOUNDARY_GROUPS)
    return mesh
