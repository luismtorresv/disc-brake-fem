"""Dataclass `Mesh` (contrato entre los generadores de malla y el solver, PLAN.md §2) y
utilidades vectorizadas sobre aristas de borde."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

_MESH_KEYS = {"backend", "refine_factor", "size_track", "size_far", "grading_distance"}


@dataclass(frozen=True)
class MeshOptions:
    """Parámetros de resolución de la malla (sección ``mesh`` del YAML).

    Atributos:
        backend: ``"gmsh"`` o ``"structured"``.
        refine_factor: divide TODOS los tamaños; 2 ⇒ malla con ~4 veces más celdas.
        size_track: tamaño de elemento junto a ``brake`` [m].
        size_far: tamaño de elemento lejos de ``brake`` [m].
        grading_distance: distancia a ``brake`` en la que el tamaño crece linealmente de
            ``size_track`` a ``size_far`` [m]. No se divide por ``refine_factor``.
    """

    backend: str = "gmsh"
    refine_factor: float = 1.0
    size_track: float = 3.0e-4
    size_far: float = 4.0e-3
    grading_distance: float = 0.01

    def __post_init__(self) -> None:
        if self.backend not in ("gmsh", "structured"):
            raise ValueError(
                f"mesh.backend debe ser 'gmsh' o 'structured'; se recibió '{self.backend}'."
            )
        for name in ("refine_factor", "size_track", "size_far", "grading_distance"):
            if not getattr(self, name) > 0:
                raise ValueError(f"mesh.{name} debe ser > 0; se recibió {getattr(self, name)}.")
        if self.size_track > self.size_far:
            raise ValueError(
                "mesh.size_track no puede ser mayor que mesh.size_far "
                f"({self.size_track} > {self.size_far})."
            )

    @property
    def h_min(self) -> float:
        """Tamaño junto a la pista, ya dividido por ``refine_factor``."""
        return self.size_track / self.refine_factor

    @property
    def h_max(self) -> float:
        """Tamaño lejano, ya dividido por ``refine_factor``."""
        return self.size_far / self.refine_factor

    @classmethod
    def from_dict(cls, d: dict) -> MeshOptions:
        """Construye las opciones a partir de la sección ``mesh`` del YAML."""
        unknown = set(d) - _MESH_KEYS
        if unknown:
            raise ValueError(f"Claves desconocidas en 'mesh': {sorted(unknown)}.")
        kw = {k: (v if k == "backend" else float(v)) for k, v in d.items()}
        return cls(**kw)


@dataclass
class Mesh:
    """Malla triangular P1 de la sección (r, z).

    Atributos:
        nodes: (N, 2) coordenadas (r, z).
        tris: (E, 3) conectividad, antihoraria (área > 0).
        region: (E,) índice de región (material) de cada triángulo.
        region_names: nombre de cada región; ``region_names[region[e]]``.
        edges: por grupo de frontera, aristas (n, 2) orientadas en sentido antihorario del
            borde (el dominio queda a la izquierda; la normal exterior es (t_z, −t_r)).
    """

    nodes: np.ndarray
    tris: np.ndarray
    region: np.ndarray
    region_names: list[str]
    edges: dict[str, np.ndarray] = field(default_factory=dict)

    @property
    def n_nodes(self) -> int:
        return self.nodes.shape[0]

    @property
    def n_elements(self) -> int:
        return self.tris.shape[0]

    def areas(self) -> np.ndarray:
        """Área con signo de cada triángulo, A = ½[(r2−r1)(z3−z1) − (r3−r1)(z2−z1)]."""
        p = self.nodes[self.tris]
        d1 = p[:, 1] - p[:, 0]
        d2 = p[:, 2] - p[:, 0]
        return 0.5 * (d1[:, 0] * d2[:, 1] - d2[:, 0] * d1[:, 1])

    def r_centroid(self) -> np.ndarray:
        """Radio del centroide r̄ = (r1 + r2 + r3)/3 de cada triángulo."""
        return self.nodes[self.tris, 0].mean(axis=1)

    def validate(self, groups: tuple[str, ...]) -> None:
        """Comprueba el contrato: triángulos antihorarios, nodos usados, y que los grupos
        ``groups`` particionen exactamente las aristas de borde. Lanza ``ValueError``."""
        if self.nodes.ndim != 2 or self.nodes.shape[1] != 2:
            raise ValueError("Mesh.nodes debe tener forma (N, 2).")
        if self.tris.ndim != 2 or self.tris.shape[1] != 3:
            raise ValueError("Mesh.tris debe tener forma (E, 3).")
        if self.region.shape != (self.n_elements,):
            raise ValueError("Mesh.region debe tener un valor por triángulo.")
        if np.any(self.areas() <= 0):
            raise ValueError("Hay triángulos con área ≤ 0 (orientación horaria o degenerados).")
        if np.setdiff1d(np.arange(self.n_nodes), self.tris.ravel()).size:
            raise ValueError("Hay nodos que no pertenecen a ningún triángulo.")
        # Un grupo puede estar vacío (p. ej. T7: `brake` cubre toda la cara superior y
        # top_inner/top_outer no tienen aristas), pero la clave debe existir.
        missing = [g for g in groups if g not in self.edges]
        if missing:
            raise ValueError(f"Grupos de frontera ausentes: {missing}.")
        extra = set(self.edges) - set(groups)
        if extra:
            raise ValueError(f"Grupos de frontera no reconocidos: {sorted(extra)}.")
        bnd = boundary_edges(self.tris)
        all_keys = np.concatenate(
            [_edge_keys(np.reshape(self.edges[g], (-1, 2)), self.n_nodes) for g in groups]
        )
        if np.unique(all_keys).size != all_keys.size:
            raise ValueError("Una arista de borde aparece en más de un grupo (o repetida).")
        if not np.array_equal(np.sort(all_keys), np.sort(_edge_keys(bnd, self.n_nodes))):
            raise ValueError("Los grupos de frontera no cubren exactamente el borde de la malla.")


def _edge_keys(edges: np.ndarray, n_nodes: int) -> np.ndarray:
    """Clave entera única de cada arista, independiente de su orientación."""
    e = np.sort(np.asarray(edges, dtype=np.int64), axis=1)
    return e[:, 0] * n_nodes + e[:, 1]


def boundary_edges(tris: np.ndarray) -> np.ndarray:
    """Aristas que pertenecen a un solo triángulo, orientadas como en ese triángulo.

    Si los triángulos son antihorarios, las aristas quedan antihorarias respecto al borde
    (dominio a la izquierda).
    """
    e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    n = int(tris.max()) + 1
    _, idx, counts = np.unique(_edge_keys(e, n), return_index=True, return_counts=True)
    return e[idx[counts == 1]]


def orient_edges(tris: np.ndarray, edges: np.ndarray, n_nodes: int) -> np.ndarray:
    """Reorienta ``edges`` como aparecen en ``boundary_edges(tris)``.

    Lanza ``ValueError`` si alguna arista no es de borde.
    """
    bnd = boundary_edges(tris)
    bkeys = _edge_keys(bnd, n_nodes)
    order = np.argsort(bkeys)
    keys = _edge_keys(edges, n_nodes)
    pos = np.searchsorted(bkeys, keys, sorter=order)
    pos = np.clip(pos, 0, bkeys.size - 1)
    found = bkeys[order[pos]] == keys
    if not np.all(found):
        raise ValueError(f"{int((~found).sum())} aristas del grupo no están en el borde.")
    return bnd[order[pos]]
