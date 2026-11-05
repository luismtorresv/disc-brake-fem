"""Generación de mallas de la sección (r, z): dataclass `Mesh` y dos backends."""

from __future__ import annotations

from ..geometry import RectangleGeometry
from .base import Mesh, MeshOptions, boundary_edges, orient_edges

__all__ = ["Mesh", "MeshOptions", "boundary_edges", "build_mesh", "orient_edges"]


def build_mesh(geo: RectangleGeometry, opts: MeshOptions) -> Mesh:
    """Genera la malla con el backend indicado en ``opts.backend``."""
    if opts.backend == "structured":
        from .structured import structured_mesh

        return structured_mesh(geo, opts)
    from .gmsh_backend import gmsh_mesh

    return gmsh_mesh(geo, opts)
