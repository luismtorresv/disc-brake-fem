"""Ensamblaje vectorizado de M, K, H y F para elementos P1 (PLAN.md §3.2–3.3).

Convenciones (CLAUDE.md):
- El factor 2π de la revolución se OMITE en todas las matrices y vectores.
- El peso de las integrales es r (``axisymmetric=True``) o 1 (modo plano, solo validación).
- Condición en ``brake``: −k ∇T·n = −q, el calor ENTRA; aporta +∫ q φ_i w ds a F.

Nada de bucles de Python sobre elementos: todo se calcula como arreglos (E, 3, 3) o (n, 2, 2)
y se ensambla con ``scipy.sparse.coo_matrix`` → CSR (las entradas repetidas se suman).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping

import numpy as np
import scipy.sparse as sp

from ..mesh import Mesh

# Tensor C_ijk = ∫ φ_i φ_j φ_k dA / A = (1 + δ_ij)(1 + δ_ik + δ_jk) / 60.
# Con r = Σ_k r_k φ_k, M^e_ij = ρc A Σ_k C_ijk r_k reproduce exactamente las fórmulas de la §3.3:
#   M_ii = ρc A [r_i/10 + (r_j + r_l)/30],  M_ij = ρc A [(r_i + r_j)/30 + r_l/60].
_I3 = np.eye(3)
_C3 = (1 + _I3[:, :, None]) * (1 + _I3[:, None, :] + _I3[None, :, :]) / 60.0

# Gauss–Legendre de 3 puntos en [0, 1] (exacto hasta grado 5).
_GL_S = np.array([0.5 - 0.5 * np.sqrt(0.6), 0.5, 0.5 + 0.5 * np.sqrt(0.6)])
_GL_W = np.array([5.0, 8.0, 5.0]) / 18.0


def element_geometry(p: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Coeficientes de las funciones de forma P1 de cada triángulo.

    Args:
        p: (E, 3, 2) coordenadas (r, z) de los nodos 1, 2, 3 de cada triángulo (antihorario).

    Returns:
        b, c: (E, 3) con b1 = z2 − z3, b2 = z3 − z1, b3 = z1 − z2 y c1 = r3 − r2, c2 = r1 − r3,
            c3 = r2 − r1; así ∇φ_i = (b_i, c_i) / (2A).
        A: (E,) área, A = ½[(r2 − r1)(z3 − z1) − (r3 − r1)(z2 − z1)].
    """
    r, z = p[..., 0], p[..., 1]
    b = np.stack([z[:, 1] - z[:, 2], z[:, 2] - z[:, 0], z[:, 0] - z[:, 1]], axis=1)
    c = np.stack([r[:, 2] - r[:, 1], r[:, 0] - r[:, 2], r[:, 1] - r[:, 0]], axis=1)
    A = 0.5 * ((r[:, 1] - r[:, 0]) * (z[:, 2] - z[:, 0]) - (r[:, 2] - r[:, 0]) * (z[:, 1] - z[:, 0]))
    return b, c, A


def element_matrices(p: np.ndarray, axisymmetric: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Matrices elementales de rigidez y masa por unidad de k y de ρc (fórmulas cerradas §3.3).

    K^e_ij = r̄ (b_i b_j + c_i c_j) / (4A)          (r̄ → 1 en modo plano)
    M^e_ij = A Σ_k C_ijk r_k                       (r_k → 1 en modo plano: A(1 + δ_ij)/12)

    Args:
        p: (E, 3, 2) coordenadas de los nodos de cada triángulo.
        axisymmetric: True ⇒ peso r; False ⇒ peso 1.

    Returns:
        Ke, Me: (E, 3, 3) cada una.
    """
    b, c, A = element_geometry(p)
    if np.any(A <= 0):
        raise ValueError("Hay triángulos con área ≤ 0: la malla debe ser antihoraria.")
    rk = p[..., 0] if axisymmetric else np.ones_like(A)[:, None].repeat(3, axis=1)
    rbar = rk.mean(axis=1)
    Ke = (rbar / (4.0 * A))[:, None, None] * (
        b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]
    )
    Me = A[:, None, None] * np.einsum("ijk,ek->eij", _C3, rk)
    return Ke, Me


def per_element(value: float | np.ndarray | Mapping[str, float], mesh: Mesh) -> np.ndarray:
    """Convierte una propiedad a un arreglo (E,) por triángulo.

    Acepta un escalar (homogéneo), un arreglo (E,) ya por elemento, o un diccionario
    ``{nombre_de_región: valor}`` que debe cubrir todas las regiones de la malla.
    """
    if isinstance(value, Mapping):
        missing = [n for n in mesh.region_names if n not in value]
        if missing:
            raise ValueError(f"Falta el valor de la propiedad para las regiones {missing}.")
        table = np.array([float(value[n]) for n in mesh.region_names])
        return table[mesh.region]
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return np.full(mesh.n_elements, float(arr))
    if arr.shape != (mesh.n_elements,):
        raise ValueError(
            f"La propiedad por elemento debe tener forma ({mesh.n_elements},); tiene {arr.shape}."
        )
    return arr


def _scatter(rows_idx: np.ndarray, vals: np.ndarray, n: int) -> sp.csr_matrix:
    """Ensambla bloques elementales (m, q, q) con conectividad (m, q) en una CSR n × n."""
    q = rows_idx.shape[1]
    rows = np.repeat(rows_idx, q, axis=1).ravel()
    cols = np.tile(rows_idx, (1, q)).ravel()
    return sp.coo_matrix((vals.ravel(), (rows, cols)), shape=(n, n)).tocsr()


def assemble_stiffness_mass(
    mesh: Mesh,
    k: float | np.ndarray | Mapping[str, float],
    rho_c: float | np.ndarray | Mapping[str, float],
    axisymmetric: bool = True,
) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    """Ensambla K_ij = ∫ k ∇φ_i·∇φ_j w dA y M_ij = ∫ ρc φ_i φ_j w dA (w = r o 1, sin 2π).

    ``k`` y ``rho_c`` (= ρ·c_p) son constantes por triángulo: escalar, arreglo (E,) o
    diccionario por región (ver :func:`per_element`).
    """
    Ke, Me = element_matrices(mesh.nodes[mesh.tris], axisymmetric)
    k_e = per_element(k, mesh)[:, None, None]
    rc_e = per_element(rho_c, mesh)[:, None, None]
    K = _scatter(mesh.tris, k_e * Ke, mesh.n_nodes)
    M = _scatter(mesh.tris, rc_e * Me, mesh.n_nodes)
    return K, M


def edge_matrices(p: np.ndarray, axisymmetric: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Matriz de masa de arista por unidad de h (fórmula cerrada §3.3).

    H^e = L [[r_a/4 + r_b/12, (r_a + r_b)/12], [(r_a + r_b)/12, r_a/12 + r_b/4]]
    (modo plano: L [[1/3, 1/6], [1/6, 1/3]]).

    Args:
        p: (n, 2, 2) coordenadas de los extremos a, b de cada arista.

    Returns:
        He: (n, 2, 2); L: (n,) longitudes.
    """
    L = np.linalg.norm(p[:, 1] - p[:, 0], axis=1)
    if axisymmetric:
        ra, rb = p[:, 0, 0], p[:, 1, 0]
        off = (ra + rb) / 12.0
        He = np.stack(
            [np.stack([ra / 4 + rb / 12, off], -1), np.stack([off, ra / 12 + rb / 4], -1)], 1
        )
    else:
        He = np.broadcast_to(np.array([[1 / 3, 1 / 6], [1 / 6, 1 / 3]]), (L.size, 2, 2))
    return L[:, None, None] * He, L


def assemble_robin(
    mesh: Mesh, h: Mapping[str, float], axisymmetric: bool = True
) -> sp.csr_matrix:
    """Ensambla H_ij = Σ_grupos ∫ h φ_i φ_j w ds sobre los grupos de ``h`` ({grupo: h})."""
    blocks, conn = [], []
    for group, hg in h.items():
        e = mesh.edges[group]
        He, _ = edge_matrices(mesh.nodes[e], axisymmetric)
        blocks.append(float(hg) * He)
        conn.append(e)
    if not blocks:
        return sp.csr_matrix((mesh.n_nodes, mesh.n_nodes))
    return _scatter(np.concatenate(conn), np.concatenate(blocks), mesh.n_nodes)


def assemble_edge_load(
    mesh: Mesh,
    groups: Iterable[str],
    q: Callable[[np.ndarray, np.ndarray], np.ndarray] | float,
    axisymmetric: bool = True,
) -> np.ndarray:
    """Vector F_i = Σ_grupos ∫ q(r, z) φ_i w ds con Gauss–Legendre de 3 puntos por arista.

    ``q`` es un flujo ENTRANTE [W/m²] (escalar o función vectorizada q(r, z)). Con w = r, φ
    lineal y q ∝ r, el integrando es de grado 3 y la regla es exacta.
    """
    e = np.concatenate([mesh.edges[g] for g in groups])
    pa, pb = mesh.nodes[e[:, 0]], mesh.nodes[e[:, 1]]
    L = np.linalg.norm(pb - pa, axis=1)
    # Puntos de Gauss: x(s) = (1 − s) x_a + s x_b; φ_a = 1 − s, φ_b = s.
    s = _GL_S[None, :]
    rq = (1 - s) * pa[:, [0]] + s * pb[:, [0]]
    zq = (1 - s) * pa[:, [1]] + s * pb[:, [1]]
    qv = np.broadcast_to(q(rq, zq) if callable(q) else float(q), rq.shape)
    w = rq if axisymmetric else 1.0
    g = qv * w * _GL_W[None, :] * L[:, None]
    fa = np.sum(g * (1 - s), axis=1)
    fb = np.sum(g * s, axis=1)
    F = np.zeros(mesh.n_nodes)
    np.add.at(F, e[:, 0], fa)
    np.add.at(F, e[:, 1], fb)
    return F
