"""Pruebas T2 (PLAN.md §6): matrices elementales y ensamblaje.

Las referencias se calculan con cuadratura de alto orden INDEPENDIENTE de las fórmulas cerradas:
- triángulos: Gauss colapsado (Duffy) 8×8, exacto hasta grado 14; gradientes de las funciones
  de forma obtenidos invirtiendo el jacobiano de la transformación afín (no con b_i, c_i);
- aristas: Gauss–Legendre de 10 puntos.
Esto protege contra errores en los coeficientes 1/10, 1/30, 1/60, 1/4, 1/12.
"""

from __future__ import annotations

import numpy as np
import pytest

from discfem.fem import (
    assemble_edge_load,
    assemble_robin,
    assemble_stiffness_mass,
    edge_matrices,
    element_matrices,
    per_element,
)
from discfem.geometry import RectangleGeometry
from discfem.mesh import MeshOptions, build_mesh

GEO = RectangleGeometry(r_inner=0.07, r_outer=0.15, half_thickness=0.01, r_a=0.08, r_b=0.14)
RTOL = 1e-12


def random_triangles(n: int, seed: int = 0, r_min: float = 0.0) -> np.ndarray:
    """Triángulos aleatorios antihorarios, no degenerados, en r ∈ [r_min, 0.2], z ∈ [0, 0.02]."""
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n:
        p = np.column_stack([rng.uniform(r_min, 0.2, 3), rng.uniform(0.0, 0.02, 3)])
        a2 = (p[1, 0] - p[0, 0]) * (p[2, 1] - p[0, 1]) - (p[2, 0] - p[0, 0]) * (p[1, 1] - p[0, 1])
        if abs(a2) < 1e-6:  # evita triángulos casi degenerados
            continue
        out.append(p if a2 > 0 else p[[0, 2, 1]])
    return np.array(out)


def quad_triangle_reference(p: np.ndarray, axisymmetric: bool, n: int = 8):
    """K y M de un triángulo con cuadratura de Duffy n×n e inversión del jacobiano."""
    x, w = np.polynomial.legendre.leggauss(n)
    x, w = 0.5 * (x + 1), 0.5 * w
    U, V = np.meshgrid(x, x, indexing="ij")
    WU, WV = np.meshgrid(w, w, indexing="ij")
    xi, eta = U, V * (1 - U)
    wq = (WU * WV * (1 - U)).ravel()
    lam = np.stack([1 - xi - eta, xi, eta]).reshape(3, -1)  # (3, nq)
    J = np.array([p[1] - p[0], p[2] - p[0]]).T  # columnas: dx/dξ, dx/dη
    detJ = np.linalg.det(J)
    G = np.linalg.inv(J)  # filas: ∇ξ, ∇η
    grads = np.array([-G[0] - G[1], G[0], G[1]])  # ∇λ_i
    r = lam.T @ p[:, 0]
    wt = (r if axisymmetric else np.ones_like(r)) * wq * detJ
    M = (lam * wt) @ lam.T
    K = (grads @ grads.T) * wt.sum()
    return K, M


@pytest.mark.parametrize("axisymmetric", [True, False])
def test_element_matrices_vs_quadrature(axisymmetric):
    P = random_triangles(200, seed=1)
    Ke, Me = element_matrices(P, axisymmetric)
    for k in range(P.shape[0]):
        Kq, Mq = quad_triangle_reference(P[k], axisymmetric)
        np.testing.assert_allclose(Ke[k], Kq, rtol=0, atol=RTOL * np.abs(Kq).max())
        np.testing.assert_allclose(Me[k], Mq, rtol=0, atol=RTOL * np.abs(Mq).max())


def test_element_mass_closed_formulas_from_plan():
    """Comprueba literalmente M_ii = A[r_i/10 + (r_j+r_l)/30], M_ij = A[(r_i+r_j)/30 + r_l/60]."""
    P = random_triangles(50, seed=2)
    _, Me = element_matrices(P, True)
    A = 0.5 * ((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1])
               - (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1]))
    r = P[..., 0]
    for i in range(3):
        j, l = [m for m in range(3) if m != i]
        np.testing.assert_allclose(Me[:, i, i], A * (r[:, i] / 10 + (r[:, j] + r[:, l]) / 30),
                                   rtol=RTOL)
        np.testing.assert_allclose(Me[:, i, j], A * ((r[:, i] + r[:, j]) / 30 + r[:, l] / 60),
                                   rtol=RTOL)


def test_element_matrices_with_vertex_on_axis():
    """r_i = 0 (eje de simetría, T6): fórmulas siguen exactas, sin singularidad."""
    P = random_triangles(50, seed=3)
    P[:, 0, 0] = 0.0
    # Mover un vértice puede invertir la orientación: se reordenan los horarios y se descartan
    # los casi degenerados.
    a2 = ((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1])
          - (P[:, 2, 0] - P[:, 0, 0]) * (P[:, 1, 1] - P[:, 0, 1]))
    P[a2 < 0] = P[a2 < 0][:, [0, 2, 1]]
    P = P[np.abs(a2) > 1e-6]
    assert P.shape[0] > 20
    Ke, Me = element_matrices(P, True)
    for k in range(P.shape[0]):
        Kq, Mq = quad_triangle_reference(P[k], True)
        np.testing.assert_allclose(Ke[k], Kq, rtol=0, atol=RTOL * np.abs(Kq).max())
        np.testing.assert_allclose(Me[k], Mq, rtol=0, atol=RTOL * np.abs(Mq).max())


@pytest.mark.parametrize("axisymmetric", [True, False])
def test_edge_matrices_vs_quadrature(axisymmetric):
    rng = np.random.default_rng(4)
    P = np.stack([np.column_stack([rng.uniform(0, 0.2, 2), rng.uniform(0, 0.02, 2)])
                  for _ in range(100)])
    He, L = edge_matrices(P, axisymmetric)
    x, w = np.polynomial.legendre.leggauss(10)
    s, w = 0.5 * (x + 1), 0.5 * w
    phi = np.stack([1 - s, s])
    for k in range(P.shape[0]):
        r = (1 - s) * P[k, 0, 0] + s * P[k, 1, 0]
        wt = (r if axisymmetric else 1.0) * w * L[k]
        Hq = (phi * wt) @ phi.T
        np.testing.assert_allclose(He[k], Hq, rtol=0, atol=RTOL * np.abs(Hq).max())


# --------------------------------------------------------------------- ensamblaje global

K_DISC, RHO, CP = 50.0, 7200.0, 460.0


@pytest.fixture(scope="module", params=["structured", "gmsh"])
def mesh(request):
    opts = MeshOptions(backend=request.param, size_track=1e-3, size_far=4e-3)
    return build_mesh(GEO, opts)


@pytest.mark.parametrize("axisymmetric", [True, False])
def test_global_K_M_properties(mesh, axisymmetric):
    K, M = assemble_stiffness_mass(mesh, K_DISC, RHO * CP, axisymmetric)
    one = np.ones(mesh.n_nodes)
    # K·1 = 0 (las funciones de forma suman 1 ⇒ el gradiente de una constante es 0).
    assert np.abs(K @ one).max() < RTOL * abs(K).max()
    # Simetría de K y M.
    assert abs(K - K.T).max() < RTOL * abs(K).max()
    assert abs(M - M.T).max() < RTOL * abs(M).max()
    # 1ᵀM1 = ρc ∫ w dA, con w = r o 1.
    A = mesh.areas()
    exact = RHO * CP * np.sum(A * (mesh.r_centroid() if axisymmetric else 1.0))
    assert abs(one @ (M @ one) - exact) / exact < RTOL


def test_mass_positive_definite_and_K_semidefinite():
    m = build_mesh(GEO, MeshOptions(backend="structured", size_track=2e-3, size_far=4e-3))
    K, M = assemble_stiffness_mass(m, K_DISC, RHO * CP, True)
    ev_m = np.linalg.eigvalsh(M.toarray())
    ev_k = np.linalg.eigvalsh(K.toarray())
    assert ev_m.min() > 0
    assert ev_k.min() > -RTOL * ev_k.max()
    assert np.sum(ev_k < RTOL * ev_k.max()) == 1  # un solo modo nulo: la constante


def test_robin_matrix_row_sums_equal_weighted_length(mesh):
    """Σ_j H_ij sumado sobre i = h ∫ r ds en los grupos con convección (prueba global de H)."""
    groups = {"top_inner": 40.0, "top_outer": 40.0, "r_inner": 40.0, "r_outer": 40.0}
    H = assemble_robin(mesh, groups, True)
    e, ri, ro, ra, rb = GEO.half_thickness, GEO.r_inner, GEO.r_outer, GEO.r_a, GEO.r_b
    exact = 40.0 * ((ra**2 - ri**2) / 2 + (ro**2 - rb**2) / 2 + ri * e + ro * e)
    total = np.ones(mesh.n_nodes) @ (H @ np.ones(mesh.n_nodes))
    assert abs(total - exact) / exact < RTOL
    assert abs(H - H.T).max() < RTOL * abs(H).max()


def test_edge_load_exact_for_uniform_pressure_profile(mesh):
    """q ∝ r (presión uniforme): Σ F = ∫_{r_a}^{r_b} c r · r dr = c (r_b³ − r_a³)/3, exacto."""
    c = 7.7e6
    F = assemble_edge_load(mesh, ["brake"], lambda r, z: c * r, True)
    exact = c * (GEO.r_b**3 - GEO.r_a**3) / 3
    assert abs(F.sum() - exact) / exact < RTOL
    # Solo los nodos de brake reciben carga, y es positiva (el calor entra).
    brake_nodes = np.unique(mesh.edges["brake"])
    assert np.all(F[brake_nodes] > 0)
    assert np.all(F[np.setdiff1d(np.arange(mesh.n_nodes), brake_nodes)] == 0)


def test_edge_load_matches_per_edge_quadrature(mesh):
    """Comparación nodo a nodo con Gauss de 10 puntos para q ∝ r (polinomio, exacto en ambas)."""
    c = 3.0e6
    F = assemble_edge_load(mesh, ["brake"], lambda r, z: c * r, True)
    x, w = np.polynomial.legendre.leggauss(10)
    s, w = 0.5 * (x + 1), 0.5 * w
    Fq = np.zeros(mesh.n_nodes)
    for a, b in mesh.edges["brake"]:
        pa, pb = mesh.nodes[a], mesh.nodes[b]
        r = (1 - s) * pa[0] + s * pb[0]
        g = c * r * r * w * np.linalg.norm(pb - pa)
        Fq[a] += np.sum(g * (1 - s))
        Fq[b] += np.sum(g * s)
    np.testing.assert_allclose(F, Fq, rtol=0, atol=RTOL * np.abs(Fq).max())


def test_edge_load_constant_and_convection_consistency(mesh):
    """F_conv = ∫ h T∞ φ_i r ds coincide con H·(T∞ 1), porque Σ_j φ_j = 1."""
    groups = ["top_inner", "top_outer", "r_inner", "r_outer"]
    h, T_inf = 40.0, 25.0
    H = assemble_robin(mesh, dict.fromkeys(groups, h), True)
    F = assemble_edge_load(mesh, groups, h * T_inf, True)
    np.testing.assert_allclose(F, H @ np.full(mesh.n_nodes, T_inf), rtol=0,
                               atol=RTOL * np.abs(F).max())


def test_properties_per_region():
    """Dos regiones con el mismo material ⇒ mismas matrices que el caso homogéneo; con
    materiales distintos, K es lineal en k por región."""
    m = build_mesh(GEO, MeshOptions(backend="structured", size_track=2e-3, size_far=4e-3))
    m.region = (m.r_centroid() > 0.11).astype(np.int64)
    m.region_names = ["inner", "outer"]
    K0, M0 = assemble_stiffness_mass(m, K_DISC, RHO * CP)
    K1, M1 = assemble_stiffness_mass(m, {"inner": K_DISC, "outer": K_DISC},
                                     {"inner": RHO * CP, "outer": RHO * CP})
    assert abs(K0 - K1).max() == 0 and abs(M0 - M1).max() == 0
    Ka, _ = assemble_stiffness_mass(m, (m.region == 0).astype(float), 1.0)
    Kb, _ = assemble_stiffness_mass(m, (m.region == 1).astype(float), 1.0)
    K2, _ = assemble_stiffness_mass(m, {"inner": 50.0, "outer": 20.0}, 1.0)
    assert abs(K2 - (50.0 * Ka + 20.0 * Kb)).max() < RTOL * abs(K2).max()


def test_per_element_errors():
    m = build_mesh(GEO, MeshOptions(backend="structured", size_track=2e-3, size_far=4e-3))
    with pytest.raises(ValueError, match="regiones"):
        per_element({"otra": 1.0}, m)
    with pytest.raises(ValueError, match="forma"):
        per_element(np.ones(3), m)
