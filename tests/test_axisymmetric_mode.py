"""Prueba T6 (PLAN.md §6): modo propio axisimétrico de un cilindro macizo aislado.

Cilindro r ∈ [0, R], R = 0,05, semiespesor 0,01, k = 50, ρ = 7200, c_p = 460, aislado en todo el
borde. Condición inicial T0 = 1 + 0,5 J0(λ r/R), con λ el primer cero de J1 (≈ 3,8317), de modo
que ∂T/∂r = 0 en r = R. Solución exacta: T = 1 + 0,5 J0(λ r/R) exp(−α λ² t / R²), α = k/(ρ c_p).

Es la prueba que valida M y K CON PESO r: en modo plano el modo propio sería un coseno, no J0.
El avance en el tiempo (Crank–Nicolson, θ = 1/2, Δt = 0,02 s, hasta t = 10 s) se escribe aquí
directamente; el módulo `timestepping` llega en la Fase 3.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse.linalg as spla
from scipy.special import j0, jn_zeros

from discfem.fem import assemble_stiffness_mass
from discfem.mesh.structured import tensor_mesh

R, E = 0.05, 0.01
K_, RHO, CP = 50.0, 7200.0, 460.0
ALPHA = K_ / (RHO * CP)
LAM = jn_zeros(1, 1)[0]
DT, T_END, THETA = 0.02, 10.0, 0.5


def run_mode(nx: int) -> tuple[float, float]:
    """Devuelve (error máx. en t = 10 s / 0,5, error máx. en todos los pasos / 0,5)."""
    nz = nx // 5
    mesh = tensor_mesh(np.linspace(0, R, nx + 1), np.linspace(0, E, nz + 1), 1, nx - 1)
    K, M = assemble_stiffness_mass(mesh, K_, RHO * CP, axisymmetric=True)
    r = mesh.nodes[:, 0]
    shape = 0.5 * j0(LAM * r / R)
    T = 1.0 + shape
    lhs = spla.splu((M + THETA * DT * K).tocsc())
    rhs_mat = (M - (1 - THETA) * DT * K).tocsr()
    n_steps = round(T_END / DT)
    err_all = 0.0
    energy0 = M @ T
    for n in range(1, n_steps + 1):
        T = lhs.solve(rhs_mat @ T)
        exact = 1.0 + shape * np.exp(-ALPHA * LAM**2 * (n * DT) / R**2)
        err_all = max(err_all, np.abs(T - exact).max())
    # Aislado: la energía Σ (M T)_i se conserva (1ᵀK = 0).
    assert abs((M @ T).sum() - energy0.sum()) / energy0.sum() < 1e-12
    err_end = np.abs(T - exact).max()
    return err_end / 0.5, err_all / 0.5


@pytest.fixture(scope="module")
def errors():
    return {nx: run_mode(nx) for nx in (10, 20, 40, 80)}


def test_t6_error_below_one_percent_nx20(errors):
    err_end, err_all = errors[20]
    assert err_end < 0.01
    assert err_all < 0.01


def test_t6_error_decreases_with_refinement(errors):
    e = [errors[nx][0] for nx in (10, 20, 40, 80)]
    assert all(e[i + 1] < e[i] for i in range(3)), e
