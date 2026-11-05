"""Pruebas de la Fase 3 (PLAN.md §6): T4 (balance de energía), T5 (estado uniforme), T7
(semi-infinito con solución exacta) y la malla temporal.

Ejecutar ``python tests/test_timestepping.py`` imprime la tabla de errores de T7 por esquema
que justifica el θ por defecto (ver el docstring de ``discfem.timestepping``).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import brentq

from discfem.fem import assemble_edge_load, assemble_robin, assemble_stiffness_mass
from discfem.geometry import RectangleGeometry
from discfem.mesh import MeshOptions, build_mesh
from discfem.mesh.structured import tensor_mesh
from discfem.timestepping import integrate, time_grid

K_, RHO, CP = 50.0, 7200.0, 460.0
SCHEMES = [(1.0, 0), (0.5, 0), (0.5, 2)]  # (θ, pasos de Rannacher)
SCHEME_IDS = ["theta1", "CN", "CN+Rannacher2"]


# ----------------------------------------------------------------------------- time_grid

def test_time_grid_auto_contains_tf():
    t = time_grid(60.0, t_f=5.56)
    assert t[0] == 0.0 and t[-1] == 60.0
    assert np.any(t == 5.56)
    dt = np.diff(t)
    np.testing.assert_allclose(dt[:100], 5.56 / 100, rtol=1e-12)
    assert np.sum(t <= 5.56) == 101


def test_time_grid_adjusts_dt_to_divide_tf():
    t = time_grid(10.0, t_f=5.56, dt=0.03)  # 5.56/0.03 = 185.33 → 185 pasos
    assert np.any(t == 5.56)
    np.testing.assert_allclose(np.diff(t)[:185], 5.56 / 185, rtol=1e-12)


def test_time_grid_short_last_step():
    t = time_grid(1.05, t_f=None, dt=0.1)
    np.testing.assert_allclose(t[-2:], [1.0, 1.05], rtol=1e-12)


@pytest.mark.parametrize("kw", [{"t_end": 0.0, "t_f": 1.0}, {"t_end": 1.0, "dt": "auto"},
                                {"t_end": 1.0, "t_f": 1.0, "dt": -1.0}])
def test_time_grid_errors(kw):
    with pytest.raises(ValueError):
        time_grid(**kw)


def test_integrate_rejects_grid_without_tf():
    m = tensor_mesh(np.linspace(0, 1, 3), np.linspace(0, 1, 3), 0, 2)
    K, M = assemble_stiffness_mass(m, 1.0, 1.0)
    with pytest.raises(ValueError, match="t_f"):
        integrate(M, K, None, 0.0, np.linspace(0, 1, 4), F_brake=lambda t: 0 * M @ np.ones(9),
                  t_f=0.5)


# ------------------------------------------------------------------- caso de referencia

GEO = RectangleGeometry(r_inner=0.07, r_outer=0.15, half_thickness=0.01, r_a=0.08, r_b=0.14)
CONV_GROUPS = ["top_inner", "top_outer", "r_inner", "r_outer"]
Q0, TF = 7.7e5, 5.56


@pytest.fixture(scope="module")
def ref_system():
    m = build_mesh(GEO, MeshOptions(backend="structured", size_track=1e-3, size_far=4e-3))
    K, M = assemble_stiffness_mass(m, K_, RHO * CP, True)
    H = assemble_robin(m, dict.fromkeys(CONV_GROUPS, 40.0), True)
    Fq = assemble_edge_load(m, ["brake"], 1.0, True)  # carga para q = 1 W/m²
    return m, K, M, H, Fq


def brake_load(Fq):
    return lambda t: Q0 * (1.0 - t / TF) * Fq


@pytest.mark.parametrize(("theta", "rs"), SCHEMES, ids=SCHEME_IDS)
@pytest.mark.parametrize("convection", [False, True], ids=["sin_conv", "con_conv"])
def test_t4_energy_balance(ref_system, theta, rs, convection):
    """T4: 1ᵀM(Tⁿ − T⁰) = E_brake − E_conv, cierre < 1e-10 relativo."""
    m, K, M, H, Fq = ref_system
    t = time_grid(20.0, t_f=TF)
    F_conv = H @ np.full(m.n_nodes, 25.0) if convection else None
    res = integrate(M, K, H if convection else None, 60.0 if convection else 25.0, t,
                    F_brake=brake_load(Fq), t_f=TF, F_conv=F_conv,
                    theta=theta, rannacher_steps=rs)
    assert res.balance_relative() < 1e-10
    if convection:
        assert res.energy_conv[-1] > 0  # T0 = 60 > T∞: sale calor
    else:
        assert res.energy_conv[-1] == 0


@pytest.mark.parametrize(("theta", "rs"), SCHEMES, ids=SCHEME_IDS)
def test_t4_brake_energy_matches_time_rule(ref_system, theta, rs):
    """La energía que entra es la regla θ aplicada a q0(1 − t/t_f) en [0, t_f], exactamente;
    después de t_f no entra nada (trampa 5)."""
    m, K, M, H, Fq = ref_system  # noqa: RUF059
    N = 100
    t = time_grid(20.0, t_f=TF)
    res = integrate(M, K, None, 25.0, t, F_brake=brake_load(Fq), t_f=TF,
                    theta=theta, rannacher_steps=rs)
    S = Q0 * (GEO.r_b**2 - GEO.r_a**2) / 2  # Σ F para q = q0, sin 2π
    dt = TF / N
    g = 1.0 - np.arange(N + 1) / N  # q(t_n)/q0
    th = np.where(np.arange(N) < rs, 1.0, theta)
    exact = S * dt * np.sum(th * g[1:] + (1 - th) * g[:-1])
    assert abs(res.energy_brake[-1] - exact) / exact < 1e-12
    i_f = np.searchsorted(t, TF)
    assert np.all(res.energy_brake[i_f:] == res.energy_brake[i_f])


def test_constant_flux_does_not_leak_after_tf(ref_system):
    """Flujo que NO se anula en t_f: con θ = 1/2, un paso que empieza en t_f no debe aportar
    (1 − θ)Δt·F(t_f). La energía de entrada debe ser exactamente q·t_f·∫r dr."""
    m, K, M, H, Fq = ref_system  # noqa: RUF059
    q, tf = 1e5, 1.0
    t = time_grid(3.0, t_f=tf, dt=0.1)
    res = integrate(M, K, None, 0.0, t, F_brake=lambda s: q * Fq, t_f=tf, theta=0.5,
                    rannacher_steps=0)
    exact = q * tf * (GEO.r_b**2 - GEO.r_a**2) / 2
    assert abs(res.energy_brake[-1] - exact) / exact < 1e-12


def test_factorizations_reused(ref_system):
    _, K, M, H, _ = ref_system
    t = time_grid(3 * TF, t_f=TF)  # t_end múltiplo de Δt: todos los pasos iguales
    r1 = integrate(M, K, H, 25.0, t, theta=1.0, rannacher_steps=0)
    r2 = integrate(M, K, H, 25.0, t, theta=0.5, rannacher_steps=2)
    r3 = integrate(M, K, H, 25.0, time_grid(20.0, t_f=TF), theta=1.0)
    assert r1.n_factorizations == 1
    assert r2.n_factorizations == 2
    assert list(r2.thetas[:3]) == [1.0, 1.0, 0.5]
    assert r3.n_factorizations == 2  # 20 s no es múltiplo de t_f/100: último paso más corto


@pytest.mark.parametrize(("theta", "rs"), SCHEMES, ids=SCHEME_IDS)
def test_t5_uniform_state_is_preserved(ref_system, theta, rs):
    """T5: T0 = T∞ y sin flujo ⇒ T queda uniforme, |T − T0| < 1e-10."""
    m, K, M, H, Fq = ref_system  # noqa: RUF059
    T_inf = 25.0
    res = integrate(M, K, H, T_inf, time_grid(30.0, t_f=TF), F_conv=H @ np.full(m.n_nodes, T_inf),
                    theta=theta, rannacher_steps=rs)
    assert np.abs(res.saved_T - T_inf).max() < 1e-10


# ---------------------------------------------------------------------------------- T7

E7, W7, Q07, TF7 = 0.06, 0.01, 7.7e5, 5.56


def semi_infinite_mesh(h1: float = 1.5e-4, n_layers: int = 40):
    """Columna plana de ancho 0,01 (2 celdas) y alto 0,06 (40 capas geométricas, la primera de
    espesor h1 junto a la cara calentada z = 0,06). `brake` cubre toda la cara superior."""
    q = brentq(lambda s: h1 * (s**n_layers - 1) / (s - 1) - E7, 1.0 + 1e-6, 2.0)
    d = np.concatenate([[0.0], np.cumsum(h1 * q ** np.arange(n_layers))])
    d[-1] = E7
    z = (E7 - d)[::-1]
    return tensor_mesh(np.linspace(0.0, W7, 3), z, 0, 2)


def run_t7(theta: float, rs: int, dt: float, h1: float = 1.5e-4):
    """Devuelve (error relativo máx. para t ≥ 0,2 s, error relativo en el instante de la malla
    más cercano a 3 s, residuo relativo del balance, ese instante)."""
    m = semi_infinite_mesh(h1)
    K, M = assemble_stiffness_mass(m, K_, RHO * CP, axisymmetric=False)
    Fq = assemble_edge_load(m, ["brake"], 1.0, axisymmetric=False)
    t = time_grid(TF7, t_f=TF7, dt=dt)
    res = integrate(M, K, None, 0.0, t, F_brake=lambda s: Q07 * (1 - s / TF7) * Fq, t_f=TF7,
                    theta=theta, rannacher_steps=rs)
    top = np.flatnonzero(m.nodes[:, 1] == E7)
    Ts = res.saved_T[:, top].mean(axis=1)
    exact = 2 * Q07 / np.sqrt(np.pi * K_ * RHO * CP) * (np.sqrt(t) - 2 * t**1.5 / (3 * TF7))
    mask = t >= 0.2 - 1e-12
    rel = np.abs(Ts[mask] - exact[mask]) / exact[mask]
    i3 = int(np.argmin(np.abs(t - 3.0)))
    return rel.max(), abs(Ts[i3] - exact[i3]) / exact[i3], res.balance_relative(), t[i3]


@pytest.mark.parametrize("dt", [0.005, 0.01])
def test_t7_implicit_euler(dt):
    err_max, err_3, bal, t3 = run_t7(1.0, 0, dt)
    assert abs(t3 - 3.0) < 1e-9  # con estos Δt, t = 3 s está en la malla
    assert err_max < 0.01
    assert err_3 < 0.002
    assert bal < 1e-10


@pytest.mark.parametrize("dt", [0.005, 0.01])
def test_t7_crank_nicolson(dt):
    err_max, _, bal, _ = run_t7(0.5, 0, dt)
    assert err_max < 0.03
    assert bal < 1e-10


@pytest.mark.parametrize("dt", [0.005, 0.01])
def test_t7_crank_nicolson_rannacher(dt):
    """Con 2 pasos de Rannacher se exige lo mismo que a θ = 1."""
    err_max, err_3, _, _ = run_t7(0.5, 2, dt)
    assert err_max < 0.01
    assert err_3 < 0.002


if __name__ == "__main__":
    import sys

    sys.stdout.reconfigure(encoding="utf-8", errors="replace")  # consola cp1252 de Windows
    print(f"{'Δt [s]':>8}  {'esquema':<18} {'máx. t≥0,2 s':>13} {'en t≈3 s':>9} "
          f"{'(t usado)':>9} {'balance':>9}")
    for dt in (0.005, 0.01, TF7 / 100):
        for (theta, rs), name in zip(SCHEMES, ["θ = 1", "θ = 1/2", "θ = 1/2 + 2 Rann."],
                                     strict=True):
            a, b, bal, t3 = run_t7(theta, rs, dt)
            print(f"{dt:8.4f}  {name:<18} {100 * a:12.3f}% {100 * b:8.3f}% {t3:9.4f} {bal:9.1e}")
