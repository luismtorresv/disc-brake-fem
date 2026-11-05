"""Integración temporal con el esquema θ (PLAN.md §3.5) y balance de energía (§3.6).

    [M + θΔt(K+H)] Tⁿ⁺¹ = [M − (1−θ)Δt(K+H)] Tⁿ + Δt [θ Fⁿ⁺¹ + (1−θ) Fⁿ]

θ = 1 es Euler implícito y θ = 1/2 es Crank–Nicolson. La matriz del lado izquierdo se
factoriza con ``scipy.sparse.linalg.splu`` una vez por cada par (θ, Δt) distinto y se reutiliza.
Con Δt constante hay una sola factorización; con Rannacher (n pasos de θ = 1 y luego θ = 1/2)
hay dos.

Elección del valor por defecto: θ = 1/2 con 2 pasos iniciales de Rannacher
(`DEFAULT_THETA`, `DEFAULT_RANNACHER_STEPS`). Se basa en T7: semi-infinito plano con flujo
q0(1 − t/t_f), 40 capas con la primera de 1,5e-4 m. Error relativo de la temperatura de la
superficie frente a la solución exacta (salida de ``python tests/test_timestepping.py``):

    Δt [s]   esquema             máx. t ≥ 0,2 s   en t = 3 s
    0,005    θ = 1               0,395 %          0,090 %
    0,005    θ = 1/2             0,070 %          0,023 %
    0,005    θ = 1/2 + 2 Rann.   0,068 %          0,023 %
    0,01     θ = 1               0,737 %          0,158 %
    0,01     θ = 1/2             0,595 %          0,022 %
    0,01     θ = 1/2 + 2 Rann.   0,111 %          0,023 %
    0,0556   θ = 1               3,452 %          0,773 %
    0,0556   θ = 1/2             6,010 %          0,502 %
    0,0556   θ = 1/2 + 2 Rann.   1,482 %          0,048 %

(Con Δt = 0,0556 s, la columna "en t = 3 s" se evalúa en el instante de la malla más cercano,
t = 3,0024 s.) Δt = 0,0556 s = t_f/100 es el paso por defecto del caso de referencia. Con ese Δt:
- θ = 1 es de primer orden: su error decae lento (0,77 % aún en t = 3 s);
- Crank–Nicolson puro oscila tras el escalón de flujo en t = 0 (6 % al inicio);
- con dos pasos de Euler implícito al arranque (Rannacher), el amortiguamiento de θ = 1
  elimina la oscilación y luego queda el segundo orden de CN: es el mejor en todas las filas.
Cuesta una factorización LU adicional. El caso de referencia del PLAN (§4) fija θ = 1 de forma
explícita en su YAML, y con eso T8/T9 se comparan contra los valores del prototipo.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

#: θ por defecto (Crank–Nicolson); justificación numérica en el docstring del módulo.
DEFAULT_THETA = 0.5
#: Pasos iniciales con θ = 1 (Rannacher) por defecto; ver el docstring del módulo.
DEFAULT_RANNACHER_STEPS = 2
#: Pasos de tiempo por defecto durante el frenado: Δt = t_f / 100.
DEFAULT_STEPS_PER_TF = 100


def time_grid(
    t_end: float,
    t_f: float | None = None,
    dt: float | str | None = "auto",
    steps_per_tf: int = DEFAULT_STEPS_PER_TF,
) -> np.ndarray:
    """Malla temporal 0 = t_0 < … < t_n = t_end que contiene t_f exactamente.

    - ``dt = "auto"``: Δt = t_f / ``steps_per_tf`` (requiere t_f).
    - ``dt`` numérico con t_f: se ajusta a Δt = t_f / round(t_f/Δt) para que t_f sea múltiplo
      entero de Δt (PLAN.md §3.5; trampa 5 de la §7).
    - Después de t_f se sigue con el mismo Δt; si t_end no es múltiplo, el último paso es más
      corto (eso cuesta una factorización adicional).
    """
    if not t_end > 0:
        raise ValueError(f"time.t_end debe ser > 0; se recibió {t_end}.")
    if t_f is not None and not t_f > 0:
        raise ValueError(f"t_f debe ser > 0; se recibió {t_f}.")
    if dt is None or dt == "auto":
        if t_f is None:
            raise ValueError("time.dt = 'auto' requiere un tiempo de frenado t_f.")
        dt = t_f / steps_per_tf
    dt = float(dt)
    if not dt > 0:
        raise ValueError(f"time.dt debe ser > 0; se recibió {dt}.")

    if t_f is not None and t_f < t_end:
        n_f = max(1, round(t_f / dt))
        dt = t_f / n_f
        head = np.linspace(0.0, t_f, n_f + 1)
    else:
        head = np.array([0.0])
    start = head[-1]
    if t_end <= start:
        return head
    ratio = (t_end - start) / dt
    n_rest = round(ratio)
    if abs(ratio - n_rest) <= 1e-9 * max(ratio, 1.0) and n_rest >= 1:
        # t_end es múltiplo de Δt (salvo redondeo): pasos idénticos, una sola factorización.
        tail = np.linspace(start, t_end, n_rest + 1)[1:]
    else:
        n_full = int(np.floor(ratio))
        tail = start + dt * np.arange(1, n_full + 1)
        tail = np.append(tail, t_end)  # último paso más corto
    return np.concatenate([head, tail])


@dataclass
class TransientResult:
    """Resultado de :func:`integrate`. Energías en J/rad (SIN el factor 2π; ver CLAUDE.md).

    Atributos:
        times: (n+1,) instantes de todos los pasos.
        saved_times, saved_T: instantes y campos guardados (cada ``save_every`` pasos, más el
            último).
        T_max, T_max_node: (n+1,) máximo nodal de T en cada instante y el nodo donde ocurre.
        energy_stored: (n+1,) 1ᵀM(Tⁿ − T⁰).
        energy_brake: (n+1,) energía acumulada que entra por la carga F_brake.
        energy_conv: (n+1,) energía acumulada que sale por convección, Σ Δt·1ᵀ(H T − F_conv).
        thetas: (n,) θ usado en cada paso.
        n_factorizations: número de factorizaciones LU realizadas.
    """

    times: np.ndarray
    saved_times: np.ndarray
    saved_T: np.ndarray
    T_max: np.ndarray
    T_max_node: np.ndarray
    energy_stored: np.ndarray
    energy_brake: np.ndarray
    energy_conv: np.ndarray
    thetas: np.ndarray
    n_factorizations: int
    extra: dict = field(default_factory=dict)

    @property
    def T_final(self) -> np.ndarray:
        return self.saved_T[-1]

    @property
    def balance_residual(self) -> np.ndarray:
        """Residuo absoluto del balance (§3.6) en cada instante:
        1ᵀM(Tⁿ − T⁰) − (E_brake − E_conv)."""
        return self.energy_stored - (self.energy_brake - self.energy_conv)

    def balance_relative(self) -> float:
        """max |residuo| / escala, con escala = max(|E_brake|, |E_conv|, |E_stored|)."""
        scale = max(
            np.abs(self.energy_brake).max(),
            np.abs(self.energy_conv).max(),
            np.abs(self.energy_stored).max(),
        )
        if scale == 0:
            return 0.0
        return float(np.abs(self.balance_residual).max() / scale)


def integrate(
    M: sp.spmatrix,
    K: sp.spmatrix,
    H: sp.spmatrix | None,
    T0: np.ndarray | float,
    times: np.ndarray,
    F_brake: Callable[[float], np.ndarray] | None = None,
    t_f: float | None = None,
    F_conv: np.ndarray | None = None,
    theta: float = DEFAULT_THETA,
    rannacher_steps: int = DEFAULT_RANNACHER_STEPS,
    save_every: int = 1,
) -> TransientResult:
    """Integra M Ṫ + (K + H) T = F_brake(t) + F_conv con el esquema θ.

    Args:
        M, K, H: matrices ensambladas (H puede ser None si no hay convección).
        T0: temperatura inicial (escalar o vector nodal).
        times: malla temporal (p. ej. de :func:`time_grid`); debe contener ``t_f``.
        F_brake: función t ↦ vector de carga del flujo de fricción (calor que entra).
        t_f: fin del frenado. Un paso [tⁿ, tⁿ⁺¹] con tⁿ ≥ t_f usa F_brake = 0 en ambos
            extremos: así un flujo que no se anula en t_f no "se filtra" al paso siguiente a
            través del término (1−θ)Fⁿ.
        F_conv: vector constante ∫ h T∞ φ_i w ds (None ⇒ 0).
        theta: θ ∈ [0, 1]; ver ``DEFAULT_THETA``.
        rannacher_steps: número de pasos iniciales con θ = 1 (irrelevante si θ = 1).
        save_every: guarda el campo cada tantos pasos (siempre guarda el primero y el último).
    """
    if not 0.0 <= theta <= 1.0:
        raise ValueError(f"time.theta debe estar en [0, 1]; se recibió {theta}.")
    if rannacher_steps < 0:
        raise ValueError(f"time.rannacher_steps debe ser ≥ 0; se recibió {rannacher_steps}.")
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or times.size < 2 or np.any(np.diff(times) <= 0):
        raise ValueError("La malla temporal debe ser creciente y tener al menos dos instantes.")
    tf_missing = t_f is not None and t_f < times[-1] and not np.any(
        np.abs(times - t_f) <= 1e-12 * max(t_f, 1.0)
    )
    if tf_missing:
        raise ValueError(
            f"La malla temporal no contiene t_f = {t_f} exactamente: el flujo se cortaría "
            "a mitad de paso y el balance de energía no cerraría (PLAN.md §7, trampa 5)."
        )

    n = M.shape[0]
    M = sp.csr_matrix(M)
    KH = sp.csr_matrix(K) if H is None else sp.csr_matrix(K + H)
    ones = np.ones(n)
    m_sum = M.T @ ones
    h_sum = np.zeros(n) if H is None else sp.csr_matrix(H).T @ ones
    fc = np.zeros(n) if F_conv is None else np.asarray(F_conv, dtype=float)
    fc_sum = fc.sum()

    T = np.full(n, float(T0)) if np.ndim(T0) == 0 else np.array(T0, dtype=float)
    T_init = T.copy()

    n_steps = times.size - 1
    T_max = np.empty(n_steps + 1)
    T_max_node = np.empty(n_steps + 1, dtype=np.int64)
    e_stored = np.zeros(n_steps + 1)
    e_brake = np.zeros(n_steps + 1)
    e_conv = np.zeros(n_steps + 1)
    thetas = np.empty(n_steps)
    T_max_node[0] = int(np.argmax(T))
    T_max[0] = T[T_max_node[0]]
    saved_t, saved_T = [times[0]], [T.copy()]

    # Caché de factorizaciones: (θ, Δt) → (LU, B). Dos pasos cuyo Δt difiere < 1e-9 relativo
    # (solo por redondeo de la malla temporal) comparten la factorización, y ese paso usa el Δt
    # de la entrada del caché en todo el esquema (LU, lado derecho y balance), para que el
    # paso siga siendo exactamente un paso θ.
    cache: list[tuple[float, float, object, sp.csr_matrix]] = []
    zero = np.zeros(n)
    fb_prev_t, fb_prev = None, None

    for k in range(n_steps):
        t0, t1 = times[k], times[k + 1]
        th = 1.0 if k < rannacher_steps else theta
        thetas[k] = th
        dt_grid = t1 - t0
        hit = next(
            (c for c in cache if c[0] == th and abs(c[1] - dt_grid) <= 1e-9 * dt_grid), None
        )
        if hit is None:
            hit = (
                th,
                dt_grid,
                spla.splu(sp.csc_matrix(M + th * dt_grid * KH)),
                sp.csr_matrix(M - (1.0 - th) * dt_grid * KH),
            )
            cache.append(hit)
        _, dt, lu, B = hit

        active = F_brake is not None and (t_f is None or t0 < t_f * (1 - 1e-12))
        if active:
            fb0 = fb_prev if fb_prev_t == t0 else np.asarray(F_brake(t0), dtype=float)
            fb1 = np.asarray(F_brake(t1), dtype=float)
            fb_prev_t, fb_prev = t1, fb1
        else:
            fb0 = fb1 = zero

        rhs = B @ T + dt * (th * fb1 + (1.0 - th) * fb0) + dt * fc
        T_new = lu.solve(rhs)

        e_brake[k + 1] = e_brake[k] + dt * (th * fb1.sum() + (1.0 - th) * fb0.sum())
        q_out_new = h_sum @ T_new - fc_sum
        q_out_old = h_sum @ T - fc_sum
        e_conv[k + 1] = e_conv[k] + dt * (th * q_out_new + (1.0 - th) * q_out_old)
        e_stored[k + 1] = m_sum @ (T_new - T_init)

        T = T_new
        T_max_node[k + 1] = int(np.argmax(T))
        T_max[k + 1] = T[T_max_node[k + 1]]
        if (k + 1) % save_every == 0 or k + 1 == n_steps:
            saved_t.append(t1)
            saved_T.append(T.copy())

    return TransientResult(
        times=times,
        saved_times=np.array(saved_t),
        saved_T=np.array(saved_T),
        T_max=T_max,
        T_max_node=T_max_node,
        energy_stored=e_stored,
        energy_brake=e_brake,
        energy_conv=e_conv,
        thetas=thetas,
        n_factorizations=len(cache),
    )
