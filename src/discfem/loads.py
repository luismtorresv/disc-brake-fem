"""Flujo de calor entrante en ``brake``: q_fric(r, t) [W/m²] (PLAN.md §3.4).

Los tres modelos se separan como q(r, t) = f(t) · s(r):
- ``f(t)`` (:meth:`time_factor`) es adimensional y vale 0 para t > t_f;
- ``s(r)`` (:meth:`shape`) está en W/m².
Así el solver ensambla una sola vez F_s,i = ∫_brake s(r) φ_i r ds y en cada paso usa
F_brake(t) = f(t)·F_s (ver :func:`brake_load_vector`).

Modelos:
- ``braking``: frenada con desaceleración constante, v(t) = v0 − a t, t_f = v0/a,
      q(r, t) = γ φ m a v(t) / (2π n_d n_c) · g(r) r / I_g,   I_g = ∫_{r_a}^{r_b} g(r) r² dr,
  con g = 1 (``uniform_pressure``) o g = 1/r (``uniform_wear``). Cumple 2π ∫ q r dr = γ P_f(t)
  con P_f(t) = φ m a v(t) / (n_d n_c), la potencia de fricción en UNA cara de un disco.
- ``constant_flux``: q = q0 uniforme en r, durante ``duration`` segundos (o toda la corrida).
- ``table``: q(t) uniforme en r, interpolación lineal en la tabla [[t, q], ...]; q = 0 después
  del último instante.

En ``constant_flux`` y ``table`` el valor dado ya es el flujo que ENTRA al disco: no se
multiplica por γ.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

_PROFILES = ("uniform_pressure", "uniform_wear")
_VEHICLE_KEYS = {
    "mass", "v0", "decel", "wheel_radius", "brake_fraction", "discs_per_axle", "faces_per_disc",
}


class Load:
    """Interfaz común. ``t_f`` es el fin del calentamiento (None: activo toda la corrida)."""

    t_f: float | None = None

    def time_factor(self, t):
        raise NotImplementedError

    def shape(self, r):
        raise NotImplementedError

    def q(self, r, t):
        """q(r, t) = f(t)·s(r) [W/m²]."""
        return self.time_factor(t) * self.shape(r)

    def _active(self, t):
        t = np.asarray(t, dtype=float)
        return np.ones_like(t) if self.t_f is None else (t <= self.t_f).astype(float)


@dataclass(frozen=True)
class BrakingLoad(Load):
    """Frenada con desaceleración prescrita (modelo ``braking``).

    Atributos (SI): mass m [kg], v0 [m/s], decel a [m/s²], wheel_radius R_w [m] (solo para ω),
    brake_fraction φ [-], discs_per_axle n_d, faces_per_disc n_c, gamma γ [-], r_a, r_b [m],
    profile ``uniform_pressure`` | ``uniform_wear``.
    """

    mass: float
    v0: float
    decel: float
    brake_fraction: float
    discs_per_axle: int
    faces_per_disc: int
    gamma: float
    r_a: float
    r_b: float
    profile: str = "uniform_wear"
    wheel_radius: float | None = None

    def __post_init__(self) -> None:
        for name in ("mass", "v0", "decel"):
            if not getattr(self, name) > 0:
                raise ValueError(f"load.vehicle.{name} debe ser > 0; se recibió {getattr(self, name)}.")
        if not 0 < self.brake_fraction <= 1:
            raise ValueError(
                f"load.vehicle.brake_fraction debe estar en (0, 1]; se recibió {self.brake_fraction}."
            )
        for name in ("discs_per_axle", "faces_per_disc"):
            v = getattr(self, name)
            if int(v) != v or v < 1:
                raise ValueError(f"load.vehicle.{name} debe ser un entero ≥ 1; se recibió {v}.")
        if not 0 < self.gamma < 1:
            raise ValueError(f"gamma debe estar en (0, 1); se recibió {self.gamma}.")
        if self.profile not in _PROFILES:
            raise ValueError(
                f"load.pressure_profile debe ser uno de {_PROFILES}; se recibió '{self.profile}'."
            )
        if not 0 <= self.r_a < self.r_b:
            raise ValueError(f"Pista inválida: r_a = {self.r_a}, r_b = {self.r_b}.")
        if self.wheel_radius is not None and not self.wheel_radius > 0:
            raise ValueError("load.vehicle.wheel_radius debe ser > 0.")

    # --------------------------------------------------------------- cinemática y potencia
    @property
    def t_f(self) -> float:  # type: ignore[override]
        """Tiempo de frenado t_f = v0 / a [s]."""
        return self.v0 / self.decel

    def velocity(self, t):
        """v(t) = v0 − a t para t ≤ t_f; 0 después [m/s]."""
        t = np.asarray(t, dtype=float)
        return np.where(t <= self.t_f, np.maximum(self.v0 - self.decel * t, 0.0), 0.0)

    def angular_velocity(self, t):
        """ω(t) = v(t) / R_w [rad/s] (requiere wheel_radius)."""
        if self.wheel_radius is None:
            raise ValueError("Para ω se necesita load.vehicle.wheel_radius.")
        return self.velocity(t) / self.wheel_radius

    def power_per_face(self, t):
        """P_f(t) = φ m a v(t) / (n_d n_c) [W], potencia de fricción en UNA cara (antes de γ)."""
        return (self.brake_fraction * self.mass * self.decel * self.velocity(t)
                / (self.discs_per_axle * self.faces_per_disc))

    @property
    def energy_per_face(self) -> float:
        """∫ P_f dt = φ m v0² / (2 n_d n_c) [J] (antes de γ)."""
        return (self.brake_fraction * self.mass * self.v0**2
                / (2 * self.discs_per_axle * self.faces_per_disc))

    # --------------------------------------------------------------- perfil radial
    @property
    def I_g(self) -> float:
        """I_g = ∫_{r_a}^{r_b} g(r) r² dr: (r_b³ − r_a³)/3 si g = 1; (r_b² − r_a²)/2 si g = 1/r."""
        if self.profile == "uniform_pressure":
            return (self.r_b**3 - self.r_a**3) / 3.0
        return (self.r_b**2 - self.r_a**2) / 2.0

    def g(self, r):
        """Perfil de presión g(r): 1 (presión uniforme) o 1/r (desgaste uniforme)."""
        r = np.asarray(r, dtype=float)
        return np.ones_like(r) if self.profile == "uniform_pressure" else 1.0 / r

    def time_factor(self, t):
        """f(t) = v(t)/v0 (= 1 − t/t_f en la frenada, 0 después)."""
        return self.velocity(t) / self.v0

    def shape(self, r):
        """s(r) = q(r, 0) = γ φ m a v0 / (2π n_d n_c) · g(r) r / I_g [W/m²]; 0 fuera de [r_a, r_b]."""
        r = np.asarray(r, dtype=float)
        c = (self.gamma * self.brake_fraction * self.mass * self.decel * self.v0
             / (2 * math.pi * self.discs_per_axle * self.faces_per_disc))
        tol = 1e-12 * self.r_b
        inside = (r >= self.r_a - tol) & (r <= self.r_b + tol)
        return np.where(inside, c * self.g(np.where(inside, r, self.r_b)) * r / self.I_g, 0.0)


@dataclass(frozen=True)
class ConstantFluxLoad(Load):
    """Flujo entrante constante q0 [W/m²], uniforme en r, durante ``duration`` [s] (None: siempre)."""

    q0: float
    duration: float | None = None

    def __post_init__(self) -> None:
        if not self.q0 >= 0:
            raise ValueError(f"load.q debe ser ≥ 0; se recibió {self.q0}.")
        if self.duration is not None and not self.duration > 0:
            raise ValueError(f"load.duration debe ser > 0; se recibió {self.duration}.")

    @property
    def t_f(self) -> float | None:  # type: ignore[override]
        return self.duration

    def time_factor(self, t):
        return self._active(t)

    def shape(self, r):
        return np.full_like(np.asarray(r, dtype=float), self.q0)


@dataclass(frozen=True)
class TableLoad(Load):
    """q(t) [W/m²] uniforme en r por interpolación lineal en ``times``/``values``; 0 después
    del último instante (que es t_f)."""

    times: tuple[float, ...]
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        t = np.asarray(self.times, dtype=float)
        v = np.asarray(self.values, dtype=float)
        if t.ndim != 1 or t.size < 2 or t.shape != v.shape:
            raise ValueError("load.points debe ser una lista de al menos dos pares [t, q].")
        if t[0] != 0.0:
            raise ValueError(f"La tabla de carga debe empezar en t = 0; empieza en {t[0]}.")
        if np.any(np.diff(t) <= 0):
            raise ValueError("Los instantes de la tabla de carga deben ser estrictamente crecientes.")
        if np.any(v < 0):
            raise ValueError("Los flujos de la tabla deben ser ≥ 0 (calor que entra).")

    @property
    def t_f(self) -> float:  # type: ignore[override]
        return float(self.times[-1])

    def time_factor(self, t):
        """q_tabla(t) en W/m² (la forma espacial es 1): interpolación lineal, 0 para t > t_f."""
        t = np.asarray(t, dtype=float)
        return np.interp(t, self.times, self.values) * self._active(t)

    def shape(self, r):
        return np.ones_like(np.asarray(r, dtype=float))


def make_load(cfg: dict, r_a: float, r_b: float, gamma: float | None = None) -> Load:
    """Construye la carga desde la sección ``load`` del YAML.

    Args:
        cfg: sección ``load``.
        r_a, r_b: radios de la pista (de la geometría).
        gamma: γ ya resuelto (ver ``materials.resolve_gamma``); obligatorio para ``braking``.
    """
    ltype = cfg.get("type")
    if ltype == "braking":
        unknown = set(cfg) - {"type", "vehicle", "pressure_profile", "gamma"}
        if unknown:
            raise ValueError(f"Claves desconocidas en 'load': {sorted(unknown)}.")
        veh = cfg.get("vehicle") or {}
        unknown = set(veh) - _VEHICLE_KEYS
        if unknown:
            raise ValueError(f"Claves desconocidas en 'load.vehicle': {sorted(unknown)}.")
        missing = sorted(_VEHICLE_KEYS - {"wheel_radius"} - set(veh))
        if missing:
            raise ValueError(f"Faltan claves en 'load.vehicle': {missing}.")
        if gamma is None:
            raise ValueError("El modelo 'braking' necesita gamma resuelto.")
        return BrakingLoad(
            mass=float(veh["mass"]),
            v0=float(veh["v0"]),
            decel=float(veh["decel"]),
            brake_fraction=float(veh["brake_fraction"]),
            discs_per_axle=veh["discs_per_axle"],
            faces_per_disc=veh["faces_per_disc"],
            wheel_radius=None if veh.get("wheel_radius") is None else float(veh["wheel_radius"]),
            gamma=float(gamma),
            r_a=r_a,
            r_b=r_b,
            profile=cfg.get("pressure_profile", "uniform_wear"),
        )
    if ltype == "constant_flux":
        unknown = set(cfg) - {"type", "q", "duration"}
        if unknown:
            raise ValueError(f"Claves desconocidas en 'load': {sorted(unknown)}.")
        if "q" not in cfg:
            raise ValueError("load.type = constant_flux requiere 'q' [W/m²].")
        dur = cfg.get("duration")
        return ConstantFluxLoad(q0=float(cfg["q"]), duration=None if dur is None else float(dur))
    if ltype == "table":
        unknown = set(cfg) - {"type", "points"}
        if unknown:
            raise ValueError(f"Claves desconocidas en 'load': {sorted(unknown)}.")
        pts = np.asarray(cfg.get("points", []), dtype=float)
        if pts.ndim != 2 or pts.shape[1] != 2:
            raise ValueError("load.points debe ser una lista de pares [t, q].")
        return TableLoad(times=tuple(pts[:, 0]), values=tuple(pts[:, 1]))
    raise ValueError(
        f"load.type debe ser 'braking', 'constant_flux' o 'table'; se recibió {ltype!r}."
    )


def brake_load_vector(mesh, load: Load, axisymmetric: bool = True, groups=("brake",)):
    """Devuelve F_brake(t) = f(t)·F_s, con F_s,i = Σ_grupos ∫ s(r) φ_i w ds ensamblado una vez.

    ``groups`` son los grupos con condición ``flux`` (por defecto solo ``brake``).
    """
    from .fem import assemble_edge_load

    F_s = assemble_edge_load(mesh, list(groups), lambda r, z: load.shape(r), axisymmetric)
    return lambda t: float(load.time_factor(t)) * F_s
