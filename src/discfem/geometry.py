"""Definición de la sección axisimétrica (r, z) del disco, sus grupos de frontera y sus regiones.

Dos tipos de sección (``geometry.type`` en el YAML):

- ``rectangle``: Ω = (r_i, r_o) × (0, e), con z = 0 como plano de simetría y la pista ``brake``
  en z = e, r_a < r < r_b (PLAN.md §3.1).
- ``profile`` (Fase 7): polígono (r, z) antihorario. Cada lado lleva el nombre de su grupo de
  frontera (``segments``); varios lados pueden compartir grupo. ``brake`` debe ser una cadena
  contigua de lados horizontales (z constante): el flujo q(r) y su normalización I_g suponen
  la pista en un plano z = cte. r_a y r_b se derivan de ella.

En ambos tipos, ``regions`` (opcional) asigna nombres a sub-polígonos con otro material
(p. ej. una capa bajo la pista). Lo que no cubre ninguna región es la región ``disc``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

#: Nombres de los grupos de frontera, en el orden en que se recorre el borde del rectángulo en
#: sentido antihorario empezando por el plano de simetría. Cada arista de borde pertenece a uno.
BOUNDARY_GROUPS: tuple[str, ...] = (
    "symmetry",   # z = 0
    "r_outer",    # r = r_o
    "top_outer",  # z = e, r_b..r_o
    "brake",      # z = e, r_a..r_b
    "top_inner",  # z = e, r_i..r_a
    "r_inner",    # r = r_i
)
#: Región por defecto: la parte del dominio que no cubre ninguna región declarada.
DEFAULT_REGION = "disc"

_RECT_KEYS = {"type", "r_inner", "r_outer", "half_thickness", "track", "axisymmetric", "regions"}
_PROFILE_KEYS = {"type", "points", "segments", "axisymmetric", "regions"}

Point = tuple[float, float]


# ------------------------------------------------------------------------- polígonos

def polygon_area(pts) -> float:
    """Área con signo (positiva si es antihorario), fórmula del cordón de zapato."""
    p = np.asarray(pts, dtype=float)
    r, z = p[:, 0], p[:, 1]
    return 0.5 * float(np.sum(r * np.roll(z, -1) - np.roll(r, -1) * z))


def polygon_r_moment(pts) -> float:
    """∫ r dA del polígono (antihorario) = (1/6) Σ (r_i + r_{i+1})(r_i z_{i+1} − r_{i+1} z_i).

    Es la referencia exacta de Σ A·r̄ de la malla (prueba T1 generalizada)."""
    p = np.asarray(pts, dtype=float)
    r, z = p[:, 0], p[:, 1]
    r1, z1 = np.roll(r, -1), np.roll(z, -1)
    return float(np.sum((r + r1) * (r * z1 - r1 * z)) / 6.0)


def _orient(a, b, c) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _on_segment(a, b, p, tol) -> bool:
    return (min(a[0], b[0]) - tol <= p[0] <= max(a[0], b[0]) + tol
            and min(a[1], b[1]) - tol <= p[1] <= max(a[1], b[1]) + tol)


def segments_intersect(a, b, c, d, tol: float = 0.0) -> bool:
    """¿Se tocan los segmentos [a, b] y [c, d]? (incluye contacto en un punto o solape)."""
    o1, o2, o3, o4 = _orient(a, b, c), _orient(a, b, d), _orient(c, d, a), _orient(c, d, b)
    if ((o1 > tol and o2 < -tol) or (o1 < -tol and o2 > tol)) and \
            ((o3 > tol and o4 < -tol) or (o3 < -tol and o4 > tol)):
        return True
    return any(abs(o) <= tol and _on_segment(s, t, p, tol) for o, s, t, p in
               ((o1, a, b, c), (o2, a, b, d), (o3, c, d, a), (o4, c, d, b)))


def points_in_polygon(points, poly) -> np.ndarray:
    """Máscara booleana: puntos estrictamente dentro del polígono (regla del rayo, vectorizada).
    Los puntos exactamente sobre el borde pueden caer de cualquier lado: úsese solo con puntos
    interiores (p. ej. centroides de triángulos de una malla conforme)."""
    P = np.asarray(points, dtype=float).reshape(-1, 2)
    V = np.asarray(poly, dtype=float)
    x, y = P[:, 0][:, None], P[:, 1][:, None]
    x0, y0 = V[:, 0][None, :], V[:, 1][None, :]
    x1, y1 = np.roll(V[:, 0], -1)[None, :], np.roll(V[:, 1], -1)[None, :]
    crosses = (y0 > y) != (y1 > y)
    with np.errstate(divide="ignore", invalid="ignore"):
        xin = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
    return (np.sum(crosses & (x < xin), axis=1) % 2) == 1


def point_segment_distance(points, a, b) -> np.ndarray:
    """Distancia de cada punto al segmento [a, b]."""
    P = np.asarray(points, dtype=float).reshape(-1, 2)
    a, b = np.asarray(a, float), np.asarray(b, float)
    ab = b - a
    t = np.clip(((P - a) @ ab) / float(ab @ ab), 0.0, 1.0)
    return np.linalg.norm(P - (a + t[:, None] * ab), axis=1)


def _validate_polygon(pts, what: str, scale: float) -> None:
    """Polígono simple, antihorario, sin lados nulos y con r ≥ 0."""
    p = np.asarray(pts, dtype=float)
    if p.ndim != 2 or p.shape[1] != 2 or p.shape[0] < 3:
        raise ValueError(f"{what}: se necesitan al menos 3 puntos [r, z].")
    if np.any(p[:, 0] < 0):
        raise ValueError(f"{what}: hay puntos con r < 0 (el eje de simetría es r = 0).")
    n = p.shape[0]
    lengths = np.linalg.norm(np.roll(p, -1, axis=0) - p, axis=1)
    if np.any(lengths <= 1e-12 * scale):
        raise ValueError(f"{what}: hay puntos repetidos consecutivos (lado de longitud nula).")
    area = polygon_area(p)
    if area <= 0:
        raise ValueError(
            f"{what}: el polígono debe recorrerse en sentido ANTIHORARIO (área con signo "
            f"{area:.3e} ≤ 0); invierta el orden de los puntos."
        )
    tol = 1e-12 * scale * scale
    for i in range(n):
        for j in range(i + 1, n):
            if j == i + 1 or (i == 0 and j == n - 1):
                continue  # lados adyacentes comparten un vértice
            if segments_intersect(p[i], p[(i + 1) % n], p[j], p[(j + 1) % n], tol):
                raise ValueError(
                    f"{what}: el polígono se autointersecta (lados {i} y {j}); debe ser simple."
                )


def _parse_regions(d, what="geometry.regions") -> dict[str, tuple[Point, ...]]:
    if d is None:
        return {}
    if not isinstance(d, dict):
        raise ValueError(f"'{what}' debe ser un diccionario {{nombre: [[r, z], ...]}}.")
    out = {}
    for name, pts in d.items():
        if name == DEFAULT_REGION:
            raise ValueError(
                f"'{DEFAULT_REGION}' es la región por defecto (lo no cubierto); use otro nombre."
            )
        try:
            out[str(name)] = tuple((float(r), float(z)) for r, z in pts)
        except (TypeError, ValueError):
            raise ValueError(f"La región '{name}' debe ser una lista de pares [r, z].") from None
    return out


class _SectionMixin:
    """Operaciones comunes a las dos secciones."""

    regions: dict
    axisymmetric: bool

    @property
    def region_names(self) -> list[str]:
        """``["disc", *regiones]``: índice de región de la malla → nombre."""
        return [DEFAULT_REGION, *self.regions]

    def boundary_segments(self) -> list[tuple[Point, Point, str]]:
        """Lados del borde (a, b, grupo) en sentido antihorario."""
        pts, segs = self.points, self.segments
        n = len(pts)
        return [(pts[k], pts[(k + 1) % n], segs[k]) for k in range(n)]

    def area(self) -> float:
        return polygon_area(self.points)

    def r_moment(self) -> float:
        """∫ r dA de la sección."""
        return polygon_r_moment(self.points)

    @property
    def z_brake(self) -> float:
        return next(a[1] for a, _, g in self.boundary_segments() if g == "brake")

    def _validate_regions(self) -> None:
        scale = float(np.ptp(np.asarray(self.points)[:, 0])) or 1.0
        outer = np.asarray(self.points, float)
        for name, poly in self.regions.items():
            _validate_polygon(poly, f"Región '{name}'", scale)
            P = np.asarray(poly, float)
            inside = points_in_polygon(P, outer)
            on_edge = np.zeros(len(P), dtype=bool)
            for a, b, _ in self.boundary_segments():
                on_edge |= point_segment_distance(P, a, b) <= 1e-9 * scale
            if not np.all(inside | on_edge):
                raise ValueError(f"La región '{name}' tiene vértices fuera de la sección.")


@dataclass(frozen=True)
class RectangleGeometry(_SectionMixin):
    """Sección rectangular del disco (medio espesor), en metros.

    Atributos:
        r_inner: radio interior r_i (≥ 0).
        r_outer: radio exterior r_o.
        half_thickness: semiespesor e; z = 0 es el plano de simetría.
        r_a, r_b: radios interior y exterior de la pista de fricción, r_i < r_a < r_b < r_o.
        axisymmetric: True ⇒ peso r en las integrales; False ⇒ peso 1 (solo validación).
        regions: sub-polígonos con otro material, {nombre: ((r, z), ...)}.
    """

    r_inner: float
    r_outer: float
    half_thickness: float
    r_a: float
    r_b: float
    axisymmetric: bool = True
    regions: dict = field(default_factory=dict)

    type = "rectangle"

    def __post_init__(self) -> None:
        if self.r_inner < 0:
            raise ValueError(f"geometry.r_inner debe ser ≥ 0; se recibió {self.r_inner}.")
        if not (self.r_inner < self.r_a < self.r_b < self.r_outer):
            raise ValueError(
                "La geometría debe cumplir r_inner < r_a < r_b < r_outer; se recibió "
                f"r_inner={self.r_inner}, r_a={self.r_a}, r_b={self.r_b}, "
                f"r_outer={self.r_outer}."
            )
        if self.half_thickness <= 0:
            raise ValueError(
                f"geometry.half_thickness debe ser > 0; se recibió {self.half_thickness}."
            )
        self._validate_regions()

    @property
    def corners(self) -> list[Point]:
        """Vértices del borde en sentido antihorario: (r_i,0), (r_o,0), (r_o,e), (r_b,e),
        (r_a,e), (r_i,e). La arista k va del vértice k al k+1 y es el grupo
        ``BOUNDARY_GROUPS[k]``."""
        e = self.half_thickness
        return [
            (self.r_inner, 0.0),
            (self.r_outer, 0.0),
            (self.r_outer, e),
            (self.r_b, e),
            (self.r_a, e),
            (self.r_inner, e),
        ]

    @property
    def points(self) -> list[Point]:
        return self.corners

    @property
    def segments(self) -> tuple[str, ...]:
        return BOUNDARY_GROUPS

    @classmethod
    def from_dict(cls, d: dict) -> RectangleGeometry:
        """Construye la geometría a partir de la sección ``geometry`` del YAML."""
        unknown = set(d) - _RECT_KEYS
        if unknown:
            raise ValueError(f"Claves desconocidas en 'geometry': {sorted(unknown)}.")
        if d.get("type", "rectangle") != "rectangle":
            raise ValueError(f"RectangleGeometry no admite geometry.type = '{d.get('type')}'.")
        try:
            track = d["track"]
            return cls(
                r_inner=float(d["r_inner"]),
                r_outer=float(d["r_outer"]),
                half_thickness=float(d["half_thickness"]),
                r_a=float(track["r_a"]),
                r_b=float(track["r_b"]),
                axisymmetric=bool(d.get("axisymmetric", True)),
                regions=_parse_regions(d.get("regions")),
            )
        except KeyError as exc:
            raise ValueError(f"Falta la clave {exc} en la sección 'geometry'.") from None


@dataclass(frozen=True)
class ProfileGeometry(_SectionMixin):
    """Sección poligonal (r, z) antihoraria con un grupo de frontera por lado.

    Atributos:
        points: vértices ((r, z), ...) en sentido antihorario.
        segments: grupo de cada lado; el lado k va de ``points[k]`` a ``points[k+1]``.
        axisymmetric: peso r (True) o 1 (False).
        regions: sub-polígonos con otro material.
    """

    points: tuple[Point, ...]
    segments: tuple[str, ...]
    axisymmetric: bool = True
    regions: dict = field(default_factory=dict)

    type = "profile"

    def __post_init__(self) -> None:
        n = len(self.points)
        if len(self.segments) != n:
            raise ValueError(
                f"geometry.segments debe tener un grupo por lado: {n} puntos ⇒ {n} lados, "
                f"pero se dieron {len(self.segments)} grupos."
            )
        bad = sorted(set(self.segments) - set(BOUNDARY_GROUPS))
        if bad:
            raise ValueError(f"Grupos desconocidos en geometry.segments: {bad}; "
                             f"opciones: {list(BOUNDARY_GROUPS)}.")
        p = np.asarray(self.points, float)
        scale = float(np.ptp(p[:, 0])) or 1.0
        _validate_polygon(p, "geometry.points", scale)
        brake = [k for k, g in enumerate(self.segments) if g == "brake"]
        if not brake:
            raise ValueError("El perfil necesita al menos un lado 'brake' (la pista).")
        # Contigüidad (cíclica) de los lados brake.
        starts = [k for k in brake if (k - 1) % n not in brake]
        if len(starts) != 1:
            raise ValueError("Los lados 'brake' deben formar una sola cadena contigua.")
        zs = {p[k, 1] for k in brake} | {p[(k + 1) % n, 1] for k in brake}
        if len(zs) != 1:
            raise ValueError(
                "Los lados 'brake' deben ser horizontales y estar en el mismo z (la pastilla "
                "apoya en un plano z = cte; q(r) se normaliza sobre [r_a, r_b])."
            )
        self._validate_regions()

    @property
    def _brake_r(self) -> np.ndarray:
        n = len(self.points)
        return np.array([self.points[j][0] for k, g in enumerate(self.segments) if g == "brake"
                         for j in (k, (k + 1) % n)])

    @property
    def r_a(self) -> float:
        """Radio interior de la pista (mínimo r de los lados ``brake``)."""
        return float(self._brake_r.min())

    @property
    def r_b(self) -> float:
        """Radio exterior de la pista (máximo r de los lados ``brake``)."""
        return float(self._brake_r.max())

    @property
    def corners(self) -> list[Point]:
        return list(self.points)

    @classmethod
    def from_dict(cls, d: dict) -> ProfileGeometry:
        unknown = set(d) - _PROFILE_KEYS
        if unknown:
            raise ValueError(f"Claves desconocidas en 'geometry' (profile): {sorted(unknown)}.")
        if "points" not in d or "segments" not in d:
            raise ValueError("geometry.type = profile necesita 'points' y 'segments'.")
        try:
            pts = tuple((float(r), float(z)) for r, z in d["points"])
        except (TypeError, ValueError):
            raise ValueError("geometry.points debe ser una lista de pares [r, z].") from None
        return cls(
            points=pts,
            segments=tuple(str(s) for s in d["segments"]),
            axisymmetric=bool(d.get("axisymmetric", True)),
            regions=_parse_regions(d.get("regions")),
        )


Geometry = RectangleGeometry | ProfileGeometry


def geometry_from_dict(d: dict) -> Geometry:
    """Construye la sección según ``geometry.type`` (``rectangle`` por defecto o ``profile``)."""
    if not isinstance(d, dict):
        raise ValueError("La sección 'geometry' debe ser un diccionario.")
    gtype = d.get("type", "rectangle")
    if gtype == "rectangle":
        return RectangleGeometry.from_dict(d)
    if gtype == "profile":
        return ProfileGeometry.from_dict(d)
    raise ValueError(f"geometry.type debe ser 'rectangle' o 'profile'; se recibió '{gtype}'.")
