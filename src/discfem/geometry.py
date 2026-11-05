"""Definición de la sección axisimétrica (r, z) del disco y de sus grupos de frontera.

Por ahora solo existe la sección rectangular Ω = (r_i, r_o) × (0, e), con z = 0 como plano de
simetría y la pista de fricción `brake` en z = e, r_a < r < r_b (ver CLAUDE.md y PLAN.md §3.1).
El perfil poligonal (`type: profile`) queda para la Fase 7.
"""

from __future__ import annotations

from dataclasses import dataclass

#: Nombres de los grupos de frontera, en el orden en que se recorre el borde en sentido
#: antihorario empezando por el plano de simetría. Cada arista de borde pertenece a uno solo.
BOUNDARY_GROUPS: tuple[str, ...] = (
    "symmetry",   # z = 0
    "r_outer",    # r = r_o
    "top_outer",  # z = e, r_b..r_o
    "brake",      # z = e, r_a..r_b
    "top_inner",  # z = e, r_i..r_a
    "r_inner",    # r = r_i
)

_GEOMETRY_KEYS = {"type", "r_inner", "r_outer", "half_thickness", "track", "axisymmetric"}


@dataclass(frozen=True)
class RectangleGeometry:
    """Sección rectangular del disco (medio espesor), en metros.

    Atributos:
        r_inner: radio interior r_i (≥ 0).
        r_outer: radio exterior r_o.
        half_thickness: semiespesor e; z = 0 es el plano de simetría.
        r_a, r_b: radios interior y exterior de la pista de fricción, r_i < r_a < r_b < r_o.
        axisymmetric: True ⇒ peso r en las integrales; False ⇒ peso 1 (solo validación).
    """

    r_inner: float
    r_outer: float
    half_thickness: float
    r_a: float
    r_b: float
    axisymmetric: bool = True

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

    @property
    def corners(self) -> list[tuple[float, float]]:
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

    @classmethod
    def from_dict(cls, d: dict) -> RectangleGeometry:
        """Construye la geometría a partir de la sección ``geometry`` del YAML."""
        unknown = set(d) - _GEOMETRY_KEYS
        if unknown:
            raise ValueError(f"Claves desconocidas en 'geometry': {sorted(unknown)}.")
        gtype = d.get("type", "rectangle")
        if gtype != "rectangle":
            raise ValueError(
                f"geometry.type = '{gtype}' no está implementado (solo 'rectangle'; "
                "'profile' llega en la Fase 7)."
            )
        try:
            track = d["track"]
            return cls(
                r_inner=float(d["r_inner"]),
                r_outer=float(d["r_outer"]),
                half_thickness=float(d["half_thickness"]),
                r_a=float(track["r_a"]),
                r_b=float(track["r_b"]),
                axisymmetric=bool(d.get("axisymmetric", True)),
            )
        except KeyError as exc:
            raise ValueError(f"Falta la clave {exc} en la sección 'geometry'.") from None
