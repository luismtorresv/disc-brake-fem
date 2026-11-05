"""Librería de materiales, efusividad térmica y partición de calor γ (PLAN.md §3.4 y §4).

Reglas (CLAUDE.md, regla 4): ninguna propiedad se inventa en el código. Toda entrada de la
librería debe traer ``source`` (texto no vacío) y ``verified`` (booleano). Si falta alguno, la
carga falla con un mensaje en español. Mientras haya materiales con ``verified: false``, cada
corrida debe emitir :class:`UnverifiedMaterialWarning` con la lista (ver :func:`warn_unverified`).

Partición del calor de fricción entre disco y pastilla, por efusividades e = √(k ρ c_p):

    γ = e_disco / (e_disco + e_pastilla)
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import yaml

_PROPERTY_KEYS = ("k", "rho", "cp")
_ALLOWED_KEYS = {"description", "k", "rho", "cp", "source", "verified"}

#: Rangos de sanidad (PLAN.md §4; aproximados, "de memoria; verificar"). Fuera de ellos se
#: advierte, sin abortar.
SANITY_RANGES = {
    "k": (0.01, 500.0, "W/(m·K)"),
    "rho": (100.0, 20_000.0, "kg/m³"),
    "cp": (100.0, 5_000.0, "J/(kg·K)"),
    "h": (1.0, 500.0, "W/(m²·K)"),
}


class UnverifiedMaterialWarning(UserWarning):
    """Hay materiales con ``verified: false``: los resultados no sirven para el informe."""


class PropertyRangeWarning(UserWarning):
    """Una propiedad está fuera del rango de sanidad de ``SANITY_RANGES``."""


@dataclass(frozen=True)
class Material:
    """Material homogéneo e isótropo con propiedades constantes (SI).

    Atributos:
        name: clave en la librería.
        k: conductividad [W/(m·K)].
        rho: densidad [kg/m³].
        cp: calor específico [J/(kg·K)].
        source: referencia bibliográfica (obligatoria).
        verified: True solo si los valores se contrastaron con ``source``.
        description: texto libre.
    """

    name: str
    k: float
    rho: float
    cp: float
    source: str
    verified: bool
    description: str = ""

    @property
    def rho_cp(self) -> float:
        """Capacidad térmica volumétrica ρ c_p [J/(m³·K)]."""
        return self.rho * self.cp

    @property
    def diffusivity(self) -> float:
        """Difusividad α = k / (ρ c_p) [m²/s]."""
        return self.k / self.rho_cp

    @property
    def effusivity(self) -> float:
        """Efusividad e = √(k ρ c_p) [W·s^½/(m²·K)]."""
        return math.sqrt(self.k * self.rho_cp)


def check_range(name: str, value: float, context: str = "") -> None:
    """Emite :class:`PropertyRangeWarning` si ``value`` está fuera de ``SANITY_RANGES[name]``."""
    lo, hi, unit = SANITY_RANGES[name]
    if not lo < value < hi:
        where = f" ({context})" if context else ""
        warnings.warn(
            f"{name} = {value:g} {unit}{where} está fuera del rango de sanidad "
            f"({lo:g}, {hi:g}) {unit}. Revise el valor.",
            PropertyRangeWarning,
            stacklevel=2,
        )


def _material_from_entry(name: str, entry: dict) -> Material:
    """Valida una entrada de la librería y construye el :class:`Material`."""
    if not isinstance(entry, dict):
        raise ValueError(  # noqa: TRY004 (error de configuración: ValueError uniforme)
            f"El material '{name}' debe ser un diccionario de propiedades."
        )
    unknown = set(entry) - _ALLOWED_KEYS
    if unknown:
        raise ValueError(f"Claves desconocidas en el material '{name}': {sorted(unknown)}.")
    for key in (*_PROPERTY_KEYS, "source", "verified"):
        if key not in entry:
            raise ValueError(
                f"Al material '{name}' le falta la clave obligatoria '{key}'. Toda entrada de la "
                "librería debe tener k, rho, cp, source y verified (CLAUDE.md, regla 4)."
            )
    source = entry["source"]
    if not isinstance(source, str) or not source.strip():
        raise ValueError(f"El material '{name}' tiene 'source' vacío: cite la fuente.")
    if not isinstance(entry["verified"], bool):
        raise ValueError(  # noqa: TRY004 (error de configuración: ValueError uniforme)
            f"En el material '{name}', 'verified' debe ser true o false."
        )
    props = {}
    for key in _PROPERTY_KEYS:
        try:
            props[key] = float(entry[key])
        except (TypeError, ValueError):
            raise ValueError(
                f"En el material '{name}', '{key}' debe ser un número; se recibió {entry[key]!r}."
            ) from None
        if not props[key] > 0:
            raise ValueError(f"En el material '{name}', '{key}' debe ser > 0.")
        check_range(key, props[key], f"material '{name}'")
    return Material(
        name=name,
        source=source,
        verified=entry["verified"],
        description=str(entry.get("description", "")),
        **props,
    )


def load_library(path: str | Path) -> dict[str, Material]:
    """Lee y valida la librería de materiales (YAML ``{nombre: {k, rho, cp, source, verified}}``)."""
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"No existe la librería de materiales '{path}'.")
    with path.open(encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict) or not data:
        raise ValueError(f"La librería '{path}' debe ser un diccionario no vacío de materiales.")
    return {name: _material_from_entry(name, entry) for name, entry in data.items()}


def resolve_library_path(library: str | Path, config_dir: str | Path | None = None) -> Path:
    """Ruta de la librería: tal cual si existe (relativa al directorio actual) y, si no,
    relativa al directorio del archivo de configuración."""
    p = Path(library)
    if p.is_absolute() or p.is_file() or config_dir is None:
        return p
    return Path(config_dir) / p


def apply_override(base: Material, override: dict) -> Material:
    """Material con algunas propiedades reemplazadas desde ``materials.overrides``.

    Los valores sobrescritos vienen del YAML de la corrida, no de la librería: el resultado
    queda con ``verified = False`` y una ``source`` que lo dice, salvo que el override traiga
    sus propios ``source`` y ``verified``.
    """
    unknown = set(override) - {"k", "rho", "cp", "source", "verified"}
    if unknown:
        raise ValueError(
            f"Claves no permitidas en el override de '{base.name}': {sorted(unknown)}."
        )
    changes = {k: float(v) for k, v in override.items() if k in _PROPERTY_KEYS}
    for key, value in changes.items():
        if not value > 0:
            raise ValueError(f"En el override de '{base.name}', '{key}' debe ser > 0.")
        check_range(key, value, f"override de '{base.name}'")
    desc = ", ".join(f"{k}={v:g}" for k, v in changes.items())
    source = override.get("source", f"override en la config ({desc}) sobre: {base.source}")
    verified = bool(override.get("verified", False))
    return replace(base, source=source, verified=verified, **changes)


def resolve_region_materials(
    regions: dict[str, str],
    library: dict[str, Material],
    overrides: dict[str, dict] | None = None,
) -> dict[str, Material]:
    """Asigna un :class:`Material` a cada región, aplicando ``overrides`` por región."""
    out = {}
    for region, mat_name in regions.items():
        if mat_name not in library:
            raise ValueError(
                f"La región '{region}' usa el material '{mat_name}', que no está en la "
                f"librería (disponibles: {sorted(library)})."
            )
        out[region] = library[mat_name]
    for region, ov in (overrides or {}).items():
        if region not in out:
            raise ValueError(f"Override para la región '{region}', que no existe en 'regions'.")
        out[region] = apply_override(out[region], ov or {})
    return out


def unverified(materials: list[Material] | dict[str, Material]) -> list[Material]:
    """Materiales con ``verified: false`` (sin repetir), en orden de aparición."""
    items = materials.values() if isinstance(materials, dict) else materials
    seen, out = set(), []
    for m in items:
        key = (m.name, m.source)
        if not m.verified and key not in seen:
            seen.add(key)
            out.append(m)
    return out


def warn_unverified(materials: list[Material] | dict[str, Material]) -> list[str]:
    """Emite :class:`UnverifiedMaterialWarning` con la lista de materiales no verificados y
    la devuelve (nombres). Debe llamarse al inicio de cada corrida."""
    bad = unverified(materials)
    names = [m.name for m in bad]
    if bad:
        detail = "; ".join(f"{m.name} (source: {m.source})" for m in bad)
        warnings.warn(
            f"Materiales NO verificados: {detail}. Los resultados no sirven para el informe "
            "hasta citar fuentes y poner verified: true.",
            UnverifiedMaterialWarning,
            stacklevel=2,
        )
    return names


def partition_coefficient(disc: Material, pad: Material) -> float:
    """γ = e_disco / (e_disco + e_pastilla), fracción del calor de fricción que entra al disco."""
    return disc.effusivity / (disc.effusivity + pad.effusivity)


def brake_regions(mesh) -> list[str]:
    """Nombres de las regiones de los triángulos que tienen una arista en ``brake``."""
    brake = np.sort(mesh.edges["brake"], axis=1)
    n = mesh.n_nodes
    bkeys = brake[:, 0] * n + brake[:, 1]
    t = mesh.tris
    tri_keys = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        e = np.sort(t[:, [a, b]], axis=1)
        tri_keys.append(e[:, 0] * n + e[:, 1])
    touches = np.isin(np.stack(tri_keys, axis=1), bkeys).any(axis=1)
    return sorted({mesh.region_names[r] for r in np.unique(mesh.region[touches])})


def resolve_gamma(
    gamma: str | float,
    mesh=None,
    region_materials: dict[str, Material] | None = None,
    pad: Material | None = None,
) -> float:
    """γ de la config: un número en (0, 1), o ``"auto"`` desde el material de la región
    adyacente a ``brake`` y el de la pastilla."""
    if gamma == "auto":
        if mesh is None or region_materials is None or pad is None:
            raise ValueError("gamma: auto requiere la malla, los materiales por región y 'pad'.")
        regions = brake_regions(mesh)
        mats = {region_materials[r] for r in regions}
        if len(mats) != 1:
            raise ValueError(
                "gamma: auto requiere un solo material bajo 'brake'; las regiones adyacentes "
                f"{regions} tienen materiales distintos. Indique gamma numérico."
            )
        return partition_coefficient(mats.pop(), pad)
    try:
        g = float(gamma)
    except (TypeError, ValueError):
        raise ValueError(f"gamma debe ser 'auto' o un número en (0, 1); se recibió {gamma!r}.")
    if not 0.0 < g < 1.0:
        raise ValueError(f"gamma debe estar en (0, 1); se recibió {g}.")
    return g
