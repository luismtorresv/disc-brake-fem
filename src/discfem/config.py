"""Carga y validación del YAML de una corrida (PLAN.md §4).

``load_config(path, overrides)`` lee el YAML, aplica ``overrides`` con claves con puntos
(p. ej. ``{"mesh.refine_factor": 2}``) y valida todas las secciones con mensajes en español.
Las propiedades de materiales NO se leen aquí: la librería se carga al empezar la corrida
(``solver.run``), donde también se emiten las advertencias de materiales no verificados.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .geometry import BOUNDARY_GROUPS, RectangleGeometry
from .mesh import MeshOptions
from .timestepping import DEFAULT_RANNACHER_STEPS, DEFAULT_THETA

_TOP_KEYS = {
    "case", "geometry", "mesh", "materials", "boundaries", "convection", "initial", "load",
    "time", "output",
}
_BC_TYPES = ("convection", "insulated", "flux")
_FORMATS = ("csv", "vtk", "png")


@dataclass(frozen=True)
class MaterialsConfig:
    """Sección ``materials``: librería, material por región, pastilla y overrides."""

    library: str
    regions: dict[str, str]
    pad: str | None = None
    overrides: dict[str, dict] = field(default_factory=dict)


@dataclass(frozen=True)
class TimeConfig:
    """Sección ``time``. ``dt`` es ``"auto"`` (t_f/100) o un número en segundos."""

    t_end: float
    theta: float = DEFAULT_THETA
    rannacher_steps: int = DEFAULT_RANNACHER_STEPS
    dt: float | str = "auto"


@dataclass(frozen=True)
class OutputConfig:
    """Sección ``output``: carpeta, sondas (r, z) y formatos (csv, vtk, png)."""

    dir: str
    probes: tuple[tuple[float, float], ...] = ()
    formats: tuple[str, ...] = _FORMATS


@dataclass(frozen=True)
class Config:
    """Configuración validada de una corrida.

    ``raw`` es el diccionario efectivo (YAML + overrides): se escribe como ``config.yaml`` en la
    carpeta de resultados para que la corrida pueda reproducirse solo con él.
    """

    case: str
    geometry: RectangleGeometry
    mesh: MeshOptions
    materials: MaterialsConfig
    boundaries: dict[str, str]
    h: float
    T_inf: float
    T0: float
    load: dict
    time: TimeConfig
    output: OutputConfig
    raw: dict
    config_dir: Path

    def groups_with(self, bc_type: str) -> list[str]:
        """Grupos de frontera con la condición ``bc_type``, en el orden de BOUNDARY_GROUPS."""
        return [g for g in BOUNDARY_GROUPS if self.boundaries[g] == bc_type]


def set_dotted(d: dict, key: str, value: Any) -> None:
    """``d["a"]["b"]["c"] = value`` para ``key = "a.b.c"``; crea los diccionarios intermedios."""
    parts = key.split(".")
    if not all(parts):
        raise ValueError(f"Clave de override inválida: '{key}'.")
    node = d
    for p in parts[:-1]:
        if p not in node or node[p] is None:
            node[p] = {}
        if not isinstance(node[p], dict):
            raise ValueError(f"Override '{key}': '{p}' no es una sección (es {node[p]!r}).")
        node = node[p]
    node[parts[-1]] = value


def _section(raw: dict, name: str, required: bool = True) -> dict:
    sec = raw.get(name)
    if sec is None:
        if required:
            raise ValueError(f"Falta la sección '{name}' en la configuración.")
        return {}
    if not isinstance(sec, dict):
        raise ValueError(f"La sección '{name}' debe ser un diccionario.")
    return sec


def _check_keys(sec: dict, name: str, allowed: set[str]) -> None:
    unknown = set(sec) - allowed
    if unknown:
        raise ValueError(f"Claves desconocidas en '{name}': {sorted(unknown)}.")


def _number(sec: dict, name: str, key: str, default: Any = None, positive: bool = False) -> float:
    if key not in sec:
        if default is None:
            raise ValueError(f"Falta '{name}.{key}'.")
        return default
    try:
        v = float(sec[key])
    except (TypeError, ValueError):
        raise ValueError(f"'{name}.{key}' debe ser un número; se recibió {sec[key]!r}.") from None
    if positive and not v > 0:
        raise ValueError(f"'{name}.{key}' debe ser > 0; se recibió {v}.")
    return v


def build_config(raw: dict, config_dir: str | Path = ".") -> Config:
    """Valida el diccionario ``raw`` (ya con overrides) y construye :class:`Config`."""
    if not isinstance(raw, dict):
        raise ValueError("La configuración debe ser un diccionario YAML.")
    _check_keys(raw, "raíz", _TOP_KEYS)
    case = str(raw.get("case", "caso"))

    geometry = RectangleGeometry.from_dict(_section(raw, "geometry"))
    mesh = MeshOptions.from_dict(_section(raw, "mesh", required=False))

    mat = _section(raw, "materials")
    _check_keys(mat, "materials", {"library", "regions", "pad", "overrides"})
    if "library" not in mat or "regions" not in mat:
        raise ValueError("'materials' necesita 'library' y 'regions'.")
    regions = mat["regions"]
    if not isinstance(regions, dict) or not regions:
        raise ValueError("'materials.regions' debe ser un diccionario {región: material}.")
    if set(regions) != {"disc"}:
        raise ValueError(
            f"La geometría rectangular tiene una sola región, 'disc'; se recibió {sorted(regions)}."
        )
    overrides = mat.get("overrides") or {}
    if not isinstance(overrides, dict):
        raise ValueError("'materials.overrides' debe ser un diccionario {región: {k: …}}.")
    materials = MaterialsConfig(
        library=str(mat["library"]), regions=dict(regions), pad=mat.get("pad"),
        overrides=dict(overrides),
    )

    bnd = _section(raw, "boundaries")
    _check_keys(bnd, "boundaries", set(BOUNDARY_GROUPS))
    missing = [g for g in BOUNDARY_GROUPS if g not in bnd]
    if missing:
        raise ValueError(f"Falta la condición de los grupos {missing} en 'boundaries'.")
    for g, bc in bnd.items():
        if bc not in _BC_TYPES:
            raise ValueError(f"boundaries.{g} = {bc!r}; debe ser uno de {_BC_TYPES}.")

    has_conv = any(v == "convection" for v in bnd.values())
    conv = _section(raw, "convection", required=has_conv)
    _check_keys(conv, "convection", {"h", "T_inf"})
    if has_conv:
        h = _number(conv, "convection", "h", positive=True)
        T_inf = _number(conv, "convection", "T_inf")
    else:
        h, T_inf = 0.0, 0.0

    init = _section(raw, "initial")
    _check_keys(init, "initial", {"T0"})
    T0 = _number(init, "initial", "T0")

    load = _section(raw, "load", required=any(v == "flux" for v in bnd.values()))
    if load and load.get("type") == "braking" and materials.pad is None \
            and load.get("gamma", "auto") == "auto":
        raise ValueError("gamma: auto necesita 'materials.pad' para calcular γ.")

    tm = _section(raw, "time")
    _check_keys(tm, "time", {"theta", "rannacher_steps", "dt", "t_end"})
    theta = _number(tm, "time", "theta", default=DEFAULT_THETA)
    if not 0 <= theta <= 1:
        raise ValueError(f"time.theta debe estar en [0, 1]; se recibió {theta}.")
    rs = tm.get("rannacher_steps", DEFAULT_RANNACHER_STEPS)
    if not isinstance(rs, int) or rs < 0:
        raise ValueError(f"time.rannacher_steps debe ser un entero ≥ 0; se recibió {rs!r}.")
    dt = tm.get("dt", "auto")
    if dt != "auto":
        dt = _number(tm, "time", "dt", positive=True)
    time = TimeConfig(t_end=_number(tm, "time", "t_end", positive=True), theta=theta,
                      rannacher_steps=rs, dt=dt)

    out = _section(raw, "output", required=False)
    _check_keys(out, "output", {"dir", "probes", "formats"})
    probes = out.get("probes") or []
    try:
        probes = tuple((float(p[0]), float(p[1])) for p in probes if len(p) == 2)
    except (TypeError, ValueError, IndexError):
        raise ValueError("output.probes debe ser una lista de pares [r, z].") from None
    if len(probes) != len(out.get("probes") or []):
        raise ValueError("output.probes debe ser una lista de pares [r, z].")
    formats = tuple(out.get("formats", _FORMATS) or ())
    bad = [f for f in formats if f not in _FORMATS]
    if bad:
        raise ValueError(f"output.formats contiene {bad}; opciones: {_FORMATS}.")
    output = OutputConfig(dir=str(out.get("dir", f"results/{case}")), probes=probes,
                          formats=formats)

    return Config(
        case=case, geometry=geometry, mesh=mesh, materials=materials, boundaries=dict(bnd),
        h=h, T_inf=T_inf, T0=T0, load=dict(load), time=time, output=output,
        raw=raw, config_dir=Path(config_dir),
    )


def load_config(
    source: str | Path | dict, overrides: dict[str, Any] | None = None
) -> Config:
    """Lee el YAML (o toma un diccionario), aplica ``overrides`` con claves con puntos y valida.

    Args:
        source: ruta al YAML o diccionario ya leído.
        overrides: ``{"sección.clave": valor}``; p. ej. ``{"mesh.refine_factor": 2}``.
    """
    if isinstance(source, dict):
        raw, config_dir = copy.deepcopy(source), Path(".")
    else:
        path = Path(source)
        if not path.is_file():
            raise ValueError(f"No existe el archivo de configuración '{path}'.")
        with path.open(encoding="utf-8") as fh:
            raw = yaml.safe_load(fh)
        config_dir = path.parent
    if not isinstance(raw, dict):
        raise ValueError("La configuración debe ser un diccionario YAML.")
    for key, value in (overrides or {}).items():
        set_dotted(raw, key, value)
    return build_config(raw, config_dir)
