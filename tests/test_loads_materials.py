"""Pruebas de la Fase 4: T3 (PLAN.md §6) para el flujo de fricción, y la librería de materiales.

T3: 2π ∫_{r_a}^{r_b} q_fric r dr = γ P_f(t) para ambos perfiles, y coincidencia con las fórmulas
cerradas de la §3.4, a 1e-12 relativo.
"""

from __future__ import annotations

import math
import re
import warnings
from pathlib import Path

import numpy as np
import pytest
import yaml

from discfem.geometry import RectangleGeometry
from discfem.loads import (
    BrakingLoad,
    ConstantFluxLoad,
    TableLoad,
    brake_load_vector,
    make_load,
)
from discfem.materials import (
    SOURCE_MISSING,
    Material,
    PropertyRangeWarning,
    UnverifiedMaterialWarning,
    apply_override,
    brake_regions,
    check_range,
    load_library,
    partition_coefficient,
    resolve_gamma,
    resolve_region_materials,
    warn_unverified,
)
from discfem.mesh import MeshOptions, build_mesh

ROOT = Path(__file__).resolve().parents[1]
LIBRARY = ROOT / "configs" / "materiales.yaml"
REF_CFG = ROOT / "configs" / "ejemplo_frenado.yaml"
RTOL = 1e-12
GEO = RectangleGeometry(r_inner=0.07, r_outer=0.15, half_thickness=0.01, r_a=0.08, r_b=0.14)
VEHICLE = {"mass": 1500.0, "v0": 27.78, "decel": 5.0, "wheel_radius": 0.30,
           "brake_fraction": 0.7, "discs_per_axle": 2, "faces_per_disc": 2}
GAMMA_REF = 0.869  # PLAN.md §6, con los materiales supuestos


def braking(profile: str, gamma: float = GAMMA_REF) -> BrakingLoad:
    cfg = {"type": "braking", "vehicle": VEHICLE, "pressure_profile": profile}
    return make_load(cfg, GEO.r_a, GEO.r_b, gamma)


def radial_integral(f, a: float, b: float, n: int = 20) -> float:
    """∫_a^b f(r) dr con Gauss–Legendre de n puntos (exacto para polinomios de grado < 2n)."""
    x, w = np.polynomial.legendre.leggauss(n)
    r = 0.5 * (b - a) * x + 0.5 * (a + b)
    return 0.5 * (b - a) * np.sum(w * f(r))


TIMES = [0.0, 1.0, 2.5, 5.0, 27.78 / 5.0]


# -------------------------------------------------------------------------------- T3

@pytest.mark.parametrize("profile", ["uniform_pressure", "uniform_wear"])
@pytest.mark.parametrize("t", TIMES)
def test_t3_power_consistency(profile, t):
    """2π ∫ q r dr = γ P_f(t)."""
    ld = braking(profile)
    lhs = 2 * math.pi * radial_integral(lambda r: ld.q(r, t) * r, GEO.r_a, GEO.r_b)
    rhs = ld.gamma * ld.power_per_face(t)
    if rhs == 0:
        assert lhs == 0
    else:
        assert abs(lhs - rhs) / rhs < RTOL


@pytest.mark.parametrize("t", TIMES[:-1])
def test_t3_closed_forms(t):
    """Fórmulas cerradas de la §3.4 para ambos perfiles."""
    m, a, phi, nd, nc = 1500.0, 5.0, 0.7, 2, 2
    v = 27.78 - a * t
    r = np.linspace(GEO.r_a, GEO.r_b, 7)
    ra, rb, g = GEO.r_a, GEO.r_b, GAMMA_REF
    q_p = 3 * g * phi * m * a * v * r / (2 * math.pi * nd * nc * (rb**3 - ra**3))
    q_w = g * phi * m * a * v / (math.pi * nd * nc * (rb**2 - ra**2))
    np.testing.assert_allclose(braking("uniform_pressure").q(r, t), q_p, rtol=RTOL)
    np.testing.assert_allclose(braking("uniform_wear").q(r, t), q_w, rtol=RTOL)


@pytest.mark.parametrize("profile", ["uniform_pressure", "uniform_wear"])
@pytest.mark.parametrize("backend", ["structured", "gmsh"])
def test_t3_fem_load_vector_power(profile, backend):
    """Con la carga ensamblada: 2π Σ_i F_brake,i(t) = γ P_f(t) (Gauss de 3 puntos es exacto
    para ambos perfiles y los nodos r_a, r_b están en la malla)."""
    mesh = build_mesh(GEO, MeshOptions(backend=backend, size_track=1e-3, size_far=4e-3))
    ld = braking(profile)
    F = brake_load_vector(mesh, ld, axisymmetric=True)
    for t in (0.0, 2.0, 4.0):
        assert abs(2 * math.pi * F(t).sum() - ld.gamma * ld.power_per_face(t)) \
            / (ld.gamma * ld.power_per_face(t)) < RTOL


def test_braking_reference_numbers():
    """Valores de la §6 (T8/T9): t_f = 5,56 s; q0 = 764 kW/m²; energía por cara = 101 kJ."""
    ld = braking("uniform_wear", gamma=0.8692006698254295)
    assert ld.t_f == pytest.approx(5.556, rel=1e-12)
    assert ld.q(0.11, 0.0) == pytest.approx(764e3, rel=1e-3)
    assert ld.energy_per_face == pytest.approx(101e3, rel=3e-3)
    # ∫ P_f dt = energía por cara (P_f lineal en t: trapecio exacto).
    t = np.linspace(0, ld.t_f, 11)
    assert np.trapezoid(ld.power_per_face(t), t) == pytest.approx(ld.energy_per_face, rel=RTOL)


def test_braking_zero_after_tf_and_kinematics():
    ld = braking("uniform_wear")
    assert ld.time_factor(ld.t_f) == 0.0
    assert np.all(ld.q(np.array([0.09, 0.12]), ld.t_f + 1.0) == 0.0)
    assert ld.time_factor(0.0) == 1.0
    assert ld.angular_velocity(0.0) == pytest.approx(27.78 / 0.30)
    assert np.all(ld.shape(np.array([0.07, 0.145])) == 0.0)  # fuera de la pista


def test_constant_and_table_loads():
    c = make_load({"type": "constant_flux", "q": 2e5, "duration": 3.0}, GEO.r_a, GEO.r_b)
    assert isinstance(c, ConstantFluxLoad) and c.t_f == 3.0
    np.testing.assert_array_equal(c.q(np.array([0.09, 0.1]), 3.0), [2e5, 2e5])
    assert c.q(0.09, 3.0001) == 0.0
    c_inf = make_load({"type": "constant_flux", "q": 1e5}, GEO.r_a, GEO.r_b)
    assert c_inf.t_f is None and c_inf.q(0.1, 1e6) == 1e5

    tb = make_load({"type": "table", "points": [[0, 0], [1, 1e5], [3, 5e4]]}, GEO.r_a, GEO.r_b)
    assert isinstance(tb, TableLoad) and tb.t_f == 3.0
    np.testing.assert_allclose(tb.q(0.1, np.array([0.5, 1.0, 2.0, 3.0, 3.5])),
                               [5e4, 1e5, 7.5e4, 5e4, 0.0], rtol=RTOL)


@pytest.mark.parametrize(
    ("cfg", "match"),
    [
        ({"type": "braking", "vehicle": {**VEHICLE, "decel": 0.0}}, "decel"),
        ({"type": "braking", "vehicle": {**VEHICLE, "brake_fraction": 1.5}}, "brake_fraction"),
        ({"type": "braking", "vehicle": VEHICLE, "pressure_profile": "parabolic"}, "pressure"),
        ({"type": "braking", "vehicle": {k: v for k, v in VEHICLE.items() if k != "mass"}},
         "Faltan"),
        ({"type": "braking", "vehicle": {**VEHICLE, "masa": 1.0}}, "desconocidas"),
        ({"type": "table", "points": [[0.5, 1.0], [1.0, 2.0]]}, "t = 0"),
        ({"type": "table", "points": [[0, 1.0], [0, 2.0]]}, "crecientes"),
        ({"type": "constant_flux"}, "requiere 'q'"),
        ({"type": "laser"}, "load.type"),
    ],
)
def test_load_config_errors(cfg, match):
    with pytest.raises(ValueError, match=match):
        make_load(cfg, GEO.r_a, GEO.r_b, GAMMA_REF)


# ------------------------------------------------------------------------- materiales

def test_library_reference_and_gamma():
    lib = load_library(LIBRARY)
    disc, pad = lib["gray_cast_iron"], lib["organic_pad"]
    assert disc.effusivity == pytest.approx(math.sqrt(50 * 7200 * 460), rel=RTOL)
    assert partition_coefficient(disc, pad) == pytest.approx(0.869, abs=5e-4)


def test_gamma_auto_from_region_adjacent_to_brake():
    cfg = yaml.safe_load(REF_CFG.read_text(encoding="utf-8"))
    lib = load_library(LIBRARY)
    mats = resolve_region_materials(cfg["materials"]["regions"], lib,
                                    cfg["materials"]["overrides"])
    mesh = build_mesh(GEO, MeshOptions(backend="structured", size_track=2e-3, size_far=4e-3))
    assert brake_regions(mesh) == ["disc"]
    g = resolve_gamma("auto", mesh, mats, lib[cfg["materials"]["pad"]])
    assert g == pytest.approx(0.8692006698254295, rel=RTOL)
    assert resolve_gamma(0.8, mesh, mats, None) == 0.8


def test_gamma_auto_rejects_two_materials_under_brake():
    lib = load_library(LIBRARY)
    mesh = build_mesh(GEO, MeshOptions(backend="structured", size_track=2e-3, size_far=4e-3))
    mesh.region = (mesh.r_centroid() > 0.11).astype(np.int64)
    mesh.region_names = ["a", "b"]
    mats = {"a": lib["gray_cast_iron"], "b": lib["organic_pad"]}
    with pytest.raises(ValueError, match="un solo material"):
        resolve_gamma("auto", mesh, mats, lib["organic_pad"])
    # Mismo material en ambas regiones: γ sigue siendo único.
    mats = {"a": lib["gray_cast_iron"], "b": lib["gray_cast_iron"]}
    assert resolve_gamma("auto", mesh, mats, lib["organic_pad"]) == pytest.approx(0.8692, abs=1e-4)


@pytest.mark.parametrize("bad", [1.0, 0.0, -0.2, "mucho"])
def test_gamma_numeric_validation(bad):
    with pytest.raises(ValueError, match="gamma"):
        resolve_gamma(bad)


def _write_lib(tmp_path, entry: dict) -> Path:
    p = tmp_path / "lib.yaml"
    p.write_text(yaml.safe_dump({"x": entry}), encoding="utf-8")
    return p


BASE = {"k": 50.0, "rho": 7200.0, "cp": 460.0, "source": "Libro, tabla 1", "verified": True}


@pytest.mark.parametrize("missing", ["k", "rho", "cp"])
def test_missing_required_key_fails_in_spanish(tmp_path, missing):
    entry = {k: v for k, v in BASE.items() if k != missing}
    with pytest.raises(ValueError, match=f"le falta la clave obligatoria '{missing}'"):
        load_library(_write_lib(tmp_path, entry))


@pytest.mark.parametrize(
    "entry",
    [
        {"k": 10.0, "rho": 1000.0, "cp": 1000.0},                      # sin source ni verified
        {"k": 10.0, "rho": 1000.0, "cp": 1000.0, "source": "  "},      # source vacía
        {"k": 10.0, "rho": 1000.0, "cp": 1000.0, "verified": False},   # sin source
    ],
)
def test_source_and_verified_are_optional(tmp_path, entry):
    """Sin source/verified el material carga, queda no verificado y la corrida lo advierte."""
    m = load_library(_write_lib(tmp_path, entry))["x"]
    assert m.source == SOURCE_MISSING and m.verified is False
    with pytest.warns(UnverifiedMaterialWarning, match=re.escape(f"x (source: {SOURCE_MISSING})")):
        assert warn_unverified([m]) == ["x"]


def test_source_without_verified_defaults_to_unverified(tmp_path):
    m = load_library(_write_lib(tmp_path, {k: v for k, v in BASE.items() if k != "verified"}))["x"]
    assert m.source == "Libro, tabla 1" and m.verified is False


@pytest.mark.parametrize(
    ("entry", "match"),
    [
        ({**BASE, "source": "  "}, "verified: true pero no declara 'source'"),
        ({k: v for k, v in BASE.items() if k != "source"}, "verified: true pero no declara"),
        ({**BASE, "source": 3}, "debe ser un texto"),
        ({**BASE, "verified": "si"}, "true o false"),
        ({**BASE, "k": -1.0}, "> 0"),
        ({**BASE, "conductividad": 1.0}, "desconocidas"),
    ],
)
def test_invalid_entries(tmp_path, entry, match):
    with pytest.raises(ValueError, match=match):
        load_library(_write_lib(tmp_path, entry))


def test_unverified_warning_lists_materials():
    lib = load_library(LIBRARY)
    with pytest.warns(UnverifiedMaterialWarning, match="gray_cast_iron.*organic_pad"):
        names = warn_unverified([lib["gray_cast_iron"], lib["organic_pad"],
                                 lib["gray_cast_iron"]])
    assert names == ["gray_cast_iron", "organic_pad"]
    ok = Material("ok", 1.0, 1000.0, 1000.0, "fuente", True)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert warn_unverified({"disc": ok}) == []


def test_range_warnings(tmp_path):
    with pytest.warns(PropertyRangeWarning, match="k = 1000"):
        load_library(_write_lib(tmp_path, {**BASE, "k": 1000.0}))
    with pytest.warns(PropertyRangeWarning, match="h = 2000"):
        check_range("h", 2000.0, "convección")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        load_library(_write_lib(tmp_path, BASE))
        check_range("h", 40.0)


def test_override_marks_unverified():
    base = Material("gci", 50.0, 7200.0, 460.0, "Libro", True)
    ov = apply_override(base, {"k": 48.0})
    assert ov.k == 48.0 and ov.rho == 7200.0 and not ov.verified
    assert "k=48" in ov.source and "Libro" in ov.source
    ov2 = apply_override(base, {"k": 48.0, "source": "Otro libro", "verified": True})
    assert ov2.verified and ov2.source == "Otro libro"
    with pytest.raises(ValueError, match="no permitidas"):
        apply_override(base, {"emisividad": 0.5})


def test_region_materials_errors():
    lib = load_library(LIBRARY)
    with pytest.raises(ValueError, match="no está en la librería"):
        resolve_region_materials({"disc": "acero"}, lib)
    with pytest.raises(ValueError, match="no existe"):
        resolve_region_materials({"disc": "gray_cast_iron"}, lib, {"pad": {"k": 1.0}})
