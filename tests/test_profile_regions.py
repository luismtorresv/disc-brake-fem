"""Pruebas de la Fase 7: perfil poligonal, regiones con materiales distintos, presets,
materiales definidos por el usuario y convección por grupo.

Criterios del PLAN (Fase 7): el caso rectangular sigue dando EXACTAMENTE los mismos
resultados (referencia guardada antes de los cambios, `tests/data/golden_rectangle.json`), y
dos regiones con material idéntico reproducen el caso homogéneo a 1e-12.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pytest
import yaml

from discfem.fem import assemble_stiffness_mass
from discfem.geometry import (
    BOUNDARY_GROUPS,
    RectangleGeometry,
    geometry_from_dict,
    point_segment_distance,
    polygon_area,
)
from discfem.loads import brake_load_vector, make_load
from discfem.materials import brake_regions, load_library, load_presets
from discfem.mesh import MeshOptions, boundary_edges, build_mesh
from discfem.solver import run
from discfem.timestepping import integrate, time_grid

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "configs" / "ejemplo_frenado.yaml"
PROFILE_CFG = ROOT / "configs" / "ejemplo_perfil.yaml"
GOLDEN = json.loads((ROOT / "tests" / "data" / "golden_rectangle.json").read_text("utf-8"))

HUB = {
    "type": "profile",
    "points": [[0.07, 0], [0.15, 0], [0.15, 0.01], [0.14, 0.01], [0.08, 0.01], [0.075, 0.01],
               [0.075, 0.016], [0.07, 0.016]],
    "segments": ["symmetry", "r_outer", "top_outer", "brake", "top_inner", "top_inner",
                 "top_inner", "r_inner"],
}
RECT = {"type": "rectangle", "r_inner": 0.07, "r_outer": 0.15, "half_thickness": 0.01,
        "track": {"r_a": 0.08, "r_b": 0.14}}
LAYER = {**RECT, "regions": {"capa": [[0.08, 0.009], [0.14, 0.009], [0.14, 0.01], [0.08, 0.01]]}}
CUSTOM_SOFT = {"k": 25.0, "rho": 7200.0, "cp": 460.0,
               "source": "valor de prueba (no es un material real)", "verified": False}


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


def _sha(a, dtype=None) -> str:
    a = np.ascontiguousarray(a if dtype is None else np.asarray(a, dtype=dtype))
    return hashlib.sha256(a.tobytes()).hexdigest()


# --------------------------------------------------------------- rectángulo sin cambios

@pytest.mark.parametrize("backend", ["gmsh", "structured"])
def test_rectangle_unchanged(backend):
    """El rectángulo da EXACTAMENTE lo mismo que antes de la Fase 7 (commit c930b37):
    malla, aristas por grupo, campo final y todas las entradas del summary de entonces."""
    gold = GOLDEN["cases"][backend]
    res = run(REF, {"mesh.backend": backend}, write=False)
    m = res.mesh
    assert _sha(m.nodes) == gold["nodes_sha256"]
    assert _sha(m.tris, np.int64) == gold["tris_sha256"]
    for g, h in gold["edges_sha256"].items():
        assert _sha(m.edges[g], np.int64) == h, g
    assert _sha(res.transient.T_final) == gold["T_final_sha256"]
    for key, value in gold["summary"].items():
        assert res.summary[key] == value, key


def test_profile_with_rectangle_corners_equals_rectangle():
    rect = RectangleGeometry.from_dict(RECT)
    prof = geometry_from_dict({"type": "profile", "points": [list(p) for p in rect.corners],
                               "segments": list(BOUNDARY_GROUPS)})
    m1 = build_mesh(rect, MeshOptions())
    m2 = build_mesh(prof, MeshOptions())
    assert np.array_equal(m1.nodes, m2.nodes) and np.array_equal(m1.tris, m2.tris)
    assert all(np.array_equal(m1.edges[g], m2.edges[g]) for g in BOUNDARY_GROUPS)
    assert (prof.r_a, prof.r_b) == (rect.r_a, rect.r_b)


# --------------------------------------------------------------- perfil poligonal

@pytest.fixture(scope="module")
def hub_mesh():
    geo = geometry_from_dict(HUB)
    return geo, build_mesh(geo, MeshOptions())


def test_profile_weighted_area_exact(hub_mesh):
    geo, m = hub_mesh
    A = m.areas()
    assert np.all(A > 0)
    assert abs(A.sum() - polygon_area(HUB["points"])) / A.sum() < 1e-12
    assert abs((A * m.r_centroid()).sum() - geo.r_moment()) / geo.r_moment() < 1e-12


def test_profile_groups_lie_on_their_segments(hub_mesh):
    geo, m = hub_mesh
    keys = [set(map(tuple, np.sort(m.edges[g], axis=1))) for g in BOUNDARY_GROUPS]
    assert sum(len(k) for k in keys) == len(set().union(*keys))  # disjuntos
    assert set().union(*keys) == set(map(tuple, np.sort(boundary_edges(m.tris), axis=1)))
    for g in BOUNDARY_GROUPS:
        segs = [(a, b) for a, b, gg in geo.boundary_segments() if gg == g]
        pts = m.nodes[np.unique(m.edges[g])]
        d = np.min([point_segment_distance(pts, a, b) for a, b in segs], axis=0)
        assert np.all(d < 1e-12), g
    # Longitud total de cada grupo = suma de sus lados.
    for g in BOUNDARY_GROUPS:
        seg = m.nodes[m.edges[g]]
        L = np.linalg.norm(seg[:, 1] - seg[:, 0], axis=1).sum()
        L_exact = sum(math.dist(a, b) for a, b, gg in geo.boundary_segments() if gg == g)
        assert abs(L - L_exact) <= 1e-12 * L_exact, g


def test_profile_brake_spans_track_and_refined(hub_mesh):
    geo, m = hub_mesh
    b = m.nodes[np.unique(m.edges["brake"])]
    assert b[:, 0].min() == geo.r_a == 0.08 and b[:, 0].max() == geo.r_b == 0.14
    assert np.all(b[:, 1] == 0.01)
    A = m.areas()
    touch = np.isin(m.tris, np.unique(m.edges["brake"])).any(axis=1)
    deep = m.nodes[m.tris].mean(axis=1)[:, 1] < 0.002
    assert A[touch].mean() < 0.1 * A[deep].mean()


def test_profile_power_consistency(hub_mesh):
    """2π Σ F_brake = γ P_f también en el perfil (pista horizontal con nodos en r_a, r_b)."""
    geo, m = hub_mesh
    cfg = yaml.safe_load(REF.read_text("utf-8"))["load"]
    ld = make_load(cfg, geo.r_a, geo.r_b, 0.869)
    F = brake_load_vector(m, ld)
    for t in (0.0, 3.0):
        assert 2 * math.pi * F(t).sum() == pytest.approx(0.869 * ld.power_per_face(t), rel=1e-12)


def test_profile_example_runs():
    res = run(PROFILE_CFG, write=False)
    s = res.summary
    assert s["geometry"]["type"] == "profile"
    assert s["energy"]["balance_relative_residual"] < 1e-10
    assert s["materials"]["disc"]["name"] == "gray_cast_iron_limpert"  # vía preset 'car'
    assert s["convection"]["r_outer"]["h"] == 60.0 and s["convection"]["top_inner"]["h"] == 40.0


def test_structured_rejects_profile_and_regions():
    with pytest.raises(ValueError, match="structured"):
        build_mesh(geometry_from_dict(HUB), MeshOptions(backend="structured"))
    with pytest.raises(ValueError, match="structured"):
        build_mesh(geometry_from_dict(LAYER), MeshOptions(backend="structured"))


@pytest.mark.parametrize(
    ("change", "match"),
    [
        ({"points": HUB["points"][::-1], "segments": HUB["segments"][::-1]}, "ANTIHORARIO"),
        ({"segments": HUB["segments"][:-1]}, "un grupo por lado"),
        ({"segments": ["lateral", *HUB["segments"][1:]]}, "desconocidos"),
        ({"points": [[0.07, 0], [0.15, 0], [0.15, 0.01], [0.14, 0.012], [0.08, 0.01],
                     [0.075, 0.01], [0.075, 0.016], [0.07, 0.016]]}, "horizontales"),
        ({"segments": ["symmetry", "r_outer", "brake", "top_outer", "brake", "top_inner",
                       "top_inner", "r_inner"]}, "cadena contigua"),
        ({"segments": ["symmetry", "r_outer", "top_outer", "top_outer", "top_inner",
                       "top_inner", "top_inner", "r_inner"]}, "brake"),
        ({"points": [[-0.01, 0], *HUB["points"][1:]]}, "r < 0"),
        ({"points": [[0.07, 0], [0.15, 0], [0.15, 0.01], [0.06, 0.005], [0.08, 0.01],
                     [0.075, 0.01], [0.075, 0.016], [0.07, 0.016]]}, "autointersecta"),
        ({"regions": {"disc": [[0.08, 0.009], [0.09, 0.009], [0.09, 0.01]]}}, "por defecto"),
        ({"regions": {"x": [[0.08, 0.009], [0.2, 0.009], [0.2, 0.01], [0.08, 0.01]]}},
         "fuera de la sección"),
    ],
)
def test_profile_validation_errors(change, match):
    with pytest.raises(ValueError, match=match):
        geometry_from_dict({**HUB, **change})


# --------------------------------------------------------------- regiones

@pytest.fixture(scope="module")
def layer_mesh():
    geo = geometry_from_dict(LAYER)
    return geo, build_mesh(geo, MeshOptions())


def test_regions_areas_and_conformity(layer_mesh):
    _, m = layer_mesh
    assert m.region_names == ["disc", "capa"]
    A = m.areas()
    capa = A[m.region == 1].sum()
    assert abs(capa - polygon_area(LAYER["regions"]["capa"])) < 1e-12 * capa
    assert abs(A.sum() - 0.01 * 0.08) < 1e-12 * A.sum()
    assert brake_regions(m) == ["capa"]
    # Conformidad: la interfaz z = 0.009, r ∈ [0.08, 0.14] está cubierta exactamente por
    # aristas de la malla (su longitud total es 0.06 m), así que ningún triángulo la cruza.
    e = np.concatenate([m.tris[:, [0, 1]], m.tris[:, [1, 2]], m.tris[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    p, q = m.nodes[e[:, 0]], m.nodes[e[:, 1]]
    on = (np.abs(p[:, 1] - 0.009) < 1e-15) & (np.abs(q[:, 1] - 0.009) < 1e-15)
    on &= (np.minimum(p[:, 0], q[:, 0]) >= 0.08) & (np.maximum(p[:, 0], q[:, 0]) <= 0.14)
    assert abs(np.linalg.norm(q[on] - p[on], axis=1).sum() - 0.06) < 1e-12


def test_regions_overlap_detected():
    over = {**RECT, "regions": {"a": [[0.08, 0.008], [0.12, 0.008], [0.12, 0.01], [0.08, 0.01]],
                                "b": [[0.10, 0.008], [0.14, 0.008], [0.14, 0.01], [0.10, 0.01]]}}
    with pytest.raises(ValueError, match="solapa|área"):
        build_mesh(geometry_from_dict(over), MeshOptions())


def test_two_regions_same_material_reproduce_homogeneous(layer_mesh):
    """Criterio del PLAN: misma malla, material idéntico en ambas regiones ⇒ el mismo campo
    que el caso homogéneo (una sola región) a 1e-12."""
    _, m = layer_mesh
    k, rc = 50.0, 7200.0 * 460.0
    K1, M1 = assemble_stiffness_mass(m, {"disc": k, "capa": k}, {"disc": rc, "capa": rc})
    K0, M0 = assemble_stiffness_mass(m, k, rc)
    cfg = yaml.safe_load(REF.read_text("utf-8"))["load"]
    geo = geometry_from_dict(LAYER)
    ld = make_load(cfg, geo.r_a, geo.r_b, 0.869)
    F = brake_load_vector(m, ld)
    t = time_grid(20.0, t_f=ld.t_f)
    r1 = integrate(M1, K1, None, 25.0, t, F_brake=F, t_f=ld.t_f, theta=1.0, rannacher_steps=0)
    r0 = integrate(M0, K0, None, 25.0, t, F_brake=F, t_f=ld.t_f, theta=1.0, rannacher_steps=0)
    rise = r0.T_max.max() - 25.0
    assert np.abs(r1.saved_T - r0.saved_T).max() < 1e-12 * rise


def test_layer_same_material_close_to_rectangle():
    """Contra la malla SIN capa la diferencia es solo de discretización (mallas distintas)."""
    a = run(REF, {"geometry": LAYER, "materials.regions.capa": "gray_cast_iron"}, write=False)
    b = run(REF, write=False)
    assert a.summary["T_max"]["value"] == pytest.approx(b.summary["T_max"]["value"], rel=5e-3)


def test_layer_less_conductive_heats_more_and_sets_gamma():
    """Capa de k menor bajo la pista ⇒ T_max mayor; γ automático usa el material de la capa."""
    base = {"geometry": LAYER, "materials.regions.capa": "gray_cast_iron"}
    soft = {"geometry": LAYER, "materials.custom": {"blando": CUSTOM_SOFT},
            "materials.regions.capa": "blando"}
    a = run(REF, base, write=False)
    b = run(REF, soft, write=False)
    assert b.summary["T_max"]["value"] > a.summary["T_max"]["value"] + 1.0
    assert b.gamma < a.gamma  # menor efusividad del lado del disco
    assert b.materials["capa"].name == "blando"
    assert b.summary["energy"]["balance_relative_residual"] < 1e-10


# --------------------------------------------------------------- materiales y config

def test_library_new_format_and_presets(tmp_path):
    lib = ROOT / "configs" / "materiales.yaml"
    mats = load_library(lib)
    presets = load_presets(lib, mats)
    assert presets["car"]["disc"] == "gray_cast_iron_limpert"
    # El formato plano anterior sigue funcionando (p. ej. carpetas de resultados viejas).
    flat = tmp_path / "flat.yaml"
    flat.write_text(yaml.safe_dump({"m": {"k": 1.0, "rho": 1000.0, "cp": 1000.0,
                                          "source": "x", "verified": True}}), "utf-8")
    assert list(load_library(flat)) == ["m"]
    assert load_presets(flat, load_library(flat)) == {}


def test_preset_and_explicit_priority():
    raw = yaml.safe_load(REF.read_text("utf-8"))
    raw["materials"] = {"library": "configs/materiales.yaml", "preset": "car"}
    res = run(raw, {"time.t_end": 6.0, "mesh.backend": "structured"}, write=False)
    assert res.materials["disc"].name == "gray_cast_iron_limpert"
    assert res.pad.name == "pad_limpert"
    raw["materials"]["custom"] = {"blando": CUSTOM_SOFT}
    raw["materials"]["regions"] = {"disc": "blando"}
    res = run(raw, {"time.t_end": 6.0, "mesh.backend": "structured"}, write=False)
    assert res.materials["disc"].name == "blando" and res.pad.name == "pad_limpert"


@pytest.mark.parametrize(
    ("change", "match"),
    [
        ({"materials.preset": "camion"}, "no existe"),
        ({"materials.custom": {"gray_cast_iron": CUSTOM_SOFT}}, "ya existe"),
        ({"materials.custom": {"x": {k: v for k, v in CUSTOM_SOFT.items() if k != "source"}},
          "materials.regions.disc": "x"}, "source"),
    ],
)
def test_material_errors(change, match):
    with pytest.raises(ValueError, match=match):
        run(REF, {**change, "mesh.backend": "structured", "time.t_end": 6.0}, write=False)


@pytest.mark.parametrize(
    ("change", "match"),
    [
        ({"geometry": LAYER}, r"capa"),  # falta el material de la región 'capa'
        ({"materials.regions.otra": "gray_cast_iron"}, "no tiene"),
        ({"convection.per_group": {"symmetry": {"h": 10.0}}}, "sin condición"),
        ({"convection.per_group": {"r_outer": {"hh": 10.0}}}, "desconocidas"),
    ],
)
def test_config_errors_fase7(change, match):
    from discfem.config import load_config

    with pytest.raises(ValueError, match=match):
        load_config(REF, change)


def test_per_group_convection():
    """per_group igual al global ⇒ resultado idéntico; h distinto en r_outer ⇒ cambia la
    energía que sale; T∞ distinto por grupo con T0 = T∞ común ya no es estacionario."""
    base = run(REF, {"mesh.backend": "structured"}, write=False).summary
    same = run(REF, {"mesh.backend": "structured",
                     "convection.per_group": {"r_outer": {"h": 40.0, "T_inf": 25.0}}},
               write=False).summary
    assert same["T_max"] == base["T_max"] and same["energy"] == base["energy"]
    more = run(REF, {"mesh.backend": "structured",
                     "convection.per_group": {"r_outer": {"h": 400.0}}}, write=False).summary
    assert more["energy"]["convection_out_J"] > base["energy"]["convection_out_J"]
    assert more["T_final"]["max"] < base["T_final"]["max"]
    # Sin flujo y T0 = T∞ en todos los grupos: uniforme (T5) también con per_group.
    calm = run(REF, {"mesh.backend": "structured", "boundaries.brake": "insulated",
                     "convection.per_group": {"r_outer": {"h": 400.0}}, "time.dt": 0.5,
                     "time.t_end": 10.0}, write=False)
    assert np.abs(calm.transient.saved_T - 25.0).max() < 1e-10
    hot_air = run(REF, {"mesh.backend": "structured", "boundaries.brake": "insulated",
                        "convection.per_group": {"r_outer": {"T_inf": 80.0}}, "time.dt": 0.5,
                        "time.t_end": 10.0}, write=False)
    assert hot_air.summary["T_final"]["max"] > 25.0 + 1e-3


def test_copy_of_raw_not_mutated_by_overrides():
    raw = yaml.safe_load(REF.read_text("utf-8"))
    before = copy.deepcopy(raw)
    run(raw, {"geometry": LAYER, "materials.regions.capa": "gray_cast_iron",
              "time.t_end": 6.0}, write=False)
    assert raw == before


# --------------------------------------------------------------- catálogo con fuente

def test_limpert_materials_match_table_3_1A():
    """Limpert, Brake Design and Safety, 2.ª ed., Tabla 3-1A (SI): k en N·m/(h·K·m) / 3600.
    Se comprueba también contra la tabla imperial 3-1B y la difusividad tabulada."""
    lib = load_library(ROOT / "configs" / "materiales.yaml")
    d, p = lib["gray_cast_iron_limpert"], lib["pad_limpert"]
    assert (d.k, d.rho, d.cp) == (174465 / 3600, 7228.0, 419.0)
    assert (p.k, p.rho, p.cp) == (4362 / 3600, 2595.0, 1465.0)
    assert d.verified and p.verified and "Tabla 3-1A" in d.source
    btu_hftF = 1.730735  # W/(m·K) por BTU/(h·ft·°F)
    assert d.k == pytest.approx(28 * btu_hftF, rel=1e-3)      # Tabla 3-1B
    assert p.k == pytest.approx(0.7 * btu_hftF, rel=1e-3)
    assert d.cp == pytest.approx(0.10 * 4186.8, rel=1e-3)
    assert d.diffusivity * 3600 == pytest.approx(0.0576, rel=1e-3)  # a tabulada, m²/h
    assert p.diffusivity * 3600 == pytest.approx(0.0011, abs=5e-5)  # 2 cifras en la tabla


def test_auto_example_uses_only_verified_materials():
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # ninguna advertencia de materiales no verificados
        res = run(ROOT / "configs" / "ejemplo_auto.yaml", {"mesh.refine_factor": 1.0},
                  write=False)
    s = res.summary
    assert s["unverified_materials"] == []
    assert s["load"]["gamma"] == pytest.approx(0.8495, abs=1e-4)
    assert s["energy"]["balance_relative_residual"] < 1e-10


def test_stainless_410_entry_from_datasheet():
    """Ficha Chilexpo 'Aleación 410' (UNS S41000): ρ = 0,266 lb/in³, c_p = 0,46 kJ/(kg·K)
    (la ficha rotula mal la unidad), k = 26,9 W/(m·K) a 100 °C. No verificado (§ nota)."""
    m = load_library(ROOT / "configs" / "materiales.yaml")["stainless_steel_410"]
    assert m.rho == pytest.approx(0.266 * 27679.9047102, rel=1e-12)
    assert (m.k, m.cp) == (26.9, 460.0)
    assert m.k == pytest.approx(187 * 0.144227889, rel=5e-3)  # 187 BTU·in/(h·ft²·°F)
    assert not m.verified and "S41000" in m.source


def test_sintered_pad_entry_goo2014():
    """Goo y Na (2014), Tabla 1 'Friction pin': ρ = 5,4 g/cm³, c = 418, k = 24 (no verificado)."""
    m = load_library(ROOT / "configs" / "materiales.yaml")["sintered_pad_goo2014"]
    assert (m.k, m.rho, m.cp) == (24.0, 5400.0, 418.0)
    assert not m.verified and "10.9725/kstle.2014.30.3.177" in m.source


def test_sintered_pad_cu_zhou2024_entry():
    """k y ρ digitalizados de Zhou et al. (2024), Fig. 6, S800; c de Goo y Na (2014)."""
    lib = load_library(ROOT / "configs" / "materiales.yaml")
    m = lib["sintered_pad_cu_zhou2024"]
    assert (m.k, m.rho) == (45.1, 5520.0)
    assert m.cp == lib["sintered_pad_goo2014"].cp == 418.0
    assert not m.verified
    assert "10.3390/met14091048" in m.source and "10.9725/kstle.2014.30.3.177" in m.source


def test_motorcycle_preset():
    """Preset de moto: disco AISI 410 + pastilla sinterizada de cobre (Zhou 2024 / Goo 2014);
    ambos se advierten como no verificados y γ = e_410 / (e_410 + e_pastilla)."""
    raw = yaml.safe_load(REF.read_text("utf-8"))
    raw["materials"] = {"library": "configs/materiales.yaml", "preset": "motorcycle"}
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        res = run(raw, {"mesh.backend": "structured", "time.t_end": 6.0}, write=False)
    assert res.materials["disc"].name == "stainless_steel_410"
    assert res.pad.name == "sintered_pad_cu_zhou2024"
    assert res.summary["unverified_materials"] == ["stainless_steel_410",
                                                   "sintered_pad_cu_zhou2024"]
    assert any("sintered_pad_cu_zhou2024" in str(w.message) for w in rec)
    d, p = res.materials["disc"], res.pad
    assert res.gamma == pytest.approx(d.effusivity / (d.effusivity + p.effusivity), rel=1e-12)
    assert res.gamma == pytest.approx(0.4834, abs=1e-4)


def test_preset_without_pad_needs_explicit_pad(tmp_path):
    """Un preset sin pastilla (permitido en la librería) exige materials.pad con gamma: auto."""
    lib = tmp_path / "lib.yaml"
    entry = {"k": 26.9, "rho": 7800.0, "cp": 460.0, "source": "prueba", "verified": True}
    lib.write_text(yaml.safe_dump({"materials": {"acero": entry, "pastilla": {**entry, "k": 1.0}},
                                   "presets": {"solo_disco": {"disc": "acero"}}}), "utf-8")
    raw = yaml.safe_load(REF.read_text("utf-8"))
    raw["materials"] = {"library": str(lib), "preset": "solo_disco"}
    ov = {"mesh.backend": "structured", "time.t_end": 6.0}
    with pytest.raises(ValueError, match="preset 'solo_disco' no la define"):
        run(raw, ov, write=False)
    raw["materials"]["pad"] = "pastilla"
    assert run(raw, ov, write=False).pad.name == "pastilla"
