"""Pruebas de la Fase 5: T8 y T9 (PLAN.md §6), configuración, salidas, reproducibilidad y CLI.

T8: caso de referencia (`configs/ejemplo_frenado.yaml`, desgaste uniforme, Δt = t_f/100,
    θ = 1 según el YAML): T_max = 102,5 ± 0,5 °C en t = 3,4 ± 0,1 s; T_max(60 s) = 72,3 ± 0,5 °C.
T9: ΔT_final ≤ γ E_cara / (ρ c π (r_o² − r_i²) e)  (cota 48,1 K; prototipo 47,3 K).
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

from discfem.config import load_config, set_dotted
from discfem.geometry import RectangleGeometry
from discfem.materials import UnverifiedMaterialWarning
from discfem.mesh import MeshOptions, build_mesh
from discfem.postprocess import interpolation_matrix
from discfem.solver import run

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "configs" / "ejemplo_frenado.yaml"


@pytest.fixture(scope="module")
def ref_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("ref")
    with pytest.warns(UnverifiedMaterialWarning):
        res = run(REF, {"output.dir": str(out)})
    return res


# --------------------------------------------------------------------------- T8 / T9

def test_t8_reference_case(ref_run):
    s = ref_run.summary
    assert s["mesh"]["n_nodes"] >= 700  # resolución ≥ 40×16 estructurada (§6)
    assert s["time"]["theta"] == 1.0 and s["load"]["pressure_profile"] == "uniform_wear"
    assert s["time"]["dt"] == pytest.approx(s["time"]["t_f"] / 100, rel=1e-9)
    assert abs(s["T_max"]["value"] - 102.5) <= 0.5
    assert abs(s["T_max"]["time"] - 3.4) <= 0.1
    assert abs(s["T_final"]["max"] - 72.3) <= 0.5
    assert s["T_final"]["time"] == 60.0


def test_t9_energy_bound(ref_run):
    s, cfg = ref_run.summary, ref_run.config
    disc = ref_run.materials["disc"]
    g = cfg.geometry
    bound = (s["load"]["gamma"] * s["load"]["energy_per_face_J"]
             / (disc.rho_cp * math.pi * (g.r_outer**2 - g.r_inner**2) * g.half_thickness))
    assert bound == pytest.approx(48.1, abs=0.05)
    dT_final = s["T_final"]["max"] - cfg.T0
    assert dT_final <= bound
    # La cota es sobre el ΔT MEDIO (energía/capacidad, sin pérdidas): también se cumple.
    M_cap = disc.rho_cp * math.pi * (g.r_outer**2 - g.r_inner**2) * g.half_thickness
    assert s["energy"]["stored_J"] / M_cap <= bound


def test_energy_consistency_reference(ref_run):
    """Balance cerrado y energía de entrada = γ E_cara·(N−1)/N: la regla θ = 1 evalúa la carga
    lineal al final de cada paso (ver reporte F3)."""
    s = ref_run.summary
    e = s["energy"]
    assert e["balance_relative_residual"] < 1e-10
    assert e["stored_J"] == pytest.approx(e["brake_in_J"] - e["convection_out_J"], rel=1e-10)
    N = 100
    assert e["brake_in_J"] == pytest.approx(s["load"]["gamma_energy_per_face_J"] * (N - 1) / N,
                                            rel=1e-10)


# --------------------------------------------------------------------------- salidas

EXPECTED_FILES = ["config.yaml", "materiales.yaml", "summary.json", "tmax.csv", "probes.csv",
                  "mesh.png", "field.png", "tmax.png", "depth_profiles.png", "field.pvd",
                  "field_0000.vtu"]


def test_outputs_written(ref_run):
    out = ref_run.output_dir
    for name in EXPECTED_FILES:
        assert (out / name).stat().st_size > 0, name
    tmax = np.loadtxt(out / "tmax.csv", delimiter=",", skiprows=1)
    assert tmax.shape == (ref_run.transient.times.size, 5)
    assert np.max(tmax[:, 1]) == ref_run.summary["T_max"]["value"]
    probes = np.loadtxt(out / "probes.csv", delimiter=",", skiprows=1)
    assert probes.shape[1] == 1 + len(ref_run.config.output.probes)
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert summary["unverified_materials"] == ["gray_cast_iron", "organic_pad"]
    assert {"T_max", "energy", "mesh", "compute_time_s", "probes"} <= set(summary)


def _strip_time(s: dict) -> dict:
    return {k: v for k, v in s.items() if k != "compute_time_s"}


def test_rerun_gives_same_summary(ref_run, tmp_path):
    """Correr dos veces da el mismo summary.json (salvo el tiempo de cómputo)."""
    with pytest.warns(UnverifiedMaterialWarning):
        res2 = run(REF, {"output.dir": str(tmp_path)})
    s1 = json.loads((ref_run.output_dir / "summary.json").read_text(encoding="utf-8"))
    s2 = json.loads((res2.output_dir / "summary.json").read_text(encoding="utf-8"))
    assert _strip_time(s1) == _strip_time(s2)


def test_copied_config_reproduces_run(ref_run, tmp_path):
    """Regla 5: la corrida se reproduce solo con el config.yaml copiado (+ su materiales.yaml)."""
    copied = ref_run.output_dir / "config.yaml"
    with pytest.warns(UnverifiedMaterialWarning):
        res2 = run(copied, {"output.dir": str(tmp_path)})
    assert _strip_time(res2.summary) == _strip_time(ref_run.summary)


def test_overrides_with_dotted_keys(tmp_path):
    with pytest.warns(UnverifiedMaterialWarning):
        res = run(REF, {"output.dir": str(tmp_path), "output.formats": [],
                        "mesh.backend": "structured", "mesh.refine_factor": 0.5,
                        "materials.overrides.disc.k": 40.0, "time.t_end": 10.0})
    assert res.summary["mesh"]["backend"] == "structured"
    assert res.materials["disc"].k == 40.0 and not res.materials["disc"].verified
    assert res.summary["T_final"]["time"] == 10.0
    saved = yaml.safe_load((tmp_path / "config.yaml").read_text(encoding="utf-8"))
    assert saved["mesh"]["refine_factor"] == 0.5
    assert saved["materials"]["overrides"]["disc"]["k"] == 40.0
    assert not (tmp_path / "tmax.csv").exists()  # formats: []


def test_run_without_writing():
    with pytest.warns(UnverifiedMaterialWarning):
        res = run(REF, {"mesh.backend": "structured", "time.t_end": 6.0}, write=False)
    assert res.output_dir is None


# --------------------------------------------------------------------------- sondas

def test_interpolation_exact_for_linear_field():
    geo = RectangleGeometry(0.07, 0.15, 0.01, 0.08, 0.14)
    mesh = build_mesh(geo, MeshOptions(backend="gmsh", size_track=1e-3, size_far=4e-3))
    f = 3.0 + 200.0 * mesh.nodes[:, 0] - 50.0 * mesh.nodes[:, 1]
    pts = np.array([[0.11, 0.01], [0.0712, 0.0], [0.15, 0.005], [0.123, 0.0047]])
    P = interpolation_matrix(mesh, pts)
    np.testing.assert_allclose(P @ f, 3.0 + 200.0 * pts[:, 0] - 50.0 * pts[:, 1], rtol=1e-12)
    with pytest.raises(ValueError, match="fuera del dominio"):
        interpolation_matrix(mesh, [[0.2, 0.0]])


# --------------------------------------------------------------------------- configuración

def ref_raw() -> dict:
    return yaml.safe_load(REF.read_text(encoding="utf-8"))


@pytest.mark.parametrize(
    ("key", "value", "match"),
    [
        ("boundaries.brake", "radiation", "debe ser uno de"),
        ("boundaries.lateral", "insulated", "desconocidas"),
        ("time.theta", 1.5, "theta"),
        ("time.rannacher_steps", -1, "rannacher"),
        ("time.t_end", 0.0, "t_end"),
        ("output.formats", ["pdf"], "formats"),
        ("output.probes", [[0.1]], "pares"),
        ("convection.h", -5.0, "> 0"),
        ("materials.regions", {"disc": "x", "capa": "y"}, "una sola región"),
        ("geometry.track.r_b", 0.2, "r_outer"),
        ("extra", 1, "raíz"),
    ],
)
def test_config_errors_in_spanish(key, value, match):
    raw = ref_raw()
    set_dotted(raw, key, value)
    with pytest.raises(ValueError, match=match):
        load_config(raw)


def test_config_missing_boundary():
    raw = ref_raw()
    del raw["boundaries"]["symmetry"]
    with pytest.raises(ValueError, match="symmetry"):
        load_config(raw)


def test_set_dotted_creates_and_rejects():
    d = {"a": {"b": 1}}
    set_dotted(d, "a.c.d", 2)
    assert d == {"a": {"b": 1, "c": {"d": 2}}}
    with pytest.raises(ValueError, match="no es una sección"):
        set_dotted(d, "a.b.x", 3)


def test_config_defaults():
    raw = ref_raw()
    del raw["time"]["theta"], raw["time"]["rannacher_steps"]
    cfg = load_config(raw)
    assert cfg.time.theta == 0.5 and cfg.time.rannacher_steps == 2
    assert cfg.groups_with("convection") == ["r_outer", "top_outer", "top_inner", "r_inner"]
    assert cfg.groups_with("flux") == ["brake"]


# --------------------------------------------------------------------------- CLI

def test_cli_run_subprocess(tmp_path):
    """`python -m discfem run` en una consola cp1252, con --set; advertencia visible."""
    import os

    env = {**os.environ, "PYTHONIOENCODING": "cp1252"}
    proc = subprocess.run(
        [sys.executable, "-m", "discfem", "run", str(REF),
         "--set", f"output.dir={tmp_path.as_posix()}", "--set", "mesh.backend=structured",
         "--set", "output.formats=[csv]"],
        check=False, env=env, capture_output=True, cwd=ROOT,
    )
    err = proc.stderr.decode("utf-8", "replace")
    assert proc.returncode == 0, err
    assert "UnverifiedMaterialWarning" in err and "gray_cast_iron" in err
    assert (tmp_path / "summary.json").is_file()


def test_cli_run_reports_config_error(tmp_path):
    from discfem.cli import main

    assert main(["run", str(REF), "--set", "time.theta=7"]) == 2
