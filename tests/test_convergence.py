"""Pruebas de la Fase 6: fórmulas de convergencia, estudio de malla del caso de referencia,
estudio de Δt, barridos (tendencias físicas) y el comando ``study``."""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pytest
import yaml

from discfem.convergence import (
    Level,
    convergence_table,
    load_study,
    mesh_study,
    observed_order,
    parameter_sweep,
    relative_change_percent,
    run_study,
    timestep_study,
)

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "configs" / "ejemplo_frenado.yaml"
# Malla estructurada gruesa y 15 s simulados: rápido, suficiente para tendencias.
FAST = {"mesh.backend": "structured", "mesh.size_track": 1e-3, "mesh.size_far": 4e-3,
        "time.t_end": 15.0}


@pytest.fixture(autouse=True)
def _quiet_material_warnings():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        yield


# ------------------------------------------------------------------------------ fórmulas

@pytest.mark.parametrize(("p", "ratio"), [(2.0, 2.0), (1.0, 2.0), (1.5, 4.0)])
def test_observed_order_exact_for_power_law(p, ratio):
    """Q(h) = Q* + C h^p con h, h/ρ, h/ρ² ⇒ el orden observado es p."""
    h = np.array([1.0, 1 / ratio, 1 / ratio**2])
    q = 10.0 + 3.0 * h**p
    assert observed_order(*q, ratio=ratio) == pytest.approx(p, rel=1e-10)


def test_observed_order_degenerate():
    assert observed_order(1.0, 1.0, 2.0) is None
    assert observed_order(1.0, 2.0, 2.0) is None


def test_relative_change_percent():
    assert relative_change_percent(99.0, 100.0) == pytest.approx(1.0)
    assert relative_change_percent(101.0, 100.0) == pytest.approx(1.0)


def _fake(T: float, probes: list[float]) -> Level:
    s = {"T_max": {"value": T, "time": 1.0}, "T_final": {"max": T}, "T0": 25.0,
         "mesh": {"n_nodes": 1, "n_elements": 1}, "time": {"dt": 0.1, "n_steps": 10},
         "probes": [{"T_max": v} for v in probes], "compute_time_s": 0.0}
    return Level("x", 0, s, np.zeros(1), np.zeros(1))


def test_selection_requires_probes_too():
    """T_max converge en el nivel 1, pero una sonda no: se selecciona el nivel 2."""
    levels = [_fake(100.0, [50.0]), _fake(100.5, [52.0]), _fake(100.6, [52.1])]
    rows, sel = convergence_table(levels, [1, 2, 4], 25.0, tol=1.0)
    assert rows[1]["eps_percent"] < 1.0 and rows[1]["eps_probes_max_percent"] > 1.0
    assert sel == 2
    assert rows[2]["order_p"] == pytest.approx(np.log2(0.5 / 0.1))
    rows, sel = convergence_table(levels, [1, 2, 4], 25.0, tol=0.01)
    assert sel is None


def test_order_skipped_when_ratio_not_constant():
    levels = [_fake(100.0, []), _fake(100.5, []), _fake(100.6, [])]
    rows, _ = convergence_table(levels, [1, 2, 3], 25.0, tol=1.0)
    assert rows[2]["order_p"] is None


# ------------------------------------------------------------------------------ estudios

def test_mesh_study_reference_closes_below_one_percent():
    """Criterio de la Fase 6: el estudio de malla del caso de referencia cierra < 1 % desde el
    segundo nivel (también en las sondas) y ε disminuye al refinar."""
    _, rows, sel = mesh_study(REF, [1.0, 2.0, 4.0], {}, tol=1.0)
    assert sel == 1
    eps = [r["eps_percent"] for r in rows[1:]]
    assert all(e < 1.0 for e in eps) and eps[1] < eps[0]
    assert all(r["eps_probes_max_percent"] < 1.0 for r in rows[1:])
    assert rows[0]["n_nodes"] < rows[1]["n_nodes"] < rows[2]["n_nodes"]
    assert rows[2]["order_p"] is not None and rows[2]["order_p"] > 0


def test_timestep_study_first_order_for_implicit_euler():
    """Con θ = 1 (YAML de referencia) el orden observado en Δt debe ser ≈ 1."""
    _, rows = timestep_study(REF, [1, 2, 4], FAST)
    assert [r["n_steps"] for r in rows][1] == 2 * rows[0]["n_steps"]
    assert rows[2]["order_p"] == pytest.approx(1.0, abs=0.05)


@pytest.mark.parametrize("key", ["load.vehicle.decel", "load.vehicle.mass", "load.vehicle.v0"])
def test_sweeps_more_load_more_temperature(key):
    """Más desaceleración, más masa o más velocidad ⇒ más T_max (tendencia esperada)."""
    base = {"load.vehicle.decel": 5.0, "load.vehicle.mass": 1500.0, "load.vehicle.v0": 27.78}
    values = [0.8 * base[key], base[key], 1.2 * base[key]]
    _, rows = parameter_sweep(REF, key, values, FAST)
    T = [r["T_max"] for r in rows]
    assert T[0] < T[1] < T[2]


def test_sweep_decel_keeps_energy():
    """La energía que entra no depende de a (solo de m v0²); la frenada es más corta."""
    _, rows = parameter_sweep(REF, "load.vehicle.decel", [4.0, 6.0], FAST)
    assert rows[0]["brake_in_J"] == pytest.approx(rows[1]["brake_in_J"], rel=1e-12)
    assert rows[1]["t_f"] < rows[0]["t_f"]


# ------------------------------------------------------------------------------ run_study

def _study_file(tmp_path: Path, **extra) -> Path:
    st = {
        "name": "prueba",
        "base": str(REF),
        "base_overrides": FAST,
        "output_dir": str(tmp_path / "out"),
        "mesh": {"refine_factors": [0.5, 1, 2], "tolerance_percent": 1.0},
        "timestep": {"factors": [1, 2]},
        "sweeps": [{"name": "a", "key": "load.vehicle.decel", "values": [4.0, 6.0]}],
        **extra,
    }
    p = tmp_path / "study.yaml"
    p.write_text(yaml.safe_dump({"study": st}), encoding="utf-8")
    return p


def test_run_study_writes_everything(tmp_path):
    summary = run_study(_study_file(tmp_path), log=lambda *a: None)
    out = tmp_path / "out"
    for stem in ("mesh_study", "timestep_study", "sweep_a"):
        for ext in (".csv", ".md", ".pdf", ".png"):
            assert (out / f"{stem}{ext}").stat().st_size > 0, stem + ext
    for name in ("study.yaml", "base_config.yaml", "materiales.yaml", "study_summary.json"):
        assert (out / name).is_file()
    js = json.loads((out / "study_summary.json").read_text(encoding="utf-8"))
    assert js["mesh"]["selected_level"] == summary["mesh"]["selected_level"]
    # Los estudios posteriores usan la malla seleccionada.
    if summary["mesh"]["selected_level"] is not None:
        assert js["effective_overrides"]["mesh.refine_factor"] == \
            js["mesh"]["selected_refine_factor"]
    assert (out / "mesh_study.pdf").read_bytes()[:4] == b"%PDF"


@pytest.mark.parametrize(
    ("change", "match"),
    [
        ({"mesh": {"refine_factors": [2, 1]}}, "creciente"),
        ({"timestep": {"factors": [1, 1.5]}}, "enteros"),
        ({"sweeps": [{"name": "x", "values": [1]}]}, "'key'"),
        ({"base": "no_existe.yaml"}, "No existe"),
        ({"colores": 1}, "desconocidas"),
    ],
)
def test_load_study_errors(tmp_path, change, match):
    with pytest.raises(ValueError, match=match):
        load_study(_study_file(tmp_path, **change))


def test_cli_study_with_set(tmp_path):
    from discfem.cli import main

    p = _study_file(tmp_path)
    rc = main(["study", str(p), "--set", "study.sweeps=[]", "--set", "study.timestep.factors=[1]",
               "--set", "study.mesh.refine_factors=[0.5,1]"])
    assert rc == 0
    assert (tmp_path / "out" / "mesh_study.csv").is_file()
    assert not (tmp_path / "out" / "sweep_a.csv").exists()
    assert list(tmp_path.glob("tmp*.yaml")) == []  # el YAML temporal se borró
