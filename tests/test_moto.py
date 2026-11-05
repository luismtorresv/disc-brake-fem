"""Pruebas del caso de motocicleta (configs/ejemplo_moto.yaml, Yamaha MT-07 de referencia)."""

from __future__ import annotations

import math
import warnings
from pathlib import Path

import pytest

from discfem.solver import run

ROOT = Path(__file__).resolve().parents[1]
MOTO = ROOT / "configs" / "ejemplo_moto.yaml"


@pytest.fixture(scope="module")
def moto():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return run(MOTO, {"mesh.refine_factor": 1.0}, write=False)


def test_moto_geometry_and_materials(moto):
    g = moto.config.geometry
    assert (g.r_outer, g.r_inner, g.half_thickness) == (0.149, 0.114, 0.00225)
    assert (g.r_a, g.r_b) == (0.1165, 0.1465)
    assert moto.materials["disc"].name == "stainless_steel_410"
    assert moto.pad.name == "sintered_pad_cu_zhou2024"
    assert moto.summary["unverified_materials"] == ["stainless_steel_410",
                                                    "sintered_pad_cu_zhou2024"]


def test_moto_energy(moto):
    """Energía por cara: φ m v0² / (2 n_d n_c); entra γ·E (salvo los 2 pasos de Euler del
    arranque de Rannacher); balance cerrado; ΔT final ≤ cota sin pérdidas (como T9)."""
    s, cfg = moto.summary, moto.config
    veh = cfg.load["vehicle"]
    E_face = (veh["brake_fraction"] * veh["mass"] * veh["v0"] ** 2
              / (2 * veh["discs_per_axle"] * veh["faces_per_disc"]))
    assert s["load"]["energy_per_face_J"] == pytest.approx(E_face, rel=1e-12)
    assert s["energy"]["brake_in_J"] == pytest.approx(moto.gamma * E_face, rel=1e-3)
    assert s["energy"]["balance_relative_residual"] < 1e-10
    disc, g = moto.materials["disc"], cfg.geometry
    cap = disc.rho_cp * math.pi * (g.r_outer**2 - g.r_inner**2) * g.half_thickness
    bound = moto.gamma * E_face / cap
    assert s["T_final"]["max"] - cfg.T0 <= bound


def test_moto_thermally_thin(moto):
    """e²/α ≪ t_f: la temperatura es casi uniforme en el espesor (superficie frente a plano
    medio en el centro de la pista, < 1 °C en el pico)."""
    d, g = moto.materials["disc"], moto.config.geometry
    assert g.half_thickness**2 / d.diffusivity < 0.25 * moto.load.t_f
    surf, mid = (p["T_max"] for p in moto.summary["probes"])
    assert 0 < surf - mid < 1.0
