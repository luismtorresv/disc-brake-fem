"""Prueba de humo de la Fase 0: el paquete y las dependencias numéricas importan."""

import importlib

import pytest

import discfem


def test_version():
    assert discfem.__version__ == "0.0.1"


@pytest.mark.parametrize("module", ["numpy", "scipy.sparse", "yaml", "gmsh", "meshio", "matplotlib"])
def test_dependency_imports(module):
    importlib.import_module(module)
