# F0 · T1 — Inspección del repositorio y arranque del proyecto

Todas las pruebas pasaron; no se relajó nada.

## 1. Qué se pidió

Fase 0 de `PLAN.md` (§5, "Fase 0 — Inspección y arranque"): ver qué existe ya en el repositorio,
qué se puede reutilizar, qué choca con la estructura de la §2, y dejar el proyecto listo con
`pytest` corriendo (aunque sea con una prueba trivial) y `git` inicializado. `prompt.txt` estaba
vacío, así que se siguió el prompt de la Fase 0 que trae `PLAN.md`. El plan se aprobó antes de
crear archivos.

## 2. Qué se hizo

- `legacy/mesh/create_mesh.py`, `legacy/mesh/disc.msh`, `legacy/mesh/disc.svg`: movidos desde
  `mesh/` con `git mv`; quedan como referencia.
- `.gitignore`: nuevo; ignora `.venv/`, cachés, `results/` y `*.msh` (salvo los de `legacy/`).
- `pyproject.toml`: nuevo; paquete `discfem` con estructura `src/`, dependencias de CLAUDE.md
  y el extra `dev` (pytest, ruff).
- `src/discfem/__init__.py`: nuevo; docstring y `__version__`.
- `src/discfem/__main__.py`: nuevo; `python -m discfem` avisa que la CLI llega en la Fase 5 y
  sale con código 1.
- `tests/test_smoke.py`: nuevo; prueba de humo (versión del paquete + que las 6 dependencias
  se pueden importar).
- `PLAN.md`: sección ESTADO, Fase 0 marcada.
- `reportes/INDEX.md` y este reporte: nuevos.
- `CLAUDE.md`, `PLAN.md`, `Anteproyecto.tex`: ahora versionados (el `.tex` sin cambios).

## 3. Cómo funciona

En esta fase no hay algoritmo numérico. Lo que se montó es la infraestructura:

- **src-layout.** El código vive en `src/discfem/`. Por eso `import discfem` solo funciona si el
  paquete está instalado (`pip install -e .`): así las pruebas usan lo que vería un usuario, y no
  un archivo suelto del directorio actual por accidente.
- **Instalación editable** (`-e`). Instala un enlace al código fuente, así que los cambios en
  `src/` se ven sin reinstalar.
- **Prueba de humo.** Si gmsh o scipy no instalan bien (pasa a menudo en Windows), falla aquí y
  no en medio de la Fase 1.

Tabla ecuación ↔ código: no aplica en esta fase (no se implementó ninguna fórmula de la §3).

**Inventario de lo que había:**

| Elemento | Estado | Uso |
|---|---|---|
| `mesh/create_mesh.py` | script Gmsh de 46 líneas | se reutiliza la geometría y los campos `Distance`+`Threshold` |
| `mesh/disc.msh` (1132 nodos), `mesh/disc.svg` | artefactos generados | solo referencia visual |
| Código FEM, pruebas, pyproject | no existían | creados (salvo FEM) |

**Conflictos del script previo con PLAN.md §2 / CLAUDE.md** (se resuelven en la Fase 1):

1. Tiene 3 grupos físicos (`brake`, `conv`, `neumann`); CLAUDE.md exige 6 (`brake`,
   `top_inner`, `top_outer`, `r_inner`, `r_outer`, `symmetry`). `conv` junta 4 curvas.
2. Los parámetros están fijos en el código (la regla 5 exige YAML) y no hay `refine_factor`.
3. Escribe `disc.msh` en el directorio actual; el plan pide la API de Gmsh en memoria → `Mesh`.
4. `import meshio` sin usar; la variable `occ` apunta a `gmsh.model.geo` (nombre engañoso).
5. La carpeta `mesh/` de la raíz chocaba de nombre con `src/discfem/mesh/` → movida a `legacy/`.
6. Valores del campo de tamaño: SizeMin 4e-4, SizeMax 4e-3, DistMin 1e-3, DistMax = e. El YAML
   de la §4 usa `size_track: 3e-4`: no son iguales; manda el YAML.

**Anteproyecto.tex**: se confirmaron las inconsistencias que lista PLAN.md §1, en estas líneas:
50 ("isotérmica/aislada"), 226–227 (`\dd\Omega_0`), 269–271 (ventilado/perforado → "Modelo 3D"),
281 ("evolución tridimensional"). No se corrigieron (lo decidiste tú). El enlace de la línea 252
(`github.com/luismtorresv/disc-brake-fem`) es el `origin` de este mismo repositorio, así que la
duda 5 de la §1 queda resuelta: lo que había era solo el script de malla.

## 4. Decisiones

- **Mover a `legacy/` en vez de borrar o dejarlo**: lo elegiste tú. Conserva la referencia sin que
  se confunda con el paquete nuevo.
- **No corregir el `.tex` ahora**: lo elegiste tú; queda pendiente para la redacción.
- **Esqueleto mínimo**: solo `__init__` y `__main__`. No se crearon los módulos vacíos de la §2:
  cada uno nace en su fase con sus pruebas, para no dejar archivos sin probar.
- **setuptools** como backend de construcción: es el estándar y no añade dependencias. Descarté
  hatch y poetry porque no aportan nada en un proyecto de curso.
- **Versiones sin fijar** en `pyproject.toml`; el congelado va en la Fase 8
  (`requirements-lock.txt`).
- **Windows**: `python` no está en el PATH (es el alias de la Microsoft Store); se usa `py`
  (Python 3.13.3) para crear el venv.

## 5. Cómo verificarlo

```
py -m venv .venv
.venv\Scripts\python -m pip install -e ".[dev]"
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
.venv\Scripts\python -m discfem
```

Salida real (final de la instalación):

```
Successfully installed colorama-0.4.6 contourpy-1.4.0 cycler-0.12.1 discfem-0.0.1 fonttools-4.68.0 gmsh-4.15.2 iniconfig-2.3.1 kiwisolver-1.5.1 markdown-it-py-4.2.0 matplotlib-3.11.2 mdurl-0.1.2 meshio-5.3.5 numpy-2.5.4 packaging-26.3 pillow-13.0.0 pluggy-1.7.0 pygments-2.21.0 pyparsing-3.3.3 pytest-9.2.0 python-dateutil-4.0.0 pyyaml-6.0.3 rich-15.0.0 ruff-0.17.3 scipy-1.18.1 tzdata-2026.6
```

pytest, ruff y CLI:

```
.......                                                                  [100%]
7 passed in 2.25s
All checks passed!
CLI pendiente (Fase 5 de PLAN.md).
exit=1
```

## 6. Resultados

No hay resultados numéricos en esta fase. Ninguna prueba de la §6 (T1–T9) aplica todavía.
Versiones instaladas que importan para el solver: numpy 2.5.4, scipy 1.18.1, gmsh 4.15.2,
meshio 5.3.5, matplotlib 3.11.2, Python 3.13.3.

## 7. Limitaciones y riesgos

- Las versiones no están fijadas: una actualización de numpy/scipy/gmsh podría cambiar el
  comportamiento antes de la Fase 8.
- Solo se probó en Windows 10 con Python 3.13; Linux y Python 3.10 sin verificar.
- `legacy/mesh/create_mesh.py` no se ejecutó en esta sesión; su `disc.msh` se tomó tal cual.
- La rueda de gmsh trae binarios nativos; en otra máquina podría faltar alguna biblioteca del
  sistema (la prueba de humo lo detectaría).

## 8. Para revisar y para el informe

1. Por qué la estructura `src/` obliga a instalar el paquete, y qué hace `pip install -e`.
2. Qué se reutiliza del script previo (geometría y campos de tamaño) y qué tiene que cambiar
   (seis grupos de frontera, parámetros desde YAML, `refine_factor`).
3. Las cuatro correcciones pendientes del `.tex` (líneas 50, 226–227, 269–271, 281).
4. Para el informe: ninguna figura sale de esta fase; `legacy/mesh/disc.svg` sirve como imagen
   de la malla anterior si quieres comparar con la de la Fase 1.
