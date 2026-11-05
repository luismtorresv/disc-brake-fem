# F8 · T1 — Documentación (README) y reproducibilidad (requirements-lock, entorno limpio)

Todas las pruebas pasaron (240 passed), tanto en el entorno de desarrollo como en dos entornos
limpios. **No se verificó en Linux**: las instrucciones de Linux del README siguen la práctica
estándar, pero no se ejecutaron (§7).

## 1. Qué se pidió

Fase 8 de `PLAN.md` (§5, "Fase 8 — Documentación y reproducibilidad"):

- escribir `README.md` con la instalación en Windows y Linux, ejemplos de cada comando, cómo
  editar un YAML, cómo añadir un material y cómo leer las salidas;
- congelar las versiones con `pip freeze` en `requirements-lock.txt`;
- verificar en un entorno limpio que `pip install -e ".[dev]"` y `pytest -q` funcionan.

## 2. Qué se hizo

- `README.md`: nuevo. Contiene:
  - instalación (Windows con `py`, Linux/macOS, versiones congeladas);
  - uso rápido de `run`, `mesh` y `study`, con `--set` y desde Python;
  - referencia de todas las secciones del YAML;
  - cómo añadir un material, con la tabla de materiales y su estado de verificación;
  - cómo leer las salidas de `run` y `study`;
  - resumen de pruebas y validación, estructura del código y referencias de datos.
- `requirements-lock.txt`: nuevo. `pip freeze --exclude-editable` del venv de desarrollo, con
  una cabecera que indica la versión de Python, el sistema y la fecha.
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

No hay algoritmo nuevo. La reproducibilidad descansa en tres niveles:

1. **Corrida:** cada carpeta de resultados guarda la configuración efectiva y la librería de
   materiales (Fase 5). `python -m discfem run results/x/config.yaml` repite la corrida.
2. **Entorno:** `requirements-lock.txt` fija las versiones exactas. Importa porque la malla de
   Gmsh puede cambiar entre versiones, y `test_rectangle_unchanged` compara bit a bit.
3. **Instalación:** `pip install -e ".[dev]"` (versiones libres) o `pip install -r
   requirements-lock.txt` seguido de `pip install -e . --no-deps` (versiones fijas).

Tabla "elemento ↔ dónde": no aplica (no se implementó ninguna fórmula).

## 4. Decisiones

- **README en español,** como los mensajes y los reportes (CLAUDE.md, regla 6). Cada clave que
  documenta se copió de `config.py` y de la ayuda real de la CLI, no de memoria.
- **El README distingue `ejemplo_frenado.yaml` de `ejemplo_auto.yaml`.** El primero existe para
  las pruebas T8 y T9 con materiales supuestos; el segundo es el caso principal con materiales
  verificados.
- **Advertencia sobre la prueba bit a bit.** Con otras versiones de gmsh o numpy,
  `test_rectangle_unchanged` puede fallar sin que el programa esté mal. El README lo explica y
  remite al archivo de versiones congeladas.
- **Lock con `--exclude-editable`.** El propio paquete no debe aparecer fijado; se instala con
  `-e . --no-deps`.
- **Entornos de verificación en el scratchpad,** para no tocar el `.venv` del proyecto.

## 5. Cómo verificarlo

Entorno limpio (A), instalación normal. Salida real:

```
py -m venv venvA
venvA\Scripts\python -m pip install -e ".[dev]"
Successfully installed colorama-0.4.6 contourpy-1.4.0 cycler-0.12.1 discfem-0.0.1 fonttools-4.68.0 gmsh-4.15.2 iniconfig-2.3.1 kiwisolver-1.5.1 markdown-it-py-4.2.0 matplotlib-3.11.2 mdurl-0.1.2 meshio-5.3.5 numpy-2.5.4 packaging-26.3 pillow-13.0.0 pluggy-1.7.0 pygments-2.21.0 pyparsing-3.3.3 pytest-9.2.0 python-dateutil-4.0.0 pyyaml-6.0.3 rich-15.0.0 ruff-0.17.3 scipy-1.18.1 tzdata-2026.6
venvA\Scripts\python -m pytest -q
240 passed in 25.79s
A: versiones idénticas al lock
```

Entorno limpio (B), desde el lock. Salida real:

```
venvB\Scripts\python -m pip install -r requirements-lock.txt
Successfully installed PyYAML-6.0.3 Pygments-2.21.0 colorama-0.4.6 ... scipy-1.18.1 tzdata-2026.6
venvB\Scripts\python -m pip install -e . --no-deps
Successfully installed discfem-0.0.1
venvB\Scripts\python -m pytest -q
240 passed in 24.98s
venvB\Scripts\python -m pip check
No broken requirements found.
```

Ejemplos del README ejecutados en el entorno B (salida real, recortada):

```
python -m discfem mesh configs/ejemplo_frenado.yaml
Malla 'gmsh' (refine_factor = 1)
  nodos: 929   triángulos: 1618
python -m discfem run configs/ejemplo_auto.yaml --set mesh.refine_factor=4 --set time.theta=1
Caso 'ejemplo_auto_100kmh': 11051 nodos, 21148 triángulos, 1080 pasos (θ = 1, Rannacher = 2, Δt = 0.05556 s)
  T_max = 106.300 en t = 3.500 s (r = 0.1095, z = 0.0100)
run('configs/ejemplo_auto.yaml', {'load.vehicle.decel': 7.0}, write=False).summary['T_max']
{'value': 117.80599576918554, 'time': 2.182714285714286, 'node': 261, 'r': 0.10969694631213703, 'z': 0.01}
python -m discfem run configs/ejemplo_auto.yaml --set "output.probes=[[0.11, 0.01]]"
  T_max = 106.695 en t = 3.500 s (r = 0.1095, z = 0.0100)
python -m discfem run <carpeta>/config.yaml        # repetición desde la config copiada
  T_max = 106.695 en t = 3.500 s (r = 0.1095, z = 0.0100)
python -m discfem study configs/estudio_malla.yaml --set "study.mesh.refine_factors=[1,2]" ...
Salidas del estudio: ...\outS
```

## 6. Resultados

| Criterio de la Fase 8 | Resultado |
|---|---|
| README: instalación Windows y Linux, cada comando, YAML, materiales, salidas | hecho (§2) |
| `requirements-lock.txt` con `pip freeze` | hecho, 24 paquetes |
| Entorno limpio con `pip install -e ".[dev]"` + `pytest -q` | 240 passed |
| Entorno limpio desde el lock | 240 passed; `pip check` sin conflictos |
| Repetir una corrida desde su `config.yaml` | mismo T_max (106,695 °C) |

## 7. Limitaciones y riesgos

- **Solo se verificó en Windows 10 con Python 3.13.3.** Las instrucciones de Linux (venv, y las
  bibliotecas que puede necesitar la rueda de gmsh) no se ejecutaron. Tampoco se probó Python
  3.10–3.12, aunque `requires-python = ">=3.10"`.
- **El lock es específico de plataforma solo en lo que pip resuelva.** Los 24 paquetes tienen
  ruedas para Linux y macOS, pero no se comprobó.
- **`test_rectangle_unchanged` depende de las versiones** (está documentado en el README).
- **Datos pendientes que el README declara:** materiales de moto no verificados (410, pastillas
  sinterizadas) y la falta de un caso de moto con geometría propia.

## 8. Para revisar y para el informe

1. **Los tres niveles de reproducibilidad:** corrida (`config.yaml` efectivo y librería
   copiada), entorno (`requirements-lock.txt`) e instalación (`-e .`).
2. **Por qué la prueba bit a bit necesita versiones fijas:** la malla de Gmsh depende de la
   versión.
3. **Qué caso usar en el informe:** `ejemplo_auto.yaml`, con materiales verificados. El de
   referencia solo sirve para las pruebas.
4. **El README como guía de uso** para quien reciba el código: instalación, YAML, materiales y
   salidas.
5. **Para el informe:** la tabla de validación del README (§6) resume T1–T9 en una sola vista.
