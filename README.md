# discfem — calor transitorio en un disco de freno (FEM axisimétrico)

Solver de elementos finitos en Python para la ecuación de calor transitoria en la **sección
axisimétrica (r, z)** de un disco de freno sólido. Lo desarrollamos para el curso Métodos
Numéricos para EDP (EAFIT, 2026).

```
ρ c_p ∂T/∂t − k [ (1/r) ∂/∂r (r ∂T/∂r) + ∂²T/∂z² ] = 0
```

con flujo de fricción en la pista, convección en las caras expuestas y simetría en z = 0.

- **Malla:** triangular P1, con Gmsh o con una malla estructurada propia.
- **Tiempo:** esquema θ (Euler implícito, Crank–Nicolson, Rannacher).
- **Carga:** flujo de fricción derivado de la frenada (m, v0, a) y repartido entre disco y
  pastilla por efusividades.
- **Entradas:** todo se configura en archivos YAML.

La especificación completa (física, fórmulas, pruebas y fases) está en [`PLAN.md`](PLAN.md); las
reglas de trabajo, en [`CLAUDE.md`](CLAUDE.md); y los reportes de cada tarea, en
[`reportes/`](reportes/INDEX.md).

**Fuera de alcance:** discos ventilados, perforados o ranurados (requieren 3D), fuente móvil de
la pastilla, radiación, acoplamiento termomecánico, desgaste, resistencia de contacto y
propiedades dependientes de T.

---

## 1. Instalación

Requisitos: Python ≥ 3.10 (probado con 3.13) y conexión a internet para instalar las ruedas
(`gmsh` incluye sus binarios; no hace falta instalar Gmsh aparte).

### Windows (PowerShell)

```powershell
py -m venv .venv                     # en Windows `python` puede ser el alias de la Store; use `py`
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install -e ".[dev]"
python -m pytest -q                  # debe terminar con "... passed"
```

Si PowerShell bloquea `activate`, ejecute una vez
`Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`, o llame directamente a
`.venv\Scripts\python -m ...` sin activar.

### Linux / macOS

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[dev]"
python -m pytest -q
```

En Linux, la rueda de `gmsh` puede necesitar bibliotecas gráficas del sistema
(`libGLU`, `libXrender`, `libXcursor`, `libXft`, `libXinerama`). Si `import gmsh` falla,
instálelas con el gestor de paquetes (en Debian/Ubuntu:
`sudo apt install libglu1-mesa libxrender1 libxcursor1 libxft2 libxinerama1`).

### Versiones congeladas

`requirements-lock.txt` fija las versiones exactas con las que se validó el proyecto. Para
reproducir los resultados bit a bit:

```bash
pip install -r requirements-lock.txt
pip install -e . --no-deps
```

> La prueba `test_rectangle_unchanged` compara la malla y los resultados del caso rectangular
> **bit a bit** contra `tests/data/golden_rectangle.json`. Con otras versiones de `gmsh` o
> `numpy` puede fallar aunque el programa esté bien (la malla de Gmsh cambia entre versiones).
> Con `requirements-lock.txt` debe pasar.

---

## 2. Uso rápido

```bash
python -m discfem run   configs/ejemplo_auto.yaml          # corrida completa
python -m discfem mesh  configs/ejemplo_frenado.yaml       # solo malla + figura
python -m discfem study configs/estudio_malla.yaml         # convergencia de malla y Δt + barridos
python -m discfem study configs/estudio_pastilla.yaml      # barrido de pastilla
```

Cualquier clave del YAML se puede cambiar desde la línea de comandos con `--set` (repetible).
El valor se interpreta como YAML:

```bash
python -m discfem run configs/ejemplo_auto.yaml --set mesh.refine_factor=4 --set time.theta=1
python -m discfem run configs/ejemplo_auto.yaml --set "output.probes=[[0.11, 0.01]]"
python -m discfem study configs/estudio_malla.yaml --set "study.mesh.refine_factors=[1,2]"
```

Desde Python:

```python
from discfem.solver import run
res = run("configs/ejemplo_auto.yaml", {"load.vehicle.decel": 7.0}, write=False)
print(res.summary["T_max"])          # {'value': ..., 'time': ..., 'node': ..., 'r': ..., 'z': ...}
```

Archivos de configuración incluidos:

| Archivo | Qué es |
|---|---|
| `configs/ejemplo_frenado.yaml` | Caso de referencia del PLAN (materiales SUPUESTOS del prototipo, θ = 1). Define las pruebas T8 y T9; no usar para resultados. |
| `configs/ejemplo_auto.yaml` | Automóvil a 100 km/h con materiales verificados (preset `car`, Limpert), malla rf = 2, θ = 1/2 + 2 Rannacher. **Caso principal.** |
| `configs/ejemplo_perfil.yaml` | Sección poligonal (cubo interior más grueso) y h distinto en `r_outer`. |
| `configs/estudio_malla.yaml` | Estudio de malla (rf 1, 2, 4, 8), de Δt (Δt, Δt/2, Δt/4) y barridos de a, v0, m y k. |
| `configs/estudio_pastilla.yaml` | Pastillas del catálogo con disco de hierro y de acero 410, y barrido de γ. |
| `configs/materiales.yaml` | Librería de materiales y presets por tipo de vehículo. |

---

## 3. Cómo escribir un YAML de corrida

Unidades SI en todo (m, s, kg, W, K o °C). Ejemplo completo comentado:
[`configs/ejemplo_auto.yaml`](configs/ejemplo_auto.yaml). Si una clave no existe o un valor no
es válido, el programa se detiene con un mensaje en español que dice qué corregir.

### `geometry`

Sección rectangular (medio espesor; z = 0 es el plano de simetría):

```yaml
geometry:
  type: rectangle
  r_inner: 0.07            # m
  r_outer: 0.15
  half_thickness: 0.01     # semiespesor e
  track: {r_a: 0.08, r_b: 0.14}   # pista de fricción en z = e
  axisymmetric: true       # false ⇒ problema plano (solo validación)
```

Sección poligonal: vértices en sentido **antihorario**, con un grupo de frontera por lado (el
lado k va del punto k al k+1). `brake` debe ser una cadena contigua de lados horizontales.

```yaml
geometry:
  type: profile
  points: [[0.07, 0], [0.15, 0], [0.15, 0.01], [0.14, 0.01], [0.08, 0.01], [0.07, 0.01]]
  segments: [symmetry, r_outer, top_outer, brake, top_inner, r_inner]
```

Regiones con otro material (rectángulo o perfil). Lo que no cubren es la región `disc`:

```yaml
  regions:
    capa: [[0.08, 0.009], [0.14, 0.009], [0.14, 0.01], [0.08, 0.01]]   # capa de 1 mm bajo la pista
```

Grupos de frontera: `brake` (pista), `top_inner`, `top_outer`, `r_inner`, `r_outer` y `symmetry`
(z = 0). Las regiones y los perfiles requieren `mesh.backend: gmsh`.

### `mesh`

```yaml
mesh:
  backend: gmsh            # gmsh | structured (solo rectángulo sin regiones)
  refine_factor: 2.0       # divide todos los tamaños; 2 ⇒ ~4 veces más celdas
  size_track: 3.0e-4       # tamaño junto a la pista [m]
  size_far: 4.0e-3         # tamaño lejos [m]
  grading_distance: 0.01   # distancia en que el tamaño crece de size_track a size_far [m]
```

### `materials`

```yaml
materials:
  library: configs/materiales.yaml
  preset: car              # car | motorcycle: define disco y pastilla
  regions: {disc: stainless_steel_410}   # opcional; tiene prioridad sobre el preset
  pad: pad_limpert         # opcional; tiene prioridad sobre el preset
  custom:                  # materiales definidos aquí mismo (mismas reglas que la librería)
    mi_acero: {k: 30.0, rho: 7800.0, cp: 470.0, source: "Ficha X, tabla Y", verified: false}
  overrides: {disc: {k: 48.0}}   # cambia propiedades de una región (queda como no verificado)
```

Con varias regiones, `regions` debe nombrar el material de cada una (`disc`, `capa`, …).

### `boundaries`, `convection`, `initial`

```yaml
boundaries:                # convection | insulated | flux, para cada uno de los 6 grupos
  brake: flux
  top_inner: convection
  top_outer: convection
  r_inner: convection
  r_outer: convection
  symmetry: insulated
convection:
  h: 40.0                  # W/(m²·K), para todos los grupos con convección
  T_inf: 25.0
  per_group:               # opcional: h y/o T∞ distintos por grupo
    r_outer: {h: 60.0}
initial: {T0: 25.0}
```

### `load` (flujo en los grupos `flux`)

```yaml
load:
  type: braking            # braking | constant_flux | table
  vehicle: {mass: 1500.0, v0: 27.78, decel: 5.0, wheel_radius: 0.30,
            brake_fraction: 0.7, discs_per_axle: 2, faces_per_disc: 2}
  pressure_profile: uniform_wear   # uniform_wear | uniform_pressure
  gamma: auto              # auto (por efusividades disco/pastilla) | número en (0, 1)
```

Las otras dos opciones dan directamente el flujo que entra al disco (sin γ):
`{type: constant_flux, q: 2.0e5, duration: 3.0}` (W/m², s) y
`{type: table, points: [[0, 0], [1, 1.0e5], [3, 0]]}` (interpolación lineal; 0 después del
último instante).

### `time` y `output`

```yaml
time:
  theta: 0.5               # 1 = Euler implícito; 0.5 = Crank–Nicolson
  rannacher_steps: 2       # pasos iniciales con θ = 1 (amortigua el arranque brusco del flujo)
  dt: auto                 # auto = t_f / steps_per_tf (t_f queda exacto en la malla temporal)
  steps_per_tf: 100
  t_end: 60.0
output:
  dir: results/mi_caso
  probes: [[0.11, 0.01], [0.11, 0.0]]   # sondas (r, z)
  formats: [csv, vtk, png]
  save_every: 1            # cada cuántos pasos se guarda el campo completo (sondas: siempre)
```

---

## 4. Cómo añadir un material

1. Abra [`configs/materiales.yaml`](configs/materiales.yaml) y añada una entrada en `materials:`:

   ```yaml
   materials:
     mi_material:
       description: Texto libre
       k: 30.0          # W/(m·K)
       rho: 7800.0      # kg/m³
       cp: 470.0        # J/(kg·K)
       source: "Autor, título, tabla y página (o URL de la ficha)"
       verified: false  # true solo si contrastó los valores con la fuente
   ```

   `source` y `verified` son **obligatorios**. Si falta alguno, el programa se detiene. Mientras
   `verified: false`, cada corrida imprime `ADVERTENCIA (UnverifiedMaterialWarning)`. Si un
   valor está fuera de los rangos de sanidad (por ejemplo, k ∉ (0,01; 500)), también advierte.
   **No escriba valores "típicos" de memoria:** cite siempre la fuente.

2. Úselo en una corrida con `materials.regions.disc: mi_material` o `materials.pad: mi_material`,
   o créele un preset:

   ```yaml
   presets:
     mi_vehiculo: {description: "...", disc: mi_material, pad: pad_limpert}
   ```

Materiales actuales:

| Nombre | Uso | Fuente | Verificado |
|---|---|---|---|
| `gray_cast_iron_limpert` | disco (auto) | Limpert, *Brake Design and Safety*, 2.ª ed., Tabla 3-1A | sí |
| `pad_limpert` | pastilla orgánica (base asbesto) | Limpert, Tabla 3-1A | sí |
| `lining_limpert` | forro de tambor | Limpert, Tabla 3-1A | sí |
| `stainless_steel_410` | disco (moto) | Ficha Chilexpo "Aleación 410" (c_p con unidad errónea en la ficha) | no |
| `sintered_pad_cu_zhou2024` | pastilla sinterizada de cobre (moto) | k y ρ: Zhou et al., *Metals* 14:1048 (2024); c: Goo y Na (2014) | no |
| `sintered_pad_goo2014` | pastilla sinterizada (referencia) | Goo y Na, *J. KSTLE* 30(3) (2014) | no |
| `gray_cast_iron`, `organic_pad` | solo para el caso de referencia del PLAN | supuestos del prototipo | no |

Presets: `car` = `gray_cast_iron_limpert` + `pad_limpert`; `motorcycle` =
`stainless_steel_410` + `sintered_pad_cu_zhou2024`. Los detalles y las reservas de cada fuente
están en `materiales.yaml` y en los reportes `F7_T2` a `F7_T5`.

---

## 5. Cómo leer las salidas

`run` escribe en `output.dir`:

| Archivo | Contenido |
|---|---|
| `config.yaml` | Configuración **efectiva** (YAML + `--set`), con la librería apuntando a la copia local. Con esta carpeta sola se repite la corrida: `python -m discfem run results/mi_caso/config.yaml`. |
| `materiales.yaml` | Copia de la librería usada. |
| `summary.json` | T_max (valor, instante, posición), T al final, sondas, γ, energía (entrada, convección, almacenada, en J **con 2π**) y cierre del balance, malla, paso de tiempo, materiales con su fuente, materiales no verificados y tiempo de cómputo. Dos corridas iguales dan el mismo archivo salvo `compute_time_s`. |
| `tmax.csv` | t, T_max, nodo y (r, z) del máximo en cada paso. |
| `probes.csv` | T en cada sonda en cada paso. |
| `field.png` | Campo T(r, z) en 4 instantes (t_f/4, pico, t_f, t_end), reflejado en z para ver el espesor completo. |
| `tmax.png` | T_max(t) y las sondas; línea vertical en t_f. |
| `depth_profiles.png` | Perfiles T(z) en el centro de la pista. |
| `mesh.png` | Malla, grupos de frontera y regiones. |
| `field_XXXX.vtu`, `field.pvd` | Los mismos 4 instantes para ParaView (abrir `field.pvd`). |

El **balance de energía** debe cerrar a ~1e-14. Si no, algo está mal. Esa cifra no mide la
precisión del resultado: para eso están los estudios de malla y de Δt.

`study` escribe en `study.output_dir`: tablas `*.csv` y `*.md`; figuras `*.pdf` (vectoriales,
para el informe) y `*.png`; `study_summary.json`; y copias de la config base, del estudio y de la
librería. En el estudio de malla, `eps_percent` es el cambio relativo de T_max entre niveles y
`order_p` el orden observado.

---

## 6. Pruebas y validación

```bash
python -m pytest -q          # 240 pruebas, ~20 s
python -m ruff check src tests
```

| Prueba | Qué valida | Resultado |
|---|---|---|
| T1 | Malla: Σ A·r̄ = ∫ r dA exacto, grupos de frontera, orientación, refinamiento | 1e-16 |
| T2 | Matrices elementales frente a cuadratura de alto orden | 1e-12 |
| T3 | Flujo de fricción: 2π ∫ q r dr = γ P_f | 1e-16 |
| T4 / T5 | Balance de energía / estado uniforme | ~1e-14 / 9e-14 |
| T6 | Modo propio axisimétrico (Bessel J0), nx = 10…80 | 0,86 → 0,03 % |
| T7 | Medio semi-infinito, flujo lineal en t | θ = 1: 0,4 %; CN + Rannacher: 0,07 % |
| T8 / T9 | Caso de referencia y cota energética | 102,30 °C en 3,39 s; 47,3 ≤ 48,1 K |

Estudios: la malla converge con ε < 1 % desde el segundo nivel (0,18 → 0,067 → 0,013 %), y en Δt
el orden observado es 1,00 para θ = 1. Detalles en `reportes/`.

---

## 7. Estructura del código

```
src/discfem/
  geometry.py       secciones rectangle/profile, regiones, utilidades de polígonos
  mesh/             Mesh (contrato), structured.py, gmsh_backend.py
  fem/assembly.py   M, K, H, F vectorizados (peso r o 1; sin el factor 2π)
  loads.py          q_fric(r, t): braking | constant_flux | table
  materials.py      librería, presets, efusividad y γ, advertencias
  timestepping.py   esquema θ, Rannacher, malla temporal con t_f exacto, balance
  config.py         validación del YAML y overrides con claves con puntos
  solver.py         run(config, overrides) -> Result
  postprocess.py    sondas, CSV, JSON, figuras, VTK
  convergence.py    estudios de malla y Δt, barridos, figura γ
  cli.py            python -m discfem run | mesh | study
```

## Referencias de datos

- R. Limpert, *Brake Design and Safety*, 2.ª ed., SAE. Tabla 3-1A.
- Y. Zhou et al., "Effect of Sintering Temperature on the Microstructure and Mechanical and
  Tribological Properties of Copper Matrix Composite for Brake Pads", *Metals* 14:1048 (2024),
  doi:10.3390/met14091048.
- B.-C. Goo, I.-K. Na, "Topology Optimization of Railway Brake Pad by Contact Analysis",
  *J. Korean Soc. Tribol. Lubr. Eng.* 30(3):177–182 (2014), doi:10.9725/kstle.2014.30.3.177.
- Chilexpo SpA, ficha técnica "Aleación 410" (UNS S41000).
