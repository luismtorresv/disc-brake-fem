# F5 · T1 — Configuración, solver, salidas y CLI (`run`), con T8 y T9

Todas las pruebas pasaron (174 passed) sin relajar tolerancias. Dos cosas para leer primero:

- **Interpretación de "mismo summary.json".** `summary.json` incluye el tiempo de cómputo, como
  pide el PLAN, así que dos corridas no dan un archivo idéntico. Lo interpreté como "idéntico
  salvo `compute_time_s`", y así lo prueba la suite (§4).
- **Sesgo en los valores de referencia de T8.** Los valores de T8 llevan incorporado el déficit
  de energía de $\theta = 1$. Con un esquema más exacto, $T_{\max}(60\ \mathrm{s}) = 72{,}8028$ °C,
  fuera de la tolerancia por 0,003 °C (§6 y §7).

## 1. Qué se pidió

Fase 5 de `PLAN.md` (§5, "Fase 5 — Configuración, CLI y salidas"):

- `config.py` (dataclasses + validación con mensajes en español), `solver.py`
  (`run(config, overrides=None) -> Result`), `postprocess.py` y `cli.py`.
- Salidas por corrida: `config.yaml` copiado, `tmax.csv`, `probes.csv`, `summary.json`, figuras
  PNG (campo en 4 instantes reflejado en $z$, $T_{\max}(t)$, perfiles en profundidad) y VTK.
- `overrides` con claves con puntos.
- Pruebas T8 y T9.

Criterios de aceptación: `python -m discfem run configs/ejemplo_frenado.yaml` produce todas las
salidas, y correr dos veces da el mismo `summary.json`.

## 2. Qué se hizo

- `src/discfem/config.py`: nuevo. `Config`, `MaterialsConfig`, `TimeConfig`, `OutputConfig`,
  `set_dotted`, `build_config` y `load_config`. Rechaza claves desconocidas y valores inválidos
  con mensajes en español.
- `src/discfem/solver.py`: nuevo. `run`, `Result` y `make_summary`. Llama a `warn_unverified`
  al inicio de cada corrida (pendiente de la Fase 4).
- `src/discfem/postprocess.py`: ampliado. `interpolation_matrix` (sondas), `snapshot_indices`,
  `plot_field_snapshots`, `plot_tmax`, `plot_depth_profiles`, `write_vtk` y `write_outputs`.
- `src/discfem/cli.py`: reescrito. `run` con `--set clave=valor` repetible; advertencias en
  formato `ADVERTENCIA (Categoría): …`; los errores de configuración salen con código 2 y
  mensaje en español. `mesh` sigue validando solo `geometry`, `mesh` y `output`.
- `src/discfem/loads.py`: `brake_load_vector` acepta los grupos con `flux` (por defecto
  `brake`).
- `pyproject.toml`: `[tool.ruff.lint] extend-ignore = ["TRY004"]`, con el motivo en un
  comentario. Se quitaron los dos `noqa` de la Fase 4.
- `tests/test_solver.py`: nuevo. T8, T9, energía, salidas, reproducibilidad, overrides, sondas,
  configuración y CLI (25 casos).
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

`run` encadena lo construido en las fases anteriores:

1. **Configuración.** `load_config` lee el YAML, aplica los overrides (`"mesh.refine_factor": 2`
   ⇒ `raw["mesh"]["refine_factor"] = 2`, creando las secciones intermedias si faltan) y valida
   cada sección.
2. **Materiales.** Se carga la librería (relativa al directorio actual o, si no está ahí, al de
   la config), se aplican los overrides por región y se emite `UnverifiedMaterialWarning` si
   corresponde. También se avisa si $h$ está fuera del rango de sanidad.
3. **Malla, $\gamma$ y carga.** $\gamma$ automático por efusividades (Fase 4). La carga se separa
   en $f(t)\,\mathbf{F}_s$, con $\mathbf{F}_s$ ensamblado una vez.
4. **Ensamblaje.** $\mathbf{K}$ y $\mathbf{M}$ con $k$ y $\rho c$ por región, $\mathbf{H}$ sobre
   los grupos `convection`, y $\mathbf{F}_{\text{conv}} = \mathbf{H}\,(T_\infty \mathbf{1})$. Esto
   usa la identidad $\int h\, T_\infty\, \varphi_i\, r\, ds = (\mathbf{H}\, T_\infty \mathbf{1})_i$,
   probada en la Fase 2.
5. **Tiempo.** `time_grid` contiene $t_f$ exactamente; `integrate` aplica el esquema $\theta$
   (Fase 3).
6. **Sondas.** `interpolation_matrix` construye una matriz dispersa $\mathbf{P}$
   ($n_{\text{sondas}} \times N$): para cada sonda, las coordenadas baricéntricas del triángulo
   que la contiene. Así $T_{\text{sonda}}(t) = \mathbf{P}\,\mathbf{T}(t)$ es la misma
   interpolación P1 que usa el método, y es exacta para campos lineales (probado).
7. **Resumen y salidas.** Las energías se reportan en J con el $2\pi$ explícito (CLAUDE.md).

**Reproducibilidad (regla 5).** `config.yaml` en la carpeta de resultados no es una copia
literal del YAML original, sino la configuración **efectiva** (original + overrides). Su
`materials.library` apunta a `materiales.yaml`, que también se copia a la carpeta. Con esa
carpeta sola se repite la corrida, y hay una prueba que lo verifica.

**Figuras.** Se dibujan 4 instantes: $t_f/4$, el instante de $T_{\max}$, $t_f$ y
$t_{\text{end}}$. El campo se refleja respecto a $z = 0$, de modo que se ve el espesor completo
del disco (las dos caras calentadas), con la misma escala de color en todos los paneles. Los
perfiles en profundidad se toman en el centro de la pista. Los VTK son esos mismos 4 instantes,
con un `field.pvd` para abrirlos como serie temporal en ParaView.

| Elemento (PLAN.md) | Código |
|---|---|
| `run(config, overrides=None) -> Result` | `solver.py:66`; `Result`, `:40` |
| Overrides con puntos | `set_dotted`, `config.py:86`; `load_config`, `:225` |
| Validación con mensajes en español | `build_config`, `config.py:132` |
| Advertencia de materiales no verificados en cada corrida | `solver.py:96` |
| $F_{\text{conv},i} = \int h\, T_\infty\, \varphi_i\, r\, ds$ | `solver.py:118` |
| Malla temporal con $t_f$, esquema $\theta$ | `solver.py:123-128` |
| Energías con $2\pi$ explícito | `make_summary`, `solver.py:152`; `two_pi`, `:157` |
| Sondas (interpolación P1) | `interpolation_matrix`, `postprocess.py:89`; `solver.py:133` |
| config.yaml efectivo + materiales.yaml, CSV, JSON, PNG, VTK | `write_outputs`, `postprocess.py:222` |
| Campo reflejado en $z$, 4 instantes | `plot_field_snapshots`, `postprocess.py:130`; `snapshot_indices`, `:120` |
| CLI `run --set` | `cmd_run`, `cli.py:77`; `_parse_set`, `:25` |

## 4. Decisiones

- **"Mismo `summary.json`".** El PLAN pide, en `summary.json`, el tiempo de cómputo y también
  que dos corridas den el mismo archivo. Ambas cosas no pueden cumplirse literalmente. Dejé
  `compute_time_s` en `summary.json`, como lo lista el PLAN, y la prueba compara todo lo demás.
  Comprobado: el `diff` entre dos corridas solo muestra esa línea (§5). La alternativa era
  sacarlo a un `timing.json`; si la prefieres, es un cambio pequeño.
- **`config.yaml` efectivo en vez de copia literal.** Una copia literal perdería los overrides
  y la corrida no se podría reproducir solo con ella. Lo malo es que se pierden los comentarios
  del YAML original.
- **Se guardan todos los pasos de tiempo.** Las sondas y $T_{\max}(t)$ los necesitan con
  resolución completa. Para 929 nodos y 1080 pasos son ~8 MB en memoria. No se expuso
  `save_every` en el YAML; queda para la Fase 6 si las mallas finas lo piden.
- **`mesh` no exige la configuración completa.** Al principio lo cambié para que usara
  `load_config` y rompí la prueba de la Fase 1, que usa un YAML mínimo. Para "solo malla +
  figura" no hace falta validar materiales ni carga, así que se revirtió a validar solo
  `geometry`, `mesh` y `output`, ahora también con `--set`.
- **Una sola región `disc` en la geometría rectangular.** `build_config` rechaza otras regiones
  con un mensaje claro. Varias regiones llegan en la Fase 7.
- **TRY004 desactivada en `pyproject.toml`** en vez de `noqa` repetidos. Todos los errores de
  configuración son `ValueError` y la CLI los atrapa en un único sitio (código de salida 2).
- **$\theta$ del caso de referencia.** No se tocó: el YAML del PLAN fija `theta: 1.0` y
  `rannacher_steps: 0`, y con eso se prueban T8 y T9. El valor por defecto del código
  ($\theta = 1/2$ + 2 Rannacher, Fase 3) solo aplica si el YAML no dice nada.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
.venv\Scripts\python -m discfem run configs/ejemplo_frenado.yaml
```

Salida real de pytest y ruff:

```
174 passed in 8.92s
All checks passed!
```

Salida real de `run` (stderr + stdout):

```
ADVERTENCIA (UnverifiedMaterialWarning): Materiales NO verificados: gray_cast_iron (source: PENDIENTE: citar fuente); organic_pad (source: PENDIENTE: citar fuente). Los resultados no sirven para el informe hasta citar fuentes y poner verified: true.
Caso 'ejemplo_frenado_100kmh': 929 nodos, 1618 triángulos, 1080 pasos (θ = 1, Rannacher = 0, Δt = 0.05556 s)
  γ = 0.8692   t_f = 5.556 s
  T_max = 102.297 en t = 3.389 s (r = 0.1091, z = 0.0100)
  T_max(t_end = 60 s) = 72.334
  energía [J, con 2π]: entra 87160.4, sale por convección 2606.0, almacenada 84554.4; cierre del balance 4.9e-14
  sonda (r=0.11, z=0.01): T_max = 102.292 en t = 3.389 s; T_final = 71.887
  sonda (r=0.11, z=0.005): T_max = 87.791 en t = 5.778 s; T_final = 71.888
  sonda (r=0.11, z=0): T_max = 86.881 en t = 7.112 s; T_final = 71.892
  Materiales NO verificados: gray_cast_iron, organic_pad
  tiempo de cómputo: 0.20 s   salidas: results\ejemplo_frenado
```

Archivos producidos (salida real de `ls results/ejemplo_frenado`):

```
config.yaml depth_profiles.png field.png field.pvd field_0000.vtu field_0001.vtu field_0002.vtu field_0003.vtu materiales.yaml mesh.png probes.csv summary.json tmax.csv tmax.png
```

Dos corridas seguidas, `diff` de los `summary.json` (salida real):

```
93c93
<   "compute_time_s": 0.20455219999894325
---
>   "compute_time_s": 0.2078867000000173
```

Comparación de esquemas y mallas (salida real del script exploratorio previo a la fase, que usa
las mismas funciones que `run`):

```
gmsh rf=1 N=929 θ=1.0 R=0: Tmax=102.30 en t=3.389; Tmax(60)=72.33; ΔTmedio final=46.17 cota=48.08; ΔT(max nodo) final=47.33 min=44.17  0.18s
gmsh rf=1 N=929 θ=0.5 R=2: Tmax=102.73 en t=3.445; Tmax(60)=72.80; ΔTmedio final=46.63 cota=48.08; ΔT(max nodo) final=47.80 min=44.61  0.17s
gmsh rf=2 N=3028 θ=1.0 R=0: Tmax=102.48 en t=3.389; Tmax(60)=72.34; ΔTmedio final=46.17 cota=48.08; ΔT(max nodo) final=47.34 min=44.17  0.43s
gmsh rf=2 N=3028 θ=0.5 R=2: Tmax=102.91 en t=3.389; Tmax(60)=72.80; ΔTmedio final=46.63 cota=48.08; ΔT(max nodo) final=47.80 min=44.61  0.45s
structured rf=1 N=2190 θ=1.0 R=0: Tmax=102.29 en t=3.389; Tmax(60)=72.33; ΔTmedio final=46.17 cota=48.08; ΔT(max nodo) final=47.33 min=44.17  0.14s
structured rf=1 N=2190 θ=0.5 R=2: Tmax=102.72 en t=3.445; Tmax(60)=72.80; ΔTmedio final=46.63 cota=48.08; ΔT(max nodo) final=47.80 min=44.61  0.14s
```

## 6. Resultados

La cota de T9 supone que toda la energía que entra a una cara se reparte uniformemente en el
semidisco, sin pérdidas:

$$
\Delta T_{\text{final}} \le \frac{\gamma\, E_{\text{cara}}}{\rho c\, \pi\,(r_o^2 - r_i^2)\, e}
$$

| Prueba | Criterio (§6) | Medido aquí (YAML de referencia, $\theta = 1$, Gmsh 929 nodos) | Prototipo |
|---|---|---|---|
| T8: $T_{\max}$ | $102{,}5 \pm 0{,}5$ °C | 102,297 °C | 102,5 |
| T8: instante de $T_{\max}$ | $3{,}4 \pm 0{,}1$ s | 3,389 s | 3,4 |
| T8: $T_{\max}(60\ \mathrm{s})$ | $72{,}3 \pm 0{,}5$ °C | 72,334 °C | 72,3 |
| T9: cota $\gamma E_{\text{cara}} / (\rho c\, \pi (r_o^2 - r_i^2)\, e)$ | 48,1 K | 48,08 K | 48,1 |
| T9: $\Delta T_{\text{final}} = T_{\max}(60) - T_0$ | $\le$ cota | 47,33 K | 47,3 |
| $\Delta T$ medio final (energía almacenada / capacidad) | $\le$ cota (adicional) | 46,17 K | — |
| Balance de energía | (control) | $4{,}9\times10^{-14}$ | — |

Lectura física:

- Entran 87,2 kJ por cara ($\gamma \cdot 101{,}3\ \mathrm{kJ} \cdot 0{,}99$), salen 2,6 kJ por
  convección en 60 s y quedan almacenados 84,6 kJ.
- La cota de T9 supone cero pérdidas y temperatura uniforme. La diferencia con el $\Delta T$
  medio ($48{,}08 - 46{,}17 = 1{,}9$ K) se reparte así: ~1 % por la regla $\theta = 1$ (0,46 K) y
  el resto por la convección.
- Las sondas muestran el retraso de la difusión: el pico llega a la superficie en 3,39 s, al
  plano medio en 7,11 s, y a los 60 s el espesor está prácticamente uniforme (71,89 °C).

**Sesgo de los valores de referencia.** Con $\theta = 1/2$ + 2 Rannacher (el esquema más exacto
según T7), la misma malla da $T_{\max} = 102{,}727$ °C y $T_{\max}(60\ \mathrm{s}) = 72{,}8028$
°C. Este último valor queda FUERA de $72{,}3 \pm 0{,}5$ por 0,003 °C. La diferencia es casi toda
energía: $\theta = 1$ con $\Delta t = t_f/100$ mete un 1 % menos (Fase 3). El $\Delta T$ medio
final pasa de 46,17 a 46,63 K, y $46{,}17/46{,}63 = 0{,}990$. Así que 72,3 °C no es "la" respuesta
del modelo continuo, sino la respuesta con $\theta = 1$ y $\Delta t = t_f/100$. La respuesta del
modelo continuo está más cerca de 72,8 °C. No se ajustó nada.

## 7. Limitaciones y riesgos

- **T8 es sensible al esquema.** Si el YAML de referencia se cambiara a $\theta = 1/2$ +
  Rannacher, T8 fallaría por 0,003 °C en $T_{\max}(60\ \mathrm{s})$ (72,8028 frente al límite de
  72,8). No sería un error del código, sino que la referencia del PLAN está calculada con
  $\theta = 1$. Recomiendo, para el informe, reportar la referencia junto con el esquema que la
  produce, o recalcularla con $\Delta t$ menor.
- **Materiales sin verificar.** Todos los números de §6 valen solo con los supuestos del
  prototipo.
- **La resolución de T8 no se estudió aquí.** Con refine_factor = 2, $T_{\max}$ sube a 102,48 °C
  (+0,18). El estudio de malla es de la Fase 6.
- **Memoria.** Se guardan todos los pasos. Con mallas de ~40 000 nodos son ~350 MB. Hará falta
  `save_every` en la Fase 6 si se usan mallas así.
- **`interpolation_matrix` recorre las sondas con un bucle de Python** (vectorizado sobre los
  triángulos). Con pocas sondas está bien; para perfiles de 101 puntos tarda poco (medido
  dentro de los 0,2 s de la corrida), pero no escala a miles de puntos.
- **El `config.yaml` efectivo pierde los comentarios** del YAML original.

## 8. Para revisar y para el informe

1. **El flujo completo de `run`:** config → materiales → malla → $\gamma$ y carga →
   $\mathbf{M}, \mathbf{K}, \mathbf{H}, \mathbf{F}$ → $t_f$ en la malla → esquema $\theta$ →
   sondas → resumen. Cada paso corresponde a una fase anterior y a su prueba.
2. **Por qué el `config.yaml` de resultados es el efectivo** y por qué se copia la librería: la
   regla 5 (una corrida se reproduce solo con su carpeta) está probada.
3. **Cómo leer el balance de energía del resumen.** Entra = almacenada + convección, con $2\pi$.
   Cierra a $10^{-14}$ por construcción; no mide exactitud.
4. **El sesgo de $\theta = 1$ en T8 (§6).** Si alguien pregunta "¿por qué 72,3 y no 72,8?", la
   respuesta es: la referencia se calculó con Euler implícito y $\Delta t = t_f/100$, que integra
   1 % menos de energía.
5. **Para el informe:**
   - `results/ejemplo_frenado/field.png`: campo en 4 instantes, espesor completo.
   - `tmax.png`: $T_{\max}$ y sondas en superficie, a media profundidad y en el plano medio.
   - `depth_profiles.png`: penetración del calor.
   - La tabla de §6, como resultados del caso de referencia y validación T8/T9.
