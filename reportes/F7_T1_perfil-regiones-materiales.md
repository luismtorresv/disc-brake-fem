# F7 · T1 — Perfil poligonal, regiones de material, presets por vehículo, materiales propios y convección por grupo

Todas las pruebas pasaron (229 passed). La fase **no está cerrada del todo**: el catálogo de
materiales reales (acero inoxidable para motos, pastilla sinterizada, hierro gris con fuente)
y el preset `motorcycle` quedan pendientes de la fuente que me vas a pasar (Limpert, *Brake
Design and Safety*). No se inventó ninguna propiedad.

Hubo dos ajustes de pruebas, explicados en §4:

- Una prueba de la Fase 5 buscaba un mensaje de error antiguo. Se actualizó el texto; el
  comportamiento sigue siendo el mismo.
- Una prueba nueva tenía un umbral mío sin justificación ("> 100 nodos"). Se reemplazó por la
  afirmación exacta que corresponde.

> **Nota (2026-11-05, F8_T2).** Los valores del perfil en §5 y §6 se calcularon con los
> materiales supuestos del prototipo. Con los materiales actuales del preset `car` (Limpert) el
> perfil da $T_{\max} = 106{,}71$ °C y $T(60\ \mathrm{s}) = 73{,}9$ °C; ver
> `F8_T2_guia-uso-e-insumos-informe.md`. La conclusión cualitativa no cambia.

## 1. Qué se pidió

Fase 7 de `PLAN.md` (§5, "Fase 7 — Generalización", en modo plan primero):

- (a) `geometry.type: profile`, un polígono $(r, z)$ antihorario con grupos de frontera por
  tramos;
- (b) regiones con materiales distintos;
- (c) nada de $k(T)$.

Criterios: el caso rectangular sigue dando exactamente los mismos resultados, y dos regiones
con material idéntico reproducen el caso homogéneo a $10^{-12}$.

Decisiones tuyas en esta sesión:

- regiones definidas como polígonos generales;
- usar materiales reales (hierro fundido gris para autos y acero inoxidable para motos), con un
  catálogo por tipo de vehículo y la posibilidad de que el usuario describa sus propios
  materiales;
- un override hipotético no basta como ejemplo;
- $h$ y $T_\infty$ por grupo de frontera;
- me pasarás la fuente (Limpert) para los valores del catálogo.

## 2. Qué se hizo

- `src/discfem/geometry.py`: reescrito.
  - `ProfileGeometry` (polígono + un grupo por lado) y `RectangleGeometry`, que ahora admite
    `regions`.
  - `geometry_from_dict`.
  - Utilidades de polígonos: área, $\int r\, dA$ exacto, intersección de segmentos, punto en
    polígono y distancia punto–segmento.
  - Validaciones en español.
- `src/discfem/mesh/gmsh_backend.py`: dos caminos.
  - `geo`, sin regiones: es el mismo código de antes generalizado a $n$ lados.
  - OCC + `fragment`, con regiones: clasificación geométrica de triángulos y aristas, y
    control de áreas por región.
- `src/discfem/mesh/structured.py`: rechaza perfiles y regiones con un mensaje claro.
- `src/discfem/materials.py`:
  - la librería admite el formato con secciones (`materials:` + `presets:`) y sigue aceptando
    el formato plano;
  - `load_presets` y `add_custom_materials`.
- `src/discfem/config.py`:
  - geometría general;
  - `materials.regions` debe nombrar exactamente las regiones de la geometría;
  - nuevas claves `materials.preset`, `materials.custom` y `convection.per_group`.
- `src/discfem/solver.py`: resuelve presets y materiales propios, y añade
  `assemble_convection` ($h$ y $T_\infty$ por grupo). El `summary` gana las secciones `geometry`
  y `convection`.
- `src/discfem/postprocess.py`:
  - `plot_mesh` sombrea las regiones;
  - el campo se refleja solo si `symmetry` está en $z = 0$;
  - los perfiles en profundidad omiten los puntos fuera de un perfil no convexo
    (`interpolation_matrix(..., allow_outside=True)`).
- `src/discfem/cli.py`: `mesh` acepta perfiles e imprime $\int r\, dA$ exacto y el área de cada
  región.
- `configs/materiales.yaml`: formato con secciones y preset `car` (mismos materiales, todavía
  no verificados). `motorcycle` está pendiente.
- `configs/ejemplo_perfil.yaml`: nuevo. Disco con cubo interior más grueso, preset `car` y $h$
  mayor en `r_outer`.
- `tests/data/golden_rectangle.json`: nuevo. Referencia del rectángulo capturada **antes** de
  cualquier cambio (commit `c930b37`).
- `tests/test_profile_regions.py`: nuevo (35 casos).
- `tests/test_solver.py`: un mensaje esperado actualizado (§4).
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

**Perfil poligonal.** En la sección axisimétrica, cualquier forma de revolución se describe
con su contorno en el plano $(r, z)$. El YAML da los vértices en sentido antihorario, de modo que
el dominio queda a la izquierda al recorrerlos, igual que en la convención de normales
exteriores de la Fase 1. También da el grupo de cada lado. Las validaciones son:

- polígono simple (sin autointersecciones), antihorario y con $r \ge 0$;
- un grupo por lado, con los nombres de CLAUDE.md;
- `brake` debe ser una cadena contigua de lados horizontales.

La última condición es necesaria porque la carga de la Fase 4,

$$
q(r, t) = \frac{\gamma\, \varphi\, m\, a\, v(t)}{2\pi\, n_d\, n_c}\cdot\frac{g(r)\, r}{I_g},
$$

reparte la potencia sobre $r \in [r_a, r_b]$ en un plano $z$ = cte. $r_a$ y $r_b$ se toman de
los extremos de esa cadena.

La prueba T1 se generaliza con la fórmula exacta del momento del polígono, que es la referencia
de $\sum A\,\bar r$ de la malla:

$$
\int_\Omega r\, dA = \frac{1}{6}\sum_i (r_i + r_{i+1})\,(r_i\, z_{i+1} - r_{i+1}\, z_i)
$$

**Malla sin regiones (kernel `geo`).** Es el código de las fases 1–6 con $n$ lados en vez de 6.
Para el rectángulo, la secuencia de llamadas a Gmsh es idéntica, por eso la malla resulta
idéntica bit a bit, y eso se verifica contra la referencia guardada antes de los cambios.

**Malla con regiones (kernel OCC).** Cada región es un polígono dentro de la sección. La
operación booleana `fragment` parte la sección en superficies que comparten las aristas de las
interfaces. Así, la malla es **conforme**: ningún triángulo cruza una interfaz entre
materiales, y $k$ y $\rho c$ son constantes en cada triángulo, como supone el ensamblaje de la
Fase 2.

Después de mallar, en vez de seguir los números de entidad que Gmsh reasigna al fragmentar,
se clasifica por geometría:

- cada triángulo va a la región que contiene su centroide (el centroide nunca está sobre una
  interfaz);
- cada arista de borde va al lado del perfil sobre el que están sus dos extremos.

El control es comparar el área mallada de cada región con la del polígono a $10^{-9}$. Si una
región se solapa con otra o sale del dominio, falla con un mensaje claro.

**Materiales.**

- `materials.preset: car` toma el disco y la pastilla de la sección `presets` de la librería.
  Lo escrito explícitamente en `regions` o `pad` tiene prioridad.
- `materials.custom` permite definir materiales en el YAML de la corrida, con las mismas
  validaciones que la librería (`source` y `verified` obligatorios, rangos de sanidad).
- $\gamma$ automático usa el material bajo `brake`. Con una capa bajo la pista, es el de la capa.

**Convección por grupo.** La condición de Robin
$-k\, \nabla T \cdot \mathbf{n} = h_g\,(T - T_{\infty,g})$ puede tener $h$ y $T_\infty$ distintos
en cada grupo $g$:

$$
\mathbf{H} = \sum_g \mathbf{H}_g(h_g), \qquad
\mathbf{F}_{\text{conv}} = \sum_g \mathbf{H}_g\,(T_{\infty,g}\, \mathbf{1}),
$$

usando la identidad $\int h\, T_\infty\, \varphi_i\, r\, ds = (\mathbf{H}\, T_\infty \mathbf{1})_i$
de la Fase 2. Si todos los grupos comparten $h$ y $T_\infty$, se ensambla en una sola llamada,
como antes, para no cambiar el orden de las sumas en punto flotante.

| Elemento | Código |
|---|---|
| $\int r\, dA$ exacto del polígono | `polygon_r_moment`, `geometry.py:50` |
| Validación del polígono (simple, antihorario, $r \ge 0$) | `_validate_polygon`, `geometry.py:103` |
| Perfil, $r_a$ y $r_b$ desde `brake` | `ProfileGeometry`, `geometry.py:276`; `r_a`/`r_b`, `:329-334` |
| Rectángulo con regiones | `RectangleGeometry`, `geometry.py:192` |
| Malla sin regiones (idéntica para el rectángulo) | `_mesh_geo_kernel`, `gmsh_backend.py:137` |
| Malla con regiones (OCC + fragment) | `_mesh_occ_kernel`, `gmsh_backend.py:172` |
| Clasificación de aristas y regiones | `_classify_boundary` (`:97`), `_classify_regions` (`:118`) |
| Presets y materiales propios | `load_presets` (`materials.py:170`), `add_custom_materials` (`:191`); `solver.py:94` |
| $\mathbf{H}$ y $\mathbf{F}_{\text{conv}}$ por grupo | `assemble_convection`, `solver.py:156`; `per_group`, `config.py:198` |
| Regiones de material = regiones de la geometría | `config.py:166` |

## 4. Decisiones

- **Dos kernels de Gmsh.** Con un solo camino OCC, la malla del rectángulo habría cambiado, y
  con ella todos los resultados validados (T8, T9, estudio de malla). El criterio del PLAN
  ("exactamente los mismos resultados") obliga a mantener el camino `geo` para el caso sin
  regiones.
- **Clasificación geométrica después de mallar**, en lugar de grupos físicos, en el camino OCC.
  Es independiente de cómo `fragment` renumera las entidades, y funciona igual para
  rectángulos y perfiles.
- **`brake` horizontal y contigua.** Una pista inclinada o partida necesitaría otra
  normalización de $q(r)$. Fuera de alcance.
- **La malla estructurada no se generaliza.** Una rejilla tensorial no sigue polígonos
  arbitrarios. El error dice qué hacer (usar `gmsh`).
- **Materiales: no invento ninguno.** Busqué fuentes primarias (§7). Las únicas legibles eran
  insuficientes o de otro grado de acero. Por eso el catálogo espera la fuente que me vas a
  pasar. El preset `car` usa los materiales del prototipo (todavía no verificados) y
  `motorcycle` queda comentado como pendiente.
- **Pruebas ajustadas (no se relajó ninguna tolerancia):**
  - `test_solver.py` esperaba el mensaje "una sola región". La configuración inválida se sigue
    rechazando, pero ahora con el mensaje general "regiones que la geometría no tiene".
  - En `test_regions_areas_and_conformity` había escrito "> 100 nodos en la interfaz", un
    número sin justificar. Con $h \approx 6{,}7\times10^{-4}$ m a 1 mm de la pista, 0,06 m de
    interfaz dan ~90 nodos (hay 91). Lo reemplacé por la afirmación exacta: las aristas de la
    malla sobre $z = 0{,}009$ suman 0,06 m a $10^{-12}$, es decir, la interfaz está completamente
    mallada.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
.venv\Scripts\python -m discfem run configs/ejemplo_perfil.yaml
```

Salida real:

```
All checks passed!
229 passed in 17.76s
```

Comparación del rectángulo contra la referencia guardada antes de los cambios. Salida real del
script de comprobación; lo mismo verifica `test_rectangle_unchanged`:

```
gmsh nodes True tris True edges True T_final True summary True
structured nodes True tris True edges True T_final True summary True
```

Exploración de los dos casos nuevos (salida real):

```
perfil 932 0.09s r_a,r_b 0.08 0.14 sumAr err 1.5029140178618032e-16 {'symmetry': 20, 'r_outer': 3, 'top_outer': 8, 'brake': 196, 'top_inner': 11, 'r_inner': 4} ['disc'] [1620]
capa 998 0.19s r_a,r_b 0.08 0.14 sumAr err 1.5400599040987277e-16 {'symmetry': 20, 'r_outer': 3, 'top_outer': 8, 'brake': 196, 'top_inner': 8, 'r_inner': 3} ['disc', 'capa'] [756, 1000]
capa=gci: Tmax 102.26864462915952 gamma 0.8692006698254295 bal 2.6225331786492112e-15
perfil: Tmax 102.29683315674657 T60 70.62176631708277 bal 5.58599567052282e-14
```

`python -m discfem run configs/ejemplo_perfil.yaml` (salida real):

```
ADVERTENCIA (UnverifiedMaterialWarning): Materiales NO verificados: gray_cast_iron (source: PENDIENTE: citar fuente); organic_pad (source: PENDIENTE: citar fuente). Los resultados no sirven para el informe hasta citar fuentes y poner verified: true.
Caso 'ejemplo_perfil_cubo': 3047 nodos, 5610 triángulos, 1080 pasos (θ = 0.5, Rannacher = 2, Δt = 0.05556 s)
  γ = 0.8692   t_f = 5.556 s
  T_max = 102.925 en t = 3.389 s (r = 0.1092, z = 0.0100)
  T_max(t_end = 60 s) = 70.961
  energía [J, con 2π]: entra 88023.2, sale por convección 3413.3, almacenada 84609.9; cierre del balance 2.2e-14
  sonda (r=0.11, z=0.01): T_max = 102.924 en t = 3.389 s; T_final = 70.928
  sonda (r=0.11, z=0): T_max = 87.501 en t = 7.001 s; T_final = 70.929
  sonda (r=0.0725, z=0.016): T_max = 68.675 en t = 52.393 s; T_final = 68.597
  Materiales NO verificados: gray_cast_iron, organic_pad
  tiempo de cómputo: 0.51 s   salidas: results\ejemplo_perfil
```

## 6. Resultados

| Criterio | Resultado |
|---|---|
| Rectángulo idéntico (malla, aristas, $T$ final, summary), gmsh y estructurada | idéntico bit a bit |
| Perfil con las 6 esquinas del rectángulo = rectángulo | misma malla y mismas aristas |
| Perfil con cubo: $\sum A\,\bar r$ frente a $\int r\, dA$ exacto | $1{,}5\times10^{-16}$ |
| Perfil: $2\pi\sum_i F_{\text{brake},i}$ frente a $\gamma P_f$ | $< 10^{-12}$ |
| Capa: área mallada de la región frente al polígono | $< 10^{-12}$ |
| Dos regiones con el mismo material frente a una sola, misma malla (criterio del PLAN) | $< 10^{-12}$ relativo al aumento de $T$ |
| Capa con el mismo material frente a la malla sin capa | 102,269 frente a 102,297 °C (0,03 %, solo discretización) |
| Capa menos conductora ($k = 25$, material de prueba) | $T_{\max}$ sube más de 1 °C; $\gamma$ baja (usa el material de la capa) |
| `per_group` igual al global | resultados idénticos |

Lectura física del ejemplo de perfil:

- El cubo más grueso añade masa térmica lejos de la pista. $T_{\max}$ casi no cambia (102,93 °C
  con CN + Rannacher, igual que el rectángulo con el mismo esquema), porque el pico es
  superficial y ocurre antes de que el calor llegue al cubo.
- La temperatura final baja de ~72,8 a 70,96 °C: misma energía repartida en más volumen, y más
  convección por $h = 60$ en `r_outer`.
- La sonda del cubo llega a su máximo a los 52 s: así de lento es el transporte radial por
  conducción.

## 7. Limitaciones y riesgos

- **Catálogo de materiales reales PENDIENTE.** Lo que encontré y por qué no lo usé:
  - Hierro gris: Jogineedi et al., SAE 2021-01-5106, Tabla 4. Son mediciones de rotores ASTM A48
    a 25 °C ($k$ = 68–78 $\mathrm{W/(m\,K)}$, $c_p$ = 0,63–0,66 $\mathrm{J/(g\,K)}$), muy
    distintas de los valores de manual, y sin la densidad en el texto extraíble.
  - Acero: la ficha de Cleveland-Cliffs es del **410S** (ferrítico, no templable), no del
    410/420 martensítico de los discos de moto.
  - Pastillas: los artículos FEM consultados difieren un orden de magnitud entre sí.

  Espero la fuente (Limpert) para cargar valores con tabla y página en `source`.
- **Superficie de fricción horizontal.** Una pista cónica o inclinada no está soportada.
- **La malla con regiones no reproduce bit a bit la malla sin regiones.** Es inevitable: las
  interfaces añaden aristas. La comparación válida ($10^{-12}$) es sobre la misma malla.
- **La clasificación de aristas usa una tolerancia geométrica** ($10^{-9}$ × tamaño). Perfiles
  con lados de menos de ~$10^{-10}$ m no están contemplados.
- **La figura del campo en el perfil** deja espacio en blanco a la izquierda (aspecto igual con
  paneles altos). Es cosmético.

## 8. Para revisar y para el informe

1. **Cómo un polígono $(r, z)$ describe un sólido de revolución** y por qué $\int r\, dA$ (no el
   área) es la prueba de que la malla está completa.
2. **Qué significa una malla conforme en las interfaces** y por qué hace falta. $k$ y $\rho c$
   son constantes por triángulo; un triángulo que cruzara la interfaz promediaría dos materiales.
3. **Por qué el rectángulo debía dar el mismo resultado bit a bit** (no romper lo validado) y
   cómo se logró: mismo kernel y misma secuencia de llamadas.
4. **El criterio "dos regiones iguales = homogéneo a $10^{-12}$"** prueba el ensamblaje por
   regiones sobre la misma malla. Contra otra malla solo cabe esperar el error de discretización
   (0,03 %).
5. **Para el informe:**
   - `results/ejemplo_perfil/field.png` y `mesh.png`: ejemplo de perfil no rectangular.
   - La malla de la capa (`mesh.png` de una corrida con `geometry.regions`) como ilustración
     de las regiones.
   - Cuando haya materiales con fuente: una comparación auto (hierro gris) frente a moto
     (acero inoxidable).
