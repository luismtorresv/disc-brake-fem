# F1 · T1 — Malla: dataclass `Mesh`, generador estructurado, generador Gmsh y `discfem mesh`

Todas las pruebas pasaron (51 passed) sin relajar tolerancias. Para Gmsh, el criterio "la celda
más pequeña toca `brake`" de T1 se cambió por uno equivalente para mallas no estructuradas; se
explica en §4 y §7.

## 1. Qué se pidió

Fase 1 de `PLAN.md` (§5, "Fase 1 — Malla"): el dataclass `Mesh` (contrato de la §2), un
generador estructurado (`mesh/structured.py`) y uno con Gmsh (`mesh/gmsh_backend.py`). Requisitos:

- (a) nodos exactamente en r_a y r_b;
- (b) en z, graduación con la celda más pequeña tocando `brake`;
- (c) `refine_factor` divide todos los tamaños;
- (d) campos `Distance` + `Threshold` de Gmsh y los seis grupos físicos;
- (e) triángulos antihorarios.

Además, el comando `python -m discfem mesh config.yaml` con figura y las pruebas T1 de la §6.

## 2. Qué se hizo

- `src/discfem/geometry.py`: nuevo. `BOUNDARY_GROUPS` (los 6 nombres) y `RectangleGeometry`
  (validación r_i < r_a < r_b < r_o, e > 0; lectura de la sección `geometry` del YAML).
- `src/discfem/mesh/base.py`: nuevo. `MeshOptions` (sección `mesh` del YAML), `Mesh`
  (`areas`, `r_centroid`, `validate`) y utilidades vectorizadas `boundary_edges`, `orient_edges`.
- `src/discfem/mesh/structured.py`: nuevo. Rejilla tensorial graduada y partida en triángulos.
- `src/discfem/mesh/gmsh_backend.py`: nuevo. Port de `legacy/mesh/create_mesh.py` a la API en
  memoria de Gmsh, con 6 grupos físicos y tamaños del YAML.
- `src/discfem/mesh/__init__.py`: nuevo. `build_mesh(geo, opts)` elige el backend.
- `src/discfem/postprocess.py`: nuevo. `plot_mesh` (sección completa + ampliación, grupos en color).
- `src/discfem/cli.py`: nuevo. Subcomandos `run | mesh | study`; solo `mesh` está implementado.
- `src/discfem/__main__.py`: modificado. Delegación a `cli.main`.
- `configs/ejemplo_frenado.yaml`: nuevo. El caso de referencia de la §4, tal cual.
- `tests/test_mesh.py`: nuevo. Pruebas T1 para ambos backends + validaciones + CLI (44 casos).
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

**Por qué importa el peso r.** En la sección axisimétrica, cada punto (r, z) representa un anillo
de radio r. Por eso todas las integrales de la §3.2 llevan el factor r (el 2π se omite). Para la
malla, la prueba natural es que Σ_e A_e · r̄_e = ∫_Ω r dr dz = e (r_o² − r_i²)/2. Aquí A_e es el
área del triángulo y r̄_e el radio de su centroide. La igualdad es exacta porque r es lineal en
cada triángulo. Si un triángulo sobra, falta o se solapa, la suma falla.

**Contrato `Mesh`.** `nodes` (N×2), `tris` (E×3, antihorarios), `region`, `region_names` y
`edges[grupo]` (n×2). Las aristas de borde se guardan orientadas con el dominio a la izquierda.
Así la normal exterior es n = (t_z, −t_r)/|t|, con t = x_b − x_a el vector de la arista. Eso hará
falta para el término de frontera en la Fase 2. `Mesh.validate` comprueba, sin bucles sobre
elementos, que:

- todas las áreas son > 0;
- todos los nodos se usan;
- los seis grupos son disjuntos y su unión es exactamente el borde.

Una arista es de borde si pertenece a un solo triángulo. Para encontrarlas, a cada arista se le
da la clave `min(a,b)·N + max(a,b)` y se cuentan las repeticiones con `np.unique`.

**Tamaño objetivo.** Los dos backends usan la misma ley. Si d es la distancia a la pista:

h(d) = h_min + (h_max − h_min) · min(d / D, 1)

- h_min = `size_track` / `refine_factor`;
- h_max = `size_far` / `refine_factor`;
- D = `grading_distance`, que no se divide porque es una distancia y no un tamaño.

En Gmsh, esto es el campo `Threshold` con DistMin = 0 y DistMax = D sobre un campo `Distance` a
la curva `brake`.

**Generador estructurado.** `graded_spacing` avanza desde la pista con d_{k+1} = d_k + h(d_k)
hasta pasar la longitud L. Después escala todo por L/d_n ≤ 1. Ninguna celda supera su tamaño
objetivo, y las celdas crecen al alejarse de la pista. Se construye así:

- **En r, por tramos:**
  - [r_i, r_a]: graduada desde r_a hacia adentro;
  - [r_a, r_b]: uniforme con ⌈(r_b − r_a)/h_min⌉ celdas;
  - [r_b, r_o]: graduada desde r_b hacia afuera.

  r_a y r_b se asignan explícitamente, así que son nodos exactos (requisito a).
- **En z:** se gradúa desde z = e hacia z = 0, así que la capa más delgada está bajo `brake`
  (requisito b; trampa 1 de la §7).
- **Triángulos:** cada rectángulo (i, j) se parte en (n00, n10, n11) y (n00, n11, n01). Ambos son
  antihorarios por construcción.
- **Grupos:** las aristas salen directamente de los índices de la rejilla.

**Generador Gmsh.** Se crean los 6 vértices de `RectangleGeometry.corners` en sentido
antihorario, 6 líneas (la línea k es el grupo `BOUNDARY_GROUPS[k]`), la superficie, los grupos
físicos y los campos. Se desactiva el tamaño heredado de puntos, bordes y curvatura para que solo
mande el campo, y se usa un solo hilo para que la malla sea reproducible. Después se leen los
nodos y los triángulos de la API y se compactan los índices. Los triángulos horarios se voltean
(vectorizado) y las aristas de cada grupo se reorientan con `orient_edges`. Como (r_a, e) y
(r_b, e) son vértices de la geometría, son nodos exactos.

| Fórmula / requisito (PLAN.md) | Código |
|---|---|
| A = ½[(r2−r1)(z3−z1) − (r3−r1)(z2−z1)] (§3.3) | `Mesh.areas`, `mesh/base.py:93-98` |
| r̄ = (r1+r2+r3)/3 (§3.3) | `Mesh.r_centroid`, `mesh/base.py:100` |
| Σ A r̄ = e(r_o² − r_i²)/2 (T1) | `test_weighted_area_exact`, `tests/test_mesh.py` |
| `refine_factor` divide los tamaños (Fase 1 c) | `MeshOptions.h_min/h_max`, `mesh/base.py:47-54` |
| h(d) lineal entre h_min y h_max | `graded_spacing`, `mesh/structured.py:17-33`; `Threshold`, `mesh/gmsh_backend.py:45-53` |
| Nodos en r_a, r_b (Fase 1 a) | `_r_grid`, `mesh/structured.py:36-48`; `RectangleGeometry.corners`, `geometry.py:60` |
| Celda más fina en z = e (Fase 1 b) | `_z_grid`, `mesh/structured.py:51-56` |
| Triángulos antihorarios (Fase 1 e) | `structured.py:77`; volteo en `gmsh_backend.py:85-90` |
| Cada arista de borde en un solo grupo | `Mesh.validate`, `mesh/base.py:104-128` |
| Normal exterior n = (t_z, −t_r) | `boundary_edges` / `orient_edges`, `mesh/base.py:137-165` |

## 4. Decisiones

- **La misma ley h(d) en ambos backends.** Así `size_track`, `size_far` y `grading_distance`
  significan lo mismo con `structured` y con `gmsh`. Descarté una razón geométrica fija (por
  ejemplo 1,2) porque habría sido un parámetro más fuera del YAML. DistMin = 0 en vez del 1e-3 del
  script previo, porque el YAML de la §4 no tiene ese parámetro.
- **Prueba "celda más pequeña toca `brake`" para Gmsh.** En una malla no estructurada, el
  triángulo de área mínima puede quedar en la segunda fila. Lo medí: con refine_factor = 1 su
  centroide está a 1,44·h_min de la pista y su área es solo un 4 % menor que la del menor
  triángulo que toca `brake`. Exigir que la toque fallaría por azar del mallador, no por un
  error. Usé dos condiciones que sí detectan la graduación invertida:
  - el centroide del triángulo mínimo está a menos de 2·h_min de `brake`;
  - el área media de los triángulos que tocan `brake` es menor que 1/10 de la de los que tocan
    `symmetry` (medido: 0,010).

  El backend estructurado sí cumple el criterio estricto, y además se verifica que los espesores
  de capa en z crecen de forma monótona.
- **Orientación de las aristas de borde fijada en `Mesh`.** El PLAN no la pedía. La añadí porque
  la Fase 2 necesita la normal exterior y es más barato garantizarla aquí que deducirla después.
- **CLI mínima ahora.** `config.py` (dataclasses y validación) es de la Fase 5. Para no
  adelantarla, `cmd_mesh` lee solo `geometry`, `mesh` y `output` con `from_dict` de cada sección.
  La copia de `config.yaml` a la carpeta de resultados ya se hace (regla 5).
- **`plot_mesh` en `postprocess.py`**, que es su lugar en la §2, en vez de crear otro módulo.
- **Error de codificación en Windows.** La consola usa cp1252 y `print("Σ A·r̄")` lanzaba
  `UnicodeEncodeError`. La prueba de la CLI dentro de pytest no lo vio, porque pytest captura en
  UTF-8. Se corrigió reconfigurando stdout/stderr a UTF-8 en `cli.main`. Se añadió una prueba
  en subproceso con `PYTHONIOENCODING=cp1252`: sin la corrección falla y con ella pasa (las dos
  corridas están en §5).

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
.venv\Scripts\python -m discfem mesh configs/ejemplo_frenado.yaml
```

Salida real de pytest y ruff:

```
All checks passed!
...................................................                      [100%]
51 passed in 2.47s
```

Prueba de regresión de codificación, con la corrección quitada temporalmente:

```
FAILED tests/test_mesh.py::test_cli_mesh - AssertionError: Traceback (most recent call last):
1 failed in 1.90s
```

`python -m discfem mesh configs/ejemplo_frenado.yaml` (Gmsh), salida real:

```
Malla 'gmsh' (refine_factor = 1)
  nodos: 929   triángulos: 1618
  área mín/máx: 3.159e-08 / 6.928e-06 m²
  Σ A·r̄ = 8.800000000000e-05 m³ (exacto e(r_o²−r_i²)/2 = 8.800000000000e-05)
  symmetry  : 20 aristas
  r_outer   : 3 aristas
  top_outer : 8 aristas
  brake     : 196 aristas
  top_inner : 8 aristas
  r_inner   : 3 aristas
  figura: results\ejemplo_frenado\mesh.png
```

La misma configuración con `backend: structured`, salida real:

```
Malla 'structured' (refine_factor = 1)
  nodos: 2190   triángulos: 3924
  área mín/máx: 2.673e-08 / 4.117e-06 m²
  Σ A·r̄ = 8.800000000000e-05 m³ (exacto e(r_o²−r_i²)/2 = 8.800000000000e-05)
  symmetry  : 218 aristas
  r_outer   : 9 aristas
  top_outer : 9 aristas
  brake     : 200 aristas
  top_inner : 9 aristas
  r_inner   : 9 aristas
  figura: results\ejemplo_frenado_structured\mesh.png
```

Mediciones con un script exploratorio (fuera de las pruebas), salida real:

```
structured 1 2190 3924 0.01s sumA rel err 2.71e-16 sumAr rel 0.0     ... argmin touches brake True  meanA brake/sym 0.066
structured 2 7361 13824 0.01s sumA rel err 2.71e-16 sumAr rel 1.54e-16 ... argmin touches brake True  meanA brake/sym 0.063
gmsh 1 929 1618 0.14s       sumA rel err 2.71e-16 sumAr rel 1.54e-16 ... argmin touches brake False meanA brake/sym 0.010
gmsh 2 3028 5580 0.14s      sumA rel err 1.36e-16 sumAr rel 0.0     ... argmin touches brake True  meanA brake/sym 0.0071
```

```
refine_factor  nodos  dist(argmin)/h_min
0.5            270    0.26
1              929    1.44
2              3028   0.47
4              11051  0.47
```

## 6. Resultados

| Criterio de la Fase 1 / T1 | structured | gmsh |
|---|---|---|
| Σ A r̄ vs e(r_o² − r_i²)/2 = 8,8e-5 m³ (tol. 1e-12 rel.) | 0 – 1,5e-16 | 0 – 1,5e-16 |
| Σ A vs e(r_o − r_i) = 8e-4 m² (tol. 1e-12 rel.) | 2,7e-16 | 1,4e-16 – 2,7e-16 |
| Cada arista de borde en un solo grupo | sí (exacto) | sí (exacto) |
| `brake` abarca exactamente [r_a, r_b], z = e | sí (igualdad `==`) | sí (igualdad `==`) |
| Triángulos antihorarios | sí | sí |
| Celda más pequeña junto a `brake` | toca `brake` | a 1,44·h_min (< 2·h_min) |
| Celdas con refine_factor 2 / 1 | 3,52 | 3,45 |
| Nodos con la config de referencia | 2190 | 929 |

La §6 no da valores numéricos de referencia para T1 más allá de las tolerancias; todas se
cumplen. Las dos mallas superan la resolución mínima que exige T8 (≈ 700 nodos). La figura
`results/ejemplo_frenado/mesh.png` muestra la malla fina bajo la pista y gruesa en el fondo,
como pide la Fase 1.

## 7. Limitaciones y riesgos

- **Criterio de T1 para Gmsh.** "La celda más pequeña toca `brake`" se cumple de forma estricta
  solo en la malla estructurada. Para Gmsh se usa el criterio de §4. Si quieres el literal,
  habría que restringirlo al backend estructurado.
- **Celdas alargadas en la malla estructurada.** Bajo la pista, Δr = h_min se mantiene en toda la
  profundidad, mientras Δz crece hasta ~2e-3 m: hay celdas con relación de aspecto ~7:1. Son
  triángulos rectángulos (ángulo máximo 90°), así que no debería afectar la convergencia de P1,
  pero gastan nodos: 2190 frente a 929 de Gmsh. No está verificado con la solución del calor;
  se verá en T8 y en la Fase 6.
- **Gmsh no es exactamente simétrico ni uniforme.** `brake` tiene 196 aristas, no 200. Una
  versión distinta de Gmsh puede dar otra malla; la reproducibilidad está probada solo dentro de
  la misma versión (4.15.2).
- **Bucle de Python en `graded_spacing`.** Es un bucle sobre puntos de una rejilla 1D (decenas de
  iteraciones), no sobre elementos, así que no contradice la regla 3.
- **Solo geometría rectangular.** `type: profile` y las regiones múltiples son de la Fase 7.
  `region` es 0 en todos los triángulos.
- **r_i = 0 (cilindro macizo de T6).** La geometría lo admite, pero no se probó. T6 necesitará
  además una malla uniforme con nz = nx/5, que hoy solo se obtiene con size_track = size_far.

## 8. Para revisar y para el informe

1. **Por qué Σ A r̄ es la prueba correcta en axisimetría.** Es ∫ r dA, el volumen del anillo
   dividido por 2π. Una malla con huecos o solapes no la cumple. Comparar con Σ A, que es la
   prueba del caso plano.
2. **La ley h(d) y por qué la celda más fina va en z = e.** Ahí el gradiente térmico es máximo
   durante el frenado. Con la graduación al revés, el prototipo daba un error de ~1,4 % que no
   bajaba al refinar (§7, trampa 1).
3. **Orientación antihoraria y normal exterior n = (t_z, −t_r).** En la Fase 2, el signo de la
   carga en `brake` y la matriz H dependen de esto.
4. **Por qué hacen falta nodos exactos en r_a y r_b.** Si `brake` no coincide con la pista, la
   energía que entra cambia (§7, trampa 6). Lo prueban `test_group_spans_exactly` y
   `test_nodes_exactly_at_track_ends`.
5. **Figura para el informe:** `results/ejemplo_frenado/mesh.png` (Gmsh, con ampliación junto a
   r_a), y la versión estructurada para comparar. La tabla de §6 sirve como tabla de verificación
   de la malla.
