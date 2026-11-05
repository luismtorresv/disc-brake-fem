# Plan de acción — solver FEM axisimétrico del disco de freno

Este plan está pensado para ejecutarse con Claude Code, **una fase por sesión**. Cada fase tiene un
objetivo, un prompt para pegar tal cual, y criterios de aceptación verificables. La física y las
fórmulas exactas están en la sección 3; las pruebas con sus tolerancias, en la sección 6.

Las tolerancias y los valores de referencia de la sección 6 salen de un prototipo mío de ~100
líneas (numpy/scipy) que implementa justo las fórmulas de la sección 3. No es el código del
proyecto; solo sirve para que las pruebas que se le piden a Claude Code sean alcanzables.

---

## 1. Alcance y decisiones tomadas

**Lo que el programa hace:** resuelve la ecuación de calor transitoria en la sección (r, z) del
disco con calentamiento por fricción en `brake`, convección en las superficies expuestas y
simetría en z = 0. Permite cambiar por archivo: geometría (dimensiones y, después, perfil
poligonal), resolución de malla, materiales y propiedades, parámetros de frenado, condiciones de
frontera, esquema temporal.

**Lo que NO hace (queda como trabajo futuro en el informe):** discos ventilados/perforados (3D),
fuente móvil, radiación, termomecánica. La versión extensa del anteproyecto los describe con 3D;
el programa no los cubre y el informe debe decirlo explícitamente.

**Decisiones que asumí y que pueden cambiar:**
- YAML como formato de entrada; librería de materiales en un YAML aparte.
- Malla con Gmsh (API de Python, sin pasar por archivo) y un generador estructurado en numpy de
  respaldo, útil para pruebas y para correr sin Gmsh.
- Elementos P1 triangulares, solver con `scipy.sparse` y factorización LU reutilizada.
- Una sola frenada seguida de enfriamiento (t_end ≥ t_f). Frenadas repetidas quedan opcionales.

### Inconsistencias que encontré entre las dos versiones del anteproyecto

Las fórmulas son idénticas en ambas. Antes de implementar, conviene corregir esto en el .tex:

1. Versión corta, resumen: dice "condición isotérmica/aislada en el plano de simetría". Es solo
   aislada (Γ_N, flujo nulo); el programa no implementa Dirichlet en z = 0.
2. Versión corta, resultados esperados: "evolución tridimensional". El modelo es 2D axisimétrico
   en (r, z, t); la frase confunde con 3D.
3. Versión corta, matrices del sistema: aparece `dΩ_0`, que no está definido. Debe ser `dr dz`.
4. Ambas versiones: la tabla de mallas manda ventilado y perforado a "Modelo 3D". Ese caso no lo
   cubre este programa.
5. La versión corta dice que el código de la malla está en
   `github.com/luismtorresv/disc-brake-fem`. No pude ver ese repositorio, así que no sé qué hay
   ya hecho; la Fase 0 empieza por inspeccionarlo.

---

## 2. Arquitectura y estructura de carpetas

```
disc-brake-fem/
├─ CLAUDE.md                  contexto y reglas para Claude Code
├─ PLAN.md                    este archivo (con sección ESTADO al final)
├─ pyproject.toml
├─ README.md
├─ configs/
│   ├─ ejemplo_frenado.yaml        caso de referencia (sección 4)
│   ├─ validacion_semiinfinito.yaml
│   ├─ estudio_malla.yaml
│   └─ materiales.yaml             librería de materiales
├─ src/discfem/
│   ├─ config.py        dataclasses + carga y validación del YAML
│   ├─ materials.py     librería, efusividad, γ
│   ├─ geometry.py      definición de la sección y de los grupos de frontera
│   ├─ mesh/
│   │   ├─ base.py          dataclass Mesh
│   │   ├─ structured.py    generador numpy
│   │   └─ gmsh_backend.py  generador con campos de tamaño
│   ├─ fem/assembly.py  M, K, H, F vectorizados
│   ├─ loads.py         q_fric(r,t): frenado | flujo constante | tabla
│   ├─ timestepping.py  esquema θ, Rannacher, malla temporal con t_f
│   ├─ solver.py        orquesta todo; devuelve un objeto Result
│   ├─ postprocess.py   Tmax(t), sondas, balance de energía, figuras, VTK/CSV
│   ├─ convergence.py   estudio de malla y de Δt
│   └─ cli.py           python -m discfem run | mesh | study
├─ tests/
└─ results/                   ignorada por git
```

Dataclass `Mesh` (contrato entre malla y solver):

```
nodes:   (N, 2) float     coordenadas (r, z)
tris:    (E, 3) int       conectividad, antihoraria (área > 0)
region:  (E,)   int       índice de región (material)
region_names: list[str]
edges:   dict[str, (n, 2) int]   aristas por grupo de frontera
```

---

## 3. Especificación física y numérica (fórmulas exactas)

### 3.1 Problema

```
ρ c_p ∂T/∂t − k [ (1/r) ∂/∂r ( r ∂T/∂r ) + ∂²T/∂z² ] = 0      en Ω × (0, t_end]
−k ∇T·n = −q_fric(r,t)        en brake
−k ∇T·n = h (T − T∞)          en grupos con convección
−k ∇T·n = 0                   en grupos aislados (symmetry por defecto)
T(r,z,0) = T0
```

Con propiedades por región: reemplazar k, ρ c_p por k(región), (ρ c_p)(región), constantes en cada
elemento.

### 3.2 Sistema semidiscreto

```
M Ṫ + (K + H) T = F(t)
M_ij = ∫ ρc φi φj r dr dz      K_ij = ∫ k ∇φi·∇φj r dr dz      H_ij = ∫_conv h φi φj r ds
F_i  = ∫_brake q_fric φi r ds + ∫_conv h T∞ φi r ds
```

Modo plano (solo validación): cambiar el peso `r` por `1` en todas las integrales.

### 3.3 Matrices elementales P1 exactas (sin cuadratura)

Triángulo con nodos 1, 2, 3 de coordenadas (r_i, z_i), antihorario:

```
b1 = z2 − z3   b2 = z3 − z1   b3 = z1 − z2
c1 = r3 − r2   c2 = r1 − r3   c3 = r2 − r1
A  = ½ [ (r2−r1)(z3−z1) − (r3−r1)(z2−z1) ]  > 0        ∇φi = (bi, ci) / (2A)
r̄  = (r1 + r2 + r3) / 3

K^e_ij = k · r̄ · (bi bj + ci cj) / (4A)                 (exacto: el gradiente es constante)
M^e_ii = ρc · A · [ r_i/10 + (r_j + r_l)/30 ]            (j, l = los otros dos nodos)
M^e_ij = ρc · A · [ (r_i + r_j)/30 + r_l/60 ]            (i ≠ j; l = el tercer nodo)
```

Para el modo plano, sustituir cada `r` por 1.

Arista (a, b) de longitud L, con radios r_a, r_b:

```
H^e = h · L · [ [ r_a/4 + r_b/12 ,  (r_a + r_b)/12 ],
                [ (r_a + r_b)/12 ,  r_a/12 + r_b/4 ] ]
F^e_i = ∫ q(r) r φi ds   →  Gauss–Legendre de 3 puntos (exacto hasta grado 5)
```

Con 3 puntos, el integrando de presión uniforme (q ∝ r, peso r, φ lineal) queda exacto.

### 3.4 Flujo de fricción (modelo `braking`)

```
v(t) = v0 − a t         t_f = v0 / a        q_fric(r,t) = 0 para t > t_f
P_f(t) = φ m a v(t) / (n_d n_c)             potencia por UNA cara de un disco
I_g = ∫_{r_a}^{r_b} g(r) r² dr              g = 1 (presión uniforme) | g = 1/r (desgaste uniforme)
q_fric(r,t) = γ φ m a v(t) / (2π n_d n_c) · g(r) r / I_g
γ = √(kρc)_disco / ( √(kρc)_disco + √(kρc)_pastilla )     (efusividades)
```

Casos cerrados, útiles como prueba:

```
presión uniforme:  q = 3 γ φ m a v(t) r / [ 2π n_d n_c (r_b³ − r_a³) ]
desgaste uniforme: q =   γ φ m a v(t)   / [  π n_d n_c (r_b² − r_a²) ]     (constante en r)
consistencia:      2π ∫ q_fric r dr = γ P_f(t)
```

Parámetros: `m` masa, `v0` velocidad inicial, `a` desaceleración, `R_w` radio de rueda (solo
aparece si se pide ω o la presión), `φ` fracción de frenado del eje, `n_d` discos por eje,
`n_c = 2` caras por disco, `γ` partición, `r_a, r_b` pista, `g(r)` perfil de presión.
`μ` y el ángulo de la pastilla θ_p se cancelan con desaceleración prescrita; no son entradas.
Si algún día se prescribe la presión en vez de la desaceleración, `μ` sí entra.

### 3.5 Esquema temporal

```
[M + θΔt(K+H)] Tⁿ⁺¹ = [M − (1−θ)Δt(K+H)] Tⁿ + Δt [ θ Fⁿ⁺¹ + (1−θ) Fⁿ ]
```

- θ = 1 Euler implícito; θ = 1/2 Crank–Nicolson.
- Δt constante: una sola factorización LU (`scipy.sparse.linalg.splu`) reutilizada. Con
  Rannacher (n pasos de θ = 1 y luego θ = 1/2) son dos factorizaciones.
- Δt por defecto: t_f / 100, ajustado para que t_f sea múltiplo entero de Δt.
- El flujo arranca de golpe en t = 0; Crank–Nicolson puede oscilar al inicio. En mi prototipo
  dio errores de hasta ~2 % en los primeros instantes frente a ~0,5 % con θ = 1.

### 3.6 Balance de energía (control de calidad)

Como `1ᵀK = 0`, sumando el sistema se obtiene `d/dt (1ᵀ M T) = 1ᵀ F_brake − Q_conv`, con
`Q_conv = 1ᵀ(H T − F_conv)`. Con el esquema θ, el cierre es exacto a precisión de máquina si se
integra con los mismos pesos θ. Todo se mide **sin** el factor 2π (ver CLAUDE.md).

---

## 4. Formato de entrada

`configs/ejemplo_frenado.yaml` (el caso de referencia; los valores de materiales y vehículo son
**supuestos míos, no verificados**):

```yaml
case: ejemplo_frenado_100kmh

geometry:
  type: rectangle            # rectangle | profile (Fase 7)
  r_inner: 0.07              # m
  r_outer: 0.15
  half_thickness: 0.01       # semiespesor e
  track: {r_a: 0.08, r_b: 0.14}
  axisymmetric: true         # false => peso 1 (solo validación)

mesh:
  backend: gmsh              # gmsh | structured
  refine_factor: 1.0         # divide TODOS los tamaños; 2 => malla ~4x más celdas
  size_track: 3.0e-4         # tamaño junto a brake [m]
  size_far: 4.0e-3           # tamaño lejos [m]
  grading_distance: 0.01     # distancia de transición [m]

materials:
  library: configs/materiales.yaml
  regions: {disc: gray_cast_iron}
  pad: organic_pad           # solo para γ
  overrides: {}              # ej.: {disc: {k: 48.0}}

boundaries:                  # condición por grupo
  brake: flux
  top_inner: convection
  top_outer: convection
  r_inner: convection
  r_outer: convection
  symmetry: insulated

convection: {h: 40.0, T_inf: 25.0}
initial: {T0: 25.0}

load:
  type: braking              # braking | constant_flux | table
  vehicle:
    mass: 1500.0
    v0: 27.78                # 100 km/h
    decel: 5.0
    wheel_radius: 0.30
    brake_fraction: 0.7
    discs_per_axle: 2
    faces_per_disc: 2
  pressure_profile: uniform_wear   # uniform_pressure | uniform_wear
  gamma: auto                      # auto | número en (0,1)

time:
  theta: 1.0
  rannacher_steps: 0
  dt: auto                   # t_f/100
  t_end: 60.0

output:
  dir: results/ejemplo_frenado
  probes: [[0.11, 0.01], [0.11, 0.005], [0.11, 0.0]]   # (r, z)
  formats: [csv, vtk, png]
```

`configs/materiales.yaml`: cada entrada **obliga** a `source` y `verified`:

```yaml
gray_cast_iron:
  description: Hierro fundido gris (disco)
  k: 50.0          # W/(m·K)   SUPUESTO DEL PROTOTIPO
  rho: 7200.0      # kg/m³     SUPUESTO DEL PROTOTIPO
  cp: 460.0        # J/(kg·K)  SUPUESTO DEL PROTOTIPO
  source: "PENDIENTE: citar fuente"
  verified: false
organic_pad:
  description: Pastilla orgánica (solo efusividad para γ)
  k: 1.5
  rho: 2500.0
  cp: 1000.0
  source: "PENDIENTE: citar fuente"
  verified: false
```

Rangos de sanidad (aproximados, de memoria; verificar): advertir, sin abortar, si k ∉ (0,01; 500)
W/(m·K), ρ ∉ (100; 20 000) kg/m³, c_p ∉ (100; 5000) J/(kg·K), h ∉ (1; 500) W/(m²·K).

---

## 5. Fases

### Fase 0 — Inspección y arranque (modo plan; sin escribir código)

**Objetivo:** saber qué existe ya y preparar el proyecto.

**Prompt:**
> Lee `CLAUDE.md` y `PLAN.md`. Inspecciona este repositorio: qué archivos hay, si ya existe código
> de malla o de FEM, qué dependencias usa y si corre. No modifiques nada todavía. Entrégame (1)
> un resumen de lo reutilizable, (2) los conflictos con la estructura de la sección 2 y (3) un
> plan de migración mínimo. Si el repositorio está vacío, di eso y propón el `pyproject.toml`.
> Espera mi confirmación antes de crear archivos.

**Criterios:** resumen entregado; decisión tomada entre reutilizar o empezar de cero; `pytest`
instalado y corriendo (aunque sea con una prueba trivial); `git` inicializado.

### Fase 1 — Malla

**Objetivo:** `Mesh` correcta, con grupos de frontera, refinable por parámetro.

**Prompt:**
> Implementa la Fase 1 de `PLAN.md`: el dataclass `Mesh` (sección 2), el generador estructurado
> (`mesh/structured.py`) y el de Gmsh (`mesh/gmsh_backend.py`). Requisitos: (a) la malla debe
> tener nodos exactamente en r_a y r_b; en el generador estructurado construye la rejilla en r
> por tramos [r_i, r_a], [r_a, r_b], [r_b, r_o]; (b) en z, graduación geométrica con la celda
> MÁS pequeña tocando `brake` (z = e); verifícalo con una prueba, porque es un error fácil de
> cometer; (c) `refine_factor` divide todos los tamaños; (d) en Gmsh usa un campo `Distance` +
> `Threshold` sobre la curva `brake`, y los grupos físicos `brake`, `top_inner`, `top_outer`,
> `r_inner`, `r_outer`, `symmetry`; (e) todos los triángulos con orientación antihoraria. Añade
> `python -m discfem mesh config.yaml` que dibuje la malla coloreando cada grupo de frontera.
> Escribe las pruebas T1 de la sección 6 y hazlas pasar.

**Criterios:** T1 pasa con ambos backends; la figura muestra la malla fina bajo la pista y gruesa
al fondo; ambos backends dan áreas totales iguales a e·(r_o − r_i) a 1e-12.

### Fase 2 — Núcleo FEM

**Objetivo:** M, K, H, F vectorizados y verificados.

**Prompt:**
> Implementa `fem/assembly.py` según la sección 3.3 de `PLAN.md`: matrices elementales cerradas
> para M y K con peso r (y variante con peso 1), matriz H y vector F por aristas con Gauss de 3
> puntos. Todo vectorizado con numpy y ensamblado con `scipy.sparse` (COO → CSR); sin bucles
> sobre elementos. Propiedades por región (arreglos por elemento). Escribe las pruebas T2 y T6 y
> hazlas pasar. En particular T2 compara las matrices elementales con una cuadratura de alto
> orden sobre triángulos aleatorios: es la prueba que protege contra errores en los coeficientes
> 1/10, 1/30, 1/60.

**Criterios:** T2 y T6 pasan. T6 es la que valida M y K con peso r; la prueba de energía por sí
sola solo valida sumas de entradas.

### Fase 3 — Integración temporal y validación analítica

**Objetivo:** esquema θ con LU reutilizada, Rannacher y t_f exacto; validar contra Duhamel.

**Prompt:**
> Implementa `timestepping.py` (sección 3.5): esquema θ con factorización reutilizada, opción de
> n pasos iniciales de Euler implícito (Rannacher), Δt automático = t_f/100 ajustado para que t_f
> sea múltiplo de Δt, y balance de energía según 3.6. Escribe T4, T5 y T7. En T7, compara θ = 1,
> θ = 1/2 y θ = 1/2 con 2 pasos Rannacher y reporta en una tabla el error máximo relativo para
> t ≥ 0,2 s. Propón el θ por defecto con base en esos números y justifícalo en un comentario del
> código, no en tu opinión.

**Criterios:** T4, T5 y T7 pasan; hay una tabla de errores por esquema.

### Fase 4 — Cargas y materiales

**Prompt:**
> Implementa `loads.py` (modelo `braking` de la sección 3.4, `constant_flux` y `table` con
> interpolación lineal) y `materials.py` (librería con `source` y `verified` obligatorios,
> efusividad, γ automático desde el material de la región adyacente a `brake`, advertencias de
> rango). Escribe T3. Si falta `source` en un material, el programa debe fallar con un mensaje en
> español; si `verified: false`, debe advertir al inicio de cada corrida.

**Criterios:** T3 pasa; γ del caso de referencia ≈ 0,87 con los materiales supuestos.

### Fase 5 — Configuración, CLI y salidas

**Prompt:**
> Implementa `config.py` (dataclasses + validación con mensajes de error claros en español),
> `solver.py` (`run(config, overrides=None) -> Result`), `postprocess.py` y `cli.py`. Salidas por
> corrida: `config.yaml` copiado, `tmax.csv` (t, T_max, nodo), `probes.csv`, `summary.json` (T_max,
> instante, balance de energía, nodos, elementos, tiempo de cómputo, materiales no verificados),
> figuras PNG (campo T(r,z) en 4 instantes reflejando z para ver el espesor completo, T_max(t),
> perfiles de T en profundidad) y VTK. `overrides` acepta claves con puntos
> (`"mesh.refine_factor": 2`) para barridos desde Python. Escribe T8 y T9.

**Criterios:** T8 y T9 pasan; `python -m discfem run configs/ejemplo_frenado.yaml` produce todas las
salidas; correr dos veces da el mismo `summary.json`.

### Fase 6 — Estudios de convergencia y paramétricos

**Prompt:**
> Implementa `convergence.py` y `python -m discfem study`. Estudio de malla: niveles con
> `refine_factor` = 1, 2, 4, 8, mismo Δt; tabla con nodos, elementos, T_max, instante de T_max,
> ε_malla = |Tmax(i+1) − Tmax(i)| / Tmax(i+1) × 100 %, orden observado p = log2(|T_h − T_h/2| /
> |T_h/2 − T_h/4|) y tiempo de cómputo; selecciona el primer nivel con ε_malla < 1 % también en
> las sondas. Estudio de Δt aparte (Δt, Δt/2, Δt/4). Barridos sobre a, v0, m y materiales con
> tabla y figura comparativa. Las figuras salen listas para el informe (PDF vectorial).

**Criterios:** el estudio de malla del caso de referencia cierra por debajo de 1 % desde el
segundo nivel (en mi prototipo bajó de 0,12 % a 0,01 %); los barridos reproducen la tendencia
esperada (más desaceleración o más masa → más T_max).

### Fase 7 — Generalización (opcional)

**Prompt (modo plan primero):**
> Propón e implementa, sin romper lo anterior: (a) `geometry.type: profile` con un polígono (r, z)
> antihorario y asignación de grupos de frontera por tramos; (b) regiones con materiales
> distintos (por ejemplo una capa fina bajo la pista con otro material); (c) nada de k(T) por
> ahora. Dime primero qué cambios requiere el generador de Gmsh y qué pruebas añadirías.

**Criterios:** el caso rectangular sigue dando exactamente los mismos resultados; hay una prueba
con dos regiones y un material idéntico en ambas que reproduce el caso homogéneo a 1e-12.

### Fase 8 — Documentación y reproducibilidad

**Prompt:**
> Escribe `README.md` (instalación en Windows y Linux, ejemplos de cada comando, cómo editar un
> YAML, cómo añadir un material, cómo leer las salidas) y congela las versiones con
> `pip freeze` en `requirements-lock.txt`. Verifica en un entorno limpio que
> `pip install -e ".[dev]"` y `pytest -q` funcionan.

---

## 6. Batería de pruebas

Los valores "medido" son de mi prototipo; úsalos como referencia de que la tolerancia es
alcanzable, no como valores a copiar sin entender.

| ID | Qué verifica | Criterio |
|---|---|---|
| T1 | Malla: Σ A r̄ = e (r_o² − r_i²)/2; cada arista de borde en un solo grupo; `brake` abarca exactamente [r_a, r_b]; triángulos antihorarios; la celda más pequeña toca `brake` | Σ A r̄ a 1e-12 relativo; el resto exacto |
| T2 | Matrices elementales vs cuadratura de alto orden (triángulos aleatorios); K·1 = 0; M simétrica y definida positiva; 1ᵀM1 = ρc ∫ r dA | 1e-12 relativo |
| T3 | q_fric: 2π ∫ q r dr = γ P_f(t) para ambos perfiles; coincide con las fórmulas cerradas de 3.4 | 1e-12 relativo |
| T4 | Balance de energía con convección apagada y encendida (3.6) | cierre < 1e-10 relativo (medido ~4e-15) |
| T5 | T0 = T∞ y sin flujo: T queda uniforme | \|T − T0\| < 1e-10 |
| T6 | Modo propio axisimétrico: cilindro r ∈ [0, R] con R = 0,05, semiespesor 0,01, k = 50, ρ = 7200, c_p = 460, aislado por todos lados, malla uniforme con nz = nx/5, θ = 1/2, Δt = 0,02, t = 10 s. T0 = 1 + 0,5 J0(λ r/R) con λ = primer cero de J1 (3,8317); solución exacta 1 + 0,5 J0(λ r/R) exp(−α λ² t / R²) | error máx./amplitud inicial < 1 % con nx = 20 (medido 0,31 %) y decreciente al refinar nx = 10, 20, 40, 80 (medido 0,86, 0,31, 0,10, 0,03 %) |
| T7 | Semi-infinito, modo plano: ancho 0,01 con 2 celdas en r, e = 0,06 con 40 capas en z, aislado salvo el flujo en la cara superior; k = 50, ρ = 7200, c_p = 460, q(t) = q0 (1 − t/t_f) con q0 = 7,7e5 y t_f = 5,56; primera celda 1e-4 a 2e-4 junto a la cara calentada; Δt = 0,005 a 0,01; T_s exacta = 2 q0 /√(π kρc) · [√t − 2 t^{3/2}/(3 t_f)] en el nodo superior | θ = 1: error máx. < 1 % para t ≥ 0,2 s (medido 0,5 %) y < 0,2 % en t = 3 s (medido 0,007–0,06 %). θ = 1/2: < 3 % (medido 1–2 %) |
| T8 | Caso de referencia (`ejemplo_frenado.yaml`), desgaste uniforme, Δt = t_f/100, resolución equivalente a una malla estructurada de 40 × 16 celdas (≈ 700 nodos) o más fina | T_max = 102,5 ± 0,5 °C en t = 3,4 ± 0,1 s; T_max(60 s) = 72,3 ± 0,5 °C |
| T9 | Cota energética del caso de referencia: ΔT_final ≤ γ E_cara / (ρ c π (r_o² − r_i²) e) | cota 48,1 K; medido 47,3 K |

Los valores de T8 y T9 solo valen con los **materiales y el vehículo supuestos** del ejemplo
(γ = 0,869; q0 = 764 kW/m²; energía por cara = 101 kJ; t_f = 5,56 s). Si cambian propiedades,
hay que recalcularlos, no ajustarlos.

---

## 7. Trampas conocidas

1. **Graduación al revés.** Mi primer prototipo puso las celdas gruesas junto a la pista y dio
   un error de ~1,4 % que no bajaba al refinar. T1 lo detecta.
2. **Factor 2π.** Se omite en el ensamblaje. Un balance con potencias físicas sin 2π falla por
   un factor 6,28.
3. **Signo de la condición en `brake`.** El calor entra: la carga es +∫ q φ r ds.
4. **Crank–Nicolson con flujo escalón.** Oscila al inicio; por eso la Fase 3 compara esquemas
   en vez de fijar θ = 1/2 por costumbre.
5. **t_f fuera de la malla temporal.** Si t_f no es múltiplo de Δt, el flujo se corta a mitad de
   paso y el balance de energía no cierra.
6. **Nodos en r_a y r_b.** Si la malla no los tiene, `brake` no coincide con la pista y la
   energía de entrada cambia.
7. **Propiedades inventadas.** El prototipo usa valores supuestos; los resultados del informe
   necesitan materiales con fuente.
8. **Eje de simetría (r_i = 0).** Solo aparece en la prueba T6 (cilindro macizo). Allí el peso r
   anula el término singular 1/r, y en mi prototipo no hizo falta ningún tratamiento especial.

---

## 8. Relación con el cronograma del informe

- Semanas 10–12: Fases 0–5 (programa funcionando y validado).
- Semanas 14–16: Fase 6 (campos térmicos y convergencia de malla, objetivos 5 y 6).
- Semana 18: Fases 7–8 y redacción; el informe cita T6 y T7 como validación y la tabla de
  convergencia como resultado.

---

## ESTADO (actualizar al terminar cada fase)

- [x] Fase 0 — Inspección y arranque (2026-11-04). Se reutiliza la geometría y los campos de
  tamaño de Gmsh del script previo, movido a `legacy/mesh/`. pyproject + esqueleto + prueba de
  humo (7 passed). En Windows `python` no está en el PATH: usar `py -m venv .venv` y
  `.venv\Scripts\python`. Reporte: `reportes/F0_T1_inspeccion-arranque.md`.
- [x] Fase 1 — Malla (2026-11-05). `Mesh`, backends `structured` y `gmsh` con la misma ley
  h(d) (Threshold con DistMin = 0), seis grupos con aristas orientadas antihorario, CLI
  `mesh` con figura. T1 pasa en ambos backends (51 passed). Para Gmsh, "celda más pequeña toca
  brake" se prueba como "a < 2·h_min de brake + área media 10× menor que en el fondo" (ver
  reporte). Reporte: `reportes/F1_T1_malla.md`.
- [x] Fase 2 — Núcleo FEM (2026-11-05). `fem/assembly.py`: M, K (peso r o 1), H y F (Gauss 3
  puntos) vectorizados, propiedades por región. T2 pasa a 1e-12; T6 da 0,86 / 0,31 / 0,10 /
  0,030 % para nx = 10/20/40/80 (igual al prototipo). T6 usa un bucle CN propio hasta que exista
  `timestepping.py`. 74 passed. Reporte: `reportes/F2_T1_ensamblaje-fem.md`.
- [x] Fase 3 — Integración temporal y validación (2026-11-05). `timestepping.py`: esquema θ
  con LU por (θ, Δt), Rannacher, `time_grid` con t_f exacto, balance por componentes. T4
  (cierre ≤ 2,8e-14), T5 (8,9e-14) y T7 pasan; T6 migrada a `integrate`. θ por defecto: 1/2 +
  2 pasos Rannacher (tabla T7 con Δt = t_f/100: 1,5 % frente a 3,5 % con θ = 1 y 6 % con CN).
  Ojo: θ = 1 con Δt = t_f/100 mete un 1 % menos de energía (regla del extremo derecho).
  101 passed. Reporte: `reportes/F3_T1_integracion-temporal.md`.
- [x] Fase 4 — Cargas y materiales (2026-11-05). `loads.py` (braking, constant_flux, table;
  q = f(t)·s(r), F_brake ensamblado una vez) y `materials.py` (source/verified obligatorios,
  overrides ⇒ no verificado, γ auto desde la región bajo brake, advertencias). T3 pasa a 1e-12;
  γ = 0,8692, q0 = 764 kW/m², E_cara = 101 kJ. Pendiente para la Fase 5: llamar a
  `warn_unverified` al inicio de `run`. 149 passed. Reporte: `reportes/F4_T1_cargas-materiales.md`.
- [x] Fase 5 — Configuración, CLI y salidas (2026-11-05). `config.py`, `solver.run`,
  `postprocess` (sondas P1, CSV, JSON, PNG, VTK) y CLI `run --set`. config.yaml de resultados
  = config efectiva + copia de materiales.yaml (reproducible, probado). T8: 102,30 °C en
  3,389 s y 72,33 °C a 60 s; T9: 47,33 K ≤ 48,08 K. `summary.json` idéntico entre corridas
  salvo `compute_time_s`. Ojo: la referencia de T8 refleja θ = 1; con CN + Rannacher,
  T(60 s) = 72,8028 °C, fuera de 72,3 ± 0,5 por 0,003 °C. 174 passed.
  Reporte: `reportes/F5_T1_config-solver-salidas.md`.
- [x] Fase 6 — Convergencia y estudios paramétricos (2026-11-05). `convergence.py` +
  `python -m discfem study configs/estudio_malla.yaml` (27 s). Malla: ε = 0,18 → 0,067 →
  0,013 % (rf 1→8), nivel seleccionado rf = 2. Δt: p = 1,001 (θ = 1); Richardson da
  T(60 s) ≈ 72,81 °C, igual que CN + Rannacher. Barridos a, v0, m crecientes; k del disco solo
  como sensibilidad ±20 % (no hay materiales con fuente). 194 passed.
  Reporte: `reportes/F6_T1_convergencia-barridos.md`.
- [~] Fase 7 — Generalización (2026-11-05). Hecho: `geometry.type: profile` (polígono +
  grupo por lado), `regions` por polígonos (Gmsh OCC + fragment, malla conforme), presets por
  vehículo, `materials.custom`, `convection.per_group`. Criterios del PLAN cumplidos: el
  rectángulo da resultados idénticos bit a bit (referencia `tests/data/golden_rectangle.json`)
  y dos regiones con el mismo material = homogéneo a 1e-12. 229 passed.
  Catálogo: hierro fundido, pastilla y forro de Limpert (Tabla 3-1A, verificados) y preset
  `car`; `configs/ejemplo_auto.yaml` sin materiales no verificados (T_max = 106,70 °C). 231 passed.
  Acero AISI 410 (ficha Chilexpo; verified: false por error de unidad en c_p y densidad
  dudosa) y pastilla sinterizada de cobre híbrida (k = 45,1 y ρ = 5520 medidos por Zhou et al.
  2024, Fig. 6, S800; c = 418 de Goo y Na 2014; verified: false; ambas de trenes). Preset
  `motorcycle` = 410 + esa pastilla (γ = 0,483). Goo y Na ya no está en el preset.
  Barrido de pastilla (`configs/estudio_pastilla.yaml`, figura `gamma_map.pdf`): T − T0 es
  exactamente proporcional a γ (dispersión ~1e-14); con el 410 el pico va de 80,9 a 119,4 °C
  según la pastilla. 240 passed.
  Caso de moto `configs/ejemplo_moto.yaml` (MT-07: 2 × Ø 298 × 4,5 mm, 264 kg; supuestos
  marcados) y `estudio_moto.yaml`: T_max = 76,7 °C, malla p = 1,9, Δt p = 2,1. 243 passed.
  PENDIENTE: contrastar 410 y pastillas sinterizadas con fuentes medidas.
  Reportes: `reportes/F7_T1_perfil-regiones-materiales.md`, `reportes/F7_T2_catalogo-materiales-limpert.md`,
  `reportes/F7_T3_acero-410-preset-moto.md`, `reportes/F7_T4_pastilla-sinterizada-preset-moto.md`,
  `reportes/F7_T5_pastilla-cobre-zhou2024.md`, `reportes/F7_T6_barrido-pastilla.md`, `reportes/F7_T7_caso-moto.md`.
- [x] Fase 8 — Documentación (2026-11-05). `README.md` (instalación Windows/Linux, comandos,
  YAML, materiales, salidas, validación) y `requirements-lock.txt` (pip freeze). Verificado en
  dos entornos limpios (Windows, Python 3.13): `pip install -e ".[dev]"` → 240 passed; desde el
  lock → 240 passed, `pip check` sin conflictos. Linux no verificado.
  Reporte: `reportes/F8_T1_documentacion-reproducibilidad.md`.
