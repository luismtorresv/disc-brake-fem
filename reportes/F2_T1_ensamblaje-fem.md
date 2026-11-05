# F2 · T1 — Núcleo FEM: ensamblaje vectorizado de M, K, H y F

Todas las pruebas pasaron (74 passed) sin relajar tolerancias. Hubo una falla intermedia en una
prueba nueva: el error estaba en la propia prueba, no en el código, y se explica en §4.

## 1. Qué se pidió

Fase 2 de `PLAN.md` (§5, "Fase 2 — Núcleo FEM"): implementar `fem/assembly.py` según la §3.3.

- Matrices elementales cerradas para M y K, con peso r y con peso 1.
- Matriz H y vector F por aristas, F con Gauss de 3 puntos.
- Todo vectorizado y ensamblado con `scipy.sparse` (COO → CSR), con propiedades por región.
- Pruebas T2 (matrices contra cuadratura de alto orden en triángulos aleatorios, K·1 = 0,
  M simétrica y definida positiva, 1ᵀM1 = ρc ∫ r dA) y T6 (modo propio axisimétrico del
  cilindro).

## 2. Qué se hizo

- `src/discfem/fem/assembly.py`: nuevo. `element_geometry`, `element_matrices`, `per_element`,
  `assemble_stiffness_mass`, `edge_matrices`, `assemble_robin`, `assemble_edge_load`.
- `src/discfem/fem/__init__.py`: nuevo. Exporta las funciones anteriores.
- `src/discfem/mesh/structured.py`: modificado. Se separó `tensor_mesh(r, z, ia, ib)`, que
  construye una malla a partir de coordenadas 1D explícitas. T6 la necesita (cilindro con
  r_i = 0 y malla uniforme nx × nz), y T7 la usará en la Fase 3. `structured_mesh` la llama sin
  cambiar su resultado: las 51 pruebas de T1 siguen pasando.
- `tests/test_assembly.py`: nuevo. T2 y pruebas de H, F, regiones y errores (21 casos).
- `tests/test_axisymmetric_mode.py`: nuevo. T6 (2 casos sobre nx = 10, 20, 40, 80).
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

**De la ecuación al sistema.** Se multiplica la ecuación de calor de la §3.1 por una función de
prueba φ_i y se integra sobre el volumen del disco. En coordenadas cilíndricas,
dV = r dr dθ dz. Como nada depende de θ, la integral en θ da un factor 2π común a todos los
términos, que se cancela (CLAUDE.md: se omite). Queda el peso r en cada integral. Después se
integra por partes el término de difusión, y aparecen las integrales de frontera. Con
T ≈ Σ_j T_j φ_j se llega a la §3.2:

M Ṫ + (K + H) T = F(t)

- M_ij = ∫ ρc φ_i φ_j r dA: masa o capacidad térmica. Aparece por la derivada temporal; es lo
  nuevo frente al caso estacionario visto en clase.
- K_ij = ∫ k ∇φ_i·∇φ_j r dA: rigidez o conducción. Es la misma del caso estacionario, pero
  con peso r.
- H_ij = ∫ h φ_i φ_j r ds: convección en las aristas de los grupos `convection`.
- F_i = ∫ q φ_i r ds + ∫ h T∞ φ_i r ds: flujo entrante en `brake` más el aporte de T∞.

**Matrices elementales (fórmulas cerradas, sin cuadratura).** En un triángulo P1, ∇φ_i =
(b_i, c_i)/(2A) es constante, así que el integrando de K es la constante (b_i b_j + c_i c_j)/(4A²)
por r. Como r es lineal, ∫ r dA = A·r̄, y por eso

K^e_ij = k r̄ (b_i b_j + c_i c_j) / (4A)

En M, el integrando φ_i φ_j r es cúbico. Se escribe r = Σ_k r_k φ_k y se usa la fórmula exacta
∫ φ_i φ_j φ_k dA = A·C_ijk, con C_ijk = (1+δ_ij)(1+δ_ik+δ_jk)/60. Eso da 6/60 = 1/10 si los tres
índices son iguales, 2/60 = 1/30 si coinciden dos y 1/60 si los tres son distintos. Así:

M^e_ij = ρc A Σ_k C_ijk r_k

Desarrollada, la fórmula es exactamente la de la §3.3: M_ii = ρc A[r_i/10 + (r_j+r_l)/30] y
M_ij = ρc A[(r_i+r_j)/30 + r_l/60]. Hay una prueba que compara literalmente con esas
expresiones. En modo plano, r_k → 1 da la masa clásica A(1+δ_ij)/12.

**Aristas.** Sobre una arista (a, b) de longitud L, con φ_a = 1 − s, φ_b = s y
r = (1−s) r_a + s r_b, la integral exacta ∫ h φ_i φ_j r ds da la matriz 2×2 de la §3.3. Por
ejemplo, H_aa = hL(r_a/4 + r_b/12). Para F se usa Gauss–Legendre de 3 puntos por arista: es
exacto hasta grado 5, y q ∝ r por r por φ es de grado 3.

**Ensamblaje vectorizado.** `element_matrices` devuelve arreglos (E, 3, 3) de una sola vez.
`_scatter` repite la conectividad para generar los índices (fila, columna) de las 9 entradas de
cada elemento, y `coo_matrix(...).tocsr()` suma los duplicados. Eso es exactamente el
ensamblaje Σ_e de clase, sin bucles de Python. F usa `np.add.at`. Las propiedades k y ρc pueden
ser un escalar, un arreglo (E,) o un diccionario {región: valor}; `per_element` las convierte a
(E,) con `region`.

| Fórmula (PLAN.md) | Código |
|---|---|
| b_i, c_i, A (§3.3) | `element_geometry`, `fem/assembly.py:32-47` |
| K^e_ij = k r̄ (b_i b_j + c_i c_j)/(4A) (§3.3) | `element_matrices`, `fem/assembly.py:68-70`; k en `:118` |
| M^e_ii, M^e_ij con 1/10, 1/30, 1/60 (§3.3) | tensor `_C3`, `fem/assembly.py:25`; `Me`, `:71` |
| Modo plano: r → 1 (§3.2) | `rk = ones` en `element_matrices`; rama `else` de `edge_matrices`, `:144` |
| H^e = hL[[r_a/4 + r_b/12, (r_a+r_b)/12], …] (§3.3) | `edge_matrices`, `fem/assembly.py:124-145`; h en `assemble_robin`, `:148` |
| F^e_i = ∫ q r φ_i ds, Gauss 3 puntos (§3.3) | `_GL_S/_GL_W`, `:28-29`; `assemble_edge_load`, `:163-190` |
| Signo en `brake`: +∫ q φ r ds (CLAUDE.md) | `assemble_edge_load` suma q tal cual (q > 0 = calor que entra) |
| Ensamblaje Σ_e (COO → CSR) | `_scatter`, `fem/assembly.py:97-102` |
| Propiedades por región | `per_element`, `fem/assembly.py:75-94` |

## 4. Decisiones

- **El tensor C_ijk en vez de escribir las 9 entradas a mano.** Es una sola línea y no tiene
  casos especiales. Para que no se pierda la conexión con el PLAN, una prueba
  (`test_element_mass_closed_formulas_from_plan`) verifica literalmente las fórmulas con 1/10,
  1/30 y 1/60.
- **Referencia de T2 independiente.** La cuadratura de referencia calcula los gradientes
  invirtiendo el jacobiano de la transformación afín, no con b_i y c_i. Usa Gauss colapsado de
  Duffy 8×8 (exacto hasta grado 14). Si la referencia usara las mismas b_i y c_i, un error de
  signo en ellas pasaría inadvertido.
- **T6 sin `timestepping.py`.** T6 necesita avanzar en el tiempo, pero el módulo es de la
  Fase 3. Para no adelantarla, la prueba escribe el bucle de Crank–Nicolson directamente
  (5 líneas con `splu`). Cuando exista `timestepping.py`, conviene cambiar T6 para que lo use.
- **F_conv sin función propia.** Como Σ_j φ_j = 1, ∫ h T∞ φ_i r ds = (H·T∞1)_i, y hay una prueba
  que lo verifica. También se puede obtener con `assemble_edge_load(…, h·T∞)`. La Fase 5 elegirá
  cuál usar.
- **`tensor_mesh` en `structured.py`** en vez de un generador nuevo. Es la misma rejilla; solo se
  separó la parte que recibe las coordenadas.
- **Falla intermedia (error en la prueba, no en el código).**
  `test_element_matrices_with_vertex_on_axis` ponía r = 0 en un vértice de triángulos aleatorios
  ya orientados, y eso invertía la orientación de algunos. `element_matrices` los rechazó con
  `ValueError`, que es lo correcto. Se corrigió la prueba reorientando después de modificar y
  descartando los casi degenerados. No se tocó ninguna tolerancia.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
```

Salida real:

```
74 passed in 3.14s
All checks passed!
```

Primera corrida, antes de corregir la prueba (salida real, recortada a las últimas líneas):

```
>           raise ValueError("Hay tri\xe1ngulos con \xe1rea ≤ 0: la malla debe ser antihoraria.")
E           ValueError: Hay tri\xe1ngulos con \xe1rea ≤ 0: la malla debe ser antihoraria.
FAILED tests/test_assembly.py::test_element_matrices_with_vertex_on_axis - Va...
1 failed, 73 passed in 3.75s
```

Tabla de T6, generada con un script que llama a `run_mode` de la prueba. Salida real:

```
  nx   nz  err(t=10)/0.5 [%]  max err todos t [%]
  10    2             0.8643               1.6783
  20    4             0.3069               0.6145  razón 2.82
  40    8             0.0989               0.2021  razón 3.10
  80   16             0.0302               0.0627  razón 3.27
```

Tiempo de ensamblaje (M, K, H y F) con mallas Gmsh del caso de referencia. Salida real:

```
refine 1: 929 nodos, 1618 triángulos, ensamblaje 2.1 ms
refine 4: 11051 nodos, 21148 triángulos, ensamblaje 19.4 ms
refine 8: 42433 nodos, 82956 triángulos, ensamblaje 82.0 ms
```

## 6. Resultados

**T2** (tolerancia 1e-12 relativa al máximo de cada matriz; todas pasan):

- K^e y M^e, con peso r y con peso 1, contra cuadratura de alto orden en 200 triángulos
  aleatorios, y en 50 más con un vértice sobre el eje r = 0.
- H^e contra Gauss de 10 puntos en 100 aristas aleatorias (peso r y peso 1).
- K·1 = 0, K y M simétricas, y 1ᵀM1 = ρc Σ A r̄, en las mallas structured y gmsh.
- M definida positiva: el menor autovalor es > 0. K semidefinida, con un solo modo nulo (la
  constante).
- Σ F para q ∝ r igual a c(r_b³ − r_a³)/3, y coincidencia nodo a nodo con Gauss de 10 puntos.
- Σ H igual a h ∫ r ds sobre los grupos con convección, y F_conv = H·T∞1.

**T6** frente a la §6:

| nx | error en t = 10 s (medido aquí) | referencia del PLAN |
|---|---|---|
| 10 | 0,86 % | 0,86 % |
| 20 | 0,31 % | 0,31 % (criterio < 1 %) |
| 40 | 0,10 % | 0,10 % |
| 80 | 0,030 % | 0,03 % |

El error decrece en forma monótona y coincide con el prototipo en todas las cifras que este
reporta. La razón entre niveles sube de 2,8 a 3,3 y se acerca al 4 que corresponde a O(h²) de
P1. En mallas tan gruesas (nz = 2) todavía no se está en el régimen asintótico.

## 7. Limitaciones y riesgos

- **El error de T6 se mide en t = 10 s.** Así salen los números del PLAN. El máximo sobre todos
  los instantes es aproximadamente el doble (0,61 % con nx = 20, todavía < 1 %); la prueba
  exige ambos por debajo de 1 %. La §6 no dice cuál de los dos usa; supuse t = 10 s porque
  reproduce los valores medidos.
- **T6 y la conservación de la energía usan su propio bucle de Crank–Nicolson.** No se probó
  todavía el esquema θ general ni Δt variable (Fase 3).
- **Propiedades constantes por elemento.** k(T) queda fuera (Fase 7: "nada de k(T)").
- **`assemble_edge_load` evalúa q(r, z) en los puntos de Gauss.** Si q tiene un salto dentro
  de una arista, la regla deja de ser exacta. Con `brake` alineada a los nodos r_a y r_b
  (Fase 1) eso no pasa en el caso de referencia.
- **No se verificó el caso r_i = 0 con convección o flujo en el eje.** Ahí H y F tienen peso
  r = 0, lo cual es correcto; solo se probó el caso aislado.

## 8. Para revisar y para el informe

1. **De dónde sale el peso r y por qué se puede omitir el 2π.** Es común a todos los términos y
   se cancela, pero cualquier potencia física (W) debe multiplicarse por 2π (§7, trampa 2).
2. **Qué valida cada prueba.** T2 protege los coeficientes 1/10, 1/30 y 1/60. T6 valida que el
   peso r esté bien puesto en conjunto: la solución exacta es J0 solo si la ecuación es la
   axisimétrica. Una prueba de energía (1ᵀM1) solo valida la suma de las entradas, como advierte
   la Fase 2.
3. **K·1 = 0 implica conservación de la energía.** Con el borde aislado, Σ (M T)_i no cambia.
   T6 lo verifica a 1e-12, y es la base del balance de la §3.6 (Fase 3).
4. **Cómo se ensambla sin bucles.** Bloques elementales (E, 3, 3), índices repetidos y una
   `coo_matrix` que suma los duplicados.
5. **Para el informe:** la tabla de T6 (nx contra error, con la razón entre niveles) es la
   validación espacial del solver axisimétrico. Se puede graficar en escala log-log (error
   contra h) con una pendiente de referencia igual a 2.
