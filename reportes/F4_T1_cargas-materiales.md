# F4 · T1 — Cargas (`loads.py`) y materiales (`materials.py`)

Todas las pruebas pasaron (149 passed) sin relajar tolerancias. Los valores de los materiales
son los SUPUESTOS DEL PROTOTIPO copiados de `PLAN.md` §4, con `verified: false`. No se inventó
ninguna propiedad.

## 1. Qué se pidió

Fase 4 de `PLAN.md` (§5, "Fase 4 — Cargas y materiales"):

- `loads.py`: modelo `braking` de la §3.4, `constant_flux` y `table` con interpolación lineal.
- `materials.py`: librería con `source` y `verified` obligatorios, efusividad, $\gamma$ automático
  desde el material de la región adyacente a `brake`, y advertencias de rango.
- Si falta `source`, el programa debe fallar con un mensaje en español. Si `verified: false`,
  debe advertir al inicio de cada corrida.
- Prueba T3 (§6). Criterio de aceptación: $\gamma \approx 0{,}87$ en el caso de referencia.

## 2. Qué se hizo

- `configs/materiales.yaml`: nuevo. `gray_cast_iron` y `organic_pad` con los valores de la §4,
  `source: "PENDIENTE: citar fuente"` y `verified: false`.
- `src/discfem/materials.py`: nuevo. `Material`, `load_library`, `resolve_library_path`,
  `apply_override`, `resolve_region_materials`, `warn_unverified`, `check_range`,
  `partition_coefficient`, `brake_regions` y `resolve_gamma`.
- `src/discfem/loads.py`: nuevo. `BrakingLoad`, `ConstantFluxLoad`, `TableLoad`, `make_load` y
  `brake_load_vector`.
- `tests/test_loads_materials.py`: nuevo. T3, más pruebas de configuración, materiales,
  $\gamma$, overrides y advertencias (48 casos).
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

**De la frenada al flujo de calor.** El vehículo, de masa $m$, frena con desaceleración
constante $a$ desde $v_0$. Por eso

$$
v(t) = v_0 - a t, \qquad t_f = \frac{v_0}{a}.
$$

La potencia de frenado del vehículo es $m\, a\, v(t)$. De ella, una fracción $\varphi$ va al eje
estudiado, que se reparte entre $n_d$ discos de $n_c$ caras cada uno. La potencia en UNA cara es:

$$
P_f(t) = \frac{\varphi\, m\, a\, v(t)}{n_d\, n_c}
$$

En la interfaz disco–pastilla se genera $P_f$, pero no todo entra al disco. Una fracción
$\gamma$ va al disco y $1 - \gamma$ a la pastilla. Para dos cuerpos semi-infinitos en contacto,
el reparto lo fija la efusividad $e = \sqrt{k \rho c_p}$, que mide cuánto calor absorbe un
material ante un cambio brusco de temperatura en su superficie:

$$
\gamma = \frac{e_{\text{disco}}}{e_{\text{disco}} + e_{\text{pastilla}}}
$$

El hierro, con $e \approx 12\,869$, absorbe mucho más que la pastilla orgánica
($e \approx 1936$): $\gamma = 0{,}869$.

**Repartir en $r$.** Ese calor se distribuye sobre la pista con un perfil $g(r)$: presión
uniforme ($g = 1$) o desgaste uniforme ($g = 1/r$). El flujo local es proporcional a
presión × velocidad de deslizamiento, y la velocidad es $\omega r$, de ahí el factor $r$. Se
normaliza con $I_g = \int_{r_a}^{r_b} g(r)\, r^2\, dr$ para que el anillo completo reciba
exactamente $\gamma P_f$:

$$
q(r, t) = \frac{\gamma\, \varphi\, m\, a\, v(t)}{2\pi\, n_d\, n_c}\cdot\frac{g(r)\, r}{I_g},
\qquad 2\pi \int_{r_a}^{r_b} q(r,t)\, r\, dr = \gamma\, P_f(t)
$$

Casos cerrados (§3.4):

$$
\text{presión uniforme: } q = \frac{3\,\gamma\,\varphi\, m\, a\, v(t)\, r}{2\pi\, n_d\, n_c\,(r_b^3 - r_a^3)},
\qquad
\text{desgaste uniforme: } q = \frac{\gamma\,\varphi\, m\, a\, v(t)}{\pi\, n_d\, n_c\,(r_b^2 - r_a^2)}
$$

Con desgaste uniforme, $q$ no depende de $r$ (764 kW/m² en $t = 0$). Con presión uniforme, $q$
crece linealmente con $r$ (542 → 949 kW/m²). El coeficiente de fricción $\mu$ no aparece: con la
desaceleración prescrita se cancela (§3.4).

**Separación $q = f(t)\, s(r)$.** Los tres modelos se escriben así. `brake_load_vector`
ensambla una sola vez

$$
F_{s,i} = \int_{\Gamma_{\text{brake}}} s(r)\, \varphi_i\, r\, ds
$$

(con `assemble_edge_load` de la Fase 2) y devuelve la función
$\mathbf{F}_{\text{brake}}(t) = f(t)\, \mathbf{F}_s$, que es lo que recibe `integrate` (Fase 3).
Para `braking`, $f(t) = v(t)/v_0 = 1 - t/t_f$, que se anula en $t_f$ y después.

**Materiales.** `load_library` valida cada entrada:

- claves obligatorias `k`, `rho`, `cp`, `source` y `verified`, con mensaje en español si falta
  alguna;
- `source` no vacío y `verified` booleano;
- propiedades $> 0$;
- no se admiten claves desconocidas (atrapa erratas).

Si una propiedad está fuera de los rangos de sanidad de la §4, emite `PropertyRangeWarning` sin
abortar. `warn_unverified` emite `UnverifiedMaterialWarning` con la lista de materiales sin
verificar y sus fuentes.

`resolve_gamma("auto")` toma los triángulos que tienen una arista en `brake`, mira su región y
su material, y exige que sea uno solo. Con dos materiales distintos bajo la pista no hay un
$\gamma$ único, y pide un $\gamma$ numérico.

| Fórmula (PLAN.md §3.4) | Código |
|---|---|
| $v(t) = v_0 - a t$; $t_f = v_0/a$ | `BrakingLoad.velocity`, `loads.py:105`; `t_f`, `:101` |
| $P_f = \varphi\, m\, a\, v / (n_d n_c)$ | `BrakingLoad.power_per_face`, `loads.py:116` |
| $I_g = \int g\, r^2\, dr$ (cerrada para $g = 1$ y $g = 1/r$) | `BrakingLoad.I_g`, `loads.py:129-133` |
| $q = \dfrac{\gamma\varphi m a v}{2\pi n_d n_c}\cdot\dfrac{g\, r}{I_g}$ | `BrakingLoad.shape` (`:144-151`) × `time_factor` (`:140`); `Load.q`, `:46` |
| $q = 0$ para $t > t_f$ | `velocity` (`np.where(t <= t_f, …)`), `loads.py:105-108` |
| $\gamma = \dfrac{\sqrt{k\rho c}_d}{\sqrt{k\rho c}_d + \sqrt{k\rho c}_p}$ | `Material.effusivity`, `materials.py:77`; `partition_coefficient`, `:230` |
| $\gamma$ automático desde la región bajo `brake` | `brake_regions` (`:235`), `resolve_gamma` (`:249`) |
| `source` y `verified` obligatorios (CLAUDE.md regla 4) | `_material_from_entry`, `materials.py:95` |
| Rangos de sanidad (§4) | `SANITY_RANGES`, `check_range`, `materials.py:82` |

## 4. Decisiones

- **$q = f(t)\, s(r)$ en vez de $q(r, t)$ genérico.** Así el vector de carga se ensambla una vez
  y en cada paso solo se escala, en vez de reensamblarlo ~1000 veces. Los tres modelos del PLAN
  cumplen esta separación. Si algún día se necesita un $q$ que no se separe (por ejemplo, una
  pista que se mueve), habrá que reensamblar.
- **`constant_flux` y `table` dan el flujo que entra al disco, sin $\gamma$.** En esos modelos
  el usuario prescribe el flujo directamente (por ejemplo, T7). Multiplicar por $\gamma$ sería
  sorpresivo. Queda escrito en el docstring del módulo.
- **`table` es uniforme en $r$, empieza en $t = 0$ y su último instante es $t_f$.** Descarté leer
  la tabla desde un CSV para no añadir formato; se escribe dentro del YAML
  (`points: [[t, q], …]`).
- **Override ⇒ `verified: false`.** Un valor que se sobrescribe en `materials.overrides` sale
  del YAML de la corrida, no de la librería. Para no "lavar" un valor sin fuente, el material
  resultante queda no verificado, salvo que el override traiga su propio `source` y `verified`.
- **Advertencias con `warnings.warn` y categorías propias.** Así las pruebas las verifican con
  `pytest.warns`, y la Fase 5 puede mostrarlas o capturarlas para `summary.json`.
- **`ValueError` también para errores de tipo en la librería.** Ruff (TRY004) sugería
  `TypeError`. Lo mantuve como `ValueError` con `noqa` y el motivo, para que todos los errores
  de configuración se manejen igual.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
```

Salida real:

```
All checks passed!
149 passed in 4.06s
```

Caso de referencia (`configs/ejemplo_frenado.yaml` + `configs/materiales.yaml`, malla Gmsh),
salida real de un script que llama a las funciones de la fase, incluida la advertencia tal como
se emite:

```
...f4.py:11: UnverifiedMaterialWarning: Materiales NO verificados: gray_cast_iron (source: PENDIENTE: citar fuente); organic_pad (source: PENDIENTE: citar fuente). Los resultados no sirven para el informe hasta citar fuentes y poner verified: true.
e_disco = 12868.57, e_pastilla = 1936.49, gamma = 0.869201
[uniform_wear] t_f = 5.5560 s, P_f(0) = 36461.2 W, E_cara = 101289.4 J, q(r_a,0) = 764236.5, q(r_b,0) = 764236.5 W/m²
   t=0.0: 2π ΣF = 31692.142923 W, γ P_f = 31692.142923 W, err rel = 1.1e-16
   t=2.0: 2π ΣF = 20283.884131 W, γ P_f = 20283.884131 W, err rel = 1.8e-16
   t=4.0: 2π ΣF = 8875.625340 W, γ P_f = 8875.625340 W, err rel = 2.0e-16
[uniform_pressure] t_f = 5.5560 s, P_f(0) = 36461.2 W, E_cara = 101289.4 J, q(r_a,0) = 542361.4, q(r_b,0) = 949132.5 W/m²
   t=0.0: 2π ΣF = 31692.142923 W, γ P_f = 31692.142923 W, err rel = 1.1e-16
   t=2.0: 2π ΣF = 20283.884131 W, γ P_f = 20283.884131 W, err rel = 1.8e-16
   t=4.0: 2π ΣF = 8875.625340 W, γ P_f = 8875.625340 W, err rel = 2.0e-16
```

## 6. Resultados

| Magnitud | Medido aquí | Referencia (§6 y criterio de la Fase 4) |
|---|---|---|
| $\gamma$ (materiales supuestos) | 0,869201 | $\approx 0{,}87$; "$\gamma = 0{,}869$" |
| $t_f$ | 5,5560 s | 5,56 s |
| $q_0$ (desgaste uniforme, $t = 0$) | 764 236 W/m² | 764 kW/m² |
| Energía por cara (antes de $\gamma$) | 101 289 J | 101 kJ |
| T3: $2\pi\int q\, r\, dr$ frente a $\gamma P_f(t)$, ambos perfiles | error $\le 10^{-12}$ (prueba) | $10^{-12}$ relativo |
| T3: fórmulas cerradas de la §3.4 | coinciden a $10^{-12}$ | $10^{-12}$ relativo |
| T3 con el vector FEM: $2\pi\sum_i F_{\text{brake},i}(t)$ frente a $\gamma P_f(t)$ | $1{,}1\times10^{-16}$ – $2{,}0\times10^{-16}$ | (adicional) |

El último renglón conecta esta fase con la anterior. Gauss de 3 puntos es exacto para ambos
perfiles, y los nodos en $r_a$ y $r_b$ (Fase 1) hacen que la pista de la malla coincida con la
física. Por eso la potencia que entra al modelo FEM, multiplicada por $2\pi$, es la potencia
física a precisión de máquina.

## 7. Limitaciones y riesgos

- **Los materiales siguen sin fuente** (`verified: false`). Todos los números del caso de
  referencia ($\gamma$, T8, T9) valen solo con estos supuestos. Antes del informe hay que citar
  fuentes y recalcular, no ajustar.
- **La advertencia "al inicio de cada corrida" todavía no está conectada.** `warn_unverified`
  existe y está probada, pero la corrida (`solver.run`) es de la Fase 5, que debe llamarla
  primero.
- **El $\gamma$ por efusividades supone contacto perfecto entre cuerpos semi-infinitos** con
  propiedades constantes. La resistencia térmica de contacto está fuera de alcance (CLAUDE.md).
- **`table` y `constant_flux` son uniformes en $r$.** Si se quiere una tabla con perfil de
  presión, hay que extender `TableLoad`.
- **La ruta de la librería es relativa al directorio actual** y, si no existe ahí, al de la
  config (`resolve_library_path`). La Fase 5 debe usar esta función para que `run` funcione
  desde cualquier carpeta.

## 8. Para revisar y para el informe

1. **La cadena física completa:**
   $m\, a\, v \to \varphi \to 1/(n_d n_c) \to \gamma \to g(r)\, r / I_g$. Cada factor tiene un
   sentido, y $2\pi\int q\, r\, dr = \gamma P_f$ es la prueba de que el reparto en $r$ no crea ni
   destruye potencia.
2. **Qué es $\gamma$ y por qué sale de las efusividades.** Es la regla de reparto para dos sólidos
   semi-infinitos en contacto. Con hierro y pastilla orgánica, ~87 % del calor entra al disco.
3. **Por qué $\mu$ y el ángulo de la pastilla no son entradas:** con desaceleración prescrita se
   cancelan (§3.4).
4. **Desgaste uniforme frente a presión uniforme:** $q$ constante frente a $q \propto r$, misma
   potencia total. En la Fase 6 se podrá ver cómo cambia el campo de $T$.
5. **Para el informe:** la tabla de §6 ($\gamma$, $t_f$, $q_0$, energía por cara) sirve como
   tabla de "parámetros derivados del caso de referencia". Debe ir acompañada de la advertencia de
   que los materiales aún no están verificados.
