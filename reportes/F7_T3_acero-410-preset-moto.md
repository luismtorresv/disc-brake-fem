# F7 · T3 — Acero inoxidable AISI 410 (ficha Chilexpo) y preset `motorcycle`

Todas las pruebas pasaron (233 passed) sin relajar tolerancias. El AISI 410 quedó con
**`verified: false`**: la ficha tiene un error de unidades en el calor específico, una densidad
dudosa y su propia advertencia de no estar verificada (§3). El preset `motorcycle` define solo
el disco, porque sigue faltando una fuente para la pastilla sinterizada.

## 1. Qué se pidió

Continuación de la Fase 7 (§5 de `PLAN.md`):

- Primero dejaste en `legacy/` una ficha de **AISI 4140** (acero Cr-Mo, UNS G41400). No la usé
  porque no es el 410 y no trae calor específico.
- Luego indicaste "Usa el 410, ya monté la ficha correcta": `legacy/410.pdf`.

## 2. Qué se hizo

- `configs/materiales.yaml`:
  - nueva entrada `stainless_steel_410`, con nota sobre los problemas de la ficha;
  - nuevo preset `motorcycle` (solo disco).
- `src/discfem/solver.py`: si `gamma: auto` y ni el preset ni la config definen pastilla, falla
  con un mensaje que dice cuál preset y qué hacer.
- `tests/test_profile_regions.py`: 2 pruebas nuevas, sobre los valores del 410 y sobre que el
  preset `motorcycle` exija pastilla y advierta que el material no está verificado.
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.
- Los PDF de `legacy/` no se versionan (`.gitignore`).

## 3. Cómo funciona

Se extrajo el texto de `legacy/410.pdf` con `pypdf`, instalado en una carpeta temporal. Es la
ficha "Aleación 410" de **Chilexpo SpA**, UNS S41000 / DIN X15Cr13 / W.N. 1.4006.

Sección "Propiedades físicas", copiada tal cual:

```
Densidad 0.266 lb/in3
Calor Especifico 0 - 100°C: 0,46  Btu/lb/ºF
Conductividad Térmica  212 °F / 100 °C   187 Btu•in/hr•ft²•°F   26,9 W/m•K
```

y la advertencia final: "Esta información técnica ha sido recompilada de diversas fuentes en
línea, incluyendo ATI®, SSC®, y Outokumpu® entre otras [...] no ha sido verificada de forma
independiente por CHILEXPO SPA."

| Propiedad | Ficha | Valor cargado | Comprobación |
|---|---|---|---|
| $k$ | 26,9 W/(m·K) a 100 °C | 26,9 | 187 BTU·in/(h·ft²·°F) × 0,1442 = 26,97 W/(m·K) ✓ |
| $c_p$ | "0,46 Btu/lb/ºF" | 460 J/(kg·K) | 0,46 BTU/(lb·°F) serían 1926 J/(kg·K), imposible para un acero (~450–500). 0,46 es el valor en kJ/(kg·K), equivalente a 0,11 BTU/(lb·°F), como en la ficha de Cleveland-Cliffs |
| $\rho$ | 0,266 lb/in³ | 7362,85 kg/m³ | $= 0{,}266 \times 27\,679{,}9$. Es menor que en otras fuentes del 410 (≈ 0,28 lb/in³ ≈ 7750). Se respeta la ficha y queda anotado |

Con la pastilla de Limpert, las efusividades $e = \sqrt{k\rho c}$ son: 410 = 9545, hierro de
Limpert = 12 115, pastilla = 2146; por tanto

$$
\gamma(\text{410 + pastilla de Limpert}) = \frac{9545}{9545 + 2146} = 0{,}8164 .
$$

**Preset `motorcycle`.** Define `disc: stainless_steel_410` sin pastilla. Con `gamma: auto`, la
config debe dar `materials.pad`; si no, el error es: "gamma: auto necesita una pastilla y el
preset 'motorcycle' no la define; indique materials.pad (o un gamma numérico)."

## 4. Decisiones

- **No usé la ficha del 4140.** Es otro acero (~1 % de Cr frente a 11,5–13,5 %), no es el de los
  discos de moto y le falta $c_p$.
- **`verified: false` para el 410.** La ficha se declara no verificada, tiene un error de unidad
  y una densidad atípica. Que cada corrida lo advierta es lo honesto. Para pasarlo a
  `verified: true` conviene contrastarlo con la ficha del fabricante (ATI, Outokumpu) o con el
  ASM Handbook.
- **$c_p = 460$ J/(kg·K)**, y no 1926, aunque la etiqueta de la ficha diga BTU. No es inventar un
  valor: es el mismo número de la ficha (0,46) con la única unidad físicamente posible, y
  coincide con 0,11 BTU/(lb·°F) de otras fichas del 410. La interpretación queda escrita en el
  `source`.
- **$\rho$ tal cual (7363 kg/m³)**, aunque parezca baja. Corregirla a 7750 sería usar un valor
  que no está en la fuente que me diste.
- **El preset `motorcycle` sin pastilla por defecto.** Asignar en silencio la pastilla de auto
  de Limpert (base asbesto) a una moto, que usa pastillas sinterizadas con $k$ mucho mayor,
  falsearía $\gamma$. Por eso la config debe elegir la pastilla de forma explícita.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
```

Salida real:

```
233 passed in 18.99s
All checks passed!
```

Conversiones (salida real):

```
rho 7362.854652913201
k desde BTU in/hr ft2 F 26.970615243
0.46 BTU/lbF en J/kgK 1925.928  | 0.11 BTU/lbF 460.548
alpha 410 7.942335361248452e-06 e410 9545.049160436416 e_hierro 12114.876476877509 gamma 410+pad_limpert 0.8164237361701493
```

Comparación con el mismo caso de auto (`ejemplo_auto.yaml`), cambiando solo el disco y con
`pad_limpert` explícita (salida real):

```
gray_cast_iron_limpert   gamma=0.8495 Tmax=106.70 en t=3.500  T60=75.88  E_in=86028 J
stainless_steel_410      gamma=0.8164 Tmax=119.44 en t=2.834  T60=70.65  E_in=82678 J
```

## 6. Resultados

Es una comparación de **material** sobre la geometría y el vehículo del auto, no un caso de
moto:

| Disco | $k$ [W/(m·K)] | $\rho c$ [J/(m³·K)] | $\gamma$ | $T_{\max}$ | $t(T_{\max})$ | $T(60\ \mathrm{s})$ |
|---|---|---|---|---|---|---|
| Hierro fundido (Limpert) | 48,46 | $3{,}03\times10^{6}$ | 0,8495 | 106,70 °C | 3,50 s | 75,88 °C |
| Inoxidable AISI 410 (Chilexpo) | 26,9 | $3{,}39\times10^{6}$ | 0,8164 | 119,44 °C | 2,83 s | 70,65 °C |

- El 410 conduce casi la mitad. El calor se queda cerca de la superficie, así que el pico es
  12,7 °C mayor y llega antes.
- Al final, con el disco casi uniforme, manda la energía sobre la capacidad. Entra menos
  energía ($\gamma$ menor, −3,9 %) en un disco con más capacidad volumétrica (+12 %), y eso da
  5,2 °C menos.

## 7. Limitaciones y riesgos

- **El 410 no está verificado**, por la unidad de $c_p$, la densidad atípica y la advertencia de
  la propia ficha. Antes del informe hay que contrastarlo con una ficha de fabricante.
- **Las propiedades del 410 son a 0–100 °C** ($c_p$) y a 100 °C ($k$), y además dependen del
  tratamiento térmico. Los discos de moto van templados; la ficha no distingue el estado para
  las propiedades físicas.
- **Falta la pastilla sinterizada.** Con ella, $\gamma$ bajaría mucho, porque una pastilla
  metálica tiene $k$ de decenas de W/(m·K), y el pico calculado para la moto cambiaría. Hasta
  tener esa fuente, cualquier caso de moto debe declarar qué pastilla usa.
- **No hay un caso de moto con su geometría y vehículo propios** (disco más delgado, masa y
  reparto de una moto). Esos datos no son propiedades de material, pero tampoco tengo fuente
  para ellos; quedan por definir contigo.

## 8. Para revisar y para el informe

1. **Cómo detectar un error de unidades en una ficha:** el orden de magnitud físico ($c_p$ de
   aceros ~450–500 J/(kg·K)) y la coherencia entre columnas (BTU frente a SI).
2. **Por qué el 410 se calienta más en el pico pero menos al final:** conductividad frente a
   capacidad térmica volumétrica y $\gamma$.
3. **Por qué el preset de moto no trae pastilla:** $\gamma$ depende fuertemente de la efusividad
   de la pastilla, y asignar la de un auto ocultaría esa dependencia.
4. **Fuentes para el informe:** Limpert (Tabla 3-1A) para hierro y pastilla; ficha Chilexpo
   "Aleación 410" (con sus advertencias) para el acero de moto.
5. **Para el informe:** la tabla de §6 como comparación de materiales del disco. Cuando haya
   pastilla sinterizada y datos de moto, un caso `ejemplo_moto.yaml` completo.
