# F7 · T5 — Pastilla sinterizada de cobre (Zhou et al. 2024) en el preset `motorcycle`

Todas las pruebas pasaron (236 passed) sin relajar tolerancias. La nueva pastilla es
**híbrida**: $k$ y $\rho$ medidos por Zhou et al. (2024), y $c$ de Goo y Na (2014), porque Zhou
no mide el calor específico. Queda con `verified: false`. **Las dos fuentes son de pastillas de
tren de alta velocidad**; no existe todavía una fuente de moto (§7).

## 1. Qué se pidió

Pediste revisar `legacy/metals-14-01048.pdf` y, si era correcta, cambiar el preset y descartar la
fuente de trenes. Al revisarla encontré dos problemas: también es de trenes y no mide $c_p$. Te
los expliqué y elegiste la **opción híbrida**: $k$ y $\rho$ de Zhou, $c$ de Goo y Na.

## 2. Qué se hizo

- `configs/materiales.yaml`:
  - nueva entrada `sintered_pad_cu_zhou2024`;
  - el preset `motorcycle` ahora usa esa pastilla con el disco `stainless_steel_410`;
  - `sintered_pad_goo2014` sale del preset, pero se conserva en la librería porque es la fuente
    de $c$ de la nueva entrada (queda anotado).
- `tests/test_profile_regions.py`: nueva `test_sintered_pad_cu_zhou2024_entry`; en
  `test_motorcycle_preset` se actualizaron la pastilla esperada y $\gamma$.
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

**La fuente.** Zhou, Y.; Zhang, Y.; Zhang, X.; Liu, J.; Wang, M., "Effect of Sintering Temperature
on the Microstructure and Mechanical and Tribological Properties of Copper Matrix Composite for
Brake Pads", *Metals* 14:1048, 2024, DOI 10.3390/met14091048.

- **Material** (Tabla 1, % en peso): Cu 57, Fe 16, grafito 9, Cr₃Fe 6, Fe₃Al 6, Sn 4, MoS₂ 2.
  Es un compuesto sinterizado de matriz de cobre, sinterizado en caliente a 700–850 °C.
- **Mediciones** (sección 2.2): densidad por el método de Arquímedes y conductividad con un
  medidor TC3000E. **No hay medición de calor específico.**
- **Muestra elegida:** S800 (800 °C), la que el artículo destaca por mayor densidad,
  conductividad y desempeño en fricción.

**Lectura de la Figura 6.** Los valores solo aparecen en una gráfica rasterizada. Se renderizó la
página 7 a 800 dpi con PyMuPDF (instalado solo en una carpeta temporal) y se localizaron por
píxeles las marcas de los ejes y el centro del marcador S800. Con interpolación lineal entre dos
marcas de valores $y_1$, $y_2$ en las filas $f_1$, $f_2$, el marcador en la fila $f$ vale
$y = y_1 + (f_1 - f)\,\dfrac{y_1 - y_2}{f_2 - f_1}$:

- **$k$:** las marcas de 45 y 30 W/(m·K) están en las filas 191,0 y 453,0 (262,5 px por
  15 W/(m·K)). El marcador está en la fila 189,0:

  $$
  k = 45 + \frac{191{,}0 - 189{,}0}{262{,}5}\times 15 = 45{,}11 \approx \mathbf{45{,}1\ W/(m\,K)}.
  $$

  La barra de error del artículo es de ≈ ±3.
- **$\rho$:** las marcas de 5,538 y 5,467 g/cm³ están en las filas 94,0 y 350,0 (256 px por
  0,071). El marcador está en la fila 158,5:

  $$
  \rho = 5{,}538 - \frac{64{,}5}{256}\times 0{,}071 = \mathbf{5{,}520\ g/cm^3}.
  $$

  El eje dice "g·cm⁻¹", una errata evidente por g·cm⁻³.

**Efecto en $\gamma$.** La efusividad de la pastilla es

$$
e_p = \sqrt{45{,}1 \cdot 5520 \cdot 418} = 10\,201,
$$

mayor que la del disco 410 ($e_d = 9545$). Así,

$$
\gamma = \frac{e_d}{e_d + e_p} = \mathbf{0{,}4834}:
$$

menos de la mitad del calor de fricción entra al disco.

## 4. Decisiones

- **Una entrada híbrida y explícita.** El `source` cita las dos fuentes y dice qué propiedad sale
  de cada una y cómo se digitalizó. Así se puede auditar, y si aparece un $c_p$ medido solo hay
  que cambiar esa línea.
- **Muestra S800.** Es el material "óptimo" del propio artículo y el único que se parece a una
  pastilla comercial bien sinterizada. Las otras muestras tienen porosidades muy distintas: $k$
  va de ≈ 2 a 45 W/(m·K).
- **Se conserva `sintered_pad_goo2014` en la librería.** Pediste descartarla del preset, y así
  se hizo. No se borra porque la nueva entrada toma de ella $c$; borrarla dejaría ese valor sin
  rastro.
- **`verified: false`.** Hay dos fuentes, valores leídos de una gráfica, un $c$ no medido y una
  aplicación ferroviaria.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
```

Salida real:

```
236 passed in 18.29s
All checks passed!
```

Digitalización y $\gamma$ (salida real):

```
k: marcador S800 (np.float64(189.0), np.int64(175), np.int64(203)) ...
k: eje en x = 944 marcas [] [np.float64(68.5), np.float64(191.0), np.float64(453.0), np.float64(716.0)]
rho: marcador S800 (np.float64(158.5), np.int64(145), np.int64(172)) ...
rho: eje en x = 259 marcas [] [np.float64(94.0), np.float64(350.0)]
k 45.114285714285714 rho 5520.111328125
e_pad 10201.075237444335 gamma 410 + Zhou/Goo 0.48338848515817295
```

(68,5 es el borde superior del panel, no una marca.)

Comparación sobre `ejemplo_auto.yaml`, misma geometría y vehículo (salida real):

```
gray_cast_iron_limpert   + pad_limpert               gamma=0.8495 Tmax=106.70 en t=3.500  T60=75.88  E_in=86028 J
stainless_steel_410      + pad_limpert               gamma=0.8164 Tmax=119.44 en t=2.834  T60=70.65  E_in=82678 J
stainless_steel_410      + sintered_pad_goo2014      gamma=0.5646 Tmax=90.31 en t=2.834  T60=56.57  E_in=57179 J
stainless_steel_410      + sintered_pad_cu_zhou2024  gamma=0.4834 Tmax=80.92 en t=2.834  T60=52.03  E_in=48952 J
```

## 6. Resultados

| Disco + pastilla | $e$ pastilla | $\gamma$ | Energía al disco | $T_{\max}$ | $T(60\ \mathrm{s})$ |
|---|---|---|---|---|---|
| Hierro (Limpert) + Limpert | 2146 | 0,850 | 86,0 kJ | 106,7 °C | 75,9 °C |
| 410 + Limpert | 2146 | 0,816 | 82,7 kJ | 119,4 °C | 70,7 °C |
| 410 + Goo y Na (descartada del preset) | 7360 | 0,565 | 57,2 kJ | 90,3 °C | 56,6 °C |
| **410 + cobre Zhou/Goo (preset `motorcycle`)** | 10 201 | **0,483** | 49,0 kJ | **80,9 °C** | 52,0 °C |

Con una pastilla sinterizada de cobre de $k$ medida (45 W/(m·K)), la pastilla absorbe más de la
mitad del calor de fricción, según el modelo de efusividades. Esto concuerda con un rasgo
conocido de las pastillas sinterizadas: transmiten mucho calor hacia la pinza y el líquido de
frenos. Aquí se ve como un $\gamma$ bajo.

## 7. Limitaciones y riesgos

- **No hay datos de moto.** Las dos fuentes de la pastilla son de tren de alta velocidad. La
  composición es del mismo tipo que una pastilla sinterizada de moto, pero las proporciones
  reales de una SBS o Brembo pueden diferir.
- **$c$ no está medido.** Viene de un valor de entrada de un modelo (Goo y Na). Para un compuesto
  con 57 % de Cu, 16 % de Fe y 9 % de grafito, el orden de magnitud (~400–450 J/(kg·K)) es
  razonable, pero no está verificado.
- **Valores leídos de una gráfica.** La incertidumbre de lectura es menor que la barra de error
  experimental de $k$ (≈ ±3, unos ±7 %). Esa incertidumbre afecta $\gamma$ en ≈ ±1,5 %.
- **El modelo de $\gamma$ por efusividades supone sólidos semi-infinitos y contacto perfecto.**
  Con una pastilla tan conductora y delgada es la hipótesis más discutible del modelo de carga.
  Conviene mencionarlo en el informe.
- **El disco 410 tampoco está verificado** (reporte F7_T3).

## 8. Para revisar y para el informe

1. **De dónde sale cada número de la pastilla:** $k$ y $\rho$ medidos (Zhou, Fig. 6, S800) y $c$
   tomado de otro artículo (Goo y Na). Hay que saber decirlo y por qué.
2. **Cómo se digitalizó una gráfica** y cuál es su incertidumbre frente a la barra de error
   experimental.
3. **Por qué $\gamma \approx 0{,}48$.** La efusividad de la pastilla de cobre supera a la del
   acero 410, así que el reparto del calor favorece a la pastilla.
4. **Lo sensible que es el resultado a la pastilla:** $T_{\max}$ va de 80,9 a 119,4 °C con el
   mismo disco. Recomiendo un barrido de la pastilla en el informe.
5. **Para el informe:** la tabla de §6 (pares disco–pastilla), citando a Limpert, Chilexpo,
   Zhou et al. y Goo y Na, con las reservas de §7.
