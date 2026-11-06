# F7 · T7 — Caso de motocicleta con su geometría (Yamaha MT-07 de referencia)

Todas las pruebas pasaron (243 passed) sin relajar tolerancias. Sobre los datos del caso:

- **Con fuente:** los discos (2 × Ø 298 mm), la masa de la moto (≈ 184 kg, de fuentes
  secundarias) y el espesor (4,5 mm, de un disco de reemplazo Brembo).
- **Supuestos que aprobaste:** ancho del anillo, piloto, $\varphi$, desaceleración y $h$.
- **Materiales:** sin verificar (AISI 410 y pastilla sinterizada de cobre).

## 1. Qué se pidió

"Haz el caso de moto con su geometría". Era el pendiente de la Fase 7 (§5 de `PLAN.md`): los
barridos anteriores usaban la geometría y el vehículo del auto. Propuse la MT-07 con una tabla de
datos que separaba lo que tiene fuente de lo supuesto, y la aprobaste.

## 2. Qué se hizo

- `configs/ejemplo_moto.yaml`: nuevo. El caso de moto; cada dato está marcado como FUENTE o
  SUPUESTO en los comentarios.
- `configs/estudio_moto.yaml`: nuevo. Estudio de malla (rf 1, 2, 4), de $\Delta t$
  ($\Delta t$, $\Delta t/2$, $\Delta t/4$), barrido del ancho del anillo (el supuesto con más
  incertidumbre) y barrido de pastilla.
- `tests/test_moto.py`: nuevo. Geometría y materiales, energía y cota, y disco térmicamente
  delgado (3 pruebas).
- `README.md`: los dos archivos nuevos en la tabla de configuraciones.
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

No hubo cambios de código: el caso usa la sección rectangular y el preset `motorcycle`
existentes.

## 3. Cómo funciona

**Geometría.** Se modela solo el **anillo flotante** de acero:

- $r_o = 0{,}149$ m, de los discos de Ø 298 mm (FUENTE);
- anillo de 35 mm, así que $r_i = 0{,}114$ m (SUPUESTO);
- semiespesor $e = 2{,}25$ mm, de los 4,5 mm de espesor (FUENTE);
- pista de 30 mm centrada en el anillo, $r_a = 0{,}1165$ y $r_b = 0{,}1465$ m (SUPUESTO).

El porta-disco de aluminio se une al anillo por remaches. Su conducción se desprecia y el borde
interior del anillo tiene convección.

**Vehículo y carga.**

- $m = 184\ (\text{moto}) + 80\ (\text{piloto}) = 264$ kg; $v_0 = 100$ km/h;
  $a = 0{,}8\,g = 7{,}848\ \mathrm{m/s^2}$, así que $t_f = v_0/a = 3{,}54$ s;
- $\varphi = 0{,}8$ en el eje delantero; $n_d = 2$ discos con $n_c = 2$ caras.

La energía por cara es

$$
E_{\text{cara}} = \frac{\varphi\, m\, v_0^2}{2\, n_d\, n_c} = 20{,}4\ \mathrm{kJ}.
$$

Con el preset `motorcycle` ($\gamma = 0{,}483$), entran 9,85 kJ por cara.

**Disco térmicamente delgado.** La difusividad del 410 es
$\alpha = k/(\rho c) = 7{,}9\times10^{-6}\ \mathrm{m^2/s}$. El tiempo de difusión a través del
semiespesor,

$$
\tau = \frac{e^2}{\alpha} \approx 0{,}64\ \mathrm{s} \ll t_f = 3{,}54\ \mathrm{s},
$$

es mucho menor que el frenado. Por eso el calor atraviesa el espesor durante el frenado y la
temperatura es casi uniforme en $z$: en el pico, la superficie y el plano medio difieren en
0,5 °C. Es lo contrario del disco de auto (10 mm de semiespesor, $e^2/\alpha \approx 6$ s
$> t_f$), donde el gradiente en el espesor es grande. Por eso este caso es mucho menos exigente
con la malla.

**Mallado.** Se escalaron los tamaños al disco delgado: size_far = 1 mm y grading_distance =
3 mm, con rf = 2 (1641 nodos). El estudio de malla confirma que el resultado no depende de la
malla.

## 4. Decisiones

- **La MT-07** porque su ficha oficial da el número y el diámetro de los discos, y existe un disco
  de reemplazo con espesor publicado. El resto de los supuestos lo aprobaste, y cada uno está
  marcado en el YAML.
- **Sin radio de rueda.** Solo se usa para calcular $\omega$ y no afecta la temperatura;
  incluirlo habría sido otro supuesto innecesario.
- **$h = 40\ \mathrm{W/(m^2\,K)}$ y $T_0 = 25$ °C, como en el auto,** para que la comparación
  aísle geometría, vehículo y materiales. En una moto, el disco expuesto probablemente tiene un
  $h$ mayor; eso bajaría la temperatura final, no el pico (§6).
- **Barrido del ancho del anillo** porque es el supuesto geométrico con más incertidumbre.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m discfem run configs/ejemplo_moto.yaml
.venv\Scripts\python -m discfem study configs/estudio_moto.yaml
```

Salida real de pytest:

```
All checks passed!
243 passed in 20.85s
```

`run` (salida real, sin la línea de advertencia de materiales no verificados):

```
Caso 'ejemplo_moto_mt07_100kmh': 1641 nodos, 2967 triángulos, 1696 pasos (θ = 0.5, Rannacher = 2, Δt = 0.0354 s)
  γ = 0.4834   t_f = 3.54 s
  T_max = 76.727 en t = 3.257 s (r = 0.1313, z = 0.0022)
  T_max(t_end = 60 s) = 66.613
  energía [J, con 2π]: entra 9846.4, sale por convección 759.5, almacenada 9086.9; cierre del balance 6.6e-14
  sonda (r=0.1315, z=0.00225): T_max = 76.726 en t = 3.257 s; T_final = 66.607
  sonda (r=0.1315, z=0): T_max = 76.205 en t = 3.540 s; T_final = 66.608
  Materiales NO verificados: stainless_steel_410, sintered_pad_cu_zhou2024
  tiempo de cómputo: 0.40 s   salidas: results\ejemplo_moto
```

`study` (salida real, 9,5 s; se omiten las advertencias):

```
| rf=1 | 503 | 847 | 76.7218 | 3.25657 | — | — | — | — | 0.177081 |
| rf=2 | 1641 | 2967 | 76.7275 | 3.25657 | 0.00732929 | 0.0108716 | 0.00881481 | — | 0.369858 |
| rf=4 | 6038 | 11449 | 76.729 | 3.25657 | 0.00196248 | 0.00291093 | 0.00202125 | 1.90096 | 1.13974 |
Nivel seleccionado: rf=2 (ε < 1 %)
| Δt/1 | 0.0353976 | 1696 | 76.7275 | 3.25657 | 66.6131 | — | — | 0.366762 |
| Δt/2 | 0.0176988 | 3391 | 76.7363 | 3.23888 | 66.6194 | 0.0115696 | — | 0.621572 |
| Δt/4 | 0.00884939 | 6781 | 76.7384 | 3.23888 | 66.621 | 0.00265324 | 2.12447 | 1.12818 |
Barrido 'ancho_anillo': geometry.r_inner = [0.1155, 0.114, 0.109, 0.104]
| 0.1155 | 0.483388 | 76.8588 | 3.25657 | 68.5913 | ...
| 0.114 | 0.483388 | 76.7275 | 3.25657 | 66.6131 | ...
| 0.109 | 0.483388 | 76.692 | 3.22118 | 61.9612 | ...
| 0.104 | 0.483388 | 76.6919 | 3.22118 | 59.393 | ...
Barrido 'pastillas_moto': materials.pad = ['pad_limpert', 'sintered_pad_goo2014', 'sintered_pad_cu_zhou2024']
| pad_limpert | 0.816424 | 112.366 | 3.25657 | 95.2829 | 3.53976 | 16630.2 | ...
| sintered_pad_goo2014 | 0.56462 | 85.42 | 3.25657 | 73.606 | 3.53976 | 11501.1 | ...
| sintered_pad_cu_zhou2024 | 0.483388 | 76.7275 | 3.25657 | 66.6131 | 3.53976 | 9846.41 | ...
```

## 6. Resultados

**Caso de moto** (preset `motorcycle`): $T_{\max}$ = **76,7 °C** en $t = 3{,}26$ s (antes de
$t_f = 3{,}54$ s); $T(60\ \mathrm{s}) = 66{,}6$ °C; balance de $6{,}6\times10^{-14}$.

**Convergencia:**

- Malla: $\varepsilon$ = 0,007 % → 0,002 %, con orden observado $p = 1{,}90$ (≈ 2, el orden de
  P1).
- Tiempo: $p = 2{,}12$, el segundo orden de Crank–Nicolson + Rannacher.

El error de discretización es despreciable ($< 0{,}02$ °C).

**Sensibilidad a los supuestos:**

| Supuesto | Rango | $T_{\max}$ | $T(60\ \mathrm{s})$ |
|---|---|---|---|
| Ancho del anillo | 33,5 → 45 mm | 76,86 → 76,69 °C (±0,1) | 68,6 → 59,4 °C |
| Pastilla | Cu sinterizada → orgánica | 76,7 → 112,4 °C | 66,6 → 95,3 °C |

- El pico casi no depende del ancho del anillo: es local bajo la pista.
- La temperatura final sí depende, porque cambia la masa que se calienta.
- La pastilla domina, igual que en el auto.

**Auto frente a moto** (100 km/h, mismos $h$ y $T_0$):

| | Auto (`ejemplo_auto`) | Moto (`ejemplo_moto`) |
|---|---|---|
| Disco | hierro, 0,07–0,15 m × 20 mm | AISI 410, 0,114–0,149 m × 4,5 mm |
| Energía por cara (antes de $\gamma$) | 101,3 kJ | 20,4 kJ |
| $\gamma$ | 0,850 | 0,483 |
| $T_{\max}$ / $T(60\ \mathrm{s})$ | 106,7 / 75,9 °C | 76,7 / 66,6 °C |

La moto disipa 5 veces menos energía por cara en un disco con 8,5 veces menos volumen (por
cara). Además, la pastilla sinterizada se lleva más de la mitad del calor. Con la pastilla
orgánica de Limpert, la moto llegaría a 112 °C, más que el auto.

## 7. Limitaciones y riesgos

- **Supuestos sin fuente:** ancho del anillo, piloto, $\varphi$, desaceleración y $h$. El ancho
  del anillo tiene poco efecto en el pico; $\varphi$, la desaceleración y la masa escalan la
  energía de forma lineal (Fase 6), así que su efecto es predecible.
- **Materiales no verificados** (410 y pastilla sinterizada, reportes F7_T3 a F7_T5). Es la
  principal fuente de incertidumbre: con otra pastilla el pico cambia hasta 36 °C.
- **Se desprecia la conducción al porta-disco por los remaches.** Si fuera importante, bajaría
  la temperatura, sobre todo la final.
- **$h = 40$ es probablemente bajo** para un disco de moto expuesto al aire. Es un supuesto de
  comparación, no un valor medido.
- **Una sola frenada.** Las frenadas repetidas, donde el disco delgado acumula calor, no están
  modeladas (opcional en el PLAN).

## 8. Para revisar y para el informe

1. **Qué datos tienen fuente y cuáles son supuestos;** el YAML los marca uno por uno.
2. **El disco de moto es térmicamente delgado** ($e^2/\alpha \approx 0{,}6\ \mathrm{s} \ll t_f$):
   temperatura casi uniforme en el espesor y malla poco exigente. Es lo contrario del disco de
   auto.
3. **El pico no depende del ancho del anillo, pero la temperatura final sí** (efecto local frente
   a efecto de capacidad).
4. **La pastilla sigue siendo el factor dominante** (76,7 frente a 112,4 °C), y por eso pesa
   tanto que no esté verificada.
5. **Para el informe** (en `results/ejemplo_moto/` y `results/estudio_moto/`):
   - `field.png`, `tmax.png` y `mesh.png`;
   - las tablas de malla y $\Delta t$ como validación;
   - `sweep_ancho_anillo.pdf` y `sweep_pastillas_moto.pdf`;
   - la tabla auto frente a moto de §6.
