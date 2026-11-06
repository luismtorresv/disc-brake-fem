# F6 · T2 — Estudio de malla y barridos sobre el caso principal, y sensibilidad a h

Todas las pruebas pasaron (243 passed) sin relajar tolerancias. No se cambió ningún dato, como
pediste: $h$ sigue en $40\ \mathrm{W/(m^2\,K)}$ y el barrido es solo de sensibilidad. Los
pendientes quedaron registrados en la nueva sección "PENDIENTES" de `PLAN.md`.

## 1. Qué se pidió

Tres cosas:

- "hacer lo de $h$": cuánto depende el resultado del coeficiente de convección, que no tiene
  fuente;
- repetir el estudio de malla y los barridos de la Fase 6 (§5 de `PLAN.md`) sobre el caso
  principal, porque la Fase 6 se hizo sobre el caso de referencia (materiales supuestos,
  $\theta = 1$);
- mantener los datos como están y dejar los pendientes registrados.

## 2. Qué se hizo

- `configs/estudio_auto.yaml`: nuevo. Sobre `ejemplo_auto.yaml`: malla (rf 1, 2, 4, 8),
  $\Delta t$ ($\Delta t$, $\Delta t/2$, $\Delta t/4$) y barridos de $a$, $v_0$, $m$, $k$ del disco
  (±20 %) y $h$ (10–160 $\mathrm{W/(m^2\,K)}$).
- `configs/estudio_moto.yaml`: se añadió el barrido de $h$.
- `PLAN.md`: nueva sección **PENDIENTES** (12 puntos) y nota en ESTADO de la Fase 6.
- `README.md`: `estudio_auto.yaml` en la tabla de configuraciones, marcado como el estudio para
  el informe.
- `configs/estudio_malla.yaml` **no se modificó**: sigue reproduciendo la Fase 6 original.
- No hubo cambios de código.

## 3. Cómo funciona

Se usa la misma maquinaria de la Fase 6 (`convergence.py`). La única novedad es la elección de la
base y del barrido de $h$. $h$ entra en la condición de Robin

$$
-k\, \nabla T \cdot \mathbf{n} = h\,(T - T_\infty)
$$

de las caras expuestas (`top_inner`, `top_outer`, `r_inner`, `r_outer`) y en las matrices
$\mathbf{H}$ y $\mathbf{F}_{\text{conv}}$ (Fase 2).

**Por qué $h$ no afecta el pico.** La potencia que sale por convección es
$h\, A\, (T - T_\infty)$. Con $h = 40$, un área expuesta de unos centímetros cuadrados por radián
y $T - T_\infty$ de decenas de grados, en los 3–5 s del frenado sale una fracción ínfima de la
energía que entra. En el caso auto entran 86 kJ por cara y salen 2,6 kJ en 60 s, casi todo
después del pico. Así, el pico lo controlan la entrada de calor y la conducción, y $h$ solo actúa
en el enfriamiento posterior.

## 4. Decisiones

- **Un estudio nuevo en vez de modificar `estudio_malla.yaml`.** El original documenta la Fase 6
  y su reporte. Cambiarlo rompería la trazabilidad de F6_T1.
- **Rango de $h$ de 10 a 160 $\mathrm{W/(m^2\,K)}$**, un factor 4 a cada lado del valor supuesto
  y dentro del rango de sanidad (1–500). Cubre desde convección natural hasta una convección
  forzada moderada.
- **No se adoptó un $h$ nuevo.** Pediste mantener los datos. La correlación con la velocidad del
  vehículo queda como pendiente 3.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m discfem study configs/estudio_auto.yaml
.venv\Scripts\python -m discfem study configs/estudio_moto.yaml
```

Salida real de `estudio_auto` (30,4 s; se omiten las advertencias):

```
| rf=1 | 929 | 1618 | 106.51 | 3.50028 | — | — | — | — | 0.210614 |
| rf=2 | 3028 | 5580 | 106.695 | 3.50028 | 0.17372 | 0.226882 | 0.177304 | — | 0.463273 |
| rf=4 | 11051 | 21148 | 106.764 | 3.50028 | 0.064262 | 0.0839107 | 0.0644471 | 1.4338 | 2.03234 |
| rf=8 | 42433 | 82956 | 106.777 | 3.50028 | 0.0126571 | 0.0165265 | 0.0125648 | 2.34383 | 10.8667 |
Nivel seleccionado: rf=2 (ε < 1 %)
| Δt/1 | 0.05556 | 1080 | 106.695 | 3.50028 | 75.8795 | — | — | 0.459313 |
| Δt/2 | 0.02778 | 2160 | 106.706 | 3.50028 | 75.8871 | 0.0102878 | — | 0.751493 |
| Δt/4 | 0.01389 | 4320 | 106.709 | 3.50028 | 75.889 | 0.00257826 | 1.99643 | 1.30213 |
desaceleración 3/4/5/6/7 m/s²: T_max 97.3198 / 101.788 / 106.695 / 112.12 / 117.806; brake_in_J 86028.5 (todas)
velocidad 16.67/22.22/27.78/33.33 m/s: T_max 61.1444 / 81.6419 / 106.695 / 136.66; brake_in_J 30977.7 / 55038.4 / 86028.5 / 123836
masa 1000/1250/1500/1750/2000 kg: T_max 79.4635 / 93.0794 / 106.695 / 120.311 / 133.927
k disco 38.77/48.4625/58.155: γ 0.834677 / 0.849504 / 0.860791; T_max 111.983 / 106.695 / 103.597
```

Barrido de $h$ con todos los decimales, leído de los CSV (salida real):

```
auto
  h=   10  T_max=106.695343  T60=77.0111  ΔT60=52.011  vs h=40: +2.2 %
  h=   20  T_max=106.695336  T60=76.5851  ΔT60=51.585  vs h=40: +1.4 %
  h=   40  T_max=106.695324  T60=75.8795  ΔT60=50.879  vs h=40: +0.0 %
  h=   80  T_max=106.695298  T60=74.6747  ΔT60=49.675  vs h=40: -2.4 %
  h=  160  T_max=106.695248  T60=72.5703  ΔT60=47.570  vs h=40: -6.5 %
moto
  h=   10  T_max=76.728450  T60=68.8868  ΔT60=43.887  vs h=40: +5.5 %
  h=   20  T_max=76.728119  T60=68.1095  ΔT60=43.109  vs h=40: +3.6 %
  h=   40  T_max=76.727460  T60=66.6131  ΔT60=41.613  vs h=40: +0.0 %
  h=   80  T_max=76.726148  T60=63.8356  ΔT60=38.836  vs h=40: -6.7 %
  h=  160  T_max=76.723556  T60=59.0290  ΔT60=34.029  vs h=40: -18.2 %
```

## 6. Resultados

**Convergencia del caso principal:**

- Malla: $\varepsilon$ = 0,17 → 0,064 → 0,013 %, igual que en la Fase 6. El nivel seleccionado es
  rf = 2, el que ya usa `ejemplo_auto.yaml`.
- Tiempo: con $\theta = 1/2$ + 2 Rannacher, el orden es **$p = 2{,}00$** y el cambio de
  $\Delta t$ a $\Delta t/4$ es de 0,013 °C. En la Fase 6, con $\theta = 1$ sobre el caso de
  referencia, era $p = 1{,}00$ y 0,34 °C. El error temporal, que antes dominaba, ahora es
  despreciable.
- Error total estimado del pico: $< 0{,}1$ °C (malla) + $< 0{,}02$ °C (tiempo).

**Barridos:** las mismas tendencias de la Fase 6, con los valores de los materiales de Limpert.

- Desaceleración (3 → 7 m/s²): pico de 97,3 a 117,8 °C, con la misma energía y una temperatura
  final casi constante (75,9 °C).
- Masa: lineal exacta,

  $$
  \frac{133{,}927 - 25}{79{,}4635 - 25} = 2{,}0000 .
  $$

- Velocidad: energía $\propto v_0^2$ y pico de 61,1 a 136,7 °C.
- $k$ del disco (±20 %): pico de 112,0 a 103,6 °C. $\gamma$ también varía (0,835–0,861), porque
  depende de la efusividad del disco.

**$h$:**

| | Pico ($h$ = 10 → 160) | $\Delta T$ a 60 s, respecto de $h = 40$ |
|---|---|---|
| Auto | 106,6953 °C, cambia $< 0{,}0001$ °C | +2,2 % … −6,5 % |
| Moto | 76,728 → 76,724 °C | +5,5 % … −18,2 % |

**$h$ no influye en el pico, que es el resultado principal.** Sí influye en el enfriamiento, y más
en la moto: su disco delgado tiene 2,4 veces más área con convección por unidad de volumen
($A/V = 121$ frente a $50\ \mathrm{m^{-1}}$), y la convección actúa sobre un $\Delta T$ casi
uniforme.

## 7. Limitaciones y riesgos

- **$h$ sigue sin fuente** (pendiente 3). Para la temperatura final de la moto, la incertidumbre
  por $h$ es la mayor después de la de la pastilla.
- **Mismo $h$ en todas las caras.** En un disco real, el borde exterior y las caras ventiladas
  por el giro tienen $h$ distintos; `convection.per_group` permite modelarlo cuando haya datos.
- **El barrido de $k$ es de sensibilidad** (override marcado como no verificado), no de
  materiales reales.
- **Los materiales de moto siguen sin verificar** (pendientes 1 y 2).

## 8. Para revisar y para el informe

1. **Qué gobierna el pico y qué el enfriamiento:** el pico depende de la entrada de calor y la
   conducción ($a$, $v_0$, $m$, $k$, pastilla); el final depende de la energía total, la
   capacidad y la convección ($h$). Por eso $h$ no afecta el pico.
2. **Con el esquema $\theta = 1/2$ + Rannacher el error temporal pasa a ser despreciable**
   ($p = 2$). En la Fase 6, con $\theta = 1$, dominaba. Es un argumento para la elección del
   esquema.
3. **La convergencia de malla del caso principal** está documentada con el mismo criterio de la
   Fase 6.
4. **La moto es más sensible a $h$** porque su disco es delgado (más área por volumen).
5. **Para el informe** (en `results/estudio_auto/`, y `results/estudio_moto/` para la moto):
   - `mesh_study.pdf` y `timestep_study.pdf` (validación del caso principal);
   - `sweep_desaceleracion.pdf`, `sweep_masa.pdf` y `sweep_velocidad_inicial.pdf` (estudio
     paramétrico);
   - `sweep_conveccion_h.pdf` de ambos casos (sensibilidad a un dato sin fuente);
   - la tabla de $h$ de §6.
