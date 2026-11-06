# F7 · T6 — Barrido de pastilla para el informe

Todas las pruebas pasaron (240 passed) sin relajar tolerancias. Es una comparación de
**materiales** sobre la geometría y el vehículo del caso auto (100 km/h), no un caso de
motocicleta. Las pastillas sinterizadas y el disco 410 siguen sin verificar (reportes F7_T3 a
F7_T5).

## 1. Qué se pidió

"Haz el barrido de pastilla para el informe", como recomendé en F7_T4 y F7_T5. La temperatura
del disco depende fuertemente de la pastilla elegida, y las propiedades de las pastillas
sinterizadas no están verificadas.

## 2. Qué se hizo

- `configs/estudio_pastilla.yaml`: nuevo. Cuatro barridos y una figura resumen:
  - las tres pastillas del catálogo con el disco de hierro (Limpert) y con el disco AISI 410;
  - $\gamma$ continuo (0,40 a 0,90) con cada disco;
  - `gamma_map`: curvas $T(\gamma)$ con las pastillas reales ubicadas encima.
- `src/discfem/convergence.py`:
  - `overrides` propios de cada barrido (aquí, el disco);
  - columna `gamma` en las tablas de barrido;
  - sección `gamma_map` y función `plot_gamma_map` (PDF + PNG), con etiquetas legibles
    (`point_labels`);
  - validaciones en español.
- `tests/test_convergence.py`: 4 pruebas nuevas, incluida la linealidad exacta en $\gamma$.
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

**Cómo entra la pastilla al modelo.** La pastilla no se malla. Solo determina qué fracción del
calor de fricción entra al disco:

$$
\gamma = \frac{e_{\text{disco}}}{e_{\text{disco}} + e_{\text{pastilla}}}, \qquad e = \sqrt{k \rho c},
$$

y el flujo en la pista es $q = \gamma \cdot P_f(t) \cdot \dfrac{g(r)\, r}{2\pi I_g}$ (Fase 4). Por
eso basta barrer las pastillas del catálogo para cada disco.

**Una propiedad exacta del modelo.** En este caso $T_0 = T_\infty = 25$ °C. Entonces el problema
para el aumento de temperatura $u = T - T_0$ es lineal y homogéneo salvo por la carga, que es
proporcional a $\gamma$. En consecuencia:

$$
u(r, z, t;\, \gamma) = \gamma \cdot u(r, z, t;\, \gamma = 1)
$$

Es decir, **$T - T_0$ es exactamente proporcional a $\gamma$** para un disco dado, y el instante
del pico no depende de $\gamma$. Las curvas de la figura son rectas que pasan por
$(\gamma = 0,\ T = T_0)$, y las pastillas reales caen exactamente sobre ellas. La prueba
`test_pad_sweep_overrides_gamma_map_and_linearity` lo verifica a $10^{-10}$, y en el estudio
completo la dispersión de $(T - T_0)/\gamma$ es $\sim 10^{-14}$.

Consecuencia práctica: para un disco, una sola simulación basta para predecir el resultado con
cualquier pastilla, porque

$$
T - T_0 = \gamma\, S_{\text{disco}},
$$

donde $S_{\text{disco}}$ es la pendiente de su recta.

| Elemento | Código |
|---|---|
| Overrides por barrido | `run_study`, `convergence.py` (`{**ov, **sw["overrides"]}`) |
| Columna $\gamma$ en las tablas | `parameter_sweep` (`row["gamma"]`), `SWEEP_COLS` |
| Figura $T(\gamma)$ con pastillas | `plot_gamma_map`, `convergence.py` |
| $\gamma$ por efusividades | `partition_coefficient`, `materials.py:230` (Fase 4) |

## 4. Decisiones

- **Barrer $\gamma$ en vez de inventar pastillas intermedias.** $\gamma$ es un parámetro del
  modelo, no una propiedad de material, así que el barrido continuo no inventa datos (regla 4).
  Las únicas pastillas con nombre son las del catálogo, con su fuente.
- **Base en `ejemplo_auto.yaml`** (rf = 2, $\theta = 1/2$ + 2 Rannacher). Es el único caso con
  geometría y vehículo completos; el título de la figura lo dice.
- **Extender `study` en vez de escribir un script aparte.** Así el barrido se reproduce con un
  solo comando y queda probado.
- **Etiquetas de la figura.** Las del primer disco van abajo a la derecha y las del segundo
  arriba a la izquierda, con nombres legibles definidos en el YAML.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
.venv\Scripts\python -m discfem study configs/estudio_pastilla.yaml
```

Salida real:

```
All checks passed!
240 passed in 21.03s
```

Estudio (salida real de la primera corrida; se omiten las 11 líneas `ADVERTENCIA` de
materiales no verificados. La segunda corrida solo cambió las etiquetas de la figura):

```
Barrido 'pastillas_hierro': materials.pad = ['pad_limpert', 'sintered_pad_goo2014', 'sintered_pad_cu_zhou2024']
| value | gamma | T_max | t_of_T_max | T_final_max | t_f | brake_in_J | compute_time_s |
|---|---|---|---|---|---|---|---|
| pad_limpert | 0.849504 | 106.695 | 3.50028 | 75.8795 | 5.556 | 86028.5 | 0.518948 |
| sintered_pad_goo2014 | 0.62207 | 84.8234 | 3.50028 | 62.2577 | 5.556 | 62996.5 | 0.48187 |
| sintered_pad_cu_zhou2024 | 0.54288 | 77.2078 | 3.50028 | 57.5148 | 5.556 | 54976.9 | 0.49634 |

Barrido 'pastillas_410': materials.pad = ['pad_limpert', 'sintered_pad_goo2014', 'sintered_pad_cu_zhou2024']
| value | gamma | T_max | t_of_T_max | T_final_max | t_f | brake_in_J | compute_time_s |
|---|---|---|---|---|---|---|---|
| pad_limpert | 0.816424 | 119.438 | 2.83356 | 70.6546 | 5.556 | 82678.5 | 0.471997 |
| sintered_pad_goo2014 | 0.56462 | 90.3115 | 2.83356 | 56.5736 | 5.556 | 57178.5 | 0.465929 |
| sintered_pad_cu_zhou2024 | 0.483388 | 80.9152 | 2.83356 | 52.0312 | 5.556 | 48952.3 | 0.472145 |

Barrido 'gamma_hierro': load.gamma = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
| 0.4 | 0.4 | 63.4673 | 3.50028 | 48.9573 | ...
| 0.9 | 0.9 | 111.551 | 3.50028 | 78.9038 | ...

Barrido 'gamma_410': load.gamma = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
| 0.4 | 0.4 | 71.2693 | 2.83356 | 47.3681 | ...
| 0.9 | 0.9 | 129.106 | 2.83356 | 75.3282 | ...

Figura γ: results\estudio_pastilla\gamma_map.pdf
real	0m11.564s
```

(Las tablas de $\gamma$ completas están en `results/estudio_pastilla/sweep_gamma_*.md`.)

Linealidad en $\gamma$ (salida real):

```
hierro  T_max        (T-T0)/γ: min 96.168281048068  max 96.168281048069  dispersión relativa 3.69e-15
hierro  T_final_max  (T-T0)/γ: min 59.893143187581  max 59.893143187582  dispersión relativa 8.30e-15
410     T_max        (T-T0)/γ: min 115.673326338197  max 115.673326338198  dispersión relativa 4.30e-15
410     T_final_max  (T-T0)/γ: min 55.920171166890  max 55.920171166892  dispersión relativa 2.27e-14
```

## 6. Resultados

**Pendientes por disco** ($T - T_0 = \gamma\, S$):

| Disco | $S_{\text{pico}}$ [°C] | $S_{\text{final}}$ (60 s) [°C] | Instante del pico |
|---|---|---|---|
| Hierro fundido (Limpert) | 96,17 | 59,89 | 3,50 s |
| AISI 410 | 115,67 | 55,92 | 2,83 s |

**Pastillas reales:**

| Pastilla | $\gamma$ hierro | $T_{\max}$ hierro | $\gamma$ 410 | $T_{\max}$ 410 |
|---|---|---|---|---|
| Orgánica (Limpert) | 0,850 | 106,7 °C | 0,816 | 119,4 °C |
| Sinterizada (Goo y Na) | 0,622 | 84,8 °C | 0,565 | 90,3 °C |
| Sinterizada Cu (Zhou) | 0,543 | 77,2 °C | 0,483 | 80,9 °C |

Lectura para el informe:

1. **La pastilla mueve más la temperatura que el disco.**
   - Con el disco 410, el pico va de 80,9 a 119,4 °C según la pastilla (38,5 °C de rango).
   - Con el hierro, va de 77,2 a 106,7 °C (29,5 °C).
   - Cambiar el disco con la misma pastilla mueve el pico 3,7 °C (Zhou) o 12,7 °C (Limpert).
2. **A igual $\gamma$, el 410 tiene un pico mayor** (pendiente 115,7 frente a 96,2), porque su
   menor conductividad retiene el calor en la superficie. Pero la temperatura final es menor
   (55,9 frente a 59,9), porque tiene más capacidad térmica volumétrica. Con pastillas reales,
   el 410 además recibe un $\gamma$ menor, porque su efusividad es más baja que la del hierro.
3. **Las pastillas sinterizadas reducen $\gamma$ a 0,48–0,62.** El calor que no entra al disco
   va a la pastilla, y en la práctica a la pinza y al líquido de frenos. El modelo no sigue ese
   calor.
4. **Incertidumbre:** con la barra de error de $k$ de Zhou (±3 W/(m·K)), $\gamma$ varía
   ≈ ±1,5 %, y por la linealidad $T - T_0$ varía lo mismo (≈ ±0,8 °C en el pico con el 410).

## 7. Limitaciones y riesgos

- **No es un caso de moto.** Se usan la geometría y el vehículo del auto, y con un disco de moto
  más delgado los valores absolutos cambiarían. Las conclusiones relativas (rangos y orden de
  pastillas) sí son las que aporta el estudio.
- **Materiales no verificados.** El disco 410 y las dos pastillas sinterizadas lo están; la
  pastilla de Limpert y el hierro sí están verificados.
- **$\gamma$ por efusividades** supone contacto perfecto entre cuerpos semi-infinitos. Con
  pastillas muy conductoras es la hipótesis más discutible del modelo.
- **La linealidad exacta vale porque $T_0 = T_\infty$.** Con $T_0 \ne T_\infty$,
  $T - T_0 = \gamma\, S$ + (respuesta a la convección con $T_0 \ne T_\infty$): las rectas ya no
  pasarían por $(0, T_0)$, aunque seguirían siendo rectas.
- **Etiquetas de la figura:** en el panel derecho, una etiqueta toca un punto del otro disco
  cerca de $\gamma \approx 0{,}55$. Se lee, pero se puede pulir.

## 8. Para revisar y para el informe

1. **La pastilla solo entra por $\gamma$, y $T - T_0$ es proporcional a $\gamma$.** Es la
   explicación de por qué las curvas son rectas y por qué una simulación por disco basta.
2. **La pastilla pesa más que el disco** en el pico térmico, con los números de §6.
3. **Pico frente a final:** la conductividad controla el pico (más alto con el 410) y la
   capacidad térmica controla el final (más bajo con el 410).
4. **Cómo citar los datos:** Limpert (verificado), Chilexpo, Goo y Na, y Zhou et al. (estos
   tres sin verificar), con las reservas de los reportes F7_T3 a F7_T5.
5. **Figuras para el informe** (en `results/estudio_pastilla/`):
   - `gamma_map.pdf`: la figura principal;
   - `sweep_pastillas_hierro.pdf` y `sweep_pastillas_410.pdf`: $T_{\max}(t)$ por pastilla;
   - las tablas `sweep_*.md`.
