# F6 · T1 — Estudios de convergencia (malla y Δt) y barridos paramétricos

Todas las pruebas pasaron (194 passed) sin relajar tolerancias. El barrido de "materiales" es
una **sensibilidad** de k del disco (±20 % sobre el valor supuesto), no una comparación entre
materiales reales. No se inventaron materiales (regla 4); se explica en §4.

## 1. Qué se pidió

Fase 6 de `PLAN.md` (§5, "Fase 6 — Estudios de convergencia y paramétricos"):

- `convergence.py` y `python -m discfem study`.
- Estudio de malla con `refine_factor` = 1, 2, 4, 8 y el mismo Δt. Tabla con nodos, elementos,
  T_max, instante de T_max, ε_malla, orden observado p y tiempo de cómputo. Seleccionar el
  primer nivel con ε < 1 % también en las sondas.
- Estudio de Δt aparte (Δt, Δt/2, Δt/4).
- Barridos sobre a, v0, m y materiales, con tabla y figura comparativa.
- Figuras listas para el informe (PDF vectorial).

Criterios: el estudio de malla del caso de referencia cierra por debajo de 1 % desde el segundo
nivel (el prototipo bajó de 0,12 % a 0,01 %), y los barridos reproducen la tendencia esperada
(más desaceleración o más masa → más T_max).

## 2. Qué se hizo

- `src/discfem/convergence.py`: nuevo. Fórmulas (`relative_change_percent`, `observed_order`),
  `convergence_table`, `mesh_study`, `timestep_study`, `parameter_sweep`, tablas CSV y
  Markdown, figuras PDF + PNG, `load_study` y `run_study`.
- `configs/estudio_malla.yaml`: nuevo. Estudio del caso de referencia: malla 1/2/4/8, Δt 1/2/4
  y barridos de a, v0, m y k del disco.
- `src/discfem/timestepping.py`: `integrate(..., probe_matrix=P)` registra P·T en todos los
  pasos sin guardar los campos completos.
- `src/discfem/config.py`: nuevas claves `time.steps_per_tf` (para el estudio de Δt con
  `dt: auto`) y `output.save_every`.
- `src/discfem/solver.py` y `postprocess.py`: las sondas salen de `probe_history`, con todos
  los pasos aunque `save_every > 1`.
- `src/discfem/cli.py`: comando `study` (también con `--set`). El estado global de `warnings`
  se restaura al salir de `main`.
- `tests/test_convergence.py`: nuevo. Fórmulas, selección, criterio de la fase, orden en Δt,
  tendencias, salidas, validación y CLI (20 casos).
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

**Qué es un estudio de convergencia.** La solución discreta T_h depende del tamaño de malla h y
del paso Δt. Si el método converge, al refinar T_h se acerca a la solución exacta, y las
diferencias entre niveles sucesivos se achican. Como no se conoce la solución exacta del caso de
referencia, se comparan niveles consecutivos:

ε_i = |Q_i − Q_{i−1}| / |Q_i| × 100 %

Aquí Q es T_max, o el máximo en el tiempo de la temperatura en cada sonda, e i es el nivel de
refinamiento. Si el error se comporta como C·h^p, con tres niveles separados por un factor ρ
(ρ = 2 aquí) se puede estimar el exponente p sin conocer la solución exacta:

p = log(|Q_{i−2} − Q_{i−1}| / |Q_{i−1} − Q_i|) / log ρ

Con ρ = 2, esta es la fórmula del PLAN. Una prueba verifica que, para Q = Q* + C h^p sintética,
la función devuelve p a 1e-10.

**Estudio de malla.** Una corrida por `refine_factor`, con el mismo Δt (`dt: auto` = t_f/100 no
depende de la malla). Se selecciona el primer nivel i ≥ 1 con ε_i < 1 % en T_max **y** en todas
las sondas. Hay una prueba con datos sintéticos donde T_max ya converge pero una sonda no: se
elige el nivel siguiente.

**Estudio de Δt.** Mismo procedimiento con Δt, Δt/2, Δt/4, sobre la malla seleccionada. Con
`dt: auto` se varía `time.steps_per_tf` (100, 200, 400), así t_f sigue en la malla temporal.

**Barridos.** Una corrida por valor de una clave con puntos (por ejemplo,
`load.vehicle.decel`). Los barridos también usan la malla seleccionada.

**Memoria.** Para el nivel más fino (42 433 nodos × 1081 pasos), guardar todos los campos
costaría ~370 MB. Los estudios piden `output.save_every` enorme y leen T_max(t) y las sondas, que
`integrate` registra en cada paso con la matriz de interpolación P (Fase 5).

| Fórmula / requisito | Código |
|---|---|
| ε = \|Q_i − Q_{i−1}\| / \|Q_i\| × 100 | `relative_change_percent`, `convergence.py:42`; uso en `:108` |
| p = log2(\|T_h − T_h/2\| / \|T_h/2 − T_h/4\|) | `observed_order`, `convergence.py:47`; uso en `:117` |
| Primer nivel con ε < 1 % en T_max y sondas | `convergence_table`, `convergence.py:85-122` (`selected = i`, `:115`) |
| Malla 1, 2, 4, 8 con el mismo Δt | `mesh_study`, `convergence.py:125` |
| Δt, Δt/2, Δt/4 | `timestep_study`, `convergence.py:134` (`steps_per_tf`, `:140`) |
| Barridos | `parameter_sweep`, `convergence.py:148` |
| Tablas y figuras PDF | `write_table` (`:174`), `plot_convergence` (`:202`), `plot_sweep` (`:243`) |
| Sondas en todos los pasos sin guardar campos | `integrate(..., probe_matrix)`, `timestepping.py:165`, `:223` |

## 4. Decisiones

- **Sensibilidad a k en vez de "materiales".** La librería solo tiene un material de disco, y
  sin fuente. Inventar entradas (acero, carbono-cerámica…) violaría la regla 4. El barrido usa
  `materials.overrides.disc.k` = 40, 50 y 60 (±20 %): queda marcado como no verificado y su
  título lo dice. El código barre cualquier clave. Cuando haya materiales con fuente en la
  librería, basta con barrer `materials.regions.disc` con sus nombres; el YAML lo explica en
  un comentario.
- **ε sobre T, como define el PLAN, y además sobre T − T0.** Dividir entre T en °C hace que ε
  dependa de la escala elegida (con T0 = 25 °C, ε se reduce). Por eso se agrega la columna
  `eps_rise_percent`, más exigente (0,24 % frente a 0,18 % en el nivel 2). La selección usa la
  definición del PLAN.
- **Los estudios posteriores usan la malla seleccionada** (`use_selected_mesh: true`). Así el
  estudio de Δt mide el error temporal con una malla ya convergida. Se puede desactivar.
- **Figuras en PDF (vectorial) y PNG.** El PDF va al informe; el PNG sirve para mirarlas
  rápido.
- **Reproducibilidad.** La carpeta del estudio guarda `study.yaml`, `base_config.yaml`,
  `materiales.yaml` y `study_summary.json` con los overrides efectivos.
- **Efecto secundario corregido en la CLI.** `main()` cambiaba los filtros y el formato de
  `warnings` de forma global. Al llamarla dentro de pytest, las advertencias se filtraban a
  otras pruebas. Ahora todo ocurre dentro de `warnings.catch_warnings()`.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
.venv\Scripts\python -m discfem study configs/estudio_malla.yaml
```

Salida real de pytest y ruff:

```
All checks passed!
194 passed in 14.83s
```

Salida real del estudio completo. Se omiten las 24 líneas
`ADVERTENCIA (UnverifiedMaterialWarning)`, una por corrida; el número se contó con `grep -c`.

```
Estudio de malla: refine_factor = [1.0, 2.0, 4.0, 8.0]
| label | n_nodes | n_elements | T_max | t_of_T_max | eps_percent | eps_rise_percent | eps_probes_max_percent | order_p | compute_time_s |
|---|---|---|---|---|---|---|---|---|---|
| rf=1 | 929 | 1618 | 102.297 | 3.38916 | — | — | — | — | 0.200466 |
| rf=2 | 3028 | 5580 | 102.481 | 3.38916 | 0.179713 | 0.2377 | 0.183792 | — | 0.463059 |
| rf=4 | 11051 | 21148 | 102.549 | 3.38916 | 0.0665369 | 0.0879868 | 0.0667101 | 1.43251 | 2.1869 |
| rf=8 | 42433 | 82956 | 102.563 | 3.38916 | 0.0130987 | 0.0173206 | 0.0130111 | 2.34455 | 10.4893 |

Nivel seleccionado: rf=2 (ε < 1 %)
Estudio de Δt: factores [1, 2, 4]
| label | dt | n_steps | T_max | t_of_T_max | T_final_max | eps_percent | order_p | compute_time_s |
|---|---|---|---|---|---|---|---|---|
| Δt/1 | 0.05556 | 1080 | 102.481 | 3.38916 | 72.3351 | — | — | 0.448073 |
| Δt/2 | 0.02778 | 2160 | 102.705 | 3.41694 | 72.5746 | 0.218298 | — | 0.713316 |
| Δt/4 | 0.01389 | 4320 | 102.817 | 3.40305 | 72.6943 | 0.108934 | 1.00127 | 1.26124 |

Barrido 'desaceleracion': load.vehicle.decel = [3.0, 4.0, 5.0, 6.0, 7.0]
| value | T_max | t_of_T_max | T_final_max | t_f | brake_in_J | compute_time_s |
|---|---|---|---|---|---|---|
| 3 | 92.9131 | 6.7598 | 72.3945 | 9.26 | 87160.4 | 0.341426 |
| 4 | 97.4479 | 4.65315 | 72.3571 | 6.945 | 87160.4 | 0.397006 |
| 5 | 102.481 | 3.38916 | 72.3351 | 5.556 | 87160.4 | 0.456444 |
| 6 | 107.945 | 2.6391 | 72.3205 | 4.63 | 87160.4 | 0.496123 |
| 7 | 113.552 | 2.14303 | 72.3102 | 3.96857 | 87160.4 | 0.5599 |

Barrido 'velocidad_inicial': load.vehicle.v0 = [16.67, 22.22, 27.78, 33.33]
| value | T_max | t_of_T_max | T_final_max | t_f | brake_in_J | compute_time_s |
|---|---|---|---|---|---|---|
| 16.67 | 59.5357 | 1.73368 | 42.0322 | 3.334 | 31385.3 | 0.632434 |
| 22.22 | 78.9647 | 2.48864 | 55.2724 | 4.444 | 55762.6 | 0.51369 |
| 27.78 | 102.481 | 3.38916 | 72.3351 | 5.556 | 87160.4 | 0.447235 |
| 33.33 | 130.444 | 4.39956 | 93.1633 | 6.666 | 125466 | 0.398797 |

Barrido 'masa': load.vehicle.mass = [1000.0, 1250.0, 1500.0, 1750.0, 2000.0]
| value | T_max | t_of_T_max | T_final_max | t_f | brake_in_J | compute_time_s |
|---|---|---|---|---|---|---|
| 1000 | 76.654 | 3.38916 | 56.5567 | 5.556 | 58106.9 | 0.450274 |
| 1250 | 89.5675 | 3.38916 | 64.4459 | 5.556 | 72633.6 | 0.451205 |
| 1500 | 102.481 | 3.38916 | 72.3351 | 5.556 | 87160.4 | 0.444635 |
| 1750 | 115.394 | 3.38916 | 80.2243 | 5.556 | 101687 | 0.443769 |
| 2000 | 128.308 | 3.38916 | 88.1135 | 5.556 | 116214 | 0.44212 |

Barrido 'conductividad_disco': materials.overrides.disc.k = [40.0, 50.0, 60.0]
| value | T_max | t_of_T_max | T_final_max | t_f | brake_in_J | compute_time_s |
|---|---|---|---|---|---|---|
| 40 | 108.061 | 3.11136 | 72.0304 | 5.556 | 85835.2 | 0.448607 |
| 50 | 102.481 | 3.38916 | 72.3351 | 5.556 | 87160.4 | 0.440565 |
| 60 | 99.0921 | 3.66696 | 72.595 | 5.556 | 88165.1 | 0.43785 |

Salidas del estudio: results\estudio_malla

real	0m27.038s
```

## 6. Resultados

**Estudio de malla** (criterio: < 1 % desde el segundo nivel; prototipo: 0,12 % → 0,01 %):

| Nivel | Nodos | T_max [°C] | ε [%] | p |
|---|---|---|---|---|
| rf = 1 | 929 | 102,297 | — | — |
| rf = 2 | 3 028 | 102,481 | 0,180 | — |
| rf = 4 | 11 051 | 102,549 | 0,067 | 1,43 |
| rf = 8 | 42 433 | 102,563 | 0,013 | 2,34 |

- Cumple el criterio: ε < 1 % desde el nivel 2 en T_max y en las sondas, y se selecciona
  rf = 2. Los valores (0,18 → 0,067 → 0,013 %) son del mismo orden que los del prototipo, y el
  último coincide (0,01 %).
- El orden observado es irregular (1,43 y 2,34). Lo esperado para P1 en temperatura nodal es
  ≈ 2. Una causa probable (no verificada) es que las mallas de Gmsh no están anidadas: cada
  nivel es otra triangulación y el error no decrece de forma estrictamente geométrica.
- El instante de T_max es el mismo en todos los niveles (3,389 s) porque está cuantizado por
  Δt = 0,0556 s.

**Estudio de Δt** (malla rf = 2): el orden observado es **p = 1,001**, exactamente el primer
orden de Euler implícito (θ = 1, fijado en el YAML de referencia). La prueba
`test_timestep_study_first_order_for_implicit_euler` lo exige a ±0,05.

- **Con la configuración por defecto, el error temporal domina al espacial.** El salto de Δt a
  Δt/4 cambia T_max en 0,33 °C. El de rf = 2 a rf = 8 la cambia en 0,08 °C.
- **Extrapolación de Richardson** (con p = 1: T* ≈ 2·T_{Δt/4} − T_{Δt/2}): da T_max ≈ 102,93 °C y
  T(60 s) ≈ 72,81 °C. Este último coincide con el 72,80 °C que dio Crank–Nicolson + Rannacher
  en la Fase 5. Es una confirmación independiente de que la referencia de T8 (72,3 °C) lleva
  el sesgo de θ = 1 con Δt = t_f/100.

**Barridos.** Las tendencias esperadas se cumplen:

- **Desaceleración** (3 → 7 m/s²): T_max sube de 92,9 a 113,6 °C. La energía es la misma
  (87,2 kJ en las 5 corridas); solo se concentra en menos tiempo, porque t_f baja de 9,3 a
  4,0 s. La temperatura final no cambia (72,3–72,4 °C).
- **Masa** (1000 → 2000 kg): T_max sube de 76,7 a 128,3 °C. Es exactamente lineal en el
  incremento: (128,308 − 25)/(76,654 − 25) = 2,0000. Así debe ser, porque el problema es
  lineal y q ∝ m.
- **Velocidad** (60 → 120 km/h): T_max sube de 59,5 a 130,4 °C. La energía crece como v0²
  (31,4 → 125,5 kJ, factor 4,0).
- **k del disco** (40 → 60): T_max baja de 108,1 a 99,1 °C, porque un disco más conductor
  aleja el calor de la superficie. La energía de entrada cambia un poco (85,8 → 88,2 kJ)
  porque γ depende de la efusividad, que depende de k.

## 7. Limitaciones y riesgos

- **El barrido de materiales es solo de sensibilidad** y usa valores supuestos. Comparar
  materiales reales requiere entradas con fuente en `materiales.yaml`.
- **Orden espacial irregular** (1,43 / 2,34). No se probó con mallas estructuradas anidadas,
  que darían un p más limpio. Con la malla estructurada, ε da 0,19 → 0,059 → 0,017 % (corrida
  exploratoria, no incluida en el estudio).
- **El criterio de selección (ε < 1 %) es poco exigente para este caso.** Ya el nivel 1 está a
  ~0,26 % del más fino. La tolerancia es configurable (`tolerance_percent`).
- **Una advertencia de materiales por corrida:** 24 líneas en un estudio. Es ruidoso pero
  cumple la regla "cada corrida advierte".
- **El estudio completo tarda ~27 s.** Las pruebas usan versiones reducidas, salvo la de malla
  con rf = 1, 2, 4 en Gmsh (~3 s).

## 8. Para revisar y para el informe

1. **Cómo se estima el orden de convergencia sin solución exacta**, con tres niveles y la
   razón ρ, y por qué ε entre niveles consecutivos no es el error, sino una estimación de él.
2. **El estudio de Δt confirma p = 1 para θ = 1 y muestra que el error temporal domina** con
   Δt = t_f/100. Es un argumento concreto para usar θ = 1/2 + Rannacher, o Δt menor, en los
   resultados finales.
3. **La extrapolación de Richardson (72,81 °C) coincide con CN + Rannacher (72,80 °C).** Así se
   explica, con dos métodos independientes, la diferencia con la referencia de T8 (72,3 °C).
4. **Las tendencias de los barridos tienen explicación física:** misma energía en menos tiempo
   (a), linealidad exacta (m), energía ∝ v0² (v0), conducción contra acumulación superficial
   (k).
5. **Para el informe** (todo en `results/estudio_malla/`):
   - `mesh_study.pdf` y la tabla `mesh_study.md`: objetivo de convergencia de malla.
   - `timestep_study.pdf`: validación del orden temporal.
   - `sweep_desaceleracion.pdf`, `sweep_masa.pdf`, `sweep_velocidad_inicial.pdf`: estudio
     paramétrico.
   - `sweep_conductividad_disco.pdf`, rotulada como sensibilidad con valores supuestos.
