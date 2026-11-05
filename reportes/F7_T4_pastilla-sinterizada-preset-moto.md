# F7 · T4 — Pastilla sinterizada (Goo y Na 2014) en el preset `motorcycle`

Todas las pruebas pasaron (235 passed) sin relajar tolerancias. La pastilla quedó con
**`verified: false`**: sus valores son datos de entrada de un modelo FEM ferroviario, no
mediciones (§3). La prueba del preset `motorcycle` cambió porque el preset ahora trae pastilla.
El comportamiento "preset sin pastilla ⇒ error claro" se sigue probando con una librería de
prueba (§4).

## 1. Qué se pidió

Continuación de la Fase 7 (§5 de `PLAN.md`):

- Dejaste el **ASM Metals Handbook, Vol. 7** en `legacy/` porque no encontraste datos oficiales
  de SBS ni de Brembo.
- Al revisarlo, el ASM no trae propiedades térmicas de pastillas sinterizadas (§3).
- Elegiste la **opción A**: cargar la pastilla de Goo y Na (2014) como `sintered_pad_goo2014`
  con `verified: false` y, siendo la más cercana a la realidad disponible, usarla en el preset
  `motorcycle`.

## 2. Qué se hizo

- `configs/materiales.yaml`:
  - nueva entrada `sintered_pad_goo2014`, con nota en la cabecera y las reservas en `source`;
  - el preset `motorcycle` = `stainless_steel_410` + `sintered_pad_goo2014`.
- `tests/test_profile_regions.py`:
  - `test_sintered_pad_entry_goo2014` (valores y fuente);
  - `test_motorcycle_preset` (disco, pastilla, advertencias y γ);
  - `test_preset_without_pad_needs_explicit_pad`, que reemplaza a la prueba anterior del preset
    de moto sin pastilla.
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

**Qué hay en el ASM Vol. 7.** Se extrajo el texto completo con `pypdf` y se revisó:

- el artículo "Friction Powder Metallurgy Materials" (pp. 2609–2616 del PDF);
- la sección "Copper-Base Friction Materials", Tabla 14 (p. 2175);
- las 124 apariciones de "thermal conductivity", "specific heat" o "diffusivity".

Hay composiciones típicas (base cobre: Cu 65–80 %, Sn 7–9 %, Fe 4–7 %, Pb 5–10 %, grafito 3–8 %,
SiO₂ 2–4 %) y coeficientes de fricción. **No hay ρ, c_p ni k del material de fricción.** El dato
más cercano es k del bronce macizo con 5 % de Sn, ≈ 71 W/(m·K) (Tabla 4, p. 2091), pero es la
matriz sin grafito, sin cerámicos y sin porosidad, así que no representa la pastilla.

**La fuente usada.** Goo, B.-C. y Na, I.-K., "Topology Optimization of Railway Brake Pad by
Contact Analysis", *J. Korean Soc. Tribol. Lubr. Eng.* 30(3):177–182, 2014,
DOI 10.9725/kstle.2014.30.3.177. Tabla 1, columna "Friction pin" (pines sinterizados de la
pastilla del KTX):

| Propiedad | Valor | Cargado |
|---|---|---|
| Densidad | 5,4 g/cm³ | 5400 kg/m³ |
| Calor específico | 418 J/(kg·°C) | 418 J/(kg·K) |
| Conductividad | 24 W/(m·°C) | 24 W/(m·K) |

**Por qué γ cambia tanto.** La fracción de calor que entra al disco es
γ = e_d / (e_d + e_p), con la efusividad e = √(k ρ c). Para la pastilla de Limpert,
e_p = 2146. Para la sinterizada, e_p = √(24 · 5400 · 418) = 7360, del mismo orden que la del
disco 410 (e_d = 9545). Así, γ pasa de 0,816 a **0,565**: la pastilla metálica absorbe el 44 %
del calor de fricción.

## 4. Decisiones

- **`verified: false`.** El artículo usa esos valores como entrada de su modelo en ABAQUS: no
  dice que los midiera ni cita de dónde salen. Además, su c = 418 es idéntico al del disco de
  acero de la misma tabla, lo que sugiere que lo copiaron, y la aplicación es ferroviaria. Todo
  queda escrito en `source` y en la cabecera de `materiales.yaml`, y cada corrida lo advierte.
- **Uso en el preset `motorcycle`.** Lo autorizaste como la aproximación más cercana disponible
  a una pastilla sinterizada de base metálica.
- **Cambio de prueba (no se relajó nada).** `test_motorcycle_preset_requires_explicit_pad`
  probaba que el preset de moto sin pastilla daba error. Como ahora el preset trae pastilla, esa
  premisa ya no vale. El comportamiento (un preset sin pastilla exige `materials.pad`) se sigue
  probando con una librería temporal en `test_preset_without_pad_needs_explicit_pad`.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
```

Salida real:

```
235 passed in 19.20s
All checks passed!
```

Efusividades y γ (salida real):

```
e410 9545.049160436416 e_goo 7360.217388093914 gamma 410+goo 0.5646198557730651 gamma 410+limpert 0.8164237361701493
```

Comparación sobre `ejemplo_auto.yaml`, misma geometría y vehículo, cambiando disco y pastilla
(salida real):

```
gray_cast_iron_limpert   + pad_limpert           gamma=0.8495 Tmax=106.70 en t=3.500  T60=75.88  E_in=86028 J
stainless_steel_410      + pad_limpert           gamma=0.8164 Tmax=119.44 en t=2.834  T60=70.65  E_in=82678 J
stainless_steel_410      + sintered_pad_goo2014  gamma=0.5646 Tmax=90.31 en t=2.834  T60=56.57  E_in=57179 J
```

## 6. Resultados

| Disco + pastilla | γ | Energía al disco | T_max | T(60 s) |
|---|---|---|---|---|
| Hierro (Limpert) + pastilla Limpert | 0,850 | 86,0 kJ | 106,7 °C | 75,9 °C |
| AISI 410 + pastilla Limpert | 0,816 | 82,7 kJ | 119,4 °C | 70,7 °C |
| AISI 410 + sinterizada (preset `motorcycle`) | 0,565 | 57,2 kJ | 90,3 °C | 56,6 °C |

La elección de la pastilla pesa más que la del disco: pasar a la sinterizada quita un 31 % de la
energía que entra al disco. Por eso el resultado de la moto es muy sensible a una propiedad no
verificada. Esto debe decirse en el informe.

Advertencia de interpretación: son comparaciones **de materiales** sobre la geometría y el
vehículo del auto, no un caso de moto real (§7).

## 7. Limitaciones y riesgos

- **Ni la pastilla ni el disco de moto están verificados.** La pastilla viene de valores de
  entrada de un FEM ferroviario; el disco, de una ficha con un error de unidades. Cualquier
  resultado de moto debe presentarse con esa reserva.
- **γ por efusividades supone contacto perfecto entre dos sólidos semi-infinitos.** Con una
  pastilla muy conductora y delgada (la sinterizada de moto mide ~4–5 mm más la placa), esa
  hipótesis es más discutible: la pastilla se "satura" antes y γ efectivo aumentaría con el
  tiempo. El modelo no lo captura (CLAUDE.md: la resistencia de contacto está fuera de alcance).
- **Recomendación para el informe:** un barrido de la pastilla entre Limpert (k ≈ 1,2) y
  Goo y Na (k = 24) para mostrar la sensibilidad. Es la opción B que propuse; no está hecha.
- **Sigue faltando un caso de moto con su geometría y vehículo** (disco delgado, masa, reparto).

## 8. Para revisar y para el informe

1. **Por qué una pastilla conductora reduce el calor que entra al disco:** γ depende de la
   razón de efusividades, y la sinterizada (e ≈ 7360) se acerca a la del acero (e ≈ 9545).
2. **Por qué la pastilla pesa más que el disco** en estos resultados (31 % frente a 4 % de
   energía).
3. **Qué fuentes se revisaron y por qué se descartaron:** ASM Vol. 7 (sin datos térmicos de
   pastillas), fichas comerciales (sin números) y Goo y Na (valores de un modelo, no medidos).
4. **Cita para el informe:** Goo y Na (2014), con la reserva explícita de que sus valores son
   datos de entrada de un modelo FEM.
5. **Para el informe:** la tabla de §6 como comparación de pares disco–pastilla, con la
   advertencia de §7.
