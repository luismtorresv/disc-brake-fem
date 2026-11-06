# F7 · T2 — Catálogo de materiales con fuente (Limpert) y ejemplo de automóvil verificado

Todas las pruebas pasaron (231 passed) sin relajar tolerancias. **Pendiente:** Limpert no trae
propiedades de acero inoxidable (discos de moto) ni de pastillas sinterizadas, así que el preset
`motorcycle` sigue sin crearse (§7).

## 1. Qué se pidió

Continuación de la Fase 7 (§5 de `PLAN.md`) con tu indicación de usar materiales reales de
discos (hierro fundido para autos, acero inoxidable para motos), elegibles por tipo de vehículo.
Me diste como fuente *Brake Design and Safety* de R. Limpert y me pediste buscar ahí los datos
(dejaste el PDF en `legacy/`).

## 2. Qué se hizo

- `configs/materiales.yaml`:
  - nuevas entradas `gray_cast_iron_limpert`, `pad_limpert` y `lining_limpert`, con
    `verified: true` y `source` con tabla y página;
  - el preset `car` ahora usa `gray_cast_iron_limpert` + `pad_limpert`;
  - los materiales del prototipo se conservan sin cambios.
- `configs/ejemplo_auto.yaml`: nuevo. El caso del automóvil a 100 km/h con `preset: car`, la
  malla del nivel seleccionado (rf = 2) y $\theta = 1/2$ + 2 Rannacher.
- `tests/test_profile_regions.py`:
  - 2 pruebas nuevas: valores de Limpert, y ejemplo sin materiales no verificados;
  - las expectativas del preset `car` se actualizaron al nuevo material.
- `.gitignore`: `legacy/*.pdf`. El libro tiene derechos de autor y no debe versionarse.
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

Se extrajo el texto del PDF con `pypdf`, instalado solo en una carpeta temporal y no en el venv
del proyecto. Luego se buscó "thermal conductivity", "specific heat" y "Table". La única tabla
de propiedades térmicas del libro es la **Tabla 3-1A, "Brake Design Values, SI Units"** (p. 123
impresa, p. 133 del PDF). Su versión imperial es la 3-1B. El texto que la introduce dice:
"Typical material properties for drum or disc, and asbestos-based lining and pad material are
listed in Table 3-1 below. Semi-metallic pad materials will have increased values for ρ and k."

Valores de la tabla SI y su conversión. El libro da $k$ en $\mathrm{N\,m/(h\,K\,m)}$, es decir,
J por hora; para pasar a $\mathrm{W/(m\,K)}$ se divide por 3600:

$$
k\ [\mathrm{W/(m\,K)}] = \frac{k\ [\mathrm{N\,m/(h\,K\,m)}]}{3600}
$$

| Material | $\rho$ [kg/m³] | $c$ [J/(kg·K)] | $k$ libro [N·m/(h·K·m)] | $k$ [W/(m·K)] | $\alpha$ libro [m²/h] |
|---|---|---|---|---|---|
| Disco (tambor o disco) | 7228 | 419 | 174 465 | 48,46 | 0,0576 |
| Pastilla | 2595 | 1465 | 4362 | 1,212 | 0,0011 |
| Forro (tambor) | 2034 | 1256 | 4174 | 1,159 | 0,00163 |

**Verificación** (por eso `verified: true`):

1. **Contra la tabla imperial 3-1B:**
   - disco: 28 BTU/(h·°F·ft) = 48,46 W/(m·K) y 0,10 BTU/(lb·°F) = 418,7 J/(kg·K);
   - pastilla: 0,7 BTU/(h·°F·ft) = 1,21 W/(m·K), 162 lb/ft³ = 2595 kg/m³ y 0,35 BTU/(lb·°F) =
     1465 J/(kg·K).
2. **Contra la difusividad tabulada $\alpha = k/(\rho c)$:**
   - disco: $1{,}6002\times10^{-5}\ \mathrm{m^2/s}$ calculado frente a
     $0{,}0576\ \mathrm{m^2/h} = 1{,}6000\times10^{-5}\ \mathrm{m^2/s}$;
   - pastilla: $0{,}00115\ \mathrm{m^2/h}$, que es el 0,0011 de la tabla redondeado a dos cifras.

**Discrepancia dentro del libro:** para el disco, la tabla imperial da 455 lb/ft³ ≈ 7288 kg/m³,
y la SI da 7228 (0,8 % de diferencia). Se usa el valor SI y queda anotado en la cabecera de
`materiales.yaml`.

| Elemento | Dónde |
|---|---|
| Entradas con fuente | `configs/materiales.yaml`, sección `materials` (`*_limpert`) |
| Preset por vehículo | `configs/materiales.yaml`, `presets.car`; resolución en `solver.py:94` |
| Prueba de los valores | `test_limpert_materials_match_table_3_1A`, `tests/test_profile_regions.py` |

## 4. Decisiones

- **Entradas nuevas en vez de modificar las del prototipo.** El caso de referencia y T8/T9 están
  definidos con los supuestos del prototipo, y el PLAN dice que, si cambian las propiedades, se
  recalculan, no se ajustan. Por eso esas entradas siguen igual, marcadas como "no usar para el
  informe", y los resultados nuevos van en `ejemplo_auto.yaml`.
- **`verified: true` para las entradas de Limpert.** Los valores se contrastaron con la tabla
  SI, la imperial y la difusividad. Si prefieres revisarlos tú en la p. 123 antes de usarlos en
  el informe, basta con volver a ponerlos en `false`.
- **El nombre `pad_limpert`, sin "orgánica" ni "sinterizada".** El libro dice que son pastillas
  a base de asbesto, que ya no se usan, aunque sus propiedades sirven como orden de magnitud de
  una pastilla orgánica. El `description` y el `source` lo dicen.
- **Sin preset `motorcycle`.** No hay datos del acero inoxidable ni de la pastilla sinterizada en
  la fuente que me diste, y no voy a inventarlos (regla 4).

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
.venv\Scripts\python -m discfem run configs/ejemplo_auto.yaml
```

Salida real:

```
231 passed in 19.72s
All checks passed!
```

Conversión y comprobaciones cruzadas (salida real):

```
k disco 48.4625 k pastilla 1.2116666666666667 k forro 1.1594444444444445
alpha disco 1.600197719555217e-05 tabla 1.6e-05
alpha pastilla 3.1871916107154523e-07 tabla 3.055555555555556e-07
imperial: 28 BTU/hFft = 48.46058  0.7 -> 1.2115144999999998  0.10 BTU/lbF = 418.68000000000006  455 lb/ft3 = 7288.4175  162 -> 2594.997
gamma Limpert 0.8495038368213138
```

`python -m discfem run configs/ejemplo_auto.yaml` (salida real; ya no aparece la advertencia de
materiales no verificados):

```
Caso 'ejemplo_auto_100kmh': 3028 nodos, 5580 triángulos, 1080 pasos (θ = 0.5, Rannacher = 2, Δt = 0.05556 s)
  γ = 0.8495   t_f = 5.556 s
  T_max = 106.695 en t = 3.500 s (r = 0.1095, z = 0.0100)
  T_max(t_end = 60 s) = 75.879
  energía [J, con 2π]: entra 86028.5, sale por convección 2825.6, almacenada 83202.9; cierre del balance 2.5e-14
  sonda (r=0.11, z=0.01): T_max = 106.695 en t = 3.500 s; T_final = 75.431
  sonda (r=0.11, z=0.005): T_max = 92.726 en t = 5.723 s; T_final = 75.431
  sonda (r=0.11, z=0): T_max = 91.694 en t = 6.889 s; T_final = 75.432
  tiempo de cómputo: 0.50 s   salidas: results\ejemplo_auto
```

## 6. Resultados

| Magnitud (rf = 2, $\theta = 1/2$ + 2 Rannacher) | Supuestos del prototipo | Limpert (verificado) |
|---|---|---|
| $k$ / $\rho$ / $c$ del disco | 50 / 7200 / 460 | 48,46 / 7228 / 419 |
| $\rho c$ del disco [J/(m³·K)] | $3{,}31\times10^{6}$ | $3{,}03\times10^{6}$ |
| $\gamma$ | 0,8692 | 0,8495 |
| $T_{\max}$ | 102,91 °C (Fase 5, exploración) | 106,70 °C en $t = 3{,}50$ s |
| $T_{\max}(60\ \mathrm{s})$ | 72,80 °C (Fase 5, rf = 1) | 75,88 °C |

Con los valores de Limpert, el disco se calienta más (+3,8 °C de pico y +3 °C al final) aunque
$\gamma$ sea un poco menor. La razón es que su capacidad térmica volumétrica es un 8,5 % menor
($c = 419$ frente a 460). Al final, con el disco casi uniforme,

$$
\Delta T \approx \frac{E}{\rho c\, V},
$$

y ese 8,5 % explica casi todo el cambio:

$$
\frac{75{,}88 - 25}{72{,}80 - 25} = 1{,}064,
$$

combinado con un 2,3 % menos de energía por $\gamma$.

## 7. Limitaciones y riesgos

- **Motocicleta pendiente.** Limpert no tiene acero inoxidable ni pastilla sinterizada. Fuera
  de la Tabla 3-1, el libro solo menciona motocicletas a propósito de la norma FMVSS 122. Para
  crear el preset `motorcycle` hace falta otra fuente, por ejemplo la ficha del fabricante del
  acero AISI 410/420 templado o un manual (ASM Handbook). Lo que encontré en la web era del 410S
  ferrítico o de fuentes secundarias (reporte F7_T1, §7).
- **La pastilla de Limpert es de base asbesto.** Sirve como valor de diseño clásico. Las
  pastillas semimetálicas tienen $\rho$ y $k$ mayores, y eso subiría la efusividad de la pastilla
  y bajaría $\gamma$.
- **Propiedades constantes.** Limpert da valores de diseño a una sola temperatura, y la
  dependencia con $T$ está fuera de alcance (Fase 7: "nada de $k(T)$").
- **El caso de referencia (T8/T9) sigue con los supuestos del prototipo** a propósito. El
  informe debe usar `ejemplo_auto.yaml` (o similares) para los resultados físicos.

## 8. Para revisar y para el informe

1. **Unidades del libro:** $k$ en $\mathrm{N\,m/(h\,K\,m)}$ → ÷3600 → $\mathrm{W/(m\,K)}$. Es un
   error fácil (factor 3600), y la prueba lo verifica contra la tabla imperial.
2. **Por qué cambian los resultados con los materiales de Limpert:** $\rho c$ menor ⇒ más
   $\Delta T$ para la misma energía; $\gamma$ menor ⇒ un poco menos de energía. El primer efecto
   domina.
3. **La cita para el informe:** Limpert, R. *Brake Design and Safety*, 2.ª ed., SAE, Tabla 3-1A,
   p. 123. Ya está en el `.bib` del anteproyecto como `limpert`.
4. **Qué falta para la comparación auto frente a moto:** una fuente para el acero inoxidable
   martensítico y para la pastilla sinterizada.
5. **Para el informe:** `results/ejemplo_auto/` (`field.png`, `tmax.png`,
   `depth_profiles.png`) como resultado principal con materiales verificados, y la tabla de §6
   (prototipo frente a Limpert) como análisis de sensibilidad a los materiales.
