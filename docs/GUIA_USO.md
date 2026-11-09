# Guía de uso de discfem

Esta guía explica, sin formalismos, qué hace el programa, cómo se usa y cómo se adapta a otros
discos, materiales y condiciones. Los detalles de cada clave del YAML están en el
[`README.md`](../README.md), y la física y las fórmulas exactas, en la sección 3 de
[`PLAN.md`](../PLAN.md).

---

## 1. Qué problema resuelve

Cuando un vehículo frena, la pastilla roza contra el disco y la energía cinética se convierte en
calor. `discfem` calcula **cómo se calienta el disco durante la frenada y cómo se enfría
después**: la temperatura en cada punto del disco y en cada instante.

Para no tener que resolver un problema 3D, aprovecha dos simetrías:

- **Simetría de revolución.** El disco es un sólido de revolución y, en promedio, se calienta
  igual en todos los ángulos (la pastilla pasa muchas veces por segundo). Basta con estudiar un
  corte radial: el plano $(r, z)$, donde $r$ es la distancia al eje y $z$ la dirección del espesor.
- **Simetría entre las dos caras.** Las dos caras del disco reciben el mismo calor, así que basta
  con la mitad del espesor. El plano medio ($z = 0$) se trata como una pared por la que no pasa
  calor.

El resultado es un rectángulo (o un polígono) en el plano $(r, z)$, que es lo que se malla.

```
 z ↑        top_inner          brake (pista)          top_outer
 e ┌────────────────────┬═══════════════════════┬───────────────┐
   │                    r_a                    r_b              │
   │ r_inner                                           r_outer  │
   │                                                            │
 0 └──────────────────────── symmetry ──────────────────────────┘
  r_i                                                          r_o  → r
```

La ecuación que se resuelve es la de calor transitoria en coordenadas cilíndricas:

$$
\rho c_p \frac{\partial T}{\partial t} - k\left[\frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial T}{\partial r}\right) + \frac{\partial^2 T}{\partial z^2}\right] = 0
$$

- **Entrada de calor:** en la pista (`brake`) entra un flujo de calor $q_{\text{fric}}(r, t)$ que
  decrece mientras el vehículo pierde velocidad y se anula cuando se detiene ($t = t_f$).
- **Salida de calor:** en las caras expuestas al aire, convección $h\,(T - T_\infty)$.
- **Simetría:** en $z = 0$ no pasa calor.

---

## 2. Cómo se usa (lo mínimo)

```bash
pip install -e ".[dev]"                                # una sola vez (ver README §1)
python -m discfem run   configs/ejemplo_auto.yaml      # una simulación completa
python -m discfem mesh  configs/ejemplo_auto.yaml      # solo la malla y su figura
python -m discfem study configs/estudio_auto.yaml      # convergencia + barridos de parámetros
```

Todo lo que define un caso vive en **un archivo YAML**: geometría, malla, materiales,
condiciones de frontera, frenada, paso de tiempo y salidas. Para crear un caso nuevo, lo normal es
copiar un YAML de `configs/` y cambiar valores. No hace falta tocar el código.

Para probar variantes sin editar el archivo, cualquier clave se cambia desde la línea de comandos:

```bash
python -m discfem run configs/ejemplo_auto.yaml --set load.vehicle.v0=33.33 --set convection.h=80
```

Si una clave está mal escrita o un valor no tiene sentido (un radio negativo, una pista fuera
del disco, un material sin fuente), el programa se detiene con un mensaje en español que dice qué
corregir.

---

## 3. Qué pasa por dentro cuando se ejecuta `run`

```
 config.yaml ──► validación ──► geometría ──► malla ──► materiales ──► ensamblaje ──► tiempo ──► salidas
                 (config.py)   (geometry.py)  (mesh/)   (materials.py) (fem/assembly) (timestepping) (postprocess)
```

1. **Lectura y validación** (`config.py`). Se carga el YAML, se aplican los `--set` y se revisa
   todo antes de calcular nada.
2. **Geometría** (`geometry.py`). Se construye el contorno de la sección (rectángulo o polígono)
   y se le pone nombre a cada tramo del borde (`brake`, `r_outer`, etc.).
3. **Malla** (`mesh/`). Se divide la sección en triángulos. La malla es fina justo debajo de la
   pista, donde la temperatura cambia rápido, y gruesa lejos de ella. Hay dos generadores: Gmsh
   (general) y uno estructurado propio (solo rectángulos, útil para pruebas). El parámetro
   `refine_factor` hace toda la malla más fina o más gruesa de una vez.
4. **Materiales** (`materials.py`). Se leen $k$, $\rho$ y $c_p$ de la librería
   `configs/materiales.yaml` y se asignan a cada región del disco. También se calcula qué
   fracción del calor de fricción se queda en el disco ($\gamma$) y qué fracción se va a la
   pastilla.
5. **Ensamblaje** (`fem/assembly.py`). Con elementos finitos lineales (P1) se construyen tres
   matrices y un vector:
   - $\mathbf{M}$ (capacidad térmica: cuánta energía hace falta para calentar cada zona),
   - $\mathbf{K}$ (conducción: cómo se mueve el calor dentro del disco),
   - $\mathbf{H}$ (convección: cuánto calor sale por las caras expuestas),
   - $\mathbf{F}$ (lo que entra: el flujo de fricción y el aire a $T_\infty$).

   Todas las integrales llevan el peso $r$, que es lo que convierte el problema plano en uno
   axisimétrico (un anillo lejos del eje tiene más volumen que uno cerca). El factor $2\pi$ de la
   vuelta completa se omite en el ensamblaje y solo se añade al reportar energías físicas.
6. **Avance en el tiempo** (`timestepping.py`). El problema queda como
   $\mathbf{M}\,\dot{\mathbf{T}} + (\mathbf{K}+\mathbf{H})\,\mathbf{T} = \mathbf{F}(t)$
   y se avanza paso a paso con el **esquema $\theta$**. Por defecto se usa Crank–Nicolson
   ($\theta = 1/2$) con dos pasos iniciales de Euler implícito (Rannacher), que evitan
   oscilaciones al arrancar la frenada de golpe. La matriz se factoriza una vez y se reutiliza en
   todos los pasos, por eso una corrida típica tarda menos de un segundo. El instante de parada
   $t_f$ siempre cae exactamente en un paso de tiempo.
7. **Salidas** (`postprocess.py`). Se guardan temperaturas, figuras, un resumen y un balance de
   energía (ver §5).

---

## 4. Cómo generaliza el programa

La idea central: **el código no sabe nada del disco concreto**. Todo lo particular (forma,
material, tipo de frenada, condiciones de borde) entra por el YAML. Esto es lo que se puede
cambiar.

### 4.1 Geometría

| Qué se quiere | Cómo |
|---|---|
| Otro disco sólido (otros radios, otro espesor, otra pista) | `geometry.type: rectangle` con `r_inner`, `r_outer`, `half_thickness` y `track: {r_a, r_b}` |
| Una sección no rectangular (cubo más grueso, escalones, chaflanes) | `geometry.type: profile` con una lista de vértices $(r, z)$ en sentido antihorario y un nombre de frontera para cada lado |
| Zonas de otro material dentro del disco (una capa bajo la pista, un inserto) | `geometry.regions` con un polígono por región; Gmsh genera una malla que respeta las interfaces |
| Un problema plano en vez de axisimétrico (solo para validar) | `axisymmetric: false` (cambia el peso $r$ por $1$) |

Ejemplos incluidos: `ejemplo_auto.yaml` (disco de auto de 20 mm), `ejemplo_moto.yaml` (disco
flotante de moto de 4,5 mm) y `ejemplo_perfil.yaml` (sección con un cubo interior más grueso).

**Lo que no cubre:** discos ventilados, perforados o ranurados. No tienen simetría de revolución y
requieren un modelo 3D.

### 4.2 Materiales

Los materiales viven en una librería aparte, `configs/materiales.yaml`. Cada material tiene $k$,
$\rho$, $c_p$ y, opcionalmente, **la fuente** de donde salieron (`source`) y si se contrastaron
con ella (`verified`). Si se omiten, el material queda "sin fuente declarada" y no verificado (útil
para pruebas). Si un material no está verificado, cada corrida lo advierte. `verified: true` exige
una `source`.

Hay varias formas de elegir materiales en un caso, de la más cómoda a la más fina:

1. **Preset por vehículo:** `materials.preset: car` o `motorcycle` define disco y pastilla de una
   vez.
2. **Por región:** `materials.regions: {disc: stainless_steel_410, capa: ...}`.
3. **Material propio dentro del caso:** `materials.custom`, con las mismas reglas que la librería.
4. **Retoque de una propiedad:** `materials.overrides: {disc: {k: 48.0}}`, útil para estudios de
   sensibilidad (queda marcado como no verificado).

La **pastilla** no se malla: solo se usa para calcular $\gamma$, la fracción del calor que entra
al disco, a partir de las efusividades $\sqrt{k\rho c}$ del disco y de la pastilla. Cambiar de
pastilla cambia $\gamma$ y, como el problema es lineal, el aumento de temperatura cambia en la
misma proporción.

Materiales disponibles hoy: hierro fundido gris y pastilla orgánica de Limpert (verificados),
acero AISI 410 y dos pastillas sinterizadas (con fuente, sin verificar). La tabla completa está
en el README §4.

### 4.3 Condiciones de frontera

El borde de la sección tiene seis grupos con nombre fijo: `brake`, `top_inner`, `top_outer`,
`r_inner`, `r_outer` y `symmetry`. A cada uno se le asigna en el YAML una de tres condiciones:

| Condición | Significado | Uso típico |
|---|---|---|
| `flux` | entra un flujo de calor $q(r, t)$ | la pista de fricción |
| `convection` | sale calor al aire: $h\,(T - T_\infty)$ | caras y bordes expuestos |
| `insulated` | no pasa calor | plano de simetría, o una cara en contacto con algo aislante |

El coeficiente $h$ y la temperatura del aire $T_\infty$ pueden ser los mismos para todo
(`convection: {h, T_inf}`) o distintos por grupo (`convection.per_group`), por ejemplo un $h$
mayor en el borde exterior, que gira más rápido.

### 4.4 La carga (cómo entra el calor)

Tres opciones en `load.type`:

- **`braking`** (la principal): el flujo se calcula a partir de la frenada (masa, velocidad
  inicial, desaceleración, fracción de frenado del eje, discos por eje). Se elige además cómo se
  reparte el calor en la pista (`uniform_wear` o `uniform_pressure`) y $\gamma$ (`auto` desde los
  materiales, o un número fijo).
- **`constant_flux`**: un flujo constante $q$ durante un tiempo dado.
- **`table`**: un flujo que sigue una tabla $(t, q)$ con interpolación lineal; sirve para
  cualquier historia de carga medida o inventada.

### 4.5 Tiempo y malla

- `time.theta` elige el esquema ($\theta = 1$: Euler implícito; $\theta = 1/2$: Crank–Nicolson),
  `rannacher_steps` los pasos iniciales amortiguados, `dt` el paso (o `auto`, que da
  $\Delta t = t_f/100$) y `t_end` hasta cuándo se simula el enfriamiento.
- `mesh.refine_factor` controla la finura de toda la malla; `size_track`, `size_far` y
  `grading_distance` controlan cómo pasa de fina (bajo la pista) a gruesa.

### 4.6 Estudios automáticos

`python -m discfem study` toma un caso base y repite la simulación cambiando algo:

- **Estudio de malla:** refina la malla ($\text{rf} = 1, 2, 4, 8$), mide cuánto cambia la
  temperatura y elige el primer nivel con cambio $< 1\,\%$.
- **Estudio de $\Delta t$:** divide el paso de tiempo y estima el orden de convergencia.
- **Barridos:** cualquier clave del YAML con una lista de valores (masa, velocidad, $h$,
  pastilla, material...). Produce tablas y figuras PDF listas para un informe.

Ver `configs/estudio_auto.yaml` como plantilla.

---

## 5. Qué produce y cómo leerlo

Cada corrida escribe en `output.dir` (por ejemplo `results/ejemplo_auto/`):

- `summary.json`: temperatura máxima (valor, instante y posición), temperaturas finales, energía
  que entró, salió y se almacenó, materiales usados con su fuente.
- `tmax.csv` y `probes.csv`: la evolución en el tiempo del máximo y de puntos elegidos (sondas).
- `field.png`, `tmax.png`, `depth_profiles.png`, `mesh.png`: figuras del campo, de la evolución y
  de la malla.
- `field.pvd` y `field_*.vtu`: para abrir en ParaView.
- `config.yaml` y `materiales.yaml`: copias exactas de lo que se usó. Con esa carpeta sola se
  puede repetir la corrida.

**Un control de calidad que conviene mirar siempre:** `energy.balance_relative_residual` en
`summary.json`. Compara la energía que entró con la que salió más la almacenada. Debe ser del
orden de $10^{-14}$. Si no lo es, algo está mal configurado. Ojo: que cierre no significa que el
resultado sea preciso; para eso están los estudios de malla y de $\Delta t$.

---

## 6. Límites del modelo y errores comunes

**Lo que el modelo supone** (y conviene decir al presentar resultados):

- El calor de la pastilla se promedia en el ángulo (no hay "mancha caliente" que gira).
- Propiedades constantes con la temperatura; sin radiación, sin deformación térmica, sin
  desgaste.
- Contacto perfecto disco–pastilla para repartir el calor ($\gamma$ por efusividades).
- Una sola frenada seguida de enfriamiento (no frenadas repetidas).
- El calor que se va a la pastilla no se sigue.

**Errores comunes:**

| Síntoma                                                       | Causa probable                                                               |
| ------------------------------------------------------------- | ---------------------------------------------------------------------------- |
| "el polígono debe recorrerse en sentido ANTIHORARIO"          | los vértices de `points` (o de una región) están al revés; invertir el orden |
| "el polígono se autointersecta"                               | dos lados del perfil se cruzan                                               |
| Error al usar `profile` o `regions` con la malla estructurada | requieren `mesh.backend: gmsh`                                               |
| Advertencia `UnverifiedMaterialWarning` en cada corrida       | algún material tiene `verified: false`; es informativa, no un error          |
| El material no se encuentra                                   | el nombre no está en la librería (o le falta `k`, `rho` o `cp`)              |
| `test_rectangle_unchanged` falla en otra máquina              | otra versión de Gmsh o numpy; instalar con `requirements-lock.txt`           |
| Corridas grandes ocupan mucho disco                           | subir `output.save_every`                                                    |
