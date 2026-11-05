# disc-brake-fem — contexto para Claude Code

Proyecto del curso Métodos Numéricos para EDP (EAFIT, 2026). Solver FEM en Python para la
ecuación de calor transitoria en la **sección axisimétrica (r, z)** de un disco de freno.
La fuente de verdad del plan de trabajo es `PLAN.md`; la de la física es la sección 3 de ese
archivo. Si algo del código contradice `PLAN.md`, avisa antes de cambiar nada.

## Alcance

Dentro:
- Dominio 2D Ω = (r_i, r_o) × (0, e) (medio espesor; z = 0 es plano de simetría), y más adelante
  perfiles poligonales (r, z) y regiones con materiales distintos.
- Malla triangular P1 (Gmsh + generador estructurado de respaldo), refinable por parámetro.
- Materiales y propiedades (k, ρ, c_p) y parámetros de frenado configurables por archivo YAML.
- Esquema θ en el tiempo, estudios de convergencia, validación analítica.

Fuera (no implementar sin que se pida):
- Discos ventilados, perforados o ranurados (requieren 3D). Se documentan como trabajo futuro.
- Fuente móvil de la pastilla (el modelo promedia en θ), radiación, acoplamiento termomecánico,
  desgaste, resistencia térmica de contacto.

## Convenciones (no negociables)

- Unidades SI en todo el código: m, s, kg, W, K. Temperaturas en °C o K, pero una sola por
  caso; la ecuación solo usa diferencias.
- Nombres de frontera: `brake` (pista de fricción, z = e, r_a < r < r_b), `top_inner`
  (r_i..r_a en z = e), `top_outer` (r_b..r_o en z = e), `r_inner`, `r_outer`, `symmetry` (z = 0).
  Cada arista de borde pertenece a exactamente un grupo. La condición de cada grupo se asigna
  en el YAML (`convection` | `insulated` | `flux`).
- Normal **exterior**. Condición en `brake`: `-k ∇T·n = -q_fric`, es decir calor que ENTRA.
  En el vector de carga aporta `+∫ q_fric φ_i r ds`.
- El factor 2π de la revolución se **omite** en todo el ensamblaje (M, K, H, F). Cualquier
  balance de energía que compare con potencias físicas debe incluir el 2π explícitamente.
- Peso de las integrales: `r` (axisimétrico) o `1` (plano, solo para validación). Es una
  opción del solver, nunca un valor fijo en el código.
- Esquema θ: `[M + θΔt(K+H)] Tⁿ⁺¹ = [M − (1−θ)Δt(K+H)] Tⁿ + Δt[θFⁿ⁺¹ + (1−θ)Fⁿ]`.
  `q_fric` se anula para t > t_f: la malla temporal debe contener t_f exactamente.

## Reglas de trabajo

1. Una fase de `PLAN.md` por sesión. No avances a la siguiente hasta que `pytest -q` pase.
2. Escribe las pruebas de la fase **antes** o junto con el código, y pega la salida de pytest.
3. Ensamblaje vectorizado con numpy/scipy; sin bucles de Python sobre elementos.
4. Nunca inventes propiedades de materiales. Toda entrada de `materiales.yaml` lleva `source`
   y `verified`. Mientras `verified: false`, cada corrida imprime una advertencia con la lista.
5. Los parámetros viven en el YAML, no en el código. Una corrida debe poder reproducirse con
   solo su `config.yaml` (el programa copia la config a la carpeta de resultados).
6. Identificadores en inglés; docstrings, comentarios y mensajes de error para el usuario en
   español.
7. Un commit por fase, con mensaje descriptivo. No reescribas historia ni hagas push sin que
   se pida.
8. Si una prueba falla, no relajes la tolerancia para hacerla pasar: explica la causa.

## Comandos

```
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
python -m pytest -q
python -m discfem run configs/ejemplo_frenado.yaml
python -m discfem mesh configs/ejemplo_frenado.yaml     # solo malla + figura
python -m discfem study configs/estudio_malla.yaml      # convergencia
```

## Dependencias

Python ≥ 3.10; numpy, scipy, pyyaml, gmsh (rueda de pip), meshio (exportar VTK), matplotlib;
pytest y ruff en desarrollo.

## Definición de terminado (por fase)

Código + pruebas pasando + docstrings + sección "ESTADO" de `PLAN.md` actualizada + commit.


## Reportes (obligatorio)

Al terminar CADA tarea (no solo cada fase), escribe un reporte en `reportes/` antes de dar la
tarea por cerrada. Nombre: `reportes/F<fase>_T<n>_<nombre-corto>.md` (ej. `F2_T1_matrices-elementales.md`).
Añade además una línea al índice `reportes/INDEX.md` (fecha, fase, tarea, enlace, estado).

Idioma: español. Público: un estudiante que domina FEM estacionario visto en clase pero no el caso
transitorio ni la geometría axisimétrica; explica esos puntos sin dar nada por sabido.

Cada reporte tiene exactamente estas secciones:
1. **Qué se pidió**: la tarea, citando la sección de PLAN.md que aplica.
2. **Qué se hizo**: lista de archivos creados o modificados, con ruta y una línea por archivo.
3. **Cómo funciona**: explicación clara, en prosa, del algoritmo. Incluye una tabla
   "ecuación ↔ código" (símbolo o fórmula de PLAN.md ↔ función y línea).
4. **Decisiones**: qué elegiste, qué alternativas descartaste y por qué.
5. **Cómo verificarlo**: comandos exactos y la salida REAL pegada (pytest, números). Si no
   ejecutaste algo, escribe "no ejecutado"; nunca resumas una salida que no viste.
6. **Resultados**: valores numéricos obtenidos y comparación con la referencia de la sección 6.
7. **Limitaciones y riesgos**: qué no cubre, qué podría fallar, qué supuestos quedan sin verificar.
8. **Para revisar y para el informe**: 3 a 5 puntos que el usuario debe entender para poder
   defender este código, y qué figura o tabla de esta tarea sirve para el informe final.

Reglas del reporte:
- No afirmes que algo funciona sin una ejecución que lo muestre en la sección 5.
- Si una prueba falló o relajaste algo, dilo en la primera línea del reporte.
- Si hay una fórmula, escríbela una vez en texto plano y explica cada símbolo.
- Al cerrar la tarea, en el chat entrega un resumen de máximo 6 líneas y la ruta del reporte.