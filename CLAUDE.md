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
4. Nunca inventes propiedades de materiales. Toda entrada de `materiales.yaml` que se use en el
   informe lleva `source` y `verified`. Para uso general ambos son opcionales en el código (por
   omisión: `source: "sin fuente declarada"`, `verified: false`); `verified: true` exige `source`.
   Mientras `verified: false`, cada corrida imprime una advertencia con la lista.
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


## Reportes

Ya no se escriben reportes por tarea. `reportes/` conserva solo los que explican el código y su
validación (índice en `reportes/INDEX.md`); no añadas nuevos salvo que se pida.
