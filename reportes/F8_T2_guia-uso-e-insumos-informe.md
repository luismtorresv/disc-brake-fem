# F8_T2 — Guía de uso y documento de insumos para el informe

Todas las pruebas pasan (243 passed); no se modificó código. Se regeneró `results/ejemplo_perfil/`,
que estaba desactualizado (ver §4).

## 1. Qué se pidió

Dos documentos en Markdown (PLAN §5, Fase 8 — Documentación):

1. Una explicación descriptiva, no formal, para quien quiera usar el programa: cómo funciona y cómo
   generaliza a otros materiales, condiciones de frontera y geometrías, siguiendo
   `legacy/criterios_proyecto.pdf` (criterio de Generalidad: "diferentes geometrías, materiales y
   condiciones de frontera" a partir de archivos de entrada).
2. Un documento con la información más significativa para redactar el informe final, siguiendo
   `legacy/criterios_proyecto.pdf` (estructura de artículo, 3500–5000 palabras) y
   `legacy/rubrica_informe.pdf`.

## 2. Qué se hizo

- `docs/GUIA_USO.md`: guía de uso descriptiva (problema, comandos, flujo interno, generalización
  por eje, salidas, límites y errores comunes).
- `docs/INSUMOS_INFORME.md`: insumos del informe por sección del artículo, con números y rutas de
  figuras, tabla criterio ↔ material, advertencias de datos, presupuesto de palabras y preguntas
  de sustentación.
- `results/ejemplo_perfil/` (ignorado por git): regenerado con la configuración actual.
- `reportes/INDEX.md`: una línea nueva.

## 3. Cómo funciona

No hay algoritmo nuevo. Los documentos se construyeron leyendo `PLAN.md` (§3 y ESTADO), `README.md`,
las secciones 6–8 de los reportes F3_T1, F6_T1, F6_T2, F7_T1, F7_T6 y F7_T7, las tablas de
`results/estudio_*/` y los `summary.json` de `results/ejemplo_*/`. Cada número del documento de
insumos lleva la ruta de la que sale.

| Concepto en los documentos | Dónde está en el código |
|---|---|
| Esquema $\theta$ con Rannacher y $t_f$ exacto | `src/discfem/timestepping.py` |
| $\mathbf{M}$, $\mathbf{K}$, $\mathbf{H}$, $\mathbf{F}$ con peso $r$ (sin $2\pi$) | `src/discfem/fem/assembly.py` |
| $q_{\text{fric}}$ y $\gamma$ por efusividades | `src/discfem/loads.py`, `src/discfem/materials.py` |
| `rectangle`, `profile`, `regions` y validación del polígono | `src/discfem/geometry.py` (`_validate_polygon`, l. 103) |
| Restricción de la malla estructurada a rectángulos | `src/discfem/mesh/structured.py` l. 64 |
| Claves del YAML y `--set` | `src/discfem/config.py` |

## 4. Decisiones

- **Markdown en `docs/`**, a pedido del usuario (se había empezado un documento de Claude Docs y se
  abandonó; quedó un borrador vacío en claude.ai).
- **La guía no repite el README**: explica el "por qué" y remite al README para cada clave.
- **El documento de insumos se ordena por secciones del artículo**, no por fases, porque así se
  escribe el informe; y abre con la tabla de criterios para ver dónde se gana cada punto.
- **Se regeneró `results/ejemplo_perfil/`**: su `summary.json` tenía $\gamma = 0{,}869$ (materiales
  supuestos del prototipo), porque se generó antes de que el preset `car` pasara a los materiales
  de Limpert. Citar el valor viejo (102,93 °C) habría sido incorrecto.
- Los mensajes de error citados en la guía se copiaron del código, no de memoria.

## 5. Cómo verificarlo

```
.venv/Scripts/python -m pytest -q
```
```
........................................................................ [ 88%]
...........................                                              [100%]
243 passed in 25.31s
```

```
.venv/Scripts/python -m discfem run configs/ejemplo_perfil.yaml
```
```
  γ = 0.8495   t_f = 5.556 s
  T_max = 106.707 en t = 3.500 s (r = 0.1094, z = 0.0100)
  T_max(t_end = 60 s) = 73.898
  energía [J, con 2π]: entra 86028.5, sale por convección 3670.1, almacenada 82358.4; cierre del balance 6.1e-15
  sonda (r=0.11, z=0.01): T_max = 106.706 en t = 3.500 s; T_final = 73.864
  sonda (r=0.11, z=0): T_max = 91.686 en t = 6.834 s; T_final = 73.865
  sonda (r=0.0725, z=0.016): T_max = 71.584 en t = 49.282 s; T_final = 71.410
  tiempo de cómputo: 0.51 s   salidas: results\ejemplo_perfil
```

Los comandos de ejemplo de la guía (`run`, `mesh`, `study`, `--set`) no se ejecutaron en esta tarea
salvo el anterior; están probados en F5_T1, F6_T1 y F8_T1.

## 6. Resultados

- Perfil con cubo, materiales de Limpert: $T_{\max} = 106{,}71$ °C (igual al rectángulo,
  106,70 °C), $T(60\ \mathrm{s}) = 73{,}9$ °C frente a 75,9 °C del rectángulo. El reporte F7_T1 daba 102,93 / 70,96 °C con los
  materiales supuestos; la conclusión cualitativa (pico igual, final menor) se mantiene.
- Los demás números del documento de insumos son los ya reportados (F3_T1, F6_T1, F6_T2, F7_T6,
  F7_T7) y coinciden con los archivos de `results/` revisados.

## 7. Limitaciones y riesgos

- Las referencias bibliográficas sugeridas nuevas (Rannacher 1984, SciPy, NumPy, Ashby) no se
  verificaron contra la fuente; el documento lo advierte.
- El reporte F7_T1 conserva los valores del perfil anteriores al cambio de materiales; en F8_T3 se
  le añadió una nota al inicio que remite a los valores nuevos.
- Las figuras sugeridas en F3_T1 (T7 por esquema, error de T6 en log-log) siguen sin generarse.
- Queda un borrador vacío de Claude Docs ("discfem — Guía de uso del programa") en claude.ai.

## 8. Para revisar y para el informe

1. La tabla de §0 de `docs/INSUMOS_INFORME.md`: qué criterio de la rúbrica cubre cada resultado.
2. La tabla de generalidad (§3.4 de los insumos): es el argumento del criterio de Generalidad.
3. La idea "el error numérico ($< 0{,}1$ °C) es mucho menor que la incertidumbre de los datos
   (pastilla: hasta 38 °C)", que debería ser la conclusión central.
4. Las advertencias de §7 de los insumos deben aparecer en el informe.
5. Figura para el informe: `results/ejemplo_perfil/field.png` regenerada, como ejemplo de
   geometría no rectangular.
