# F8_T3 — Notación matemática en LaTeX en reportes y documentos

No se modificó código ni se corrieron pruebas nuevas: es un cambio de formato en los `.md`. La
compilación de las fórmulas en Obsidian **no se verificó** visualmente desde aquí (§5).

## 1. Qué se pidió

Reemplazar el texto matemático de los reportes y documentos Markdown hechos hasta ahora por LaTeX,
para que Obsidian lo renderice (complemento de la Fase 8 de `PLAN.md`, documentación).

## 2. Qué se hizo

- `docs/GUIA_USO.md`, `docs/INSUMOS_INFORME.md`: toda la notación en LaTeX, conservando los
  cambios de formato de tablas que hiciste en Obsidian.
- `reportes/F1_T1` a `reportes/F7_T7`: notación en LaTeX en prosa, tablas y fórmulas.
- `reportes/F8_T1`, `reportes/F8_T2`: los pocos símbolos que tenían.
- `reportes/F7_T1`: además, una nota al inicio con los valores actuales del perfil (materiales de
  Limpert), que F8_T2 había dejado como pendiente.
- `reportes/F0_T1`: sin cambios (solo tiene valores de configuración y una cita literal del `.tex`).
- `reportes/INDEX.md`: una línea nueva.

## 3. Cómo funciona

Convenciones aplicadas (las que entiende Obsidian, que usa MathJax):

| Caso | Forma |
|---|---|
| Símbolo en la prosa | `$\theta$`, `$\Delta t$`, `$T_{\max}$`, `$q_{\text{fric}}$` |
| Fórmula aparte (antes en bloque de código) | bloque `$$ … $$`; varias líneas con `\begin{aligned}` |
| Matrices y vectores | negrita: `$\mathbf{M}$`, `$\mathbf{T}^{n+1}$`, `$\mathbf{1}^\top$` |
| Valor absoluto dentro de tablas | `\lvert … \rvert` (una barra vertical rompería la tabla) |
| Coma decimal dentro de math | `0{,}5` (sin llaves, LaTeX añade un espacio tras la coma) |
| Potencias de diez | `$10^{-14}$`, `$2{,}5\times10^{-14}$` en vez de `1e-14` |

Qué **no** se convirtió, a propósito:

- las salidas reales de programas (bloques de código en §5 de cada reporte), porque son copia
  literal de lo que imprimió el programa;
- nombres de archivos, claves del YAML e identificadores de código (van en `` `…` ``);
- unidades sueltas en la prosa (m/s², kg/m³) y citas textuales de fuentes.

Además se añadieron, donde el texto las mencionaba sin escribirlas, las fórmulas explícitas: la
solución de Duhamel de T7 (F3), la matriz $\mathbf{H}^e$ (F2), los casos cerrados de $q$ (F4), la
cota de T9 (F5), la extrapolación de Richardson (F6) y la interpolación usada al digitalizar la
Figura 6 de Zhou (F7_T5). Todas salen de `PLAN.md` §3 y §6 o del propio reporte.

## 4. Decisiones

- **`$…$` y `$$…$$`, no `\(…\)`**: es la sintaxis nativa de Obsidian y de GitHub.
- **No tocar las salidas de los programas**: la regla de los reportes exige pegar la salida real;
  pasarla a LaTeX la falsearía.
- **No reescribir F7_T1 con los valores nuevos**: el reporte documenta lo que se midió en su
  momento; se añadió una nota en vez de cambiar los números.
- **`README.md` y `PLAN.md` no se convirtieron**: no los pediste y el `PLAN.md` es tu fuente de
  verdad. Si los quieres en LaTeX, es la misma operación.

## 5. Cómo verificarlo

Se revisaron todos los archivos con un script (en la carpeta temporal del trabajo) que quita los
bloques de código y el LaTeX, y lista lo que queda. Comprueba tres cosas: que no queden `$`
impares en ninguna línea, que ningún `$` tenga espacio justo por dentro (Obsidian no lo
renderiza) y que no haya `|` dentro de math en filas de tabla.

Salida real de la primera pasada (resumida a los problemas encontrados):

```
reportes\F7_T7_caso-moto.md issues: 1
    ('espacio', 'T_{\\max} = ')
```

Se corrigió ese caso. En la segunda pasada, los únicos restos de símbolos Unicode fuera de LaTeX
son unidades, la cita del libro de Limpert ("increased values for ρ and k") y los valores de
configuración de F0_T1, que se dejaron así a propósito.

**No ejecutado:** abrir los archivos en Obsidian y mirar el render. Si alguna fórmula no compila,
lo más probable es `\begin{smallmatrix}` en la tabla de F2_T1 §3 (MathJax lo admite, pero es la
construcción menos común que se usó).

## 6. Resultados

19 archivos revisados, 18 modificados. No cambió ningún número.

## 7. Limitaciones y riesgos

- El render no se verificó en Obsidian (§5).
- En GitHub, las fórmulas dentro de tablas a veces no se renderizan; en Obsidian sí.
- El `INDEX.md` y las salidas de programa siguen en texto plano.

## 8. Para revisar y para el informe

1. Las fórmulas de los reportes se pueden copiar casi tal cual al `.tex` del informe (cambiando
   `{,}` según el paquete de idioma que uses).
2. Las fórmulas nuevas (Duhamel, Richardson, cota de T9, $\mathbf{H}^e$) son las que la rúbrica
   pide ver explícitas en la metodología.
3. Revisa en Obsidian F2_T1 §3 (matrices elementales) y F3_T1 §3 (esquema $\theta$ y balance): son
   las fórmulas más largas.
