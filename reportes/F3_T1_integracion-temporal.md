# F3 · T1 — Integración temporal (esquema θ, Rannacher, t_f exacto), balance de energía y validación analítica

Todas las pruebas pasaron (101 passed) sin relajar tolerancias. Durante el desarrollo fallaron
dos cosas: una era un error mío en el caché de factorizaciones, la otra una expectativa
equivocada en una prueba; ambas se explican en §4. Hay una discrepancia con el prototipo: aquí
Crank–Nicolson da errores menores en T7 (§6).

## 1. Qué se pidió

Fase 3 de `PLAN.md` (§5, "Fase 3 — Integración temporal y validación analítica"):

- `timestepping.py` según la §3.5: esquema θ con la factorización reutilizada, n pasos
  iniciales de Euler implícito (Rannacher), Δt automático = t_f/100 ajustado para que t_f sea
  múltiplo de Δt, y balance de energía según la §3.6.
- Pruebas T4 (balance), T5 (estado uniforme) y T7 (semi-infinito).
- En T7, comparar θ = 1, θ = 1/2 y θ = 1/2 con 2 pasos Rannacher en una tabla de error máximo
  relativo para t ≥ 0,2 s, y proponer el θ por defecto "con base en esos números", justificado
  en un comentario del código.

## 2. Qué se hizo

- `src/discfem/timestepping.py`: nuevo. `time_grid`, `integrate` y `TransientResult`, con la
  justificación numérica del θ por defecto en el docstring del módulo.
- `src/discfem/mesh/base.py`: modificado. `Mesh.validate` admite grupos vacíos (la clave debe
  existir y la partición del borde se sigue exigiendo).
- `src/discfem/mesh/structured.py`: modificado. `tensor_mesh` acepta ia = 0 e ib = n_r − 1, de
  modo que `brake` puede cubrir toda la cara superior. T7 lo necesita.
- `tests/test_timestepping.py`: nuevo. T4, T5, T7, la malla temporal, la reutilización de LU y
  la no filtración de flujo después de t_f (25 casos). Ejecutado como script, imprime la tabla
  de T7.
- `tests/test_axisymmetric_mode.py`: modificado. T6 ahora usa `integrate` en vez de su propio
  bucle (pendiente que quedó de la Fase 2).
- `PLAN.md` (ESTADO), `reportes/INDEX.md`, este reporte.

## 3. Cómo funciona

**Del sistema semidiscreto a un sistema por paso.** La Fase 2 dejó M Ṫ + (K + H) T = F(t), un
sistema de EDO, una por nodo. Para avanzar de tⁿ a tⁿ⁺¹ = tⁿ + Δt, el esquema θ evalúa todo lo
que no es derivada temporal como un promedio ponderado entre el instante viejo y el nuevo:

[M + θΔt(K+H)] Tⁿ⁺¹ = [M − (1−θ)Δt(K+H)] Tⁿ + Δt [θ Fⁿ⁺¹ + (1−θ) Fⁿ]

- Tⁿ es conocido y Tⁿ⁺¹ es la incógnita; Δt es el paso.
- θ = 1 (Euler implícito) usa solo el instante nuevo. Es de primer orden en Δt y no oscila
  nunca (L-estable).
- θ = 1/2 (Crank–Nicolson, CN) es la regla del trapecio: segundo orden en Δt, pero amortigua
  poco los modos rápidos. Ante un salto brusco, como el flujo que arranca de golpe en t = 0,
  puede oscilar.
- Rannacher: los primeros n pasos con θ = 1 amortiguan el salto, y después se sigue con CN.

Cada paso es un sistema lineal "estacionario" como los de clase, con la matriz
A = M + θΔt(K+H). Si θ y Δt no cambian, A tampoco cambia. Por eso se factoriza una vez
(`splu`, LU dispersa) y en cada paso solo se hacen dos sustituciones triangulares, que son muy
baratas. `integrate` guarda una factorización por cada par (θ, Δt) distinto: una con Δt
constante, dos con Rannacher, y una más si el último paso es más corto.

**Por qué t_f debe estar en la malla temporal.** El flujo de fricción se apaga en t_f. Si t_f
cayera dentro de un paso, ese paso mezclaría "con flujo" y "sin flujo", y la energía que entra
dependería de Δt de forma arbitraria. `time_grid` construye el tramo [0, t_f] con exactamente
round(t_f/Δt) pasos iguales (Δt = t_f/100 por defecto), y sigue con el mismo Δt hasta t_end.
`integrate` rechaza una malla que no contenga t_f. Además, un paso que empieza en tⁿ ≥ t_f usa
F_brake = 0 en ambos extremos. Sin esto, con θ < 1, el término (1−θ)Fⁿ con Fⁿ = F(t_f) metería
energía después del frenado si el flujo no se anula en t_f (por ejemplo, un `constant_flux`).
Hay una prueba de esto.

**Balance de energía (§3.6).** Se premultiplica el esquema por 1ᵀ (sumar todas las
ecuaciones). Como 1ᵀK = 0 (Fase 2), la conducción no crea ni destruye energía y queda:

1ᵀM(Tⁿ⁺¹ − Tⁿ) = Δt[θ(1ᵀF_brakeⁿ⁺¹ − Q_convⁿ⁺¹) + (1−θ)(1ᵀF_brakeⁿ − Q_convⁿ)]

con Q_conv = 1ᵀ(H T − F_conv), que es el calor que sale por convección. `integrate` acumula
E_brake (lo que entra) y E_conv (lo que sale) con exactamente esos pesos θ, y calcula
E_stored = 1ᵀM(Tⁿ − T⁰). El residuo E_stored − (E_brake − E_conv) solo puede deberse al
redondeo y a la exactitud de la solución LU. Las energías están en J/rad, sin el 2π.

| Fórmula (PLAN.md) | Código |
|---|---|
| Esquema θ (§3.5) | `integrate`, `timestepping.py:149`; matriz A en `:236`; lado derecho en `:250` |
| θ = 1 en los primeros n pasos (Rannacher) | `timestepping.py:226` |
| Una LU por (θ, Δt), reutilizada | caché en `integrate`, `timestepping.py:220-240` |
| Δt = t_f/100, t_f múltiplo de Δt (§3.5) | `time_grid`, `timestepping.py:52`; ajuste en `:79` |
| q_fric = 0 para t > t_f; t_f en la malla (CLAUDE.md, trampa 5) | `tf_missing`, `:184`; `active`, `:242` |
| d/dt(1ᵀMT) = 1ᵀF_brake − Q_conv (§3.6) | `e_brake`, `e_conv`, `e_stored`, `timestepping.py:253-257`; `balance_residual`, `:132` |
| θ por defecto | `DEFAULT_THETA`, `DEFAULT_RANNACHER_STEPS`, `timestepping.py:45-47` y docstring del módulo |

## 4. Decisiones

- **θ por defecto = 1/2 con 2 pasos Rannacher.** La Fase 3 pide decidir con números, no por
  costumbre. Medí T7 también con Δt = t_f/100 = 0,0556 s, que es el paso por defecto real del
  caso de referencia y es 5 a 10 veces mayor que el de T7. Con ese Δt, CN + Rannacher tiene el
  menor error en todas las filas (tabla en §6). θ = 1 queda en 3,5 % y CN puro en 6 % al inicio.
  El YAML de referencia (`configs/ejemplo_frenado.yaml`, §4 del PLAN) fija `theta: 1.0` de forma
  explícita, así que no cambia nada para T8 y T9. **Esto lo decides tú:** si prefieres θ = 1
  por defecto, es cambiar una constante.
- **Caché tolerante a redondeo (error mío, corregido).** La primera versión usaba la clave
  (θ, Δt con 12 cifras). Los Δt de `linspace` difieren en unos ulps, así que con Δt = 0,005
  salían 3 factorizaciones en vez de 1. Ahora, dos pasos con el mismo θ cuyo Δt difiere menos de
  1e-9 relativo comparten la LU, y ese paso usa el Δt de la LU también en el lado derecho y en el
  balance, para seguir siendo exactamente un paso θ. Un error del mismo tipo en `time_grid` (con
  t_f = t_end no se ajustaba el último paso) también se corrigió.
- **Expectativa errónea en una prueba.** `test_factorizations_reused` esperaba 1 LU con
  t_end = 20 s, pero 20 no es múltiplo de 0,0556, así que el último paso es más corto y hacen
  falta 2. Era la prueba la que estaba mal: se cambió a t_end = 3·t_f para el caso de 1 LU, y se
  mantuvo el de 20 s esperando 2.
- **Balance por componente en vez de un solo total.** Separar E_brake y E_conv permite comparar
  con la potencia física (Fase 4: 2π ∫ q r dr = γ P_f) y ver qué fracción sale por convección.
- **Grupos de frontera vacíos.** T7 necesita `brake` sobre toda la cara superior. Descarté
  crear una geometría especial: es más simple permitir `top_inner` y `top_outer` vacíos en
  `tensor_mesh`. `Mesh.validate` sigue exigiendo que los seis grupos existan y particionen el
  borde.
- **Discrepancia con el prototipo en CN.** Aquí CN puro da 0,07–0,6 % para t ≥ 0,2 s con
  Δt = 0,005–0,01, y el PLAN reporta 1–2 %. No ajusté nada para acercarme. Una posible causa
  (no verificada) es que el prototipo arrancara con F⁰ = 0 en vez de F⁰ = q0; ese F⁰ = 0
  equivale a otra discretización del escalón y oscila más.

## 5. Cómo verificarlo

```
.venv\Scripts\python -m pytest -q
.venv\Scripts\python -m ruff check src tests
.venv\Scripts\python tests\test_timestepping.py
```

Salida real de pytest y ruff:

```
All checks passed!
.............................                                            [100%]
101 passed in 3.85s
```

Tabla T7, salida real de `python tests\test_timestepping.py`:

```
  Δt [s]  esquema             máx. t≥0,2 s  en t≈3 s (t usado)   balance
  0.0050  θ = 1                     0.395%    0.090%    3.0000   6.5e-14
  0.0050  θ = 1/2                   0.070%    0.023%    3.0000   4.6e-14
  0.0050  θ = 1/2 + 2 Rann.         0.068%    0.023%    3.0000   4.5e-14
  0.0100  θ = 1                     0.737%    0.158%    3.0000   6.4e-14
  0.0100  θ = 1/2                   0.595%    0.022%    3.0000   5.8e-14
  0.0100  θ = 1/2 + 2 Rann.         0.111%    0.023%    3.0000   5.8e-14
  0.0556  θ = 1                     3.452%    0.773%    3.0024   3.2e-14
  0.0556  θ = 1/2                   6.010%    0.502%    3.0024   3.3e-14
  0.0556  θ = 1/2 + 2 Rann.         1.482%    0.048%    3.0024   4.0e-14
```

T4 y T5 sobre la malla estructurada del caso de referencia (511 nodos; flujo q0(1 − t/t_f) con
q0 = 7,7e5 y t_f = 5,56; t_end = 20 s). Salida real de un script que reutiliza las funciones de
la prueba:

```
malla: 511 nodos; pasos: 360
conv=False theta1         cierre=7.41e-15  E_brake=1.3987e+04 J/rad  E_conv=0.0000e+00  E_stored=1.3987e+04  LU=2
conv=False CN             cierre=3.22e-15  E_brake=1.4128e+04 J/rad  E_conv=0.0000e+00  E_stored=1.4128e+04  LU=2
conv=False CN+Rannacher2  cierre=1.53e-14  E_brake=1.4125e+04 J/rad  E_conv=0.0000e+00  E_stored=1.4125e+04  LU=3
conv=True  theta1         cierre=3.90e-15  E_brake=1.3987e+04 J/rad  E_conv=2.2346e+02  E_stored=1.3763e+04  LU=2
conv=True  CN             cierre=3.73e-15  E_brake=1.4128e+04 J/rad  E_conv=2.2418e+02  E_stored=1.3904e+04  LU=2
conv=True  CN+Rannacher2  cierre=2.83e-14  E_brake=1.4125e+04 J/rad  E_conv=2.2417e+02  E_stored=1.3901e+04  LU=3
T5 max|T-T0| (θ=1): 8.881784197001252e-14
```

(LU = 2 o 3 aquí porque 20 s no es múltiplo de Δt: el último paso es más corto.)

## 6. Resultados

| Prueba | Criterio de la §6 | Medido aquí | Prototipo (§6) |
|---|---|---|---|
| T4 sin convección | cierre < 1e-10 | 3,2e-15 – 1,5e-14 | ~4e-15 |
| T4 con convección | cierre < 1e-10 | 3,7e-15 – 2,8e-14 | ~4e-15 |
| T5 | \|T − T0\| < 1e-10 | 8,9e-14 | — |
| T7, θ = 1, Δt = 0,005 | < 1 % (t ≥ 0,2 s) y < 0,2 % en t = 3 s | 0,395 % y 0,090 % | 0,5 % y 0,007–0,06 % |
| T7, θ = 1, Δt = 0,01 | ídem | 0,737 % y 0,158 % | — |
| T7, θ = 1/2 | < 3 % | 0,070 % (Δt = 0,005); 0,595 % (Δt = 0,01) | 1–2 % |
| T7, θ = 1/2 + 2 Rann. | (se le exige lo mismo que a θ = 1) | 0,068 % y 0,023 %; 0,111 % y 0,023 % | — |
| T6 (ahora con `integrate`) | < 1 % con nx = 20 | sin cambios: 0,31 % | 0,31 % |

**La energía que entra depende del esquema.** Con flujo q0(1 − t/t_f) y N = 100 pasos, la regla
θ = 1 evalúa la carga al final de cada paso: E = E_exacta·(N−1)/N, que es 1 % menos. Medido:
13 987 frente a 14 128 J/rad, y la prueba `test_t4_brake_energy_matches_time_rule` lo verifica
a 1e-12 para los tres esquemas. CN integra exacto una carga lineal en el tiempo (14 128 J/rad
= q0·t_f/2·(r_b² − r_a²)/2).

## 7. Limitaciones y riesgos

- **θ = 1 con el Δt por defecto subestima en 1 % la energía de frenado.** Además, en T7 tiene
  un 3,5 % de error temprano. El caso de referencia (T8, T9) usa θ = 1 porque así lo fija el
  YAML del PLAN. Si T8 no cae en ±0,5 °C en la Fase 5, esta es la primera causa a revisar.
- **El criterio "t = 3 s" de T7 con Δt = t_f/100** se evalúa en t = 3,0024 s, el instante de
  la malla más cercano. No se interpoló.
- **No hay Δt adaptativo ni frenadas repetidas.** Varios t_f exigirían que `time_grid` aceptara
  una lista de instantes de quiebre; hoy solo admite uno.
- **Se guardan todos los campos por defecto (`save_every = 1`).** Con 42 000 nodos y ~1100
  pasos serían ~370 MB. Para corridas grandes hay que subir `save_every`; la Fase 5 decidirá qué
  se guarda.
- **La discrepancia de CN con el prototipo (§4) no está explicada con certeza.**
- **El balance usa `dt` del caché.** Si dos pasos difieren en más de 1e-9 relativo, cada uno
  usa su propio Δt, y eso es lo correcto. Si una malla temporal de un usuario tuviera pasos que
  difieren en menos de eso pero de forma intencional (algo improbable), se tratarían como
  iguales.

## 8. Para revisar y para el informe

1. **El esquema θ convierte el problema transitorio en una sucesión de problemas
   "estacionarios"** con la matriz M + θΔt(K+H). Por eso se factoriza una vez y se reutiliza.
   Qué cambia entre θ = 1 y θ = 1/2: orden de precisión contra amortiguamiento.
2. **Por qué CN oscila con el flujo escalón y cómo lo arregla Rannacher.** Se ve en la fila
   Δt = 0,0556 de la tabla: 6,0 % → 1,5 %.
3. **Por qué t_f debe ser un nodo de la malla temporal** y por qué un paso que empieza en t_f no
   debe usar F(t_f) (trampa 5). `test_constant_flux_does_not_leak_after_tf` lo demuestra.
4. **El balance de energía discreto cierra a precisión de máquina** porque 1ᵀK = 0 y porque las
   energías se integran con los mismos pesos θ. Sirve como control de calidad de cada corrida,
   no como prueba de exactitud (θ = 1 cierra a 1e-15 aunque mete un 1 % menos de energía).
5. **Para el informe:** la tabla T7 (Δt × esquema) es la validación temporal y la justificación
   del esquema elegido. Una figura sugerida es T_s(t) exacta frente a los tres esquemas con
   Δt = t_f/100, ampliada en 0–1 s para mostrar la oscilación de CN; aún no está generada.
