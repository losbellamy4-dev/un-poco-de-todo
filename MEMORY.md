# MEMORY — estado del laboratorio

Última actualización: 2026-10-12 (sesión 1, en curso)

## Objetivo elegido y por qué
Sesión 1 partió de un repo vacío. Criterio: resultados **verificables
objetivamente**, alcanzables con 4 CPUs y sin API de LLM, con valor real
(publicables / acreditables). Elegido: **encontrar códigos lineales que mejoren
las cotas inferiores de las tablas de Grassl (codetables.de)**.
- Verificación exacta e independiente (motor C propio + GAP/GUAVA).
- ~107 000 huecos abiertos (cota inferior < superior) en q ∈ {2,3,4,5,7,8,9}.
- Un récord se envía a Markus Grassl (codetables.de) y queda acreditado.

Alternativas en cola: esquemas de multiplicación de matrices (flip graphs,
catálogo de Sedoglavic), covering designs (LJCR no responde desde aquí),
covering arrays (tablas de Colbourn tras login), códigos esféricos.

## Infraestructura (codes/) — ver codes/README.md
- `data/bounds_q{q}.tsv`: tablas completas descargadas 2026-10-12 (`scripts/scrape_bounds.py`).
- `src/mindist.c` → `libmindist.so` (+ `make` en src/):
  - `md_compute`: distancia mínima Brouwer–Zimmermann, OpenMP, ~2e8 palabras/s (4 hilos).
  - `md_minwords` / `md_words_range`: lista proyectiva de palabras de peso en [d, dhi].
  - `md_extend` / `md_extend2`: columnas extra que suben d en 1 o 2.
- `src/qtsearch.c`: muestreo aleatorio QT (modos histograma/extensión).
- `src/qtx.c`: búsqueda local QT con objetivo = log nº esperado de extensiones con j columnas.
- `scripts/qt_driver.py`, `scripts/qtx_driver.py`: planificación + ejecución paralela.
- `scripts/verify.py`: verificación exacta (C) + GAP/GUAVA → `results/records.jsonl`.
- `scripts/impact.py`: qué entradas de la tabla mejora un código (acortar/perforar/subcódigo).
- GAP: `apt-get install -y gap gap-guava` (no viene instalado; reinstalar cada sesión).

## Lecciones (evidencia)
- Códigos aleatorios: 4–6 por debajo de la tabla. QT de 1 generador igualan el récord
  con frecuencia 0.1–10 % en GF(7) n≤70 pero **no lo superan** (~2M+ muestras).
- Extensión por una columna funciona muy bien en códigos QT (palabras mínimas en
  pocas órbitas de tamaño l), p.ej. [58,8,39]_7 → [59,8,40]_7 en el 80 % de intentos,
  pero en esos puntos la tabla ya iguala. En los puntos "planos" (lb(n+1)=lb(n))
  los códigos que igualan tienen 300–600 palabras mínimas: demasiadas.
- La búsqueda local (qtclimb/qtx) encuentra códigos con d = lb−1 y sólo 1 órbita de
  palabras mínimas (p.ej. [32,11,15]_7 con A_15 = 16).
- `pkill -f "patrón"` mata tu propio shell si el patrón aparece en el comando:
  usar `pkill -f "^python3 qt_driver"`.

## Resultados de las búsquedas clásicas (sesión 1) — todo negativo
- QT aleatorio/escalada en GF(3..9): iguala récords, nunca los supera.
- Extensión 1/2 columnas y Construction X (r≤4, aux MDS) sobre códigos QT: 0 éxitos.
- Construction X sobre los propios códigos récord de la tabla (GF(7), 287 reconstruidos,
  `xtable.py`): 0 éxitos → Grassl ya aplica X sistemáticamente. Mercado eficiente.

## Giro: códigos cuánticos (tabla QECC de qubits, `data/qbounds_q4.tsv`)
- Cotas inferiores débiles en n≥35: mesetas largas por relleno ("ExtendCode" de un
  código menor), p.ej. [[40,8,8]] = [[34,8,8]] + 6 qubits; [[41,1]]: lb 11 en 12 pasos.
- Motor `src/qmd.c` (`qmd_distance`): distancia simpléctica de códigos aditivos GF(4),
  n≤64, máscaras de 64 bits, B-Z sobre posiciones de qubit, exclusión de subcódigo
  (distancia impura N(S)\S). Validado contra fuerza bruta (300 códigos), [[5,1,3]], Steane.
- Familias: grafos circulantes/bicirculantes ([[n,0]], `qsearch.py`) y "aditivos
  cíclicos" S={(c, c·s): c∈<g>} con s palindrómico, deg g = k (`qcyclic.py`, `qbatch.py`).
  En n=52 los grafos circulantes aleatorios igualan el récord (d=14) en ~3 % de intentos.

## Resultados verificados
(ninguno todavía)

## Siguientes pasos (prioridad)
1. Cuántico: resultados de `runs/qs_*` (grafos k=0) y `runs/qb1.log` (cíclicos k>0).
   Añadir escalada (d, -A_d) y familias de grafos multi-circulantes.
2. Verificador independiente para códigos cuánticos (segunda implementación).
3. Clásico (baja prioridad): multi-twisted, 2-generadores; parser de construcciones
   ampliado (QR codes, SubfieldSubcode, Construction B).
