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

## Resultados verificados
(ninguno todavía)

## Siguientes pasos (prioridad)
1. Ver resultados de `runs/qtx1` (dobles planos, j=2). Si 0: caza directa (j=0)
   en GF(8)/GF(9)/GF(4) con m hasta 13 (los grupos de Aydin aún encuentran
   récords QT en F4, F8, F9).
2. Extensiones con j≥3 columnas (búsqueda greedy/recocido de V k×j).
3. Familias más ricas: multi-twisted (λ distinto por bloque), 2-generadores.
4. Parsear construcciones de codetables (QC + punct/short) para extender los
   propios códigos récord de la tabla.
