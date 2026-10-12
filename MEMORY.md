# MEMORY — estado del laboratorio

Última actualización: 2026-10-12 (fin de la sesión 1). Lee también `journal/`.

## Objetivo y criterio
Resultados **verificables objetivamente**, alcanzables con 4 CPUs y sin API de LLM.
Sesión 1 eligió: superar cotas inferiores de las tablas de códigos de Grassl
(codetables.de: lineales [n,k,d]_q y cuánticos [[n,k,d]]). Resultado: **0 récords**,
infraestructura sólida y un benchmark publicable (`bench/`).

## Lección principal (evidencia de ~10 CPU-horas)
Donde verificar d es barato (≤1e10 palabras), las tablas son un mercado eficiente:
QT aleatorio, escalada, extensiones, Construction X, grafos circulantes... igualan
los récords pero no los superan. Donde la tabla es débil (n grande), probar d por
enumeración es inviable (1e15+). ⇒ La próxima ventaja debe venir de **cotas de
distancia demostrables sin enumerar** o de familias algebraicas nuevas, no de más
muestreo.

## Infraestructura (ver codes/README.md y bench/README.md)
- `codes/src/mindist.c` (+`qmd.c`) → `libmindist.so`; `make` en codes/src.
  - `md_compute` B-Z GF(q≤9) ~2e8 palabras/s; `md_isd` rechazo probabilístico;
    `md_minwords`/`md_words_range`; `md_extend`/`md_extend2`/`md_avoid` (Construction X).
  - `qmd_distance`: distancia simpléctica de códigos aditivos GF(4), n≤64, con exclusión
    del estabilizador (distancia impura).
- Búsquedas: `qtsearch.c`, `qtclimb.c`, `qtx.c`, `qtX.c` + `qt_driver.py`, `qtx_driver.py`;
  cuánticas `qsearch.py` (grafos), `qcyclic.py`/`qbatch.py` (aditivos cíclicos),
  `qqc.py` (QC hermíticos autoortogonales).
- Tablas: `codes/data/bounds_q*.tsv`, `codes/data/qbounds_q4.tsv`; construcciones
  cacheadas en `codes/data/constructions/` y `qconstructions/`; `table_codes.py`
  reconstruye los códigos récord desde su descripción Magma.
- Verificación: `verify.py` (C + GAP), `qverify.py` (qmd + imagen binaria [3n] + GAP),
  `bench/verify_submission.py`.
- GAP: `apt-get install -y gap gap-guava` (reinstalar cada sesión; ~2 min).

## Hechos útiles
- Códigos QT: palabras mínimas en órbitas de tamaño l ⇒ extensiones fáciles, pero
  sólo donde la tabla ya iguala.
- Tablas clásicas cerradas bajo suma directa, (u|u+v) y Construction X sobre sus
  propios códigos (comprobado GF(7)).
- Tabla cuántica k=0, n=117–124: débil (18–21) pero ya incluye la Construction X
  cuántica del mejor cíclico; mejorar exige d≥20 demostrable en [119,60]_4.
- Entradas cuánticas con "Construction type: GV" no tienen código explícito.
- Red: LJCR (proxy 502) y coveringrepository.com (Cloudflare) bloqueados; Zenodo,
  codetables.de, fmm.univ-lille.fr accesibles.
- `pkill -f patrón` mata tu shell si el patrón está en el comando: usar `^python3 x`.
- Los contadores "matched" en modo target son engañosos: en modo target la función
  devuelve la primera palabra < T, no el mínimo.

- **Tope estructural** (aprendido tarde): si C ⊆ g·R^m con <g> hermítico autoortogonal,
  C^⊥h ⊇ (<g>^⊥h)^m ⇒ d ≤ d(<g>^⊥h) (código cíclico corto). `qqc.py`/`qqc2.py` y `qcx.py`
  sufren este tope. Análogo: familias "grafo + subcódigo" están topadas por la distancia
  del grafo. Antes de lanzar una familia, buscar el subcódigo del dual que la limita.
- Familia sin ese tope: doble circulante unitario {(c, c·u): c∈<g>}, u = b/b̄* (`qdc.py`);
  su dual contiene el autodual {(y, y·u)} ⇒ tope = distancia de ese código autodual.
- Motor cuántico 128 qubits: `qmd128.c` (`qdist128`, `quantum_distance_any`).
- `bch_scan.py`: todos los códigos constacíclicos tipo BCH (todo n, λ, paso b) con
  distancia *certificada* por la cota BCH: 15 773 pares (n,k) en q=3,4,5,7,8,9 → 0 superan
  la tabla. Grassl ya domina las familias algebraicas clásicas.
- `src/arcsearch.c` (tabú sobre conjuntos de puntos en PG(k−1,q)): ni siquiera alcanza
  los códigos de Kohnert ([57,5,45]_8). Para k pequeño, Kohnert (grupos de automorfismos
  prescritos + reducción de retículos) es muy superior a la búsqueda local.

## Resultados verificados
(ninguno)

## Siguientes pasos (prioridad)
1. Revisar `runs/qtX2*`, `runs/qb1.log`, `runs/qqc1*` (búsquedas que siguieron corriendo).
2. Cotas demostrables más finas que BCH (Hartmann–Tzeng, Roos, espectrales para QC) en n
   grande; o el método de Kohnert (automorfismos prescritos + LLL/ILP) implementado aquí
   para grupos no probados (es el que produce los récords de k pequeño).
3. Publicar/difundir `bench/` (artefacto o repo) para sistemas tipo AlphaEvolve.
4. Dominios alternativos solo si hay acceso a SOTA actual (covering designs vía Zenodo
   snapshot es arriesgado: puede estar desactualizado).
