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

## Resultados verificados
(ninguno)

## Siguientes pasos (prioridad)
1. Revisar `runs/qtX2*`, `runs/qb1.log`, `runs/qqc1*` (búsquedas que siguieron corriendo).
2. **Cotas demostrables** para códigos QC/QT grandes (Semenov–Trifonov / spectral bounds,
   BCH/HT/Roos para constacíclicos en GF(4), GF(8), GF(9)): buscar códigos cuya cota
   *probada* supere entradas débiles de n grande (cadenas de perforación).
3. Publicar/difundir `bench/` (artefacto o repo) para sistemas tipo AlphaEvolve.
4. Dominios alternativos solo si hay acceso a SOTA actual (covering designs vía Zenodo
   snapshot es arriesgado: puede estar desactualizado).
