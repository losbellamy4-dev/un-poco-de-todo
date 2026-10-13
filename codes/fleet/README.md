# Flota de sesiones paralelas (sesión 2)

La sesión coordinadora reparte el trabajo en sesiones cloud independientes (4 CPUs cada una).
Las sesiones de la flota **no hacen commit ni push**: informan a la coordinadora con
`send_message`, y ella verifica e integra los resultados en este repositorio.

## Preparación (cada sesión)
```bash
cd codes/src && make -s libmindist.so && cd ../..
python3 -c "import highspy" || pip install -q highspy
# GAP solo si hace falta verificación cruzada: apt-get install -y gap gap-guava
```

## Fragmentos de cómputo (s01–s30)
`codes/fleet/shards/sNN.txt`: las 1612 metas `q k n d` (d = cota inferior de la tabla + 1,
tomada de `codes/data/bounds_q*.tsv`), repartidas en round-robin.
```bash
nohup codes/fleet/run_shard.sh s07 500 4 > runs/s07.out 2>&1 &
```
Se tarda unas 2 h. Cada línea de `runs/fleet_s07.log` es `(q, k, n, d) <grupos probados> <salida>`.
Un `*** FOUND [n,k,d]_q` significa que se encontró un código con d ≥ cota de la tabla + 1 (la
distancia la recalcula el motor C). La matriz generadora queda en `codes/results/candidates.jsonl`.

## Formato del informe (send_message a la coordinadora)
```
FLEET-REPORT <id>
status: done | partial (<hechos>/<total>)
found: <número>   (0 si no hay)
<por cada hallazgo: q n k d, y la línea JSON completa de candidates.jsonl>
notes: <≤150 palabras: qué se hizo, tiempos, qué aprendiste>
```
Antes de informar de un hallazgo, recalcula d con
`python3 codes/scripts/verify.py` o `gf.min_distance` y compáralo con
`codes/data/bounds_q{q}.tsv` (columna lower). Ojo: las líneas de calibración
`[9,3,7]_8` y `[60,5,48]_8` de candidates.jsonl son antiguas y NO son hallazgos.
