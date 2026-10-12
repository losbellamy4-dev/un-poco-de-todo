# codes — caza de récords de códigos lineales [n,k,d]_q

Objetivo: encontrar códigos lineales con distancia mínima mayor que la cota
inferior publicada en las tablas de Markus Grassl (https://www.codetables.de).

## Uso rápido
```bash
cd codes/src && make                       # libmindist.so, qtsearch, qtclimb, qtx
apt-get install -y gap gap-guava           # verificador independiente
cd ../scripts
python3 test_mindist.py && python3 test_ext.py && python3 test_ext2.py   # tests
python3 qtx_driver.py --q 7 --slack -1 --secs 30 --out ../../runs/x.jsonl
python3 verify.py                          # verifica results/candidates.jsonl
python3 impact.py 7 58 8 40                # entradas de la tabla que mejoraría
```

## Método
1. **Distancia mínima**: Brouwer–Zimmermann con cadena de conjuntos de
   información disjuntos; modo `target` para rechazo rápido. Cuerpos GF(q) con
   polinomios de Conway (mismo que Magma/GAP).
2. **Familias**: códigos cuasi-twisted de 1 generador (g f_1, …, g f_m) en
   (F[x]/(x^l − λ))^m, g | x^l − λ, k = l − deg g, f_1 = 1.
3. **Extensión**: si las palabras de peso bajo de un código son pocas, añadir j
   columnas V (k×j) puede subir d en j. Objetivo de búsqueda local:
   F = j·k·ln q + Σ_w A_w · ln P(Bin(j, 1−1/q) ≥ T−w).
4. **Verificación**: distancia exacta con el motor C y, por separado, con
   GAP/GUAVA; sólo entonces se registra en `results/records.jsonl`.

## Ficheros
- `data/bounds_q{q}.tsv` — tablas (n, k, lower, upper) a 2026-10-12.
- `results/candidates.jsonl` — salidas brutas de búsqueda (sin verificar).
- `results/records.jsonl` — récords verificados (matriz generadora incluida).
