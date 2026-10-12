# Instrucciones para sesiones futuras (léelo primero)

Este repositorio es un laboratorio de investigación autónoma. Cada sesión debe
continuar el trabajo anterior, no empezar de cero.

1. **Lee `MEMORY.md`** antes de hacer nada: estado actual, resultados verificados,
   lecciones aprendidas y la cola priorizada de siguientes pasos.
2. Revisa la última entrada de `journal/` para el detalle de la sesión previa.
3. Al terminar (o cada ~hora de trabajo), actualiza `MEMORY.md` y añade/actualiza
   tu entrada en `journal/`, y haz commit + push. El contenedor es efímero:
   lo que no está en git se pierde.

## Convenciones
- Todo resultado "extraordinario" se verifica con un método independiente antes
  de registrarse en `codes/results/` (ver `codes/README.md`).
- Los cálculos largos se lanzan en segundo plano y escriben logs en `runs/`
  (ignorado por git salvo resúmenes).
- No incluir identificadores de modelo en commits ni ficheros.

## Proyectos
- `codes/` — búsqueda de nuevos récords de códigos lineales [n,k,d]_q frente a
  las tablas de Grassl (codetables.de). Motor C de distancia mínima
  (Brouwer–Zimmermann) + búsqueda de códigos cuasi-cíclicos/cuasi-twisted.
