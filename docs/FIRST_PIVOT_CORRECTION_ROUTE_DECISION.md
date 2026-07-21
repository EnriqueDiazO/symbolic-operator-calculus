# Decisión de ruta para la corrección completa

## Resultado

**C: insufficient evidence**

La decisión usa únicamente la matriz de interfaces y las fuentes verificadas.
No expresa una preferencia heurística y no afirma compactitud, pertenencia a una
álgebra, admisibilidad Mellin ni Fredholmness.

## Ruta A — Mellin PDO

| Requisito mínimo | Evidencia actual | Estado |
|---|---|---|
| `C_2^{norm}=Op(c)+K` | no existe construcción de `c` ni resto `K` | no demostrado |
| `c` en una clase Mellin admisible | `R_1` está clasificado como regularizador, pero la palabra completa no | no demostrado |
| reglas de semiproducto | `Cauchy × R_1`, `R_1 × multiplicación` y `R_1 × dilatación` carecen de regla aplicable | faltante central |
| control de `R_1` con factores laterales | `R_1` permanece atómico y no puede moverse | no demostrado |
| no degeneración del símbolo | no hay símbolo de la corrección | no aplicable todavía |

KKL 2014, Thm. 5.8, solo puede entrar después de probar pertenencia de símbolo y
sus condiciones de borde/fibra. Los trabajos 2008 y 2024 están localizados, pero
el índice aún no contiene un teorema verificado que cierre las interfaces
necesarias. Por tanto la salida requerida
`A: supported by verified closure rules` no está sustentada.

## Ruta B — Álgebra cuspídea de 2025

| Requisito mínimo | Evidencia actual | Estado |
|---|---|---|
| todos los factores son generadores de la álgebra | solo hay etiquetas semánticas; `U_k`, `\widehat G_k`, `R_1` y los bloques concretos no han sido identificados como generadores | no demostrado |
| cierre del producto completo | Thm. 6.2 exige primero pertenencia a `A` | no demostrado |
| ideal compacto identificado | Lemas 3.1–3.2 y Cor. 3.3 verifican modelos explícitos; falta la identificación con los bloques de Fase N | adaptación abierta |
| homomorfismo de símbolos aplicable | verificado para la álgebra de la fuente | hipótesis no transferidas |
| criterio de Fredholm disponible | Thm. 6.2 lo da para elementos de `A` | fuera de alcance hasta probar pertenencia |
| compatibilidad con shifts y `R_1` | no existe una regla verificada que inserte el regularizador Mellin en esa álgebra | faltante central |

Por tanto la salida requerida
`B: supported by verified algebra membership` tampoco está sustentada.

## Por qué se selecciona C

La expansión C3 muestra que en cada uno de los 16 términos aparecen las mismas
interfaces críticas: un Wiener--Hopf localizado con una pieza auxiliar, una
pieza de Cauchy con `R_1`, y `R_1` con un factor derecho normalizado. No existe
una cadena certificada que atraviese ese núcleo. Las fuentes verificadas ofrecen
criterios **condicionados a** membresía, pero no prueban esa membresía para la
corrección.

La evidencia es por ello insuficiente para decidir entre A y B sin un nuevo
lema. El resultado formal sí permite mantener ambas rutas como candidatas, pero
ninguna está apoyada por todas las reglas de cierre requeridas.

## Próximo lema matemático sugerido

El siguiente objetivo debería ser un **lema de cierre del núcleo mixto**, con
hipótesis explícitas sobre pesos, dilataciones, coeficientes y localización:

> Para cada elección `Q_1 in {U_1^{-1}P^+, P^-}` y
> `B_1 in {U_1, \widehat G_1}`, clasificar rigurosamente
> `Q_1 R_1 B_1 W_{1,2}^{+}` módulo compactos, sin reordenar factores.

Una versión útil debe concluir exactamente una de estas alternativas:

1. `Op(c)+K` con `c` en una clase Mellin identificada y con estimaciones de
   semiproducto; o
2. pertenencia demostrada a la álgebra cuspídea (o a una extensión explícita)
   compatible con el homomorfismo de símbolos.

Este núcleo aparece en los 16 términos y separa directamente las rutas A y B.
Solo después tendría sentido acoplar el bloque izquierdo y estudiar condiciones
de no degeneración; no se anticipa una conclusión de Fredholmness.
