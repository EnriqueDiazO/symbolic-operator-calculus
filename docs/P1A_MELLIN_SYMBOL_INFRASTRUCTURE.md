# P1-A: infraestructura tipada para símbolos escalares de Mellin

P1-A comienza con declaraciones explícitas de papeles y dominios. Un
`MellinVariable` combina un `sp.Symbol` ya existente con uno de los papeles
`OUTPUT_SPATIAL`, `INPUT_SPATIAL`, `RELATIVE_MULTIPLICATIVE`, `FREQUENCY` o
`PARAMETER`. El nombre del símbolo nunca determina ese papel y las factories no
aceptan strings ni crean símbolos implícitamente.

`MellinSymbolDomain` exige exactamente una frecuencia y compone, sin
proyectarlos a expresiones desnudas, los objetos P0-B `AssumptionContext`,
`ComplexDomain` y `SingularSet`. Sus operaciones solo añaden hipótesis,
singularidades, exclusiones o intersecciones. Un resultado `UNDETERMINED` no se
promueve a compatibilidad demostrada y evidencia externa solo se registra como
`EVIDENCE_SUPPLIED`.

Esta infraestructura describe datos escalares. No construye una transformada,
integral, convolución, pseudodiferencial u operador de Mellin; tampoco afirma
acotación, pertenencia a una clase funcional, compacidad o Fredholmness.
