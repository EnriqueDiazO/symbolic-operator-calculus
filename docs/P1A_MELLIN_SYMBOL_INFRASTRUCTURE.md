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

## Expresiones y símbolos

`MellinExpression` valida que cada símbolo libre tenga un papel declarado,
que los parámetros no se presenten como variables y que las variables mudas
se registren explícitamente. Suma, producto, negación, sustitución,
diferenciación escalar y conjugación formal conservan el dominio y combinan
hipótesis, singularidades y evidencia usando el estatus más débil. `.as_expr()`
es una proyección explícita que no permite reconstruir esos metadatos.

`MellinSymbol` es un tipo separado que clasifica dependencias por papeles y
símbolos libres: `FREQUENCY_ONLY`, `SPACE_FREQUENCY`,
`RELATIVE_FREQUENCY`, `PARAMETRIC_CONSTANT` o `CONSTANT`. La igualdad
estructural, la consulta escalar de SymPy, `ConditionalIdentity` y
`ApproximateEquality` permanecen operaciones distintas.

Los ejemplos incluidos registran las condiciones de
`coth(pi*(lambda+i*kappa))`, el cociente exponencial sobre `sinh`, el
multiplicador `gamma**(i*lambda)` con `gamma > 0`, un símbolo dependiente del
espacio y una expresión relativa. Las dos identidades diagonales son
`ConditionalIdentity` escalares; una comprobación de residual cero no las
convierte en identidades operatoriales.

El renderizado usa las etiquetas “Mellin scalar symbol”, “Formal Mellin
expression”, “Conditional scalar identity”, “Evidence supplied” y
“Symbolically checked under assumptions”. No usa “Mellin operator”, “Fredholm
symbol” ni “Exact operator identity”.
