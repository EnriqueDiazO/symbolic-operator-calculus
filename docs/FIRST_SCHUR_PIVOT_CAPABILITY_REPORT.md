# Informe de capacidad: primer pivote de Schur normalizado

Fecha: 2026-07-20.

Este informe describe lo que certifica la implementación formal, lo que se
introduce como regla suministrada y lo que permanece como problema analítico.
No contiene una afirmación de Fredholmness ni de pertenencia a una álgebra de
operadores.

## Resultado principal

La paquetería construye de manera determinista, en su AST no conmutativo, la
identidad formal

```tex
N_{2}^{(1)}
= N_{2}
- Z_{2}^{-1} A_{2,1} T_{1,-} R_{1} Z_{1}^{-1} A_{1,2} T_{2,-}.
```

También registra, como relación distinta y no certificada analíticamente, la
equivalencia módulo compactos

```tex
N_{2}^{(1)}
\simeq N_{2}
+ Z_{2}^{-1}
  (\widetilde V_{\alpha_2}-G_2) W_{2,1}^{-}
  T_{1,-} R_1 Z_1^{-1}
  (\widetilde V_{\alpha_1}-G_1) W_{1,2}^{+} T_{2,-}.
```

`R_1` permanece como un átomo `formal_regularizer`; no se calcula ni se
expande un símbolo inverso.

## Cálculos certificados formalmente

La certificación disponible es **estructural dentro del cálculo formal**, no
analítica. La implementación verifica:

- la construcción ordenada de
  `A22^(1) = A22 - A21 A11^(-1) A12`;
- la distribución finita de `Z2^(-1)` a la izquierda y `T2,-` a la derecha;
- la sustitución por subpalabra contigua, sin conmutar ni reordenar factores;
- el signo negativo de la corrección exacta;
- el reconocimiento de la palabra `Z2^(-1) A22 T2,-` por el átomo `N2`;
- el signo positivo que resulta al introducir la regla con signo negativo para
  `A21`;
- la expansión opcional y exclusiva de las dos diferencias lineales, con
  coeficientes `(+, -, -, +)`;
- la reconstrucción exacta de la expansión de cuatro términos desde el
  metadato factorizado;
- que cada término de la corrección contiene una sola ocurrencia atómica de
  `R1`;
- que las salidas y su LaTeX son deterministas.

Estas comprobaciones no verifican dominios, acotación, compactness, cierre de
álgebras ni propiedades de Fredholm.

## Transformaciones que son reglas suministradas

Las siguientes relaciones se almacenan explícitamente como hipótesis o reglas;
la paquetería no demuestra su contenido analítico:

1. `A11^(-1) = T1,- R1 Z1^(-1)` es una `FORMAL_SUBSTITUTION`.
2. `Z2^(-1) A22 T2,- = N2` es la regla suministrada que define/reconoce el
   término diagonal normalizado.
3. `A21 ≃ -(Vtilde_alpha2-G2) Wminus_21` es una
   `MOD_COMPACT_EQUIVALENCE` sin certificación interna.
4. `A12 ≃ (Vtilde_alpha1-G1) Wplus_12` es una
   `MOD_COMPACT_EQUIVALENCE` sin certificación interna.
5. Propagar las dos últimas reglas al producto completo produce una declaración
   módulo compactos; no constituye una prueba automática de estabilidad del
   ideal compacto bajo todas las composiciones mostradas.

La traza utiliza tipos distintos para `EXACT_EQUALITY`,
`FORMAL_SUBSTITUTION`, `MOD_COMPACT_EQUIVALENCE` y
`ANALYTIC_PROOF_OBLIGATION`.

## Clasificación declarada de factores

La tabla es metadato suministrado. Una fila no prueba pertenencia a la clase
indicada.

| Factor | Clase operatorial conocida | Estatus | Fuente matemática |
|---|---|---|---|
| $Z_2^{-1}$ | MultiplicationOperator | exacto | definición del transporte |
| $\widetilde V_{\alpha_2}$ | Multiplication × Dilation | exacto | transporte del shift |
| $G_2$ | MultiplicationOperator | exacto | coeficiente branchwise |
| $W_{2,1}^{-}$ | LocalizedWienerHopfOperator | módulo compactos | localización cuspídea 2025 |
| $T_{1,-}$ | AuxiliaryShiftOperator | invertible | resultado branchwise |
| $R_1$ | MellinPDORegularizer | regularizador módulo compactos | KKL 2014, Thm. 5.8 |
| $Z_1^{-1}$ | MultiplicationOperator | exacto | normalización diagonal |
| $\widetilde V_{\alpha_1}$ | Multiplication × Dilation | exacto | transporte del shift |
| $G_1$ | MultiplicationOperator | exacto | coeficiente branchwise |
| $W_{1,2}^{+}$ | LocalizedWienerHopfOperator | módulo compactos | localización cuspídea 2025 |
| $T_{2,-}$ | AuxiliaryShiftOperator | invertible | resultado branchwise |

## Afirmaciones analíticas abiertas

1. Justificar la sustitución Wiener--Hopf módulo compactos.
2. Demostrar que la corrección completa pertenece a una álgebra cerrada
   apropiada.
3. Demostrar cierre bajo los productos que aparecen.
4. Decidir si la salida es un Mellin PDO admisible o un elemento de la álgebra
   cuspídea de 2025.
5. Solo después, aplicar un criterio de Fredholmness sin presuponer su
   conclusión.

En particular, no están demostradas la compactness del residuo, la acotación de
cada composición, la admisibilidad de un símbolo Mellin, la pertenencia a la
álgebra cuspídea ni la Fredholmness de `N2^(1)`.

## Siguiente expresión que debe analizarse

El objeto analítico inmediato es la corrección completa

```tex
\mathcal C_2
:=
Z_2^{-1}
(\widetilde V_{\alpha_2}-G_2)
W_{2,1}^{-}
T_{1,-}
R_1
Z_1^{-1}
(\widetilde V_{\alpha_1}-G_1)
W_{1,2}^{+}
T_{2,-}.
```

El siguiente cálculo debe determinar un espacio común, el ideal compacto
pertinente y una propiedad de cierre que permita clasificar `\mathcal C_2`.
Solo entonces tendría sentido analizar `N_2+\mathcal C_2` con un criterio de
Fredholm.

## Recomendación

**Conclusión: C) todavía no puede decidirse.**

La traza actual no contiene evidencia suficiente para elegir rigurosamente
entre:

- A) reducir toda la corrección a un Mellin PDO admisible; o
- B) trabajar directamente en la álgebra cuspídea de 2025.

La presencia simultánea de dilataciones transportadas, dos factores
Wiener--Hopf localizados, shifts auxiliares y un regularizador Mellin hace
plausibles ambos marcos, pero la paquetería no dispone de un teorema de cierre
que decida cuál absorbe el producto completo. La recomendación operativa es
analizar primero `\mathcal C_2` en ambos lenguajes y escoger A o B únicamente
cuando se haya demostrado una representación o un cierre explícito.

## Artefactos reproducibles

- `notebooks/first_schur_pivot_normalization.ipynb` ejecuta la traza completa y
  muestra las dos salidas principales, la tabla y las obligaciones.
- `docs/FIRST_SCHUR_PIVOT_DERIVATION.tex` es un fragmento sin preámbulo y sin
  macros privadas.
- `docs/PHASE_N_AUDIT_FIRST_SCHUR_PIVOT.md` documenta la capacidad previa y la
  ampliación mínima justificada.
