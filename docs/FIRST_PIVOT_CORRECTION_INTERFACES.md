# Matriz de interfaces de la corrección completa

## Alcance

La tabla cuenta interfaces **contiguas** en los 16 términos C3. Un conteo no
nulo solo indica presencia estructural. Ninguna fila se promueve a regla por la
clase declarada de sus factores.


| Interfaz | Nº de términos | Regla existente | Fuente | Estado |
|---|---:|---|---|---|
| multiplicación × Wiener--Hopf | 16 | no implemented semiproduct rule | paper exact normalization; closure rule absent | SOURCE_VERIFIED_RULE_MISSING |
| dilatación × Wiener--Hopf | 16 | no implemented closure rule | paper exact normalization; closure rule absent | SOURCE_VERIFIED_RULE_MISSING |
| Wiener--Hopf × dilatación | 16 | exact auxiliary expansion only | T_{k,-} definition; semiproduct rule absent | SOURCE_VERIFIED_RULE_MISSING |
| Wiener--Hopf × Cauchy factor | 16 | none | no rule found in the implemented API | NO_RULE_FOUND |
| Cauchy factor × R_1 | 16 | none verified | KKL 2014 source verification pending | SOURCE_UNVERIFIED |
| Wiener--Hopf × R_1 | 0 | not contiguous in C3 | structural inspection | NOT_APPLICABLE |
| R_1 × multiplicación | 8 | none verified | Mellin semiproduct source verification pending | SOURCE_UNVERIFIED |
| R_1 × dilatación | 8 | none verified | Mellin-dilation source verification pending | SOURCE_UNVERIFIED |
| R_1 × Wiener--Hopf | 0 | not contiguous in C3 | structural inspection | NOT_APPLICABLE |
| Wiener--Hopf × auxiliar | 8 | T_{k,-}=U_k^{-1}P^+ + P^- | exact supplied auxiliary identity; interface rule absent | SOURCE_VERIFIED_RULE_MISSING |
| productos con cutoffs/localización explícitos | 0 | no explicit cutoff atom in C3 | localized status is metadata on W factors | NOT_APPLICABLE |


## Interfaces abiertas más recurrentes

Todas estas interfaces aparecen en 16 de 16 términos y constituyen el cuello de
botella estructural:

1. factor izquierdo normalizado × Wiener--Hopf localizado;
2. Wiener--Hopf izquierdo × primera pieza del auxiliar;
3. factor de Cauchy izquierdo × `R_1`;
4. `R_1` × factor derecho normalizado;
5. factor derecho normalizado × Wiener--Hopf localizado;
6. Wiener--Hopf derecho × última pieza del auxiliar.

Las filas `R_1 × multiplicación` y `R_1 × dilatación` se reparten 8/8 según el
factor derecho sea `\widehat G_1` o `U_1`. Las interfaces directas
Wiener--Hopf × `R_1` y `R_1` × Wiener--Hopf no ocurren en C3, de modo que no
pueden usarse para saltar factores intermedios.

## Motifs que sí se reconocen


| Motif | Estado | Frecuencia |
|---|---|---:|
| `left_auxiliary_piece` | EXACT_RECOGNIZED | 16 |
| `mellin_regularizer_core` | UNRESOLVED | 16 |
| `mixed_cauchy_dilation` | EXACT_RECOGNIZED | 16 |
| `normalized_left_shift` | SUPPLIED_MOD_COMPACT | 16 |
| `normalized_right_shift` | SUPPLIED_MOD_COMPACT | 16 |
| `right_auxiliary_piece` | EXACT_RECOGNIZED | 16 |


El motif `mellin_regularizer_core` permanece `UNRESOLVED`: únicamente reconoce
la posición de `R_1` y sus vecinos. Los motivos `normalized_*_shift` heredan el
nivel módulo compactos del modelo de Fase N. Las piezas auxiliares y los
subproductos contiguos `U_k^{-1}P^+` son exactos por las definiciones
suministradas.

## Obligaciones inducidas

- Probar reglas para multiplicación/dilatación junto a los Wiener--Hopf
  localizados, con hipótesis y fuente verificadas.
- Probar cómo cada pieza de Cauchy a la izquierda actúa sobre `R_1`.
- Probar cómo `R_1` se compone a la derecha con `U_1` y `\widehat G_1`.
- Identificar cutoffs/localizaciones explícitos antes de aplicar los resultados
  cuspídeos de 2025; la metadata `LocalizedWienerHopf*` no los suministra.
- Probar cierre de la palabra completa en una clase Mellin o en una álgebra
  cuspídea antes de invocar cualquier criterio de Fredholmness.
