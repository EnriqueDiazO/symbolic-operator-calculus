# Mapa de referencias para la corrección completa

## Contrato de evidencia

La consulta se hizo con `litctl` en el índice local, en modo de solo lectura, en
el commit `0275d2e8332c000c5677062b3471a16105751fa7`. Para cada afirmación
utilizada se exigió clave BibTeX resuelta, registro de teorema y lectura de las
páginas pertinentes. Una fuente temática o un resultado de búsqueda no se
convirtió en regla.

## Resultados verificados en su marco de origen

| Clave | Resultado | Evidencia | Hipótesis esenciales | Conclusión de la fuente | Aplicación a la corrección |
|---|---|---|---|---|---|
| `KKL2014Regularization` | Thm. 5.8 | impresa 206; PDF 18; checksum `7ed69d096f9e31983cb9b4701e8b486fd15736d82163144f27f9bae244ce41f6` | símbolo en `E-tilde(R+,V(R))`; para suficiencia, no anulación de todos los valores de borde y fibra de (5.46) | criterio de Fredholm para `Op(a)` y regularizadores `Op(b)+K` con valores recíprocos | `unproved_adaptation`: no se ha construido el símbolo de la corrección, probado su clase ni las condiciones de no anulación |
| `Karlovich2025Cusps` | Lemma 3.1 | artículo 102, pp. 8–10; PDF 8–10; checksum `301f8435688416e6ccd88be1a44f85e12372c9a557bdc1dbea301318b846be23` | curva/peso/índices de (3.7)–(3.9) | compactitud de tres errores transportados explícitos | `unproved_adaptation`: los bloques `W_{2,1}^{-}`, `W_{1,2}^{+}` no han sido identificados con esos errores |
| `Karlovich2025Cusps` | Lemma 3.2 | artículo 102, pp. 11–12; PDF 11–12; mismo checksum | hipótesis de la Sección 3 y casos (3.16)–(3.18) | compactitud de los errores conjugados explícitos | `unproved_adaptation`: falta la identificación bloque por bloque |
| `Karlovich2025Cusps` | Cor. 3.3 | artículo 102, pp. 12–13; PDF 12–13; mismo checksum | curva cuspídea (1.1)–(1.2), relación de pesos y acotación del operador de Cauchy | matrices transportadas de Cauchy módulo compactos | `unproved_adaptation`: no certifica automáticamente las sustituciones de Fase N |
| `Karlovich2025Cusps` | Thm. 6.1 | artículo 102, p. 24; derivación PDF 20–24; mismo checksum | los operadores son los generadores transportados de la álgebra definida en la fuente | símbolos matriciales explícitos de los generadores | `unproved_adaptation`: los factores de `C_2` no están identificados con esos generadores |
| `Karlovich2025Cusps` | Thm. 6.2 | artículo 102, p. 25; PDF 25; mismo checksum | curva, peso y `p` de la fuente; operador perteneciente a la álgebra `A` de (3.1) | homomorfismo de álgebra y criterio determinantal para operadores de `A` | `unproved_adaptation`: pertenencia de la palabra completa y compatibilidad de shifts/`R_1` abiertas |

Ninguna fila anterior es una aplicación directa ni una especialización ya
completada. `litctl use audit` mantiene abiertas las tres adaptaciones: el
regularizador Mellin, la localización compacta de los bloques y el criterio de
la álgebra cuspídea mixta.

## Fuentes consultadas sin resultado canónico aplicable

| Fuente | Estado del índice | Resultado de la consulta | Consecuencia aquí |
|---|---|---|---|
| `Karlovich2017HasemanSO`, Thms. 5.6, 6.2 y 7.1 | `ambiguous`, sin enlace ni texto | no existen registros de teorema; no se leyó el capítulo correcto | `ANALYTIC_PROOF_OBLIGATION: source verification pending` |
| artículo 2014, *Fredholmness and index of simplest singular integral operators with two slowly oscillating shifts*, DOI `10.7153/oam-08-52` | PDF local leído en el índice, pero sin clave BibTeX canónica | Thm. 5.2, pp. impresas 952–954/PDF 18–20, y prueba de Thm. 1.1 en PDF 20 respaldan los auxiliares en su propio marco | no se crea una regla certificada hasta resolver la omisión bibliográfica y comprobar las hipótesis branchwise |
| `Karlovich2008Nonlocal` | resolución manual, PDF y texto disponibles; cero registros verificados | búsquedas localizan resultados Mellin y shifts, pero no hay registro de teorema auditado | fuente prometedora; `source verification pending` |
| `KarlovichMonsivais2024MellinPDO` | resolución manual por título/DOI, PDF y texto disponibles; cero registros verificados | búsquedas localizan semiproductos y álgebras Mellin, pero no hay registro de teorema auditado | no justifica todavía ninguna interfaz `R_1 × ...` ni cierre |
| `Karlovich2025Cusps`, Thm. 7.1 | registro con enunciado, prueba fuente no disponible bajo clave canónica | estado `unverified` | no se usa la reducción real-lineal/compleja |

## Correspondencia fuente–interfaz

| Interfaz u obligación | Evidencia disponible | Estado operativo |
|---|---|---|
| sustituciones Wiener--Hopf de Fase N | Cor. 3.3 y sus lemas son análogos verificables en el modelo 2025, pero falta identificar los bloques concretos | regla suministrada `MOD_COMPACT_EQUIVALENCE`, `UNCERTIFIED` |
| `Cauchy factor × R_1` | KKL 5.8 no es una regla de semiproducto | `SOURCE_UNVERIFIED` |
| `R_1 × U_1` y `R_1 × \widehat G_1` | 2008/2024 sugieren tecnología Mellin; ningún teorema aplicable está indexado | `SOURCE_UNVERIFIED` |
| cierre en una clase Mellin | KKL 5.8 exige membresía antes de aplicar el criterio | obligación abierta |
| cierre en la álgebra cuspídea | Karlovich 2025, Thm. 6.2, solo aplica a elementos de `A` | obligación abierta |
| Fredholmness | ambos caminos exigen hipótesis posteriores no comprobadas | fuera del alcance de esta fase |

## Conclusión bibliográfica

La literatura verificada ofrece dos marcos de destino plausibles, pero no una
regla central que absorba la palabra no conmutativa completa. No se añade al
motor ninguna regla analítica nueva. La clasificación correcta de toda
aplicación externa en esta fase es `unproved_adaptation` o
`source verification pending`.
