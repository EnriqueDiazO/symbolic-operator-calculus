# Revisión adversarial de la cancelación de dilataciones en la corrección de Schur

Fecha de revisión: 2026-07-21.

## Dictamen

```text
BLOCKED_WITH_PRECISE_OBLIGATION
```

La identidad abstracta

\[
U_\gamma\operatorname{Op}_{M,\delta}(r)U_{\gamma^{-1}}
=\operatorname{Op}_{M,\delta}(r^\gamma),
\qquad
r^\gamma(x,\lambda)=r(\gamma x,\lambda),
\]

es exacta con las convenciones del artículo. También es correcto propagar
representaciones de los tres factores módulo compactos a su producto, siempre
que todos los operadores sean acotados sobre la misma realización ponderada.

Lo que no está demostrado es que los factores reales que rodean a
\(\mathcal R_1\) en la corrección del artículo tengan las formas
\(W(b_{2,1})U_\gamma\) y
\(U_{\gamma^{-1}}W(b_{1,2})\). Los bloques contienen multiplicadores de
transporte, coeficientes, localizadores, factores de Cauchy y los operadores
auxiliares \(\mathcal T_{1,-}\) y \(Z_1^{-1}\). No existe una regla certificada
que permita atravesarlos o eliminarlos para hacer adyacente el par de
dilataciones. Hay además una diferencia material entre el regularizador
Mellin \(\mathcal R_1\) y el regularizador completo
\(\mathcal A_{1,1}^{(-1)}\), así como un conflicto de signo en la notación de
la corrección.

En consecuencia, esta fase no certifica la cancelación en el pivote real, la
pertenencia del sándwich final a un álgebra mixta ni una conclusión de
Fredholm.

## Alcance y criterio de los veredictos

Se auditan por separado dos niveles:

- **Identidad abstracta**: el sándwich ya viene escrito, sobre un mismo
  espacio, como
  \(W(b_{2,1})U_\gamma\operatorname{Op}_{M,\delta}(r)
  U_{\gamma^{-1}}W(b_{1,2})\).
- **Aplicación real**: se parte de los bloques y del regularizador que aparecen
  efectivamente en `sections/06_shur_reduction.tex`.

Los estados significan:

| Estado | Criterio |
|---|---|
| `PASS` | La afirmación está demostrada con las convenciones y objetos exactos relevantes. |
| `PASS_WITH_ASSUMPTION` | El paso es correcto sólo bajo una hipótesis indicada y conservada en la salida. |
| `OPEN` | No hay contradicción conocida, pero falta una demostración aplicable al objeto real. |
| `FAIL` | La afirmación o identificación propuesta no coincide con las fórmulas reales o usa un paso inválido. |

La presencia de `FAIL` impide declarar la fase completada.

## Convenciones fijadas antes de auditar

Sea

\[
X_\delta^p:=L^p(\mathbb R_+,w_\delta),
\qquad
\|f\|_{X_\delta^p}^p
=\int_0^\infty |x^\delta f(x)|^p\,dx,
\]

y póngase

\[
\kappa_\delta:=\delta+\frac1p,
\qquad
(V_\gamma f)(x):=f(\gamma x),
\qquad
U_\gamma:=\gamma^{\kappa_\delta}V_\gamma.
\]

Entonces

\[
\|V_\gamma f\|_{X_\delta^p}
=\gamma^{-\kappa_\delta}\|f\|_{X_\delta^p},
\qquad
\|U_\gamma f\|_{X_\delta^p}=\|f\|_{X_\delta^p},
\qquad
U_\gamma^{-1}=U_{\gamma^{-1}}.
\]

La isometría

\[
(\Phi_\delta f)(x)=x^{\kappa_\delta}f(x)
\]

lleva \(X_\delta^p\) a
\(L^p(\mathbb R_+,dx/x)\). Para evitar un error de tipos, en este informe
se define la realización ponderada del Mellin PDO por

\[
\operatorname{Op}_{M,\delta}(r)
:=\Phi_\delta^{-1}\operatorname{Op}_M(r)\Phi_\delta.
\]

La cuantización de Mellin es la izquierda:

\[
[\operatorname{Op}_M(r)f](x)
=\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
r(x,\lambda)\left(\frac{x}{y}\right)^{i\lambda}
f(y)\frac{dy}{y}\,d\lambda.
\]

Con este orden,

\[
U_\gamma\operatorname{Op}_{M,\delta}(r)U_{\gamma^{-1}}
=\operatorname{Op}_{M,\delta}\bigl(r(\gamma\,\cdot,\lambda)\bigr).
\]

El orden inverso produce \(r(x/\gamma,\lambda)\). Estas dos fórmulas no
son intercambiables.

Para \(m=2\), el artículo no introduce literalmente
\(\gamma_{k,j}\). Introduce

\[
\beta_{k,j}=\alpha_j^{-1}\circ\alpha_k,
\qquad
\beta_{k,j}(x)=\frac{\gamma_k}{\gamma_j}x.
\]

En lo que sigue,

\[
\gamma:=\frac{\gamma_2}{\gamma_1},
\qquad
\beta_{2,1}(x)=\gamma x,
\qquad
\beta_{1,2}(x)=\gamma^{-1}x.
\]

## Lema abstracto que sí supera la revisión

Supóngase que todos los operadores siguientes actúan acotadamente en el
mismo \(X_\delta^p\), que la normalización de ambas dilataciones coincide y
que

\[
\begin{aligned}
A_{2,1}&\equiv W(b_{2,1})U_\gamma,\\
R_{1,1}&\equiv \operatorname{Op}_{M,\delta}(r_{1,1}),\\
A_{1,2}&\equiv U_{\gamma^{-1}}W(b_{1,2})
\end{aligned}
\qquad \pmod{\mathcal K(X_\delta^p)}.
\]

Supóngase además que
\(r_{1,1}\) y
\(r_{1,1}^\gamma(x,\lambda):=r_{1,1}(\gamma x,\lambda)\)
pertenecen a una clase que da operadores Mellin acotados. Entonces

\[
\begin{aligned}
A_{2,1}R_{1,1}A_{1,2}
&\equiv
W(b_{2,1})U_\gamma
\operatorname{Op}_{M,\delta}(r_{1,1})
U_{\gamma^{-1}}W(b_{1,2})\\
&=
W(b_{2,1})
\operatorname{Op}_{M,\delta}(r_{1,1}^\gamma)
W(b_{1,2})
\pmod{\mathcal K(X_\delta^p)}.
\end{aligned}
\]

La segunda línea usa una `EXACT_IDENTITY`. La primera conserva la relación
más débil de las tres representaciones de entrada. Si las tres son exactas,
todo el resultado es exacto. Si la representación del regularizador es sólo
una hipótesis, el resultado completo permanece
`CONDITIONAL_ON_REGULARIZER_REPRESENTATION`.

Esta proposición no demuestra sus propias hipótesis de factorización para los
bloques del artículo.

## Matriz de los diez controles

| Nº | Control | Identidad abstracta | Aplicación real | Resultado adversarial |
|---:|---|---|---|---|
| 1 | Escala radial | `PASS` | `OPEN` | El orden \(U_\gamma\operatorname{Op}U_{\gamma^{-1}}\) da \(r(\gamma x,\lambda)\); falta exponer ese orden en los factores reales. |
| 2 | Factor ponderado | `PASS` | `PASS_WITH_ASSUMPTION` | \(\kappa_\delta=\delta+1/p\). Los escalares se cancelan sólo si el par usa el mismo \(p,\delta\) y la misma normalización; \(\rho_k\) no es parte de esa normalización. |
| 3 | Signo de Fourier | `PASS` | `PASS_WITH_ASSUMPTION` | Con \(e^{-it\lambda}\), \(V_aW(b)V_a^{-1}=W(b(\cdot/a))\). La regla exacta se verificó para el núcleo auxiliar, no para todo el bloque de Schur. |
| 4 | Truncación Wiener--Hopf | `PASS` | `PASS` | Una dilatación positiva conserva \((0,\infty)\), por lo que extensión por cero y restricción a la semirrecta no cambian el signo ni la escala. Esto no trata \(\chi_\infty\). |
| 5 | Orden de composición | `PASS` | `FAIL` | En el lema abstracto el par inverso es adyacente. En \(\mathcal B_{2,1}^{-}\mathcal R_1\mathcal B_{1,2}^{+}\) intervienen \(\mathcal T_{1,-}\), \(Z_1^{-1}\), \(\widetilde V_k-G_k\), proyecciones y cortes; no se ha probado la factorización requerida ni es lícito conmutarlos. |
| 6 | Estabilidad módulo compactos al multiplicar | `PASS` | `PASS` | \(\mathcal K(X)\) es un ideal bilateral y los factores usados son acotados. Esto propaga equivalencias, pero no autoriza permutar factores. |
| 7 | Estabilidad de símbolos | `PASS` | `OPEN` | \(\widetilde{\mathcal E}\) es estable bajo \(x\mapsto\gamma x\), pero no está probada la pertenencia del sándwich WH--Mellin--WH ni de la ampliación con dilataciones. |
| 8 | Efecto de \(\chi_\infty\) | `PASS_WITH_ASSUMPTION` | `OPEN` | La identidad desnuda no contiene cortes. En el bloque real aparecen \(\chi_\infty(\gamma_kx)\) y \(\chi_\infty(\gamma_jy)\). Hace falta conservarlos o consumir explícitamente una equivalencia compacta aplicable a la palabra completa. |
| 9 | Espacio del regularizador | `PASS_WITH_ASSUMPTION` | `FAIL` | \(\mathcal R_1=\Phi_\delta^{-1}\operatorname{Op}(r_{1,1})\Phi_\delta\) es Mellin; el regularizador de \(\mathcal A_{1,1}\) es \(\mathcal T_{1,-}\mathcal R_1Z_1^{-1}\). Identificarlos borra factores no triviales. |
| 10 | Caso matricial | `PASS_WITH_ASSUMPTION` | `OPEN` | Para matrices finitas, los residuos compactos se suman entrada por entrada. No hay identificación Mellin de las entradas cuspídeas ni cálculo Fredholm de una álgebra matricial mixta aplicable al pivote. |

## Guardas exigidas a la implementación

La revisión del código encontró y cerró varias rutas de sobrecertificación. La
API final falla cerradamente salvo que:

- el operador y las dos dilataciones estén declarados acotados;
- las escalas sean inversas y coincidan espacio, (p), peso, normalización y
  convención de acción;
- los contextos de hipótesis no contengan una contradicción explícita;
- Wiener--Hopf use la convención de Fourier del repositorio y Mellin PDO use
  su cuantización izquierda;
- la estabilidad por escala radial tenga evidencia explícita, conservada en
  el símbolo transformado;
- no haya un localizador espacial oculto: el tipo actual no puede representar
  exactamente (chi(x)\mapsto\chi(\gamma x)), por lo que lo rechaza incluso si
  un booleano lo marca como controlado;
- todas las equivalencias módulo compactos nombren un mismo ideal.

Un resultado bloqueado conserva la palabra de Schur original, las relaciones
y evidencia de entrada y sus capas lógicas, pero no crea una relación
semántica de salida. Si la representación Mellin del regularizador sólo se
asume, la conclusión completa se conserva como `FormalIdentity`; nunca se
promueve a `ExactIdentity` ni a `ModCompactEquivalence` no condicionada.

## Desarrollo de los controles

### 1. Escala radial

Desde la definición integral, para una función de prueba \(f\),

\[
\begin{aligned}
[D_\gamma\operatorname{Op}_M(r)D_{\gamma^{-1}}f](x)
&=\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
r(\gamma x,\lambda)
\left(\frac{\gamma x}{y}\right)^{i\lambda}
f(y/\gamma)\frac{dy}{y}\,d\lambda\\
&=\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
r(\gamma x,\lambda)
\left(\frac{x}{z}\right)^{i\lambda}
f(z)\frac{dz}{z}\,d\lambda,
\end{aligned}
\]

con \(y=\gamma z\). Por tanto, la escala correcta para la cadena escrita en
el encargo es \(\gamma x\). El documento
`docs/LEFT_DILATION_MELLIN_COVARIANCE_DERIVATION.md:13-36,54-98`
demuestra la versión de orden inverso, con \(x/\gamma\); ambas son coherentes.

El fallo de aplicación no está en esta cuenta, sino en que el artículo aún
no produce el par adyacente alrededor de \(\mathcal R_1\).

### 2. Factor ponderado

El artículo usa \(\int|x^\delta f(x)|^pdx\), no
\(\int|f(x)|^px^\delta dx\). El cambio \(y=\gamma x\) da

\[
\|\gamma^\nu f(\gamma\,\cdot)\|_{X_\delta^p}^p
=\gamma^{p\nu-p\delta-1}\|f\|_{X_\delta^p}^p,
\]

y fija \(\nu=\kappa_\delta\). Aunque la covariancia por conjugación también
vale para la dilatación no normalizada, no se pueden mezclar \(V_\gamma\) y
\(U_\gamma\) sin registrar los escalares. Tampoco puede absorberse en ese
escalar el multiplicador no constante

\[
\rho_k(x)=\frac{\gamma_kx-ic_k}{x-ic_k}
\]

que forma parte de \(\widetilde V_k=\rho_kV_k\).

### 3. Signo de Fourier

La convención del artículo es

\[
(\mathcal Ff)(\lambda)=\int_{\mathbb R}f(t)e^{-it\lambda}\,dt.
\]

Para \(h_a(t)=a h(at)\), \(a>0\),

\[
\mathcal Fh_a(\lambda)=(\mathcal Fh)(\lambda/a).
\]

De aquí se obtiene exactamente

\[
V_aW(b)V_a^{-1}=W\bigl(b(\cdot/a)\bigr),
\]

sin signo adicional ni factor \(2\pi\). Las dos colocaciones correctas son

\[
V_aW(b)=W\bigl(b(\cdot/a)\bigr)V_a,
\qquad
W(b)V_a=V_aW\bigl(b(a\,\cdot)\bigr).
\]

El artículo verifica estas escalas en su operador auxiliar
\(\mathcal W_{k,j}=V_kW(b_{k,j})V_j^{-1}\). Esta verificación no identifica
por sí sola \(A_{k,j}\) con \(\mathcal W_{k,j}\).

### 4. Truncación Wiener--Hopf

El operador Wiener--Hopf se define por extensión desde la semirrecta,
convolución en la recta y restricción posterior. Como \(a>0\), la aplicación
\(x\mapsto ax\) preserva exactamente \(\mathbb R_+\). El cambio de variable no
produce un término de borde en cero ni cambia \(\chi_+\) por \(\chi_-\).

Este `PASS` se refiere a la truncación estructural \(\chi_+\) de Wiener--Hopf.
El localizador cuspídeo \(\chi_\infty\) es otro operador y se audita en el
punto 8.

### 5. Orden de composición

La identidad auxiliar exacta del artículo es

\[
\mathcal W_{k,j}
=V_{\beta_{k,j}}W(b_{k,j}^{(j)})
=W(b_{k,j}^{(k)})V_{\beta_{k,j}}.
\]

Para \(\gamma=\gamma_2/\gamma_1\), pueden escogerse las colocaciones

\[
\mathcal W_{2,1}=W(b_{2,1}^{(2)})V_\gamma,
\qquad
\mathcal W_{1,2}=V_{\gamma^{-1}}W(b_{1,2}^{(2)}).
\]

Eso sí produciría una cancelación al rodear inmediatamente un Mellin PDO.
Pero los bloques reales son

\[
A_{k,j}=\frac12(\widetilde V_k-G_k)T_{k,j}
\]

y, en el modelo de Schur,

\[
\begin{aligned}
\mathcal B_{2,1}^{-}
&=(\widetilde V_2-G_2)W_{2,1}^{-}\mathcal T_{1,-},\\
\mathcal B_{1,2}^{+}
&=Z_1^{-1}(\widetilde V_1-G_1)W_{1,2}^{+}.
\end{aligned}
\]

No es posible sustituir estas palabras completas por
\(W(b_{2,1})V_\gamma\) y
\(V_{\gamma^{-1}}W(b_{1,2})\) usando sólo la identidad del núcleo auxiliar.
El orden real constituye un `FAIL` de la aplicación propuesta.

También debe conservarse el signo. El pivote se define por

\[
A_{2,2}^{(1)}=A_{2,2}-A_{2,1}A_{1,1}^{(-1)}A_{1,2}.
\]

Como el modelo normalizado de \(A_{2,1}\) lleva signo negativo, el artículo
denomina \(C_{2,2}^{(1)}\) al término positivo que se suma al pivote. Por
tanto, módulo compactos,

\[
C_{2,2,\mathrm{art}}^{(1)}
\equiv-\,A_{2,1}A_{1,1}^{(-1)}A_{1,2},
\]

no al producto crudo con signo positivo usado en el enunciado inicial del
encargo.

### 6. Estabilidad módulo compactos al multiplicar

Si \(A-A'\), \(R-R'\) y \(B-B'\) son compactos y los seis operadores son
acotados en el mismo espacio, entonces

\[
ARB-A'R'B'
=(A-A')RB+A'(R-R')B+A'R'(B-B')
\]

es compacto. La misma afirmación vale para matrices finitas, entrada por
entrada. Este paso no necesita conmutatividad.

Su alcance es limitado: no demuestra ninguna de las tres equivalencias de
entrada, no cambia el orden y no convierte un operador no identificado en un
Mellin PDO.

### 7. Estabilidad de símbolos

Para \(r\in\widetilde{\mathcal E}(\mathbb R_+,V(\mathbb R))\), el cambio
\(x\mapsto x/\gamma\) conserva la norma en \(C_b(\mathbb R_+,V)\), la
oscilación radial lenta, la continuidad uniforme por traslaciones en
\(\lambda\), las colas uniformes de \(\partial_\lambda r\) y las fibras de
extremo. Esto está demostrado en
`docs/RADIAL_SCALING_INVARIANCE_OF_MELLIN_SYMBOLS.md:3-105`.

La estabilidad para \(r(\gamma x,\lambda)\) sigue aplicando el mismo resultado
con la constante \(\gamma^{-1}\). Por tanto, la escala del símbolo Mellin
supera la revisión.

No se sigue que

\[
W(b_{2,1})\operatorname{Op}_{M,\delta}(r^\gamma)W(b_{1,2})
\]

pertenezca a \(\widetilde{\mathcal E}\), a la álgebra de BKM, a la de
Duduchava o a otra álgebra con cálculo Fredholm. La estabilidad radial y el
cierre mixto son obligaciones distintas.

### 8. Efecto de \(\chi_\infty\)

El artículo fija \(\chi_\infty=0\) en \([0,1]\) y
\(\chi_\infty=1\) en \([2,\infty)\), y usa

\[
W_{k,j}^{\infty}
=M_{\chi_\infty}W(b_{k,j})M_{\chi_\infty}.
\]

Bajo conjugación,

\[
V_aM_{\chi_\infty}V_a^{-1}
=M_{\chi_\infty(a\,\cdot)}.
\]

En consecuencia, el núcleo localizado conjugado contiene un corte de salida
\(\chi_\infty(\gamma_kx)\) y uno de entrada
\(\chi_\infty(\gamma_jy)\). No desaparecen por la covariancia de \(W(b)\).

El artículo sí dispone de una equivalencia compacta entre el bloque
off-diagonal y su localización. Esa equivalencia puede propagarse por factores
acotados. Sin embargo, una derivación de cancelación debe elegir y documentar
una de estas dos rutas:

1. conservar ambos cortes con sus escalas exactas durante toda la
   factorización; o
2. reemplazar el bloque localizado por el no localizado mediante la
   equivalencia compacta concreta y registrar el residuo resultante.

La identidad abstracta no realiza ninguna de las dos rutas para
\(\mathcal B_{2,1}^{-}\) y \(\mathcal B_{1,2}^{+}\); por eso la aplicación
permanece `OPEN`. Eliminar los cortes declarando una igualdad exacta sería
`FAIL`.

### 9. Espacio del regularizador

El artículo construye primero

\[
\mathcal R_1
=\Phi_\delta^{-1}\operatorname{Op}(r_{1,1})\Phi_\delta,
\]

que regulariza al producto normalizado \(\mathcal N_1\). El regularizador que
se inserta para el bloque diagonal es, en cambio,

\[
\mathcal A_{1,1}^{(-1)}
=\mathcal T_{1,-}\mathcal R_1Z_1^{-1}.
\]

Por tanto:

- si `R11` significa \(\mathcal R_1\), su representación como Mellin PDO
  transportado es exacta una vez elegido el símbolo regularizador;
- si `R11` significa el regularizador de \(\mathcal A_{1,1}\) que aparece en
  el complemento de Schur, la afirmación
  \(R_{1,1}\simeq\operatorname{Op}_M(r_{1,1})\) omite
  \(\mathcal T_{1,-}\) y \(Z_1^{-1}\) y es falsa como identificación de los
  objetos escritos.

El artículo absorbe correctamente esos dos factores en
\(\mathcal B_{2,1}^{-}\) y \(\mathcal B_{1,2}^{+}\), pero esas palabras
normalizadas todavía no han sido factorizadas con un par de dilataciones
adyacentes. Éste es un fallo de aplicación, no de la covariancia Mellin.

### 10. Caso matricial

En una matriz finita de operadores, una matriz cuyos bloques son compactos es
compacta en la suma finita de los espacios de rama. Por ello, la propagación
de residuos compactos y la eliminación de un par inverso ya adyacente se
extienden entrada por entrada.

Eso no resuelve la situación real. Para usar un cálculo matricial de Fredholm
se requieren primero símbolos Mellin admisibles para las entradas
off-diagonales localizadas. El registro de KM 2024, Teorema 7.1, es sólo una
analogía para estrellas sin cúspide: la hipótesis de separación angular falla
en la pure cusp star. No puede suministrar el símbolo matricial faltante.

## Trazabilidad al artículo auxiliar

Raíz de solo lectura:

```text
/home/enriquedo/PersonalProjects/Papers/DiazOcampoHasemanBVP2026
```

| Hecho usado | Ruta y líneas | Fuerza |
|---|---|---|
| Peso \(w_\delta=x^\delta\) | `sections/03_transport_and_kernels.tex:241-262` | Exacta |
| Norma real, \(\Phi_\delta\) y \(\kappa_\delta\) | `sections/03_transport_and_kernels.tex:1064-1108` | Exacta |
| \(V_\gamma\), \(U_\gamma\), isometría y multiplicador \(\gamma^{i\lambda}\) | `sections/03_transport_and_kernels.tex:1368-1534` | Exacta |
| Cuantización Mellin izquierda | `sections/03_transport_and_kernels.tex:870-969` | Definición exacta |
| Localizador \(\chi_\infty\) y residuo compacto | `sections/03_transport_and_kernels.tex:734-773` | Módulo compactos |
| Fourier \(e^{-it\lambda}\) y definición de \(W(b)\) | `sections/04_target_theorem53.tex:17-48` | Definición exacta |
| Símbolos off-diagonales y sus signos | `sections/04_target_theorem53.tex:84-179` | Exacta en el modelo WH |
| Localización WH de \(T_{k,j}\) | `sections/04_target_theorem53.tex:228-237,288-338` | Módulo compactos |
| Convención \(V_kf=f(\gamma_k\cdot)\) y \(\widetilde V_k=\rho_kV_k\) | `sections/05_haseman_operator_with_shift.tex:80-160` | Exacta |
| Bloques reales \(A_{k,j}\) | `sections/05_haseman_operator_with_shift.tex:301-320` | Exacta |
| Modelo localizado de los bloques | `sections/05_haseman_operator_with_shift.tex:1668-1726` | Módulo compactos |
| Factorización relativa de \(\mathcal W_{k,j}\) | `sections/05_haseman_operator_with_shift.tex:1790-2166` | Exacta, sólo para el auxiliar no localizado |
| Cortes reescalados tras conjugación | `sections/05_haseman_operator_with_shift.tex:2172-2278` | Exacta |
| Símbolo Mellin candidato del bloque localizado | `sections/05_haseman_operator_with_shift.tex:2282-2472` | Verificación de clase pendiente |
| Signos de \(A_{1,2}\) y \(A_{2,1}\) para \(m=2\) | `sections/06_shur_reduction.tex:44-95` | Módulo compactos |
| \(\mathcal R_1\) y \(\mathcal A_{1,1}^{(-1)}\) | `sections/06_shur_reduction.tex:109-195` | Representaciones exactas; propiedad regularizadora módulo compactos |
| Definición y signo del pivote | `sections/06_shur_reduction.tex:281-348` | Modelo módulo compactos |
| \(C_{2,2}^{(1)}=\mathcal B_{2,1}^{-}\mathcal R_1\mathcal B_{1,2}^{+}\) | `sections/06_shur_reduction.tex:359-432` | Exacta para la corrección modelo |

## Puerta de evidencia bibliográfica

Raíz de solo lectura:

```text
/home/enriquedo/PersonalProjects/haseman-literature-index
```

La auditoría bibliográfica concluye `Caso D`: las identidades de covariancia
son derivaciones directas exactas, pero ninguna fuente verificada proporciona
una misma álgebra que contenga los símbolos WH cuspídeos concretos, Mellin PDO
radiales, dilataciones relativas y un cálculo Fredholm inverse-closed. Véanse
`reports/MIXED_WIENER_HOPF_MELLIN_DILATION_AUDIT.md:5-38,286-379,403-465`.

### Resultados que sí pueden usarse, con sus límites

| Fuente verificada | Registro | Lo que aporta y lo que no aporta |
|---|---|---|
| BKM 2026, Teoremas 8.3 y 8.6 | `catalog/theorem_index.yaml:3-89` | Inverse-closedness y criterio local sólo para miembros de su \(\mathfrak D\). Los símbolos WH actuales fallan las uniones de esa clase y no hay generador de dilatación. |
| Duduchava 1987, Lemas 2.4, 2.6 y Teorema 3.2 | `catalog/theorem_index.yaml:90-222` | Conmutador Mellin--WH compacto y cálculo de su álgebra exacta; no cubre un factor de dilatación no trivial ni un Mellin PDO radial general. |
| Karlovich 2006, Teoremas 11.11, 12.1 y 12.3 | `catalog/theorem_index.yaml:315-450` | Cálculo de Mellin PDO, compactos e inversas dentro de la clase; no es una interfaz WH--dilatación. |
| Karlovich 2009, Teoremas 4.1 y 4.3 | `catalog/theorem_index.yaml:451-554` | Compactitud y Fredholm para Mellin PDO admisibles; exige pertenencia previa del símbolo. |
| Karlovich 2025, resultados 3.1--3.3, 4.2--4.4, 5.1, 6.1--6.2 | `catalog/theorem_index.yaml:855-1243` | Modelo cuspídeo y álgebras ambiente dentro de su alcance; no incorpora Mellin PDO radiales generales ni dilataciones constantes en la álgebra de semirrecta. |
| Karlovich--Monsiváis 2024, Teoremas 4.3, 5.1 y 7.1 | `catalog/theorem_index.yaml:1321-1481` | Semiproducto y cálculo Mellin; 7.1 sólo sirve como analogía sin cúspide, no como especialización a la pure cusp star. |
| Karlovich--Ramírez de Arellano 2004, Teoremas 5.7 y 6.4 | `catalog/theorem_index.yaml:1482-1576` | Pertenencia de núcleos exactos con singularidad fija; falta identificarlos, por igualdad de núcleos o equivalencia compacta, con \(W(b)V_\gamma\). Su \(M_\gamma\) es multiplicación por una potencia, no una dilatación. |
| Karlovich--Rosales-Méndez 2023, resultados verificados 3.1--3.3, 4.1--4.3, 9.2 y 10.3 | `catalog/theorem_index.yaml:1577-1921,2242-2293` | Antecedentes de shifts y álgebras en estrellas de rayos; no cierran el producto cuspídeo WH--Mellin--dilatación. |
| KKL 2014, Lema 5.7 y Teorema 5.8 | `catalog/theorem_index.yaml:2294-2373` | Si el símbolo está en \(\widetilde E\) y satisface todas las condiciones de no anulación, el regularizador es \(\operatorname{Op}(b)+K\). No prueba esas premisas para el bloque completo. |

Los usos aplicables conservan abiertas precisamente las interfaces necesarias:
`catalog/uses.yaml:394-462` para las entradas matriciales Mellin y
`catalog/uses.yaml:492-570` para BKM, Duduchava y la extensión por
dilataciones. La ruta alternativa de singularidades fijas también permanece
abierta en `catalog/uses.yaml:588-615`.

### Resultados excluidos de esta revisión

No se usan como teoremas certificados:

| Fuente o registros | Líneas | Motivo de exclusión |
|---|---|---|
| Karlovich 2025, resultados 7.1 y 7.2 | `catalog/theorem_index.yaml:1244-1320` | `proof_read: false` y `status: unverified`. |
| Karlovich--Rosales-Méndez 2023, resultados 5.1--5.4 y 6.1--6.3 | `catalog/theorem_index.yaml:1922-2241` | `proof_read: false` y `status: unverified`; sólo se mencionan, cuando procede, como analogía. |
| Karlovich--Rosales 2022 QC | `catalog/papers.yaml:1019-1037` | `verification_status: unverified` y sin registro theorem-centric certificado. |

Las fórmulas de dilatación de
`reports/MIXED_WIENER_HOPF_MELLIN_DILATION_AUDIT.md:286-348`
son derivaciones directas del informe, no teoremas externos; el propio informe
lo declara en `:484-486`.

## Obligación precisa que bloquea la cancelación real

La interfaz correcta no es el producto esquemático con un regularizador
desnudo, sino la identidad de modelo

\[
C_{2,2,\mathrm{art}}^{(1)}
=\mathcal B_{2,1}^{-}\mathcal R_1\mathcal B_{1,2}^{+},
\]

con

\[
\begin{aligned}
\mathcal B_{2,1}^{-}
&=(\widetilde V_2-G_2)W_{2,1}^{-}\mathcal T_{1,-},\\
\mathcal B_{1,2}^{+}
&=Z_1^{-1}(\widetilde V_1-G_1)W_{1,2}^{+},\\
\mathcal R_1
&=\Phi_\delta^{-1}\operatorname{Op}(r_{1,1})\Phi_\delta.
\end{aligned}
\]

Para desbloquear la fase debe probarse, sobre el mismo \(X_\delta^p\), una
descomposición exacta o módulo compactos de **esas dos palabras completas**.
La prueba debe:

1. expandir o factorizar \(\mathcal B_{2,1}^{-}\) y
   \(\mathcal B_{1,2}^{+}\) sin cambiar el orden de
   \(\mathcal T_{1,-}\), \(Z_1^{-1}\), las proyecciones, \(G_k\),
   \(\rho_k\) y los bloques WH;
2. identificar, para cada sumando no compacto, un par
   \(U_g\,\mathcal R_1\,U_{g^{-1}}\) realmente adyacente, o demostrar que el
   sumando restante es compacto por un resultado aplicable;
3. conservar los signos de \(A_{2,1}\) y de
   \(C_{2,2,\mathrm{art}}^{(1)}\);
4. conservar \(\chi_\infty(\gamma_kx)\) y
   \(\chi_\infty(\gamma_jy)\), o eliminar cada par de cortes usando una
   equivalencia compacta explícita y registrar su residuo;
5. verificar que el símbolo resultante permanece en la clase Mellin
   declarada;
6. probar por separado que la suma final pertenece a una misma álgebra cuyo
   cociente admita regularizadores y un símbolo de Fredholm.

Una factorización suficiente, pero actualmente no demostrada, sería

\[
\mathcal B_{2,1}^{-}\equiv L_{2,1}U_\gamma,
\qquad
\mathcal B_{1,2}^{+}\equiv U_{\gamma^{-1}}L_{1,2}
\pmod{\mathcal K(X_\delta^p)},
\]

con todos los factores omitidos absorbidos en \(L_{2,1}\) y \(L_{1,2}\) por
identidades o semiproductos certificados, no por notación. También sería
admisible una descomposición finita por sumandos con distintas razones de
escala, siempre que cada par y cada residuo se controle por separado. No basta
con demostrar la identidad para
\(V_2W(b)V_1^{-1}\), pues ése es sólo el núcleo auxiliar.

## Clasificación final

### Demostrado exactamente

- normalización, inversa y composición de \(U_\gamma\);
- escala Fourier
  \(V_\gamma W(b)V_{\gamma^{-1}}=W(b(\cdot/\gamma))\);
- covariancia Mellin
  \(U_\gamma\operatorname{Op}_{M,\delta}(r)U_{\gamma^{-1}}
  =\operatorname{Op}_{M,\delta}(r(\gamma\cdot,\lambda))\);
- factorización relativa del auxiliar no localizado
  \(\mathcal W_{k,j}\);
- representación central
  \(\mathcal R_1=\Phi_\delta^{-1}\operatorname{Op}(r_{1,1})\Phi_\delta\);
- identidad del modelo
  \(C_{2,2,\mathrm{art}}^{(1)}
  =\mathcal B_{2,1}^{-}\mathcal R_1\mathcal B_{1,2}^{+}\).

### Demostrado módulo compactos

- sustitución de \(T_{k,j}\) por el modelo WH localizado;
- forma normalizada de \(A_{1,2}\) y \(A_{2,1}\) en el caso \(m=2\);
- propagación de equivalencias ya demostradas por multiplicación con
  operadores acotados.

### Condicional

- cancelación del par de dilataciones en un sándwich que ya tenga las
  factorizaciones externas requeridas;
- identificación de un regularizador suministrado a la API con
  \(\operatorname{Op}_{M,\delta}(r)\), salvo cuando el objeto sea
  explícitamente \(\mathcal R_1\);
- eliminación de los cortes cuspídeos dentro de la palabra completa.

### Todavía abierto

- las factorizaciones de \(\mathcal B_{2,1}^{-}\) y
  \(\mathcal B_{1,2}^{+}\) que harían adyacente el par inverso;
- la pertenencia del sándwich WH--Mellin--WH a una álgebra mixta;
- la pertenencia del producto real y del pivote a un álgebra de Fredholm
  inverse-closed;
- el símbolo de la corrección y las condiciones de no anulación del pivote;
- cualquier conclusión de Fredholm o de índice para el pivote completo.

Por estas razones, el resultado de fase permanece
`BLOCKED_WITH_PRECISE_OBLIGATION`.
