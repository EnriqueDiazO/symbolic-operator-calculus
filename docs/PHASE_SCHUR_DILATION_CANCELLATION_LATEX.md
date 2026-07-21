conditional draft — not yet ready for insertion into the article

# Lema condicional de cancelación de dilataciones en la corrección de Schur

El siguiente fragmento es autocontenido, pero deliberadamente condicional. Su
conclusión sólo se activa después de demostrar las factorizaciones externas de
los bloques que se enumeran entre las hipótesis. No afirma que esas
factorizaciones valgan para los bloques reales del artículo.

## Fragmento LaTeX

```latex
\begin{lemma}[Cancelación condicional de dilataciones relativas]
\label{lem:conditional-relative-dilation-cancellation}
Sea $1<p<\infty$ y sea
\[
  -\frac1p<\delta<1-\frac1p,
  \qquad
  \nu:=\delta+\frac1p.
\]
Considérese el espacio
\[
  X_\delta:=L^p(\mathbb R_+,w_\delta),
  \qquad
  w_\delta(x):=x^\delta,
\]
con la convención
\[
  \|f\|_{X_\delta}^p
  :=
  \int_0^\infty |x^\delta f(x)|^p\,dx.
\]
Para $\gamma>0$, defínase la dilatación ponderada normalizada por
\[
  \bigl(V_\gamma^{(\nu)}f\bigr)(x)
  :=
  \gamma^\nu f(\gamma x),
  \qquad x>0.
\]
Sea
\[
  d\mu(x):=\frac{dx}{x},
  \qquad
  (\Phi_\delta f)(x):=x^\nu f(x),
\]
de modo que
\[
  \Phi_\delta:X_\delta\longrightarrow L^p(\mathbb R_+,d\mu)
\]
es un isomorfismo isométrico. Para un símbolo $r=r(x,\lambda)$, defínase
formalmente
\[
  \bigl(\operatorname{Op}_M(r)g\bigr)(x)
  :=
  \frac1{2\pi}
  \int_{\mathbb R}\int_0^\infty
  r(x,\lambda)
  \left(\frac{x}{y}\right)^{i\lambda}
  g(y)\,\frac{dy}{y}\,d\lambda
\]
sobre funciones de prueba, y defínase su realización ponderada por
\[
  \operatorname{Op}_{M,\delta}(r)
  :=
  \Phi_\delta^{-1}\operatorname{Op}_M(r)\Phi_\delta.
\]
Para fijar también la notación Wiener--Hopf, úsese la convención
\[
  (\mathcal Ff)(\lambda)
  :=
  \int_{\mathbb R}f(t)e^{-it\lambda}\,dt,
  \qquad
  (\mathcal F^{-1}g)(t)
  :=
  \frac1{2\pi}\int_{\mathbb R}g(\lambda)e^{it\lambda}\,d\lambda.
\]
Si $e_+$ denota extensión por cero desde $\mathbb R_+$ a $\mathbb R$,
$r_+$ denota restricción a $\mathbb R_+$ y $M_b$ es multiplicación por un
multiplicador de Fourier $b$, escríbase
\[
  W(b):=r_+\mathcal F^{-1}M_b\mathcal F e_+,
\]
siempre que este operador sea acotado en $X_\delta$.

Supóngase lo siguiente.
\begin{enumerate}
\item El operador $\operatorname{Op}_M(r_{1,1})$ se extiende a un operador
acotado en $L^p(\mathbb R_+,d\mu)$.
\item El símbolo reescalado
\[
  r_{1,1}^{\gamma}(x,\lambda)
  :=
  r_{1,1}(\gamma x,\lambda)
\]
pertenece a la misma clase admisible que $r_{1,1}$ y
$\operatorname{Op}_M(r_{1,1}^{\gamma})$ también es acotado. Si se denota la
variable espacial por $r$, esta convención se escribe
\[
  r_{1,1}^{\gamma}(r,\lambda)
  =
  r_{1,1}(\gamma r,\lambda),
\]
o, esquemáticamente,
\[
  r^{\gamma}(r,\lambda)=r(\gamma r,\lambda).
\]
\item Los operadores $W(b_{2,1})$, $W(b_{1,2})$, $A_{2,1}$,
$A_{1,2}$ y $R_{1,1}$ son acotados en $X_\delta$, y todos emplean los
mismos parámetros $p$, $\delta$, $\nu$ y la misma convención de dilatación.
\item Existen operadores compactos $K_{2,1}$, $K_{1,2}$ y $K_R$ en
$X_\delta$ tales que, preservando exactamente el orden escrito,
\begin{align*}
  A_{2,1}
  &=
  W(b_{2,1})V_\gamma^{(\nu)}+K_{2,1},
  \\
  A_{1,2}
  &=
  V_{\gamma^{-1}}^{(\nu)}W(b_{1,2})+K_{1,2},
  \\
  R_{1,1}
  &=
  \operatorname{Op}_{M,\delta}(r_{1,1})+K_R.
\end{align*}
En el caso exacto se toma $K_{2,1}=K_{1,2}=K_R=0$.
\end{enumerate}

Entonces se tiene la identidad exacta de covariancia
\[
  V_\gamma^{(\nu)}
  \operatorname{Op}_{M,\delta}(r_{1,1})
  V_{\gamma^{-1}}^{(\nu)}
  =
  \operatorname{Op}_{M,\delta}(r_{1,1}^{\gamma}),
\]
y, en consecuencia,
\[
  A_{2,1}R_{1,1}A_{1,2}
  \equiv
  W(b_{2,1})
  \operatorname{Op}_{M,\delta}(r_{1,1}^{\gamma})
  W(b_{1,2})
  \pmod{\mathcal K(X_\delta)}.
\]
Si las tres representaciones de la hipótesis 4 son exactas, la última
congruencia es una igualdad exacta.
\end{lemma}

\begin{proof}
Primero se verifica la normalización ponderada. Para $f\in X_\delta$ y
$y=\gamma x$,
\begin{align*}
  \|V_\gamma^{(\nu)}f\|_{X_\delta}^p
  &=
  \int_0^\infty
  \left|x^\delta\gamma^\nu f(\gamma x)\right|^p\,dx
  \\
  &=
  \gamma^{p\nu-p\delta-1}
  \int_0^\infty |y^\delta f(y)|^p\,dy
  \\
  &=
  \|f\|_{X_\delta}^p,
\end{align*}
porque $p\nu-p\delta-1=0$. Por tanto,
$V_\gamma^{(\nu)}$ es una isometría invertible y
\[
  \bigl(V_\gamma^{(\nu)}\bigr)^{-1}
  =
  V_{\gamma^{-1}}^{(\nu)}.
\]

Sea
\[
  (D_\gamma g)(x):=g(\gamma x)
\]
en $L^p(\mathbb R_+,d\mu)$. Para una función de prueba $g$,
\begin{align*}
  \bigl(
    \Phi_\delta
    V_\gamma^{(\nu)}
    \Phi_\delta^{-1}g
  \bigr)(x)
  &=
  x^\nu\gamma^\nu(\gamma x)^{-\nu}g(\gamma x)
  \\
  &=
  g(\gamma x)
  =
  (D_\gamma g)(x).
\end{align*}
Así,
\[
  \Phi_\delta V_\gamma^{(\nu)}\Phi_\delta^{-1}=D_\gamma,
  \qquad
  \Phi_\delta V_{\gamma^{-1}}^{(\nu)}\Phi_\delta^{-1}
  =D_{\gamma^{-1}}.
\]

Aplicando ahora la definición integral del Mellin PDO,
\begin{align*}
 &\bigl(
   D_\gamma\operatorname{Op}_M(r_{1,1})D_{\gamma^{-1}}g
  \bigr)(x)
  \\
 &\quad=
  \bigl(
    \operatorname{Op}_M(r_{1,1})D_{\gamma^{-1}}g
  \bigr)(\gamma x)
  \\
 &\quad=
  \frac1{2\pi}
  \int_{\mathbb R}\int_0^\infty
  r_{1,1}(\gamma x,\lambda)
  \left(\frac{\gamma x}{y}\right)^{i\lambda}
  g(y/\gamma)\,\frac{dy}{y}\,d\lambda.
\end{align*}
Con el cambio de variables $y=\gamma z$ se tiene $dy/y=dz/z$ y
\begin{align*}
 &\bigl(
   D_\gamma\operatorname{Op}_M(r_{1,1})D_{\gamma^{-1}}g
  \bigr)(x)
  \\
 &\quad=
  \frac1{2\pi}
  \int_{\mathbb R}\int_0^\infty
  r_{1,1}(\gamma x,\lambda)
  \left(\frac{x}{z}\right)^{i\lambda}
  g(z)\,\frac{dz}{z}\,d\lambda
  \\
 &\quad=
  \bigl(\operatorname{Op}_M(r_{1,1}^{\gamma})g\bigr)(x).
\end{align*}
La igualdad vale primero sobre el núcleo de funciones de prueba usado para
definir las integrales; las hipótesis de acotación permiten extenderla a todo
$L^p(\mathbb R_+,d\mu)$. Al conjugar por $\Phi_\delta^{-1}$ y
$\Phi_\delta$ se obtiene
\[
  V_\gamma^{(\nu)}
  \operatorname{Op}_{M,\delta}(r_{1,1})
  V_{\gamma^{-1}}^{(\nu)}
  =
  \operatorname{Op}_{M,\delta}(r_{1,1}^{\gamma}).
\]

Para propagar las equivalencias módulo compactos, póngase
\begin{align*}
  A'_{2,1}&:=W(b_{2,1})V_\gamma^{(\nu)},
  \\
  R'_{1,1}&:=\operatorname{Op}_{M,\delta}(r_{1,1}),
  \\
  A'_{1,2}&:=V_{\gamma^{-1}}^{(\nu)}W(b_{1,2}).
\end{align*}
La identidad telescópica
\begin{align*}
 &A_{2,1}R_{1,1}A_{1,2}
  -A'_{2,1}R'_{1,1}A'_{1,2}
  \\
 &\quad=
  (A_{2,1}-A'_{2,1})R_{1,1}A_{1,2}
  \\
 &\qquad
  +A'_{2,1}(R_{1,1}-R'_{1,1})A_{1,2}
  \\
 &\qquad
  +A'_{2,1}R'_{1,1}(A_{1,2}-A'_{1,2})
\end{align*}
muestra que la diferencia es compacta, porque $\mathcal K(X_\delta)$ es un
ideal bilateral y todos los factores restantes son acotados. Finalmente,
\begin{align*}
  A'_{2,1}R'_{1,1}A'_{1,2}
  &=
  W(b_{2,1})
  V_\gamma^{(\nu)}
  \operatorname{Op}_{M,\delta}(r_{1,1})
  V_{\gamma^{-1}}^{(\nu)}
  W(b_{1,2})
  \\
  &=
  W(b_{2,1})
  \operatorname{Op}_{M,\delta}(r_{1,1}^{\gamma})
  W(b_{1,2}),
\end{align*}
lo que prueba la congruencia. Si los tres residuos compactos son cero, todos
los pasos son igualdades exactas.
\end{proof}
```

## Estado lógico del lema

```text
dilation_covariance: EXACT_IDENTITY
dilation_cancellation: EXACT_IDENTITY
initial_block_representations: EQUALITY_MODULO_COMPACTS_OR_EXACT_BY_HYPOTHESIS
regularizer_mellin_representation: CONDITIONAL_ON_REGULARIZER_REPRESENTATION
resulting_schur_product_relation: WEAKEST_INPUT_RELATION
mixed_algebra_membership: UNPROVED_MIXED_ALGEBRA_MEMBERSHIP
fredholm_consequence: UNPROVED
```

La igualdad exacta del núcleo central no fortalece las relaciones de entrada.
Si cualquiera de las tres representaciones iniciales vale únicamente módulo
compactos, el resultado completo también vale únicamente módulo compactos.

## Advertencia sobre los bloques reales del artículo

El lema anterior **no aplica todavía** a los bloques que aparecen realmente en
la reducción del artículo. Las fórmulas auditadas son

\[
A_{2,1}
\simeq
-(\widetilde V_{\alpha_2}-G_2)W_{2,1}^{-},
\qquad
A_{1,2}
\simeq
(\widetilde V_{\alpha_1}-G_1)W_{1,2}^{+},
\]

y el regularizador de Schur elegido es

\[
A_{1,1}^{(-1)}
=
\mathcal T_{1,-}
\Phi_\delta^{-1}\operatorname{Op}_M(r_{1,1})
\Phi_\delta Z_1^{-1}.
\]

Así, el término real se reorganiza como

\[
C_{2,2,\mathrm{artículo}}^{(1)}
=
\mathcal B_{2,1}^{-}\mathcal R_1\mathcal B_{1,2}^{+},
\]

con

\[
\mathcal B_{2,1}^{-}
=
(\widetilde V_{\alpha_2}-G_2)
W_{2,1}^{-}\mathcal T_{1,-},
\]

\[
\mathcal B_{1,2}^{+}
=
Z_1^{-1}
(\widetilde V_{\alpha_1}-G_1)W_{1,2}^{+}.
\]

Todavía falta demostrar que estos dos factores, incluidos sus coeficientes,
normalizadores y localizadores `chi_infinity`, admiten las factorizaciones con
dilataciones inversas requeridas por el lema. La identidad exacta demostrada
para el núcleo auxiliar no localizado
`V_k W(b_{k,j}) V_j^{-1}` no cubre automáticamente estos factores.

Tampoco se ha demostrado que el sándwich final

\[
W(b_{2,1})
\operatorname{Op}_{M,\delta}(r_{1,1}^{\gamma})
W(b_{1,2})
\]

pertenezca a una clase mixta cerrada que permita aplicar un criterio de
Fredholm. Este borrador no prueba Fredholmness, no calcula un índice y no
resuelve el pivote de Schur.

## Nota editorial sobre el conflicto de signo

Debe distinguirse el producto crudo

\[
P_{\mathrm{raw}}
:=
A_{2,1}A_{1,1}^{(-1)}A_{1,2}
\]

del objeto llamado `C_{2,2}^{(1)}` en el artículo. El pivote se define por

\[
A_{2,2}^{(1)}
=
A_{2,2}-P_{\mathrm{raw}}.
\]

Como el modelo normalizado de `A_{2,1}` contiene un signo negativo, el artículo
denomina

\[
C_{2,2,\mathrm{artículo}}^{(1)}
:=
(\widetilde V_{\alpha_2}-G_2)
W_{2,1}^{-}
A_{1,1}^{(-1)}
(\widetilde V_{\alpha_1}-G_1)
W_{1,2}^{+},
\]

de modo que

\[
C_{2,2,\mathrm{artículo}}^{(1)}
\simeq
-P_{\mathrm{raw}}.
\]

El lema se formula para `P_raw`. Si se traduce su conclusión a la notación
editorial del artículo, debe introducirse el signo menos correspondiente. No
deben usarse ambos objetos bajo el mismo identificador.

## Decisión

```text
conditional_lemma: READY_AS_CONDITIONAL_DRAFT
application_to_actual_article_blocks: BLOCKED_WITH_PRECISE_OBLIGATION
mixed_algebra_membership: UNPROVED_MIXED_ALGEBRA_MEMBERSHIP
fredholm_or_index_claim: FORBIDDEN_AT_THIS_STAGE
```
