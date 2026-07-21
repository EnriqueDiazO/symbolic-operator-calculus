# Cancelación de dilataciones en la primera corrección de Schur

## Alcance y veredicto

Este documento separa dos afirmaciones que no deben confundirse.

1. El lema abstracto

   \[
   \mathsf V_\gamma\operatorname{Op}_M(a)
   \mathsf V_{\gamma^{-1}}
   =\operatorname{Op}_M(a^\gamma),
   \qquad
   a^\gamma(r,\lambda)=a(\gamma r,\lambda),
   \]

   es una **identidad operatorial exacta** bajo las convenciones del artículo.

2. El artículo no ha reducido sus bloques reales a

   \[
   A_{2,1}\simeq W(b_{2,1})\mathsf V_\gamma,
   \qquad
   A_{1,2}\simeq
   \mathsf V_{\gamma^{-1}}W(b_{1,2}).
   \]

   Por ello, el lema no puede aplicarse todavía a la corrección real. La
   decisión para el artículo es `BLOCKED_WITH_PRECISE_OBLIGATION`, no una
   certificación del pivote.

Las cuatro capas lógicas son:

| Capa | Estado |
|---|---|
| Conjugación de la Mellin PDO por el par inverso | `EXACT_IDENTITY` |
| Sustitución de los bloques por modelos, si se prueban los modelos requeridos | `EQUALITY_MODULO_COMPACTS` |
| Representación Mellin del regularizador abstracto suministrado a la API | `CONDITIONAL_ON_REGULARIZER_REPRESENTATION` |
| Pertenencia del sándwich final a una álgebra Fredholm | `UNPROVED_MIXED_ALGEBRA_MEMBERSHIP` |

## 1. Espacio ponderado y convención de dilatación

El artículo usa el peso multiplicativo

\[
w_\delta(x)=x^\delta,
\qquad
\|f\|_{p,\delta}^p
=\int_0^\infty |x^\delta f(x)|^p\,dx,
\]

no la convención alternativa
\(\int |f(x)|^p x^\delta\,dx\). Véanse
`sections/03_transport_and_kernels.tex:241-262`, `:311-329` y
`:1084-1099` en el repositorio del artículo.

Para registrar sin ambigüedad tanto la dilatación desnuda como la normalizada,
fijamos una potencia real \(\nu\) y definimos

\[
[\mathsf V_\gamma^{(\nu)}f](x)
:=\gamma^\nu f(\gamma x),
\qquad \gamma>0.
\]

La notación del artículo es

\[
V_\gamma=\mathsf V_\gamma^{(0)},
\qquad
U_\gamma=\mathsf V_\gamma^{(\kappa_\delta)},
\qquad
\kappa_\delta:=\delta+\frac1p.
\]

Esta distinción es necesaria: las factorizaciones Wiener--Hopf relativas de la
Sección 5 usan \(V_\gamma\), mientras que la isometría del espacio ponderado es
\(U_\gamma\).

### 1.1 Cálculo explícito de la norma

Con \(y=\gamma x\),

\[
\begin{aligned}
\|\mathsf V_\gamma^{(\nu)}f\|_{p,\delta}^p
&=\int_0^\infty
  |x^\delta\gamma^\nu f(\gamma x)|^p\,dx\\
&=\gamma^{p\nu-p\delta-1}
  \int_0^\infty |y^\delta f(y)|^p\,dy.
\end{aligned}
\]

Por tanto,

\[
\|\mathsf V_\gamma^{(\nu)}\|
=\gamma^{\nu-\delta-1/p}.
\]

La única potencia que hace isométrica esta familia para todo \(\gamma>0\) es

\[
\boxed{\nu=\kappa_\delta=\delta+\frac1p.}
\]

Esto reproduce la derivación del artículo en
`sections/03_transport_and_kernels.tex:1371-1467`; no se ha supuesto el
exponente.

### 1.2 Inversa, identidad y composición

Para la misma potencia \(\nu\),

\[
\mathsf V_1^{(\nu)}=I,
\qquad
\mathsf V_\gamma^{(\nu)}
\mathsf V_\eta^{(\nu)}
=\mathsf V_{\gamma\eta}^{(\nu)},
\qquad
(\mathsf V_\gamma^{(\nu)})^{-1}
=\mathsf V_{\gamma^{-1}}^{(\nu)}.
\]

En particular, el factor de normalización se cancela exactamente en una
conjugación por dilataciones inversas. La API exige que ambos factores usen la
misma \(\nu\), el mismo \(p\), el mismo peso, el mismo espacio y la misma
convención de acción.

### 1.3 Transporte a la medida de Haar

Con

\[
[\Phi_\delta f](x)=x^{\kappa_\delta}f(x),
\]

se tiene una isometría de
\(L^p(\mathbb R_+,w_\delta)\) a
\(L^p(\mathbb R_+,d\mu)\), \(d\mu(x)=dx/x\), y

\[
\Phi_\delta U_\gamma\Phi_\delta^{-1}=D_\gamma,
\qquad
[D_\gamma g](x)=g(\gamma x).
\]

Así, las identidades siguientes pueden demostrarse en el espacio de Haar y
transportarse de vuelta sin escalares residuales.

## 2. Covariancia de multiplicaciones

Sea \([M_af](x)=a(x)f(x)\). Entonces, punto por punto,

\[
\begin{aligned}
[\mathsf V_\gamma^{(\nu)}M_a
  \mathsf V_{\gamma^{-1}}^{(\nu)}f](x)
&=\gamma^\nu a(\gamma x)
  \gamma^{-\nu}f(x)\\
&=a(\gamma x)f(x).
\end{aligned}
\]

Por consiguiente,

\[
\boxed{
\mathsf V_\gamma^{(\nu)}M_a
\mathsf V_{\gamma^{-1}}^{(\nu)}
=M_{a^\gamma},
\qquad a^\gamma(x)=a(\gamma x).
}
\]

La escala \(\gamma x\), y no \(x/\gamma\), viene fijada por el orden escrito.
El orden inverso produce
\(\mathsf V_{\gamma^{-1}}M_a\mathsf V_\gamma=M_{a(x/\gamma)}\).

## 3. Covariancia Wiener--Hopf

### 3.1 Convención de Fourier y truncación

El artículo fija

\[
(\mathcal Ff)(\lambda)
=\int_{\mathbb R}f(t)e^{-it\lambda}\,dt,
\qquad
(\mathcal F^{-1}g)(t)
=\frac1{2\pi}\int_{\mathbb R}g(\lambda)e^{it\lambda}\,d\lambda,
\]

y

\[
W(b)=r_+\mathcal F^{-1}b\mathcal F e_+,
\]

donde \(e_+\) extiende por cero y \(r_+\) restringe a la semirrecta. Véase
`sections/04_target_theorem53.tex:17-48`.

Para un núcleo de convolución inicialmente integrable,

\[
[W(b)f](x)=\int_0^\infty h(x-y)f(y)\,dy,
\qquad h=\mathcal F^{-1}b.
\]

La positividad de \(\gamma\) conserva \((0,\infty)\), de modo que extensión,
restricción y cambio de escala son compatibles. Para multiplicadores
generales, el cálculo sobre una clase densa se extiende por acotación; no se
reemplaza esta definición por una convención de transformada de SymPy.

### 3.2 Cálculo del núcleo

Aplicando los factores desde la derecha y poniendo \(y=\gamma z\),

\[
\begin{aligned}
[\mathsf V_\gamma^{(\nu)}W(b)
  \mathsf V_{\gamma^{-1}}^{(\nu)}f](x)
&=\gamma^\nu\int_0^\infty
  h(\gamma x-y)\gamma^{-\nu}f(y/\gamma)\,dy\\
&=\int_0^\infty \gamma h(\gamma(x-z))f(z)\,dz.
\end{aligned}
\]

El núcleo nuevo es \(h_\gamma(t)=\gamma h(\gamma t)\). Con la convención
anterior,

\[
\begin{aligned}
(\mathcal Fh_\gamma)(\lambda)
&=\int_{\mathbb R}\gamma h(\gamma t)e^{-it\lambda}\,dt\\
&=\int_{\mathbb R}h(u)e^{-iu\lambda/\gamma}\,du
=b(\lambda/\gamma).
\end{aligned}
\]

No aparece otro Jacobiano, signo ni factor \(2\pi\). Así,

\[
\boxed{
\mathsf V_\gamma^{(\nu)}W(b)
\mathsf V_{\gamma^{-1}}^{(\nu)}
=W(b_\gamma),
\qquad b_\gamma(\lambda)=b(\lambda/\gamma).
}
\]

### 3.3 Formas con la dilatación a izquierda o derecha

Multiplicando la identidad exacta por la dilatación correspondiente se obtiene

\[
\boxed{
\mathsf V_\gamma W(b)
=W\bigl(b(\cdot/\gamma)\bigr)\mathsf V_\gamma,
}
\]

y, de manera equivalente,

\[
\boxed{
W(b)\mathsf V_\gamma
=\mathsf V_\gamma W\bigl(b(\gamma\,\cdot)\bigr).
}
\]

Por ello no es lícito mover \(\mathsf V_\gamma\) a través de \(W(b)\) sin
reescalar el símbolo.

### 3.4 Localización

Si el factor real es \(M_{\chi_\infty}W(b)M_{\chi_\infty}\), entonces

\[
\mathsf V_\gamma M_{\chi_\infty}
\mathsf V_{\gamma^{-1}}
=M_{\chi_\infty(\gamma\,\cdot)}.
\]

La covariancia del Wiener--Hopf desnudo no elimina estos cortes. Cualquier
regla que los atraviese debe conservar los cortes reescalados o aportar una
equivalencia módulo compactos específica. La API rechaza una localización no
controlada.

## 4. Convolución de Mellin

En \(L^p(\mathbb R_+,d\mu)\), una convolución de Mellin es el caso de símbolo
sin dependencia radial:

\[
M^0(c)=\operatorname{Op}_M(c(\lambda)).
\]

La fórmula de la sección siguiente no modifica \(\lambda\), y por tanto

\[
\boxed{
D_\gamma M^0(c)D_{\gamma^{-1}}=M^0(c)
}
\]

exactamente. La misma identidad vale para \(U_\gamma\) en el espacio
ponderado después de conjugar por \(\Phi_\delta\).

## 5. Mellin PDO: derivación desde la integral

La cuantización del artículo es

\[
[\operatorname{Op}_M(a)f](x)
=\frac1{2\pi}
 \int_{\mathbb R}\int_0^\infty
 a(x,\lambda)
 \left(\frac{x}{y}\right)^{i\lambda}
 f(y)\,\frac{dy}{y}\,d\lambda.
\]

Tomamos primero \(f\in C_c^\infty(\mathbb R_+)\) y truncamos simétricamente
la integral en \(\lambda\). El cambio de variable siguiente se hace en la
integral interior; no requiere intercambiar integraciones.

Como \([D_{\gamma^{-1}}f](y)=f(y/\gamma)\),

\[
\begin{aligned}
&[D_\gamma\operatorname{Op}_M(a)D_{\gamma^{-1}}f](x)\\
&=\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
 a(\gamma x,\lambda)
 \left(\frac{\gamma x}{y}\right)^{i\lambda}
 f(y/\gamma)\,\frac{dy}{y}\,d\lambda.
\end{aligned}
\]

Con \(y=\gamma z\), se tiene \(dy/y=dz/z\) y

\[
\left(\frac{\gamma x}{\gamma z}\right)^{i\lambda}
=\left(\frac{x}{z}\right)^{i\lambda}.
\]

Luego

\[
\begin{aligned}
[D_\gamma\operatorname{Op}_M(a)D_{\gamma^{-1}}f](x)
&=\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
 a(\gamma x,\lambda)
 \left(\frac{x}{z}\right)^{i\lambda}
 f(z)\,\frac{dz}{z}\,d\lambda\\
&=[\operatorname{Op}_M(a^\gamma)f](x),
\end{aligned}
\]

donde

\[
\boxed{a^\gamma(x,\lambda)=a(\gamma x,\lambda).}
\]

Si ambos lados son acotados, la igualdad sobre la clase densa tiene una única
extensión acotada. Al volver al espacio ponderado, los escalares de
\(U_\gamma\) y \(U_{\gamma^{-1}}\) se cancelan. Por tanto,

\[
\boxed{
U_\gamma\,\Phi_\delta^{-1}\operatorname{Op}_M(a)\Phi_\delta
U_{\gamma^{-1}}
=\Phi_\delta^{-1}\operatorname{Op}_M(a^\gamma)\Phi_\delta.
}
\]

La fórmula existente de Fase R,
\(D_\gamma^{-1}\operatorname{Op}_M(a)D_\gamma
=\operatorname{Op}_M(a(x/\gamma,\lambda))\), corresponde al orden inverso y
es consistente con esta derivación.

## 6. Estabilidad de la clase de símbolos

Sea
\(S_\gamma a(x,\lambda)=a(\gamma x,\lambda)\), \(\gamma>0\). Para las
clases escalares \(\mathcal E\) y
\(\widetilde{\mathcal E}(\mathbb R_+,V(\mathbb R))\) usadas por el artículo,
la estabilidad se prueba directamente:

1. \(x\mapsto\gamma x\) es un homeomorfismo de \(\mathbb R_+\), y

   \[
   \sup_{x>0}\|a(\gamma x,\cdot)\|_V
   =\sup_{z>0}\|a(z,\cdot)\|_V.
   \]

2. La continuidad con valores en \(V\) se conserva por composición.

3. Los módulos de oscilación radial satisfacen

   \[
   \operatorname{cm}^{\infty}_t(S_\gamma a)
   =\operatorname{cm}^{\infty}_{\gamma t}(a),
   \]

   de modo que se conservan los límites en \(0\) y \(\infty\).

4. Las traslaciones en \(\lambda\) conmutan con \(S_\gamma\), y sus supremos
   uniformes no cambian tras \(z=\gamma x\).

5. Las colas uniformes de \(\partial_\lambda a\) son idénticas después del
   mismo cambio radial.

6. La escala radial induce un automorfismo del álgebra lentamente oscilante y
   conserva cada fibra sobre \(0\) e \(\infty\).

Por tanto, para estas clases escalares,

\[
a\in\widetilde{\mathcal E}
\Longrightarrow
S_\gamma a\in\widetilde{\mathcal E}.
\]

Esta prueba no implica que productos con el multiplicador oscilatorio
\(\gamma^{i\lambda}\) pertenezcan a la misma clase, ni establece una versión
matricial. La implementación sólo conserva el nombre de la clase cuando se
aporta evidencia explícita de estabilidad radial.

## 7. Lema condicional abstracto

### Conditional dilation cancellation in the first Schur correction

Sean \(X=L^p(\mathbb R_+,w_\delta)\), \(\gamma>0\) y una potencia común
\(\nu\). Supóngase que:

1. todos los factores son acotados en el mismo espacio \(X\);
2. con relación \(\mathrel{\mathcal R}_L\), exacta o módulo compactos,

   \[
   A_{2,1}\mathrel{\mathcal R}_L
   W(b_{2,1})\mathsf V_\gamma^{(\nu)};
   \]

3. con relación \(\mathrel{\mathcal R}_R\), exacta o módulo compactos,

   \[
   A_{1,2}\mathrel{\mathcal R}_R
   \mathsf V_{\gamma^{-1}}^{(\nu)}W(b_{1,2});
   \]

4. el regularizador posee, de manera certificada o como hipótesis explícita,

   \[
   R_{1,1}\mathrel{\mathcal R}_M
   \Phi_\delta^{-1}\operatorname{Op}_M(r_{1,1})\Phi_\delta;
   \]

5. la clase declarada de \(r_{1,1}\) es estable bajo
   \(r\mapsto\gamma r\);
6. no hay entre las dos dilataciones un corte, proyección, multiplicador u
   otro factor además de la Mellin PDO representada;
7. si alguna relación de entrada es módulo compactos, se usa explícitamente
   que los compactos forman un ideal bilateral en \(\mathcal B(X)\).

Entonces la parte de modelo satisface exactamente

\[
\begin{aligned}
&W(b_{2,1})\mathsf V_\gamma^{(\nu)}
 \Phi_\delta^{-1}\operatorname{Op}_M(r_{1,1})\Phi_\delta
 \mathsf V_{\gamma^{-1}}^{(\nu)}W(b_{1,2})\\
&\qquad=
W(b_{2,1})
\Phi_\delta^{-1}\operatorname{Op}_M(r_{1,1}^\gamma)\Phi_\delta
W(b_{1,2}),
\end{aligned}
\]

con

\[
r_{1,1}^\gamma(r,\lambda)=r_{1,1}(\gamma r,\lambda).
\]

Al sustituir los modelos en el producto de Schur, la relación de salida es la
más débil de \(\mathcal R_L,\mathcal R_M,\mathcal R_R\). En particular, si
alguna es módulo compactos, la conclusión global es sólo

\[
A_{2,1}R_{1,1}A_{1,2}
\simeq
W(b_{2,1})
\Phi_\delta^{-1}\operatorname{Op}_M(r_{1,1}^\gamma)\Phi_\delta
W(b_{1,2}).
\]

Si la representación Mellin del regularizador es una hipótesis, toda esta
conclusión queda marcada
`CONDITIONAL_ON_REGULARIZER_REPRESENTATION`, aunque la cancelación interna sea
exacta.

### Prueba de la propagación módulo compactos

Escribamos

\[
A_{2,1}=L+K_L,
\qquad
R_{1,1}=M+K_M,
\qquad
A_{1,2}=R+K_R,
\]

donde sólo se incluyen los \(K\) correspondientes a relaciones módulo
compactos. Al expandir
\((L+K_L)(M+K_M)(R+K_R)-LMR\), cada sumando contiene al menos un factor
compacto y los demás son acotados. Por la propiedad de ideal bilateral y por
estabilidad bajo sumas finitas, el residual es compacto. Esta prueba no
convierte ninguna equivalencia de entrada en igualdad exacta.

## 8. Por qué el lema no se aplica todavía al artículo

El artículo demuestra solamente

\[
A_{1,2}\simeq
(\widetilde V_{\alpha_1}-G_1)W_{1,2}^+,
\qquad
A_{2,1}\simeq
-(\widetilde V_{\alpha_2}-G_2)W_{2,1}^-,
\]

en `sections/06_shur_reduction.tex:44-95`. El regularizador del bloque usado
en Schur es

\[
A_{1,1}^{(-1)}
=T_{1,-}\mathcal R_1Z_1^{-1},
\qquad
\mathcal R_1
=\Phi_\delta^{-1}\operatorname{Op}_M(r_{1,1})\Phi_\delta,
\]

en `sections/06_shur_reduction.tex:161-195`. Por tanto, la palabra real
contiene

\[
(\widetilde V_{\alpha_2}-G_2)W_{2,1}^-
T_{1,-}\mathcal R_1Z_1^{-1}
(\widetilde V_{\alpha_1}-G_1)W_{1,2}^+,
\]

no el patrón contiguo
\(W V_\gamma\operatorname{Op}_M(r)V_{\gamma^{-1}}W\). Los cortes
\(\chi_\infty\), los coeficientes, las proyecciones de \(T_{1,-}\) y
\(Z_1^{-1}\) tampoco pueden cruzarse sin reglas adicionales.

Además, el artículo define

\[
C_{2,2,\mathrm{artículo}}^{(1)}
\simeq-A_{2,1}A_{1,1}^{(-1)}A_{1,2},
\]

mientras que el encargo denomina \(C_{2,2}^{(1)}\) al producto sin el signo
menos. Ambos objetos se conservan separados.

## 9. Obligaciones abiertas

1. Probar una normalización de los bloques reales que produzca exactamente el
   par adyacente
   \(W(b_{2,1})\mathsf V_\gamma\) y
   \(\mathsf V_{\gamma^{-1}}W(b_{1,2})\), o demostrar otra regla que controle
   todos los factores interpuestos.
2. Controlar los localizadores \(\chi_\infty\) bajo escala radial, con una
   igualdad exacta que los conserve o una equivalencia módulo compactos
   certificada.
3. Identificar qué parte del regularizador completo
   \(T_{1,-}\mathcal R_1Z_1^{-1}\) queda entre el par de dilataciones.
4. Probar la pertenencia de
   \(W(b_{2,1})\operatorname{Op}_M(r_{1,1}^\gamma)W(b_{1,2})\) a una álgebra
   cerrada adecuada.
5. Construir después, y no antes, el símbolo Fredholm de la corrección y el
   cálculo del pivote.
6. Formular y probar por separado cualquier extensión matricial.

Nada de lo anterior se resuelve por la identidad exacta de covariancia.
