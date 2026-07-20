# Mathematical Conventions

Este documento fija las convenciones matematicas del MVP. El codigo posterior
debe tratar estas reglas como parte de la especificacion, no como comentarios
decorativos.

## Convencion de Fourier

La transformada de Fourier se define por

\[
(\mathcal F f)(\lambda)
=
\int_{\mathbb R} f(t)e^{-it\lambda}\,dt.
\]

La transformada inversa se define por

\[
(\mathcal F^{-1} g)(t)
=
\frac{1}{2\pi}
\int_{\mathbb R} g(\lambda)e^{it\lambda}\,d\lambda.
\]

En particular, con \(d=c_2-c_1>0\),

\[
b^+_{1,2}(\lambda)
=
\chi_+(\lambda)e^{-d\lambda},
\qquad
b^-_{2,1}(\lambda)
=
\chi_-(\lambda)e^{d\lambda}.
\]

Bajo esta convencion,

\[
K^+_{1,2}(t)
=
\frac{1}{2\pi}
\int_0^\infty e^{-d\lambda}e^{it\lambda}\,d\lambda
=
\frac{1}{2\pi(d-it)},
\]

y

\[
K^-_{2,1}(t)
=
\frac{1}{2\pi}
\int_{-\infty}^{0} e^{d\lambda}e^{it\lambda}\,d\lambda
=
\frac{1}{2\pi(d+it)}.
\]

SymPy puede devolver estas integrales como `Piecewise` por condiciones de
convergencia. Si se usa `conds="none"`, debe hacerse solo despues de declarar
simbolicamente \(d>0\) y de documentar que se estan descartando condiciones que
ya estan cubiertas por la hipotesis del MVP.

## Orden de composicion

El producto de operadores se lee de izquierda a derecha como composicion que
actua sobre funciones desde la derecha:

\[
(ABCf)(x) = (A(B(Cf)))(x).
\]

Por tanto, el factor mas a la derecha actua primero. El orden escrito de los
factores es parte de la estructura matematica y no puede ser reordenado.

La expresion principal es

\[
(\widetilde V_{\alpha_2}-G_2)
W^-_{2,1}
R_{1,1}
(\widetilde V_{\alpha_1}-G_1)
W^+_{1,2}.
\]

Su expansion obligatoria contiene, en este orden, los cuatro terminos:

\[
\widetilde V_{\alpha_2}
W^-_{2,1}
R_{1,1}
\widetilde V_{\alpha_1}
W^+_{1,2},
\]

\[
-
\widetilde V_{\alpha_2}
W^-_{2,1}
R_{1,1}
G_1
W^+_{1,2},
\]

\[
-
G_2
W^-_{2,1}
R_{1,1}
\widetilde V_{\alpha_1}
W^+_{1,2},
\]

y

\[
+
G_2
W^-_{2,1}
R_{1,1}
G_1
W^+_{1,2}.
\]

Los signos son \(+, -, -, +\).

El orden anterior debe preservarse en la representacion estructural del motor.
No debe inferirse del orden canonico que una biblioteca simbolica pueda elegir
para imprimir o almacenar sumas.

## Identidad y cero estructurales

`Product(())` es la identidad estructural de composicion: representa la palabra
vacia del monoide de operadores. `LinearCombination(())` es el cero aditivo
estructural. Por ejemplo, para un operando valido \(f\),

\[
\operatorname{Product}(())f=f,
\qquad
\operatorname{LinearCombination}(())f=0.
\]

El objeto publico `I` es, en cambio, un atomo explicito de identidad usado en
formulas. Por ello, `Product(())` y `Product((I,))` son estructuralmente
distintos, aunque sus acciones sean equivalentes. Del mismo modo,
`Term(1, Product(()))` y `Term(1, I)` conservan estructuras distintas y la
segunda forma se normaliza a `Term(1, Product((I,)))`.

La igualdad del AST es igualdad estructural, no equivalencia algebraica
general. En particular, no se eliminan automaticamente factores `I` de
`Product((I, A))`, `Product((A, I))` ni `Product((I, I))`.

## Coeficientes del AST

Los coeficientes del AST son numeros concretos de Python de tipo `int`,
`float`, `complex` o `fractions.Fraction`. Los booleanos no son coeficientes,
aunque `bool` sea una subclase de `int`. Todo coeficiente admitido debe ser
finito; se rechazan `nan`, los infinitos y los complejos que tengan una parte
real o imaginaria no finita.

Los ceros de todos los tipos admitidos son validos y conservan su tipo dentro
de un `Term`; el filtrado estructural vigente los elimina al construir una
`LinearCombination`. Los complejos finitos son validos. Si su parte imaginaria
es cero, se presentan como numeros reales; en otro caso se renderizan como un
factor escalar LaTeX agrupado antes del producto de operadores.

Los escalares y expresiones SymPy no forman parte del contrato de coeficientes
del AST. Esta frontera mantiene coeficientes numericos concretos y hashables;
el soporte para parametros simbolicos requiere una fase arquitectonica
posterior.

## No conmutatividad

Todos los operadores del MVP son no conmutativos. En general,

\[
AB \ne BA.
\]

El motor no debe simplificar un producto de operadores usando conmutatividad a
menos que exista una regla matematica explicita, local y testeada. Para el MVP
no se propone ninguna regla de conmutacion entre operadores.

## Acciones basicas

Para un operador de multiplicacion,

\[
(G_k f)(x)=G_k(x)f(x).
\]

Para el shift transportado
\(\widetilde V_{\alpha_k}=\rho_k V_{\alpha_k}\), con
\(\alpha_k(x)=\gamma_k x\),

\[
(\widetilde V_{\alpha_k} f)(x)
=
\rho_k(x) f(\gamma_k x).
\]

Un operador integral formal actua como

\[
(Kf)(x)
=
\int_0^\infty K(x,y)f(y)\,dy.
\]

Los operadores localizados de Wiener-Hopf se representan mediante

\[
L^+_{1,2}(x,y)
=
\chi_\infty(x)K^+_{1,2}(x-y)\chi_\infty(y),
\]

y

\[
L^-_{2,1}(x,y)
=
\chi_\infty(x)K^-_{2,1}(x-y)\chi_\infty(y).
\]

## Regularizador formal

El simbolo

\[
R_{1,1}:=A_{1,1}^{(-1)}
\]

denota un regularizador formal modulo operadores compactos. No es un inverso
exacto en el MVP y no posee automaticamente un kernel.

La expresion integral

\[
(R_{1,1}f)(u)
=
\int_0^\infty R_{1,1}(u,v)f(v)\,dv.
\]
solo se puede construir cuando el usuario aporta explicitamente una
`KernelRepresentation`. Este objeto identifica el regularizador, la expresion
del kernel, sus variables, el dominio de integracion, su estatus semantico y
las hipotesis o evidencia externas. Sin ese objeto, las APIs de acciones y de
Schur lanzan `KernelRepresentationRequiredError`; no se crea automaticamente
una funcion simbolica `R11(u, v)`.

Una `KernelRepresentation` formal o asumida no prueba que el kernel exista como
funcion ordinaria. Incluso con una representacion certificada externamente, una
`Integral` de SymPy no prueba convergencia, acotacion ni existencia del operador.
No se permiten cancelaciones automaticas como \(A_{1,1}R_{1,1}=I\) ni
\(R_{1,1}A_{1,1}=I\).

## Tipos de relacion semantica

El modulo `semantics.py` distingue cuatro tipos inmutables y no convertibles:

- `ExactIdentity`, con evidencia explicita;
- `FormalIdentity`, que no afirma validez analitica o exacta;
- `ModCompactEquivalence`, certificada solo si el usuario aporta evidencia;
- `ApproximateEquality`, con tolerancia, norma o criterio y residual explicitos.

En particular, el MVP distingue

\[
A=B
\]

para igualdad exacta, y

\[
A\simeq B
\]

para equivalencia modulo compactos.

Una equivalencia modulo compactos no puede ser usada automaticamente como una
identidad exacta. `ModCompactRelation` y `ModCompactSchurRelation` conservan sus
endpoints historicos, pero exponen una `ModCompactEquivalence` no certificada
por defecto. Los campos de espacio, ideal, residual y evidencia son datos
aportados; el codigo no infiere que un residual sea compacto.

## Invariantes de operadores Wiener–Hopf relativos

El dominio Wiener–Hopf relativo usa indices enteros estrictamente positivos y
solo representa pares off-diagonal, de modo que \(k\ne j\). Sus escalas son
expresiones SymPy declaradas reales, positivas y finitas. La dilatacion relativa

\[
\beta_{k,j}(x)=\frac{\gamma_k}{\gamma_j}x
\]

almacena las dos escalas de procedencia y valida algebraicamente ese cociente.

Para una misma factorizacion L1 se admiten tres productos canonicos ordenados:

\[
V_{\alpha_k}W(b_{k,j})V_{\alpha_j}^{-1},\qquad
V_{\beta_{k,j}}W(b_{k,j}^{\mathrm L}),\qquad
W(b_{k,j}^{\mathrm R})V_{\beta_{k,j}}.
\]

El orden de estos factores, su par \((k,j)\), sus escalas y la procedencia de
los tres simbolos forman parte de la estructura validada. En
`RelativeWienerHopfIdentity`, `exact_relations` contiene dos objetos
`ExactIdentity` inmutables, cada uno con evidencia de la forma canonica y de la
procedencia L1 validadas. No se usa un booleano `exact`. Estas relaciones no
afirman por si solas que los tres productos hayan sido aplicados
independientemente a una funcion. Esa semantica ejecutable y la verificacion
algebraica independiente de correspondencias quedan fuera de A3.1 y
corresponden a A3.2/A3.3.

## Semántica ejecutable de productos Wiener–Hopf relativos

Las dilataciones relativas no llevan peso. Para una expresion \(q\), sus
acciones directa, inversa y relativa son, respectivamente,

\[
(V_{\alpha_k}q)(x)=q(\gamma_kx),\qquad
(V_{\alpha_j}^{-1}q)(x)=q(x/\gamma_j),\qquad
(V_{\beta_{k,j}}q)(x)=q((\gamma_k/\gamma_j)x).
\]

Un factor Wiener–Hopf usa el kernel de convolucion asociado a su placement:

\[
(W(h)q)(x)=\int_0^\infty h(x-y)q(y)\,dy.
\]

Los productos se aplican factor por factor desde la derecha. Cada una de las
tres formas canonicas conserva una accion directa propia. La normalizacion
realiza el cambio de variable acotado que corresponde a esa forma y obtiene en
los tres casos

\[
K_{k,j}^{\mathrm{rel}}(x,y)
=\gamma_jh_{k,j}(\gamma_kx-\gamma_jy).
\]

En el producto original, el factor \(\gamma_j\) surge del cambio de variable
inducido por \(V_{\alpha_j}^{-1}\); no forma parte de la accion atomica de una
dilatacion. La igualdad de los tres kernels se comprueba algebraicamente y se
registra como evidencia ejecutable derivada. Esta evidencia es distinta de
`RelativeWienerHopfIdentity.exact_relations`, que registra identidades
estructurales tipadas y no depende de las variables elegidas para aplicar los
productos.

## Correspondencias Fourier de símbolos relativos

Bajo la convención Fourier del proyecto, para una escala positiva \(a\),

\[
h_a(t)=a h(at)
\quad\Longleftrightarrow\quad
b_a(\lambda)=b(\lambda/a).
\]

Por tanto, las colocaciones izquierda y derecha del operador relativo cumplen

\[
b^{\mathrm L}_{k,j}(\lambda)
=b_{k,j}(\lambda/\gamma_j),
\qquad
b^{\mathrm R}_{k,j}(\lambda)
=b_{k,j}(\lambda/\gamma_k).
\]

Estas correspondencias se verifican reconstruyendo los dos símbolos desde el
símbolo original y los dos kernels escalados desde el kernel original. Los
resultados calculados se comparan algebraicamente con los campos almacenados
por L1. La evidencia `correspondences_verified` es inmutable y derivada; no es
un alias de esos campos.

Las tres evidencias tienen significados distintos: `identity.exact_relations`
registra la coherencia estructural tipada de los productos, `actions_verified` certifica la
igualdad de sus acciones normalizadas y `correspondences_verified` certifica
las correspondencias Fourier reconstruidas independientemente.

## Intervalos de integracion

Salvo que se indique lo contrario, los operadores integrales del MVP integran
sobre

\[
(0,\infty).
\]

Las integrales de Fourier usan \(\mathbb R\), \((0,\infty)\) o
\((-\infty,0)\) segun la funcion caracteristica \(\chi_+\) o \(\chi_-\).

## Variables mudas

La composicion debe generar variables mudas frescas. Las reglas minimas son:

- no reutilizar una variable libre como variable ligada;
- no reutilizar una variable ya ligada en la misma expresion;
- aceptar un conjunto explicito de nombres prohibidos;
- preferir nombres legibles como \(u\), \(v\), \(w\), \(z\), y despues
  variantes indexadas si esos nombres ya estan ocupados;
- verificar por tests que no hay captura accidental de \(x\), \(y\), \(u\),
  \(v\) u otras variables presentes.

Para operadores integrales \(A\), \(B\), \(C\), la composicion debe producir
estructuralmente

\[
K_{ABC}(x,y)
=
\int_0^\infty\int_0^\infty
K_A(x,u)K_B(u,v)K_C(v,y)\,dv\,du,
\]

respetando el orden \(K_A K_B K_C\).

## Kernel combinado

Se definen

\[
M_{1,2}(v,y)
=
\rho_1(v)L^+_{1,2}(\gamma_1v,y)
-
G_1(v)L^+_{1,2}(v,y),
\]

y

\[
M_{2,1}(x,u)
=
\rho_2(x)L^-_{2,1}(\gamma_2x,u)
-
G_2(x)L^-_{2,1}(x,u).
\]

Si, y solo si, se ha aportado una `KernelRepresentation` de \(R_{1,1}\),
el kernel combinado formal del MVP se representa como

\[
C_{2,2}^{(1)}(x,y)
=
\int_0^\infty
\int_0^\infty
M_{2,1}(x,u)
R_{1,1}(u,v)
M_{1,2}(v,y)
\,dv\,du.
\]

Despues puede construirse formalmente

\[
(C_{2,2}^{(1)}f)(x)
=
\int_0^\infty C_{2,2}^{(1)}(x,y)f(y)\,dy.
\]

Ambas integrales quedan envueltas en `KernelAnnotatedExpression`, que conserva
la representacion, el estatus y las hipotesis aportadas. La reduccion algebraica
de Schur \(A_{22}-A_{21}R_{11}A_{12}\) sigue existiendo sin kernel; lo que se
prohibe es su materializacion integral automatica. Ninguna de estas expresiones
demuestra compacidad, convergencia, acotacion o propiedades Fredholm.

## LaTeX

LaTeX es solo una proyeccion de objetos simbolicos internos. Ninguna expresion
matematica del motor debe depender de una cadena LaTeX para su identidad,
expansion, composicion o comparacion estructural.
