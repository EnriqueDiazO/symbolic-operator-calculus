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
exacto en el MVP.

Cuando se usa su kernel,

\[
(R_{1,1}f)(u)
=
\int_0^\infty R_{1,1}(u,v)f(v)\,dv.
\]

El kernel \(R_{1,1}(u,v)\) permanece formal. No se permiten cancelaciones
automaticas como \(A_{1,1}R_{1,1}=I\) ni \(R_{1,1}A_{1,1}=I\).

## Igualdad exacta y equivalencia modulo compactos

El MVP distingue dos tipos de relacion:

\[
A=B
\]

para igualdad exacta, y

\[
A\simeq B
\]

para equivalencia modulo compactos.

Una equivalencia modulo compactos no puede ser usada automaticamente como una
identidad exacta. Cualquier transformacion que dependa de \(\simeq\) debe estar
marcada estructuralmente como tal.

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

El kernel combinado del MVP es

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

Despues,

\[
(C_{2,2}^{(1)}f)(x)
=
\int_0^\infty C_{2,2}^{(1)}(x,y)f(y)\,dy.
\]

## LaTeX

LaTeX es solo una proyeccion de objetos simbolicos internos. Ninguna expresion
matematica del motor debe depender de una cadena LaTeX para su identidad,
expansion, composicion o comparacion estructural.
