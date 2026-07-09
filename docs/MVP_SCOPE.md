# MVP Scope

Fecha de auditoria: 2026-07-09.

Este repositorio inicia un proyecto de investigacion matematica computacional
para experimentar con calculo simbolico de composiciones de operadores
integrales no conmutativos. El primer MVP no sera un CAS general: solo debe
modelar, expandir, componer y renderizar el caso concreto de
\(C_{2,2}^{(1)}\).

## Auditoria inicial

- Repositorio: vacio salvo por `.git`.
- Rama actual: `main`.
- Commits: ninguno todavia.
- Python disponible: `python3` es Python 3.10.12.
- `python` no esta configurado directamente en pyenv.
- `pip`: 22.0.2 para Python 3.10.
- SymPy local: 1.14.0.
- pytest local: 9.0.3.

No se hizo commit, no se creo paquete Python y no se implemento el motor.

## Fuentes revisadas

- SymPy 1.14.0 documenta simbolos no conmutativos mediante
  `Symbol(..., commutative=False)` y preservacion del orden de factores no
  conmutativos:
  https://docs.sympy.org/latest/tutorials/intro-tutorial/manipulation.html
- La guia de assumptions de SymPy distingue simbolos conmutativos y no
  conmutativos:
  https://docs.sympy.org/latest/guides/assumptions.html
- SymPy permite integrales no evaluadas con `Integral` y evaluacion posterior
  con `.doit()`:
  https://docs.sympy.org/latest/tutorials/intro-tutorial/calculus.html
- La documentacion de integrales cubre integrales definidas con limites:
  https://docs.sympy.org/latest/modules/integrals/integrals.html
- SymPy incluye salida LaTeX dentro de su sistema de impresion:
  https://docs.sympy.org/latest/tutorials/intro-tutorial/printing.html
- PyPI y la documentacion publica muestran SymPy 1.14.0 como version actual:
  https://pypi.org/project/sympy/

## Capacidades observadas de SymPy

Pruebas exploratorias locales, sin escribir codigo de proyecto:

- `A*B == B*A` es `False` cuando `A` y `B` son simbolos no conmutativos.
- `expand((A - B)*(C - D))` preserva el orden relativo de factores no
  conmutativos.
- SymPy puede ordenar los terminos de una suma expandida de forma canonica; por
  tanto, el orden de los cuatro terminos del MVP no debe depender de la
  representacion `Add` cruda.
- `Integral(K(x, y)*f(y), (y, 0, oo))` representa una integral formal sobre
  \((0,\infty)\).
- Para \(d>0\), `integrate(exp(-d*l)*exp(I*t*l), (l, 0, oo))` devuelve una
  expresion `Piecewise` por defecto.
- Con la misma hipotesis y `conds="none"`, SymPy devuelve una forma equivalente
  a \(1/(d-it)\); analogamente, la integral sobre \((-\infty,0)\) devuelve una
  forma equivalente a \(1/(d+it)\).

Decision: SymPy se usara para expresiones escalares, funciones simbolicas,
integrales formales, simplificaciones locales y LaTeX. El orden de operadores
y el orden de terminos de la expansion no se confiaran a `Mul` o `Add` como
representacion principal del motor; se propone una representacion estructural
propia y pequena para productos no conmutativos y combinaciones lineales.

## Objetivo del MVP

Representar la expresion

\[
(\widetilde V_{\alpha_2}-G_2)
W^-_{2,1}
R_{1,1}
(\widetilde V_{\alpha_1}-G_1)
W^+_{1,2},
\qquad R_{1,1}=A_{1,1}^{(-1)},
\]

expandirla en cuatro terminos con signos \(+, -, -, +\), preservar el orden
exacto de factores, convertir la accion sobre una funcion en integrales
anidadas y construir el kernel combinado

\[
C_{2,2}^{(1)}(x,y)
=
\int_0^\infty\int_0^\infty
M_{2,1}(x,u) R_{1,1}(u,v) M_{1,2}(v,y)\,dv\,du.
\]

## Fuera de alcance

- CAS general.
- Teoria completa del articulo.
- Calculo del kernel formal \(R_{1,1}\).
- Ideal completo de operadores compactos.
- Interfaz grafica, API web o base de datos.
- Mathematica, Wolfram Engine, NCAlgebra basado en Wolfram o software
  propietario.
- SageMath, GAP, Singular, OSCAR, PyLops, NumPy, SciPy y Matplotlib en el MVP.
- Reglas globales de simplificacion o monkey patching de SymPy.

## Dependencias propuestas

Dependencias de ejecucion:

- Python >= 3.10.
- SymPy >= 1.14, < 2.

Dependencias de desarrollo:

- pytest >= 9, < 10.

No se propone anadir NumPy, SciPy, Matplotlib ni herramientas web. `mpmath`
entrara solo como dependencia transitiva de SymPy.

## `pyproject.toml` propuesto

Este archivo se propone para la siguiente fase, pero no se crea todavia.

```toml
[build-system]
requires = ["setuptools>=69"]
build-backend = "setuptools.build_meta"

[project]
name = "symbolic-operator-calculus"
version = "0.1.0"
description = "Minimal symbolic experiments for noncommutative integral-operator compositions."
readme = "README.md"
requires-python = ">=3.10"
dependencies = [
    "sympy>=1.14,<2",
]

[project.optional-dependencies]
dev = [
    "pytest>=9,<10",
]

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-q"
```

## Estructura minima propuesta

```text
docs/
  MATHEMATICAL_CONVENTIONS.md
  MVP_SCOPE.md
src/
  symbolic_operator_calculus/
    __init__.py
    operators.py
    actions.py
    composition.py
    fourier.py
    kernels.py
    latex.py
    relations.py
tests/
  test_actions.py
  test_combined_kernel.py
  test_composition.py
  test_fourier.py
  test_latex.py
  test_noncommutative_expansion.py
  test_relations.py
```

Responsabilidades propuestas:

- `operators.py`: atomos de operador, producto ordenado y combinacion lineal.
- `actions.py`: reglas de accion de multiplicacion, shift transportado e
  integral formal.
- `kernels.py`: kernels formales \(L^+\), \(L^-\), \(R_{1,1}\), \(M_{1,2}\),
  \(M_{2,1}\) y \(C_{2,2}^{(1)}\).
- `composition.py`: generacion segura de variables mudas y composicion de
  kernels.
- `fourier.py`: convencion de Fourier y calculo simbolico de \(K^+\) y
  \(K^-\).
- `relations.py`: igualdad exacta y equivalencia modulo compactos como tipos
  separados.
- `latex.py`: renderizado LaTeX como proyeccion, no como representacion interna.

## Estructura de tests propuesta

Los tests deben comparar estructura matematica, no solo strings LaTeX.

- `test_noncommutative_expansion.py`: tests 1 a 4.
- `test_actions.py`: tests 5 a 7 y partes de 10 a 11.
- `test_composition.py`: tests 8 a 10.
- `test_combined_kernel.py`: tests 11 a 12.
- `test_relations.py`: tests 13 a 14.
- `test_fourier.py`: tests 15 a 16.
- `test_latex.py`: test 17.

## Plan de implementacion posterior

1. Crear `pyproject.toml`, `src/` y `tests/` sin logica matematica compleja.
2. Implementar atomos de operador y producto no conmutativo estructural.
3. Implementar expansion distributiva solo para sumas/productos del MVP.
4. Implementar reglas de accion para \(G_k\), \(\widetilde V_{\alpha_k}\),
   operadores integrales y \(R_{1,1}\).
5. Implementar generador de variables mudas con avoidance set explicito.
6. Implementar composicion de kernels y accion sobre \(f\).
7. Implementar \(M_{1,2}\), \(M_{2,1}\), \(C_{2,2}^{(1)}\) y su accion.
8. Implementar calculo Fourier con hipotesis \(d>0\) y documentacion del
   `Piecewise` que SymPy devuelve por defecto.
9. Implementar tipos de relacion exacta y modulo compactos.
10. Implementar renderizado LaTeX.
11. Ejecutar los 17 tests obligatorios.

## Riesgos matematicos

- Confundir \(A_{1,1}^{(-1)}\) con un inverso exacto. Mitigacion: tipo
  `FormalRegularizer` sin reglas automaticas de cancelacion.
- Usar equivalencias modulo compactos como igualdades exactas. Mitigacion:
  tipos de relacion separados y sin conversion implicita.
- Permitir reordenamiento de operadores al delegar demasiado en SymPy.
  Mitigacion: producto de operadores como tupla ordenada propia.
- Captura accidental de variables mudas. Mitigacion: generador de nombres con
  conjunto de variables prohibidas y tests especificos.
- Interpretar shifts sobre la variable equivocada. Mitigacion: reglas de accion
  locales y tests sobre aparicion de \(\gamma_k x\).

## Riesgos tecnicos

- SymPy puede normalizar factores conmutativos dentro de `Mul`. Esto es
  aceptable para kernels escalares, pero no para productos de operadores.
- SymPy puede ordenar terminos dentro de `Add`. Esto no debe definir el orden
  matematico esperado de la expansion \(+, -, -, +\).
- Algunas integrales de Fourier aparecen como `Piecewise` por condiciones de
  convergencia. El MVP debe conservar la salida original o aplicar una regla
  minima bajo hipotesis explicita \(d>0\).
- Comparar LaTeX directamente seria fragil. Los tests de LaTeX solo verificaran
  que las expresiones principales se pueden renderizar y contienen notacion
  esperada; la correccion matematica se verificara estructuralmente.
- Python 3.10.12 es suficiente para el MVP, pero conviene evitar sintaxis que
  requiera 3.11 o 3.12 si el proyecto fija `>=3.10`.

## Criterios de aceptacion

- La expresion principal se expande en cuatro terminos.
- Los signos son exactamente \(+, -, -, +\).
- El orden de factores coincide con el orden matematico especificado.
- La accion de cada tipo de operador coincide con su regla declarada.
- Las composiciones generan variables mudas frescas.
- \(M_{2,1} R_{1,1} M_{1,2}\) aparece en ese orden en el kernel combinado.
- \(R_{1,1}\) permanece formal.
- Igualdad exacta y equivalencia modulo compactos son relaciones distintas.
- Las transformadas inversas producen \(K^+\) y \(K^-\) bajo \(d>0\).
- Todas las expresiones principales exportan a LaTeX.
