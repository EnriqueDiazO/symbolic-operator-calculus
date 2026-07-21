# Auditoría de fuentes: cancelación de dilataciones en la corrección de Schur

Fecha de auditoría: 2026-07-21.

## 1. Alcance y decisión

Esta auditoría contrasta la forma deseada

\[
A_{2,1}\simeq W(b_{2,1})V_\gamma,
\qquad
A_{1,2}\simeq V_{\gamma^{-1}}W(b_{1,2}),
\qquad
R_{1,1}\simeq \operatorname{Op}_M(r_{1,1})
\]

con las fórmulas que aparecen realmente en el artículo auxiliar de solo
lectura

```text
/home/enriquedo/PersonalProjects/Papers/DiazOcampoHasemanBVP2026
```

Se inspeccionaron prioritariamente:

```text
sections/03_transport_and_kernels.tex
sections/04_target_theorem53.tex
sections/05_haseman_operator_with_shift.tex
sections/06_shur_reduction.tex
NotationMacros2.tex
```

La decisión de esta auditoría es:

```text
BLOCKED_WITH_PRECISE_OBLIGATION
```

Las convenciones integrales del artículo permiten derivar una cancelación de
dilataciones para un sándwich abstracto una vez que existen factores
adyacentes `V_gamma Op_M(r) V_gamma_inverse`. Sin embargo, el artículo no
demuestra que los bloques reales del complemento de Schur admitan las dos
factorizaciones externas que producirían esos factores adyacentes. En
particular:

1. las identidades `W(b) V_gamma` y `V_gamma_inverse W(b)` se demuestran para
   unos núcleos auxiliares `mathcal W_{k,j}`, no para `A_{k,j}`;
2. los bloques `A_{k,j}` contienen coeficientes, factores de transporte y dos
   localizadores `chi_infinity`;
3. el regularizador que entra en Schur es
   `A_{1,1}^{(-1)} = T_{1,-} R_1 Z_1^{-1}`, no un Mellin PDO desnudo;
4. el signo de `C_{2,2}^{(1)}` en el artículo no coincide con el uso de ese
   nombre para el producto crudo `A_{2,1} R_{1,1} A_{1,2}`.

No se declara resuelto el pivote ni se afirma pertenencia a un álgebra mixta de
Fredholm.

## 2. Estados semánticos usados

Esta auditoría separa estrictamente los siguientes estados:

| Estado | Significado en este documento |
|---|---|
| `EXACT_IDENTITY` | Igualdad algebraica u operatorial que el artículo presenta como exacta, o que se deriva directamente de sus definiciones integrales sin reordenar factores. |
| `EQUALITY_MODULO_COMPACTS` | Igualdad de clases en el cociente por compactos; no se promueve a igualdad exacta. |
| `CONDITIONAL_ON_REGULARIZER_REPRESENTATION` | Paso válido sólo después de fijar qué operador denota `R_{1,1}` y suministrar la representación Mellin compatible con ese operador y espacio. |
| `UNPROVED_MIXED_ALGEBRA_MEMBERSHIP` | No hay prueba citada en el artículo de que el producto final Wiener--Hopf--Mellin--Wiener--Hopf pertenezca a la clase requerida. |

Las afirmaciones de esta auditoría describen lo escrito en el artículo. No
constituyen por sí mismas una verificación independiente de los teoremas
externos que el artículo cita.

## 3. Tabla de extracción y acciones

| Concepto | Artículo | Estado | Acción requerida en el repositorio principal |
|---|---|---|---|
| `Upsilon_k` | `sections/03_transport_and_kernels.tex:268-337`: `(Upsilon_k f)(x)=(x-ic_k)^{-1}f((x-ic_k)^{-1})`, con inversa explícita e isometría. | `EXACT_IDENTITY` | Registrar la fórmula, dominio, codominio y convención de peso sin cambiar el exponente. |
| Espacio transportado | `sections/03_transport_and_kernels.tex:339-429`: `sigma: L^p(Gamma,w_tilde) -> L_m^p(R_+,w_delta)` es isométrico. | `EXACT_IDENTITY` | Usar el producto de `m` copias de la misma realización ponderada. |
| Peso | `sections/03_transport_and_kernels.tex:61-130,241-262,1084-1099`: `w_delta(x)=x^delta` y la norma es `int |x^delta f(x)|^p dx`. | `EXACT_IDENTITY` | Modelar el peso multiplicativo, equivalente a densidad `x^(p delta)`, no a densidad `x^delta`. |
| Normalización | `sections/03_transport_and_kernels.tex:1064-1108,1368-1469`: `kappa_delta=delta+1/p`, `U_gamma=gamma^kappa_delta V_gamma`. | `EXACT_IDENTITY` | Distinguir `V_gamma` de `U_gamma`; no borrar el escalar de normalización. |
| Convención de dilatación | `sections/05_haseman_operator_with_shift.tex:86-121`: `(V_{alpha_k}f)(x)=f(gamma_k x)`. | `EXACT_IDENTITY` | Fijar la convención pullback antes de toda regla de covariancia. |
| Dilatación relativa | `sections/05_haseman_operator_with_shift.tex:1876-1892`: `beta_{k,j}=alpha_j^{-1} o alpha_k`, `beta_{k,j}(x)=(gamma_k/gamma_j)x`. | `EXACT_IDENTITY` | Si se usa `gamma_{k,j}`, documentarlo como alias derivado `gamma_k/gamma_j`; el artículo no usa ese nombre. |
| Núcleo WH conjugado | `sections/05_haseman_operator_with_shift.tex:1790-2166`: `mathcal W_{k,j}=V_k W(b_{k,j})V_j^{-1}=V_beta W(b^(j))=W(b^(k))V_beta`. | `EXACT_IDENTITY` | Reutilizar sólo para `mathcal W_{k,j}` o demostrar primero que el factor real de Schur se reduce a ese núcleo. |
| `chi_infinity` | `sections/03_transport_and_kernels.tex:747-762` y `sections/04_target_theorem53.tex:228-237,304-338`: localización simétrica de salida y entrada. | `EQUALITY_MODULO_COMPACTS` para sustituir `T_{k,j}` | Conservar ambos cortes y abrir una obligación explícita para su covariancia bajo dilatación. |
| Bloque exacto `A_{k,j}` | `sections/05_haseman_operator_with_shift.tex:301-320`: `A_{k,j}=1/2 (Vtilde_k-G_k)T_{k,j}` para `k != j`. | `EXACT_IDENTITY` | Mantener el factor de fila a la izquierda y `T_{k,j}` a la derecha. |
| Modelo WH de `A_{k,j}` | `sections/05_haseman_operator_with_shift.tex:1668-1726`; caso `m=2` en `sections/06_shur_reduction.tex:44-95`. | `EQUALITY_MODULO_COMPACTS` | No promoverlo a igualdad exacta ni reordenar el factor `(Vtilde_k-G_k)`. |
| Símbolo `r_{1,1}` | `sections/06_shur_reduction.tex:134-160`: existe bajo las hipótesis de Fredholm y regulariza `Op(n_1)` módulo compactos. | `CONDITIONAL_ON_REGULARIZER_REPRESENTATION` | Mantenerlo existencial; no inventar kernel ni fórmula cerrada. |
| Operador `mathcal R_1` | `sections/06_shur_reduction.tex:161-171`: `R_1=Phi_delta^{-1} Op(r_{1,1}) Phi_delta`. | `EXACT_IDENTITY` una vez elegido `r_{1,1}` | Representarlo como Mellin PDO transportado al espacio ponderado. |
| Regularizador de Schur | `sections/06_shur_reduction.tex:173-195`: `A_{1,1}^{(-1)}=T_{1,-} R_1 Z_1^{-1}`. | `EXACT_IDENTITY` de la representación elegida; propiedades inversas sólo módulo compactos | No identificar `A_{1,1}^{(-1)}` con `Op(r_{1,1})`. |
| Corrección del artículo | `sections/06_shur_reduction.tex:299-348`: el pivote resta `A21 A11^(-1) A12`, pero `C_{2,2}^{(1)}` se define con el factor positivo resultante del signo negativo de `A21`. | `EQUALITY_MODULO_COMPACTS` respecto al producto de los bloques reales | Conservar por separado `P_raw=A21 A11^(-1) A12` y `C_article`, con `C_article ~= -P_raw`. |
| Representación Mellin del bloque localizado | `sections/05_haseman_operator_with_shift.tex:2282-2472`: símbolo candidato y verificación de clase pendiente. | No certificada | No usarla como teorema ni como evidencia de cierre. |
| Producto mixto final | No hay teorema en las secciones auditadas que cierre el sándwich final en una álgebra de Fredholm. | `UNPROVED_MIXED_ALGEBRA_MEMBERSHIP` | Emitir obligación y detener cualquier afirmación de Fredholm o índice. |

## 4. Transporte y peso reales

La rama se parametriza por

\[
\eta_k(x)=(x-ic_k)^{-1},
\qquad
|dt|=|x-ic_k|^{-2}\,dx.
\]

Fuente:
`/home/enriquedo/PersonalProjects/Papers/DiazOcampoHasemanBVP2026/sections/03_transport_and_kernels.tex:9-27`.

El peso de rama es

\[
w_k(t)=|1+ic_kt|^{1-2/p-\mu}|t|^\mu,
\]

y la convención ponderada es

\[
\|f\|_{L^p(\Gamma_k,w_k)}=\|w_kf\|_{L^p(\Gamma_k)}.
\]

Fuentes:
`sections/03_transport_and_kernels.tex:61-101` y
`sections/03_transport_and_kernels.tex:105-130` del artículo auxiliar.

Se define

\[
\delta=1-\frac2p-\mu,
\qquad
w_\delta(x)=x^\delta.
\]

Fuente: `sections/03_transport_and_kernels.tex:241-262`.

Por tanto, la norma real en la semirrecta es

\[
\|f\|_{L^p(\mathbb R_+,w_\delta)}^p
=
\int_0^\infty |x^\delta f(x)|^p\,dx.
\]

La igualdad aparece explícitamente en
`sections/03_transport_and_kernels.tex:1084-1099`. El motor no debe
interpretar `w_delta` como una densidad `x^delta dx`; la densidad equivalente
es `x^(p delta) dx`.

El transporte de rama es

\[
\Upsilon_k:L^p(\Gamma_k,w_k)\longrightarrow
L^p(\mathbb R_+,w_\delta),
\]

\[
(\Upsilon_kf)(x)
=
(x-ic_k)^{-1}f\!\left((x-ic_k)^{-1}\right),
\]

con inversa

\[
(\Upsilon_k^{-1}\varphi)(t)
=
t^{-1}\varphi(t^{-1}+ic_k).
\]

El artículo prueba la isometría y la sobreyectividad en
`sections/03_transport_and_kernels.tex:268-337`. Globalmente,

\[
\sigma f=(\Upsilon_1f_1,\ldots,\Upsilon_mf_m)
\]

es un isomorfismo isométrico hacia
`L_m^p(R_+,w_delta)`; véase
`sections/03_transport_and_kernels.tex:339-429`.

## 5. Dilatación ponderada y parámetro relativo

El artículo fija

\[
\alpha_k(x)=\gamma_kx,
\qquad
(V_{\alpha_k}f)(x)=f(\gamma_kx),
\qquad
\gamma_k>0.
\]

Fuente: `sections/05_haseman_operator_with_shift.tex:86-121`.

El shift que resulta del transporte desde la curva no es sólo `V_{alpha_k}`:

\[
\widetilde V_{\alpha_k}
=
\rho_kV_{\alpha_k},
\qquad
\rho_k(x)=\frac{\gamma_kx-ic_k}{x-ic_k}.
\]

Además, `rho_k(0)=1` y `rho_k(infinity)=gamma_k`. Fuente:
`sections/05_haseman_operator_with_shift.tex:123-160`.

Con

\[
\kappa_\delta=\delta+\frac1p,
\]

la norma real da

\[
\|V_\gamma f\|^p
=
\gamma^{-p\delta-1}\|f\|^p
=
\gamma^{-p\kappa_\delta}\|f\|^p.
\]

Por ello, la dilatación isométrica es

\[
U_\gamma=\gamma^{\kappa_\delta}V_\gamma.
\]

El artículo prueba esta normalización en
`sections/03_transport_and_kernels.tex:1368-1469`, después de definir
`kappa_delta` en `sections/03_transport_and_kernels.tex:1064-1108`.

No aparece un parámetro llamado `gamma_{k,j}`. La notación real es

\[
\beta_{k,j}:=\alpha_j^{-1}\circ\alpha_k,
\qquad
\beta_{k,j}(x)=\frac{\gamma_k}{\gamma_j}x,
\]

y el orden se verifica mediante

\[
V_{\alpha_k}V_{\alpha_j}^{-1}
=
V_{\alpha_j^{-1}\circ\alpha_k}.
\]

Fuente: `sections/05_haseman_operator_with_shift.tex:1876-1892`.

Si la API usa `gamma_{k,j}`, su definición respaldada como alias derivado es

\[
\gamma_{k,j}:=\frac{\gamma_k}{\gamma_j}.
\]

Para `m=2`, si `gamma=gamma_2/gamma_1`, entonces
`gamma_{2,1}=gamma` y `gamma_{1,2}=gamma^{-1}`.

## 6. Lo que sí factoriza exactamente como Wiener--Hopf y dilatación

La convención de Fourier es

\[
(\mathcal Ff)(\lambda)
=
\int_{\mathbb R}f(t)e^{-it\lambda}\,dt,
\]

y `W(b)=chi_+ F^{-1} b F chi_+ I`; véase
`sections/04_target_theorem53.tex:17-48`.

La regla de escala que usa el artículo es

\[
\mathcal F[a h(a\,\cdot)](\lambda)
=
(\mathcal Fh)(\lambda/a),
\qquad a>0.
\]

Fuente: `sections/05_haseman_operator_with_shift.tex:1894-1900`.

Se define el núcleo auxiliar no localizado

\[
\mathcal W_{k,j}
:=
V_{\alpha_k}W(b_{k,j})V_{\alpha_j}^{-1}
\]

en `sections/05_haseman_operator_with_shift.tex:1790-1806`. Para este operador,
el artículo demuestra la identidad exacta

\[
\mathcal W_{k,j}
=
V_{\beta_{k,j}}W\!\left(b_{k,j}^{(j)}\right)
=
W\!\left(b_{k,j}^{(k)}\right)V_{\beta_{k,j}},
\]

con

\[
b_{k,j}^{(j)}(\lambda)
=
b_{k,j}(\lambda/\gamma_j),
\qquad
b_{k,j}^{(k)}(\lambda)
=
b_{k,j}(\lambda/\gamma_k).
\]

Fuentes:

- `sections/05_haseman_operator_with_shift.tex:1902-2023` para la colocación
  izquierda de la dilatación;
- `sections/05_haseman_operator_with_shift.tex:2025-2132` para la colocación
  derecha;
- `sections/05_haseman_operator_with_shift.tex:2134-2166` para la identidad
  conjunta y la correspondencia de símbolos.

En particular, con `gamma=gamma_2/gamma_1`, se obtiene exactamente

\[
\mathcal W_{2,1}
=
W\!\left(b_{2,1}^{(2)}\right)V_\gamma,
\]

\[
\mathcal W_{1,2}
=
V_{\gamma^{-1}}W\!\left(b_{1,2}^{(2)}\right),
\]

donde ambos símbolos que rodearían una cancelación central usan la escala de
la rama 2:

\[
b_{2,1}^{(2)}(\lambda)=b_{2,1}(\lambda/\gamma_2),
\qquad
b_{1,2}^{(2)}(\lambda)=b_{1,2}(\lambda/\gamma_2).
\]

Estas fórmulas tienen estado `EXACT_IDENTITY`, pero su sujeto es
`mathcal W_{k,j}`, no `A_{k,j}`.

## 7. Localización: posición de `chi_infinity` y límite de la regla exacta

El corte satisface

\[
0\leq\chi_\infty\leq1,
\qquad
\chi_\infty=0\text{ en }[0,1],
\qquad
\chi_\infty=1\text{ en }[2,\infty).
\]

Como `x=infinity` corresponde a la cúspide, localiza el modelo cerca de ese
punto. Fuente: `sections/03_transport_and_kernels.tex:747-762`.

Para `k != j`, el artículo prueba

\[
T_{k,j}
\simeq
\chi_\infty W(b_{k,j})\chi_\infty I.
\]

El primer corte localiza la salida y el segundo la entrada. Fuentes:

- `sections/04_target_theorem53.tex:228-237`;
- `sections/04_target_theorem53.tex:304-338`.

El artículo sólo obtiene compactitud del conmutador entre el multiplicador por
`chi_infinity` y el Wiener--Hopf correspondiente en
`sections/04_target_theorem53.tex:318-338`. Eso no es una regla general para
cruzar el corte a través de una dilatación.

Al conjugar el bloque localizado aparecen explícitamente

\[
\chi_\infty(\gamma_kx)
\quad\text{y}\quad
\chi_\infty(\gamma_jy),
\]

además de `rho_k(x)` y `rho_j(y)^{-1}`. Véase
`sections/05_haseman_operator_with_shift.tex:2172-2278`. Por tanto, la
factorización exacta del núcleo no localizado no puede trasplantarse al bloque
localizado sin una prueba adicional de compactitud o una identidad que retenga
los cortes reescalados.

## 8. Forma y orden reales de `A12` y `A21`

Antes del modelo Wiener--Hopf, el artículo da la identidad exacta

\[
\mathcal A_{k,j}
=
\frac12
(\widetilde V_{\alpha_k}-G_kI)T_{k,j},
\qquad k\ne j.
\]

Fuente: `sections/05_haseman_operator_with_shift.tex:301-320`.

Al sustituir el modelo localizado sólo se obtiene

\[
\mathcal A_{k,j}
\simeq
\frac12
(\widetilde V_{\alpha_k}-G_kI)
\chi_\infty W(b_{k,j})\chi_\infty I.
\]

Fuente: `sections/05_haseman_operator_with_shift.tex:1668-1726`. El factor de
la fila `k` permanece a la izquierda; no se conmuta con Wiener--Hopf.

Para dos ramas, `sections/06_shur_reduction.tex:44-95` fija

\[
A_{1,2}
\simeq
(\widetilde V_{\alpha_1}-G_1)W_{1,2}^{+},
\]

\[
A_{2,1}
\simeq
-(\widetilde V_{\alpha_2}-G_2)W_{2,1}^{-},
\]

donde

\[
W_{1,2}^{+}
=
\chi_\infty
W\!\left(\chi_+(\cdot)e^{(c_1-c_2)(\cdot)}\right)
\chi_\infty I,
\]

\[
W_{2,1}^{-}
=
\chi_\infty
W\!\left(\chi_-(\cdot)e^{(c_2-c_1)(\cdot)}\right)
\chi_\infty I.
\]

El orden de actuación, de derecha a izquierda, es:

1. corte de entrada;
2. Wiener--Hopf;
3. corte de salida;
4. factor `(Vtilde_k-G_k)`.

No hay en estas fórmulas una dilatación relativa aislada a la derecha de
`A21` ni a la izquierda de `A12`. Identificar estos bloques con
`mathcal W_{2,1}` y `mathcal W_{1,2}` sería una inferencia no demostrada.

## 9. Regularizador real y estatus Mellin

La notación del artículo distingue:

1. `r_{1,1}`, símbolo Mellin regularizador;
2. `mathcal R_1`, regularizador del producto diagonal normalizado `mathcal N_1`;
3. `A_{1,1}^{(-1)}`, regularizador del bloque diagonal que entra en Schur.

Los macros relevantes están en `NotationMacros2.tex:147-152,216-219`. Aunque
existe un macro para `mathcal R_{1,1}`, las secciones auditadas usan
`mathcal R_1` y `A_{1,1}^{(-1)}`.

Primero,

\[
\mathcal N_1
=
\widehat{\mathcal A}_{1,1}\mathcal T_{1,-},
\qquad
\Phi_\delta\mathcal N_1\Phi_\delta^{-1}
=
\operatorname{Op}(n_1).
\]

La segunda igualdad se declara exacta en
`sections/05_haseman_operator_with_shift.tex:792-900`; su reutilización en la
reducción está en `sections/06_shur_reduction.tex:109-132`.

Bajo las hipótesis de Fredholm, el artículo invoca un teorema de regularización
y obtiene existencialmente

\[
r_{1,1}\in\widetilde{\mathcal E}
\]

tal que

\[
\operatorname{Op}(n_1)\operatorname{Op}(r_{1,1})\simeq I,
\qquad
\operatorname{Op}(r_{1,1})\operatorname{Op}(n_1)\simeq I.
\]

Fuente: `sections/06_shur_reduction.tex:134-160`. No se proporciona una fórmula
cerrada ni un kernel para `r_{1,1}`.

Se define exactamente

\[
\mathcal R_1
:=
\Phi_\delta^{-1}\operatorname{Op}(r_{1,1})\Phi_\delta.
\]

Fuente: `sections/06_shur_reduction.tex:161-171`.

Pero el regularizador usado por el complemento de Schur es

\[
\boxed{
A_{1,1}^{(-1)}
=
\mathcal T_{1,-}\mathcal R_1Z_1^{-1}
=
\mathcal T_{1,-}
\Phi_\delta^{-1}\operatorname{Op}(r_{1,1})
\Phi_\delta Z_1^{-1}.
}
\]

Fuente: `sections/06_shur_reduction.tex:173-195`. Sus relaciones de
regularización son módulo compactos, como se verifica en
`sections/06_shur_reduction.tex:196-260`. El superíndice `(-1)` no significa
inversa exacta; la misma advertencia aparece en
`sections/05_haseman_operator_with_shift.tex:1336-1398`.

Consecuencias semánticas:

- si `R11` del código significa `mathcal R_1`, su representación Mellin
  ponderada está fijada exactamente una vez elegido el símbolo existencial;
- si `R11` significa el regularizador del bloque `A11` que aparece en Schur,
  entonces `R11=A_{1,1}^{(-1)}` y no es `Op_M(r_{1,1})` desnudo;
- ninguna de las dos interpretaciones autoriza inventar un kernel ordinario;
- la fórmula `R11 ~= Op_M(r11)` es ambigua e insuficiente mientras no declare
  el espacio, las conjugaciones por `Phi_delta` y los factores
  `T_{1,-}`, `Z_1^{-1}`.

## 10. Conflicto de signo en `C_{2,2}^{(1)}`

El artículo define el pivote

\[
A_{2,2}^{(1)}
:=
A_{2,2}-A_{2,1}A_{1,1}^{(-1)}A_{1,2}.
\]

Fuente: `sections/06_shur_reduction.tex:281-305`.

Sea

\[
P_{\mathrm{raw}}
:=
A_{2,1}A_{1,1}^{(-1)}A_{1,2}.
\]

Como el modelo normalizado de `A21` lleva un signo negativo y el de `A12` no,
`sections/06_shur_reduction.tex:310-348` da

\[
P_{\mathrm{raw}}
\simeq
-Q,
\]

donde

\[
Q
=
(\widetilde V_{\alpha_2}-G_2)
W_{2,1}^{-}
A_{1,1}^{(-1)}
(\widetilde V_{\alpha_1}-G_1)
W_{1,2}^{+}.
\]

El artículo llama precisamente `C_{2,2}^{(1)}` a este factor positivo `Q`.
Por tanto,

\[
\boxed{
C_{2,2,\mathrm{artículo}}^{(1)}
\simeq
-P_{\mathrm{raw}}.
}
\]

Si una API llama `C_{2,2}^{(1)}` al producto crudo
`A21 R11 A12`, esa API y el artículo usan nombres iguales para objetos de signo
opuesto. Deben conservarse dos identificadores semánticos distintos. El pivote
queda

\[
A_{2,2}^{(1)}
\simeq
A_{2,2}+C_{2,2,\mathrm{artículo}}^{(1)}.
\]

## 11. Factorización real de la corrección

Al sustituir el regularizador elegido, el artículo conserva exactamente el
orden

\[
\begin{aligned}
C_{2,2,\mathrm{artículo}}^{(1)}
={}&
(\widetilde V_{\alpha_2}-G_2)
W_{2,1}^{-}
\mathcal T_{1,-}
\mathcal R_1
Z_1^{-1}
(\widetilde V_{\alpha_1}-G_1)
W_{1,2}^{+}.
\end{aligned}
\]

Fuente: `sections/06_shur_reduction.tex:359-380`.

Luego define

\[
\mathcal B_{2,1}^{-}
:=
(\widetilde V_{\alpha_2}-G_2)
W_{2,1}^{-}\mathcal T_{1,-},
\]

\[
\mathcal B_{1,2}^{+}
:=
Z_1^{-1}
(\widetilde V_{\alpha_1}-G_1)
W_{1,2}^{+},
\]

y obtiene

\[
C_{2,2,\mathrm{artículo}}^{(1)}
=
\mathcal B_{2,1}^{-}\mathcal R_1\mathcal B_{1,2}^{+}.
\]

Fuentes: `sections/06_shur_reduction.tex:382-432`.

Ésta es la interfaz correcta para una futura cancelación. La obligación no es
simplemente volver a factorizar `A21` y `A12`: debe tratar los factores reales
`mathcal B_{2,1}^{-}` y `mathcal B_{1,2}^{+}`, o demostrar una normalización
equivalente que explique explícitamente dónde quedan `T_{1,-}` y `Z_1^{-1}`.

## 12. Covariancia Mellin derivable de las convenciones

La definición del Mellin PDO es

\[
(\operatorname{Op}(a)f)(x)
=
\frac1{2\pi}
\int_{\mathbb R}\int_0^\infty
a(x,\lambda)
\left(\frac{x}{y}\right)^{i\lambda}
f(y)\frac{dy}{y}\,d\lambda.
\]

Fuente: `sections/03_transport_and_kernels.tex:950-965`.

Con `(V_gamma f)(x)=f(gamma x)`, una sustitución directa `y=gamma z` da,
inicialmente sobre funciones de prueba,

\[
\boxed{
V_\gamma\operatorname{Op}(a)V_{\gamma^{-1}}
=
\operatorname{Op}(a^\gamma),
\qquad
a^\gamma(x,\lambda)=a(\gamma x,\lambda).
}
\]

No aparece un Jacobiano adicional porque `dy/y=dz/z`. La misma fórmula vale
para `U_gamma` porque los escalares normalizadores de la conjugación se
cancelan. En la dirección inversa,

\[
V_{\gamma^{-1}}\operatorname{Op}(a)V_\gamma
=
\operatorname{Op}\!\left(a(x/\gamma,\lambda)\right).
\]

Esta regla tiene estado `EXACT_IDENTITY` cuando el operador está definido en
un dominio común y la identidad sobre funciones de prueba se extiende por
acotación. No crea por sí sola los factores adyacentes en la corrección real.

## 13. Representación Mellin pendiente de los bloques localizados

Para el bloque conjugado y localizado, el artículo obtiene una amplitud exacta
y define sólo un símbolo candidato; véase
`sections/05_haseman_operator_with_shift.tex:2282-2420`. La identificación con
`Op(a_{k,j})` se formula bajo la condición de que aplique inversión Mellin en
`sections/05_haseman_operator_with_shift.tex:2422-2452`.

El propio artículo declara en
`sections/05_haseman_operator_with_shift.tex:2454-2472` que todavía deben
probarse las estimaciones en una clase admisible de símbolos. Por tanto:

```text
full_localized_off_diagonal_mellin_representation: UNPROVED
```

No puede usarse ese pasaje para justificar cierre, acotación o composición en
el álgebra mixta deseada.

## 14. Obligación precisa que bloquea la cancelación

Para cambiar la decisión a una cancelación certificada hace falta demostrar,
en el mismo espacio `L^p(R_+,w_delta)` y con una única convención de dilatación,
un lema que cubra todos los puntos siguientes.

### 14.1 Factorizaciones externas del término real

Debe probarse, exactamente o módulo compactos con evidencia trazable, una
factorización de los factores que aparecen realmente en

\[
C_{2,2,\mathrm{artículo}}^{(1)}
=
\mathcal B_{2,1}^{-}\mathcal R_1\mathcal B_{1,2}^{+},
\]

de la forma

\[
\mathcal B_{2,1}^{-}
\simeq
W(\widehat b_{2,1})D_\gamma,
\qquad
\mathcal B_{1,2}^{+}
\simeq
D_{\gamma^{-1}}W(\widehat b_{1,2}),
\]

o una variante explícitamente equivalente. Aquí `D_gamma` debe declararse
como `V_gamma` o `U_gamma`; si se usa la versión normalizada, deben conservarse
todos los escalares `gamma^kappa_delta`.

La prueba debe explicar, sin conmutaciones implícitas:

1. cómo se tratan `(Vtilde_2-G_2)` y `(Vtilde_1-G_1)`;
2. cómo se absorben o conservan `T_{1,-}` y `Z_1^{-1}`;
3. qué símbolos reescalados aparecen, con la convención
   `b(lambda/gamma_k)` o `b(lambda/gamma_j)` correcta;
4. cómo se controlan `chi_infinity(gamma_k x)` y
   `chi_infinity(gamma_j y)`;
5. por qué cada residuo descartado es compacto en el espacio ponderado real.

Demostrar únicamente la identidad para `mathcal W_{k,j}` no satisface esta
obligación.

### 14.2 Representación central inequívoca

Debe fijarse si el operador central es:

```text
mathcal R_1 = Phi_delta^{-1} Op(r_1,1) Phi_delta
```

o

```text
A_1,1^(-1) = T_1,- mathcal R_1 Z_1^(-1).
```

Si se usa la primera forma, las factorizaciones externas deben ser las de
`mathcal B_{2,1}^{-}` y `mathcal B_{1,2}^{+}`. Si se usa la segunda, no puede
afirmarse `R11 ~= Op_M(r11)` sin transportar también los factores auxiliares.

### 14.3 Aplicación de covariancia y relación resultante

Sólo después de 14.1 y 14.2 puede aplicarse

\[
D_\gamma\mathcal R_1D_{\gamma^{-1}}
=
\Phi_\delta^{-1}
\operatorname{Op}\!\left(r_{1,1}(\gamma x,\lambda)\right)
\Phi_\delta
\]

en la normalización correspondiente. La relación de salida debe ser la más
débil de las relaciones de entrada. Si cualquiera de las factorizaciones
externas sólo vale módulo compactos, el resultado completo sólo puede tener
estado `EQUALITY_MODULO_COMPACTS`.

### 14.4 Cierre del producto mixto

Incluso después de cancelar las dilataciones, queda por probar que

\[
W(\widehat b_{2,1})
\operatorname{Op}_M(r_{1,1}^{\gamma})
W(\widehat b_{1,2})
\]

pertenece a la clase operatorial requerida y admite el cálculo de Fredholm que
se quiera aplicar. Hasta disponer de ese teorema:

```text
mixed_algebra_membership: UNPROVED_MIXED_ALGEBRA_MEMBERSHIP
```

## 15. Resultado semántico final de la auditoría

```text
transport_and_weight_conventions: EXACT_IDENTITY
relative_dilation_core_factorization: EXACT_IDENTITY
actual_A12_A21_wh_models: EQUALITY_MODULO_COMPACTS
actual_A12_A21_relative_dilation_factorization: UNPROVED
weighted_normalized_product_mellin_representation: EXACT_IDENTITY
chosen_mathcal_R1_mellin_representation: CONDITIONAL_ON_REGULARIZER_REPRESENTATION
A11_regularizer_as_bare_mellin_pdo: NOT_SUPPORTED
localized_off_diagonal_mellin_representation: UNPROVED
dilation_cancellation_after_adjacency_is_established: EXACT_IDENTITY
initial_block_representations_needed_for_adjacency: UNPROVED
mixed_algebra_membership: UNPROVED_MIXED_ALGEBRA_MEMBERSHIP
phase_decision: BLOCKED_WITH_PRECISE_OBLIGATION
```

La fase puede implementar tipos, hipótesis y obligaciones que representen este
resultado condicional. No debe emitir `CERTIFIED_CANCELLATION` para la
corrección de Schur real hasta satisfacer las obligaciones de la sección 14.
