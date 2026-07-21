# Exact right-dilation composition in Mellin quantization

## Result

With the paper's normalized dilation and Mellin convention, the correct option
is **A**:

\[
 \operatorname{Op}(a)D_\gamma
 =\operatorname{Op}_{\mathrm{right}\text{-}\gamma}(a,d_\gamma),
 \qquad
 d_\gamma(\lambda)=\gamma^{i\lambda}.
\]

On the dense test class where the oscillatory integral is read directly, the
right-hand side is exactly the Kohn--Nirenberg quantization with amplitude
\(a(x,\lambda)d_\gamma(\lambda)\). The factorized notation is retained because
that product need not belong to the standard Mellin symbol class.

## 1. Spaces and normalization

Let

\[
 X=L^p(\mathbb R_+,d\mu),
 \qquad d\mu(x)=\frac{dx}{x},
 \qquad 1<p<\infty,
\]

and let \(X_\delta=L^p(\mathbb R_+,w_\delta)\), where

\[
 \kappa_\delta=\delta+\frac1p,
 \qquad
 (\Phi_\delta f)(x)=x^{\kappa_\delta}f(x).
\]

The paper's normalized dilation is

\[
 U_\gamma=\gamma^{\kappa_\delta}V_\gamma,
 \qquad
 (V_\gamma f)(x)=f(\gamma x),
 \qquad \gamma>0.
\]

Conjugating it gives, exactly,

\[
\begin{aligned}
 (\Phi_\delta U_\gamma\Phi_\delta^{-1}g)(x)
 &=x^{\kappa_\delta}\gamma^{\kappa_\delta}
   (\gamma x)^{-\kappa_\delta}g(\gamma x)\\
 &=g(\gamma x).
\end{aligned}
\]

Write \(D_\gamma g(x):=g(\gamma x)\). There is no residual weight factor.
In contrast, the unnormalized \(V_\gamma\) would produce
\(\gamma^{-\kappa_\delta}D_\gamma\).

## 2. Mellin transform and sign

The conventions are

\[
 (\mathcal Mg)(\lambda)
 =\int_0^\infty g(x)x^{-i\lambda}\frac{dx}{x},
 \qquad
 (\mathcal M^{-1}h)(x)
 =\frac1{2\pi}\int_{\mathbb R}h(\lambda)x^{i\lambda}\,d\lambda.
\]

For \(g\) in a Mellin test class and \(z=\gamma x\),

\[
\begin{aligned}
 \mathcal M(D_\gamma g)(\lambda)
 &=\int_0^\infty g(\gamma x)x^{-i\lambda}\frac{dx}{x}\\
 &=\int_0^\infty g(z)(z/\gamma)^{-i\lambda}\frac{dz}{z}\\
 &=\gamma^{i\lambda}(\mathcal Mg)(\lambda).
\end{aligned}
\]

Thus

\[
 d_\gamma(\lambda)=\gamma^{+i\lambda}.
\]

The plus sign is also checked using
\(g_0(x)=e^{-(\log x)^2}\). Since

\[
 (\mathcal Mg_0)(\lambda)=\sqrt\pi e^{-\lambda^2/4},
\]

translation by \(\log\gamma\) in logarithmic coordinates gives

\[
 \mathcal M[g_0(\gamma\,\cdot)](\lambda)
 =e^{+i\lambda\log\gamma}\sqrt\pi e^{-\lambda^2/4}.
\]

This rules out \(\gamma^{-i\lambda}\).

## 3. Quantization and ordered composition

The paper uses the left/Kohn--Nirenberg Mellin quantization

\[
 [\operatorname{Op}(a)g](x)
 =\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
 a(x,\lambda)\left(\frac{x}{y}\right)^{i\lambda}
 g(y)\frac{dy}{y}\,d\lambda.
\]

Take initially \(g\in C_c^\infty(\mathbb R_+)\), and interpret the outer
integral through symmetric frequency truncations. For each finite truncation,
the inner integral is absolutely integrable. Therefore the following uses only
the substitution \(z=\gamma y\) inside that integral; it does not exchange the
two integrations:

\[
\begin{aligned}
 [\operatorname{Op}(a)D_\gamma g](x)
 &=\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
 a(x,\lambda)\left(\frac{x}{y}\right)^{i\lambda}
 g(\gamma y)\frac{dy}{y}\,d\lambda\\
 &=\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
 a(x,\lambda)
 \left(\frac{\gamma x}{z}\right)^{i\lambda}
 g(z)\frac{dz}{z}\,d\lambda\\
 &=\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
 a(x,\lambda)\gamma^{i\lambda}
 \left(\frac{x}{z}\right)^{i\lambda}
 g(z)\frac{dz}{z}\,d\lambda.
\end{aligned}
\]

Hence the order is

\[
 (a,d_\gamma),
\]

not \((d_\gamma,a)\), and no commutation of operators has occurred.

## 4. Extension by continuity

Suppose \(\operatorname{Op}(a)\in\mathcal B(X)\). Since \(D_\gamma\) is an
invertible isometry on \(X\),

\[
 \|\operatorname{Op}(a)D_\gamma\|_{\mathcal B(X)}
 \leq \|\operatorname{Op}(a)\|_{\mathcal B(X)}.
\]

The identity on the dense test class therefore has a unique bounded extension.
We denote that extension by

\[
 \operatorname{Op}_{\mathrm{right}\text{-}\gamma}(a,d_\gamma).
\]

This definition does not use a standard-symbol boundedness theorem for
\(a d_\gamma\). It records the exact factorized operator already known to be
bounded from its two ordered factors.

Returning to the weighted space gives

\[
 R_1U_1
 =\Phi_\delta^{-1}
   \operatorname{Op}_{\mathrm{right}\text{-}\gamma_1}
   (r_1,d_{\gamma_1})
   \Phi_\delta.
\]

The AST continues to store the left-hand side as the ordered product
`Product((R1, U1))`; \(R_1\) is not expanded.

## 5. Symbol-class boundary

If \(\gamma\ne1\), then \(d_\gamma\notin V(\mathbb R)\). Consequently, the
exact operator identity does not prove

\[
 a(x,\lambda)d_\gamma(\lambda)
 \in\widetilde{\mathcal E}(\mathbb R_+,V(\mathbb R)).
\]

The Phase Q representation deliberately preserves the ordered pair
\((a,d_\gamma)\). No Fredholm, algebra-closure, or Wiener--Hopf conclusion is
drawn.
