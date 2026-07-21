# Exact left-dilation covariance in Mellin quantization

## Result

For \(\gamma>0\), define

\[
 a_\gamma(x,\lambda):=a(x/\gamma,\lambda),
 \qquad
 d_{\gamma^{-1}}(\lambda)=\gamma^{-i\lambda}.
\]

With the paper's Mellin convention, both candidate formulas are exact:

\[
 \boxed{
 D_\gamma^{-1}\operatorname{Op}(a)
 =\operatorname{Op}(a_\gamma)D_{\gamma^{-1}}
 =\operatorname{Op}_{\mathrm{right}\text{-}\gamma^{-1}}
   (a_\gamma,d_{\gamma^{-1}})
 }
\]

and

\[
 \boxed{
 D_\gamma^{-1}\operatorname{Op}(a)D_\gamma
 =\operatorname{Op}(a_\gamma).
 }
\]

No scalar normalization remains after the weighted conjugation. The first
formula is retained as an ordered factorized symbol when \(\gamma\ne1\); it
does not assert that \(a_\gamma d_{\gamma^{-1}}\) belongs to the standard
Mellin class.

## Definitions and dense test class

On \(X=L^p(\mathbb R_+,d\mu)\), \(d\mu(x)=dx/x\),

\[
 [\operatorname{Op}(a)f](x)
 =\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
 a(x,\lambda)\left(\frac{x}{y}\right)^{i\lambda}
 f(y)\frac{dy}{y}\,d\lambda,
\]

and \(D_\gamma f(x)=f(\gamma x)\). Initially take
\(f\in C_c^\infty(\mathbb R_+)\) and truncate the outer frequency integral
symmetrically. The inner integral is then absolutely integrable, and the
following calculation does not exchange integrations.

## Direct derivation of the left action

Evaluate the Mellin PDO at \(x/\gamma\):

\[
\begin{aligned}
 [D_\gamma^{-1}\operatorname{Op}(a)f](x)
 &= [\operatorname{Op}(a)f](x/\gamma)\\
 &=\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
 a(x/\gamma,\lambda)
 \left(\frac{x/\gamma}{y}\right)^{i\lambda}
 f(y)\frac{dy}{y}\,d\lambda\\
 &=\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
 a_\gamma(x,\lambda)\gamma^{-i\lambda}
 \left(\frac{x}{y}\right)^{i\lambda}
 f(y)\frac{dy}{y}\,d\lambda.
\end{aligned}
\]

Thus the radial argument is \(x/\gamma\), the oscillatory factor is on the
right, and its exponent is \(-i\lambda\). By the exact Phase Q right-action
formula applied to \(D_{\gamma^{-1}}\), the last amplitude is precisely the
ordered operator \(\operatorname{Op}(a_\gamma)D_{\gamma^{-1}}\).

No operator was commuted.

## Exact conjugation

Right multiplication by \(D_\gamma\) now gives

\[
\begin{aligned}
 D_\gamma^{-1}\operatorname{Op}(a)D_\gamma
 &=\operatorname{Op}(a_\gamma)
   D_{\gamma^{-1}}D_\gamma\\
 &=\operatorname{Op}(a_\gamma).
\end{aligned}
\]

Equivalently, the frequency factors cancel exactly:

\[
 d_{\gamma^{-1}}(\lambda)d_\gamma(\lambda)
 =\gamma^{-i\lambda}\gamma^{i\lambda}=1.
\]

## Independent sign and radial checks

Let \(g_0(x)=e^{-(\log x)^2}\). In logarithmic coordinates,
\(D_{\gamma^{-1}}g_0\) is the translate
\(e^{-(u-\log\gamma)^2}\). Its Fourier/Mellin transform is

\[
 \gamma^{-i\lambda}\sqrt\pi e^{-\lambda^2/4},
\]

which confirms the negative sign. For a spatial multiplication symbol
\(a=a(x)\), direct evaluation gives

\[
 D_\gamma^{-1}M_aD_\gamma=M_{a(\cdot/\gamma)},
\]

confirming \(x/\gamma\), rather than \(\gamma x\).

## Bounded extension

If \(a\in\widetilde{\mathcal E}\), radial scaling preserves that class and
\(\operatorname{Op}(a_\gamma)\) is bounded. Since both dilations are
invertible isometries, the identities on the dense test class have unique
bounded extensions to \(X\).

On the weighted space,

\[
 \Phi_\delta U_\gamma^{\pm1}\Phi_\delta^{-1}
 =D_{\gamma^{\pm1}},
\]

so conjugating the two identities by \(\Phi_\delta^{-1}\) introduces no
additional scalar. These are exact operator identities, independently of the
later modulo-compact Cauchy semiproducts.
