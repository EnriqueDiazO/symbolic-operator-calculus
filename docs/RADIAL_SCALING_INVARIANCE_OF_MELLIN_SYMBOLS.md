# Radial-scaling invariance of Mellin symbols

## Statement

Let \(\gamma>0\) and

\[
 a_\gamma(x,\lambda):=a(x/\gamma,\lambda).
\]

If
\(a\in\widetilde{\mathcal E}(\mathbb R_+,V(\mathbb R))\), then

\[
 a_\gamma\in
 \widetilde{\mathcal E}(\mathbb R_+,V(\mathbb R)).
\]

This is a direct verification of every contract in the definition used by the
paper.

## Boundedness and continuity into \(V\)

The map \(x\mapsto x/\gamma\) is a homeomorphism of \(\mathbb R_+\). Hence
\(x\mapsto a(x/\gamma,\cdot)\) is continuous into \(V(\mathbb R)\), and

\[
 \sup_{x>0}\|a_\gamma(x,\cdot)\|_V
 =\sup_{z>0}\|a(z,\cdot)\|_V.
\]

Thus the `C_b(R+,V(R))` norm is unchanged.

## Slow radial oscillation

For the paper's multiplicative oscillation seminorm,

\[
\begin{aligned}
 \operatorname{cm}^{\infty}_t(a_\gamma)
 &=\sup_{x_1,x_2\in[t,2t]}
 \|a(x_1/\gamma,\cdot)-a(x_2/\gamma,\cdot)\|_\infty\\
 &=\operatorname{cm}^{\infty}_{t/\gamma}(a).
\end{aligned}
\]

Since \(t/\gamma\to0\) when \(t\to0\), and \(t/\gamma\to\infty\) when
\(t\to\infty\), both slowly oscillating endpoint conditions are preserved.

## Covariable translations

Translation acts only on \(\lambda\), and therefore

\[
 (a_\gamma)^h(x,\lambda)=a^h(x/\gamma,\lambda).
\]

The substitution \(z=x/\gamma\) gives exactly

\[
 \sup_{x>0}\|a_\gamma(x,\cdot)-(a_\gamma)^h(x,\cdot)\|_V
 =\sup_{z>0}\|a(z,\cdot)-a^h(z,\cdot)\|_V.
\]

Hence the defining uniform translation limit remains zero.

## Uniform derivative tails

Since

\[
 \partial_\lambda a_\gamma(x,\lambda)
 =\partial_\lambda a(x/\gamma,\lambda),
\]

one also has

\[
\begin{aligned}
 &\sup_{x>0}\int_{\mathbb R\setminus[-M,M]}
 |\partial_\lambda a_\gamma(x,\lambda)|\,d\lambda\\
 &\qquad=
 \sup_{z>0}\int_{\mathbb R\setminus[-M,M]}
 |\partial_\lambda a(z,\lambda)|\,d\lambda.
\end{aligned}
\]

The strengthened `E-tilde` tail condition is therefore unchanged.

## Endpoint limits and fibers

Any ordinary radial endpoint limit is preserved because scaling by a fixed
positive constant leaves sequences approaching zero or infinity at the same
endpoint. More generally, radial scaling defines an automorphism of the
slowly oscillating coefficient algebra. Its transpose maps each maximal-ideal
fiber over \(0\) or \(\infty\) into the same fiber, because the automorphism
fixes the endpoint values of continuous functions on the compactified
half-line. Thus scaling does not create or exchange endpoint fibers.

## Status and scope

The result is `ANALYTICALLY_PROVED`. It proves standard membership of
\(a_\gamma\), but not membership of
\(a_\gamma d_{\gamma^{-1}}\) when \(\gamma\ne1\). The latter remains an
ordered factorized representation.
