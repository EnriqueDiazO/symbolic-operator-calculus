# Oscillatory dilation symbol: classification and factorized boundary

## Classification result

For \(\gamma>0\), set

\[
 d_\gamma(\lambda)=\gamma^{i\lambda}
 =e^{i\lambda\log\gamma}.
\]

There are two disjoint cases.

| Case | Result |
|---|---|
| \(\gamma=1\) | \(d_\gamma=1\), so its variation is zero and its membership in \(V(\mathbb R)\) is certified. |
| \(\gamma\ne1\) | \(d_\gamma\notin V(\mathbb R)\). This says nothing against boundedness of the associated dilation operator. |

The project uses the standard bounded-variation symbol class: its elements are
absolutely continuous functions whose total variation

\[
 \operatorname{Var}(v;\mathbb R)
 =\int_{\mathbb R}|v'(\lambda)|\,d\lambda
\]

is finite, with the corresponding supremum-plus-variation norm.

## Proof of nonmembership for a nontrivial dilation

Although \(|d_\gamma(\lambda)|=1\),

\[
 d_\gamma'(\lambda)
 =i\log(\gamma)d_\gamma(\lambda),
 \qquad
 |d_\gamma'(\lambda)|=|\log\gamma|.
\]

Consequently, for every \(T>0\),

\[
 \operatorname{Var}(d_\gamma;[-T,T])
 =2T|\log\gamma|.
\]

When \(\gamma\ne1\), this tends to infinity with \(T\). Hence the total
variation on \(\mathbb R\) is infinite.

The obstruction is also visible at both endpoints. With
\(c=\log\gamma\ne0\), choose sequences tending to either endpoint along which
\(c\lambda\) is an even multiple of \(\pi\), and sequences along which it is an
odd multiple. The corresponding values of \(d_\gamma\) are \(1\) and \(-1\).
Thus the oscillation persists and no limit exists at \(+\infty\) or
\(-\infty\).

This proof excludes only standard \(V(\mathbb R)\)-membership. The normalized
operator \(D_\gamma f(x)=f(\gamma x)\) is nevertheless an invertible isometry
on \(L^p(\mathbb R_+,dx/x)\).

## Declarative factorized class

Phase Q uses the notation

\[
 \widetilde{\mathcal E}^{\,\mathrm{right}\text{-}\gamma}
 :=\{(a,d_\gamma):a\in
 \widetilde{\mathcal E}(\mathbb R_+,V(\mathbb R))\}
\]

only as a typed collection of ordered pairs. It is not declared to be an
algebra or a norm-closed class. The following facts alone are established:

- \((a,d_\gamma)\) represents the already bounded operator
  \(\operatorname{Op}(a)D_\gamma\);
- right composition is exact by direct substitution in the quantization;
- the order \((a,d_\gamma)\) is retained;
- boundedness follows from boundedness of \(\operatorname{Op}(a)\) and the
  isometry of \(D_\gamma\);
- none of this implies \(a d_\gamma\in\widetilde{\mathcal E}\).

Phase Q originally recorded seven open obligations. Phase R discharges the
specific left-composition obligation by exact covariance. Six remain: sums
with distinct frequencies, products of factorized elements, involution, a
Fredholm symbol, norm closure, and later interaction with localized
Wiener--Hopf factors. The last item remains only a future obligation and is not
studied in either calculation.

## Proven damping lemma

Let \(h\in V(\mathbb R)\cap L^1(\mathbb R)\), with \(h\) absolutely
continuous. Then \(d_\gamma h\) is absolutely continuous and

\[
 (d_\gamma h)'
 =d_\gamma\bigl(h'+i\log(\gamma)h\bigr)
 \quad\text{a.e.}
\]

Therefore

\[
 \boxed{
 \operatorname{Var}(d_\gamma h)
 \leq
 \operatorname{Var}(h)
 +|\log\gamma|\,\|h\|_{L^1}
 }.
\]

Also \(|d_\gamma h|=|h|\), so zero limits of \(h\) at either endpoint are
preserved. In particular, the finite endpoint limits supplied by bounded
variation must be zero when \(h\in L^1\).

This is an independently proved elementary lemma. It is not applied to any
localized Wiener--Hopf factor here: an applicable \(L^1\) damping factor has
not been identified. In particular, the argument never uses the false
inclusion \(V(\mathbb R)\cap C_0(\mathbb R)\subset L^1(\mathbb R)\).
