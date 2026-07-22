# P1-C LaTeX draft

**Status:** certified bilateral draft under the article's stated Fredholm and
weighted-space hypotheses. The diagonal formula is source-complete. This
draft is **not** evidence for common-algebra membership, a single Mellin-PDO
representation of the full word, or a Schur-correction symbol.

```latex
% P1-C certified bilateral diagonal-regularizer draft.
% Scope: the diagonal operator only.  Do not cite this passage for a common
% algebra, a single Mellin-PDO representation, or a Schur symbol.

\begin{proposition}[Bilateral factored regularizer of the first diagonal]
Let
\[
  X=L^p(\mathbb R_+,w_\delta),
  \qquad
  \mathcal A_{1,1}
   =\widetilde V_{\alpha_1}P^+ +G_1P^-,
\]
and suppose the Fredholm and symbol hypotheses used for the regularization of
\(\mathcal N_1\) hold.  Put
\[
 \mathcal H_1:=\mathcal N_1,
 \qquad
 \mathcal N_1
   =\widehat{\mathcal A}_{1,1}\mathcal T_{1,-},
 \qquad
 \mathcal T_{1,-}=U_1^{-1}P^+ +P^-.
\]
Here \(\mathcal H_1\) is only an alias for the notation \(\mathcal N_1\)
used in the preceding construction.  Let
\[
 \mathcal R_1
   =\Phi^{-1}\operatorname{Op}(r_{1,1})\Phi
\]
be a two-sided regularizer of \(\mathcal H_1\) on \(X\), so that
\[
 \mathcal H_1\mathcal R_1-I\in\mathcal K(X),
 \qquad
 \mathcal R_1\mathcal H_1-I\in\mathcal K(X).
\]
If \(Z_1=M_{\zeta_1}\) is the invertible multiplier from the exact
factorization
\[
 \mathcal A_{1,1}=Z_1\widehat{\mathcal A}_{1,1},
\]
then the ordered operator
\[
 \mathcal B_1
   :=\mathcal T_{1,-}\mathcal R_1Z_1^{-1}
\]
is a two-sided regularizer of \(\mathcal A_{1,1}\):
\[
 \mathcal A_{1,1}\mathcal B_1-I\in\mathcal K(X),
 \qquad
 \mathcal B_1\mathcal A_{1,1}-I\in\mathcal K(X).
\]
\end{proposition}

\begin{proof}
The two exact definitions give the oriented link identity
\[
 Z_1^{-1}\mathcal A_{1,1}\mathcal T_{1,-}
 =\widehat{\mathcal A}_{1,1}\mathcal T_{1,-}
 =\mathcal H_1.
\]
Left multiplication by \(Z_1\), exact cancellation, and association yield
\[
\begin{aligned}
 \mathcal A_{1,1}\mathcal B_1-I
 &=\mathcal A_{1,1}\mathcal T_{1,-}\mathcal R_1Z_1^{-1}-I\\
 &=Z_1\bigl(\mathcal H_1\mathcal R_1-I\bigr)Z_1^{-1}
 \in\mathcal K(X).
\end{aligned}
\]
For the opposite product, right-multiply the link by the exact inverse
\(\mathcal T_{1,-}^{-1}\).  Then
\[
 Z_1^{-1}\mathcal A_{1,1}
 =\mathcal H_1\mathcal T_{1,-}^{-1},
\]
and hence, independently,
\[
\begin{aligned}
 \mathcal B_1\mathcal A_{1,1}-I
 &=\mathcal T_{1,-}\mathcal R_1Z_1^{-1}\mathcal A_{1,1}-I\\
 &=\mathcal T_{1,-}
   \bigl(\mathcal R_1\mathcal H_1-I\bigr)
   \mathcal T_{1,-}^{-1}
 \in\mathcal K(X).
\end{aligned}
\]
Both inclusions use the same bilateral ideal \(\mathcal K(X)\) and bounded
exactly invertible exterior factors.  No commutation of factors is used.
\end{proof}

\begin{remark}[Certification boundary]
The proposition does not assert that
\(\mathcal T_{1,-}\mathcal R_1Z_1^{-1}\) is one Mellin
pseudodifferential operator or that it belongs to a common operator algebra
with the off-diagonal Schur factors.  Consequently no symbol of the complete
Schur correction is asserted here.
\end{remark}
```

## Semantic ledger

- Exact: the definitions of \(Z_1\), \(Z_1^{-1}\),
  \(\mathcal T_{1,-}\), \(\mathcal T_{1,-}^{-1}\), the transport formulas,
  and the link identity.
- Modulo compact operators: the two core-regularizer relations and the two
  complete-product conclusions.
- Conditional: applicability retains the Fredholm, symbol-class, weight, and
  nonvanishing hypotheses stated by the article.
- Not proved: full-word common-algebra membership, a single Mellin-PDO model,
  cutoff replacement, the WH--Mellin--WH sandwich calculus, and a Schur symbol.
