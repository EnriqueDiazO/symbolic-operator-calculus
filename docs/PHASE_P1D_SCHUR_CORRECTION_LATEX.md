# P1-D LaTeX draft: first Schur correction

The following fragment deliberately separates the exact model assembly from
the modulo-compact pivot relation and from the blocked algebra/symbol claims.

```latex
\paragraph{Exact assembly of the first model correction.}
Put
\[
 \mathcal L_{2,1}
 :=(\widetilde V_{\alpha_2}-\mathcal G_2)\mathcal W^-_{2,1},
 \qquad
 \mathcal R_{1,2}
 :=(\widetilde V_{\alpha_1}-\mathcal G_1)\mathcal W^+_{1,2}.
\]
Using the two-sided factored regularizer of the first diagonal block,
\[
 \mathcal A_{1,1}^{(-1)}
 =\mathcal T_{1,-}\mathcal R_1Z_1^{-1},
\]
the first model correction is the exact ordered product
\[
\begin{aligned}
 \mathcal C_{2,2}^{(1)}
 &=\mathcal L_{2,1}\mathcal A_{1,1}^{(-1)}\mathcal R_{1,2}\\
 &=(\widetilde V_{\alpha_2}-\mathcal G_2)
   \mathcal W^-_{2,1}\mathcal T_{1,-}\mathcal R_1Z_1^{-1}
   (\widetilde V_{\alpha_1}-\mathcal G_1)\mathcal W^+_{1,2}.
\end{aligned}
\tag{P1-D.1}
\]
Equivalently, with the article definitions
\[
 \mathcal B^-_{2,1}=\mathcal L_{2,1}\mathcal T_{1,-},
 \qquad
 \mathcal B^+_{1,2}=Z_1^{-1}\mathcal R_{1,2},
\]
one has the exact factorization
\[
 \mathcal C_{2,2}^{(1)}
 =\mathcal B^-_{2,1}\mathcal R_1\mathcal B^+_{1,2}.
\tag{P1-D.2}
\]
No commutation, distribution, cutoff removal, or dilation cancellation is
used in (P1-D.1)--(P1-D.2).

\paragraph{Modulo-compact relation for the pivot.}
The source convention is
\[
 \mathcal A_{2,2}^{(1)}
 :=\mathcal A_{2,2}
   -\mathcal A_{2,1}\mathcal A_{1,1}^{(-1)}\mathcal A_{1,2}.
\]
Since
\[
 \mathcal A_{2,1}\simeq-\mathcal L_{2,1},
 \qquad
 \mathcal A_{1,2}\simeq+\mathcal R_{1,2}
 \pmod{\mathcal K(X)},
\]
the normalized model satisfies only
\[
 \mathcal A_{2,2}^{(1)}
 \simeq \mathcal A_{2,2}+\mathcal C_{2,2}^{(1)}
 \pmod{\mathcal K(X)}.
\tag{P1-D.3}
\]
The positive sign in (P1-D.3) is the product of the external Schur minus and
the minus in the normalized model of \(\mathcal A_{2,1}\).

\paragraph{Algebra-membership audit.}
All four cutoff multipliers in
\(\mathcal W^-_{2,1}\) and \(\mathcal W^+_{1,2}\) are retained. The audited
sources do not provide a single inverse-closed algebra that simultaneously
contains these localized Wiener--Hopf blocks, the general radial Mellin
pseudodifferential regularizer \(\mathcal R_1\), the nontrivial dilation
interfaces, and their exact ordered products. Accordingly the present
classification is
\[
 \operatorname{status}_{\rm alg}
 =\texttt{BLOCKED\_BY\_COMMON\_ALGEBRA\_SOURCE\_GAP}.
\tag{P1-D.4}
\]

\paragraph{Symbol gate.}
Because common-algebra membership and a multiplicative symbol homomorphism
for the mixed word have not been proved, no Fredholm symbol is assigned to
\(\mathcal C_{2,2}^{(1)}\). The result is
\[
 \operatorname{status}_{\rm sym}
 =\texttt{BLOCKED\_BY\_ALGEBRA\_MEMBERSHIP}.
\tag{P1-D.5}
\]
Equation (P1-D.5) is a blocking status, not a candidate symbolic formula.
Likewise, (P1-D.1) alone does not prove Fredholmness of the first pivot.
```

## Classification summary

| Statement | Strength |
|---|---|
| Equations (P1-D.1) and (P1-D.2) for the model correction | exact |
| Off-diagonal block substitutions and pivot equation (P1-D.3) | modulo compact operators in the same \(\mathcal K(X)\) |
| Existence of the P1-C regularizer symbol under the article hypotheses | conditional on those stated Fredholm hypotheses |
| Common mixed-algebra membership | blocked by a source gap |
| Complete correction symbol | not demonstrated; blocked by algebra membership |
| Fredholmness of the first pivot | not demonstrated |

