# P1-D source audit: first Schur correction

## Decision

The active article supports the exact model correction

\[
 \mathcal C_{2,2}^{(1)}
 =\mathcal L_{2,1}\,
   (\mathcal T_{1,-}\mathcal R_1 Z_1^{-1})\,
   \mathcal R_{1,2},
 \tag{1}
\]

where

\[
 \mathcal L_{2,1}
 =(\widetilde V_{\alpha_2}-\mathcal G_2)\mathcal W^-_{2,1},
 \qquad
 \mathcal R_{1,2}
 =(\widetilde V_{\alpha_1}-\mathcal G_1)\mathcal W^+_{1,2}.
\]

P1-C independently certifies the central parenthesis as a two-sided
regularizer of \(\mathcal A_{1,1}\) modulo the common ideal
\(\mathcal K(X)\). Therefore the strongest assembly classification is

```text
CERTIFIED_EXACT_SCHUR_ASSEMBLY
```

for the model correction itself. The passage from the original Schur pivot
to the model pivot remains modulo compact operators. No source located in
the article or audited literature proves that the complete seven-factor word
belongs to one inverse-closed algebra carrying a Fredholm symbol.

## Authoritative article locations

All article paths below are under
`/home/enriquedo/PersonalProjects/Papers/DiazOcampoHasemanBVP2026` and were
inspected read-only.

| Object | Source | Strength |
|---|---|---|
| Scalar and matrix modulo-compact convention on \(X=L^p(\mathbb R_+,w_\delta)\) | `sections/06_shur_reduction.tex:10-12` | definition |
| Diagonal blocks \(\mathcal A_{1,1}\), \(\mathcal A_{2,2}\) | `sections/06_shur_reduction.tex:27-42` | exact definitions |
| Raw localized off-diagonal blocks | `sections/06_shur_reduction.tex:44-59` | modulo compact |
| Normalized signs and \(\mathcal W^+_{1,2}\), \(\mathcal W^-_{2,1}\) | `sections/06_shur_reduction.tex:60-95` | exact definitions inside modulo-compact block models |
| P1-C diagonal regularizer \(\mathcal T_{1,-}\mathcal R_1Z_1^{-1}\) | `sections/06_shur_reduction.tex:109-260` | two-sided modulo compact under stated hypotheses |
| Raw Schur pivot \(\mathcal A_{2,2}-\mathcal A_{2,1}\mathcal A_{1,1}^{(-1)}\mathcal A_{1,2}\) | `sections/06_shur_reduction.tex:281-305` | exact definition |
| Positive normalized correction in the pivot | `sections/06_shur_reduction.tex:310-354` | modulo compact as part of the pivot |
| Seven-factor correction word | `sections/06_shur_reduction.tex:359-380` | exact identity for the model correction |
| Article factors \(\mathcal B^-_{2,1}\), \(\mathcal B^+_{1,2}\) | `sections/06_shur_reduction.tex:382-399` | exact definitions |
| \(\mathcal C_{2,2}^{(1)}=\mathcal B^-_{2,1}\mathcal R_1\mathcal B^+_{1,2}\) | `sections/06_shur_reduction.tex:401-456` | exact |
| Conditional Fredholm equivalence with the pivot | `sections/06_shur_reduction.tex:460-607` | modulo compact; no pivot Fredholmness conclusion |

The article's names `B^-` and `B^+` already contain the interfaces
`T_{1,-}` and `Z_1^{-1}`:

\[
 \mathcal B^-_{2,1}=\mathcal L_{2,1}\mathcal T_{1,-},
 \qquad
 \mathcal B^+_{1,2}=Z_1^{-1}\mathcal R_{1,2}.
\]

Thus the exact source form

\[
 \mathcal B^-_{2,1}\mathcal R_1\mathcal B^+_{1,2}
\]

and the P1-D grouped form (1) expand to the same ordered word. This avoids
the interface duplication diagnosed in P1-B.

## Exact expanded word

Only definition substitution and association are used:

\[
\begin{aligned}
 \mathcal C_{2,2}^{(1)}={}&
 (\widetilde V_{\alpha_2}-\mathcal G_2)
 \mathcal W^-_{2,1}
 \mathcal T_{1,-}
 \mathcal R_1
 Z_1^{-1}
 (\widetilde V_{\alpha_1}-\mathcal G_1)
 \mathcal W^+_{1,2}.
\end{aligned}
\tag{2}
\]

The implementation must retain both differences as noncommutative linear
combinations and must not distribute (2). It must also keep the P1-C
regularizer as a certified structured object in the grouped representation;
expansion is an evidence ledger, not a replacement for that certificate.

## Domains and branch orientation

Let \(X_j\) denote branch \(j\)'s copy of the same weighted scalar space.
Then

\[
 \mathcal R_{1,2}:X_2\to X_1,
 \quad
 \mathcal A_{1,1}^{(-1)}:X_1\to X_1,
 \quad
 \mathcal L_{2,1}:X_1\to X_2.
\]

Consequently (1) is an endomorphism of \(X_2\). Branch compatibility is an
independent typed gate; equality of the underlying weighted-space parameters
does not erase the source and target branch labels.

## Sign audit

The exact algebraic pivot is

\[
 \mathcal A_{2,2}^{(1)}
 :=\mathcal A_{2,2}
   -\mathcal A_{2,1}\mathcal A_{1,1}^{(-1)}\mathcal A_{1,2}.
\]

The normalized block relations are

\[
 \mathcal A_{1,2}\simeq +\mathcal R_{1,2},
 \qquad
 \mathcal A_{2,1}\simeq -\mathcal L_{2,1}.
\]

Hence the two minus signs combine and

\[
 \mathcal A_{2,2}^{(1)}
 \simeq \mathcal A_{2,2}+\mathcal C_{2,2}^{(1)}.
\]

The raw product, its external Schur sign, and the positive model correction
are distinct data. The source proves a conditional Fredholm equivalence; it
does not prove the assembled pivot Fredholm.

## Cutoffs and dilations

Each localized Wiener--Hopf block is exactly

\[
 \mathcal W^\pm=M_{\chi_\infty}W(b^\pm)M_{\chi_\infty}.
\]

All four cutoff multipliers remain in (2). No cutoff replacement is needed
to establish the exact assembly. The localization theorem at
`sections/03_transport_and_kernels.tex:747-764` concerns replacement of the
original off-diagonal operator by its doubly localized model modulo compact
operators; it does not delete either cutoff from the model word.

The article proves exact dilation covariance and the orientation
\(\beta_{k,j}=\alpha_j^{-1}\circ\alpha_k\) in
`sections/05_haseman_operator_with_shift.tex:1876-1892` and
`:2140-2278`. Under conjugation, the cutoffs acquire branch scales rather
than becoming identical. Inside (2), the inverse dilation contained in the
first summand of \(\mathcal T_{1,-}\) is not adjacent to the dilation inside
\(\widetilde V_{\alpha_1}\); intervening projections, the Mellin
regularizer, multipliers, and linear combinations forbid exact cancellation.

## Article boundary

`sections/06_fredholm_program.tex:4-31` explicitly leaves construction of the
reduced algebra and Fredholm symbol to future work. The active article no
longer makes the older, unsupported full-pivot symbolic claim. P1-D follows
that operatorial route: it certifies (1)-(2), records the modulo-compact pivot
relation, and blocks algebra/symbol promotion.

