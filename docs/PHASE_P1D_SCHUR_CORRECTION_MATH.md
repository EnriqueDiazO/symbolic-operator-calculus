# P1-D mathematics: exact first Schur correction

## Objects

Fix branch copies

\[
 X_j=L^p(\mathbb R_+,w_\delta),\qquad j=1,2,
\]

while retaining their branch labels. Define

\[
 D_k:=\widetilde V_{\alpha_k}-\mathcal G_k,
 \quad
 L_{2,1}:=D_2\mathcal W^-_{2,1}:X_1\to X_2,
 \quad
 R_{1,2}:=D_1\mathcal W^+_{1,2}:X_2\to X_1.
\]

P1-C supplies the certified structured regularizer

\[
 B_1:=\mathcal T_{1,-}\mathcal R_1Z_1^{-1}:X_1\to X_1
\]

with both

\[
 \mathcal A_{1,1}B_1-I\in\mathcal K(X_1),
 \qquad
 B_1\mathcal A_{1,1}-I\in\mathcal K(X_1).
\]

## Three exact views of one correction

The grouped P1-C insertion is

\[
 C_{2,2}^{(1)}=L_{2,1}B_1R_{1,2}.
 \tag{G}
\]

The article's block grouping is

\[
 C_{2,2}^{(1)}=B^-_{2,1}\mathcal R_1B^+_{1,2},
 \quad
 B^-_{2,1}=L_{2,1}\mathcal T_{1,-},
 \quad
 B^+_{1,2}=Z_1^{-1}R_{1,2}.
 \tag{S}
\]

Substitution, with no commutation or distribution, gives the same expanded
word from (G) and (S):

\[
 D_2\mathcal W^-_{2,1}\mathcal T_{1,-}\mathcal R_1Z_1^{-1}
 D_1\mathcal W^+_{1,2}.
 \tag{E}
\]

The assembly certificate records (G), (S), and (E) separately. Equality of
their typed ASTs is required; matching display names is not evidence.

## Exactness versus modulo compactness

The equations (G)-(E) concern the model correction and are exact under the
source definitions and the certified P1-C candidate representation. The
original off-diagonal block substitutions have a different strength:

\[
 \mathcal A_{2,1}\simeq-L_{2,1},
 \qquad
 \mathcal A_{1,2}\simeq+R_{1,2}
 \pmod{\mathcal K}.
\]

Since bounded left and right multiplication preserve the bilateral compact
ideal,

\[
 \mathcal A_{2,1}B_1\mathcal A_{1,2}
 \simeq -C_{2,2}^{(1)}.
\]

Applying the external Schur minus yields

\[
 \mathcal A_{2,2}-\mathcal A_{2,1}B_1\mathcal A_{1,2}
 \simeq \mathcal A_{2,2}+C_{2,2}^{(1)}.
 \tag{P}
\]

Thus `CERTIFIED_EXACT_SCHUR_ASSEMBLY` describes (G)-(E), while (P) is a
`ModCompactEquivalence` in the same typed ideal. Neither relation proves
Fredholmness of the right-hand side.

## Cutoff ledger

Each \(\mathcal W^\pm\) contains a left and a right
\(M_{\chi_\infty}\). The four actual audit records have status
`NO_REPLACEMENT_NEEDED` because all four multipliers remain in (E).

If a later derivation requests replacement, exact covariance can transport a
cutoff:

\[
 V_\gamma M_\chi V_{\gamma^{-1}}
 =M_{\chi_\gamma},\qquad \chi_\gamma(x)=\chi(\gamma x).
\]

This identity changes the cutoff scale and is not a theorem that the old and
new cutoffs are equal or compactly equivalent. A requested replacement must
therefore fail unless a separate exact or modulo-compact certificate is
provided.

## Relative-dilation ledger

Relative dilations are audited in the original ordered tree. A reciprocal
pair may be cancelled only if it is an actual adjacent pair, or reduced by a
specific certified covariance trace. The implementation never distributes
the differences \(D_k\), never chooses a summand of \(\mathcal T_{1,-}\),
and never commutes intervening factors to manufacture adjacency.

For the article word, the possible \(U_1^{-1}\)/\(U_1\) occurrence is
separated by a projection, \(\mathcal R_1\), \(Z_1^{-1}\), multipliers, and
linear combinations. Its audit is

```text
BLOCKED_BY_INTERVENING_FACTORS
```

Other branch dilations without a typed inverse partner receive
`NO_RELATIVE_DILATION_PAIR`. These statuses do not weaken the exact assembly:
they block only a proposed simplification or symbol derivation that needs the
pair resolved.

## Algebra and symbol result

The correction is a bounded endomorphism of \(X_2\) by typed composition.
That fact is not membership in a named inverse-closed algebra. The actual
factor-coverage matrix leaves the general radial Mellin regularizer,
dilation-containing interfaces, mixed product closure, and multiplicative
symbol theorem without one common source. Therefore:

```text
assembly_status = CERTIFIED_EXACT_SCHUR_ASSEMBLY
algebra_status  = BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP
symbol_status   = BLOCKED_BY_ALGEBRA_MEMBERSHIP
fredholm_status = NOT_CERTIFIED
```

A synthetic fully sourced extension can exercise the successful symbol API,
but it does not change the article instance.

## Fail-closed inference order

The public constructor validates in this order:

1. typed branch/domain/codomain and weighted-space compatibility;
2. exact definitions of both exterior factors and source blocks;
3. P1-C bilateral central-regularizer status and common compact ideal;
4. equality of the grouped, source-factorized, and expanded ASTs;
5. sign and off-diagonal modulo-compact evidence;
6. cutoff and dilation audit records;
7. common-algebra membership;
8. only then, symbol prerequisites.

An earlier failure is not masked by a later attractive formula. In
particular, a candidate symbol is absent whenever algebra membership is not
certified.

