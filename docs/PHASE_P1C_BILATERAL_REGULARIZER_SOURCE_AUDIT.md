# P1-C source audit: bilateral diagonal regularizer

## Decision

The article itself supports the strongest admissible P1-C classification,
subject to its stated Fredholm hypotheses:

`CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER`.

The certified candidate is the ordered word

\[
 B_1:=\mathcal T_{1,-}\mathcal R_1Z_1^{-1}.
\]

This statement does **not** identify \(B_1\) with \(\mathcal R_1\), turn the
whole word into one Mellin pseudodifferential operator, prove common-algebra
membership, or construct a Schur-correction symbol.

## Authoritative notation and the P1-C alias

The article calls the Mellin core \(\mathcal N_1\), not \(\mathcal H_1\).
Therefore this phase uses the declared alias

\[
 \mathcal H_1:=\mathcal N_1
\]

only at the P1-C API boundary. No occurrence of `H_1` or `mathcal H_1` was
found in the article source.

## Factor-by-factor audit

| Object | Article definition or theorem | Space and status |
|---|---|---|
| \(\mathcal A_{1,1}\) | `sections/06_shur_reduction.tex:27-42`: \(\widetilde V_{\alpha_1}P^++G_1P^-\) | Bounded scalar operator on \(X=L^p(\mathbb R_+,w_\delta)\) under the article assumptions |
| \(Z_1=M_{\zeta_1}\) | `sections/05_haseman_operator_with_shift.tex:339-408`: \(\zeta_1(x)=\gamma_1^{-\kappa}(\gamma_1x-ic_1)/(x-ic_1)\) and \(\widetilde V_1=Z_1U_1\) | Exact multiplication interface on \(X\) |
| \(Z_1^{-1}\) | Same passage: \(\zeta_1\) extends continuously, has positive nonzero endpoint values, and is bounded away from zero for \(\gamma_1>0\), \(0<\kappa<1\) | Exact bounded inverse multiplication by \(\zeta_1^{-1}\) on \(X\) |
| \(\widehat{\mathcal A}_{1,1}\) | `sections/05_haseman_operator_with_shift.tex:429-507`: \(U_1P^++\widehat G_1P^-\), with \(\mathcal A_{1,1}=Z_1\widehat{\mathcal A}_{1,1}\) | Exact; the proof explicitly does not commute \(Z_1\) through a projection |
| \(\mathcal T_{1,-}\) | `sections/05_haseman_operator_with_shift.tex:515-529`: \(U_1^{-1}P^++P^-\) | Exact noncommutative sum on \(X\), not a bare shift |
| \(\mathcal T_{1,-}^{-1}\) | `sections/05_haseman_operator_with_shift.tex:532-782`, especially `770-778`: \(\mathcal T_{1,+}\Phi^{-1}\operatorname{Op}(\theta_1^{-1})\Phi\) | Exact bounded inverse on \(X\); both auxiliary operators have index zero |
| \(\mathcal H_1:=\mathcal N_1\) | `sections/05_haseman_operator_with_shift.tex:792-900` and `sections/06_shur_reduction.tex:109-132`: \(\mathcal N_1=\widehat{\mathcal A}_{1,1}\mathcal T_{1,-}\), \(\Phi\mathcal N_1\Phi^{-1}=\operatorname{Op}(n_1)\) | Exact core and exact Mellin transport; the API alias is local to P1-C |
| \(\mathcal R_1\) | `sections/06_shur_reduction.tex:134-171`: \(\mathcal R_1=\Phi^{-1}\operatorname{Op}(r_{1,1})\Phi\) | Two-sided regularizer of \(\mathcal N_1\) modulo \(\mathcal K(X)\) under the Fredholm hypotheses at `101-107` |
| \(B_1\) | `sections/05_haseman_operator_with_shift.tex:1336-1398`; specialized at `sections/06_shur_reduction.tex:173-195` | Ordered candidate \(\mathcal T_{1,-}\mathcal R_1Z_1^{-1}\) |

## Exact link identity

The source gives two exact identities:

\[
 \mathcal A_{1,1}=Z_1\widehat{\mathcal A}_{1,1},\qquad
 \mathcal N_1=\widehat{\mathcal A}_{1,1}\mathcal T_{1,-}.
\]

Consequently, using the exact inverse of \(Z_1\),

\[
 Z_1^{-1}\mathcal A_{1,1}\mathcal T_{1,-}
 =\mathcal N_1=\mathcal H_1.
\]

This is `ARTICLE_INTERNAL_EXACT`. It has no compact ideal because it is an
operator equality. The equivalent orientation
\(\mathcal A_{1,1}\mathcal T_{1,-}=Z_1\mathcal H_1\) follows by exact left
multiplication and cancellation, not by name matching or commutation.

## Mellin-core regularizer

Under the article's Fredholm assumptions, `sections/06_shur_reduction.tex:134-159`
records both Mellin-side relations

\[
 \operatorname{Op}(n_1)\operatorname{Op}(r_{1,1})\simeq I,
 \qquad
 \operatorname{Op}(r_{1,1})\operatorname{Op}(n_1)\simeq I.
\]

After exact transport by \(\Phi\), lines `161-171` state that \(\mathcal R_1\)
is a two-sided regularizer of \(\mathcal N_1\):

\[
 \mathcal H_1\mathcal R_1-I\in\mathcal K(X),\qquad
 \mathcal R_1\mathcal H_1-I\in\mathcal K(X).
\]

The general proposition states these two hypotheses explicitly at
`sections/05_haseman_operator_with_shift.tex:1361-1376`.

## Independent full-product checks in the article

The general proof at `sections/05_haseman_operator_with_shift.tex:1400-1537`
sets \(B=\mathcal T_{k,-}\mathcal R_kZ_k^{-1}\) and computes independently

\[
 \mathcal A_{k,k}B-I
 =Z_k(\mathcal N_k\mathcal R_k-I)Z_k^{-1}\in\mathcal K(X)
\]

at lines `1473-1493`, and

\[
 B\mathcal A_{k,k}-I
 =\mathcal T_{k,-}(\mathcal R_k\mathcal N_k-I)
  \mathcal T_{k,-}^{-1}\in\mathcal K(X)
\]

at lines `1495-1529`. The specialization repeats both products in
`sections/06_shur_reduction.tex:196-260`.

## Compact ideal and typing

The common ideal is

\[
 \mathcal K(X),\qquad X=L^p(\mathbb R_+,w_\delta).
\]

`sections/06_shur_reduction.tex:10-12` fixes modulo-compact equivalence on the
scalar component space (or the corresponding matrix space where applicable).
Every factor in the diagonal proof acts on the same scalar \(X\). The Mellin
realization uses \(\Phi\) only inside an exact conjugation; it does not change
the space attached to the transported operator \(\mathcal R_1\).

## Existing implementation map

| Concept | Existing class | State before P1-C | P1-C action |
|---|---|---|---|
| \(\mathcal A_{1,1}\) | `OperatorAtom` plus typed block models | No diagonal regularizer model | Add a typed diagonal operator wrapper |
| \(\mathcal T_{1,-}\) | `AuxiliaryOperator` | Exact definition retained, no inverse certificate | Reuse and attach a typed exact-inverse interface |
| \(\mathcal H_1\) | `MellinPseudodifferentialOperator` and transport metadata | No generic transported core wrapper | Add a core model retaining article alias and Mellin data |
| \(\mathcal R_1\) | `RegularizerOperator`, `RegularizerMellinRepresentation`, `TransportedMellinCore` | Representation exists; sidedness absent | Reuse and add an oriented core-regularizer certificate |
| \(Z_1\), \(Z_1^{-1}\) | `MultiplicationOperator` | Typed, but no paired inverse proof | Reuse and add a nonvanishing exact-inverse interface |
| Modulo compact | `ModCompactEquivalence` | Available | Reuse unchanged with a typed bilateral ideal descriptor |
| Exact equality | `ExactIdentity` | Available | Reuse unchanged |
| Noncommutative word | `Product`/`LinearCombination` | Available | Reuse unchanged; never commute or auto-expand |
| Bilateral regularizer | None | Missing | Add an explicit two-product certificate |

`src/symbolic_operator_calculus/expressions.py` does not exist. The repository's
authoritative noncommutative expression layer is
`src/symbolic_operator_calculus/operators.py`; Mellin scalar expressions live
in `src/symbolic_operator_calculus/mellin/expressions.py`.

## Source integrity

The article repository was inspected read-only. No article, bibliography, or
literature-index file is modified by P1-C.
