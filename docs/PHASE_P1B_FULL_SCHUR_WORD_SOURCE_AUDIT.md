# P1-B source audit: full first-Schur word

## Scope and verdict

This audit is intentionally literal.  The requested input word

\[
  \mathcal B_{2,1}^{-}\mathcal T_{1,-}\mathcal R_1
  Z_1^{-1}\mathcal B_{1,2}^{+}
\]

is **not** the word denoted by the same symbols in the article.  In the
article, `B^-_{2,1}` already ends in `T_{1,-}` and `B^+_{1,2}` already starts
with `Z_1^{-1}`.  Consequently, expanding the requested five-factor word
duplicates both auxiliary interfaces.  This is a source-level conflict, not
an algebraic simplification opportunity, and P1-B must fail closed with
`BLOCKED_BY_REGULARIZER_INTERFACE` for that literal input.

The article itself instead proves the exact model identity

\[
  C_{2,2}^{(1)}
  =\mathcal B_{2,1}^{-}\mathcal R_1\mathcal B_{1,2}^{+}.
\]

No statement below changes either convention.

## Authoritative locations

All paths in this section are under
`/home/enriquedo/PersonalProjects/Papers/DiazOcampoHasemanBVP2026`.

| Object | Source | Strength | Literal content |
|---|---|---|---|
| Weighted space and Mellin transport | `sections/03_transport_and_kernels.tex:1064-1108` | exact | `Phi_delta: L^p(R_+,w_delta) -> L^p(R_+,dmu)`, `Phi_delta f=x^(delta+1/p)f` |
| Half-line Cauchy factors | `sections/03_transport_and_kernels.tex:1114-1197` | exact | `P^pm=(I pm S_R+)/2`; after `Phi_delta` they are Mellin convolutions |
| Cutoff | `sections/03_transport_and_kernels.tex:747-762` | definition plus compact localization | `chi_infty=0` on `[0,1]`, `=1` on `[2,infinity)` and `T_kj-chi T_kj chi I` is compact |
| Branch dilation | `sections/05_haseman_operator_with_shift.tex:86-121` | exact | `alpha_k(x)=gamma_k x`, `gamma_k>0`, `V_k f(x)=f(gamma_k x)`; orientation preserving |
| Transported shift | `sections/05_haseman_operator_with_shift.tex:123-160` | exact | `Vtilde_k=M_{rho_k}V_k`, `rho_k=(gamma_k x-i c_k)/(x-i c_k)` |
| Normalized dilation and multiplier | `sections/05_haseman_operator_with_shift.tex:339-408` | exact | `U_k=gamma_k^(delta+1/p)V_k`, `U_k^-1=U_{gamma_k^-1}`, `Vtilde_k=Z_kU_k` |
| `Z_k` | `sections/05_haseman_operator_with_shift.tex:357-388` | exact | multiplication by `zeta_k=gamma_k^(-delta-1/p)rho_k`; bounded and invertible |
| Auxiliary `T_{k,-}` | `sections/05_haseman_operator_with_shift.tex:515-529` | exact definition | `T_{k,-}=U_k^-1 P^+ + P^-` |
| Full diagonal regularizer | `sections/05_haseman_operator_with_shift.tex:1336-1397` | exact representation; regularizing property modulo compacts | `A_kk^(-1):=T_{k,-} R_k Z_k^-1` |
| Normalized off-diagonal models | `sections/06_shur_reduction.tex:44-95` | modulo compacts from actual entries | `A12 ~= (Vtilde_1-G_1)W^+_12`, `A21 ~= -(Vtilde_2-G_2)W^-_21` |
| Mellin core `R_1` | `sections/06_shur_reduction.tex:109-171` | exact definition after existential regularizer symbol | `R_1=Phi_delta^-1 Op_M(r_11) Phi_delta` |
| Full first-block regularizer | `sections/06_shur_reduction.tex:173-195` | exact representation | `A11^(-1)=T_{1,-} Phi^-1 Op_M(r_11) Phi Z_1^-1` |
| Expanded model correction | `sections/06_shur_reduction.tex:359-380` | exact for `C`; contribution to pivot modulo compacts | seven ordered factors displayed below |
| Article `B^-` | `sections/06_shur_reduction.tex:382-390` | exact definition | `(Vtilde_2-G_2)W^-_21 T_{1,-}` |
| Article `B^+` | `sections/06_shur_reduction.tex:391-399` | exact definition | `Z_1^-1(Vtilde_1-G_1)W^+_12` |
| Article three-factor form | `sections/06_shur_reduction.tex:401-456` | exact for model correction | `C=B^- R_1 B^+` |

## Factor-by-factor article word

The unsimplified article word is

\[
  (\widetilde V_{\alpha_2}-\mathcal G_2)\,
  \underbrace{M_{\chi_\infty}W(b^-_{2,1})M_{\chi_\infty}}_{W^-_{2,1}}\,
  (U_1^{-1}P^+ + P^-)\,
  \Phi_\delta^{-1}\operatorname{Op}_M(r_{1,1})\Phi_\delta\,
  Z_1^{-1}\,
  (\widetilde V_{\alpha_1}-\mathcal G_1)\,
  \underbrace{M_{\chi_\infty}W(b^+_{1,2})M_{\chi_\infty}}_{W^+_{1,2}}.
\]

Every factor acts on `L^p(R_+,w_delta)` except the central
`Op_M(r_11)`, which acts on `L^p(R_+,dmu)` between the explicitly displayed
`Phi_delta^-1` and `Phi_delta`.  Each scalar block has the same scalar
domain/codomain; the off-diagonal matrix compatibility is `A21: X_1 -> X_2`,
`A11^(-1): X_1 -> X_1`, and `A12: X_2 -> X_1`, with the branch copies
isomorphic to the same weighted half-line space.

| Pos. | Factor | Kind and retained metadata | Dilation/cutoff information |
|---:|---|---|---|
| 1 | `Vtilde_2-G_2` | ordered linear combination of transported shift and coefficient multiplication | `Vtilde_2=M_{rho_2}V_{gamma_2}`; no commutation is implicit |
| 2 | `W^-_21` | localized Wiener--Hopf operator | left and right `M_chi`; `b^-_21=lambda<0 exponential` |
| 3 | `T_{1,-}` | ordered linear combination | `U_1^-1 P^+ + P^-`; the inverse dilation occurs only in the first summand and precedes `P^+` |
| 4 | `R_1` | transported Mellin PDO | `Phi^-1 Op_M(r_11) Phi`; it is not the full diagonal regularizer |
| 5 | `Z_1^-1` | invertible coefficient multiplication | multiplication by `zeta_1^-1`; cannot be absorbed across adjacent factors without a rule |
| 6 | `Vtilde_1-G_1` | ordered linear combination | `Vtilde_1=M_{rho_1}V_{gamma_1}` |
| 7 | `W^+_12` | localized Wiener--Hopf operator | left and right `M_chi`; `b^+_12=lambda>0 exponential` |

Thus all four cutoff multipliers in the correction remain present.  In the
article's direct conjugation computation they become
`chi_infty(gamma_k x)` and `chi_infty(gamma_j y)`; see
`sections/05_haseman_operator_with_shift.tex:2172-2278`.  The relative scale
is `beta_{k,j}(x)=(gamma_k/gamma_j)x`, with this orientation proved in
`sections/05_haseman_operator_with_shift.tex:1876-1892`.

## Literal requested word versus article word

Expanding only definitions, without cancelling or commuting, gives

\[
\begin{aligned}
 &\mathcal B^-_{2,1}\mathcal T_{1,-}\mathcal R_1
 Z_1^{-1}\mathcal B^+_{1,2}\\
 &=(\widetilde V_{\alpha_2}-\mathcal G_2)W^-_{2,1}
 \mathcal T_{1,-}\,\mathcal T_{1,-}\,
 \Phi_\delta^{-1}\operatorname{Op}_M(r_{1,1})\Phi_\delta\,
 Z_1^{-1}\,Z_1^{-1}
 (\widetilde V_{\alpha_1}-\mathcal G_1)W^+_{1,2}.
\end{aligned}
\]

This differs structurally from the exact article word.  Neither
`T_{1,-}^2=T_{1,-}` nor `Z_1^{-2}=Z_1^{-1}` is stated or true in general.
The five-argument API must therefore retain the original five atoms, record
the expanded duplication, and reject the article-labelled instance by typed
factor provenance rather than by string matching.

## Previous-phase certification table

| Element | Certified | Conditional | Not proved |
|---|:---:|:---:|:---:|
| Multiplication covariance under a typed dilation | yes |  |  |
| Wiener--Hopf covariance with rescaled symbol | yes |  |  |
| Mellin-PDO covariance with rescaled radial symbol | yes |  |  |
| Cancellation of an explicitly adjacent inverse pair | yes |  |  |
| Chosen `R_1=Phi^-1 Op_M(r_11)Phi` interface | yes, as the article's exact definition | existence depends on the stated Fredholm hypotheses |  |
| Application of cancellation to the real pivot word |  |  | yes |

## Sign convention

The algebraic Schur pivot is
`A22^(1)=A22-A21 A11^(-1) A12` at
`sections/06_shur_reduction.tex:281-305`.  Because the normalized model of
`A21` carries a leading minus sign (`:86-95`), substitution produces the
positive displayed model correction (`:310-348`).  Hence the article's
`C22^(1)` is the positive model expression, whereas the raw product
`A21 A11^(-1) A12` is its negative modulo compact operators.  These are two
different typed sign layers.

## Audit classification

The earliest fatal issue for the requested five-factor call is
`BLOCKED_BY_REGULARIZER_INTERFACE`.  Even after using the article's corrected
three-factor interface, extracting inverse relative dilations remains blocked
by the ordered sums, projections, multipliers, and retained cutoffs.  No
Fredholm or pivot conclusion follows from this audit.
