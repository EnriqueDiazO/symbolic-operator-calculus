# P1-B cutoff control

## Exact covariance

Let the normalized dilation act by

\[
  (V_\gamma f)(x)=\gamma^\nu f(\gamma x),\qquad
  (V_{\gamma^{-1}}f)(x)=\gamma^{-\nu}f(x/\gamma),
\]

on one fixed weighted space, with `gamma>0`.  For a multiplication cutoff
`M_chi`, direct evaluation gives

\[
\begin{aligned}
 (V_\gamma M_\chi V_{\gamma^{-1}}f)(x)
 &=\gamma^\nu\chi(\gamma x)\gamma^{-\nu}f(x)\\
 &=\chi(\gamma x)f(x).
\end{aligned}
\]

Hence

\[
  V_\gamma M_\chi V_{\gamma^{-1}}=M_{\chi_\gamma},
  \qquad \chi_\gamma(x)=\chi(\gamma x),
\]

exactly.  The normalization factors cancel, but the radial scale does not.
For the article's cutoff, if `chi=0` on `[0,1]` and `chi=1` on
`[2,infinity)`, then `chi_gamma=0` on `[0,1/gamma]` and is one on
`[2/gamma,infinity)`.  Its coordinate system, support description, endpoint
role, and new radial scale must be retained as data.

The article independently displays these scaled multipliers in the actual
localized kernel: `chi_infinity(gamma_k x)` and
`chi_infinity(gamma_j y)` in
`sections/05_haseman_operator_with_shift.tex:2172-2278`.

## What the verified sources prove

`sections/03_transport_and_kernels.tex:734-773`, specializing the verified
Karlovich 2025 cusp localization result, proves

\[
  T_{k,j}-M_{\chi_\infty}T_{k,j}M_{\chi_\infty}\in\mathcal K.
\]

The verified literature-index report
`reports/MIXED_LOCAL_CUSP_ALGEBRA_ROUTE_DECISION.md` records Karlovich 2025,
Corollary 3.3 as supporting the present two-cutoff Wiener--Hopf models modulo
the original compact localization.  That result requires retaining both
cuts.  The same report concludes `MIXED_LOCAL_ALGEBRA_NOT_AVAILABLE` for the
larger word.  DLS1995 Lemma 2.9 concerns a different equal-order cusp
equivalence and is only a specialization/analogy for the present coordinates.
The unavailable 2025 mixed-localization primary is not inferred.

None of these verified results states, for the mixed P1-B products,

\[
  (M_{\chi(\gamma\cdot)}-M_\chi)K\in\mathcal K,
  \qquad
  K(M_{\chi(\gamma\cdot)}-M_\chi)\in\mathcal K.
\]

Compactness cannot be deduced merely from a picture of separated supports.

## Conditional replacement lemma

The implementation may record, but must not prove, the following schema.

**Cutoff replacement modulo compact operators (conditional).**  Let `K` be a
bounded operator on the stated weighted space and let `chi_1`, `chi_2` be
typed cutoffs in the same coordinate system.  If a separately supplied
analytic hypothesis proves the relevant one-sided products, then:

- left: `(M_chi1-M_chi2)K` compact implies
  `M_chi1 K ~= M_chi2 K`;
- right: `K(M_chi1-M_chi2)` compact implies
  `K M_chi1 ~= K M_chi2`;
- both: both statements may be recorded only when both are evidenced.

The claim has independent fields `status=ASSUMED|CERTIFIED|UNPROVED`,
`ideal`, and `side=LEFT|RIGHT|BOTH`.  `ASSUMED` produces only a formal,
conditional relation.  `CERTIFIED` requires nonempty source and evidence and
may produce a modulo-compact relation.  `UNPROVED` produces no rewrite.

## Classification

For the actual P1-B cores, cutoff replacement is
`BLOCKED_BY_CUTOFF_CONTROL`.  This is a secondary obstruction to the earlier
five-factor regularizer-interface conflict.  The API must never equate the
original and scaled cutoff automatically.
