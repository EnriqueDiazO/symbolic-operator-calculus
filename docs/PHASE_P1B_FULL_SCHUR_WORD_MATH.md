# P1-B mathematics: full first-Schur word

## Result

The literal five-factor request is classified

```text
BLOCKED_BY_REGULARIZER_INTERFACE
```

This is not a failure of the abstract dilation-covariance identity.  It is a
source/interface mismatch that occurs earlier: the article's
`B^-_{2,1}` already contains `T_{1,-}`, and its `B^+_{1,2}` already contains
`Z_1^{-1}`.  The requested call therefore duplicates those interfaces.

The implementation nevertheless certifies exact and modulo-compact synthetic
full-word factorizations.  Those cases demonstrate that the API can expose
and cancel an actual adjacent inverse-dilation pair without claiming that the
article's word has that shape.

## Three operator levels

The analysis keeps three expressions separate.

1. The authoritative input AST is always

   \[
   B^-_{2,1}\,T_{1,-}\,R_1\,Z_1^{-1}\,B^+_{1,2}.
   \]

2. An expanded evidentiary AST records supplied block and regularizer
   representations.  For the article-labelled input it contains two
   `T_{1,-}` and two `Z_1^{-1}` factors.

3. A factorized AST exists only after every typed interface, domain,
   orientation, cutoff claim, and relation strength has passed validation.

Blocked results have no `factorized_word` and no semantic relation.  This
prevents an attractive candidate formula from being mistaken for a derived
identity.

## Exact cutoff covariance

For

\[
 (V_\gamma f)(x)=\gamma^\nu f(\gamma x),\qquad
 (V_{\gamma^{-1}}f)(x)=\gamma^{-\nu}f(x/\gamma),
\]

direct calculation gives

\[
 V_\gamma M_\chi V_{\gamma^{-1}}
 =M_{\chi_\gamma},
 \qquad \chi_\gamma(x)=\chi(\gamma x).
\]

`SpatialCutoffOperator` retains the symbol, support description, regions on
which the cutoff is zero/one, coordinate system, and radial scale.
`conjugate_cutoff_by_dilation` returns a new cutoff and an `ExactIdentity`;
it does not return a boolean and does not equate the new cutoff with the old
one.

A `CutoffEquivalenceClaim` is a separate object.  `UNPROVED` enables no
rewrite, `ASSUMED` produces only a `FormalIdentity`, and `CERTIFIED` with
evidence produces a `ModCompactEquivalence` in the named ideal.  A two-sided
localization needs either a `BOTH` claim or separate `LEFT` and `RIGHT`
claims.  This makes side-specific compactness an input rather than an
inference from support geometry.

## Typed real factors

The full-word layer adds only metadata absent from the existing operator
model:

- `SpatialCutoffOperator` and `LocalizedOperator` store the exact ordered
  product `M_left K M_right`;
- `CauchyProjectionOperator` retains polarity and its exact
  Mellin-convolution model;
- `CoefficientMultiplicationOperator` distinguishes coefficient roles;
- `TransportedShiftOperator` stores the exact ordered factorization
  `M_rho V_gamma`;
- `AuxiliaryOperator` can store the genuine linear combination
  `U_1^{-1}P^+ +P^-`;
- `TransportedMellinCore` stores
  `Phi_delta^{-1} Op_M(r_11) Phi_delta` and never identifies this central
  core with `T_{1,-}R_1Z_1^{-1}`;
- `FullWordBlock` augments the existing `SchurBlockRepresentation` with typed
  interface provenance and cutoff claims.

Existing `WienerHopfOperator`, `WeightedDilationOperator`,
`MellinPseudodifferentialOperator`, `RegularizerMellinRepresentation`,
`Product`, `LinearCombination`, and semantic relation types are reused.

## The only successful core pattern

The exact synthetic path accepts the contiguous supplied models

\[
 B^-\,T_-=W_LV_\gamma,
 \qquad
 Z^{-1}B^+=V_{\gamma^{-1}}W_R,
 \qquad
 R=\Phi^{-1}\operatorname{Op}_M(r)\Phi.
\]

The expanded AST is therefore

\[
 W_LV_\gamma\Phi^{-1}\operatorname{Op}_M(r)\Phi
 V_{\gamma^{-1}}W_R.
\]

On the weighted realization used by the previous phase, the contiguous
middle subproduct is delegated to
`cancel_relative_dilations_in_schur_correction`.  The certified covariance
is

\[
 V_\gamma\operatorname{Op}_M(r)V_{\gamma^{-1}}
 =\operatorname{Op}_M(r^\gamma),
 \qquad r^\gamma(x,\lambda)=r(\gamma x,\lambda).
\]

No search, commutation, or distribution is performed to create the pattern.
If a projection or extra multiplication separates the pair, the result is
`BLOCKED_BY_OPERATOR_ORDER`.  If a localized exterior lacks complete claims,
the result is `BLOCKED_BY_CUTOFF_CONTROL`.  An assumed or untransported
regularizer representation gives `BLOCKED_BY_REGULARIZER_INTERFACE`.

## Relation strength

- all component relations exact: `FULL_WORD_EXACT_FACTORIZATION` with an
  `ExactIdentity`;
- any certified modulo-compact block/regularizer input:
  `FULL_WORD_MOD_COMPACTS_FACTORIZATION`, only in one common named ideal;
- complete cutoff claims present:
  `FACTORIZATION_AFTER_CUTOFF_EQUIVALENCE`; assumed claims keep a
  `FormalIdentity`, while certified claims give only a modulo-compact
  relation;
- any missing prerequisite: one of the three blocking statuses, with no
  derived relation.

Compact-ideal strength, hypotheses, evidence, and the original five-factor
order are retained in the result.

## Sign layers

`RawSchurCorrection` stores the intrinsically positive ordered product
`A21 A11^(-1) A12`.  `SignedSchurContribution` applies exactly one external
sign.  `SchurPivot` adds that signed term to `A22`.  The article's positive
model `C22^(1)` results from combining the outer Schur minus with the minus
already present in the normalized model for `A21`; it is not silently used as
the raw-product convention.

## What remains open

The implementation always retains these analytic obligations:

```text
left_core_factorization
right_core_factorization
cutoff_replacement_mod_compacts
regularizer_mellin_representation
wh_mellin_wh_sandwich_membership
fredholm_algebra_membership
schur_correction_symbol
```

The synthetic tests validate the inference engine, not the article-level
premises.  The article correction, full pivot symbol, Fredholm-algebra
membership, and Fredholm conclusion remain open.
