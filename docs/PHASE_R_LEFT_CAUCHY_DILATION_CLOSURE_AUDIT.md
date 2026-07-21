# Phase R audit: left Cauchy–dilation closure

Audit date: 2026-07-21. This report was created before changing Phase R
production code.

## Required checkpoints

```text
pwd / Git root
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

branch
p1a-mellin-symbol-infrastructure

HEAD
4f5ce22bbfbf55299c5b498b59fb0e015dcab4ee
docs(mellin): analyze regularizer-dilation closure

immediate parent
4e6a8ab533c525beb95214138ecbf66b99b34d0e
feat(mellin): model exact right-dilation composition

upstream
origin/p1a-mellin-symbol-infrastructure

divergence (upstream-only, local-only)
0  0

git status --short
?? docs/FULL_CODE_AUDIT.md
```

The protected file remained untracked and byte-identical:

```text
sha256  5d6797bf0ee751dd164f42a7af7438b5226c9e0c20ff79ef80ce78ea7d7eb4b1
```

Before and after the read-only consultation:

- literature index: branch `main`, HEAD
  `3d5840f08a115fd620f1a59a66ad886ef07259c3`, clean;
- paper: branch `main`, HEAD
  `fbd2a8a8124a1f064640f85a439529ec15f3665f`, clean.

No external HEAD or worktree changed.

## Baseline

```text
pytest
937 passed in 45.59s

python -m compileall src tests
completed successfully

ruff check .
All checks passed!

git diff --check
<no output>
```

## Existing public API

The repository already provides:

- the authoritative noncommutative `OperatorAtom`/`Product` AST with atomic
  `R1`, `U1`, `U1_inverse`, `P_plus`, `P_minus`, and `Ghat1`;
- `StandardMellinSymbol`, `OscillatoryMellinFactor`, and the ordered
  `FactorizedOscillatoryMellinSymbol` from Phase Q;
- exact `RightDilationCompositionTrace` for
  `Product((R1, U1))`;
- disjoint `ExactIdentity`, `ModCompactEquivalence`, `OperatorRule`, and
  `AnalyticProofObligation` types;
- the Phase P/Q closure graph and deterministic LaTeX/YAML renderers.

The Phase Q factorized-symbol validator is specialized to `R1 U1`, and there
is no typed record for radial scaling, left covariance, or the four truncated
left cores. A minimal generalization of that wrapper and one Phase R trace are
therefore needed. No second AST or general symbolic-composition engine is
needed.

## Verified paper definitions

The paper uses

\[
 (\mathcal Mf)(\lambda)=\int_0^\infty f(x)x^{-i\lambda}\frac{dx}{x}
\]

and

\[
 [\operatorname{Op}(a)f](x)
 =\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
 a(x,\lambda)(x/y)^{i\lambda}f(y)\frac{dy}{y}\,d\lambda.
\]

The weighted isometry and normalized dilation satisfy

\[
 (\Phi_\delta f)(x)=x^{\delta+1/p}f(x),
 \qquad
 \Phi_\delta U_\gamma\Phi_\delta^{-1}=D_\gamma,
 \qquad
 (D_\gamma f)(x)=f(\gamma x).
\]

The half-line Cauchy factors are

\[
 P^\pm=\tfrac12(I\pm S_{\mathbb R_+}),
 \qquad
 \Phi_\delta P^\pm\Phi_\delta^{-1}
 =\operatorname{Op}(p^\pm),
\]

where

\[
 p^\pm(\lambda)=\tfrac12\left(1\pm
 \coth(\pi(\lambda+i\kappa_\delta))\right)\in V(\mathbb R).
\]

They are not complementary projections:

\[
 P^-P^+=P^+P^-\ne0.
\]

The paper also supplies

\[
 R_1=\Phi_\delta^{-1}\operatorname{Op}(r_1)\Phi_\delta,
 \quad r_1\in\widetilde{\mathcal E},
 \qquad
 \widehat G_1(x)=\zeta_1(x)^{-1}G_1(x).
\]

Under its explicit hypothesis `G_1 in SO`, the lambda-independent symbol
`Ghat_1(x)` belongs to the same strengthened Mellin class.

## Symbol-class contracts to check under radial scaling

The paper defines `E-tilde` by the following cumulative requirements:

1. bounded continuous mapping from the radial variable into `V(R)`;
2. slow oscillation on multiplicative intervals `[t,2t]` at zero and infinity;
3. uniform continuity under covariable translations in the `V` norm;
4. uniform vanishing of the derivative tails outside `[-M,M]`.

Phase R must verify each contract after `x -> x/gamma`, rather than invoke
informal invariance.

## Literature evidence

The mandated CLI workflow (`resolve`, `search`, `pages`, `theorem show`, and
`use audit`) was executed from the literature-index repository.

- `KKL2014TwoShifts`, Theorem 3.3, printed p. 942/PDF p. 8, checksum
  `3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec`:
  if both symbols lie in `E(R+,V(R))`, their ordered Mellin-PDO product is the
  pointwise product symbol modulo compact operators.
- `KKL2014TwoShifts`, Lemma 3.4, printed p. 942/PDF p. 8, same checksum:
  an exact triple product is available under separated dependence. It is not
  needed for the two-factor Cauchy semiproduct and cannot be used with a
  non-`V` dilation multiplier.
- `Karlovich2017HasemanSO`, Theorem 3.5, printed p. 469/PDF pp. 7–9, checksum
  `9f92e37599da89f88cbdcdb85d70ef69645aecdc987013f2b54fb878915e9091`,
  was audited as required. It is a later Fredholm criterion and is not applied
  anywhere in Phase R.

The paper definitions prove that each `p^pm` is a spatially constant admissible
Mellin symbol. Consequently Theorem 3.3 is a direct specialization for
`P^pm Op(a)` once membership of `a` has been proved.

## Gate A conclusion

All checkpoints, external-stability checks, conventions, and baseline tests
pass. Direct evaluation at `x/gamma` predicts

\[
 D_\gamma^{-1}\operatorname{Op}(a)
 =\operatorname{Op}(a_\gamma)D_{\gamma^{-1}},
 \qquad
 a_\gamma(x,\lambda)=a(x/\gamma,\lambda),
\]

and hence

\[
 D_\gamma^{-1}\operatorname{Op}(a)D_\gamma
 =\operatorname{Op}(a_\gamma).
\]

These formulas still require a direct dense-core proof, sign check, bounded
extension, and semantic tests. Phase R may proceed with the minimal extension.
