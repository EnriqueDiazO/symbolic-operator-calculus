# Phase Q audit: exact right-dilation composition

Audit date: 2026-07-21. This file was created before changing Phase Q
production code.

## Symbolic repository checkpoint

```text
pwd
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git root
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

branch
p1a-mellin-symbol-infrastructure

HEAD
c60aedb1d1fa59b985251b2caae8695585774044

upstream
origin/p1a-mellin-symbol-infrastructure

divergence (upstream, local)
0  0

git status --short
?? docs/FULL_CODE_AUDIT.md
```

The Phase P commits are:

```text
77d9e06 feat(analysis): integrate literature-backed closure obligations
c60aedb docs(schur): select minimal first-pivot closure lemma
```

The protected audit remained untracked and byte-identical:

```text
sha256  5d6797bf0ee751dd164f42a7af7438b5226c9e0c20ff79ef80ce78ea7d7eb4b1
```

## External checkpoints

Before and after the read-only consultation:

- literature index: branch `main`, HEAD
  `3d5840f08a115fd620f1a59a66ad886ef07259c3`, clean;
- paper: branch `main`, HEAD
  `fbd2a8a8124a1f064640f85a439529ec15f3665f`, clean.

No external HEAD changed.

## Baseline

```text
pytest
918 passed in 46.46s

python -m compileall src tests
completed successfully

ruff check .
All checks passed!

git diff --check
<no output>
```

## Existing public API

The repository already provides:

- the noncommutative `OperatorAtom` and `Product` AST, including atomic `R1`
  and `U1`;
- `MellinExpression`, `MellinSymbol`, explicit variable roles, domains,
  assumptions, singularities, and deterministic renderers;
- `build_dilation_multiplier`, which constructs
  \(\gamma^{i\lambda}\) under the explicit branch hypothesis \(\gamma>0\)
  but intentionally makes no \(V(\mathbb R)\)-membership claim;
- disjoint `ExactIdentity`, `ModCompactEquivalence`, `OperatorRule`, and
  `AnalyticProofObligation` types;
- the Phase P closure graph and result/interface metadata.

The API does not implement the action of a general Mellin PDO and does not
currently distinguish a standard Mellin symbol from a standard symbol followed
by an oscillatory dilation factor. This limitation must be documented before
adding the narrow semantic wrapper required by Phase Q. No second operator AST
is needed.

## Verified paper conventions

The paper defines

\[
 (\mathcal Mf)(\lambda)
 =\int_0^\infty f(x)x^{-i\lambda}\frac{dx}{x},
 \qquad
 (\mathcal M^{-1}g)(x)
 =\frac1{2\pi}\int_{\mathbb R}g(\lambda)x^{i\lambda}\,d\lambda,
\]

and the left/Kohn--Nirenberg Mellin quantization

\[
 [\operatorname{Op}(a)f](x)
 =\frac1{2\pi}\int_{\mathbb R}\int_0^\infty
 a(x,\lambda)\left(\frac{x}{y}\right)^{i\lambda}
 f(y)\frac{dy}{y}\,d\lambda.
\]

On the weighted space,

\[
 \kappa_\delta=\delta+\frac1p,
 \qquad
 (\Phi_\delta f)(x)=x^{\kappa_\delta}f(x),
 \qquad
 U_\gamma=\gamma^{\kappa_\delta}V_\gamma,
 \quad
 (V_\gamma f)(x)=f(\gamma x).
\]

The paper proves exactly

\[
 (\Phi_\delta U_\gamma\Phi_\delta^{-1}g)(x)=g(\gamma x)
\]

and, with the displayed Mellin convention,

\[
 \mathcal M[g(\gamma\,\cdot)](\lambda)
 =\gamma^{i\lambda}(\mathcal Mg)(\lambda).
\]

Thus the normalized dilation carries no residual scalar weight factor and the
sign is \(+i\lambda\). The unnormalized \(V_\gamma\) would retain the separate
factor \(\gamma^{-\kappa_\delta}\).

## Source evidence for the coefficient case

The literature CLI was used according to its evidence contract.

- `KKL2014TwoShifts`, Theorem 3.3, printed p. 942/PDF p. 8, checksum
  `3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec`:
  \(\operatorname{Op}(a)\operatorname{Op}(b)\simeq
  \operatorname{Op}(ab)\) for \(a,b\in\mathcal E(\mathbb R_+,V(\mathbb R))\).
- The same source, Lemma 3.4, printed p. 942/PDF p. 8: an exact triple
  product under separated dependence; it is not needed for the two-factor
  coefficient semiproduct.
- `KKL2014Regularization`, Theorem 5.8, printed p. 206/PDF p. 18, checksum
  `7ed69d096f9e31983cb9b4701e8b486fd15736d82163144f27f9bae244ce41f6`:
  the supplied regularizer symbol lies in the same strengthened Mellin class
  once the theorem's hypotheses have been established.

The paper proves separately that, under \(G_1\in SO(\mathbb R_+)\), the
lambda-independent symbol \((x,\lambda)\mapsto\widehat G_1(x)\) lies in
\(\widetilde{\mathcal E}(\mathbb R_+,V(\mathbb R))\). Since \(\Phi_\delta\)
and \(M_{\widehat G_1}\) are multiplication operators, conjugation leaves
that coefficient operator unchanged. These facts make Theorem 3.3 a concrete
specialization candidate, to be recorded separately from the dilation case.

## Gate A conclusion

All checkpoints and the baseline pass, the exact conventions are known, and
the protected file is unchanged. Phase Q may proceed with a small declarative
extension. Direct substitution \(z=\gamma y\) in the inner integral predicts
option A:

\[
 \operatorname{Op}(a)U_\gamma
 =\operatorname{Op}(a\,d_\gamma),
 \qquad d_\gamma(\lambda)=\gamma^{i\lambda},
\]

after conjugation to \(L^p(\mathbb R_+,d\mu)\). This is an operator identity
to be proved on a dense test class and extended by continuity; it is not a
claim that \(a\,d_\gamma\) belongs to the standard symbol class.
