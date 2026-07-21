# Phase P integrated audit: minimal closure lemma

Audit date: 2026-07-20. This report was completed before changing production
code for Phase P.

## 1. Symbolic repository checkpoint

```text
pwd
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git rev-parse --show-toplevel
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git branch --show-current
p1a-mellin-symbol-infrastructure

git rev-parse HEAD
b39c96c4fc28898918c536bb3081344dd3e0b403

git rev-parse --abbrev-ref --symbolic-full-name @{upstream}
origin/p1a-mellin-symbol-infrastructure

git rev-list --left-right --count @{upstream}...HEAD
0  0

git status --short
?? docs/FULL_CODE_AUDIT.md

git diff --cached --name-only
<empty>
```

The local HEAD, branch, upstream name, and allowed dirty entry match the gate.
The upstream had independently advanced to the same HEAD, so the observed
divergence is `0/0`, rather than the historical `0/6`; no push was performed in
this phase.

The protected file remains byte-identical:

```text
sha256  5d6797bf0ee751dd164f42a7af7438b5226c9e0c20ff79ef80ce78ea7d7eb4b1
size    57033 bytes
mtime   2026-07-16 19:32:24.724846953 -0600
mode    664
```

## 2. Phase N and O history

```text
4f37d85 feat(schur): derive normalized first pivot formally
9d7db86 feat(semantics): track relation kinds and proof obligations
b7436bd feat(schur): export normalized pivot derivation
ca73104 feat(schur): normalize complete first-pivot correction
1ab9192 feat(analysis): catalog noncommutative correction terms
b39c96c docs(schur): analyze complete first-pivot correction
```

## 3. Baseline

Commands were executed with the repository `.venv/bin` at the front of `PATH`:

```text
pytest
903 passed in 44.05s

python -m compileall src tests
Listing 'src'...
Listing 'src/symbolic_operator_calculus'...
Listing 'src/symbolic_operator_calculus/mellin'...
Listing 'src/symbolic_operator_calculus.egg-info'...
Listing 'tests'...

ruff check .
All checks passed!

git diff --check
<no output>
```

Gate A is open.

## 4. External repository stability

The literature index was inspected read-only at branch `main`, HEAD
`3d5840f08a115fd620f1a59a66ad886ef07259c3`, with an empty
`git status --short`. Its HEAD and clean status were unchanged after `resolve`,
`theorem show`, `search`, `pages`, `theorem validate`, `use audit`, and
`related first-schur-pivot`.

The paper was inspected read-only at branch `main`, HEAD
`fbd2a8a8124a1f064640f85a439529ec15f3665f`, with an empty
`git status --short`. Its HEAD and clean status were unchanged after inspection.

`litctl theorem validate` reported `Theorem catalog is valid.` The use audit
keeps every closure-relevant transfer open except the source-level/direct use of
Karlovich 2017, Theorem 3.5, for an already admissible Mellin PDO. In
particular, semiproduct, fixed-singularity, cusp localization, mixed-algebra,
and regularizer transfers retain unchecked hypotheses.

## 5. Verified evidence consumed

All five canonical keys resolved to local PDFs. The relevant checksums are:

| Key | Checksum | Verified records inspected |
|---|---|---|
| `Karlovich2017HasemanSO` | `9f92e37599da89f88cbdcdb85d70ef69645aecdc987013f2b54fb878915e9091` | Thms. 3.5, 5.3, 5.6, 6.2, 7.1 |
| `KKL2014TwoShifts` | `3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec` | Thm. 1.1; Lemma 2.7; Thm. 3.3; Lemma 3.4; Thm. 3.6; Lemmas 4.4–4.5; Thm. 5.2 |
| `Karlovich2009MellinPDO` | `b37ae39fe0e49d5efd67341382b8d1a6df44b05579cffae7a18dd312ba485c99` | Thm. 4.3, scalar scope verified |
| `KKL2014Regularization` | `7ed69d096f9e31983cb9b4701e8b486fd15736d82163144f27f9bae244ce41f6` | Thm. 5.8 |
| `Karlovich2025Cusps` | `301f8435688416e6ccd88be1a44f85e12372c9a557bdc1dbea301318b846be23` | Lemmas 3.1–3.2; Cor. 3.3; Thms. 5.1, 6.1, 6.2 |

The CLI resolution modes were: manual match for
`Karlovich2017HasemanSO`, exact DOI for `KKL2014TwoShifts`, manual match for
`Karlovich2009MellinPDO`, exact title for `KKL2014Regularization`, and exact
DOI for `Karlovich2025Cusps`. Every enumerated result returned status
`verified` from `litctl theorem show`.

`litctl related first-schur-pivot` classified Karlovich 2017, Theorem 3.5,
as direct only for an operator already known to be an admissible Mellin PDO;
the remaining Karlovich 2017 results are analogies. It retained the listed
KKL 2014, Karlovich 2009, regularization, and 2025 cusp uses as pending
adaptations. The output did not supply a closure rule for any of the four
present words.

The exact hypotheses, pages, conclusions, application kinds, and limitations
will be exported in `FIRST_PIVOT_EVIDENCE_MANIFEST.yaml`. No theorem was inferred
from a search result or prior report.

## 6. Paper definitions relevant to closure

The paper works on `L^p(R_+,w_delta)` with `1<p<infinity` and
`-1/p<delta<1-1/p`. The isometry

```text
Phi_delta f(x) = x^(delta+1/p) f(x)
```

maps this space onto `L^p(R_+,dmu)`. Under it, the Cauchy factors `P^±` are
Mellin convolutions with `V(R)` symbols
`(1 ± coth(pi(lambda+i(delta+1/p))))/2`.

`U_1` is the normalized constant dilation
`gamma_1^(delta+1/p)V_{gamma_1}` and is an invertible isometry. Its Mellin
multiplier is `gamma_1^(i lambda)`, which the paper explicitly states is not in
`V(R)` when `gamma_1 != 1`.

Assuming `G_1 in SO(R_+)`, the normalized coefficient `Ghat_1` is a
multiplication operator with coefficient in `SO(R_+)`. The regularizer is

```text
R_1 = Phi_delta^(-1) Op(r_1) Phi_delta,
```

where `r_1` belongs to the strengthened Mellin class
`E-tilde(R_+,V(R))`; it is a two-sided regularizer modulo compacts for the
normalized diagonal product.

The exact normalized right block is

```text
W_{1,2}^+ = chi_infinity W(chi_+(.) exp((c_1-c_2)(.))) chi_infinity I.
```

It acts between the branch-2 and branch-1 copies of the weighted half-line
space. Its unlocalized kernel comes from
`(tau-x+i(c_1-c_2))^(-1)/(pi i)`, and the relation with the transported
off-diagonal Cauchy block is only modulo compact operators.

## 7. The central non-identification

`KKL2014TwoShifts`, Lemmas 4.4–4.5, concern the fixed-singularity operator

```text
(R_y f)(t) = (pi i)^(-1) integral_0^infinity
             (t/tau)^(1/y-1/p) f(tau)/(tau+t) dtau,
```

with Mellin multiplier `1/sinh(pi(lambda+i/y))`. This is neither the kernel nor
the localization structure of `W_{1,2}^+`. Therefore those lemmas are an
`UNPROVED_ADAPTATION`, not an identification or an applicable closure rule.

Likewise, Karlovich 2017, Theorem 5.6, gives an exact Mellin PDO for a shifted
off-diagonal block on a logarithmic star. The present geometry is a pure cusp,
uses a power-weighted transport, and produces a localized Wiener–Hopf block.
The result is methodological `ANALOGY`, not a specialization already proved.

## 8. Existing public API and minimal extension

The existing `OperatorAtom`, `Product`, `Term`, `LinearCombination`, and
`compose_ordered` types are sufficient to represent the four words without a
new AST. Phase O already provides factor signatures, relation kinds, interface
statuses, source-aware operator rules, and proof-obligation records.

The minimum reusable extension is declarative metadata over existing `Product`
objects:

1. four stable closure-word records;
2. typed evidence references containing the mandatory provenance fields;
3. a dependency graph whose nodes reuse analytic proof-obligation semantics;
4. candidate-lemma and decision records that cannot represent a candidate as
   proved unless all blocking obligations are discharged.

No symbolic product, membership, compactness, or Fredholm inference engine is
justified. `R_1` must remain an atom, and no code may create a Mellin symbol for
`W_{1,2}^+`.

## 9. Audit conclusion

All gate conditions are satisfied. Phase P may proceed with the narrow
declarative extension above. The expected outcome remains conditional: the
evidence makes H1 formulable as a minimal candidate once the two `B_1` factors
are represented in an admissible Mellin calculus, but it does not certify H1,
H2, or H3. The final decision must be derived from the explicit obligation
graph rather than assumed here.
