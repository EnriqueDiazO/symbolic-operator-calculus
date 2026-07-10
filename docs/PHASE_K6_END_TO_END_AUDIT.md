# Phase K6B audit: end-to-end closure of the first reduced m=2 block

Audit date: 2026-07-10.

## 1. Repository state and baseline

```text
pwd
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git rev-parse --show-toplevel
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git branch --show-current
main

git log -1 --oneline
001d719 feat: add compact action for first reduced Schur model

git status --short
<clean>

git diff --cached --name-only
<empty>

upstream
origin/main
```

`origin` uses
`git@github-personal:EnriqueDiazO/symbolic-operator-calculus.git` for fetch and
push.

Baseline:

```text
PYTHONPATH=src pytest -q
265 passed in 28.28s

python3 -m compileall -q src tests
passed

ruff check .
All checks passed!

git diff --check
passed (no output)
```

## 2. Complete chain map

| Step | Input | Output | Module | Source of truth | Rewritten information | Exact conversion? | Automatic/manual data |
|---|---|---|---|---|---|---|---|
| 1 | `a22_first_schur_reduction()` | `FirstSchurReduction` | `blocks.py` / `relations.py` | Exact block metadata | m=2 indices and exterior `-1` are encoded once | No | Automatic constructor |
| 2 | `a22_first_schur_mod_compact_relation()` | `ModCompactSchurRelation` | `blocks.py` | Exact reduction plus expanded model constructors | No formula duplicated | No | Automatic |
| 3 | `a22_first_schur_model()` | Eight-term `LinearCombination` | `blocks.py` | `a22_exact_operator()` followed by the correction | Terms are concatenated, not re-entered | No | Automatic |
| 4 | `a22_first_schur_correction()` | Four-term positive correction | `blocks.py` | Schur sign, A21 leading sign, `main_expression()` | No product rewritten in K4B | No | Automatic |
| 5 | `factor_first_schur_correction()` | `FirstSchurCorrectionFactorization` | `kernels.py` | The K4B correction | Splits products around their unique `R11`; does not spell factors manually | No | Automatic; optional invalid input can be supplied for validation |
| 6 | Factorization | LEFT / `R11` / RIGHT | `kernels.py` | Ordered AST prefixes and suffixes | None | No | Automatic Cartesian validation |
| 7 | LEFT | scalar `M21` | `kernels.py` + `actions.py` | Apply LEFT, extract `KernelCombination`, project locally | No manual M21 formula | No | Caller supplies variables and optional rules |
| 8 | RIGHT | scalar `M12` | `kernels.py` + `actions.py` | Apply RIGHT, extract `KernelCombination`, project locally | No manual M12 formula | No | Caller supplies variables and optional rules |
| 9 | Derived M21/M12 | nested `C22` kernel | `kernels.py` | Derived side kernels plus formal `R11(u,v)` | Only the compact nesting is constructed | No | Dummy variables automatic or explicit |
| 10 | `C22`, input function | `Integral(C22(x,y) f(y),dy)` | `kernels.py` | `combined_kernel_c22()` | None | No | Function class/input variable required |
| 11 | `f(x)` and rules | `FirstSchurCompactModelAction` | `schur.py` | `a22_exact_operator()` and `apply_combined_kernel_c22()` | No A22/M/C formula duplicated | No | Rules and variables optional |

The chain is operational and traceable. The only intentionally local scalar
bridge not derived by applying an AST side is the unevaluated formal center
`R11(u,v)` used inside the compact kernel.

## 3. Exact level

`a22_first_schur_reduction()` returns:

```text
FirstSchurReduction(
  diagonal=ExactBlock("A",2,2),
  left=ExactBlock("A",2,1),
  regularizer=FormalRegularizer(ExactBlock("A",1,1), R11),
  right=ExactBlock("A",1,2),
  offdiagonal_product_coefficient=-1)
```

The probe confirmed that `diagonal`, `left`, and `right` are all
`ExactBlock`, not `OperatorAtom`. `FirstSchurReduction` has no `apply`,
`as_expr`, or `rewrite`. No public dispatcher accepts it. The exact meaning
remains only `A22 - A21 R11 A12` metadata.

## 4. Modulo-compact relation

`a22_first_schur_mod_compact_relation()` returns a frozen
`ModCompactSchurRelation` whose `exact` is the first reduction and whose
`model` is the expanded eight-term `LinearCombination`. Probes established:

```text
relation.exact == a22_first_schur_reduction()  True
relation.model == a22_first_schur_model()      True
relation.exact == relation.model               False
has rewrite/as_expr                            False
```

### Relation/action classification

This is **case B**. `apply_a22_first_schur_model_compact()` stores the relation
but constructs the action through the semantically specialized sources
`a22_exact_operator()` and `apply_combined_kernel_c22()`, rather than traversing
`relation.model` directly. K6A tests prove expanded/compact equivalence in both
formal and explicit modes, and separately assert that the stored relation is
the public relation constructor.

The frozen result type validates the type of `relation`, not semantic equality
of arbitrary manually supplied `diagonal` and `correction` fields with its
model. The public factory is synchronized; direct dataclass construction can
create inconsistent metadata. This is a bounded validation gap, not an exact
conversion.

## 5. Expanded model

The coefficient probe produced exactly:

```text
(Fraction(1,2), Fraction(1,2), Fraction(1,2), Fraction(-1,2),
 1, -1, -1, 1)
```

`model.terms[:4] == a22_exact_operator().terms` and
`model.terms[4:] == a22_first_schur_correction().terms` are both true. The
first four products are the ordered diagonal expansion; the final four retain
the noncommutative correction order. Construction concatenates AST tuples and
never inspects a SymPy `Add`.

## 6. Schur sign

The production computation remains:

```text
Schur exterior sign            -1
leading A21_WH sign            -1
correction sign product        +1
correction coefficients        (1,-1,-1,1)
```

`correction_sign` multiplies `main_expression()` and therefore controls the
returned correction. Neither K5B nor K6A adds a separate global correction
sign. K5B derives side signs from the correction coefficients, while K6A
delegates to the compact K5B action.

## 7. Structural factorization

`factor_first_schur_correction()` calls
`a22_first_schur_correction()` when no explicit validation input is provided.
It locates the unique `R11` in each term, preserves prefix/suffix order,
requires the ordered Cartesian product, and reconstructs the correction before
returning.

The result is:

```text
LEFT  coefficients (1,-1)
  Vtilde_alpha2 Wminus_21
  G2 Wminus_21

CENTER
  R11

RIGHT coefficients (1,-1)
  Vtilde_alpha1 Wplus_12
  G1 Wplus_12
```

The reconstruction `LEFT * R11 * RIGHT` is structurally equal to the K4B
correction. No second hard-coded factorization exists in production.

## 8. M21 and M12

`m21_kernel()` obtains LEFT from the factorization and passes it to
`_kernel_expression_from_factor()`. That helper applies the ordered
combination, calls `extract_applied_kernels()`, and uses
`KernelCombination.as_expr()` only as the final scalar projection. `m12_kernel`
does the same with RIGHT.

No occurrences of `rho1`, `rho2`, `G1_scalar`, `G2_scalar`, `gamma1`,
`gamma2`, `Lplus_12_kernel`, or `Lminus_21_kernel` remain in `kernels.py` as
manual production constants. The formal formulas emerge from the selected
action mapping. Tests also mutate the factor signs and verify that M21/M12
follow them.

## 9. Combined C22 kernel

`combined_kernel_c22()` calls the derived `m21_kernel()` and `m12_kernel()` and
constructs:

```text
Integral(
  M21(x,u) * Integral(
    R11(u,v) * M12(v,y),
    (v,0,oo)),
  (u,0,oo))
```

M21 is outside, `R11*M12` is in the inner `v` integral, and the outer integral
is in `u`. Formal expansion has coefficients `(1,-1,-1,1)` and no spurious
`1/2`, `2`, or `-2`. There is no parallel four-term formula.

There is one residual nominal duplication: `R11_kernel = Function("R11")` is
declared in both `actions.py` and `kernels.py`. M21/M12 are action-derived, but
the compact center is not extracted by applying the `R11` atom. Both uses are
formal and tests bind them to the same public `R11` semantics, but the string
name could theoretically drift.

## 10. Formal and Fourier-explicit paths

Both paths invoke exactly the same public chain; only the rules mapping differs.

### Formal (`rules=None`)

The compact result contains `Lplus_12`, `Lminus_21`, and `R11`. It has four
ordered diagonal terms, one correction input integral, inert `PrincipalValue`,
free symbols `gamma1,gamma2,x`, and bound symbols `y,u,v` under ordinary names.

### Fourier explicit

With `explicit_wiener_hopf_rules(decay=d)`, the compact result contains
normalized Fourier kernels and `chi_infinity`, has no `Lplus_12` or
`Lminus_21`, and continues to contain formal `R11`. Its free symbols are
`d,gamma1,gamma2,x`; bound symbols remain `y,u,v`.

This is not a completely evaluated kernel: only the off-diagonal Wiener-Hopf
factors are Fourier explicit.

## 11. Compact model action

`apply_a22_first_schur_model_compact()` accepts only an undefined input
function applied as `f(variable)`. It resolves one rules mapping and passes it
to both branches:

```text
diagonal  <- apply_linear_combination_ordered(a22_exact_operator(), ...)
correction<- apply_combined_kernel_c22(...)
```

The frozen `FirstSchurCompactModelAction` keeps `relation`, ordered `diagonal`,
and one outer correction `Integral` separate. The diagonal retains
`PrincipalValue`. No specialized eight-term action formula exists in
`schur.py`; the expanded route remains the generic linear-combination engine.

## 12. Expanded/compact equivalence

K6A parameterizes the end-to-end comparison over formal and explicit rules.
It compares:

```text
apply_linear_combination_ordered(a22_first_schur_model(), ...)
```

against the compact result. It establishes:

1. the four compact diagonal `AppliedTerm` objects equal the first four
   expanded terms;
2. extraction from the expanded correction yields four kernels with
   coefficients `(1,-1,-1,1)`;
3. expanding the compact `M21*R11*M12` integrand reproduces those four kernels
   term by term using algebraic comparison;
4. the same formal `R11` remains in both modes.

No invalid structural equality is requested between eight terms and the
`4 + Integral` representation.

## 13. Variable safety

The input variable is chosen in `schur.py` with `fresh_symbol()` while avoiding
the external variable and operand symbols. `combined_kernel_c22()` then chooses
outer and middle variables with the same shared policy, adding each chosen
symbol to the avoid set.

Observed defaults:

| External variable | Free external symbol retained | Bound symbols | Correction input |
|---|---|---|---|
| `x` | `x` | `y,u,v` | `y` |
| `y` | `y` | `u,v,w` | `u` |
| `u` | `u` | `y,v,w` | `y` |
| `v` | `v` | `y,u,w` | `y` |

All three internal variables are distinct and disjoint from final free symbols.
Explicit valid `input=y, outer=u, middle=v` is preserved. Explicit
`input=external`, `outer=input`, and `middle=outer` collisions raise
`CompactSchurActionError` or `ValueError` as appropriate.

The diagonal principal-value integrals and compact correction integrals have
separate binding scopes. `collect_bound_symbols()` sees their union; safe
substitution tests and the external-name probes show no capture.

## 14. `FirstSchurCompactModelAction.as_expr()`

The method is exactly:

```python
return self.diagonal.as_expr() + self.correction
```

It is a local scalar projection. No production code calls this method to
reconstruct an AST, choose operator order, rewrite a relation, or define exact
equality. Searches found its uses only in K6A tests for final projection,
free/bound-symbol inspection, and coexistence checks. Ordered meaning remains
in the `diagonal` field and the correction stays a separate field before
projection.

## 15. Dependency map

```text
operators     -> (none)
fourier       -> (none)
substitution  -> (none)
actions       -> fourier, operators, substitution
relations     -> operators
blocks        -> operators, relations
kernels       -> actions, blocks, operators, substitution
schur         -> actions, blocks, kernels, operators, relations, substitution
package init  -> all public modules
```

There are no import cycles, delayed imports, or runtime import tricks. No lower
layer imports `schur.py`. `kernels.py` contains the specific m=2 correction
factorization because that structure is required to derive kernels; it does not
contain unrelated block metadata or compact diagonal action logic. `schur.py`
delegates rather than repeating block or kernel formulas.

## 16. Public surface for the first reduced Schur block

| Layer | Public API |
|---|---|
| Exact | `FirstSchurReduction`, `a22_first_schur_reduction` |
| Expanded model | `a22_first_schur_model`, `a22_first_schur_correction` |
| Relation | `ModCompactSchurRelation`, `a22_first_schur_mod_compact_relation` |
| Factorization | `FirstSchurCorrectionFactorization`, `SchurFactorizationError`, `factor_first_schur_correction` |
| Kernels | `KernelTerm`, `KernelCombination`, `m21_kernel`, `m12_kernel`, `c22_integrand`, `combined_kernel_c22`, `apply_combined_kernel_c22` |
| Compact model action | `FirstSchurCompactModelAction`, `CompactSchurActionError`, `apply_a22_first_schur_model_compact` |

The names distinguish exact reduction, model, relation, kernels, and compact
model action. `a22_exact_operator()` denotes exact A22, while
`a22_first_schur_reduction()` denotes exact reduced A22. The only mild
historical ambiguity is `combined_kernel_c22`: its name does not say “first
Schur” or “model”, although its module context and docstring identify its role.
No current API falsely calls the compact model action exact.

## 17. End-to-end test coverage

No single test explicitly names and asserts every one of the eleven audit steps
from exact reduction through relation, factorization, derived sides, kernels,
and compact action.

The longest composed K6A test starts from `a22_first_schur_model()`, applies its
eight terms, constructs the compact action, extracts the correction kernels,
and compares them with the expanded compact C22 in both rule modes. Internally,
the compact C22 calls M21/M12, which call factorization, which starts from the
K4B correction. Separate K4B tests bind relation exact/model to the public
constructors, and K5B tests explicitly traverse correction through
factorization and C22.

Thus behavioral coverage spans the complete chain across a small set of
composable tests, but there is no single narrative traceability test beginning
with `a22_first_schur_reduction()` and ending with the compact action.

## 18. Classified technical debt

| Severity | Finding | Classification rationale |
|---|---|---|
| CRITICAL | None | No mathematical inconsistency, forbidden exact conversion, broken mode, or capture defect was found. |
| HIGH | None | Both production paths and expanded/compact equivalence are covered. |
| MEDIUM | Formal `R11_kernel` is declared separately in `actions.py` and `kernels.py` | Residual nominal second source; currently coherent and formal, but a name drift would not be prevented by construction. |
| LOW | `FirstSchurCompactModelAction` validates relation type but not semantic consistency of manually supplied fields | Public factory is correct and tested; only direct dataclass construction can create mismatch. |
| LOW | No one test enumerates the entire eleven-step semantic chain | Coverage exists across K4B/K5B/K6A; this is traceability clarity, not missing behavior. |
| LOW | `combined_kernel_c22` has a historically broad name | Context is unambiguous in current m=2 scope; no functional risk. |
| NO ES DEUDA | Separate `ModCompactRelation` and `ModCompactSchurRelation` | Their validated endpoint types and semantics differ intentionally. |
| NO ES DEUDA | Compact action is case B rather than traversing `relation.model` | It preserves compactness and has formal/explicit equivalence tests. |
| NO ES DEUDA | `as_expr()` creates a SymPy sum | It is only a final scalar projection and is not an AST source. |
| NO ES DEUDA | `R11` remains formal in explicit mode | This is the stated mathematical scope. |

## 19. Closure assessment

The m=2 first reduced block is coherently closed at the requested level:

- exact metadata remains non-applicable;
- modulo-compact semantics remain explicit;
- the expanded model preserves exact ordering and coefficients;
- the Schur sign controls the correction;
- M21/M12 derive from the correction AST;
- formal and Fourier-explicit C22 use the same pipeline;
- the compact action combines exact A22 action with compact model correction;
- expanded and compact representations agree termwise;
- dummy-variable selection is capture-safe.

No critical second source of truth remains. The duplicated formal R11 scalar
name is a medium-strength traceability issue, not a competing formula for the
reduction.

## 20. Required recommendation

**OPTION A: the m=2 case is coherently closed. The next phase should focus on
representation, explanation, or export of the derivation, not new algebra.**

There is no small missing mathematical connection requiring option B, no
architectural inconsistency requiring option C, and no mathematical blockade
supporting option D. Recursive Schur reduction, arbitrary `m`, and general
operator matrices remain deliberately out of scope.
