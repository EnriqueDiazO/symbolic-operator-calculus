# Phase K5 audit: Schur reduction to combined kernel

## 1. Repository state and baseline

Audit date: 2026-07-09.

```text
pwd
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git rev-parse --show-toplevel
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git branch --show-current
main

git log -1 --oneline
44a3fe3 feat: model first Schur reduction modulo compacts

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

Baseline validations:

```text
PYTHONPATH=src pytest -q
233 passed in 2.87s

python3 -m compileall -q src tests
passed

ruff check .
All checks passed!

git diff --check
passed (no output)
```

## 2. Module map and real APIs

| Responsibility | Module | Real API |
|---|---|---|
| Noncommutative AST and four-term correction | `operators.py` | `OperatorAtom`, `Product`, `Term`, `LinearCombination`, `expand_ordered`, `main_expression` |
| Atomic and composed actions | `actions.py` | `mvp_atomic_rules`, `explicit_wiener_hopf_rules`, `apply_atom`, `apply_product`, `apply_linear_combination_ordered`, `apply_linear_combination` |
| Exact/off-diagonal/Schur block constructors | `blocks.py` | `a22_first_schur_reduction`, `a22_first_schur_mod_compact_relation`, `a22_first_schur_model`, `a22_first_schur_correction` |
| Semantic metadata and relations | `relations.py` | `ExactBlock`, `FormalRegularizer`, `FirstSchurReduction`, `WienerHopfModel`, `ModCompactRelation`, `ModCompactSchurRelation` |
| Normalized Fourier kernels | `fourier.py` | `kplus_kernel`, `kminus_kernel`, `localized_lplus_kernel`, `localized_lminus_kernel` |
| Kernel extraction and combined kernel | `kernels.py` | `extract_integral_kernel`, `extract_applied_kernels`, `m12_kernel`, `m21_kernel`, `c22_integrand`, `combined_kernel_c22`, `apply_combined_kernel_c22` |
| Capture-safe substitution | `substitution.py` | `substitute_free_variable`, `collect_bound_symbols`, `fresh_symbol`, `alpha_rename_integral` |

`main_expression()` lives in `operators.py`. It is the structural source for
the positive four-term correction. Ordered application occurs in `actions.py`.
Nested-integral kernel extraction and the separate scalar construction of
`M12`, `M21`, and `C22` occur in `kernels.py`.

The final compact action is built by:

```python
apply_combined_kernel_c22(x, f, y)
```

which returns `Integral(combined_kernel_c22(x, y) * f(y), (y, 0, oo))`.

## 3. Traceability chain

| Step | Input API | Output | Automatic? | Manual information / duplication |
|---|---|---|---|---|
| 1 | `a22_first_schur_reduction()` | `FirstSchurReduction` | Yes | Exact block indices and exterior `-1` are encoded by the m=2 constructor. |
| 2 | `a22_first_schur_mod_compact_relation()` | `ModCompactSchurRelation` | Yes | Calls the exact and model constructors; no implicit rewrite. |
| 3 | `a22_first_schur_model()` | Eight-term `LinearCombination` | Yes | Concatenates `a22_exact_operator().terms` and correction terms. |
| 4 | `a22_first_schur_correction()` | Four-term `LinearCombination` | Yes | Derives a scalar sign and scales `main_expression()`. |
| 5 | `main_expression()` | Four ordered terms | Yes | This is the structural source of products and base signs. |
| 6 | `expand_ordered(...)` | Four ordered terms | Yes but redundant for this already-expanded expression | No additional formula. |
| 7 | `apply_linear_combination_ordered(...)` | `AppliedLinearCombination` | Yes | Caller supplies operand, external variable, and formal or explicit rules. |
| 8 | The applied terms | Three nested integrals per term | Yes | Dummy variables are selected by the action engine. |
| 9 | `extract_applied_kernels(...)` | Four-term `KernelCombination` | Yes | Caller supplies `x`, function class `f`, and free kernel input `y`. |
| 10 | `m21_kernel(x, u)` | Scalar two-term `M21` | No connection from step 9 | Formula is manually encoded in `kernels.py`. |
| 11 | `m12_kernel(v, y)` | Scalar two-term `M12` | No connection from step 9 | Formula is manually encoded in `kernels.py`. |
| 12 | `combined_kernel_c22(x, y)` | Factorized nested integral | Automatic only from the manual M constructors | Does not consume `KernelCombination` or an applied correction. |
| 13 | `apply_combined_kernel_c22(x, f, y)` | Final input integral | Yes from step 12 | Caller again supplies variables and function class. |

The chain is real and reusable from K4B through kernel extraction (steps
1-9). It is broken as an API chain between steps 9 and 10-12. The H/I scalar
factorization is mathematically consistent with the extracted kernels, but it
is a parallel construction rather than a derivation from the K4B object.

## 4. `main_expression()` audit

The implementation is:

```python
correction_sign * main_expression()
```

The probe produced:

```text
a22_first_schur_correction() == main_expression()  True
a22_first_schur_correction() is main_expression()  False
```

There is structural equality, not object identity. `main_expression()` creates
a new immutable value on each call, and scalar multiplication creates another
`LinearCombination`. Products and coefficients originate directly from
`main_expression`; K4B does not spell out any of its four products again.

The resulting coefficients are `(1, -1, -1, 1)` and factor order is preserved.

## 5. Correction sign audit

`a22_first_schur_correction()` computes:

```python
schur_sign = a22_first_schur_reduction().offdiagonal_product_coefficient  # -1
leading_a21_sign = a21_wh_model().expression.terms[0].coefficient         # -1
correction_sign = schur_sign * leading_a21_sign                           # +1
return correction_sign * main_expression()
```

This is **case A**: the derived sign controls the returned expression because
it is used as the scalar multiplier. It is not merely asserted or discarded.
The positive signs already inside `main_expression()` are not a second K4B
expansion; they define the positive `B21 R11 B12` source which is then scaled.

There remains a small coupling assumption: the first coefficient of
`a21_wh_model()` is treated as the global leading sign. This is valid for its
current normalized two-term shape and is protected by K1/K4B tests.

## 6. Origin of `M21` and `M12`

The real formulas in `kernels.py` are manually encoded:

```text
M12(v,y) = rho1(v) Lplus_12(gamma1*v,y) - G1(v) Lplus_12(v,y)
M21(x,u) = rho2(x) Lminus_21(gamma2*x,u) - G2(x) Lminus_21(x,u)
```

They are not constructed from `a12_wh_model`, `a21_wh_model`,
`main_expression`, `AppliedLinearCombination`, or `KernelCombination`.
Consequently their minus signs are independent constants in scalar formulas.
The tests prove algebraic agreement with the four action-derived kernels, but
the AST is not their source at runtime.

This is a genuine second source of truth for factor formulas and signs.

## 7. Origin and shape of `C22`

`combined_kernel_c22()` is classification **D**: it manually composes the
separately constructed `M21` and `M12` around the formal `R11` kernel. It does
not consume the correction, its action, or its extracted kernels.

Its real structure is:

```text
Integral(
    M21(x,u) * Integral(
        R11(u,v) * M12(v,y),
        (v,0,oo)
    ),
    (u,0,oo)
)
```

Thus the scalar order is `M21`, `R11`, `M12`; `v` is integrated inside and
`u` outside. With ordinary variables `x,y`, defaults are `u,v`. The formal
formula has no spurious `2`, `-2`, or `1/2`.

The action-derived route uses the action engine's fresh-symbol policy and has
tests for capture avoidance. In contrast, `combined_kernel_c22()` defaults
directly to `Symbol("u")` and `Symbol("v")`; it does not call `fresh_symbol`
and does not validate collisions with caller-supplied `x` or `y`. Callers can
pass explicit `outer_variable` and `middle_variable`, but freshness is manual.
This is a real alpha-renaming/capture robustness gap in the compact constructor.

The extracted formal kernel sum is not structurally equal to the factorized
`combined_kernel_c22()` object. Phase I tests expand the scalar product and
prove termwise mathematical agreement instead. No production API performs
that factorization.

## 8. Formal versus explicit Fourier pipeline

### Formal route

The following works:

```python
correction = a22_first_schur_correction()
applied = apply_linear_combination_ordered(
    correction, f(x), x, mvp_atomic_rules()
)
kernels = extract_applied_kernels(applied, x, f, y)
```

It yields four kernels, coefficients `(1,-1,-1,1)`, and bound variables
`u,v`. The separate `combined_kernel_c22()` yields the equivalent formal
factorized formula, but there is no automatic transition between them.

### Explicit route

Replacing the rules with:

```python
explicit_wiener_hopf_rules(decay=d)
```

also applies and extracts all four correction kernels. Each extracted kernel
contains `chi_infinity` and normalized Fourier denominators, contains no formal
`Lplus_12`/`Lminus_21`, and preserves bound `u,v`. For example, its first term
has the structure

```text
Integral(
  chi(u) chi(gamma2*x) rho2(x) /
    (2*pi*(d + i*(gamma2*x-u)))
  * Integral(
      R11(u,v) chi(y) chi(gamma1*v) rho1(v) /
        (2*pi*(d - i*(gamma1*v-y))),
      (v,0,oo)),
  (u,0,oo))
```

The other three terms replace shifted multipliers according to the same AST.

However, `m12_kernel`, `m21_kernel`, and `combined_kernel_c22` hard-code formal
`Lplus_12` and `Lminus_21` and accept no rule mapping. Therefore there is no
compact explicit `C22` end-to-end. Classification for the requested complete
chain is **neither mode is connected end-to-end**: formal has parallel compact
construction; explicit stops after four extracted kernels.

The exact break is `KernelCombination -> M21/M12/combined_kernel_c22`.

## 9. Complete reduced-block action

`a22_first_schur_model()` can already be applied as eight ordered terms. A
probe confirmed:

```text
AppliedTerm count: 8
PrincipalValue present: yes
Integral atoms in SymPy projection: 10
free symbols: gamma1, gamma2, x
bound symbols: u, v, y
```

There is no API that returns the semantically compact form directly. It can be
constructed without duplicating formulas by combining existing APIs:

```python
apply_linear_combination(a22_exact_operator(), f(x), x, rules)
+ apply_combined_kernel_c22(x, f, y, outer_variable=u, middle_variable=v)
```

This is classification **B**: existing APIs suffice, although no named helper
encapsulates the composition. Such a helper would be small, but creating it is
outside this audit.

Calling `extract_applied_kernels` on the complete eight-term applied model
fails at the identity contribution with `KernelExtractionError`, because that
API intentionally expects every term to have a final input integral. Kernel
extraction must be restricted to the correction.

## 10. `PrincipalValue` coexistence

The ordered `AppliedLinearCombination` safely holds both diagonal
`PrincipalValue` terms and ordinary nested correction integrals. The SymPy
projection also holds both. Tests cover shift substitution, fresh bound
variables, free symbols, alpha-renaming, and preservation of the inert
`PrincipalValue` wrapper.

No accidental evaluation was observed. `AppliedLinearCombination` preserves
conceptual order; conversion with `as_expr()` produces a SymPy `Add`, whose
display/order is controlled by SymPy and must not be used as the ordered source
of truth.

The compact probe (`A22` action plus `apply_combined_kernel_c22`) retained
`PrincipalValue`, had five integral atoms, free symbols `gamma1,gamma2,x`, and
bound symbols `u,v,y`. The only material variable-safety caveat is the compact
kernel constructor's manual default `u,v` selection described above.

## 11. Semantic relation types

| Type | Semantics | Implicit conversion or equality behavior |
|---|---|---|
| `ExactBlock` | Reference to an exact indexed block | None |
| `WienerHopfModel` | Normalized AST model for one exact block | None |
| `ModCompactRelation` | Specific exact-block/WH-model equivalence modulo compacta | None |
| `FormalRegularizer` | Associates a formal atom with its exact target block | None; no left/right orientation |
| `FirstSchurReduction` | Exact metadata for `A22 - A21 R11 A12` | None; exact blocks remain non-applicable |
| `ModCompactSchurRelation` | Specific reduced-block/model equivalence modulo compacta | None |

The types have distinct semantics. The two modulo-compact relation classes
share a pattern but deliberately validate different exact/model types; this is
not harmful duplication at the current small scope. Neither behaves as exact
equality, exposes `rewrite`, nor converts sides implicitly. They should remain
specific until a broader use case justifies generalization.

## 12. Public API surface

The current public surface, grouped by role, is:

- **AST:** `OperatorAtom`, `Product`, `Term`, `LinearCombination`,
  `expand_ordered`, `main_expression`, `G1`, `G2`, `I`, `S_Rplus`, `R11`,
  `Vtilde_alpha1`, `Vtilde_alpha2`, `Wplus_12`, `Wminus_21`.
- **Actions:** `AtomicAction`, `AppliedTerm`, `AppliedLinearCombination`,
  `IdentityAction`, `MultiplicationAction`, `TransportedShiftAction`,
  `IntegralKernelAction`, `FormalRegularizerAction`,
  `PrincipalValueIntegralAction`, `PrincipalValue`, `mvp_atomic_rules`,
  `explicit_wiener_hopf_rules`, `apply_atom`, `apply_product`,
  `apply_linear_combination_ordered`, `apply_linear_combination`, and the
  action error types.
- **Fourier:** `InverseFourierIntegral`, `FourierEvaluationError`,
  `frequency_variable`, `time_variable`, `positive_decay_symbol`,
  `forward_exponential`, `inverse_exponential`, `inverse_prefactor`,
  `bplus_symbol`, `bminus_symbol`, `positive_inverse_integral`,
  `negative_inverse_integral`, `kplus_kernel`, `kminus_kernel`,
  `localized_lplus_kernel`, `localized_lminus_kernel`.
- **Blocks:** `pplus_operator`, `pminus_operator`, `a11_exact_operator`,
  `a22_exact_operator`, `a12_wh_model`, `a21_wh_model`,
  `a12_mod_compact_relation`, `a21_mod_compact_relation`,
  `a11_formal_regularizer`.
- **Relations:** `ExactBlock`, `WienerHopfModel`, `ModCompactRelation`,
  `FormalRegularizer`, `FirstSchurReduction`, `ModCompactSchurRelation`.
- **Kernels:** `KernelTerm`, `KernelCombination`, `KernelExtractionError`,
  `extract_integral_kernel`, `extract_applied_kernels`, `m12_kernel`,
  `m21_kernel`, `c22_integrand`, `combined_kernel_c22`,
  `apply_combined_kernel_c22`.
- **Schur:** `a22_first_schur_reduction`,
  `a22_first_schur_correction`, `a22_first_schur_model`,
  `a22_first_schur_mod_compact_relation`.
- **Substitution utilities:** `SafeSubstitutionError`, `fresh_symbol`,
  `collect_bound_symbols`, `substitute_free_variable`.

No circular export or import was found. Naming is mostly explicit. The notable
ambiguities are historical: `main_expression` does not reveal that it is the
positive first Schur correction, and `combined_kernel_c22` does not expose in
its name whether it is formal or explicit. `m12_kernel`/`m21_kernel` overlap
conceptually with action-derived kernels but are independent implementations.

## 13. Dependency audit

Direct internal imports are:

```text
operators       -> (none)
fourier         -> (none)
substitution    -> (none)
actions         -> fourier, operators, substitution
relations       -> operators
blocks          -> operators, relations
kernels         -> actions, operators
package __init__-> actions, blocks, fourier, kernels, operators, relations,
                  substitution
```

There are no cycles, deferred imports, or runtime import tricks. Low-level
`operators.py` and semantic `relations.py` do not import SymPy. `kernels.py`
depends on the action result type, which matches its extraction role. The
`actions -> fourier` dependency was introduced intentionally so explicit rules
reuse Fourier kernels; it does not create a reverse dependency.

No module of lower structural level imports `blocks.py` or Schur metadata.

## 14. Test coverage by phase

| Phase | Main tests | Protected invariants |
|---|---|---|
| H | `test_full_application.py`, `test_noncommutative_expansion.py` | Four products, signs, action order, shifts, three nested integrals, ordered results |
| I | `test_kernels.py` | Kernel extraction, bound/free variables, M formulas, factorized C22, final compact action |
| J | `test_fourier.py` | Fourier convention, normalized kernels, localization, original-vs-normalized factors |
| K1 | `test_exact_coefficients.py`, `test_relations.py`, `test_offdiagonal_blocks.py` | `Fraction`, exact/model distinction, normalized A12/A21 coefficients |
| K2 | `test_explicit_wiener_hopf.py` | Formal default, explicit normalized actions, shifted cutoffs, full A12/A21 actions |
| K3A | `test_diagonal_blocks.py` | Identity, structural S_Rplus, exact projections, exact A11/A22 products |
| K3B | `test_principal_value_action.py` | Principal value, singular denominator, shifts, freshness, complete diagonal actions |
| K4A | `test_formal_regularizer.py` | Formal metadata, unchanged R11 action, no inverse/cancellation/orientation |
| K4B | `test_first_schur_reduction.py` | Exact reduction, sign derivation, reuse of `main_expression`, eight-term model, relation and application |

Tests overwhelmingly compare AST or SymPy structure. No mathematical invariant
depends on fragile string equality. A few tests parse source files with `ast`
to enforce the architectural absence of SymPy imports or float literals; these
are appropriate source-level constraints, not string representations of math.

There is no test starting specifically from
`a22_first_schur_correction()` and ending at `combined_kernel_c22()`. H/I start
from `main_expression`; K4B proves the correction equals that expression and
can be applied, but the two facts are not composed into one end-to-end test.
There is also no explicit-to-compact-C22 test because no such production path
exists.

## 15. Confirmed technical debt

1. **Primary disconnect:** `KernelCombination` cannot be factored or converted
   into `M21 R11 M12`; `combined_kernel_c22` is parallel to the action-derived
   route.
2. **Duplicated formulas/signs:** `m12_kernel` and `m21_kernel` manually repeat
   multipliers, shifts, formal kernels, and minus signs already implied by the
   AST/actions.
3. **Formal-only compact constructor:** `combined_kernel_c22` cannot consume
   `explicit_wiener_hopf_rules`, although the four-term extraction route can.
4. **Freshness gap:** compact C22 defaults to fixed `u,v` rather than using the
   action engine's collision avoidance.
5. **No composed regression:** no test crosses the K4B constructor boundary
   through action/extraction to compact C22.
6. **Historical naming:** `main_expression` is less specific than its current
   mathematical role, though renaming is not justified within this audit.

`correction_sign` is not dead code and is not technical debt of case C.

## 16. Required recommendation for K5B

**OPTION A: K5B should connect semantically
`a22_first_schur_correction()` to `C22`.**

The next phase should make the action-derived correction (or its extracted
`KernelCombination`) the source consumed by the compact combined-kernel path,
while preserving ordered AST provenance and supporting formal and explicit
rules. It should also centralize fresh dummy selection. It should not begin
recursive Schur reduction, general matrices, or arbitrary `m`.

This recommendation is preferred over option B because the existing pieces are
mathematically consistent and well tested; the defect is the missing semantic
connection itself, not a prerequisite mathematical error. Option C is false
because the pipeline is not fully connected, and option D is unsupported by
the audit.
