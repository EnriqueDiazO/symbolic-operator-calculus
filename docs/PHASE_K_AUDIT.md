# Phase K0 Audit: Aij Pre-Implementation

This report audits the current repository before implementing the blocks
`A_{k,j}`. It does not implement new operators, actions, relations, or tests.

## 1. Executive Summary

The project currently has a solid m=2 MVP pipeline:

- a pure structural noncommutative operator AST;
- atomic actions into SymPy scalar expressions and formal integrals;
- safe free-variable substitution under `Integral`;
- ordered application of the four Schur-correction terms;
- kernel extraction for those applied terms;
- formal `M12`, `M21`, `C22`, and `C22 f`;
- evaluated normalized Fourier kernels `Kplus` and `Kminus`.

The most important blocker for the true `A_{k,j}` phase is exact scalar
support. The AST accepts `int | float | complex` coefficients but rejects both
`Fraction(1, 2)` and `sympy.Rational(1, 2)`. Therefore it cannot represent
`1/2 (Vtilde_alpha_k - G_k) T_{k,j}` exactly today. It can only use the
inexact float `0.5`.

The second important blocker is mathematical relation typing. The code has no
`ExactEquality`, no `ModCompactRelation`, and no equivalent type. Calling the
normalized Wiener-Hopf model simply `A12` or `A21` would risk confusing an
exact block with an equivalence modulo compact operators.

The off-diagonal normalized m=2 Wiener-Hopf models are largely available today
if written as already-expanded `LinearCombination` objects:

- `A12_WH = Vtilde_alpha1 * Wplus_12 - G1 * Wplus_12`;
- `A21_WH = -Vtilde_alpha2 * Wminus_21 + G2 * Wminus_21`.

They apply correctly to `f(x)` using the current engine. They use formal
`Lplus_12` and `Lminus_21` by default. Explicit Fourier kernels can be used
without modifying code by passing a local `rules` mapping with replacement
`IntegralKernelAction` callables.

The diagonal blocks `A11` and `A22` are not ready: there are no atoms or
actions for identity, `S_Rplus`, `Pplus`, or `Pminus`, and exact `1/2` is not
available.

Recommended path: implement the next real phase in small steps. K1 should stay
with m=2, add exact coefficient support and relation wrappers, and expose
off-diagonal Wiener-Hopf models with safe names. K2 should connect formal
Wiener-Hopf kernels to explicit Fourier kernels through a small kernel-rule
API. K3 should add diagonal projection machinery.

## 2. Git State

Commands requested:

```text
pwd
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git rev-parse --show-toplevel
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git branch --show-current
main

git status --short
<empty>

git log -10 --oneline --decorate
a8588a7 (HEAD -> main, origin/main) feat: evaluate normalized Wiener-Hopf Fourier kernels
6cc07f0 feat: derive Schur correction integrals and combined kernel
d1cae8f feat: support safe integral composition and variable substitution
579450d feat: add symbolic operator actions and scalar algebra
91fb1ba feat: establish noncommutative operator AST

git rev-parse --abbrev-ref --symbolic-full-name @{u}
origin/main

git remote -v
origin git@github-personal:EnriqueDiazO/symbolic-operator-calculus.git (fetch)
origin git@github-personal:EnriqueDiazO/symbolic-operator-calculus.git (push)
```

Answers:

1. Last real commit: `a8588a7 feat: evaluate normalized Wiener-Hopf Fourier kernels`.
2. Included in commits: AST, actions, safe substitution, ordered application,
   kernel extraction, `M12/M21/C22`, and Fourier kernels.
3. Phases only as local changes: none at audit start.
4. H/I/J files are committed. H/I are in `6cc07f0`; J is in `a8588a7`.
5. Initial tree state: clean.
6. Upstream exists: `origin/main`.
7. Origin exists: `git@github-personal:EnriqueDiazO/symbolic-operator-calculus.git`.
8. No push was performed.

| Phase | State | Commit | Current files |
|---|---|---|---|
| AST / early operator algebra | Committed | `91fb1ba` | `operators.py`, `tests/test_operators.py`, `tests/test_noncommutative_expansion.py` |
| Actions and scalar algebra | Committed | `579450d` | `actions.py`, `tests/test_actions.py` |
| Safe integral composition / substitution | Committed | `d1cae8f` | `substitution.py`, action updates, substitution tests |
| H: ordered full application | Committed | `6cc07f0` | `AppliedTerm`, `AppliedLinearCombination`, `tests/test_full_application.py` |
| I: extracted kernels and C22 | Committed | `6cc07f0` | `kernels.py`, `tests/test_kernels.py` |
| J: Fourier Kplus/Kminus | Committed | `a8588a7` | `fourier.py`, `tests/test_fourier.py` |
| K0 audit | Local report only | none | `docs/PHASE_K_AUDIT.md` |

## 3. Architecture Map

High-level flow:

```text
OperatorAtom
     |
     v
Product / LinearCombination
     |
     v
AtomicAction mapping
     |
     v
SymPy Expr / Integral
     |
     v
AppliedLinearCombination
     |
     v
KernelCombination
     |
     v
M12, M21, C22

Fourier layer:
bplus/bminus -> unevaluated inverse Fourier Integral -> evaluated Kplus/Kminus
             -> localized Lplus/Lminus explicit kernels
```

Connections:

- `operators.py` is pure Python and does not import SymPy.
- `actions.py` bridges the AST to SymPy scalar expressions and formal
  integrals.
- `substitution.py` protects free-variable substitutions under SymPy
  `Integral`.
- `kernels.py` extracts formal kernels from applied expressions and builds
  formal `M12/M21/C22`.
- `fourier.py` computes explicit normalized Wiener-Hopf kernels, but it is not
  wired into `actions.py` by default.
- `__init__.py` exports the public MVP API.

## 4. Public API Inventory

### `operators.py`

Responsibility: structural noncommutative operator language.

Public classes:

- `OperatorAtom(name: str, kind: str = "operator")`
- `Product(factors: tuple[OperatorAtom, ...])`
- `Term(coefficient: Scalar, product: Product)`
- `LinearCombination(terms: tuple[Term, ...])`

Public functions:

- `expand_ordered(expression: Factor) -> LinearCombination`
- `main_expression() -> LinearCombination`

Important internals:

- `_is_scalar(value)`
- `_checked_coefficient(coefficient)`
- `_unsupported_distribution()`

SymPy dependency: none.

Immutable structures: all main classes are `@dataclass(frozen=True)`.

Important assumptions:

- all operator atoms are noncommutative;
- products contain only `OperatorAtom`;
- `LinearCombination` terms remain ordered;
- distributivity across `LinearCombination` is intentionally not implemented.

### `actions.py`

Responsibility: apply operator atoms/products/combinations to SymPy scalar
expressions.

Public classes:

- `AppliedTerm(coefficient, product, expression)`
- `AppliedLinearCombination(result_terms)`
- `MultiplicationAction(coefficient)`
- `TransportedShiftAction(rho, gamma)`
- `IntegralKernelAction(kernel)`
- `FormalRegularizerAction(kernel)`
- error classes: `AtomicActionError`, `MissingAtomicActionError`,
  `UnsupportedAtomicExpressionError`, `UnsupportedProductExpressionError`,
  `UnsupportedLinearCombinationExpressionError`, `ProductApplicationError`,
  `LinearCombinationApplicationError`, `MissingIntegrationVariableError`,
  `UnsafeScalarSubstitutionError`

Public functions:

- `mvp_atomic_rules() -> Mapping[OperatorAtom, AtomicAction]`
- `apply_atom(atom, operand, variable, rules, *, integration_variable=None)`
- `apply_product(product, operand, variable, rules, *, integration_variable=None)`
- `apply_linear_combination(combination, operand, variable, rules, *, integration_variable=None)`
- `apply_linear_combination_ordered(combination, operand, variable, rules, *, integration_variable=None)`

Important internals:

- `_kernel_integral(kernel, operand, variable, integration_variable)`
- `_replace_free_variable(operand, variable, replacement)`
- `_atom_needs_integration_variable(atom, rules)`
- `_fresh_integration_variable(expression, output_variable, generated_variables)`

SymPy dependency: central; actions produce `sp.Expr` and `sp.Integral`.

Immutable structures: action classes and applied-result classes are frozen
dataclasses.

Important assumptions:

- `Product` is applied right-to-left;
- `Product(())` acts like identity on an operand;
- integral atoms require an integration variable unless in a multi-factor
  product where local dummy variables are generated;
- the default rule mapping is hard-coded but passed explicitly to application
  functions.

### `substitution.py`

Responsibility: safe substitution of free SymPy variables under integrals.

Public classes:

- `SafeSubstitutionError`

Public functions:

- `collect_bound_symbols(expression: sp.Expr) -> set[sp.Symbol]`
- `fresh_symbol(avoid, *, preferred_names=("y", "u", "v", "w", "z")) -> sp.Symbol`
- `substitute_free_variable(expression, old_variable, new_expression) -> sp.Expr`
- `alpha_rename_integral(integral, rename_map) -> sp.Integral`

Important internals:

- `_substitute_free`
- `_substitute_free_in_integral`
- `_substitute_free_in_limits`
- `_rename_bound_occurrences`

SymPy dependency: direct and focused on `Integral`, symbols, and expression
trees.

Immutable structures: no dataclasses; functions return new SymPy expressions.

Important assumptions:

- only the MVP integral patterns are supported;
- alpha-renaming is local to integrals and is not a general alpha-equivalence
  engine for all SymPy constructs.

### `kernels.py`

Responsibility: extract kernels from applied expressions and build formal
Schur-correction kernels.

Public classes:

- `KernelExtractionError`
- `KernelTerm(coefficient, product, kernel)`
- `KernelCombination(terms)`

Public functions:

- `extract_integral_kernel(applied_expression, output_variable, input_function, input_variable)`
- `extract_applied_kernels(applied_combination, output_variable, input_function, input_variable)`
- `m12_kernel(left_variable, input_variable)`
- `m21_kernel(output_variable, right_variable)`
- `c22_integrand(output_variable, outer_variable, middle_variable, input_variable)`
- `combined_kernel_c22(output_variable, input_variable, *, outer_variable=None, middle_variable=None)`
- `apply_combined_kernel_c22(output_variable, input_function, input_variable, *, outer_variable=None, middle_variable=None)`

Important internals:

- `_remove_input_integral`
- `_remove_exact_input_factor`
- `_canonicalize_remaining_integrals`
- `_split_single_integral_factor`

SymPy dependency: direct; builds and inspects `sp.Integral`, `sp.Mul`, and
symbolic functions.

Immutable structures: `KernelTerm` and `KernelCombination` are frozen
dataclasses.

Important assumptions:

- extraction is for MVP nested integral shapes;
- `C22` uses formal `Lplus_12`, `Lminus_21`, `R11`, `rho`, `G`, and `gamma`;
- no relation modulo compact operators is represented.

### `fourier.py`

Responsibility: explicit Fourier convention and evaluated normalized kernels.

Public classes:

- `FourierEvaluationError`
- `InverseFourierIntegral(symbol_expression, integrand, integration_variable, lower_limit, upper_limit, decay_parameter, time_variable)`

Public functions:

- `frequency_variable() -> sp.Symbol`
- `time_variable() -> sp.Symbol`
- `positive_decay_symbol() -> sp.Symbol`
- `forward_exponential(time=None, frequency=None) -> sp.Expr`
- `inverse_exponential(time=None, frequency=None) -> sp.Expr`
- `inverse_prefactor() -> sp.Expr`
- `bplus_symbol(decay=None, frequency=None) -> sp.Expr`
- `bminus_symbol(decay=None, frequency=None) -> sp.Expr`
- `positive_inverse_integral(time=None, decay=None, frequency=None) -> InverseFourierIntegral`
- `negative_inverse_integral(time=None, decay=None, frequency=None) -> InverseFourierIntegral`
- `kplus_kernel(time=None, decay=None, frequency=None) -> sp.Expr`
- `kminus_kernel(time=None, decay=None, frequency=None) -> sp.Expr`
- `localized_lplus_kernel(output_variable, input_variable, *, decay=None) -> sp.Expr`
- `localized_lminus_kernel(output_variable, input_variable, *, decay=None) -> sp.Expr`

SymPy dependency: direct; performs symbolic integration and uses assumptions.

Immutable structures: `InverseFourierIntegral` is frozen.

Important assumptions:

- Fourier convention is project-specific:
  `exp(-I*t*lambda)` forward, `exp(+I*t*lambda)` inverse, prefactor
  `1/(2*pi)`;
- `d` is real and positive;
- `evaluate()` uses `conds="none"` only after checking `d.is_real` and
  `d.is_positive`.

### `__init__.py`

Responsibility: public API re-export.

It exports operator atoms, AST classes, action classes/functions, substitution
utilities, kernel utilities, and Fourier utilities. It contains no logic beyond
imports and `__all__`.

## 5. AST Exactness Audit

Questions:

1. Atomic operator: `OperatorAtom(name, kind="operator")`, frozen dataclass.
2. Ordered product: `Product(factors=(...))`, frozen dataclass, tuple order is
   the structure.
3. Linear combination: `LinearCombination(terms=(Term(...), ...))`, frozen
   dataclass. Each `Term` stores `coefficient` and `Product`.
4. Term order is tuple order; no sorting is used.
5. Coefficients are stored on `Term.coefficient`.
6. Accepted coefficient types: `int`, `float`, `complex`, excluding `bool`.
7. Exact `1/2`: no.

Probe results:

| Coefficient | Result |
|---|---|
| `1` | accepted |
| `0.5` | accepted as `float` |
| `1+0j` | accepted as `complex` |
| `Fraction(1, 2)` | rejected: `TypeError` |
| `sympy.Rational(1, 2)` | rejected: `TypeError` |
| `True` | rejected: `TypeError` |

Additional scalar probes:

- `Fraction(1, 2) * (Vtilde_alpha1 - G1)` is rejected.
- `sympy.Rational(1, 2) * (Vtilde_alpha1 - G1)` is rejected.
- `0.5 * (Vtilde_alpha1 - G1)` is accepted, but inexact.

Conclusion: the current AST cannot represent
`1/2 (Vtilde_alpha_k - G_k) T_{k,j}` exactly without code changes.

## 6. Indexing Audit

Current atoms:

```text
G1
G2
Vtilde_alpha1
Vtilde_alpha2
Wplus_12
Wminus_21
R11
```

Answers:

1. They are not generic parametrized instances.
2. They are hard-coded `OperatorAtom` objects with names.
3. There is no internal concept of index `k,j`.
4. `G3` can be created today as `OperatorAtom("G3")`.
5. `Wplus_13` can be created today as `OperatorAtom("Wplus_13")`.
6. A local action can be registered by passing a custom rules mapping to
   `apply_atom`, `apply_product`, or `apply_linear_combination_ordered`.
7. Obstacles to m arbitrary:
   - no typed index data on atoms;
   - no constructors enforcing `k<j` or `k>j`;
   - default `mvp_atomic_rules()` is hard-coded for m=2;
   - no symbol registry for `c_k`, `gamma_k`, `rho_k`, `G_k`, `L_{k,j}`;
   - no general sign/normalization policy for off-diagonal symbols;
   - no diagonal projection operators.

Probe:

```text
G3 = OperatorAtom("G3")
Wplus_13 = OperatorAtom("Wplus_13")
G3 * Wplus_13 -> Product((G3, Wplus_13))
apply_atom(G3, ..., mvp_atomic_rules()) -> MissingAtomicActionError
local_rules[G3] = MultiplicationAction(Function("G3")) -> works
local_rules[Wplus_13] = IntegralKernelAction(Function("Lplus_13")) -> works
```

## 7. Action-Rule Audit

Answers:

1. `mvp_atomic_rules()` is hard-coded.
2. It returns `MappingProxyType`, so the default mapping is immutable.
3. Another mapping can be constructed externally with `dict(mvp_atomic_rules())`.
4. `apply_product` accepts alternative mappings via its `rules` argument.
5. `Lplus_12` can be replaced by an explicit kernel without modifying
   `actions.py` by replacing the action for `Wplus_12` in a local rules map.
6. Yes, this can be done with a local rule:
   `IntegralKernelAction(lambda out, inn: localized_lplus_kernel(out, inn, decay=d))`.
7. `IntegralKernelAction` expects a callable with signature like
   `Callable[..., sp.Expr]`; in practice it is called as `kernel(variable,
   integration_variable)`.

Diagram:

```text
OperatorAtom
     |
     v
AtomicAction
     |
     v
SymPy expression / Integral
```

Important behavior:

- `apply_atom` applies exactly one atom.
- `apply_product` requires a `Product`, applies right-to-left, and generates
  local dummy variables for multi-factor products containing integral actions.
- `apply_linear_combination_ordered` preserves ordered applied terms and
  coefficients.
- `apply_linear_combination` projects the ordered result to a SymPy sum.

Limitation:

```text
(Vtilde_alpha1 - G1) * Wplus_12
```

raises:

```text
TypeError: Products involving LinearCombination require distributivity and are not supported in this phase.
```

Therefore A12/A21 must be built today as already-expanded
`LinearCombination` objects.

## 8. Kernel/Fourier Integration Audit

Actions layer:

- `actions.py` uses formal undefined SymPy functions:
  `Lplus_12(x, y)` and `Lminus_21(x, y)`.

Fourier layer:

- `fourier.py` computes explicit expressions:
  `chi_infinity(x) * Kplus(x-y) * chi_infinity(y)`;
  `chi_infinity(x) * Kminus(x-y) * chi_infinity(y)`.

Probe:

```text
actions.Lplus_12_kernel: Lplus_12, type UndefinedFunction
localized_lplus_kernel(x,y): chi_infinity(x)*chi_infinity(y)/(2*pi*d*(1 - I*(x - y)/d))
formal action: Integral(Lplus_12(x, y)*f(y), (y, 0, oo))
formal Lplus_12(x,y) == explicit localized Lplus: False
```

Answers:

1. `Lplus_12` in `actions.py` and `localized_lplus_kernel` in `fourier.py`
   are not the same mathematical object in the code.
2. They do not have the same type: one is an applied undefined SymPy function,
   the other is a SymPy `Mul` expression.
3. They do not share symbolic identity.
4. They live in disconnected layers.
5. They can be connected manually through a local action rule, or by manual
   expression replacement after applying a formal action.
6. There is no public kernel-substitution API.
7. Minimal change to apply A12 using Fourier explicit kernels: add a small
   kernel rule/registry API or expose helper rules such as
   `explicit_wiener_hopf_rules(decay=d)` that maps `Wplus_12` and `Wminus_21`
   to `IntegralKernelAction(localized_lplus/lminus)`.

Manual replacement probe:

```text
Integral(Lplus_12(x, y)*f(y), (y, 0, oo))
-> Integral(chi_infinity(x)*chi_infinity(y)*f(y)/(2*pi*d*(1 - I*(x - y)/d)), (y, 0, oo))
```

This works as an ad hoc probe but does not preserve formal kernel provenance.

## 9. Normalization and Sign Audit

Reference normalization:

For `k<j`:

```text
b_{k,j} = 2 chi_+ exp((c_k-c_j) lambda)
```

For m=2, `d=c2-c1>0`, so

```text
b_{1,2} = 2 chi_+ exp(-d lambda)
1/2 (Vtilde_alpha1 - G1) W(b_{1,2})
  -> (Vtilde_alpha1 - G1) Wplus_12
```

For `k>j`:

```text
b_{k,j} = -2 chi_- exp((c_k-c_j) lambda)
```

Thus

```text
b_{2,1} = -2 chi_- exp(d lambda)
1/2 (Vtilde_alpha2 - G2) W(b_{2,1})
  -> -(Vtilde_alpha2 - G2) Wminus_21
```

| Block | Original symbol | Exterior factor | Normalized block | Final sign |
|---|---|---:|---|---:|
| `A12` | `2 chi_+ exp(-d lambda)` | `1/2` | `Wplus_12 = chi_+ exp(-d lambda)` | `+` |
| `A21` | `-2 chi_- exp(d lambda)` | `1/2` | `Wminus_21 = chi_- exp(d lambda)` | `-` |

Duplication risks:

- Adding `1/2` again when using normalized `Wplus_12` would incorrectly scale
  `A12` by `1/2`.
- Keeping the original factor `2` inside a future `Wplus_12` and also using the
  normalized formula would incorrectly scale by `2`.
- For `A21`, the `-2` from the original symbol has already become the exterior
  minus in the normalized formula. A future implementation must not both put
  `-2` into the symbol and add the exterior minus.
- Current `fourier.bplus_symbol` and `fourier.bminus_symbol` are normalized:
  they return `exp(-d*lambda)` and `exp(d*lambda)`, with no factor `2` or `-2`.

## 10. A12/A21 Probes

Direct construction probe:

```text
(Vtilde_alpha1 - G1) * Wplus_12
-> TypeError: Products involving LinearCombination require distributivity...
```

Expanded construction works:

```python
A12_WH = LinearCombination((
    Term(1, Vtilde_alpha1 * Wplus_12),
    Term(-1, G1 * Wplus_12),
))

A21_WH = LinearCombination((
    Term(-1, Vtilde_alpha2 * Wminus_21),
    Term(1, G2 * Wminus_21),
))
```

Formal A12 result:

```text
coefficients: (1, -1)
term 1: rho1(x)*Integral(Lplus_12(gamma1*x, y)*f(y), (y, 0, oo))
term 2: G1(x)*Integral(Lplus_12(x, y)*f(y), (y, 0, oo))
sum:
-G1(x)*Integral(Lplus_12(x, y)*f(y), (y, 0, oo))
+rho1(x)*Integral(Lplus_12(gamma1*x, y)*f(y), (y, 0, oo))
integrals in summed expression: 2
```

Formal A21 result:

```text
coefficients: (-1, 1)
term 1: rho2(x)*Integral(Lminus_21(gamma2*x, y)*f(y), (y, 0, oo))
term 2: G2(x)*Integral(Lminus_21(x, y)*f(y), (y, 0, oo))
sum:
G2(x)*Integral(Lminus_21(x, y)*f(y), (y, 0, oo))
-rho2(x)*Integral(Lminus_21(gamma2*x, y)*f(y), (y, 0, oo))
integrals in summed expression: 2
```

What is preserved:

- ordered `AppliedTerm` entries;
- coefficients;
- originating `Product`;
- formal SymPy integral structure.

What is lost or absent:

- no relation marker saying this is modulo compact operators;
- no exact block type;
- no direct unexpanded product `(V-G)W`;
- formal `Lplus_12`/`Lminus_21` do not remember Fourier derivation.

## 11. Explicit-Kernel Probes

Without permanent code changes, explicit kernels work via local rules:

```python
d = positive_decay_symbol()
explicit_rules = dict(mvp_atomic_rules())
explicit_rules[Wplus_12] = IntegralKernelAction(
    lambda out, inn: localized_lplus_kernel(out, inn, decay=d)
)
explicit_rules[Wminus_21] = IntegralKernelAction(
    lambda out, inn: localized_lminus_kernel(out, inn, decay=d)
)
```

Explicit A12 result:

```text
term 1:
rho1(x)*Integral(
  chi_infinity(y)*chi_infinity(gamma1*x)*f(y)
  /(2*pi*d*(1 - I*(gamma1*x - y)/d)),
  (y, 0, oo)
)

term 2:
G1(x)*Integral(
  chi_infinity(x)*chi_infinity(y)*f(y)
  /(2*pi*d*(1 - I*(x - y)/d)),
  (y, 0, oo)
)

sum:
- term2 + term1
```

Explicit A21 result:

```text
term 1:
rho2(x)*Integral(
  chi_infinity(y)*chi_infinity(gamma2*x)*f(y)
  /(2*pi*d*(1 + I*(gamma2*x - y)/d)),
  (y, 0, oo)
)

term 2:
G2(x)*Integral(
  chi_infinity(x)*chi_infinity(y)*f(y)
  /(2*pi*d*(1 + I*(x - y)/d)),
  (y, 0, oo)
)

sum:
term2 - term1
```

Conclusion:

- It works using local rules.
- It requires manual rule construction today.
- A small public API would make this reliable.
- Ordered applied-term structure is preserved.
- There is still no explicit relation object tying the exact block to the
  Wiener-Hopf model.

## 12. Diagonal-Block Gap

Target:

```text
Akk = Vtilde_alpha_k Pplus + G_k Pminus
Pplus = 1/2 (I + S_Rplus)
Pminus = 1/2 (I - S_Rplus)
```

Answers:

1. OperatorAtom for `I`: no predefined one.
2. OperatorAtom for `S_Rplus`: no.
3. `Pplus` and `Pminus`: no.
4. Identity action: no.
5. Cauchy singular integral action on `R+`: no.
6. Formal representation of `Pplus = 1/2(I + S_Rplus)`:
   partial only if one creates ad hoc atoms and uses inexact `0.5`;
   exact `1/2` is impossible today.
7. Application: no, because there are no actions for `I` or `S_Rplus`.
8. What blocks A11/A22:
   - no exact scalar `1/2`;
   - no projection atoms;
   - no identity action;
   - no singular integral action;
   - no diagonal block constructor;
   - no exact-vs-mod-compact relation type.

Distinctions:

| Layer | Current state |
|---|---|
| Formal representation | Partial with ad hoc atoms and float `0.5`; not exact |
| Formal application | No |
| Integral realization | No |
| Evaluation | No |

## 13. Mathematical-Relation Gap

Answers:

1. `ExactEquality`: no.
2. `ModCompactRelation`: no.
3. Equivalent type: no.
4. The current code could confuse exact equality and equivalence modulo compact
   operators if future constructors use names like `A12` directly for WH
   models.
5. Yes, it would be dangerous to call the WH model simply `A12`.
6. Recommended internal names:
   - `ExactBlock` for exact formulas;
   - `WienerHopfModel` for normalized WH approximants;
   - `ModCompactRelation(exact, model)` for `exact ~= model`.

Recommendation:

- Avoid using `A12` as a single overloaded object.
- Prefer `A12_exact` or `ExactBlock("A", 1, 2)` for the exact block.
- Prefer `A12_wh_model` or `WienerHopfModel(1, 2, ...)` for the normalized WH
  model.
- Connect them with an explicit `ModCompactRelation`.

## 14. Generalization Strategies

### Strategy A: keep MVP exclusively m=2

Changes needed:

- Add only m=2 block helpers: `A12_WH`, `A21_WH`, later `A11`, `A22`.
- Add exact coefficient support.
- Add relation wrappers.

Risk: low.

Complexity: low.

Benefit for current article: high if the immediate target is only m=2.

Benefit for m>=3: low.

### Strategy B: create a general `A(k,j)` constructor but keep concrete m=2 actions

Changes needed:

- Introduce indexed block metadata and constructors.
- Keep existing atoms/actions for `G1`, `G2`, `Vtilde_alpha1`,
  `Vtilde_alpha2`, `Wplus_12`, `Wminus_21`.
- Add normalization/sign policy for known m=2 off-diagonal blocks.
- Add exact coefficient support and relation wrappers.

Risk: moderate.

Complexity: moderate.

Benefit for current article: high.

Benefit for m>=3: medium; the API can grow without forcing full symbolic
generalization now.

### Strategy C: generalize operators, kernels, and symbols to arbitrary indices now

Changes needed:

- Indexed atom classes or factories;
- general `G_k`, `rho_k`, `gamma_k`, `W_{k,j}`, `b_{k,j}`, `c_k`;
- sign and support logic for all `k,j`;
- generalized action registry;
- likely diagonal projection support;
- relation model.

Risk: high.

Complexity: high.

Benefit for current article: uncertain; could slow the immediate work.

Benefit for m>=3: high, but premature.

Recommendation: Strategy B. It keeps m=2 computations stable while preventing
names and APIs from painting the project into a corner.

## 15. Test Matrix

| Mathematical capability | Current test evidence | State |
|---|---|---|
| Noncommutativity | `tests/test_operators.py`, `tests/test_noncommutative_expansion.py` | Covered |
| Ordered expansion | `test_main_expression_expands_to_exactly_four_terms`, factor-order tests | Covered |
| Atomic actions | `tests/test_actions.py` multiplication, shifts, kernels, regularizer | Covered |
| Product composition | `apply_product` order and integral composition tests | Covered |
| Safe substitution | `tests/test_substitution.py` | Covered |
| Alpha-renaming | substitution capture/nested tests | Covered |
| Four applied integrals | `tests/test_full_application.py` | Covered |
| Kernels K1...K4 | `tests/test_kernels.py` I1-I10 | Covered |
| M12/M21 | `tests/test_kernels.py` I11-I14 | Covered |
| C22 | `tests/test_kernels.py` I15-I18 | Covered |
| Fourier Kplus/Kminus | `tests/test_fourier.py` J1-J17 | Covered |
| A12 | no named test; probe only with expanded LC | Partial / not formalized |
| A21 | no named test; probe only with expanded LC | Partial / not formalized |
| A11 | no tests or implementation | Missing |
| A22 | no tests or implementation | Missing |
| Modulo compact relations | no tests or implementation | Missing |

Validation commands:

```text
pytest
py_compile
ruff check .
git diff --check
```

Final validation result after creating this report:

```text
pytest: 154 passed
py_compile: passed
ruff check .: All checks passed!
git diff --check: passed
git status --short: ?? docs/PHASE_K_AUDIT.md
```

## 16. Aij Gap Analysis

| Block | Representable today | Applicable today | Formal kernel | Explicit kernel | Correct relation |
|---|---|---|---|---|---|
| `A11` | Partial | No | No | No | No |
| `A12` | Partial | Yes | Yes | Partial | No |
| `A21` | Partial | Yes | Yes | Partial | No |
| `A22` | Partial | No | No | No | No |
| `A22^(1)` | Yes | Yes | Yes | Partial | No |

Explanations:

- `A11`: one could invent ad hoc atoms, but `Pplus/Pminus` and exact `1/2` are
  missing. No action exists.
- `A12`: the normalized WH model is representable only as an expanded
  `LinearCombination`, not as `(V-G)W`. It applies formally. Explicit kernels
  work with local rules.
- `A21`: same as `A12`, with the important exterior minus encoded by
  coefficients `(-1, +1)`.
- `A22`: same diagonal blockers as `A11`.
- `A22^(1)`: the Schur correction kernel/action exists as `C22`, but relation
  typing and explicit Fourier integration into formal `C22` are not built in.

## 17. Recommended K1/K2/K3

### K1: exact off-diagonal normalized WH blocks for m=2

Scope:

- Add exact coefficient support (`Fraction` or `sympy.Rational`) to the AST.
- Add explicit relation wrappers:
  `ExactBlock`, `WienerHopfModel`, `ModCompactRelation`.
- Add m=2 normalized off-diagonal constructors:
  `A12_wh_model`, `A21_wh_model`.
- Keep them expanded internally to avoid adding general distributivity.
- Preserve signs:
  `A12_WH = +Vtilde_alpha1*Wplus_12 - G1*Wplus_12`;
  `A21_WH = -Vtilde_alpha2*Wminus_21 + G2*Wminus_21`.

Do not include:

- diagonal blocks;
- arbitrary m;
- general Schur recursion.

Tests:

- exact `1/2` accepted and preserved;
- `Fraction(1,2)` or `Rational(1,2)` round-trips;
- `A12_wh_model` coefficients `(1, -1)`;
- `A21_wh_model` coefficients `(-1, 1)`;
- relation object says model is modulo compact, not exact equality.

### K2: formal-to-explicit kernel connection

Scope:

- Add a small kernel-rule API for replacing `Wplus_12`/`Wminus_21` formal
  kernels with Fourier-localized kernels.
- Keep formal and explicit modes both available.
- Do not infer replacements from visual SymPy expression order.

Tests:

- formal A12/A21 still use `Lplus_12`, `Lminus_21`;
- explicit A12/A21 use `chi_infinity` and Fourier denominators;
- ordered `AppliedTerm` provenance remains intact;
- no duplicate factor `2` or `-2`.

### K3: diagonal blocks

Scope:

- Add atoms for identity and `S_Rplus`.
- Add exact projections `Pplus = 1/2(I + S_Rplus)` and
  `Pminus = 1/2(I - S_Rplus)`.
- Add enough actions to apply the diagonal blocks formally.
- Keep any analytic realization of `S_Rplus` tightly scoped.

Tests:

- exact `Pplus/Pminus` coefficients;
- `A11` and `A22` formal construction;
- identity action;
- singular-integral action if implemented;
- exact equality vs modulo compact relation remains explicit.

General indexing should wait until after K1/K2 stabilize. A Strategy B style
constructor can expose `A(k,j)` metadata early while still dispatching only m=2
implemented cases.
