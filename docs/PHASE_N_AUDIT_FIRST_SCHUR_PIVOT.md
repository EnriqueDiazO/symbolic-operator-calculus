# Phase N audit: normalized first Schur pivot

Audit date: 2026-07-20.

This audit was completed before changing the symbolic engine. It concerns only
`symbolic-operator-calculus`; no paper repository or paper LaTeX source was
inspected or modified.

## 1. Initial repository state

```text
pwd
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git rev-parse --show-toplevel
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git branch --show-current
p1a-mellin-symbol-infrastructure

git rev-parse HEAD
bc3ba4a671cb3c66a25665ccfd968e4b52ad0a72

git rev-parse --abbrev-ref --symbolic-full-name @{upstream}
origin/p1a-mellin-symbol-infrastructure

git status --short
?? docs/FULL_CODE_AUDIT.md
```

`docs/FULL_CODE_AUDIT.md` was present and untracked before this phase. It is
user-owned state and must not be edited, removed, staged, or included in a
phase commit.

No `AGENTS.md` file is present in the repository.

## 2. Python environment and tool versions

The shell initially had no active virtual environment:

```text
VIRTUAL_ENV=None
```

The unqualified `python` command initially failed because the active pyenv
selection did not provide that shim. The repository does contain `.venv`, so
the reproducible baseline was run with `$PWD/.venv/bin` prepended to `PATH`.
That makes the commands retain the requested spellings while selecting the
repository environment.

```text
Python 3.10.12
sys.prefix=/home/enriquedo/PersonalProjects/symbolic-operator-calculus/.venv
base_prefix=/usr
SymPy 1.14.0
pytest 9.1.1
Ruff 0.15.17
```

Ruff is supplied by the surrounding user environment; `.venv` has no private
Ruff executable. For reference, the initially resolved user-level `pytest`
outside `.venv` reported 9.0.3.

## 3. Reproducible baseline

With `PATH="$PWD/.venv/bin:$PATH"`:

```text
$ pytest
........................................................................ [  8%]
........................................................................ [ 16%]
........................................................................ [ 25%]
........................................................................ [ 33%]
........................................................................ [ 42%]
........................................................................ [ 50%]
........................................................................ [ 59%]
........................................................................ [ 67%]
........................................................................ [ 76%]
........................................................................ [ 84%]
........................................................................ [ 93%]
.........................................................                [100%]
849 passed in 43.75s

$ python -m compileall src tests
Listing 'src'...
Listing 'src/symbolic_operator_calculus'...
Listing 'src/symbolic_operator_calculus/mellin'...
Listing 'src/symbolic_operator_calculus.egg-info'...
Listing 'tests'...

$ ruff check .
All checks passed!

$ git diff --check
<no output>
```

All four effective baseline commands exited with status 0.

## 4. Relevant repository structure

```text
src/symbolic_operator_calculus/
  __init__.py                 public API
  operators.py                noncommutative AST
  blocks.py                   exact blocks and WH/Schur constructors
  relations.py                exact-block and modulo-compact metadata
  semantics.py                semantic relation types and certification status
  derivations.py              existing first-Schur kernel/action trace
  kernels.py                  correction factorization and kernels
  rendering.py                ordered LaTeX rendering
  exporting.py                complete-document TeX export
  substitution.py             capture-safe scalar SymPy substitution only
  schur.py                    compact-model action
  mellin/                     scalar Mellin domains, expressions, symbols

tests/
  test_first_schur_reduction.py
  test_first_schur_derivation_trace.py
  test_first_schur_latex_rendering.py
  test_first_schur_tex_export.py
  test_noncommutative_expansion.py
  test_relations.py
  test_semantic_safety.py
  test_substitution.py
  ...

notebooks/
  relative_wiener_hopf.ipynb
  test.ipynb
  thesis_symbolic_operator_calculus_guide.ipynb

docs/
  MATHEMATICAL_CONVENTIONS.md
  P0_SEMANTIC_SAFETY.md
  PHASE_K5_AUDIT.md
  PHASE_K6_END_TO_END_AUDIT.md
  PHASE_K7_VISUAL_EXPORT_AUDIT.md
  PHASE_L1_DILATION_WIENER_HOPF_SYMBOL.md
  PHASE_L2_RELATIVE_WH_OPERATOR_TRACE.md
  ...
```

## 5. What the public API already represents

- `OperatorAtom`, `Product`, `Term`, and `LinearCombination` form one small,
  immutable noncommutative AST. Products store factor order directly in tuples;
  operator algebra is not delegated to commutative SymPy objects.
- `Product` contains only `OperatorAtom` factors. `LinearCombination` contains
  ordered scalar-weighted products and preserves term order.
- `R11` is an atom with `kind="formal_regularizer"`; it has no implicit inverse,
  symbol, or kernel. An explicit caller-supplied kernel representation is
  required before kernel-level action.
- `ExactBlock`, `FirstSchurReduction`, and `FormalRegularizer` record the formal
  unnormalized reduction
  `A22 - A21 R11 A12` without treating exact blocks as AST atoms.
- `a12_mod_compact_relation()` and `a21_mod_compact_relation()` separately
  declare the two normalized Wiener--Hopf models. Their semantic relation is a
  `ModCompactEquivalence`, uncertified unless external evidence is supplied.
- `ExactIdentity`, `FormalIdentity`, `ModCompactEquivalence`, and
  `ApproximateEquality` are distinct types. `require_exact_identity` rejects a
  modulo-compact relation where exactness is required.
- `factor_first_schur_correction()` validates and reconstructs the existing
  four-term correction around the atomic regularizer.
- `render_product_latex` and `render_linear_combination_latex` preserve stored
  order. The specialized first-Schur renderer visibly distinguishes `=` and
  `\simeq`.
- `export_first_schur_derivation_tex()` writes an already-rendered derivation,
  but it emits a complete LaTeX document with a preamble rather than a reusable
  fragment.

## 6. What the public API already derives

The existing m=2 pipeline derives the unnormalized formal reduction and its
four-term Wiener--Hopf correction:

```text
A22^(1) = A22 - A21 R11 A12
A21 ~= -(Vtilde_alpha2 - G2) Wminus_21
A12 ~=  (Vtilde_alpha1 - G1) Wplus_12
```

It correctly derives the positive correction sign as the product of the outer
Schur minus sign and the leading minus sign of the A21 model. It can then apply
the expanded correction to operands and, only when a regularizer kernel is
explicitly supplied, build kernel/action expressions.

This is related to the requested calculation but is not the requested
normalization: the current correction uses `R11` directly and contains none of
`Z_2^{-1}`, `T_{1,-}`, `Z_1^{-1}`, or `T_{2,-}`.

## 7. What the package cannot certify

The current package does not certify any of the following:

- compactness of the residual in a supplied Wiener--Hopf equivalence;
- boundedness, algebra membership, closure under the displayed products, or
  admissibility as a Mellin PDO;
- Fredholmness or an index statement;
- that a caller-supplied regularizer is an exact inverse;
- the analytic validity of `A_{1,1}^{(-1)} = T_{1,-} R_1 Z_1^{-1}`;
- the exact normalized identity
  `Z_2^{-1} A_{2,2} T_{2,-} = N_2` unless it is supplied as an explicit formal
  rule/evidence-bearing identity;
- propagation of modulo-compact equivalence through arbitrary surrounding
  products as an analytic theorem.

Construction of a formal or semantic relation records its asserted scope. It
is not a proof of the asserted mathematical relation.

## 8. Blocking API gaps for the normalized pivot

The calculation cannot be completed solely with the current public API:

1. The required normalization and pivot atoms are not registered.
2. `substitution.py` substitutes free scalar SymPy variables inside integrals;
   it does not substitute operator atoms or ordered subwords in the operator
   AST.
3. Multiplication by a `LinearCombination` intentionally raises an error, so a
   grouped product such as `(Vtilde_alpha2 - G2) Wminus_21 ...` cannot itself be
   an AST `Product`. The package only stores its controlled four-term
   expansion and reconstructs factorization through validated metadata.
4. The existing `FirstSchurDerivationTrace` targets kernels/actions and does not
   retain the normalization, inverse substitution, `N_2` recognition, or a
   sequence of typed relation kinds.
5. There is no proof-obligation type and no factor-classification output.
6. The current exporter produces a standalone document, not the requested
   reusable no-preamble fragment.

Manually instantiating `Term` and `Product` could spell the final expanded
words, but that would not derive or validate the requested substitutions and
would not meet the traceability goal.

## 9. Minimum justified extension

The minimum extension should reuse the existing AST and add no parallel
expression tree:

- register only the missing `OperatorAtom` constants and their central LaTeX
  names;
- add dedicated normalized-first-pivot builders that construct each expression
  from `Product`, `Term`, and `LinearCombination` and validate substitutions by
  exact ordered subword replacement;
- add a small immutable trace whose steps carry an explicit enum with
  `EXACT_EQUALITY`, `FORMAL_SUBSTITUTION`,
  `MOD_COMPACT_EQUIVALENCE`, or `ANALYTIC_PROOF_OBLIGATION`;
- represent the factorized WH correction as validated metadata over its
  existing expanded `LinearCombination`, following the precedent of
  `FirstSchurCorrectionFactorization`, rather than allowing nested sums in
  `Product`;
- store factor classifications and proof obligations as presentation/semantic
  metadata without claiming that the software proves them;
- provide deterministic renderers and a reusable-fragment exporter for this
  trace.

No generic algebra-membership engine, compactness prover, Mellin inverse-symbol
calculus, or broad AST refactor is justified by this task.

## 10. Audit conclusion

The API already provides the correct safety foundation: immutable ordered
products, controlled linear combinations, atomic regularizers, and distinct
semantic relation types. It also contains the unnormalized first-Schur model
and the exact signs needed for the two off-diagonal substitutions.

It does **not** currently expose enough public structure to derive the
normalized pivot step by step. A narrow extension is required. That extension
must treat both normalized identities as supplied formal rules, keep the
Wiener--Hopf replacement modulo compacta, keep `R_1` atomic, and emit analytic
proof obligations rather than analytic conclusions.
