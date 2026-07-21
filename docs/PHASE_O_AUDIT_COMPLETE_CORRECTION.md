# Phase O audit: complete first-pivot correction

Audit date: 2026-07-20. This audit was completed before changing production
code for Phase O.

Write access for this phase is restricted to
`/home/enriquedo/PersonalProjects/symbolic-operator-calculus`. The paper and
literature-index repositories were inspected only through read-only Git
commands.

## 1. Initial Git state

```text
pwd
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git rev-parse --show-toplevel
/home/enriquedo/PersonalProjects/symbolic-operator-calculus

git branch --show-current
p1a-mellin-symbol-infrastructure

git rev-parse HEAD
b7436bd2abd30abb1dd69c69a0db94887b6aedc3

git rev-parse --abbrev-ref --symbolic-full-name @{upstream}
origin/p1a-mellin-symbol-infrastructure

git rev-list --left-right --count @{upstream}...HEAD
0	3

git status --short
?? docs/FULL_CODE_AUDIT.md

git diff --cached --name-only
<empty>
```

The branch therefore began zero commits behind and three commits ahead of its
upstream. There were no unrelated tracked changes and nothing staged.

## 2. Integrity of the pre-existing untracked file

`docs/FULL_CODE_AUDIT.md` was present before Phase O and remains untracked. Its
initial fingerprint is:

```text
sha256  5d6797bf0ee751dd164f42a7af7438b5226c9e0c20ff79ef80ce78ea7d7eb4b1
size    57033 bytes
mtime   2026-07-16 19:32:24.724846953 -0600
mode    664
```

Git porcelain v2 reported `? docs/FULL_CODE_AUDIT.md`, and the staged-file list
was empty. This file must remain byte-identical, untracked, unstaged, and absent
from every Phase O commit.

## 3. Python environment and baseline

The shell itself had no activated virtual environment (`VIRTUAL_ENV=None`).
Commands were run with `$PWD/.venv/bin` prepended to `PATH`, resolving:

```text
Python 3.10.12
sys.prefix=/home/enriquedo/PersonalProjects/symbolic-operator-calculus/.venv
SymPy 1.14.0
pytest 9.1.1
Ruff 0.15.17
```

Effective baseline:

```text
$ pytest
........................................................................ [  8%]
........................................................................ [ 16%]
........................................................................ [ 24%]
........................................................................ [ 32%]
........................................................................ [ 40%]
........................................................................ [ 49%]
........................................................................ [ 57%]
........................................................................ [ 65%]
........................................................................ [ 73%]
........................................................................ [ 81%]
........................................................................ [ 90%]
........................................................................ [ 98%]
...............                                                          [100%]
879 passed in 45.25s

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

All four commands exited with status 0. Gate A is open.

## 4. Verified Phase N commits

HEAD contains the three expected commits in order:

```text
4f37d854411eb9192aed996c9c14e6f534f670b8
feat(schur): derive normalized first pivot formally

9d7db8609592065db9e65cefc57f6f10ffa1b991
feat(semantics): track relation kinds and proof obligations

b7436bd2abd30abb1dd69c69a0db94887b6aedc3
feat(schur): export normalized pivot derivation
```

Phase N created these files:

- `docs/PHASE_N_AUDIT_FIRST_SCHUR_PIVOT.md`;
- `docs/FIRST_SCHUR_PIVOT_CAPABILITY_REPORT.md`;
- `docs/FIRST_SCHUR_PIVOT_DERIVATION.tex`;
- `notebooks/first_schur_pivot_normalization.ipynb`;
- `src/symbolic_operator_calculus/normalized_schur.py`;
- `tests/test_normalized_first_schur_pivot.py`;
- `tests/test_normalized_schur_semantics.py`;
- `tests/test_first_schur_pivot_notebook.py`;
- `tests/test_normalized_first_schur_tex_export.py`.

It also extended `__init__.py`, `operators.py`, `substitution.py`,
`semantics.py`, `rendering.py`, and `exporting.py`.

## 5. External repositories at audit time

The paper repository existed at
`/home/enriquedo/PersonalProjects/Papers/DiazOcampoHasemanBVP2026`, with HEAD
`c1c0979c07e138d71783077359003c1af6e2a794`. Its `git status --short` output
was empty.

The literature index existed at
`/home/enriquedo/PersonalProjects/haseman-literature-index`, with HEAD
`0b4dbed213169f2e6ad9a0954ad764fe97a6295c`. Its `git status --short` output
was empty.

These HEAD values and clean statuses are the **initial** read-only integrity
checkpoints. During the later literature consultation, that external
repository appeared at HEAD
`0275d2e8332c000c5677062b3471a16105751fa7`, still with empty
`git status --short`. No Phase O command wrote to it; the independent HEAD
advance is recorded instead of treating the earlier value as a final
checkpoint.

## 6. Public Phase N API

### Normalized derivation

`NormalizedFirstSchurPivotDerivation` retains exact definitions, the formal
inverse substitution, diagonal recognition, the final exact formal result,
the modulo-compact result, factorization metadata, typed semantic steps,
factor classifications, and proof obligations.

`build_normalized_first_schur_pivot_derivation()` deterministically builds the
trace and exposes the expanded WH correction through
`expanded_mod_compact_result`.

### Ordered operator algebra and substitution

- `OperatorAtom`, `Product`, `Term`, and `LinearCombination` remain the only
  operator-expression AST.
- `compose_ordered()` explicitly distributes finite linear combinations while
  retaining tuple order. It never sorts, commutes, cancels, or combines terms.
- `substitute_operator_subproduct()` replaces contiguous subwords left to
  right and distributes the replacement. It performs no analytic inference.

### Relations and proof obligations

`DerivationRelationKind` distinguishes:

- `EXACT_EQUALITY`;
- `FORMAL_SUBSTITUTION`;
- `MOD_COMPACT_EQUIVALENCE`;
- `ANALYTIC_PROOF_OBLIGATION`.

`ExactIdentity`, `FormalIdentity`, and `ModCompactEquivalence` remain distinct
types. `AnalyticProofObligation` is a pending statement with no discharge or
certification method. `SemanticDerivationStep` binds a typed relation kind to a
payload without proving it.

### Factor classification

`FactorClassification`, `OperatorClass`, and `FactorStatus` store declarative
metadata. They do not infer membership or expose proof methods. Phase N already
classifies the factors of the unnormalized WH correction, including `R1` as a
declared `MellinPDORegularizer` with a modulo-compact regularizer status.

### Rendering and export

The public API provides deterministic exact and modulo-compact LaTeX renderers,
Markdown renderers for classifications and obligations, a rendered-trace
structure, and a reusable no-preamble TeX exporter. The specialized renderer
preserves factor order and visually separates `=`, formal substitution, and
`\simeq`.

## 7. Capability assessment for Phase O

The current public AST and ordered algebra are sufficient to build every
requested finite expansion without another AST:

- C0 can be represented by validated factorization metadata over one expanded
  `LinearCombination`;
- C1 follows from `compose_ordered()` applied to the two differences;
- C2 follows from the already public `t1_minus_expression()` and
  `t2_minus_expression()` while retaining the differences as factorization
  metadata;
- C3 follows from explicit distribution of all four two-term factors and must
  contain 16 ordered terms.

The exact normalizations

```text
Zk^(-1) (Vtilde_alphak - Gk) = Uk - Ghatk
```

can be recorded as supplied `ExactIdentity` objects and applied through
contiguous structural substitution. This is exact only for those supplied
normalization rules. It does not upgrade the global pivot relation, whose
Wiener--Hopf origin remains modulo compact operators.

## 8. Missing pieces and minimum justified extension

The API does not yet contain atoms for `U1`, `U2`, `Ghat1`, `Ghat2`, or the
complete-correction label. It also does not assign stable term IDs or retain
term-level signatures, motifs, open interfaces, and proof obligations.

The minimum extension is therefore:

1. register only the missing noncommutative atoms and their LaTeX names;
2. add a dedicated complete-correction trace that reuses `Product`,
   `LinearCombination`, `compose_ordered`, semantic relation types, and the
   Phase N trace;
3. add lightweight immutable term-analysis metadata whose authoritative
   expression is still the existing `Term`/`Product` AST;
4. add an exact contiguous motif recognizer and an interface classifier that
   never commute factors and never invent rules;
5. generate catalogs, matrices, notebook outputs, and TeX from those public
   structures.

No generic rewrite engine, operator-algebra membership prover, compactness
prover, Mellin symbol inversion, or parallel expression tree is justified.

## 9. Logical boundary before implementation

The following must remain separate throughout Phase O:

- the supplied exact normalization of the two `Zk`-prefixed factors;
- the modulo-compact Wiener--Hopf substitutions inherited from Phase N;
- subsequent finite formal distribution;
- declarative class labels;
- unproved analytic closure and membership questions.

In particular, removing `Z1_inverse` or `Z2_inverse` is permitted only through
the two named supplied exact identities. `R1` may not move, expand, or acquire
an inverse-symbol formula. No term may be reordered to manufacture a motif.

## 10. Audit conclusion

All Phase A checks passed. The expected Phase N history is present, the
baseline is green, the only dirty entry is the protected untracked audit file,
and both external repositories began clean. Phase O may proceed with the
narrow extension described above.
