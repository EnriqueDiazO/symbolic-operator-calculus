# Phase K7 visual and semantic export audit

## 1. Git state

Audit date: 2026-07-10. Repository and working directory:

```text
/home/enriquedo/PersonalProjects/symbolic-operator-calculus
```

Branch and HEAD:

```text
main
79c14e2 feat: add structured trace for first Schur derivation
```

Initial status:

```text
 M notebooks/test.ipynb
 M src/symbolic_operator_calculus/__init__.py
 M src/symbolic_operator_calculus/derivations.py
 M src/symbolic_operator_calculus/kernels.py
?? src/symbolic_operator_calculus/exporting.py
?? src/symbolic_operator_calculus/rendering.py
?? tests/test_first_schur_latex_rendering.py
?? tests/test_first_schur_tex_export.py
?? tests/test_k7c_semantic_outputs.py
?? tests/test_semantic_latex_ordering.py
```

`git diff --cached --name-only` was empty. The audit did not alter any
pre-existing K7B/K7B.1/K7C file.

## 2. Baseline

The initial validation completed successfully:

```text
PYTHONPATH=src pytest -q                   314/314 passed
python3 -m compileall -q src tests         passed, no output
ruff check .                               All checks passed!
git diff --check                           passed, no output
```

## 3. Notebook execution

The repository notebook was not executed in place. A temporary kernelspec was
created under `/tmp/k7d-jupyter`, using the project's editable virtual
environment, and `jupyter nbconvert --execute` wrote the executed copy to
`/tmp/k7d-audit/executed.ipynb`. Execution finished with code 0. The original
notebook retained zero stored outputs.

The temporary output contains semantic `text/latex` payloads, rather than the
raw SymPy names that existed before K7C.

## 4. Four applied terms

The executed notebook displays exactly:

```text
Término 1; coeficiente = 1
Término 2; coeficiente = -1
Término 3; coeficiente = -1
Término 4; coeficiente = 1
```

The expressions use `L^-_{2,1}`, `L^+_{1,2}`, `R_{1,1}`, `\rho_1`,
`\rho_2`, `\gamma_1`, and `\gamma_2`. No internal scalar-function name is
present in the visible LaTeX payloads.

The expressions attached to terms 2 and 3 do not receive an extra leading
minus. Their `-1` belongs exclusively to `AppliedTerm.coefficient`.

## 5. Function versus multiplication operator

Operator-level output contains `G_{1}\,I` and `G_{2}\,I`, including the
off-diagonal models and correction expansion. Applied scalar expressions
contain `G_1(v)` and `G_2(x)` without `I`.

No visible occurrence of `G_1 I(v)` or `G_2 I(x)` was found. No operator model
subtracts a bare scalar `G_k` from `\widetilde V_{\alpha_k}`.

## 6. Plus/minus notation

The executed notebook and both exported documents use:

```latex
L^-_{2,1},\quad L^+_{1,2},\quad W^-_{2,1},\quad W^+_{1,2}.
```

Thus the sign is a superscript and the ordered index pair is a subscript.
None of the rejected variants (`Lminus_21`, `Lplus_12`, `L^{-}_{21}`, or
`L_{minus,21}`) appears in visible output.

## 7. M21

The actual `text/latex` output is ordered as:

```latex
M_{2,1}(x,u)
=\rho_2(x)L^-_{2,1}(\gamma_2x,u)
-G_2(x)L^-_{2,1}(x,u).
```

The shifted term is first, `\rho_2` precedes its kernel, and the `G_2` term is
second with a minus sign. It is not reordered through a SymPy `Add`.

## 8. M12

The actual `text/latex` output is ordered as:

```latex
M_{1,2}(v,y)
=\rho_1(v)L^+_{1,2}(\gamma_1v,y)
-G_1(v)L^+_{1,2}(v,y).
```

The shifted term and its scalar-factor order are preserved.

## 9. Off-diagonal models

The executed notebook shows `\simeq`, the correct signs, multiplication
operators, and factorization:

```latex
A_{2,1}\simeq
-\left(\widetilde V_{\alpha_2}-G_2 I\right)W^-_{2,1},
\qquad
A_{1,2}\simeq
\left(\widetilde V_{\alpha_1}-G_1 I\right)W^+_{1,2}.
```

## 10. Combined kernel C22

The compact structural formula is the first line of the displayed aligned
environment:

```latex
C_{2,2}^{(1)}(x,y)
=\int_0^\infty M_{2,1}(x,u)
  \int_0^\infty R_{1,1}(u,v)M_{1,2}(v,y)\,dv\,du.
```

The scalar expansion follows on the second line. The compact explanation
therefore precedes the large calculated expression.

## 11. Action on f

The first line displayed for the action is:

```latex
\left(C_{2,2}^{(1)}f\right)(x)
=\int_0^\infty C_{2,2}^{(1)}(x,y)f(y)\,dy.
```

The existing expanded action follows. The notebook does not present the large
integral without first identifying the acting operator.

## 12. Explicit Fourier realization

The notebook executed a second trace with
`explicit_wiener_hopf_rules(decay=d)`. In the realized lateral factors:

- `L^-_{2,1}` and `L^+_{1,2}` disappear;
- `\chi_\infty` appears at the localized variables;
- normalized denominators contain `2\pi d` and the corresponding
  `1\pm i(\cdot)/d` factor;
- `R_{1,1}` remains formal elsewhere in the trace;
- the shifted `\rho_k` term remains first and the `-G_k` term second.

No change to Fourier normalization was detected.

## 13. TeX export

Formal and explicit traces were rendered and exported to
`/tmp/k7d-audit/formal.tex` and `/tmp/k7d-audit/explicit.tex`. Each document:

- contains `\documentclass{article}`, `amsmath`, `amssymb`, and `geometry`;
- contains exactly nine sections in the renderer's order;
- contains no internal name from the audit's forbidden list;
- consumes `RenderedFirstSchurDerivation` rather than recalculating a trace.

No generated TeX or PDF file was written into the repository.

## 14. TeX compilation and visual inspection

`latexmk` 4.76, using pdfTeX, compiled both documents with exit code 0:

```text
formal.pdf     1 page, letter size
explicit.pdf   2 pages, letter size
```

There were no fatal TeX errors, but the logs and rasterized page inspection
revealed material horizontal overflow:

| Document | Line/area | Overfull amount |
|---|---|---:|
| formal | correction expansion | 70.72 pt |
| formal | combined kernel | 27.36 pt |
| formal | compact model action | 563.45 pt |
| explicit | correction expansion | 70.72 pt |
| explicit | combined kernel | 105.29 pt |
| explicit | compact model action | 641.38 pt |

The right-hand portion of the compact-model action is outside the printable
page in the formal PDF. This is not merely cosmetic: the compiled document is
mathematically correct but not fully readable on the page.

## 15. Notebook as consumer

Structural inspection found eight code cells and zero stored outputs. The
notebook:

- imports and calls public package APIs;
- does not construct `Term`, `Product`, or `Integral`;
- obtains the four terms from `a22_first_schur_correction()` and
  `apply_linear_combination_ordered()`;
- obtains C22 and its action from their public APIs;
- uses `Math` with the rendering layer for affected final output;
- contains no absolute machine path.

It is a consumer of the mathematical pipeline, not a second source of truth.

## 16. Traceability

| Visible output | Source |
|---|---|
| Four applied terms | `a22_first_schur_correction()` + `apply_linear_combination_ordered()` |
| Off-diagonal models | real mod-compact relations through the derivation trace |
| M21/M12 | ordered `KernelCombination` fields |
| C22 | `combined_kernel_c22()` and the corresponding rendered trace step |
| Action on `f` | `apply_combined_kernel_c22()` |
| TeX document | `RenderedFirstSchurDerivation` |

No independently maintained mathematical formula was found in notebook code
or in the exporter.

## 17. Classified findings

### ALTO — exported long formulas exceed the page

The correction expansion and, especially, the compact-model action are not
line-broken for document export. The maximum overfull box is 641.38 pt, and
raster inspection confirms clipped/off-page content. The exported source is
compilable, but its visual presentation is not checkpoint-quality.

### NO ES DEUDA — scalar expansion order inside C22

The compact M21/R11/M12 line preserves and explains the semantic structure.
The following fully expanded SymPy scalar expression may use commutative
ordering without changing the audited ordered M21 and M12 definitions. This is
not a second mathematical source.

### No critical or mathematical findings

No CRÍTICO issue, second source of truth, normalization change, or mathematical
error was detected.

## 18. Recommendation

**OPCIÓN B — Hace falta una corrección visual pequeña antes del checkpoint.**

Add deliberate line breaks or suitable multiline environments for the
correction expansion, combined-kernel expansion, and compact-model action,
then repeat the two-document compilation audit. Do not change the underlying
trace, scalar actions, signs, kernels, or Fourier normalization.
