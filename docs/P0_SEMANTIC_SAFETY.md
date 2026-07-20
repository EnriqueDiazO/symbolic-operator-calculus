# P0 semantic safety

This phase separates algebraic structure from analytic interpretation.  It
does not add Mellin calculus, compactness proofs, Fredholm theory, boundedness
results, or function-class membership tests.

## Relation types

`ExactIdentity`, `FormalIdentity`, `ModCompactEquivalence`, and
`ApproximateEquality` are distinct frozen dataclasses.  There is no implicit
conversion between them.  An API that requires exactness uses
`require_exact_identity`; a modulo-compact relation is rejected.

`ModCompactEquivalence` stores optional `space`, `compact_ideal`, `residual`,
and `evidence`.  With no evidence its `certification_status` is
`UNCERTIFIED`.  Supplying evidence records an external justification; it does
not make the software a compactness prover.

The compatibility types `ModCompactRelation` and `ModCompactSchurRelation`
retain their `exact` and `model` endpoint names and expose the new semantic
object as `semantic_relation`.

## Regularizers and kernel representations

`FormalRegularizer` is a `RegularizerOperator`; neither class has an implicit
kernel.  The default atomic rules therefore contain a
`FormalRegularizerAction` whose `kernel_representation` is `None`.

Calling an integral action for that regularizer raises
`KernelRepresentationRequiredError`.  Production code never creates
`sp.Function("R11")` automatically.

A caller that has an additional formal assumption or external result may
provide it explicitly:

```python
import sympy as sp
import symbolic_operator_calculus as soc

u0, v0 = sp.symbols("u0 v0")
representation = soc.KernelRepresentation(
    operator=soc.a11_formal_regularizer(),
    kernel_expression=sp.Function("R11")(u0, v0),
    output_variable=u0,
    input_variable=v0,
    integration_domain=soc.IntegrationDomain(0, sp.oo, "positive half-line"),
    semantic_status=soc.KernelRepresentationStatus.FORMAL,
    hypotheses=("R11 admits this formal kernel representation",),
)

rules = soc.mvp_atomic_rules(regularizer_kernel=representation)
```

Results involving the representation are returned as
`KernelAnnotatedExpression`.  Its `expression` field is the SymPy expression;
`kernel_representations`, `hypotheses`, and `semantic_statuses` retain the
semantic data.  Calling `as_expr()` is an explicit projection that discards
the wrapper, not a proof of convergence or operator existence.

The following APIs now require the explicit `regularizer_kernel` argument when
they reach the Schur kernel layer:

- `c22_integrand`;
- `combined_kernel_c22`;
- `apply_combined_kernel_c22`;
- `apply_a22_first_schur_model_compact`;
- `build_first_schur_derivation_trace`.

The algebraic AST expression

\[
A_{22}-A_{21}R_{11}A_{12}
\]

and its ordered four-term expansion remain available without a kernel.  Only
their automatic interpretation as ordinary integral operators is prohibited.

## Mathematical limits

A `KernelRepresentation`, even with externally supplied evidence, does not by
itself prove convergence, boundedness, compactness, Fredholmness, or the
existence of an ordinary kernel in a chosen function space.  A SymPy
`Integral` is a symbolic expression, not an analytic certificate.
