from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest
import sympy as sp

from semantic_helpers import explicit_r11_kernel_representation
from symbolic_operator_calculus import (
    ApproximateEquality,
    CertificationStatus,
    ExactBlock,
    ExactIdentity,
    ExactIdentityRequiredError,
    ExactIdentityScope,
    FormalIdentity,
    FormalRegularizerAction,
    G1,
    G2,
    IntegrationDomain,
    KernelAnnotatedExpression,
    KernelRepresentation,
    KernelRepresentationRequiredError,
    KernelRepresentationStatus,
    LinearCombination,
    ModCompactEquivalence,
    ModCompactRelation,
    Product,
    R11,
    RegularizerOperator,
    Term,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    a11_formal_regularizer,
    a12_wh_model,
    a22_first_schur_correction,
    a22_first_schur_mod_compact_relation,
    a22_first_schur_reduction,
    apply_a22_first_schur_model_compact,
    apply_atom,
    apply_combined_kernel_c22,
    apply_linear_combination,
    apply_product,
    build_first_schur_derivation_trace,
    c22_integrand,
    combined_kernel_c22,
    extract_integral_kernel,
    mvp_atomic_rules,
    render_scalar_latex,
    require_exact_identity,
)


def kernel_representation(
    *,
    operator=None,
    status=KernelRepresentationStatus.FORMAL,
    hypotheses=("caller hypothesis",),
    evidence=None,
    name="K",
):
    output, input_ = sp.symbols(f"_{name}_output _{name}_input")
    return KernelRepresentation(
        operator=a11_formal_regularizer() if operator is None else operator,
        kernel_expression=sp.Function(name)(output, input_),
        output_variable=output,
        input_variable=input_,
        integration_domain=IntegrationDomain(0, sp.oo, "positive half-line"),
        semantic_status=status,
        hypotheses=hypotheses,
        evidence=evidence,
    )


def test_semantic_relation_types_are_disjoint_and_immutable():
    exact = ExactIdentity(
        "A",
        "B",
        evidence="symbolic derivation",
        scope=ExactIdentityScope.STRUCTURAL,
    )
    formal = FormalIdentity("A", "B", justification="formal substitution")
    modulo_compact = ModCompactEquivalence("A", "B")
    approximate = ApproximateEquality(
        "A",
        "B",
        tolerance=sp.Rational(1, 1000),
        norm="operator norm",
        residual="A-B",
    )

    assert len({type(exact), type(formal), type(modulo_compact), type(approximate)}) == 4
    assert exact != formal != modulo_compact != approximate
    assert not hasattr(formal, "exact")
    assert not hasattr(modulo_compact, "exact")
    assert not hasattr(approximate, "exact")
    with pytest.raises(FrozenInstanceError):
        exact.left = "changed"


def test_exact_consumer_rejects_modulo_compact_equivalence():
    relation = ModCompactEquivalence("A", "B", evidence="external lemma")

    with pytest.raises(ExactIdentityRequiredError):
        require_exact_identity(relation)
    assert require_exact_identity(
        ExactIdentity(
            "A",
            "B",
            evidence="symbolic derivation",
            scope=ExactIdentityScope.STRUCTURAL,
        )
    ).right == "B"


@pytest.mark.parametrize(
    "non_exact",
    (
        FormalIdentity("A", "B"),
        ModCompactEquivalence("A", "B"),
        ApproximateEquality("A", "B", 1, "real norm", "A-B"),
        True,
        {"left": "A", "right": "B", "evidence": "looks exact"},
    ),
)
def test_exact_consumer_rejects_every_non_exact_and_duck_typed_value(non_exact):
    with pytest.raises(ExactIdentityRequiredError):
        require_exact_identity(non_exact)


def test_exact_identity_requires_explicit_scope_and_retains_hypotheses():
    with pytest.raises(TypeError):
        ExactIdentity("A", "B", evidence="derivation", scope="structural")

    identity = ExactIdentity(
        "A",
        "B",
        evidence="derivation",
        scope=ExactIdentityScope.SCALAR_SYMBOLIC,
        hypotheses=("x is real",),
    )

    assert identity.scope is ExactIdentityScope.SCALAR_SYMBOLIC
    assert identity.hypotheses == ("x is real",)
    with pytest.raises(TypeError):
        ExactIdentity(
            "A",
            "B",
            evidence="derivation",
            scope=ExactIdentityScope.STRUCTURAL,
            hypotheses="one hypothesis, not characters",
        )


def test_approximate_equality_requires_tolerance_norm_and_residual():
    with pytest.raises(ValueError):
        ApproximateEquality("A", "B", 0, "norm", "residual")
    with pytest.raises(ValueError):
        ApproximateEquality("A", "B", 1, None, "residual")
    with pytest.raises(ValueError):
        ApproximateEquality("A", "B", 1, "norm", None)


@pytest.mark.parametrize("tolerance", (-1, sp.nan, sp.oo, -sp.oo, sp.I))
def test_approximate_equality_rejects_invalid_or_criterion_incompatible_tolerance(
    tolerance,
):
    with pytest.raises(ValueError):
        ApproximateEquality("A", "B", tolerance, "real-valued norm", "A-B")


def test_approximate_equality_rejects_boolean_tolerance():
    with pytest.raises(TypeError):
        ApproximateEquality("A", "B", True, "real-valued norm", "A-B")


def test_regularizer_without_kernel_representation_cannot_generate_integral():
    x, y = sp.symbols("x y")
    f = sp.Function("f")
    action = mvp_atomic_rules()[R11]

    assert isinstance(action.regularizer, RegularizerOperator)
    assert action.regularizer == a11_formal_regularizer()
    assert action.kernel_representation is None
    assert not hasattr(action, "kernel")
    with pytest.raises(KernelRepresentationRequiredError, match="regularizer is formal"):
        apply_atom(
            R11,
            f(x),
            x,
            mvp_atomic_rules(),
            integration_variable=y,
        )


def test_indirect_product_and_linear_combination_cannot_materialize_r11():
    x, y = sp.symbols("x y")
    f = sp.Function("f")
    product = Product((R11,))
    combination = LinearCombination((Term(1, product),))

    with pytest.raises(KernelRepresentationRequiredError):
        apply_product(
            product,
            f(x),
            x,
            mvp_atomic_rules(),
            integration_variable=y,
        )
    with pytest.raises(KernelRepresentationRequiredError):
        apply_linear_combination(
            combination,
            f(x),
            x,
            mvp_atomic_rules(),
            integration_variable=y,
        )


@pytest.mark.parametrize(
    "unsafe_builder",
    (
        lambda x, y, u, v, f: c22_integrand(x, u, v, y),
        lambda x, y, u, v, f: combined_kernel_c22(x, y),
        lambda x, y, u, v, f: apply_combined_kernel_c22(x, f, y),
        lambda x, y, u, v, f: apply_a22_first_schur_model_compact(f(x), x),
        lambda x, y, u, v, f: build_first_schur_derivation_trace(
            f(x),
            x,
            input_variable=y,
            outer_variable=u,
            middle_variable=v,
        ),
    ),
)
def test_every_public_schur_kernel_builder_requires_representation(unsafe_builder):
    x, y, u, v = sp.symbols("x y u v")

    with pytest.raises(KernelRepresentationRequiredError):
        unsafe_builder(x, y, u, v, sp.Function("f"))


def test_production_code_does_not_create_an_automatic_r11_function():
    source_root = Path(__file__).parents[1] / "src" / "symbolic_operator_calculus"
    source = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (source_root / "actions.py", source_root / "kernels.py")
    )

    assert 'sp.Function("R11")' not in source
    assert "R11_kernel" not in source


def test_explicit_kernel_representation_allows_annotated_integral():
    x, y = sp.symbols("x y")
    f = sp.Function("f")
    representation = explicit_r11_kernel_representation()
    rules = mvp_atomic_rules(regularizer_kernel=representation)

    result = apply_atom(R11, f(x), x, rules, integration_variable=y)

    assert isinstance(result, KernelAnnotatedExpression)
    assert isinstance(result.expression, sp.Integral)
    assert result.expression == sp.Integral(
        sp.Function("R11")(x, y) * f(y),
        (y, 0, sp.oo),
    )
    assert result.kernel_representations == (representation,)
    assert result.semantic_statuses == (KernelRepresentationStatus.FORMAL,)
    assert result.hypotheses == representation.hypotheses


@pytest.mark.parametrize(
    "unsafe_builder",
    (
        lambda wrong, x, y, u, v, f: c22_integrand(
            x, u, v, y, regularizer_kernel=wrong
        ),
        lambda wrong, x, y, u, v, f: combined_kernel_c22(
            x, y, regularizer_kernel=wrong
        ),
        lambda wrong, x, y, u, v, f: apply_combined_kernel_c22(
            x, f, y, regularizer_kernel=wrong
        ),
        lambda wrong, x, y, u, v, f: apply_a22_first_schur_model_compact(
            f(x), x, regularizer_kernel=wrong
        ),
        lambda wrong, x, y, u, v, f: build_first_schur_derivation_trace(
            f(x),
            x,
            input_variable=y,
            outer_variable=u,
            middle_variable=v,
            regularizer_kernel=wrong,
        ),
    ),
)
def test_representation_for_another_regularizer_is_rejected_everywhere(
    unsafe_builder,
):
    x, y, u, v = sp.symbols("x y u v")
    f = sp.Function("f")
    other_regularizer = RegularizerOperator("other target", "other regularizer")
    wrong = kernel_representation(operator=other_regularizer, name="K_other")

    with pytest.raises(ValueError, match="represent this regularizer"):
        FormalRegularizerAction(a11_formal_regularizer(), wrong)
    with pytest.raises(ValueError, match="A11 regularizer"):
        unsafe_builder(wrong, x, y, u, v, f)


def test_externally_certified_representation_requires_evidence():
    with pytest.raises(ValueError, match="does not verify"):
        kernel_representation(
            status=KernelRepresentationStatus.EXTERNALLY_CERTIFIED,
            evidence=None,
            name="K_external",
        )
    with pytest.raises(TypeError, match="hypothesis objects"):
        kernel_representation(
            hypotheses="one hypothesis, not characters",
            name="K_bad_hypotheses",
        )


def test_composition_retains_mixed_statuses_hypotheses_evidence_and_operator():
    x, y, u = sp.symbols("x y u")
    f = sp.Function("f")
    regularizer = a11_formal_regularizer()
    formal = kernel_representation(
        status=KernelRepresentationStatus.FORMAL,
        hypotheses=("formal hypothesis",),
        name="K_formal",
    )
    assumed = kernel_representation(
        status=KernelRepresentationStatus.ASSUMED,
        hypotheses=("assumed hypothesis",),
        evidence="external assumption record",
        name="K_assumed",
    )

    first = FormalRegularizerAction(regularizer, formal).apply(
        f(x), x, integration_variable=y
    )
    composed = FormalRegularizerAction(regularizer, assumed).apply(
        first, x, integration_variable=u
    )

    assert composed.kernel_representations == (formal, assumed)
    assert composed.semantic_statuses == (
        KernelRepresentationStatus.FORMAL,
        KernelRepresentationStatus.ASSUMED,
    )
    assert composed.hypotheses == ("formal hypothesis", "assumed hypothesis")
    assert composed.evidences == (None, "external assumption record")
    assert composed.represented_operators == (regularizer, regularizer)
    assert "FORMAL" in repr(composed) and "external assumption record" in repr(
        composed
    )
    assert isinstance(hash(composed), int)
    assert not hasattr(composed, "doit")
    assert not hasattr(composed, "simplify")


def test_kernel_extraction_retains_annotations_instead_of_silently_unwrapping():
    x, y = sp.symbols("x y")
    f = sp.Function("f")
    representation = kernel_representation(name="K_extract")
    applied = FormalRegularizerAction(
        a11_formal_regularizer(), representation
    ).apply(f(x), x, integration_variable=y)

    extracted = extract_integral_kernel(applied, x, f, y)

    assert isinstance(extracted, KernelAnnotatedExpression)
    assert extracted.expression == sp.Function("K_extract")(x, y)
    assert extracted.kernel_representations == (representation,)
    assert extracted.hypotheses == representation.hypotheses


def test_explicit_sympy_projection_does_not_recover_semantic_status_automatically():
    x, y = sp.symbols("x y")
    f = sp.Function("f")
    representation = kernel_representation(name="K_projection")
    annotated = FormalRegularizerAction(
        a11_formal_regularizer(), representation
    ).apply(f(x), x, integration_variable=y)

    projected = annotated.as_expr()
    simplified = sp.simplify(projected)

    assert isinstance(projected, sp.Expr)
    assert not isinstance(projected, KernelAnnotatedExpression)
    assert not hasattr(projected, "semantic_statuses")
    assert not hasattr(simplified, "semantic_statuses")
    with pytest.raises(AttributeError):
        _ = annotated.function


def test_rendering_makes_caller_status_and_evidence_visible():
    x, y = sp.symbols("x y")
    f = sp.Function("f")
    formal = kernel_representation(name="K_render")
    external = kernel_representation(
        status=KernelRepresentationStatus.EXTERNALLY_CERTIFIED,
        evidence="external certificate",
        name="K_render",
    )
    formal_expression = FormalRegularizerAction(
        a11_formal_regularizer(), formal
    ).apply(f(x), x, integration_variable=y)
    external_expression = FormalRegularizerAction(
        a11_formal_regularizer(), external
    ).apply(f(x), x, integration_variable=y)

    formal_latex = render_scalar_latex(formal_expression)
    external_latex = render_scalar_latex(external_expression)

    assert formal_latex != external_latex
    assert formal_expression != external_expression
    assert "kernel status: formal" in formal_latex
    assert "evidence objects: 0" in formal_latex
    assert "kernel status: externally certified" in external_latex
    assert "evidence objects: 1" in external_latex


def test_explicit_schur_kernel_preserves_status_and_hypotheses():
    x, y, u, v = sp.symbols("x y u v")
    representation = explicit_r11_kernel_representation()

    result = combined_kernel_c22(
        x,
        y,
        outer_variable=u,
        middle_variable=v,
        regularizer_kernel=representation,
    )

    assert isinstance(result, KernelAnnotatedExpression)
    assert result.kernel_representations == (representation,)
    assert result.semantic_statuses == (KernelRepresentationStatus.FORMAL,)
    assert result.hypotheses == representation.hypotheses
    assert result.expression.has(sp.Function("R11")(u, v))


def test_first_schur_reduction_and_correction_keep_order_and_signs():
    reduction = a22_first_schur_reduction()
    correction = a22_first_schur_correction()

    assert reduction.diagonal == ExactBlock("A", 2, 2)
    assert reduction.left == ExactBlock("A", 2, 1)
    assert reduction.regularizer.target == ExactBlock("A", 1, 1)
    assert reduction.right == ExactBlock("A", 1, 2)
    assert reduction.offdiagonal_product_coefficient == -1
    assert tuple(term.coefficient for term in correction.terms) == (1, -1, -1, 1)
    assert all(term.product.factors[2] is R11 for term in correction.terms)
    assert tuple(term.product.factors for term in correction.terms) == (
        (Vtilde_alpha2, Wminus_21, R11, Vtilde_alpha1, Wplus_12),
        (Vtilde_alpha2, Wminus_21, R11, G1, Wplus_12),
        (G2, Wminus_21, R11, Vtilde_alpha1, Wplus_12),
        (G2, Wminus_21, R11, G1, Wplus_12),
    )


def test_modulo_compact_relation_without_evidence_is_uncertified():
    relation = ModCompactRelation(
        ExactBlock("A", 1, 2),
        a12_wh_model(),
        space="X",
        compact_ideal="K(X)",
        residual="A12 - model",
    )

    assert relation.certification_status is CertificationStatus.UNCERTIFIED
    assert relation.semantic_relation.space == "X"
    assert relation.semantic_relation.compact_ideal == "K(X)"
    assert relation.semantic_relation.residual == "A12 - model"
    assert relation.semantic_relation.evidence is None


def test_modulo_compact_evidence_is_stored_but_not_invented():
    evidence = {"source": "externally supplied theorem"}
    relation = ModCompactEquivalence("A", "B", evidence=evidence)

    assert relation.certification_status is CertificationStatus.EVIDENCE_SUPPLIED
    assert relation.evidence is evidence


def test_modulo_compact_relation_implies_neither_structure_nor_exact_identity():
    relation = a22_first_schur_mod_compact_relation()

    assert relation.exact != relation.model
    assert relation.semantic_relation.left != relation.semantic_relation.right
    assert not isinstance(relation.semantic_relation, ExactIdentity)
    with pytest.raises(ExactIdentityRequiredError):
        require_exact_identity(relation)
    with pytest.raises(ExactIdentityRequiredError):
        require_exact_identity(relation.semantic_relation)
