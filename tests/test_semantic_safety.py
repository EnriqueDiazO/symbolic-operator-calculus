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
    FormalIdentity,
    G1,
    G2,
    KernelAnnotatedExpression,
    KernelRepresentationRequiredError,
    KernelRepresentationStatus,
    ModCompactEquivalence,
    ModCompactRelation,
    R11,
    RegularizerOperator,
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
    build_first_schur_derivation_trace,
    c22_integrand,
    combined_kernel_c22,
    mvp_atomic_rules,
    require_exact_identity,
)


def test_semantic_relation_types_are_disjoint_and_immutable():
    exact = ExactIdentity("A", "B", evidence="symbolic derivation")
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
        ExactIdentity("A", "B", evidence="symbolic derivation")
    ).right == "B"


def test_approximate_equality_requires_tolerance_norm_and_residual():
    with pytest.raises(ValueError):
        ApproximateEquality("A", "B", 0, "norm", "residual")
    with pytest.raises(ValueError):
        ApproximateEquality("A", "B", 1, None, "residual")
    with pytest.raises(ValueError):
        ApproximateEquality("A", "B", 1, "norm", None)


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

    assert relation.certification_status is CertificationStatus.CERTIFIED
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
