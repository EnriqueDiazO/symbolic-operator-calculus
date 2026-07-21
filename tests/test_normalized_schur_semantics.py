import pytest

from symbolic_operator_calculus import (
    G1,
    G2,
    R1,
    T1_minus,
    T2_minus,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    Z1_inverse,
    Z2_inverse,
    AnalyticProofObligation,
    DerivationRelationKind,
    FactorClassification,
    FactorStatus,
    OperatorClass,
    SemanticDerivationStep,
    build_normalized_first_schur_pivot_derivation,
    normalized_first_schur_factor_classifications,
    normalized_first_schur_proof_obligations,
)


EXPECTED_FACTORS = (
    Z2_inverse,
    Vtilde_alpha2,
    G2,
    Wminus_21,
    T1_minus,
    R1,
    Z1_inverse,
    Vtilde_alpha1,
    G1,
    Wplus_12,
    T2_minus,
)

EXPECTED_CLASSES = (
    OperatorClass.MULTIPLICATION_OPERATOR,
    OperatorClass.MULTIPLICATION_DILATION,
    OperatorClass.MULTIPLICATION_OPERATOR,
    OperatorClass.LOCALIZED_WIENER_HOPF_OPERATOR,
    OperatorClass.AUXILIARY_SHIFT_OPERATOR,
    OperatorClass.MELLIN_PDO_REGULARIZER,
    OperatorClass.MULTIPLICATION_OPERATOR,
    OperatorClass.MULTIPLICATION_DILATION,
    OperatorClass.MULTIPLICATION_OPERATOR,
    OperatorClass.LOCALIZED_WIENER_HOPF_OPERATOR,
    OperatorClass.AUXILIARY_SHIFT_OPERATOR,
)

EXPECTED_STATUSES = (
    FactorStatus.EXACT,
    FactorStatus.EXACT,
    FactorStatus.EXACT,
    FactorStatus.MOD_COMPACT,
    FactorStatus.INVERTIBLE,
    FactorStatus.REGULARIZER_MOD_COMPACT,
    FactorStatus.EXACT,
    FactorStatus.EXACT,
    FactorStatus.EXACT,
    FactorStatus.MOD_COMPACT,
    FactorStatus.INVERTIBLE,
)

EXPECTED_SOURCES = (
    "definición del transporte",
    "transporte del shift",
    "coeficiente branchwise",
    "localización cuspídea 2025",
    "resultado branchwise",
    "KKL 2014, Thm. 5.8",
    "normalización diagonal",
    "transporte del shift",
    "coeficiente branchwise",
    "localización cuspídea 2025",
    "resultado branchwise",
)


def test_factor_classification_table_has_exact_requested_order_and_values():
    rows = normalized_first_schur_factor_classifications()

    assert tuple(row.factor for row in rows) == EXPECTED_FACTORS
    assert tuple(row.operator_class for row in rows) == EXPECTED_CLASSES
    assert tuple(row.status for row in rows) == EXPECTED_STATUSES
    assert tuple(row.mathematical_source for row in rows) == EXPECTED_SOURCES


def test_factor_metadata_is_declarative_and_does_not_certify_membership():
    row = next(
        item
        for item in normalized_first_schur_factor_classifications()
        if item.factor is R1
    )

    assert isinstance(row, FactorClassification)
    assert row.operator_class is OperatorClass.MELLIN_PDO_REGULARIZER
    assert row.status is FactorStatus.REGULARIZER_MOD_COMPACT
    assert not hasattr(row, "certify")
    assert not hasattr(row, "prove")
    assert R1.is_formal_regularizer


def test_trace_distinguishes_every_required_relation_kind():
    trace = build_normalized_first_schur_pivot_derivation()
    kinds = tuple(step.relation_kind for step in trace.semantic_steps)

    assert set(kinds) == set(DerivationRelationKind)
    assert kinds[:9] == (
        DerivationRelationKind.EXACT_EQUALITY,
        DerivationRelationKind.EXACT_EQUALITY,
        DerivationRelationKind.EXACT_EQUALITY,
        DerivationRelationKind.FORMAL_SUBSTITUTION,
        DerivationRelationKind.EXACT_EQUALITY,
        DerivationRelationKind.EXACT_EQUALITY,
        DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
        DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
        DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
    )
    assert all(
        kind is DerivationRelationKind.ANALYTIC_PROOF_OBLIGATION
        for kind in kinds[9:]
    )


def test_trace_step_payloads_reference_the_actual_relation_objects():
    trace = build_normalized_first_schur_pivot_derivation()
    by_key = {step.key: step for step in trace.semantic_steps}

    assert by_key["reduced_definition"].payload is trace.reduced_definition
    assert by_key["inverse_substitution"].payload is trace.inverse_substitution
    assert by_key["exact_result"].payload is trace.exact_result
    assert (
        by_key["normalized_mod_compact_result"].payload
        is trace.mod_compact_result
    )
    assert tuple(step.payload for step in trace.semantic_steps[9:]) == (
        trace.proof_obligations
    )


def test_all_five_analytic_proof_obligations_are_pending_and_explicit():
    obligations = normalized_first_schur_proof_obligations()

    assert tuple(item.key for item in obligations) == (
        "justify_wh_mod_compact",
        "classify_complete_correction",
        "prove_product_closure",
        "choose_analytic_framework",
        "apply_fredholm_criterion_later",
    )
    assert len(obligations) == 5
    assert all(isinstance(item, AnalyticProofObligation) for item in obligations)
    assert all(
        item.relation_kind
        is DerivationRelationKind.ANALYTIC_PROOF_OBLIGATION
        for item in obligations
    )
    assert "sin afirmar aún su conclusión" in obligations[-1].statement
    assert all(not hasattr(item, "discharged") for item in obligations)


def test_trace_reuses_deterministic_classifications_and_obligations():
    first = build_normalized_first_schur_pivot_derivation()
    second = build_normalized_first_schur_pivot_derivation()

    assert first.factor_classifications == second.factor_classifications
    assert first.proof_obligations == second.proof_obligations
    assert first.semantic_steps == second.semantic_steps


def test_semantic_metadata_validates_types_without_interpreting_evidence():
    with pytest.raises(TypeError):
        FactorClassification(
            R1,
            "MellinPDORegularizer",
            FactorStatus.REGULARIZER_MOD_COMPACT,
            "source",
        )
    with pytest.raises(TypeError):
        FactorClassification(
            R1,
            OperatorClass.MELLIN_PDO_REGULARIZER,
            "regularizer",
            "source",
        )
    with pytest.raises(ValueError):
        AnalyticProofObligation("empty", "   ")
    with pytest.raises(TypeError):
        SemanticDerivationStep("step", "EXACT_EQUALITY", object())
