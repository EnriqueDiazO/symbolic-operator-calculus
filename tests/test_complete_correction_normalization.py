from symbolic_operator_calculus import (
    C2,
    C2_normalized,
    G1,
    G2,
    Ghat1,
    Ghat2,
    N2,
    N2_first,
    R1,
    T1_minus,
    T2_minus,
    U1,
    U2,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    Z1_inverse,
    Z2_inverse,
    DerivationRelationKind,
    ExactIdentity,
    ModCompactEquivalence,
    OperatorRule,
    RuleCertificationStatus,
    build_complete_correction_derivation,
    complete_correction_expression,
    factor_normalized_complete_correction,
    normalized_complete_correction_expression,
    render_complete_correction_derivation_latex,
    render_normalized_complete_correction_latex,
    render_original_complete_correction_latex,
)


ORIGINAL_FACTORS = (
    (
        Z2_inverse,
        Vtilde_alpha2,
        Wminus_21,
        T1_minus,
        R1,
        Z1_inverse,
        Vtilde_alpha1,
        Wplus_12,
        T2_minus,
    ),
    (
        Z2_inverse,
        Vtilde_alpha2,
        Wminus_21,
        T1_minus,
        R1,
        Z1_inverse,
        G1,
        Wplus_12,
        T2_minus,
    ),
    (
        Z2_inverse,
        G2,
        Wminus_21,
        T1_minus,
        R1,
        Z1_inverse,
        Vtilde_alpha1,
        Wplus_12,
        T2_minus,
    ),
    (
        Z2_inverse,
        G2,
        Wminus_21,
        T1_minus,
        R1,
        Z1_inverse,
        G1,
        Wplus_12,
        T2_minus,
    ),
)

NORMALIZED_FACTORS = (
    (U2, Wminus_21, T1_minus, R1, U1, Wplus_12, T2_minus),
    (U2, Wminus_21, T1_minus, R1, Ghat1, Wplus_12, T2_minus),
    (Ghat2, Wminus_21, T1_minus, R1, U1, Wplus_12, T2_minus),
    (Ghat2, Wminus_21, T1_minus, R1, Ghat1, Wplus_12, T2_minus),
)


def coefficients(expression):
    return tuple(term.coefficient for term in expression.terms)


def factors(expression):
    return tuple(term.product.factors for term in expression.terms)


def test_original_complete_correction_reuses_phase_n_order_and_signs():
    expression = complete_correction_expression()

    assert coefficients(expression) == (1, -1, -1, 1)
    assert factors(expression) == ORIGINAL_FACTORS


def test_exact_normalization_removes_only_both_z_factors():
    expression = normalized_complete_correction_expression()

    assert coefficients(expression) == (1, -1, -1, 1)
    assert factors(expression) == NORMALIZED_FACTORS
    assert all(
        Z1_inverse not in term.product.factors
        and Z2_inverse not in term.product.factors
        for term in expression.terms
    )
    assert all(len(term.product.factors) == 7 for term in expression.terms)
    assert all(term.product.factors.count(R1) == 1 for term in expression.terms)


def test_normalized_factorization_reconstructs_without_commutation():
    factorization = factor_normalized_complete_correction()

    assert factorization.expanded == normalized_complete_correction_expression()
    assert coefficients(factorization.left_difference) == (1, -1)
    assert factors(factorization.left_difference) == ((U2,), (Ghat2,))
    assert factorization.bridge.factors == (T1_minus, R1)
    assert coefficients(factorization.right_difference) == (1, -1)
    assert factors(factorization.right_difference) == ((U1,), (Ghat1,))
    assert factorization.suffix.factors == (T2_minus,)


def test_global_pivot_relation_stays_mod_compact_while_z_rules_are_exact():
    trace = build_complete_correction_derivation()

    assert isinstance(trace.pivot_relation, ModCompactEquivalence)
    assert trace.pivot_relation.left is N2_first
    assert factors(trace.pivot_relation.right) == ((N2,), (C2,))
    assert isinstance(trace.left_normalization_rule.relation, ExactIdentity)
    assert isinstance(trace.right_normalization_rule.relation, ExactIdentity)
    assert isinstance(trace.normalized_identity, ExactIdentity)
    assert trace.normalized_definition.left is C2_normalized
    assert trace.normalized_definition.right == (
        normalized_complete_correction_expression()
    )
    assert tuple(step.relation_kind for step in trace.semantic_steps) == (
        DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
        DerivationRelationKind.EXACT_EQUALITY,
        DerivationRelationKind.EXACT_EQUALITY,
        DerivationRelationKind.EXACT_EQUALITY,
    )


def test_every_new_rule_has_name_logic_preconditions_source_and_status():
    trace = build_complete_correction_derivation()
    pivot_rule, left_rule, right_rule = trace.rules

    assert all(isinstance(rule, OperatorRule) for rule in trace.rules)
    assert tuple(rule.name for rule in trace.rules) == (
        "complete_correction_from_first_pivot",
        "normalize_left_branch_factor",
        "normalize_right_branch_factor",
    )
    assert pivot_rule.relation_kind is (
        DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
    )
    assert pivot_rule.certification_status is RuleCertificationStatus.UNCERTIFIED
    assert left_rule.certification_status is (
        RuleCertificationStatus.CERTIFIED_EXACT
    )
    assert right_rule.certification_status is (
        RuleCertificationStatus.CERTIFIED_EXACT
    )
    assert all(rule.preconditions for rule in trace.rules)
    assert all(rule.source for rule in trace.rules)


def test_complete_correction_latex_is_exact_and_deterministic():
    first = build_complete_correction_derivation()
    second = build_complete_correction_derivation()

    assert render_original_complete_correction_latex(first) == (
        r"\mathcal C_{2} = Z_{2}^{-1}\,"
        r"\left(\widetilde V_{\alpha_2} - G_{2}\right)\,W^-_{2,1}\,"
        r"T_{1,-}\,R_{1}\,Z_{1}^{-1}\,"
        r"\left(\widetilde V_{\alpha_1} - G_{1}\right)\,W^+_{1,2}\,T_{2,-}"
    )
    assert render_normalized_complete_correction_latex(first) == (
        r"\mathcal C_{2}^{\mathrm{norm}} = "
        r"\left(U_{2} - \widehat G_{2}\right)\,W^-_{2,1}\,T_{1,-}\,R_{1}\,"
        r"\left(U_{1} - \widehat G_{1}\right)\,W^+_{1,2}\,T_{2,-}"
    )
    rendered = render_complete_correction_derivation_latex(first)
    assert rendered == render_complete_correction_derivation_latex(second)
    assert rendered.steps[0].latex == (
        r"N_{2}^{(1)} \simeq N_{2} + \mathcal C_{2}"
    )
    assert " = " in rendered.steps[2].latex
    assert r"\simeq" not in rendered.steps[2].latex
    assert " = " in rendered.steps[3].latex
    assert r"\simeq" not in rendered.steps[3].latex


def test_r1_remains_atomic_and_no_inverse_symbol_is_computed():
    assert R1.is_formal_regularizer
    assert not hasattr(R1, "symbol")
    assert not hasattr(R1, "inverse")
    assert U1 * Wplus_12 != Wplus_12 * U1
    assert U2 * Wminus_21 != Wminus_21 * U2
