import sympy as sp

from symbolic_operator_calculus import (
    R1,
    U1,
    DerivationRelationKind,
    DilationFrequencyCase,
    ExactIdentity,
    ExactIdentityScope,
    FactorizedOscillatoryMellinSymbol,
    Product,
    RuleCertificationStatus,
    StandardMellinMembership,
    build_oscillatory_mellin_factor,
    build_right_dilation_composition_trace,
    render_right_dilation_composition_latex,
)


def test_exact_trace_keeps_r1_u1_order_and_atomic_regularizer():
    lam, gamma = sp.symbols("lambda gamma")
    trace = build_right_dilation_composition_trace(lam, gamma)

    assert trace.ast_product == Product((R1, U1))
    assert trace.ast_product.factors.count(R1) == 1
    assert R1.is_formal_regularizer
    assert isinstance(trace.factorized_symbol, FactorizedOscillatoryMellinSymbol)
    assert trace.factorized_symbol.ast_product is trace.ast_product
    assert trace.factorized_symbol.ordered_symbol_factors == (
        trace.factorized_symbol.standard_factor,
        trace.factorized_symbol.oscillatory_factor,
    )


def test_exact_identity_and_rule_have_operatorial_strength_and_definition_proof():
    lam, gamma = sp.symbols("lambda gamma")
    trace = build_right_dilation_composition_trace(lam, gamma)

    assert isinstance(trace.exact_identity, ExactIdentity)
    assert trace.exact_identity.scope is ExactIdentityScope.OPERATORIAL
    assert trace.exact_identity.left is trace.ast_product
    assert trace.exact_identity.right is trace.factorized_symbol
    assert trace.rule.name == "exact_right_dilation_mellin_composition"
    assert trace.rule.relation_kind is DerivationRelationKind.EXACT_EQUALITY
    assert trace.rule.certification_status is RuleCertificationStatus.CERTIFIED_EXACT
    assert "direct derivation" in trace.rule.source
    assert "integrations are not interchanged" in (
        trace.exact_identity.evidence.substitution
    )


def test_multiplier_sign_is_positive_and_weight_normalization_cancels():
    lam = sp.Symbol("lambda", real=True)
    gamma = sp.Symbol("gamma")
    positive_gamma = sp.Symbol("positive_gamma", positive=True)
    x = sp.Symbol("x", positive=True)
    kappa = sp.Symbol("kappa", real=True)
    trace = build_right_dilation_composition_trace(lam, gamma)

    assert trace.factorized_symbol.oscillatory_factor.expression == sp.Pow(
        gamma,
        sp.I * lam,
        evaluate=False,
    )
    normalization = (
        x**kappa
        * positive_gamma**kappa
        * (positive_gamma * x) ** (-kappa)
    )
    assert sp.powsimp(normalization, force=True) == 1
    assert "log-Gaussian" in trace.exact_identity.evidence.sign_check


def test_log_gaussian_check_selects_gamma_to_positive_i_lambda():
    lam = sp.Symbol("lambda", real=True)
    gamma = sp.Symbol("gamma", positive=True)
    u = sp.Symbol("u", real=True)
    log_gamma = sp.log(gamma)
    shifted = sp.exp(-((u + log_gamma) ** 2)) * sp.exp(-sp.I * lam * u)
    transformed = sp.integrate(shifted, (u, -sp.oo, sp.oo))
    expected = gamma ** (sp.I * lam) * sp.sqrt(sp.pi) * sp.exp(-(lam**2) / 4)

    assert sp.simplify(transformed - expected) == 0


def test_identity_and_nontrivial_frequency_cases_are_kept_disjoint():
    lam, gamma = sp.symbols("lambda gamma")
    identity = build_oscillatory_mellin_factor(
        lam,
        gamma,
        case=DilationFrequencyCase.IDENTITY,
    )
    nontrivial = build_oscillatory_mellin_factor(lam, gamma)

    assert identity.expression == 1
    assert identity.standard_v_membership is StandardMellinMembership.CERTIFIED
    assert nontrivial.expression == sp.Pow(gamma, sp.I * lam, evaluate=False)
    assert nontrivial.standard_v_membership is (
        StandardMellinMembership.EXCLUDED_NONTRIVIAL_DILATION
    )


def test_factorized_class_does_not_claim_standard_product_membership():
    lam, gamma = sp.symbols("lambda gamma")
    factorized = build_right_dilation_composition_trace(lam, gamma).factorized_symbol

    assert factorized.provisional_class == "tilde-E^(right-gamma)"
    assert factorized.standard_product_membership is (
        StandardMellinMembership.NOT_DEMONSTRATED
    )
    assert len(factorized.pending_properties) == 6
    assert all(
        item.key != "right_gamma_left_composition"
        for item in factorized.pending_properties
    )
    assert not hasattr(factorized, "as_expr")


def test_phase_q_latex_is_deterministic_and_has_no_false_membership_claim():
    lam, gamma = sp.symbols("lambda gamma")
    first = build_right_dilation_composition_trace(lam, gamma)
    second = build_right_dilation_composition_trace(lam, gamma)
    latex = render_right_dilation_composition_latex(first)

    assert latex == render_right_dilation_composition_latex(second)
    assert latex.startswith(r"R_{1}\,U_{1}\overset{\mathrm{exact}}{=}")
    assert r"\gamma_1^{i\lambda}" in latex
    assert r"\notin\widetilde" not in latex
    assert r"\text{not demonstrated}" in latex
