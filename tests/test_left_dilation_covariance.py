import sympy as sp

from symbolic_operator_calculus import (
    R1,
    U1,
    U1_inverse,
    AnalyticLemmaStatus,
    DerivationRelationKind,
    DilationFrequencyCase,
    ExactIdentityScope,
    Product,
    RuleCertificationStatus,
    StandardMellinMembership,
    build_left_dilation_covariance_trace,
    render_left_dilation_covariance_latex,
)


def test_left_covariance_keeps_order_and_has_two_exact_operator_identities():
    lam, gamma, x = sp.symbols("lambda gamma x")
    trace = build_left_dilation_covariance_trace(lam, gamma, x)

    assert trace.source_product == Product((U1_inverse, R1))
    assert trace.conjugated_product == Product((U1_inverse, R1, U1))
    assert trace.left_identity.scope is ExactIdentityScope.OPERATORIAL
    assert trace.conjugated_identity.scope is ExactIdentityScope.OPERATORIAL
    assert all(
        rule.relation_kind is DerivationRelationKind.EXACT_EQUALITY
        for rule in trace.rules
    )
    assert all(
        rule.certification_status is RuleCertificationStatus.CERTIFIED_EXACT
        for rule in trace.rules
    )


def test_inverse_frequency_radial_argument_and_exact_cancellation_are_correct():
    lam, gamma, x = sp.symbols("lambda gamma x")
    trace = build_left_dilation_covariance_trace(lam, gamma, x)

    assert trace.inverse_factor.expression == sp.Pow(
        gamma, -sp.I * lam, evaluate=False
    )
    assert trace.scaling.radial_argument == x / gamma
    assert trace.frequency_cancellation == 1
    assert trace.factorized_symbol.standard_product_membership is (
        StandardMellinMembership.NOT_DEMONSTRATED
    )


def test_log_gaussian_check_selects_gamma_to_negative_i_lambda():
    lam = sp.Symbol("lambda", real=True)
    gamma = sp.Symbol("gamma", positive=True)
    u = sp.Symbol("u", real=True)
    shifted = sp.exp(-((u - sp.log(gamma)) ** 2)) * sp.exp(-sp.I * lam * u)
    transformed = sp.integrate(shifted, (u, -sp.oo, sp.oo))
    expected = gamma ** (-sp.I * lam) * sp.sqrt(sp.pi) * sp.exp(-(lam**2) / 4)

    assert sp.simplify(transformed - expected) == 0


def test_radial_scaling_checks_every_project_contract():
    lam, gamma, x = sp.symbols("lambda gamma x")
    scaling = build_left_dilation_covariance_trace(lam, gamma, x).scaling

    assert scaling.status is AnalyticLemmaStatus.ANALYTICALLY_PROVED
    assert scaling.scaled_symbol.membership is StandardMellinMembership.CERTIFIED
    assert len(scaling.verified_contracts) == 6
    text = " ".join(scaling.verified_contracts)
    for fragment in ("C_b(V)", "continuity", "cm_t", "translation", "tail", "fibers"):
        assert fragment in text


def test_gamma_one_reduces_inverse_factor_and_product_to_standard_membership():
    lam, gamma, x = sp.symbols("lambda gamma x")
    trace = build_left_dilation_covariance_trace(
        lam,
        gamma,
        x,
        case=DilationFrequencyCase.IDENTITY,
    )

    assert trace.inverse_factor.expression == 1
    assert trace.factorized_symbol.standard_product_membership is (
        StandardMellinMembership.CERTIFIED
    )
    assert trace.frequency_cancellation == 1


def test_covariance_latex_is_deterministic_and_records_all_signs():
    lam, gamma, x = sp.symbols("lambda gamma x")
    first = build_left_dilation_covariance_trace(lam, gamma, x)
    second = build_left_dilation_covariance_trace(lam, gamma, x)
    latex = render_left_dilation_covariance_latex(first)

    assert latex == render_left_dilation_covariance_latex(second)
    assert r"a_{\gamma}(x,\lambda)=a(x/\gamma,\lambda)" in latex
    assert r"d_{\gamma^{-1}}(\lambda)=\gamma^{-i\lambda}" in latex
    assert r"D_{\gamma}^{-1}\operatorname{Op}(a)D_{\gamma}" in latex
