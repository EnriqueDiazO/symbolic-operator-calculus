from symbolic_operator_calculus import (
    Ghat1,
    R1,
    CertificationStatus,
    DerivationRelationKind,
    Product,
    RuleCertificationStatus,
    StandardMellinMembership,
    build_coefficient_semiproduct_trace,
    render_coefficient_semiproduct_latex,
)


def test_coefficient_case_preserves_order_and_is_not_an_exact_identity():
    trace = build_coefficient_semiproduct_trace()

    assert trace.ast_product == Product((R1, Ghat1))
    assert trace.relation.left is trace.ast_product
    assert trace.relation.right is trace.product_symbol
    assert trace.relation.certification_status is CertificationStatus.EVIDENCE_SUPPLIED
    assert trace.rule.relation_kind is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
    assert trace.rule.certification_status is (
        RuleCertificationStatus.CERTIFIED_MOD_COMPACT
    )


def test_all_three_coefficient_symbols_have_separately_recorded_membership():
    trace = build_coefficient_semiproduct_trace()

    for symbol in (trace.left_symbol, trace.right_symbol, trace.product_symbol):
        assert symbol.membership is StandardMellinMembership.CERTIFIED
        assert symbol.symbol_class == "tilde-E(R_+,V(R))"
    assert trace.product_symbol.name == "r_1 Ghat_1"


def test_semiproduct_evidence_has_exact_source_pages_checksum_and_hypotheses():
    evidence = build_coefficient_semiproduct_trace().relation.evidence

    assert evidence.theorem == "KKL2014TwoShifts Theorem 3.3"
    assert evidence.printed_page == "942"
    assert evidence.pdf_page == 8
    assert evidence.checksum == (
        "3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec"
    )
    assert len(evidence.paper_hypotheses) == 3


def test_coefficient_latex_is_deterministic_and_explicitly_mod_compact():
    first = build_coefficient_semiproduct_trace()
    second = build_coefficient_semiproduct_trace()
    latex = render_coefficient_semiproduct_latex(first)

    assert latex == render_coefficient_semiproduct_latex(second)
    assert r"\operatorname{Op}(r_1)M_{\widehat G_1}" in latex
    assert r"\simeq\operatorname{Op}(r_1\widehat G_1)" in latex
    assert "CERTIFIED\\_MOD\\_COMPACT" in latex
