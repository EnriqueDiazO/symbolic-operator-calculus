import sympy as sp

from symbolic_operator_calculus import (
    Ghat1,
    P_minus,
    P_plus,
    R1,
    U1,
    U1_inverse,
    CertificationStatus,
    DerivationRelationKind,
    FactorizedOscillatoryMellinSymbol,
    Product,
    RuleCertificationStatus,
    StandardMellinMembership,
    build_first_pivot_left_core_closure,
    render_first_pivot_left_cores_latex,
)


def _trace():
    return build_first_pivot_left_core_closure(*sp.symbols("lambda gamma x"))


def test_four_left_core_words_are_exact_ordered_ast_products_with_atomic_r1():
    trace = _trace()
    expected = (
        ("L-+", Product((P_minus, R1, U1))),
        ("L++", Product((U1_inverse, P_plus, R1, U1))),
        ("L-Ghat", Product((P_minus, R1, Ghat1))),
        ("L+Ghat", Product((U1_inverse, P_plus, R1, Ghat1))),
    )

    assert tuple((core.identifier, core.ast_product) for core in trace.cores) == expected
    assert all(core.ast_product.factors.count(R1) == 1 for core in trace.cores)
    assert R1.is_formal_regularizer


def test_cauchy_semiproducts_are_certified_only_modulo_compacts():
    trace = _trace()

    assert tuple(relation.left for relation in trace.cauchy_relations) == (
        Product((P_plus, R1)),
        Product((P_minus, R1)),
    )
    assert all(
        relation.certification_status is CertificationStatus.EVIDENCE_SUPPLIED
        for relation in trace.cauchy_relations
    )
    assert all(
        rule.relation_kind is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
        and rule.certification_status is RuleCertificationStatus.CERTIFIED_MOD_COMPACT
        for rule in trace.cauchy_rules
    )
    assert trace.cauchy_relations[0].evidence.printed_page == "942"


def test_two_outputs_are_standard_and_two_remain_right_factorized():
    cores = {core.identifier: core for core in _trace().cores}

    assert cores["L++"].output_kind == "STANDARD_MELLIN_PDO"
    assert cores["L-Ghat"].output_kind == "STANDARD_MELLIN_PDO"
    assert cores["L++"].output_symbol.membership is StandardMellinMembership.CERTIFIED
    assert cores["L-Ghat"].output_symbol.membership is StandardMellinMembership.CERTIFIED
    for identifier in ("L-+", "L+Ghat"):
        output = cores[identifier].output_symbol
        assert isinstance(output, FactorizedOscillatoryMellinSymbol)
        assert output.standard_product_membership is StandardMellinMembership.NOT_DEMONSTRATED
    assert cores["L-+"].output_symbol.oscillatory_factor.expression == sp.Pow(
        sp.Symbol("gamma"), sp.I * sp.Symbol("lambda"), evaluate=False
    )
    assert cores["L+Ghat"].output_symbol.oscillatory_factor.expression == sp.Pow(
        sp.Symbol("gamma"), -sp.I * sp.Symbol("lambda"), evaluate=False
    )


def test_all_final_relations_are_mod_compact_and_residues_stay_compact():
    trace = _trace()

    assert all(
        core.rule.relation_kind is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
        for core in trace.cores
    )
    assert all(core.compact_residual_terms for core in trace.cores)
    proof = " ".join(trace.compact_ideal_proof)
    assert "AKB is compact" in proof
    assert "finite sums" in proof
    assert "conjugation preserves compactness" in proof


def test_latex_is_deterministic_ordered_and_excludes_the_wiener_hopf_factor():
    first = _trace()
    second = _trace()
    latex = render_first_pivot_left_cores_latex(first)

    assert latex == render_first_pivot_left_cores_latex(second)
    assert latex.count(r"\simeq") == 4
    assert r"U_{1}^{-1}\,P^{+}\,R_{1}\,U_{1}" in latex
    assert r"p^+(\lambda)r_1(x/\gamma_1,\lambda)" in latex
    assert r"d_{\gamma_1^{-1}}" in latex
    assert "W_{1,2}" not in latex
    assert "P^+P^-=0" not in latex
