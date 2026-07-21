from symbolic_operator_calculus import (
    Ghat1,
    P_minus,
    P_plus,
    R1,
    U1,
    U1_inverse,
    Wplus_12,
    Product,
    build_first_pivot_closure_analysis,
    closure_words,
    render_closure_words_latex,
    render_product_latex,
)


def test_four_right_core_words_are_exact_existing_ast_products():
    words = closure_words()

    assert tuple(word.identifier for word in words) == (
        "K++",
        "K+Ghat",
        "K-+",
        "K-Ghat",
    )
    assert tuple(word.product.factors for word in words) == (
        (U1_inverse, P_plus, R1, U1, Wplus_12),
        (U1_inverse, P_plus, R1, Ghat1, Wplus_12),
        (P_minus, R1, U1, Wplus_12),
        (P_minus, R1, Ghat1, Wplus_12),
    )
    assert all(isinstance(word.product, Product) for word in words)


def test_every_word_keeps_one_atomic_r1_and_wplus_at_the_right_edge():
    for word in closure_words():
        assert word.product.factors.count(R1) == 1
        assert word.product.factors[-1] is Wplus_12
        assert tuple(item.factor for item in word.signature) == word.product.factors
        assert tuple(item.position for item in word.signature) == tuple(
            range(1, len(word.product.factors) + 1)
        )


def test_interfaces_are_contiguous_slices_and_do_not_reorder_r1():
    for word in closure_words():
        for interface in word.interfaces:
            assert word.product.factors[
                interface.start_position - 1 : interface.end_position
            ] == interface.factors
        assert word.interfaces[0].factors[-1] is R1
        assert word.interfaces[1].factors[0] is R1


def test_latex_is_deterministic_and_preserves_factor_order():
    first = closure_words()
    second = closure_words()

    assert render_closure_words_latex(first) == render_closure_words_latex(second)
    assert tuple(render_product_latex(word.product) for word in first) == (
        r"U_{1}^{-1}\,P^{+}\,R_{1}\,U_{1}\,W^+_{1,2}",
        r"U_{1}^{-1}\,P^{+}\,R_{1}\,\widehat G_{1}\,W^+_{1,2}",
        r"P^{-}\,R_{1}\,U_{1}\,W^+_{1,2}",
        r"P^{-}\,R_{1}\,\widehat G_{1}\,W^+_{1,2}",
    )


def test_complete_closure_analysis_is_deterministic():
    assert build_first_pivot_closure_analysis() == build_first_pivot_closure_analysis()
