from symbolic_operator_calculus import (
    G1,
    G2,
    R11,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    OperatorAtom,
    Product,
    expand_ordered,
    main_expression,
)


def factors(term):
    return term.product.factors


def test_main_expression_uses_noncommutative_ordered_products():
    expression = main_expression()
    assert isinstance(expression, Product)

    atoms = (
        Vtilde_alpha2,
        G2,
        Wminus_21,
        R11,
        Vtilde_alpha1,
        G1,
        Wplus_12,
    )
    assert all(not atom.is_commutative for atom in atoms)
    assert R11.is_formal_regularizer

    assert Vtilde_alpha2 * G2 != G2 * Vtilde_alpha2


def test_main_expression_expands_to_exactly_four_terms():
    expanded = expand_ordered(main_expression())

    assert len(expanded.terms) == 4


def test_main_expression_expanded_signs_are_in_required_order():
    expanded = expand_ordered(main_expression())

    assert tuple(term.sign for term in expanded.terms) == (1, -1, -1, 1)


def test_main_expression_expanded_factors_are_in_required_order():
    expanded = expand_ordered(main_expression())

    assert tuple(factors(term) for term in expanded.terms) == (
        (Vtilde_alpha2, Wminus_21, R11, Vtilde_alpha1, Wplus_12),
        (Vtilde_alpha2, Wminus_21, R11, G1, Wplus_12),
        (G2, Wminus_21, R11, Vtilde_alpha1, Wplus_12),
        (G2, Wminus_21, R11, G1, Wplus_12),
    )


def test_difference_product_expands_without_sorting_terms():
    A = OperatorAtom("A")
    B = OperatorAtom("B")
    C = OperatorAtom("C")
    D = OperatorAtom("D")

    expanded = expand_ordered((A - B) * (C - D))

    assert tuple(term.sign for term in expanded.terms) == (1, -1, -1, 1)
    assert tuple(factors(term) for term in expanded.terms) == (
        (A, C),
        (A, D),
        (B, C),
        (B, D),
    )


def test_reversed_products_are_structurally_distinct():
    A = OperatorAtom("A")
    B = OperatorAtom("B")
    C = OperatorAtom("C")

    assert A * B * C != C * B * A
