import sympy as sp

from symbolic_operator_calculus import (
    G1,
    G2,
    R11,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    AppliedLinearCombination,
    AppliedTerm,
    KernelAnnotatedExpression,
    apply_linear_combination,
    apply_linear_combination_ordered,
    expand_ordered,
    main_expression,
)
from semantic_helpers import regularizer_rules


REQUIRED_FACTORS = (
    (Vtilde_alpha2, Wminus_21, R11, Vtilde_alpha1, Wplus_12),
    (Vtilde_alpha2, Wminus_21, R11, G1, Wplus_12),
    (G2, Wminus_21, R11, Vtilde_alpha1, Wplus_12),
    (G2, Wminus_21, R11, G1, Wplus_12),
)


def scalar_symbols():
    x, y, u, v = sp.symbols("x y u v")
    gamma1, gamma2 = sp.symbols("gamma1 gamma2")
    return x, y, u, v, gamma1, gamma2, sp.Function("f")


def expanded_main_expression():
    return expand_ordered(main_expression())


def applied_main_expression():
    x, *_rest, f = scalar_symbols()
    return apply_linear_combination_ordered(
        expanded_main_expression(),
        f(x),
        x,
        regularizer_rules(),
    )


def expected_applied_expressions():
    x, y, u, v, gamma1, gamma2, f = scalar_symbols()
    rho1 = sp.Function("rho1")
    rho2 = sp.Function("rho2")
    g1 = sp.Function("G1")
    g2 = sp.Function("G2")
    lplus = sp.Function("Lplus_12")
    lminus = sp.Function("Lminus_21")
    r11 = sp.Function("R11")

    shifted_inner = sp.Integral(lplus(gamma1 * v, y) * f(y), (y, 0, sp.oo))
    unshifted_inner = sp.Integral(lplus(v, y) * f(y), (y, 0, sp.oo))

    return (
        rho2(x)
        * sp.Integral(
            lminus(gamma2 * x, u)
            * sp.Integral(r11(u, v) * rho1(v) * shifted_inner, (v, 0, sp.oo)),
            (u, 0, sp.oo),
        ),
        rho2(x)
        * sp.Integral(
            lminus(gamma2 * x, u)
            * sp.Integral(r11(u, v) * g1(v) * unshifted_inner, (v, 0, sp.oo)),
            (u, 0, sp.oo),
        ),
        g2(x)
        * sp.Integral(
            lminus(x, u)
            * sp.Integral(r11(u, v) * rho1(v) * shifted_inner, (v, 0, sp.oo)),
            (u, 0, sp.oo),
        ),
        g2(x)
        * sp.Integral(
            lminus(x, u)
            * sp.Integral(r11(u, v) * g1(v) * unshifted_inner, (v, 0, sp.oo)),
            (u, 0, sp.oo),
        ),
    )


def expressions(applied):
    return tuple(term.expression.expression for term in applied.result_terms)


def coefficients(applied):
    return tuple(term.coefficient for term in applied.result_terms)


def has_external_factor(expression, factor):
    if isinstance(expression, KernelAnnotatedExpression):
        expression = expression.expression
    return expression == factor or (expression.is_Mul and factor in expression.args)


def test_h1_original_expression_expands_to_four_terms():
    expanded = expanded_main_expression()

    assert len(expanded.terms) == 4


def test_h2_expanded_coefficients_stay_in_required_order():
    expanded = expanded_main_expression()

    assert tuple(term.coefficient for term in expanded.terms) == (1, -1, -1, 1)


def test_h3_expanded_products_keep_exact_factor_order():
    expanded = expanded_main_expression()

    assert tuple(term.product.factors for term in expanded.terms) == REQUIRED_FACTORS


def test_h4_application_produces_four_scalar_results():
    applied = applied_main_expression()

    assert isinstance(applied, AppliedLinearCombination)
    assert len(applied.result_terms) == 4
    assert all(isinstance(term, AppliedTerm) for term in applied.result_terms)
    assert all(
        isinstance(term.expression, KernelAnnotatedExpression)
        for term in applied.result_terms
    )


def test_h5_each_applied_result_contains_three_nested_integrals():
    applied = applied_main_expression()

    assert tuple(
        len(term.expression.atoms(sp.Integral)) for term in applied.result_terms
    ) == (3, 3, 3, 3)


def test_h6_only_first_two_terms_contain_rho2_of_external_variable():
    x, *_rest = scalar_symbols()
    rho2_x = sp.Function("rho2")(x)
    applied = applied_main_expression()

    assert tuple(term.expression.has(rho2_x) for term in applied.result_terms) == (
        True,
        True,
        False,
        False,
    )


def test_h7_only_last_two_terms_contain_external_g2_factor():
    x, *_rest = scalar_symbols()
    g2_x = sp.Function("G2")(x)
    applied = applied_main_expression()

    assert tuple(
        has_external_factor(term.expression, g2_x) for term in applied.result_terms
    ) == (False, False, True, True)


def test_h8_only_terms_one_and_three_contain_rho1_on_intermediate_variable():
    _, _, _, v, *_rest = scalar_symbols()
    rho1_v = sp.Function("rho1")(v)
    applied = applied_main_expression()

    assert tuple(term.expression.has(rho1_v) for term in applied.result_terms) == (
        True,
        False,
        True,
        False,
    )


def test_h9_only_terms_two_and_four_contain_g1_on_intermediate_variable():
    _, _, _, v, *_rest = scalar_symbols()
    g1_v = sp.Function("G1")(v)
    applied = applied_main_expression()

    assert tuple(term.expression.has(g1_v) for term in applied.result_terms) == (
        False,
        True,
        False,
        True,
    )


def test_h10_only_terms_one_and_three_shift_lplus_by_gamma1():
    _, y, _, v, gamma1, *_rest = scalar_symbols()
    shifted_lplus = sp.Function("Lplus_12")(gamma1 * v, y)
    applied = applied_main_expression()

    assert tuple(
        term.expression.has(shifted_lplus) for term in applied.result_terms
    ) == (True, False, True, False)


def test_h11_only_terms_one_and_two_shift_lminus_by_gamma2():
    x, _, u, _, _, gamma2, _ = scalar_symbols()
    shifted_lminus = sp.Function("Lminus_21")(gamma2 * x, u)
    applied = applied_main_expression()

    assert tuple(
        term.expression.has(shifted_lminus) for term in applied.result_terms
    ) == (True, True, False, False)


def test_h12_signs_remain_ordered_after_application():
    applied = applied_main_expression()

    assert coefficients(applied) == (1, -1, -1, 1)
    assert tuple(term.as_expr() for term in applied.result_terms) == (
        applied.result_terms[0].expression.expression,
        -applied.result_terms[1].expression.expression,
        -applied.result_terms[2].expression.expression,
        applied.result_terms[3].expression.expression,
    )


def test_h13_ordered_structure_is_inspectable_without_sympy_add():
    applied = applied_main_expression()

    assert not isinstance(applied, sp.Add)
    assert applied.terms is applied.result_terms
    assert tuple(term.product.factors for term in applied.result_terms) == REQUIRED_FACTORS
    assert applied.result_terms[0].product.factors == REQUIRED_FACTORS[0]


def test_h14_sympy_projection_matches_sum_of_ordered_results():
    x, *_rest, f = scalar_symbols()
    expanded = expanded_main_expression()
    applied = apply_linear_combination_ordered(expanded, f(x), x, regularizer_rules())

    expected_sum = sum(
        (term.as_expr() for term in applied.result_terms),
        sp.Integer(0),
    )
    assert expressions(applied) == expected_applied_expressions()
    assert applied.as_expr() == expected_sum
    result = apply_linear_combination(expanded, f(x), x, regularizer_rules())
    assert isinstance(result, KernelAnnotatedExpression)
    assert result.expression == expected_sum
