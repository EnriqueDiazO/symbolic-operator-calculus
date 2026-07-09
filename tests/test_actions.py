import pytest
import sympy as sp

from symbolic_operator_calculus import (
    G1,
    G2,
    R11,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    FormalRegularizerAction,
    IntegralKernelAction,
    LinearCombination,
    LinearCombinationApplicationError,
    MissingAtomicActionError,
    MissingIntegrationVariableError,
    OperatorAtom,
    Product,
    ProductApplicationError,
    UnsupportedAtomicExpressionError,
    UnsupportedProductExpressionError,
    Term,
    apply_atom,
    apply_linear_combination,
    apply_product,
    main_expression,
    mvp_atomic_rules,
)


def symbols_and_functions():
    x, y = sp.symbols("x y")
    f = sp.Function("f")
    h = sp.Function("h")
    return x, y, f, h


def lc(*terms):
    return LinearCombination(
        tuple(Term(coefficient, product) for coefficient, product in terms)
    )


def integral_count(expression):
    return len(expression.atoms(sp.Integral))


def test_multiplication_action_for_g1_on_function():
    x, _, f, _ = symbols_and_functions()

    result = apply_atom(G1, f(x), x, mvp_atomic_rules())

    assert result == sp.Function("G1")(x) * f(x)


def test_multiplication_action_for_g2_on_general_expression():
    x, _, _, h = symbols_and_functions()

    result = apply_atom(G2, h(x), x, mvp_atomic_rules())

    assert result == sp.Function("G2")(x) * h(x)


def test_transported_shift_alpha1_on_function():
    x, _, f, _ = symbols_and_functions()
    gamma1 = sp.Symbol("gamma1")

    result = apply_atom(Vtilde_alpha1, f(x), x, mvp_atomic_rules())

    assert result == sp.Function("rho1")(x) * f(gamma1 * x)


def test_transported_shift_alpha2_on_function():
    x, _, f, _ = symbols_and_functions()
    gamma2 = sp.Symbol("gamma2")

    result = apply_atom(Vtilde_alpha2, f(x), x, mvp_atomic_rules())

    assert result == sp.Function("rho2")(x) * f(gamma2 * x)


def test_transported_shift_alpha1_on_integral_operand():
    x, y, f, _ = symbols_and_functions()
    gamma1 = sp.Symbol("gamma1")
    K = sp.Function("K")
    operand = sp.Integral(K(x, y) * f(y), (y, 0, sp.oo))

    result = apply_atom(Vtilde_alpha1, operand, x, mvp_atomic_rules())

    expected = sp.Function("rho1")(x) * sp.Integral(
        K(gamma1 * x, y) * f(y),
        (y, 0, sp.oo),
    )
    assert result == expected


def test_transported_shift_actions_remain_structurally_distinct():
    x, _, f, _ = symbols_and_functions()
    rules = mvp_atomic_rules()

    alpha1_result = apply_atom(Vtilde_alpha1, f(x), x, rules)
    alpha2_result = apply_atom(Vtilde_alpha2, f(x), x, rules)

    assert alpha1_result != alpha2_result
    assert rules[Vtilde_alpha1] != rules[Vtilde_alpha2]


def test_wiener_hopf_positive_action_is_formal_integral():
    x, y, f, _ = symbols_and_functions()

    result = apply_atom(
        Wplus_12,
        f(x),
        x,
        mvp_atomic_rules(),
        integration_variable=y,
    )

    expected = sp.Integral(sp.Function("Lplus_12")(x, y) * f(y), (y, 0, sp.oo))
    assert result == expected


def test_wiener_hopf_negative_action_is_formal_integral():
    x, y, f, _ = symbols_and_functions()

    result = apply_atom(
        Wminus_21,
        f(x),
        x,
        mvp_atomic_rules(),
        integration_variable=y,
    )

    expected = sp.Integral(sp.Function("Lminus_21")(x, y) * f(y), (y, 0, sp.oo))
    assert result == expected


def test_formal_regularizer_action_is_formal_integral():
    x, y, f, _ = symbols_and_functions()

    result = apply_atom(
        R11,
        f(x),
        x,
        mvp_atomic_rules(),
        integration_variable=y,
    )

    expected = sp.Integral(sp.Function("R11")(x, y) * f(y), (y, 0, sp.oo))
    assert result == expected


def test_formal_regularizer_action_is_semantically_distinct():
    rules = mvp_atomic_rules()

    assert isinstance(rules[R11], FormalRegularizerAction)
    assert rules[R11].is_formal_regularizer
    assert isinstance(rules[Wplus_12], IntegralKernelAction)
    assert not isinstance(rules[Wplus_12], FormalRegularizerAction)


def test_unregistered_atom_fails_with_specific_error():
    x, _, f, _ = symbols_and_functions()
    unknown = OperatorAtom("Unknown")

    with pytest.raises(MissingAtomicActionError):
        apply_atom(unknown, f(x), x, mvp_atomic_rules())


def test_apply_atom_rejects_product_instead_of_composing():
    x, _, f, _ = symbols_and_functions()

    with pytest.raises(UnsupportedAtomicExpressionError):
        apply_atom(G1 * G2, f(x), x, mvp_atomic_rules())


def test_apply_product_two_atoms_uses_rightmost_factor_first():
    x, _, f, _ = symbols_and_functions()
    gamma1 = sp.Symbol("gamma1")

    result = apply_product(Vtilde_alpha1 * G1, f(x), x, mvp_atomic_rules())

    expected = sp.Function("rho1")(x) * sp.Function("G1")(gamma1 * x) * f(
        gamma1 * x
    )
    assert result == expected


def test_apply_product_three_atoms_uses_composition_order():
    x, _, f, _ = symbols_and_functions()
    gamma1 = sp.Symbol("gamma1")

    result = apply_product(G2 * Vtilde_alpha1 * G1, f(x), x, mvp_atomic_rules())

    expected = (
        sp.Function("G2")(x)
        * sp.Function("rho1")(x)
        * sp.Function("G1")(gamma1 * x)
        * f(gamma1 * x)
    )
    assert result == expected


def test_apply_product_order_is_not_reversed():
    x, _, f, _ = symbols_and_functions()

    left_then_right = apply_product(Vtilde_alpha1 * G1, f(x), x, mvp_atomic_rules())
    reversed_order = apply_product(G1 * Vtilde_alpha1, f(x), x, mvp_atomic_rules())

    assert left_then_right != reversed_order


def test_apply_product_one_factor_matches_apply_atom():
    x, _, f, _ = symbols_and_functions()

    result = apply_product(Product((G1,)), f(x), x, mvp_atomic_rules())

    assert result == apply_atom(G1, f(x), x, mvp_atomic_rules())


def test_apply_product_empty_product_returns_operand_unchanged():
    x, _, f, _ = symbols_and_functions()
    operand = f(x)

    result = apply_product(Product(()), operand, x, mvp_atomic_rules())

    assert result is operand


def test_apply_product_rejects_non_product_input():
    x, _, f, _ = symbols_and_functions()

    with pytest.raises(UnsupportedProductExpressionError):
        apply_product(G1, f(x), x, mvp_atomic_rules())


def test_apply_product_wraps_unregistered_atom_error():
    x, _, f, _ = symbols_and_functions()
    unknown = OperatorAtom("Unknown")
    product = Product((unknown,))

    with pytest.raises(ProductApplicationError) as exc_info:
        apply_product(product, f(x), x, mvp_atomic_rules())

    assert isinstance(exc_info.value.__cause__, MissingAtomicActionError)
    assert exc_info.value.product == product
    assert exc_info.value.factor == unknown
    assert exc_info.value.factor_position == 1
    assert exc_info.value.factor_count == 1


def test_apply_product_wraps_apply_atom_error():
    x, _, f, _ = symbols_and_functions()
    product = Product((Wplus_12,))

    with pytest.raises(ProductApplicationError) as exc_info:
        apply_product(product, f(x), x, mvp_atomic_rules())

    assert isinstance(exc_info.value.__cause__, MissingIntegrationVariableError)


def test_apply_product_error_context_identifies_failed_factor_position():
    x, _, f, _ = symbols_and_functions()
    unknown = OperatorAtom("Unknown")
    product = G1 * unknown * G2

    with pytest.raises(ProductApplicationError) as exc_info:
        apply_product(product, f(x), x, mvp_atomic_rules())

    assert exc_info.value.product == product
    assert exc_info.value.factor == unknown
    assert exc_info.value.factor_position == 2
    assert exc_info.value.factor_count == 3
    assert isinstance(exc_info.value.__cause__, MissingAtomicActionError)


def test_apply_product_accepts_valid_sympy_scalar_operand():
    x, _, f, _ = symbols_and_functions()
    operand = x**2 + f(x)

    result = apply_product(Product((G1,)), operand, x, mvp_atomic_rules())

    assert result == sp.Function("G1")(x) * (x**2 + f(x))


def test_apply_product_handles_nontrivial_scalar_expression():
    x, _, f, _ = symbols_and_functions()
    operand = x + f(x)

    result = apply_product(G2 * G1, operand, x, mvp_atomic_rules())

    assert result == sp.Function("G2")(x) * sp.Function("G1")(x) * (x + f(x))


def test_apply_product_allows_first_applied_factor_to_generate_simple_expression():
    x, _, f, _ = symbols_and_functions()

    result = apply_product(G1 * G2, f(x), x, mvp_atomic_rules())

    assert result == sp.Function("G1")(x) * sp.Function("G2")(x) * f(x)


def test_apply_product_composes_shift_after_integral_kernel():
    x, y, f, _ = symbols_and_functions()
    gamma1 = sp.Symbol("gamma1")
    product = Vtilde_alpha1 * Wplus_12

    result = apply_product(product, f(x), x, mvp_atomic_rules())

    expected = sp.Function("rho1")(x) * sp.Integral(
        sp.Function("Lplus_12")(gamma1 * x, y) * f(y),
        (y, 0, sp.oo),
    )
    assert result == expected
    assert integral_count(result) == 1


def test_apply_product_composes_regularizer_after_shifted_kernel():
    x, y, f, _ = symbols_and_functions()
    v = sp.Symbol("v")
    gamma1 = sp.Symbol("gamma1")
    product = R11 * Vtilde_alpha1 * Wplus_12

    result = apply_product(product, f(x), x, mvp_atomic_rules())

    expected = sp.Integral(
        sp.Function("R11")(x, v)
        * sp.Function("rho1")(v)
        * sp.Integral(
            sp.Function("Lplus_12")(gamma1 * v, y) * f(y),
            (y, 0, sp.oo),
        ),
        (v, 0, sp.oo),
    )
    assert result == expected
    assert integral_count(result) == 2


def test_mvp_atomic_rules_mapping_is_stable():
    rules = mvp_atomic_rules()

    with pytest.raises(TypeError):
        rules[G1] = rules[G1]


def test_apply_product_does_not_mutate_product():
    x, _, f, _ = symbols_and_functions()
    product = G2 * G1
    original_factors = product.factors

    apply_product(product, f(x), x, mvp_atomic_rules())

    assert product.factors == original_factors


def test_apply_product_does_not_mutate_original_operand():
    x, _, f, _ = symbols_and_functions()
    operand = x + f(x)
    original_operand = operand

    apply_product(G2 * G1, operand, x, mvp_atomic_rules())

    assert operand == original_operand


def test_main_term_product_produces_three_nested_integrals():
    x, y, f, _ = symbols_and_functions()
    u, v = sp.symbols("u v")
    gamma1 = sp.Symbol("gamma1")
    gamma2 = sp.Symbol("gamma2")
    product = main_expression().terms[0].product

    result = apply_product(product, f(x), x, mvp_atomic_rules())

    expected = sp.Function("rho2")(x) * sp.Integral(
        sp.Function("Lminus_21")(gamma2 * x, u)
        * sp.Integral(
            sp.Function("R11")(u, v)
            * sp.Function("rho1")(v)
            * sp.Integral(
                sp.Function("Lplus_12")(gamma1 * v, y) * f(y),
                (y, 0, sp.oo),
            ),
            (v, 0, sp.oo),
        ),
        (u, 0, sp.oo),
    )
    assert result == expected
    assert integral_count(result) == 3


def test_main_term_rho1_uses_regularizer_integration_variable():
    x, _, f, _ = symbols_and_functions()
    v = sp.Symbol("v")
    result = apply_product(main_expression().terms[0].product, f(x), x, mvp_atomic_rules())

    assert result.has(sp.Function("rho1")(v))
    assert not result.has(sp.Function("rho1")(x))


def test_main_term_rho2_uses_external_variable():
    x, _, f, _ = symbols_and_functions()
    u, v, y = sp.symbols("u v y")
    result = apply_product(main_expression().terms[0].product, f(x), x, mvp_atomic_rules())

    assert result.has(sp.Function("rho2")(x))
    assert not result.has(sp.Function("rho2")(u))
    assert not result.has(sp.Function("rho2")(v))
    assert not result.has(sp.Function("rho2")(y))


def test_main_term_gamma1_only_shifts_lplus_output_variable():
    x, y, f, _ = symbols_and_functions()
    v = sp.Symbol("v")
    gamma1 = sp.Symbol("gamma1")
    result = apply_product(main_expression().terms[0].product, f(x), x, mvp_atomic_rules())

    assert result.has(sp.Function("Lplus_12")(gamma1 * v, y))
    assert not result.has(sp.Function("Lplus_12")(v, y))
    assert not result.has(f(gamma1 * y))


def test_main_term_gamma2_only_shifts_lminus_output_variable():
    x, _, f, _ = symbols_and_functions()
    u = sp.Symbol("u")
    gamma2 = sp.Symbol("gamma2")
    result = apply_product(main_expression().terms[0].product, f(x), x, mvp_atomic_rules())

    assert result.has(sp.Function("Lminus_21")(gamma2 * x, u))
    assert not result.has(sp.Function("Lminus_21")(x, u))
    assert not result.has(sp.Function("Lminus_21")(gamma2 * x, gamma2 * u))


def test_main_term_bound_variables_are_distinct_and_not_captured():
    x, _, f, _ = symbols_and_functions()
    y, u, v = sp.symbols("y u v")
    result = apply_product(main_expression().terms[0].product, f(x), x, mvp_atomic_rules())

    bound_variables = {
        limit[0]
        for integral in result.atoms(sp.Integral)
        for limit in integral.limits
    }
    assert bound_variables == {y, u, v}
    assert bound_variables.isdisjoint(result.free_symbols)


def test_main_term_product_does_not_mutate_original_product():
    x, _, f, _ = symbols_and_functions()
    product = main_expression().terms[0].product
    original_factors = product.factors

    apply_product(product, f(x), x, mvp_atomic_rules())

    assert product.factors == original_factors


def test_apply_linear_combination_two_atoms():
    x, _, f, _ = symbols_and_functions()
    combination = lc((1, G1), (1, G2))

    result = apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert result == sp.Function("G1")(x) * f(x) + sp.Function("G2")(x) * f(x)


def test_apply_linear_combination_three_terms():
    x, _, f, _ = symbols_and_functions()
    combination = lc((1, G1), (-1, G2), (1, G1 * G2))

    result = apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert result == (
        sp.Function("G1")(x) * f(x)
        - sp.Function("G2")(x) * f(x)
        + sp.Function("G1")(x) * sp.Function("G2")(x) * f(x)
    )


def test_apply_linear_combination_uses_distinct_supported_coefficients():
    x, _, f, _ = symbols_and_functions()
    combination = lc((2, G1), (-3, G2))

    result = apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert tuple(term.coefficient for term in combination.terms) == (2, -3)
    assert result == (
        2 * sp.Function("G1")(x) * f(x)
        - 3 * sp.Function("G2")(x) * f(x)
    )


def test_apply_linear_combination_negative_coefficient():
    x, _, f, _ = symbols_and_functions()
    combination = lc((-3, G1))

    result = apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert result == -3 * sp.Function("G1")(x) * f(x)


def test_apply_linear_combination_single_term():
    x, _, f, _ = symbols_and_functions()
    combination = lc((1, G1))

    result = apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert result == sp.Function("G1")(x) * f(x)


def test_apply_linear_combination_empty_combination_is_zero():
    x, _, f, _ = symbols_and_functions()

    result = apply_linear_combination(
        LinearCombination(()),
        f(x),
        x,
        mvp_atomic_rules(),
    )

    assert result == sp.Integer(0)


def test_apply_linear_combination_operator_atom_term_matches_apply_atom():
    x, _, f, _ = symbols_and_functions()
    combination = lc((1, G1))

    result = apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert result == apply_atom(G1, f(x), x, mvp_atomic_rules())


def test_apply_linear_combination_product_term_matches_apply_product():
    x, _, f, _ = symbols_and_functions()
    product = G2 * G1
    combination = lc((1, product))

    result = apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert result == apply_product(product, f(x), x, mvp_atomic_rules())


def test_apply_linear_combination_mixes_atom_and_product_terms():
    x, _, f, _ = symbols_and_functions()
    combination = lc((2, G1), (-3, G2 * G1))

    result = apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert result == (
        2 * sp.Function("G1")(x) * f(x)
        - 3 * sp.Function("G2")(x) * sp.Function("G1")(x) * f(x)
    )


def test_apply_linear_combination_preserves_product_order_inside_term():
    x, _, f, _ = symbols_and_functions()
    gamma1 = sp.Symbol("gamma1")
    combination = lc((1, Vtilde_alpha1 * G1))

    result = apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert result == (
        sp.Function("rho1")(x)
        * sp.Function("G1")(gamma1 * x)
        * f(gamma1 * x)
    )


def test_apply_linear_combination_wraps_unregistered_atom_term():
    x, _, f, _ = symbols_and_functions()
    unknown = OperatorAtom("Unknown")
    combination = lc((1, unknown))

    with pytest.raises(LinearCombinationApplicationError) as exc_info:
        apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert isinstance(exc_info.value.__cause__, MissingAtomicActionError)


def test_apply_linear_combination_wraps_error_inside_product_term():
    x, _, f, _ = symbols_and_functions()
    unknown = OperatorAtom("Unknown")
    combination = lc((1, G1 * unknown))

    with pytest.raises(LinearCombinationApplicationError) as exc_info:
        apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert isinstance(exc_info.value.__cause__, ProductApplicationError)


def test_apply_linear_combination_error_context_identifies_failed_term():
    x, _, f, _ = symbols_and_functions()
    unknown = OperatorAtom("Unknown")
    combination = lc((1, G1), (-1, unknown), (1, G2))

    with pytest.raises(LinearCombinationApplicationError) as exc_info:
        apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert exc_info.value.combination == combination
    assert exc_info.value.term == combination.terms[1]
    assert exc_info.value.term_position == 2
    assert exc_info.value.term_count == 3
    assert exc_info.value.coefficient == -1


def test_apply_linear_combination_preserves_original_error_cause_chain():
    x, _, f, _ = symbols_and_functions()
    unknown = OperatorAtom("Unknown")
    combination = lc((1, G1 * unknown))

    with pytest.raises(LinearCombinationApplicationError) as exc_info:
        apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert isinstance(exc_info.value.__cause__, ProductApplicationError)
    assert isinstance(exc_info.value.__cause__.__cause__, MissingAtomicActionError)


def test_apply_linear_combination_wraps_missing_variable_for_integral_atom():
    x, _, f, _ = symbols_and_functions()
    combination = lc((1, Wplus_12))

    with pytest.raises(LinearCombinationApplicationError) as exc_info:
        apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert isinstance(exc_info.value.__cause__, MissingIntegrationVariableError)


def test_apply_linear_combination_nontrivial_scalar_operand():
    x, _, f, _ = symbols_and_functions()
    operand = x**2 + f(x)
    combination = lc((1, G1), (-1, G2))

    result = apply_linear_combination(combination, operand, x, mvp_atomic_rules())

    assert result == (
        sp.Function("G1")(x) * (x**2 + f(x))
        - sp.Function("G2")(x) * (x**2 + f(x))
    )


def test_apply_linear_combination_does_not_mutate_combination():
    x, _, f, _ = symbols_and_functions()
    combination = lc((1, G1), (-1, G2))
    original_terms = combination.terms

    apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert combination.terms == original_terms


def test_apply_linear_combination_does_not_mutate_operand():
    x, _, f, _ = symbols_and_functions()
    operand = x + f(x)
    original_operand = operand

    apply_linear_combination(lc((1, G1), (-1, G2)), operand, x, mvp_atomic_rules())

    assert operand == original_operand


def test_apply_linear_combination_preserves_apply_atom_stability():
    x, _, f, _ = symbols_and_functions()

    apply_linear_combination(lc((1, G1), (-1, G2)), f(x), x, mvp_atomic_rules())

    assert apply_atom(G1, f(x), x, mvp_atomic_rules()) == sp.Function("G1")(x) * f(x)


def test_apply_linear_combination_preserves_apply_product_stability():
    x, _, f, _ = symbols_and_functions()
    product = G2 * G1

    apply_linear_combination(lc((1, G1), (-1, G2)), f(x), x, mvp_atomic_rules())

    assert apply_product(product, f(x), x, mvp_atomic_rules()) == (
        sp.Function("G2")(x) * sp.Function("G1")(x) * f(x)
    )


def test_apply_linear_combination_applies_shifted_integral_product():
    x, y, f, _ = symbols_and_functions()
    gamma1 = sp.Symbol("gamma1")
    combination = lc((1, Vtilde_alpha1 * Wplus_12))

    result = apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    expected = sp.Function("rho1")(x) * sp.Integral(
        sp.Function("Lplus_12")(gamma1 * x, y) * f(y),
        (y, 0, sp.oo),
    )
    assert result == expected
