from fractions import Fraction

import sympy as sp

from symbolic_operator_calculus import (
    G1,
    G2,
    I,
    S_Rplus,
    Vtilde_alpha1,
    Vtilde_alpha2,
    PrincipalValue,
    PrincipalValueIntegralAction,
    Product,
    a11_exact_operator,
    a22_exact_operator,
    apply_atom,
    apply_linear_combination_ordered,
    apply_product,
    collect_bound_symbols,
    mvp_atomic_rules,
    pminus_operator,
    pplus_operator,
)


HALF = Fraction(1, 2)


def variables_and_function():
    x, y = sp.symbols("x y")
    return x, y, sp.Function("f")


def principal_value_integral(operand, x, y):
    return PrincipalValue(sp.Integral(operand / (y - x), (y, 0, sp.oo))) / (
        sp.pi * sp.I
    )


def test_s_rplus_has_a_principal_value_atomic_action():
    assert isinstance(mvp_atomic_rules()[S_Rplus], PrincipalValueIntegralAction)


def test_s_rplus_action_has_exact_coefficient_denominator_and_limits():
    x, y, f = variables_and_function()

    result = apply_atom(
        S_Rplus,
        f(x),
        x,
        mvp_atomic_rules(),
        integration_variable=y,
    )
    principal_value = next(iter(result.atoms(PrincipalValue)))
    integral = principal_value.args[0]

    assert result == principal_value_integral(f(y), x, y)
    assert isinstance(principal_value, PrincipalValue)
    assert isinstance(integral, sp.Integral)
    assert sp.simplify(result / principal_value) == 1 / (sp.pi * sp.I)
    assert integral.function == f(y) / (y - x)
    assert integral.function != f(y) / (x - y)
    assert integral.limits == ((y, 0, sp.oo),)
    assert integral.variables == [y]
    assert y not in result.free_symbols
    assert collect_bound_symbols(result) == {y}


def test_operator_identity_and_imaginary_unit_are_distinct():
    x, y, f = variables_and_function()
    identity_result = apply_atom(I, f(x), x, mvp_atomic_rules())
    singular_result = apply_atom(
        S_Rplus,
        f(x),
        x,
        mvp_atomic_rules(),
        integration_variable=y,
    )

    assert I.name == "I"
    assert identity_result == f(x)
    assert singular_result.has(sp.I)
    assert not isinstance(I, sp.Expr)


def test_s_rplus_uses_a_fresh_bound_variable_in_product_application():
    x, y, v = sp.symbols("x y v")
    f, h = sp.Function("f"), sp.Function("h")

    result = apply_product(
        Vtilde_alpha1 * S_Rplus,
        f(x) + h(y),
        x,
        mvp_atomic_rules(),
    )
    integral = next(iter(result.atoms(sp.Integral)))

    assert integral.variables == [v]
    assert y in result.free_symbols
    assert v not in result.free_symbols


def test_transported_shifts_move_the_principal_value_singularity():
    x, y, f = variables_and_function()
    gamma1, gamma2 = sp.symbols("gamma1 gamma2")
    rho1, rho2 = sp.Function("rho1"), sp.Function("rho2")

    shifted1 = apply_product(
        Vtilde_alpha1 * S_Rplus,
        f(x),
        x,
        mvp_atomic_rules(),
        integration_variable=y,
    )
    shifted2 = apply_product(
        Vtilde_alpha2 * S_Rplus,
        f(x),
        x,
        mvp_atomic_rules(),
        integration_variable=y,
    )

    assert shifted1 == rho1(x) * principal_value_integral(
        f(y), gamma1 * x, y
    )
    assert shifted2 == rho2(x) * principal_value_integral(
        f(y), gamma2 * x, y
    )


def test_projection_actions_preserve_exact_ordered_contributions():
    x, y, f = variables_and_function()
    rules = mvp_atomic_rules()

    pplus = apply_linear_combination_ordered(
        pplus_operator(), f(x), x, rules, integration_variable=y
    )
    pminus = apply_linear_combination_ordered(
        pminus_operator(), f(x), x, rules, integration_variable=y
    )

    assert tuple(term.coefficient for term in pplus.terms) == (HALF, HALF)
    assert tuple(term.coefficient for term in pminus.terms) == (HALF, -HALF)
    assert tuple(term.expression for term in pplus.terms) == (
        f(x),
        principal_value_integral(f(y), x, y),
    )
    assert tuple(term.expression for term in pminus.terms) == (
        f(x),
        principal_value_integral(f(y), x, y),
    )
    assert pplus.as_expr() == sp.Rational(1, 2) * f(x) + sp.Rational(
        1, 2
    ) * principal_value_integral(f(y), x, y)
    assert pminus.as_expr() == sp.Rational(1, 2) * f(x) - sp.Rational(
        1, 2
    ) * principal_value_integral(f(y), x, y)


def expected_diagonal_terms(x, y, f, gamma, rho, multiplier):
    return (
        rho(x) * f(gamma * x),
        rho(x) * principal_value_integral(f(y), gamma * x, y),
        multiplier(x) * f(x),
        multiplier(x) * principal_value_integral(f(y), x, y),
    )


def test_complete_a11_action_preserves_four_ordered_signed_terms():
    x, y, f = variables_and_function()
    model = a11_exact_operator()
    applied = apply_linear_combination_ordered(
        model, f(x), x, mvp_atomic_rules(), integration_variable=y
    )

    assert tuple(term.coefficient for term in applied.terms) == (
        HALF,
        HALF,
        HALF,
        -HALF,
    )
    assert tuple(term.expression for term in applied.terms) == expected_diagonal_terms(
        x,
        y,
        f,
        sp.Symbol("gamma1"),
        sp.Function("rho1"),
        sp.Function("G1"),
    )


def test_complete_a22_action_preserves_four_ordered_signed_terms():
    x, y, f = variables_and_function()
    model = a22_exact_operator()
    applied = apply_linear_combination_ordered(
        model, f(x), x, mvp_atomic_rules(), integration_variable=y
    )

    assert tuple(term.coefficient for term in applied.terms) == (
        HALF,
        HALF,
        HALF,
        -HALF,
    )
    assert tuple(term.expression for term in applied.terms) == expected_diagonal_terms(
        x,
        y,
        f,
        sp.Symbol("gamma2"),
        sp.Function("rho2"),
        sp.Function("G2"),
    )


def test_no_operator_power_or_projection_rewrite_rules_were_added():
    squared = S_Rplus * S_Rplus

    assert isinstance(squared, Product)
    assert squared.factors == (S_Rplus, S_Rplus)
    assert squared != I
    assert len(pplus_operator().terms) == 2
    assert len(pminus_operator().terms) == 2
    assert (G1, G2) != (I, S_Rplus)
