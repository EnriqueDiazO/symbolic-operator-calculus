import sympy as sp

from symbolic_operator_calculus import (
    Wminus_21,
    Wplus_12,
    a12_wh_model,
    a21_wh_model,
    apply_atom,
    apply_linear_combination_ordered,
    apply_product,
    explicit_wiener_hopf_rules,
    kminus_kernel,
    kplus_kernel,
    mvp_atomic_rules,
    positive_decay_symbol,
)


def variables_and_functions():
    x, y = sp.symbols("x y")
    return x, y, sp.Function("f"), sp.Function("chi_infinity")


def normalized_plus_integral(x, y, f, chi_infinity, d):
    return sp.Integral(
        chi_infinity(x)
        * kplus_kernel(x - y, d)
        * chi_infinity(y)
        * f(y),
        (y, 0, sp.oo),
    )


def normalized_minus_integral(x, y, f, chi_infinity, d):
    return sp.Integral(
        chi_infinity(x)
        * kminus_kernel(x - y, d)
        * chi_infinity(y)
        * f(y),
        (y, 0, sp.oo),
    )


def test_default_wiener_hopf_rules_remain_formal():
    x, y, f, _ = variables_and_functions()
    rules = mvp_atomic_rules()

    plus = apply_atom(Wplus_12, f(x), x, rules, integration_variable=y)
    minus = apply_atom(Wminus_21, f(x), x, rules, integration_variable=y)

    assert plus == sp.Integral(sp.Function("Lplus_12")(x, y) * f(y), (y, 0, sp.oo))
    assert minus == sp.Integral(
        sp.Function("Lminus_21")(x, y) * f(y), (y, 0, sp.oo)
    )


def test_explicit_wiener_hopf_atoms_use_normalized_fourier_kernels():
    x, y, f, chi_infinity = variables_and_functions()
    d = positive_decay_symbol()
    rules = explicit_wiener_hopf_rules(decay=d)

    plus = apply_atom(Wplus_12, f(x), x, rules, integration_variable=y)
    minus = apply_atom(Wminus_21, f(x), x, rules, integration_variable=y)

    assert plus == normalized_plus_integral(x, y, f, chi_infinity, d)
    assert minus == normalized_minus_integral(x, y, f, chi_infinity, d)


def test_explicit_atoms_do_not_use_original_unnormalized_kernels():
    x, y, f, chi_infinity = variables_and_functions()
    d = positive_decay_symbol()
    rules = explicit_wiener_hopf_rules(decay=d)

    plus = apply_atom(Wplus_12, f(x), x, rules, integration_variable=y)
    minus = apply_atom(Wminus_21, f(x), x, rules, integration_variable=y)
    original_plus = sp.Integral(
        chi_infinity(x) * 2 * kplus_kernel(x - y, d) * chi_infinity(y) * f(y),
        (y, 0, sp.oo),
    )
    original_minus = sp.Integral(
        chi_infinity(x) * -2 * kminus_kernel(x - y, d) * chi_infinity(y) * f(y),
        (y, 0, sp.oo),
    )

    assert plus != original_plus
    assert minus != original_minus


def test_explicit_shift_changes_kernel_and_left_cutoff_for_both_actions():
    x, y, f, chi_infinity = variables_and_functions()
    d = positive_decay_symbol()
    gamma1, gamma2 = sp.symbols("gamma1 gamma2")
    rho1, rho2 = sp.Function("rho1"), sp.Function("rho2")
    rules = explicit_wiener_hopf_rules(decay=d)

    plus = apply_product(
        a12_wh_model().expression.terms[0].product,
        f(x),
        x,
        rules,
        integration_variable=y,
    )
    minus = apply_product(
        a21_wh_model().expression.terms[0].product,
        f(x),
        x,
        rules,
        integration_variable=y,
    )

    assert plus == rho1(x) * normalized_plus_integral(
        gamma1 * x, y, f, chi_infinity, d
    )
    assert minus == rho2(x) * normalized_minus_integral(
        gamma2 * x, y, f, chi_infinity, d
    )
    assert plus.has(chi_infinity(gamma1 * x))
    assert minus.has(chi_infinity(gamma2 * x))


def test_a12_wh_model_has_normalized_coefficients_and_explicit_action():
    x, y, f, chi_infinity = variables_and_functions()
    d = positive_decay_symbol()
    gamma1 = sp.Symbol("gamma1")
    rho1, g1 = sp.Function("rho1"), sp.Function("G1")
    model = a12_wh_model()

    applied = apply_linear_combination_ordered(
        model.expression,
        f(x),
        x,
        explicit_wiener_hopf_rules(decay=d),
        integration_variable=y,
    )

    assert tuple(term.coefficient for term in model.expression.terms) == (1, -1)
    assert tuple(term.coefficient for term in applied.terms) == (1, -1)
    assert applied.terms[0].expression == rho1(x) * normalized_plus_integral(
        gamma1 * x, y, f, chi_infinity, d
    )
    assert applied.terms[1].expression == g1(x) * normalized_plus_integral(
        x, y, f, chi_infinity, d
    )


def test_a21_wh_model_has_normalized_coefficients_and_explicit_action():
    x, y, f, chi_infinity = variables_and_functions()
    d = positive_decay_symbol()
    gamma2 = sp.Symbol("gamma2")
    rho2, g2 = sp.Function("rho2"), sp.Function("G2")
    model = a21_wh_model()

    applied = apply_linear_combination_ordered(
        model.expression,
        f(x),
        x,
        explicit_wiener_hopf_rules(decay=d),
        integration_variable=y,
    )

    assert tuple(term.coefficient for term in model.expression.terms) == (-1, 1)
    assert tuple(term.coefficient for term in applied.terms) == (-1, 1)
    assert applied.terms[0].expression == rho2(x) * normalized_minus_integral(
        gamma2 * x, y, f, chi_infinity, d
    )
    assert applied.terms[1].expression == g2(x) * normalized_minus_integral(
        x, y, f, chi_infinity, d
    )
