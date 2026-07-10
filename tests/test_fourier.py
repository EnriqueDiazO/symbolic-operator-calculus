import sympy as sp

from symbolic_operator_calculus import (
    bminus_symbol,
    bplus_symbol,
    forward_exponential,
    frequency_variable,
    inverse_exponential,
    inverse_prefactor,
    kminus_kernel,
    kplus_kernel,
    localized_lminus_kernel,
    localized_lplus_kernel,
    negative_inverse_integral,
    positive_decay_symbol,
    positive_inverse_integral,
    time_variable,
)


def algebraically_equivalent(left, right):
    return sp.simplify(left - right) == 0


def test_j1_inverse_convention_uses_positive_exponential():
    lambda_ = frequency_variable()
    t = time_variable()

    assert inverse_exponential(t, lambda_) == sp.exp(sp.I * t * lambda_)
    assert forward_exponential(t, lambda_) == sp.exp(-sp.I * t * lambda_)


def test_j2_inverse_prefactor_is_one_over_two_pi():
    assert inverse_prefactor() == sp.Rational(1, 2) / sp.pi


def test_j3_decay_parameter_is_real_and_positive():
    d = positive_decay_symbol()

    assert d.is_real
    assert d.is_positive


def test_j4_positive_integral_is_constructed_before_evaluation():
    calculation = positive_inverse_integral()

    assert isinstance(calculation.integral, sp.Integral)
    assert calculation.integral.function == calculation.integrand
    assert calculation.integral != calculation.evaluate()


def test_j5_positive_integral_limits_are_zero_to_infinity():
    calculation = positive_inverse_integral()

    assert calculation.limits == (frequency_variable(), 0, sp.oo)
    assert calculation.integral.limits == ((frequency_variable(), 0, sp.oo),)


def test_j6_positive_integral_evaluates_to_expected_expression():
    d = positive_decay_symbol()
    t = time_variable()
    calculation = positive_inverse_integral(t, d)

    assert isinstance(calculation.evaluate_default(), sp.Piecewise)
    assert calculation.assumptions["decay_parameter_is_positive"] is True
    assert algebraically_equivalent(calculation.evaluate(), 1 / (d - sp.I * t))


def test_j7_kplus_evaluates_to_expected_kernel():
    d = positive_decay_symbol()
    t = time_variable()

    assert algebraically_equivalent(
        kplus_kernel(t, d),
        1 / (2 * sp.pi * (d - sp.I * t)),
    )


def test_j8_negative_integral_is_constructed_before_evaluation():
    calculation = negative_inverse_integral()

    assert isinstance(calculation.integral, sp.Integral)
    assert calculation.integral.function == calculation.integrand
    assert calculation.integral != calculation.evaluate()


def test_j9_negative_integral_limits_are_negative_infinity_to_zero():
    calculation = negative_inverse_integral()

    assert calculation.limits == (frequency_variable(), -sp.oo, 0)
    assert calculation.integral.limits == ((frequency_variable(), -sp.oo, 0),)


def test_j10_negative_integral_evaluates_to_expected_expression():
    d = positive_decay_symbol()
    t = time_variable()
    calculation = negative_inverse_integral(t, d)

    assert isinstance(calculation.evaluate_default(), sp.Piecewise)
    assert calculation.assumptions["decay_parameter_is_positive"] is True
    assert algebraically_equivalent(calculation.evaluate(), 1 / (d + sp.I * t))


def test_j11_kminus_evaluates_to_expected_kernel():
    d = positive_decay_symbol()
    t = time_variable()

    assert algebraically_equivalent(
        kminus_kernel(t, d),
        1 / (2 * sp.pi * (d + sp.I * t)),
    )


def test_j12_equivalence_checks_are_algebraic_not_string_based():
    d = positive_decay_symbol()
    t = time_variable()
    positive_result = positive_inverse_integral(t, d).evaluate()
    negative_result = negative_inverse_integral(t, d).evaluate()

    assert algebraically_equivalent(positive_result, 1 / (d - sp.I * t))
    assert algebraically_equivalent(negative_result, 1 / (d + sp.I * t))
    assert not isinstance(positive_result, str)
    assert not isinstance(negative_result, str)


def test_j13_kplus_and_kminus_are_structurally_distinct():
    d = positive_decay_symbol()
    t = time_variable()

    assert kplus_kernel(t, d) != kminus_kernel(t, d)


def test_j14_localized_lplus_uses_kplus_of_x_minus_y():
    x, y = sp.symbols("x y", real=True)
    d = positive_decay_symbol()
    chi_infinity = sp.Function("chi_infinity")

    result = localized_lplus_kernel(x, y, decay=d)

    assert result == chi_infinity(x) * kplus_kernel(x - y, d) * chi_infinity(y)


def test_j15_localized_lminus_uses_kminus_of_x_minus_y():
    x, y = sp.symbols("x y", real=True)
    d = positive_decay_symbol()
    chi_infinity = sp.Function("chi_infinity")

    result = localized_lminus_kernel(x, y, decay=d)

    assert result == chi_infinity(x) * kminus_kernel(x - y, d) * chi_infinity(y)


def test_j16_chi_infinity_remains_formal():
    x, y = sp.symbols("x y", real=True)
    lplus = localized_lplus_kernel(x, y)
    lminus = localized_lminus_kernel(x, y)
    chi_infinity = sp.Function("chi_infinity")

    assert lplus.has(chi_infinity(x), chi_infinity(y))
    assert lminus.has(chi_infinity(x), chi_infinity(y))
    assert not lplus.has(sp.Heaviside)
    assert not lminus.has(sp.Heaviside)


def test_j17_b_symbols_use_half_line_weights_without_heaviside():
    lambda_ = frequency_variable()
    d = positive_decay_symbol()

    assert bplus_symbol(d, lambda_) == sp.exp(-d * lambda_)
    assert bminus_symbol(d, lambda_) == sp.exp(d * lambda_)
    assert not bplus_symbol(d, lambda_).has(sp.Heaviside)
    assert not bminus_symbol(d, lambda_).has(sp.Heaviside)
