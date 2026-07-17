import pytest
import sympy as sp

from symbolic_operator_calculus import (
    FourierEvaluationError,
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


def test_original_b12_kernel_is_twice_the_normalized_positive_kernel():
    d = positive_decay_symbol()
    t = time_variable()
    original_b12_kernel = 2 * kplus_kernel(t, d)

    assert algebraically_equivalent(
        original_b12_kernel,
        1 / (sp.pi * (d - sp.I * t)),
    )


def test_original_b21_kernel_is_minus_twice_the_normalized_negative_kernel():
    d = positive_decay_symbol()
    t = time_variable()
    original_b21_kernel = -2 * kminus_kernel(t, d)

    assert algebraically_equivalent(
        original_b21_kernel,
        -1 / (sp.pi * (d + sp.I * t)),
    )


def test_decay_default_is_used_only_when_argument_is_none():
    lambda_ = frequency_variable()
    default = positive_decay_symbol()

    assert bplus_symbol(None, lambda_) == sp.exp(-default * lambda_)
    assert bminus_symbol(None, lambda_) == sp.exp(default * lambda_)
    assert positive_inverse_integral(decay=None).decay_parameter == default
    assert negative_inverse_integral(decay=None).decay_parameter == default
    assert default.is_real is True
    assert default.is_positive is True


@pytest.mark.parametrize(
    "decay",
    [
        1,
        2.5,
        sp.Integer(2),
        sp.Rational(3, 2),
        sp.Symbol("d", positive=True),
        sp.exp(sp.Symbol("r", real=True)),
    ],
)
def test_positive_explicit_decay_is_preserved(decay):
    lambda_ = frequency_variable()
    expected_decay = sp.sympify(decay)

    assert bplus_symbol(decay, lambda_) == sp.exp(-expected_decay * lambda_)
    assert bminus_symbol(decay, lambda_) == sp.exp(expected_decay * lambda_)
    assert positive_inverse_integral(decay=decay).decay_parameter == expected_decay
    assert negative_inverse_integral(decay=decay).decay_parameter == expected_decay


def test_explicit_positive_symbol_is_not_replaced_by_default():
    lambda_ = frequency_variable()
    explicit_decay = sp.Symbol("delta", positive=True)
    result = bplus_symbol(explicit_decay, lambda_)

    assert result.has(explicit_decay)
    assert not result.has(positive_decay_symbol())


@pytest.mark.parametrize("decay", [0, 0.0, sp.Integer(0)])
@pytest.mark.parametrize(
    "consumer",
    [
        bplus_symbol,
        bminus_symbol,
        lambda value: positive_inverse_integral(decay=value),
        lambda value: negative_inverse_integral(decay=value),
    ],
)
def test_zero_decay_is_rejected_by_public_fourier_constructors(decay, consumer):
    with pytest.raises(FourierEvaluationError):
        consumer(decay)


def test_explicit_zero_is_never_replaced_by_positive_default():
    lambda_ = frequency_variable()

    with pytest.raises(FourierEvaluationError):
        bplus_symbol(sp.Integer(0), lambda_)


@pytest.mark.parametrize("decay", [-1, -0.5, sp.Integer(-2)])
def test_negative_decay_is_rejected(decay):
    with pytest.raises(FourierEvaluationError):
        bplus_symbol(decay)


@pytest.mark.parametrize("decay", [1 + 2j, sp.I])
def test_nonreal_decay_is_rejected(decay):
    with pytest.raises(FourierEvaluationError):
        bminus_symbol(decay)


@pytest.mark.parametrize(
    "decay",
    [
        float("nan"),
        float("inf"),
        float("-inf"),
        sp.nan,
        sp.oo,
        -sp.oo,
    ],
)
def test_nonfinite_decay_is_rejected(decay):
    with pytest.raises(FourierEvaluationError):
        positive_inverse_integral(decay=decay)


@pytest.mark.parametrize(
    "decay",
    [
        sp.Symbol("d"),
        sp.Symbol("d", negative=True),
        sp.Symbol("d", zero=True),
        sp.Symbol("d", real=False),
    ],
)
def test_incompatible_or_unknown_decay_assumptions_are_rejected(decay):
    with pytest.raises(FourierEvaluationError):
        negative_inverse_integral(decay=decay)


@pytest.mark.parametrize("decay", [True, "1", [1], object()])
def test_incompatible_decay_types_are_rejected(decay):
    with pytest.raises(TypeError):
        bplus_symbol(decay)


def test_positive_explicit_decay_keeps_normalized_kernel_formulas():
    d = sp.Symbol("delta", positive=True)
    t = time_variable()

    assert algebraically_equivalent(
        kplus_kernel(t, d),
        1 / (2 * sp.pi * (d - sp.I * t)),
    )
    assert algebraically_equivalent(
        kminus_kernel(t, d),
        1 / (2 * sp.pi * (d + sp.I * t)),
    )
