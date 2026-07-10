from dataclasses import FrozenInstanceError

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    DilationConjugatedWienerHopfFactorization,
    FourierScalingError,
    RelativeWienerHopfError,
    article_wiener_hopf_kernel,
    factor_dilation_conjugated_wiener_hopf,
    forward_exponential,
    halfline_indicator_is_scale_invariant,
    inverse_exponential,
    inverse_prefactor,
    kminus_kernel,
    kplus_kernel,
    m2_relative_wiener_hopf_factorizations,
    scaled_convolution_kernel,
    scaled_fourier_symbol,
)


@pytest.fixture(scope="module")
def variables():
    gamma1, gamma2, d = sp.symbols(
        "gamma1 gamma2 d",
        real=True,
        positive=True,
    )
    x, y, t, lambda_ = sp.symbols("x y t lambda", real=True)
    return gamma1, gamma2, d, x, y, t, lambda_


@pytest.fixture(scope="module")
def m2(variables):
    gamma1, gamma2, d, *_ = variables
    return m2_relative_wiener_hopf_factorizations(gamma1, gamma2, d)


def test_project_fourier_convention_and_original_kernel_symbol_pairs(variables):
    _, _, d, _, _, t, lambda_ = variables
    h12 = article_wiener_hopf_kernel(t, -d)
    h21 = article_wiener_hopf_kernel(t, d)

    assert forward_exponential(t, lambda_) == sp.exp(-sp.I * t * lambda_)
    assert inverse_exponential(t, lambda_) == sp.exp(sp.I * t * lambda_)
    assert inverse_prefactor() == 1 / (2 * sp.pi)
    assert sp.simplify(h12 - 2 * kplus_kernel(t, d)) == 0
    assert sp.simplify(h21 - (-2) * kminus_kernel(t, d)) == 0


def test_conjugated_action_change_of_variable_produces_exact_gamma_j(variables):
    gamma1, gamma2, d, x, y, t, _ = variables
    h = article_wiener_hopf_kernel(t, -d)
    reconstructed = gamma2 * h.xreplace({t: gamma1 * x - gamma2 * y})
    expected = -gamma2 / (
        sp.pi * sp.I * (gamma1 * x - gamma2 * y + sp.I * d)
    )

    assert sp.simplify(reconstructed - expected) == 0
    assert sp.simplify(reconstructed / h.xreplace({t: gamma1 * x - gamma2 * y})) == gamma2


def test_fourier_scaling_law_is_structural_and_requires_positive_scale(variables):
    gamma1, _, _, _, _, t, lambda_ = variables
    h = sp.Function("h")
    b = sp.Function("b")

    assert scaled_convolution_kernel(h(t), t, gamma1) == gamma1 * h(gamma1 * t)
    assert scaled_fourier_symbol(b(lambda_), lambda_, gamma1) == b(lambda_ / gamma1)
    with pytest.raises(FourierScalingError):
        scaled_convolution_kernel(h(t), t, sp.Symbol("a", real=True))


def test_relative_scale_and_both_exact_kernel_reconstructions(m2):
    w12, w21 = m2

    assert isinstance(w12, DilationConjugatedWienerHopfFactorization)
    assert w12.relative_scale == w12.gamma_k / w12.gamma_j
    assert w21.relative_scale == w21.gamma_k / w21.gamma_j
    for factorization in m2:
        assert factorization.kernels_reconstruct_exactly
        assert sp.simplify(
            factorization.left_reconstructed_kernel
            - factorization.right_reconstructed_kernel
        ) == 0
        assert sp.simplify(
            factorization.left_reconstructed_kernel
            - factorization.integral_kernel
        ) == 0


def test_left_and_right_scaled_kernels_use_the_correct_dilation(m2):
    for factorization in m2:
        t = factorization.time_variable
        h = factorization.original_kernel
        assert factorization.left_convolution_kernel == (
            factorization.gamma_j * h.xreplace({t: factorization.gamma_j * t})
        )
        assert factorization.right_convolution_kernel == (
            factorization.gamma_k * h.xreplace({t: factorization.gamma_k * t})
        )


def test_positive_halfline_indicators_are_stable_only_with_positive_scale(variables):
    gamma1, _, _, _, _, _, lambda_ = variables

    assert halfline_indicator_is_scale_invariant("plus", lambda_, gamma1)
    assert halfline_indicator_is_scale_invariant("minus", lambda_, gamma1)
    with pytest.raises(RelativeWienerHopfError):
        halfline_indicator_is_scale_invariant(
            "plus",
            lambda_,
            sp.Symbol("gamma", real=True),
        )


def test_general_k_less_j_symbols_keep_original_factor_two(variables):
    gamma1, gamma2, d, *_ = variables
    result = factor_dilation_conjugated_wiener_hopf(
        1,
        3,
        gamma1,
        gamma2,
        sp.Integer(0),
        d,
    )
    lambda_ = result.frequency_variable
    chi_plus = sp.Function("chi_plus")

    assert result.original_symbol == 2 * chi_plus(lambda_) * sp.exp(-d * lambda_)
    assert result.left_symbol == 2 * chi_plus(lambda_) * sp.exp(
        -d * lambda_ / gamma2
    )
    assert result.right_symbol == 2 * chi_plus(lambda_) * sp.exp(
        -d * lambda_ / gamma1
    )


def test_general_k_greater_j_symbols_keep_original_factor_minus_two(variables):
    gamma1, gamma2, d, *_ = variables
    result = factor_dilation_conjugated_wiener_hopf(
        3,
        1,
        gamma2,
        gamma1,
        d,
        sp.Integer(0),
    )
    lambda_ = result.frequency_variable
    chi_minus = sp.Function("chi_minus")

    assert result.original_symbol == -2 * chi_minus(lambda_) * sp.exp(d * lambda_)
    assert result.left_symbol == -2 * chi_minus(lambda_) * sp.exp(
        d * lambda_ / gamma1
    )
    assert result.right_symbol == -2 * chi_minus(lambda_) * sp.exp(
        d * lambda_ / gamma2
    )


def test_all_four_m2_symbols_use_the_shift_side_denominator(m2):
    w12, w21 = m2
    lambda_ = w12.frequency_variable
    chi_plus = sp.Function("chi_plus")
    chi_minus = sp.Function("chi_minus")
    gamma1, gamma2 = w12.gamma_k, w12.gamma_j
    d = -w12.delta_c

    assert w12.left_symbol == 2 * chi_plus(lambda_) * sp.exp(-d * lambda_ / gamma2)
    assert w12.right_symbol == 2 * chi_plus(lambda_) * sp.exp(-d * lambda_ / gamma1)
    assert w21.left_symbol == -2 * chi_minus(lambda_) * sp.exp(d * lambda_ / gamma1)
    assert w21.right_symbol == -2 * chi_minus(lambda_) * sp.exp(d * lambda_ / gamma2)
    for symbol in (
        w12.left_symbol,
        w12.right_symbol,
        w21.left_symbol,
        w21.right_symbol,
    ):
        assert not symbol.has(sp.pi)
        assert not symbol.has(sp.Rational(1, 2))


def test_gamma_and_delta_assumptions_are_enforced(variables):
    gamma1, _, d, *_ = variables
    undecided = sp.Symbol("gamma", real=True)

    with pytest.raises(RelativeWienerHopfError):
        factor_dilation_conjugated_wiener_hopf(
            1,
            2,
            gamma1,
            undecided,
            0,
            d,
        )
    with pytest.raises(RelativeWienerHopfError):
        factor_dilation_conjugated_wiener_hopf(
            1,
            2,
            gamma1,
            gamma1,
            d,
            0,
        )


def test_factorization_metadata_is_immutable(m2):
    with pytest.raises(FrozenInstanceError):
        m2[0].relative_scale = sp.Integer(1)


def test_nontrivial_numeric_kernel_sanity():
    result = factor_dilation_conjugated_wiener_hopf(
        1,
        2,
        sp.Rational(3, 2),
        sp.Rational(5, 3),
        sp.Integer(1),
        sp.Integer(4),
    )
    substitutions = {
        result.output_variable: sp.Rational(7, 5),
        result.input_variable: sp.Rational(4, 5),
    }
    values = [
        complex(sp.N(expression.subs(substitutions), 16))
        for expression in (
            result.integral_kernel,
            result.left_reconstructed_kernel,
            result.right_reconstructed_kernel,
        )
    ]

    assert abs(values[0] - values[1]) < 1e-13
    assert abs(values[0] - values[2]) < 1e-13
