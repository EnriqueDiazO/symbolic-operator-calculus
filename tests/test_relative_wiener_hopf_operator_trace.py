from dataclasses import FrozenInstanceError

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    DilationOperatorModel,
    OrderedRelativeOperatorProduct,
    RelativeWienerHopfDerivationTrace,
    WienerHopfOperatorModel,
    build_relative_wiener_hopf_trace,
    render_relative_wiener_hopf_trace_latex,
)


@pytest.fixture(scope="module")
def trace():
    gamma1, gamma2, d = sp.symbols("gamma1 gamma2 d", positive=True, real=True)
    return build_relative_wiener_hopf_trace(1, 2, gamma1, gamma2, 0, d)


def test_identity_is_immutable_typed_and_exact(trace):
    assert isinstance(trace, RelativeWienerHopfDerivationTrace)
    assert trace.identity.exact
    assert isinstance(trace.original_operator, OrderedRelativeOperatorProduct)
    assert not isinstance(trace.original_operator, sp.Mul)
    with pytest.raises(FrozenInstanceError):
        trace.identity.exact = False


def test_products_preserve_noncommutative_factor_order(trace):
    original = trace.original_operator.factors
    left = trace.left_operator.factors
    right = trace.right_operator.factors
    assert len(original) == 3
    assert isinstance(original[0], DilationOperatorModel)
    assert isinstance(original[1], WienerHopfOperatorModel)
    assert isinstance(original[2], DilationOperatorModel)
    assert original[0].dilation.k == 1
    assert original[1].placement == "original"
    assert original[2].inverse and original[2].dilation.k == 2
    assert [type(factor) for factor in left] == [
        DilationOperatorModel,
        WienerHopfOperatorModel,
    ]
    assert [type(factor) for factor in right] == [
        WienerHopfOperatorModel,
        DilationOperatorModel,
    ]


def test_trace_reuses_l1_objects_and_reconstructions(trace):
    l1 = trace.factorization
    assert trace.relative_dilation.scale is l1.relative_scale
    assert trace.original_operator.factors[1].symbol is l1.original_symbol
    assert trace.left_operator.factors[1].symbol is l1.left_symbol
    assert trace.right_operator.factors[0].symbol is l1.right_symbol
    assert trace.action.kernel is l1.integral_kernel
    assert trace.reconstruction_checks == (True, True)
    assert sp.simplify(l1.left_reconstructed_kernel - trace.action.kernel) == 0
    assert sp.simplify(l1.right_reconstructed_kernel - trace.action.kernel) == 0


def test_symbol_correspondence_and_common_action(trace):
    l1 = trace.factorization
    assert sp.simplify(trace.symbol_correspondence_forward - l1.right_symbol) == 0
    assert sp.simplify(trace.symbol_correspondence_reverse - l1.left_symbol) == 0
    assert trace.action.integral.function == (
        trace.action.kernel * trace.action.input_function(trace.action.input_variable)
    )


def test_m2_reverse_case_keeps_minus_two_and_denominators():
    gamma1, gamma2, d = sp.symbols("gamma1 gamma2 d", positive=True, real=True)
    trace = build_relative_wiener_hopf_trace(2, 1, gamma2, gamma1, d, 0)
    frequency = trace.factorization.frequency_variable
    chi_minus = sp.Function("chi_minus")
    assert trace.factorization.left_symbol == (
        -2 * chi_minus(frequency) * sp.exp(d * frequency / gamma1)
    )
    assert trace.factorization.right_symbol == (
        -2 * chi_minus(frequency) * sp.exp(d * frequency / gamma2)
    )


def test_renderer_accepts_only_trace_and_shows_mathematical_notation(trace):
    rendered = render_relative_wiener_hopf_trace_latex(trace)
    latex = rendered.as_latex()
    assert len(rendered.steps) == 10
    assert r"V_{\alpha_1}" in latex
    assert r"V_{\alpha_2}^{-1}" in latex
    assert r"V_{\beta_{1,2}}" in latex
    assert r"b_{1,2}^{\mathrm L}" in latex
    assert r"b_{1,2}^{\mathrm R}" in latex
    assert r"\mathcal K_{1,2}(x,y)" in latex
    assert r"\int_0^\infty" in latex
    assert r"\simeq" not in latex
    assert r"\widetilde V" not in latex
    assert r"\rho" not in latex
    for internal_name in ("relative_scale", "left_symbol", "right_symbol", "V_beta"):
        assert internal_name not in latex
    with pytest.raises(TypeError):
        render_relative_wiener_hopf_trace_latex(trace.factorization)
