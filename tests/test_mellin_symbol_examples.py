import pytest
import sympy as sp

from symbolic_operator_calculus import (
    ConditionalIdentity,
    ConditionalVerificationStatus,
    ExactIdentity,
    MellinSymbolDependency,
    build_diagonal_mellin_symbols,
    build_dilation_multiplier,
    build_exponential_sinh_symbol,
    build_relative_frequency_symbol,
    build_space_frequency_symbol,
    render_conditional_scalar_identity_latex,
    render_conditional_scalar_identity_text,
    render_mellin_symbol_latex,
    render_mellin_symbol_text,
)


def test_diagonal_symbol_family_has_expected_formulas_and_conditions() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    family = build_diagonal_mellin_symbols(lam, kappa)
    z = sp.pi * (lam + sp.I * kappa)

    assert family.z.as_expr() == z
    assert family.s.as_expr() == sp.coth(z)
    assert family.p_minus.as_expr() == (1 - sp.coth(z)) / 2
    assert family.p_plus.as_expr() == (1 + sp.coth(z)) / 2
    for symbol in (family.z, family.s, family.p_minus, family.p_plus):
        assert symbol.dependency is MellinSymbolDependency.FREQUENCY_ONLY
        assert symbol.assumptions.contains_exact(kappa > 0)
        assert symbol.assumptions.contains_exact(kappa < 1)
    assert family.s.singularities.singularities
    assert family.s.singularities.avoidances


def test_diagonal_relations_are_checked_conditional_scalar_identities() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    family = build_diagonal_mellin_symbols(lam, kappa)
    z = sp.pi * (lam + sp.I * kappa)

    assert isinstance(family.product_identity, ConditionalIdentity)
    assert not isinstance(family.product_identity, ExactIdentity)
    assert sp.simplify(
        family.product_identity.lhs + 1 / (4 * sp.sinh(z) ** 2)
    ) == 0
    assert sp.simplify(
        family.derivative_identity.lhs + sp.pi / sp.sinh(z) ** 2
    ) == 0
    for identity in (family.product_identity, family.derivative_identity):
        assert identity.verification_status is (
            ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
        )
        assert identity.assumption_context.contains_exact(kappa > 0)
        assert identity.singular_set.singularities


def test_exponential_over_sinh_symbol_keeps_bounds_and_poles() -> None:
    lam, beta, kappa = sp.symbols("lambda beta kappa")
    symbol = build_exponential_sinh_symbol(lam, beta, kappa)
    shifted = lam + sp.I * kappa

    assert symbol.as_expr() == sp.exp(beta * shifted) / sp.sinh(sp.pi * shifted)
    for assumption in (beta > -sp.pi, beta < sp.pi, kappa > 0, kappa < 1):
        assert symbol.assumptions.contains_exact(assumption)
    assert symbol.singularities.singularities
    for attribute in ("bounded", "integrable", "variation"):
        with pytest.raises(AttributeError):
            getattr(symbol, attribute)


def test_dilation_multiplier_uses_explicit_positive_real_branch() -> None:
    lam, gamma = sp.symbols("lambda gamma")
    symbol = build_dilation_multiplier(lam, gamma)

    assert symbol.as_expr() == sp.Pow(gamma, sp.I * lam, evaluate=False)
    assert symbol.assumptions.contains_exact(gamma > 0)
    assert symbol.singularities.avoidances
    for attribute in ("finite_variation", "V", "operator"):
        with pytest.raises(AttributeError):
            getattr(symbol, attribute)


def test_space_frequency_example_is_scalar_not_a_pdo() -> None:
    lam, x, beta, kappa = sp.symbols("lambda x beta kappa")
    omega = sp.Function("omega")(x)
    symbol = build_space_frequency_symbol(lam, x, omega, beta, kappa)

    assert symbol.dependency is MellinSymbolDependency.SPACE_FREQUENCY
    assert {item.symbol for item in symbol.spatial_variables} == {x}
    assert sp.exp(sp.I * omega * lam) in sp.Mul.make_args(symbol.as_expr())
    with pytest.raises(AttributeError):
        getattr(symbol, "pdo")


def test_space_frequency_example_rejects_undeclared_omega_variables() -> None:
    lam, x, y, beta, kappa = sp.symbols("lambda x y beta kappa")

    with pytest.raises(ValueError, match="only"):
        build_space_frequency_symbol(lam, x, x + y, beta, kappa)


def test_relative_example_uses_relative_and_frequency_roles() -> None:
    lam, r = sp.symbols("lambda r")
    q = sp.Function("q")(r, lam)
    symbol = build_relative_frequency_symbol(lam, r, q)

    assert symbol.as_expr() == q
    assert symbol.dependency is MellinSymbolDependency.RELATIVE_FREQUENCY
    assert symbol.relative_variable is not None
    assert symbol.relative_variable.symbol == r
    with pytest.raises(AttributeError):
        getattr(symbol, "transform")


def test_relative_example_rejects_missing_or_unknown_dependencies() -> None:
    lam, r, unknown = sp.symbols("lambda r unknown")

    with pytest.raises(ValueError, match="depend"):
        build_relative_frequency_symbol(lam, r, lam)
    with pytest.raises(ValueError, match="undeclared"):
        build_relative_frequency_symbol(lam, r, lam + r + unknown)


def test_rendering_exposes_semantics_without_operator_claims() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    family = build_diagonal_mellin_symbols(lam, kappa)

    text = render_mellin_symbol_text(family.s)
    latex = render_mellin_symbol_latex(family.s)
    identity_text = render_conditional_scalar_identity_text(
        family.derivative_identity
    )
    identity_latex = render_conditional_scalar_identity_latex(
        family.derivative_identity
    )

    assert "Mellin scalar symbol" in text
    assert "Formal Mellin expression" in text
    assert "frequency" in text
    assert "singular" in text.lower()
    assert "Mellin scalar symbol" in latex
    assert "Conditional scalar identity" in identity_text
    assert "Symbolically checked under assumptions" in identity_text
    assert "Conditional scalar identity" in identity_latex
    combined = "\n".join((text, latex, identity_text, identity_latex))
    for forbidden in (
        "Mellin operator",
        "Exact operator identity",
        "Fredholm symbol",
    ):
        assert forbidden not in combined
