import pytest
import sympy as sp

from symbolic_operator_calculus import (
    SafeSubstitutionError,
    collect_bound_symbols,
    fresh_symbol,
    substitute_free_variable,
)


def integral_fixture():
    x, y = sp.symbols("x y")
    f = sp.Function("f")
    K = sp.Function("K")
    return x, y, f, K, sp.Integral(K(x, y) * f(y), (y, 0, sp.oo))


def test_free_symbols_distinguishes_free_and_bound_variables():
    x, y, _, _, expression = integral_fixture()

    assert expression.free_symbols == {x}
    assert y not in expression.free_symbols


def test_integral_bound_symbols_identifies_bound_variable():
    _, y, _, _, expression = integral_fixture()

    assert tuple(expression.bound_symbols) == (y,)
    assert tuple(expression.variables) == (y,)
    assert collect_bound_symbols(expression) == {y}


def test_substitute_free_variable_changes_free_x_but_not_bound_y():
    x, y, f, K, expression = integral_fixture()
    gamma = sp.Symbol("gamma")

    result = substitute_free_variable(expression, x, gamma * x)

    expected = sp.Integral(K(gamma * x, y) * f(y), (y, 0, sp.oo))
    assert result == expected
    assert expression == sp.Integral(K(x, y) * f(y), (y, 0, sp.oo))


def test_substitute_free_variable_alpha_renames_on_capture_risk():
    x, y, f, K, expression = integral_fixture()
    u = sp.Symbol("u")

    result = substitute_free_variable(expression, x, y)

    expected = sp.Integral(K(y, u) * f(u), (u, 0, sp.oo))
    assert result == expected
    assert result.free_symbols == {y}
    assert collect_bound_symbols(result) == {u}


def test_substitute_free_variable_keeps_nested_integral_variables_distinct():
    x, y, f, K, _ = integral_fixture()
    u, v = sp.symbols("u v")
    expression = sp.Integral(
        K(x, u) * sp.Integral(K(u, y) * f(y), (y, 0, sp.oo)),
        (u, 0, sp.oo),
    )

    result = substitute_free_variable(expression, x, u)

    expected = sp.Integral(
        K(u, v) * sp.Integral(K(v, y) * f(y), (y, 0, sp.oo)),
        (v, 0, sp.oo),
    )
    assert result == expected
    assert collect_bound_symbols(result) == {y, v}
    assert result.free_symbols == {u}


def test_fresh_symbol_is_deterministic_and_uses_indexed_fallback():
    y, u, v, w, z = sp.symbols("y u v w z")

    assert fresh_symbol({u}) == y
    assert fresh_symbol({y, u, v, w, z}) == sp.Symbol("y1")


def test_substitute_free_variable_rejects_non_sympy_expression():
    x = sp.Symbol("x")

    with pytest.raises(SafeSubstitutionError):
        substitute_free_variable("x", x, x)
