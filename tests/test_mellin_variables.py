import pytest
import sympy as sp

from symbolic_operator_calculus import (
    AssumptionContext,
    MellinVariable,
    MellinVariableRole,
    input_spatial_variable,
    mellin_frequency,
    mellin_parameter,
    output_spatial_variable,
    relative_multiplicative_variable,
)


def test_every_mellin_role_is_explicit() -> None:
    x, y, r, lam, kappa = sp.symbols("x y r lambda kappa")

    declarations = (
        output_spatial_variable(x),
        input_spatial_variable(y),
        relative_multiplicative_variable(r),
        mellin_frequency(lam),
        mellin_parameter(kappa),
    )

    assert tuple(item.role for item in declarations) == tuple(MellinVariableRole)


def test_names_do_not_infer_roles() -> None:
    lambda_named = sp.Symbol("lambda")
    x_named = sp.Symbol("x")

    assert output_spatial_variable(lambda_named).role is MellinVariableRole.OUTPUT_SPATIAL
    assert mellin_frequency(x_named).role is MellinVariableRole.FREQUENCY


def test_same_symbol_with_different_roles_is_distinct_and_hashable() -> None:
    symbol = sp.Symbol("ambiguous")
    frequency = mellin_frequency(symbol)
    spatial = output_spatial_variable(symbol)

    assert frequency != spatial
    assert len({frequency, spatial}) == 2


def test_homonymous_symbols_with_different_sympy_assumptions_are_distinct() -> None:
    real_symbol = sp.Symbol("x", real=True)
    positive_symbol = sp.Symbol("x", positive=True)

    assert mellin_parameter(real_symbol) != mellin_parameter(positive_symbol)


@pytest.mark.parametrize("invalid", ["lambda", sp.Symbol("x") + 1, 2])
def test_variable_declarations_reject_non_symbols(invalid: object) -> None:
    with pytest.raises(TypeError, match="SymPy Symbol"):
        MellinVariable(invalid, MellinVariableRole.FREQUENCY)  # type: ignore[arg-type]

    with pytest.raises(TypeError, match="SymPy Symbol"):
        mellin_frequency(invalid)  # type: ignore[arg-type]


def test_external_context_does_not_mutate_sympy_assumptions() -> None:
    lam = sp.Symbol("lambda")
    declaration = mellin_frequency(
        lam,
        assumption_context=AssumptionContext((sp.Contains(lam, sp.S.Reals),)),
    )

    assert lam.is_real is None
    assert declaration.assumption_context.contains_exact(
        sp.Contains(lam, sp.S.Reals)
    )


def test_relative_variable_has_no_implicit_positivity() -> None:
    r = sp.Symbol("r")
    relative = relative_multiplicative_variable(r)

    assert r.is_positive is None
    assert relative.assumption_context.is_empty


def test_with_assumptions_returns_a_new_declaration() -> None:
    kappa = sp.Symbol("kappa")
    original = mellin_parameter(kappa)
    narrowed = original.with_assumptions(kappa > 0)

    assert original.assumption_context.is_empty
    assert narrowed.assumption_context.contains_exact(kappa > 0)
    assert narrowed.role is MellinVariableRole.PARAMETER
