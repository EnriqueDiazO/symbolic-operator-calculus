import ast
from fractions import Fraction
from pathlib import Path

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    G1,
    G2,
    I,
    S_Rplus,
    Vtilde_alpha1,
    Vtilde_alpha2,
    LinearCombination,
    ModCompactRelation,
    OperatorAtom,
    WienerHopfModel,
    a11_exact_operator,
    a22_exact_operator,
    apply_atom,
    mvp_atomic_rules,
    pminus_operator,
    pplus_operator,
)


HALF = Fraction(1, 2)


def coefficients(expression):
    return tuple(term.coefficient for term in expression.terms)


def factors(expression):
    return tuple(term.product.factors for term in expression.terms)


def test_identity_and_s_rplus_are_structural_operator_atoms():
    assert isinstance(I, OperatorAtom)
    assert isinstance(S_Rplus, OperatorAtom)
    assert I != S_Rplus


def test_operators_module_does_not_import_sympy():
    operators_path = (
        Path(__file__).parents[1] / "src/symbolic_operator_calculus/operators.py"
    )
    module = ast.parse(operators_path.read_text(encoding="utf-8"))
    imported_roots = {
        alias.name.split(".")[0]
        for node in ast.walk(module)
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    }

    assert "sympy" not in imported_roots


def test_identity_action_returns_the_exact_input_expression():
    x = sp.Symbol("x")
    operand = sp.Integral(sp.Function("f")(x), (x, 0, 1)) + x**2

    result = apply_atom(I, operand, x, mvp_atomic_rules())

    assert result is operand


def test_s_rplus_has_an_atomic_action_after_k3b():
    assert S_Rplus in mvp_atomic_rules()


def test_pplus_is_the_exact_ordered_two_term_projection():
    projection = pplus_operator()

    assert isinstance(projection, LinearCombination)
    assert len(projection.terms) == 2
    assert coefficients(projection) == (HALF, HALF)
    assert factors(projection) == ((I,), (S_Rplus,))
    assert all(isinstance(coefficient, Fraction) for coefficient in coefficients(projection))


def test_pminus_is_the_exact_ordered_two_term_projection():
    projection = pminus_operator()

    assert isinstance(projection, LinearCombination)
    assert len(projection.terms) == 2
    assert coefficients(projection) == (HALF, -HALF)
    assert factors(projection) == ((I,), (S_Rplus,))
    assert all(isinstance(coefficient, Fraction) for coefficient in coefficients(projection))


@pytest.mark.parametrize(
    ("constructor", "shift", "multiplier"),
    (
        (a11_exact_operator, Vtilde_alpha1, G1),
        (a22_exact_operator, Vtilde_alpha2, G2),
    ),
)
def test_diagonal_exact_operator_has_expected_ordered_expansion(
    constructor,
    shift,
    multiplier,
):
    block = constructor()

    assert isinstance(block, LinearCombination)
    assert len(block.terms) == 4
    assert coefficients(block) == (HALF, HALF, HALF, -HALF)
    assert factors(block) == (
        (shift, I),
        (shift, S_Rplus),
        (multiplier, I),
        (multiplier, S_Rplus),
    )
    assert all(isinstance(coefficient, Fraction) for coefficient in coefficients(block))
    assert not isinstance(block, WienerHopfModel)
    assert not isinstance(block, ModCompactRelation)


def test_diagonal_blocks_contain_no_float_or_factor_two_coefficients():
    for block in (a11_exact_operator(), a22_exact_operator()):
        assert all(not isinstance(coefficient, float) for coefficient in coefficients(block))
        assert all(coefficient not in (2, -2) for coefficient in coefficients(block))
