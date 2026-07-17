from fractions import Fraction
from pathlib import Path
import ast

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    G1,
    LinearCombination,
    OperatorAtom,
    Product,
    Term,
    apply_linear_combination,
    expand_ordered,
    mvp_atomic_rules,
)


def term_data(combination):
    return tuple(
        (term.coefficient, term.product.factors)
        for term in combination.terms
    )


@pytest.mark.parametrize(
    "coefficient",
    (
        0,
        1,
        -1,
        0.0,
        1.5,
        -2.5,
        Fraction(0, 1),
        Fraction(1, 2),
        Fraction(-3, 4),
        0j,
        1 + 0j,
        1 + 2j,
        -1 + 2j,
        2j,
        -2j,
    ),
)
def test_finite_python_coefficients_are_accepted_and_hashable(coefficient):
    term = Term(coefficient, OperatorAtom("A"))

    assert term.coefficient == coefficient
    assert hash(term) == hash(Term(coefficient, OperatorAtom("A")))
    if isinstance(coefficient, complex):
        assert isinstance(term.coefficient, complex)


@pytest.mark.parametrize("coefficient", [True, False])
def test_bool_is_still_rejected_as_coefficient(coefficient):
    with pytest.raises(TypeError):
        Term(coefficient, OperatorAtom("A"))


@pytest.mark.parametrize(
    "coefficient",
    [float("nan"), float("inf"), float("-inf")],
)
def test_nonfinite_float_coefficients_are_rejected(coefficient):
    with pytest.raises(ValueError, match=r"Term\.coefficient.*finite"):
        Term(coefficient, OperatorAtom("A"))


@pytest.mark.parametrize(
    "coefficient",
    [
        complex(float("nan"), 0),
        complex(0, float("nan")),
        complex(float("inf"), 0),
        complex(0, float("inf")),
        complex(float("-inf"), 1),
        complex(1, float("-inf")),
    ],
)
def test_nonfinite_complex_coefficients_are_rejected(coefficient):
    with pytest.raises(ValueError, match=r"Term\.coefficient.*finite"):
        Term(coefficient, OperatorAtom("A"))


@pytest.mark.parametrize(
    "coefficient",
    [
        sp.Integer(2),
        sp.Rational(1, 2),
        sp.Float(1.5),
        sp.I,
        sp.Symbol("a"),
        sp.Symbol("a", real=True),
        sp.Symbol("a", positive=True),
    ],
)
def test_sympy_coefficients_remain_outside_the_ast_contract(coefficient):
    with pytest.raises(TypeError, match=r"Term\.coefficient"):
        Term(coefficient, OperatorAtom("A"))


def test_complex_zero_remains_stored_and_is_filtered_from_combinations():
    product = Product((G1,))
    term = Term(0j, product)

    assert term.coefficient == 0j
    assert isinstance(term.coefficient, complex)
    assert LinearCombination((term,)) == LinearCombination(())


def test_fraction_one_half_is_accepted_as_coefficient():
    A = OperatorAtom("A")

    term = Term(Fraction(1, 2), A)

    assert term.coefficient == Fraction(1, 2)


def test_fraction_one_half_remains_fraction_not_float():
    A = OperatorAtom("A")

    result = Fraction(1, 2) * A

    assert result.terms[0].coefficient == Fraction(1, 2)
    assert isinstance(result.terms[0].coefficient, Fraction)
    assert not isinstance(result.terms[0].coefficient, float)


def test_negative_fraction_remains_exact():
    A = OperatorAtom("A")

    result = Fraction(-1, 2) * A

    assert result.terms[0].coefficient == Fraction(-1, 2)
    assert isinstance(result.terms[0].coefficient, Fraction)


def test_expanding_half_difference_preserves_fraction_coefficients_in_order():
    A = OperatorAtom("A")
    B = OperatorAtom("B")

    result = expand_ordered(Fraction(1, 2) * (A - B))

    assert term_data(result) == (
        (Fraction(1, 2), (A,)),
        (Fraction(-1, 2), (B,)),
    )


def test_operators_module_still_does_not_import_sympy():
    source = Path("src/symbolic_operator_calculus/operators.py").read_text()
    tree = ast.parse(source)
    imports = [
        alias.name
        for node in tree.body
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    ]

    assert "sympy" not in imports


def test_fraction_coefficient_projects_to_sympy_rational_in_actions_layer():
    x = sp.Symbol("x")
    f = sp.Function("f")
    combination = LinearCombination((Term(Fraction(1, 2), Product((G1,))),))

    result = apply_linear_combination(combination, f(x), x, mvp_atomic_rules())

    assert result == sp.Rational(1, 2) * sp.Function("G1")(x) * f(x)
