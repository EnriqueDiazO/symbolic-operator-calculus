from fractions import Fraction
from pathlib import Path
import ast

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


def test_bool_is_still_rejected_as_coefficient():
    A = OperatorAtom("A")

    try:
        Term(True, A)
    except TypeError:
        pass
    else:
        raise AssertionError("bool coefficients must be rejected")


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
