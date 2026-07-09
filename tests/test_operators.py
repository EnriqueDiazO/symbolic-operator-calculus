import ast
from pathlib import Path

import pytest

from symbolic_operator_calculus import LinearCombination, OperatorAtom, Product, Term


def atoms():
    return (
        OperatorAtom("A"),
        OperatorAtom("B"),
        OperatorAtom("C"),
        OperatorAtom("D"),
    )


def term_data(combination):
    return tuple(
        (term.coefficient, term.product.factors)
        for term in combination.terms
    )


def test_operator_product_two_atoms_remains_product():
    A, B, _, _ = atoms()

    result = A * B

    assert isinstance(result, Product)
    assert result.factors == (A, B)


def test_operator_product_three_atoms_preserves_order():
    A, B, C, _ = atoms()

    result = A * B * C

    assert result.factors == (A, B, C)


def test_product_rejects_scalar_left_factor():
    A, _, _, _ = atoms()

    with pytest.raises(TypeError):
        Product((2, A))


def test_product_rejects_scalar_right_factor():
    A, _, _, _ = atoms()

    with pytest.raises(TypeError):
        Product((A, 2))


def test_operator_addition_creates_ordered_linear_combination():
    A, B, _, _ = atoms()

    result = A + B

    assert isinstance(result, LinearCombination)
    assert term_data(result) == ((1, (A,)), (1, (B,)))


def test_operator_subtraction_creates_negative_coefficient():
    A, B, _, _ = atoms()

    result = A - B

    assert term_data(result) == ((1, (A,)), (-1, (B,)))


def test_operator_negation_creates_negative_single_term():
    A, _, _, _ = atoms()

    result = -A

    assert term_data(result) == ((-1, (A,)),)


def test_product_negation_preserves_factor_order():
    A, B, _, _ = atoms()

    result = -(A * B)

    assert term_data(result) == ((-1, (A, B)),)


def test_left_scalar_multiplication_of_atom_creates_coefficient():
    A, _, _, _ = atoms()

    result = 2 * A

    assert isinstance(result, LinearCombination)
    assert term_data(result) == ((2, (A,)),)


def test_right_scalar_multiplication_of_atom_creates_coefficient():
    A, _, _, _ = atoms()

    result = A * 2

    assert isinstance(result, LinearCombination)
    assert term_data(result) == ((2, (A,)),)


def test_negative_scalar_multiplication_of_atom():
    A, _, _, _ = atoms()

    result = -2 * A

    assert term_data(result) == ((-2, (A,)),)


def test_left_scalar_multiplication_of_product():
    A, B, _, _ = atoms()

    result = 3 * (A * B)

    assert term_data(result) == ((3, (A, B)),)


def test_right_scalar_multiplication_of_product():
    A, B, _, _ = atoms()

    result = (A * B) * 3

    assert term_data(result) == ((3, (A, B)),)


def test_atom_plus_scaled_atom():
    A, B, _, _ = atoms()

    result = A + 2 * B

    assert term_data(result) == ((1, (A,)), (2, (B,)))


def test_compound_linear_combination_with_scaled_terms():
    A, B, C, _ = atoms()

    result = A + 2 * B - 3 * C

    assert term_data(result) == ((1, (A,)), (2, (B,)), (-3, (C,)))


def test_scalar_times_linear_combination_scales_all_coefficients():
    A, B, _, _ = atoms()

    result = 2 * (A + B)

    assert term_data(result) == ((2, (A,)), (2, (B,)))


def test_linear_combination_times_scalar_scales_all_coefficients():
    A, B, _, _ = atoms()

    result = (A - B) * 3

    assert term_data(result) == ((3, (A,)), (-3, (B,)))


def test_sum_of_two_linear_combinations_is_flat_and_ordered():
    A, B, C, D = atoms()

    result = (A + B) + (C - D)

    assert term_data(result) == (
        (1, (A,)),
        (1, (B,)),
        (1, (C,)),
        (-1, (D,)),
    )


def test_empty_linear_combination_represents_zero_operator():
    result = LinearCombination(())

    assert result.terms == ()


def test_empty_product_represents_empty_composition():
    result = Product(())

    assert result.factors == ()


def test_zero_times_atom_is_empty_linear_combination():
    A, _, _, _ = atoms()

    result = 0 * A

    assert result == LinearCombination(())


def test_atom_plus_zero_returns_single_term_combination():
    A, _, _, _ = atoms()

    result = A + 0

    assert term_data(result) == ((1, (A,)),)


def test_zero_plus_atom_returns_single_term_combination():
    A, _, _, _ = atoms()

    result = 0 + A

    assert term_data(result) == ((1, (A,)),)


def test_product_rejects_linear_combination_factor():
    A, B, _, _ = atoms()

    with pytest.raises(TypeError):
        Product((A + B, A))


def test_linear_combination_times_atom_is_not_supported():
    A, B, C, _ = atoms()

    with pytest.raises(TypeError):
        (A + B) * C


def test_atom_times_linear_combination_is_not_supported():
    A, B, C, _ = atoms()

    with pytest.raises(TypeError):
        A * (B + C)


def test_two_scaled_atoms_work_naturally():
    A, B, _, _ = atoms()

    result = 2 * A + 3 * B

    assert term_data(result) == ((2, (A,)), (3, (B,)))


def test_term_zero_coefficient_is_removed_by_linear_combination():
    A, _, _, _ = atoms()

    result = LinearCombination((Term(0, A),))

    assert result == LinearCombination(())


def test_term_hashing_and_equality_remain_structural():
    A, B, _, _ = atoms()

    first = 2 * A + 3 * B
    second = LinearCombination((Term(2, Product((A,))), Term(3, Product((B,)))))

    assert first == second
    assert hash(first) == hash(second)


def test_operator_product_equality_remains_noncommutative():
    A, B, _, _ = atoms()

    assert A * B != B * A


def test_float_and_complex_coefficients_are_supported():
    A, B, _, _ = atoms()

    result = 1.5 * A + (2 + 3j) * B

    assert term_data(result) == ((1.5, (A,)), (2 + 3j, (B,)))


def test_bool_is_not_a_supported_coefficient():
    A, _, _, _ = atoms()

    with pytest.raises(TypeError):
        True * A


def test_product_rejects_nested_product_factor():
    A, B, _, _ = atoms()

    with pytest.raises(TypeError):
        Product((A * B,))


def test_operators_module_does_not_import_sympy():
    source = Path("src/symbolic_operator_calculus/operators.py").read_text()
    tree = ast.parse(source)
    imports = [
        alias.name
        for node in tree.body
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    ]

    assert "sympy" not in imports
