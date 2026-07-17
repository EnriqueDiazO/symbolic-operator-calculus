import ast
from fractions import Fraction
from pathlib import Path

import pytest

from symbolic_operator_calculus import (
    G1,
    G2,
    I,
    R11,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    LinearCombination,
    OperatorAtom,
    Product,
    Term,
)


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


def test_operator_atom_valid_fields_remain_structural_and_hashable():
    default = OperatorAtom("A")
    explicit = OperatorAtom("A", "operator")
    wiener_hopf = OperatorAtom("Wplus_12", "wiener_hopf")
    spaced = OperatorAtom(" A ", " operator ")

    assert default == explicit
    assert default != wiener_hopf
    assert hash(default) == hash(explicit)
    assert {default, explicit, wiener_hopf} == {default, wiener_hopf}
    assert {default: "atom"}[explicit] == "atom"
    assert spaced.name == " A "
    assert spaced.kind == " operator "


@pytest.mark.parametrize("value", [None, 1, [], {}, (), b"A"])
def test_operator_atom_rejects_non_string_name(value):
    with pytest.raises(TypeError, match=r"OperatorAtom\.name"):
        OperatorAtom(value)


@pytest.mark.parametrize("value", ["", " ", "   ", "\t", "\n"])
def test_operator_atom_rejects_empty_or_whitespace_name(value):
    with pytest.raises(ValueError, match=r"OperatorAtom\.name"):
        OperatorAtom(value)


@pytest.mark.parametrize("value", [None, 1, [], {}, (), b"A"])
def test_operator_atom_rejects_non_string_kind(value):
    with pytest.raises(TypeError, match=r"OperatorAtom\.kind"):
        OperatorAtom("A", value)


@pytest.mark.parametrize("value", ["", " ", "   ", "\t", "\n"])
def test_operator_atom_rejects_empty_or_whitespace_kind(value):
    with pytest.raises(ValueError, match=r"OperatorAtom\.kind"):
        OperatorAtom("A", value)


def test_public_operator_constants_remain_hashable_and_composable():
    constants = (
        I,
        G1,
        G2,
        Vtilde_alpha1,
        Vtilde_alpha2,
        Wplus_12,
        Wminus_21,
        R11,
    )

    assert len(set(constants)) == len(constants)
    assert all({constant: constant.name}[constant] == constant.name for constant in constants)
    assert Product(constants).factors == constants
    assert tuple(Term(1, constant).product.factors for constant in constants) == tuple(
        (constant,) for constant in constants
    )


def test_product_normalizes_supported_iterables_and_remains_hashable():
    A, B, _, _ = atoms()

    from_tuple = Product((A, B))
    from_list = Product([A, B])
    from_generator = Product(atom for atom in (A, B))

    assert from_tuple == from_list == from_generator
    assert isinstance(from_list.factors, tuple)
    assert from_generator.factors == (A, B)
    assert hash(from_tuple) == hash(from_list) == hash(from_generator)
    assert {from_tuple, from_list, from_generator} == {from_tuple}
    assert {from_tuple: "product"}[from_list] == "product"


@pytest.mark.parametrize("factor", [object(), "A", None])
def test_product_rejects_non_atom_factors_with_local_error(factor):
    with pytest.raises(TypeError, match=r"Product\.factors"):
        Product((factor,))


def test_product_rejects_non_iterable_container_with_local_error():
    with pytest.raises(TypeError, match=r"Product\.factors"):
        Product(None)


@pytest.mark.parametrize("coefficient", [1, -1, 0, 0.5, 1 + 2j, Fraction(1, 2)])
def test_term_supported_coefficients_remain_hashable(coefficient):
    A, _, _, _ = atoms()

    term = Term(coefficient, Product((A,)))

    assert term.coefficient == coefficient
    assert hash(term) == hash(Term(coefficient, A))


def test_term_normalizes_operator_atom_shorthand():
    A, _, _, _ = atoms()

    from_atom = Term(1, A)
    from_product = Term(1, Product((A,)))

    assert from_atom == from_product
    assert isinstance(from_atom.product, Product)
    assert from_atom.product.factors == (A,)
    assert hash(from_atom) == hash(from_product)
    assert {from_atom: "term"}[from_product] == "term"


@pytest.mark.parametrize("product", [None, "A", object(), [OperatorAtom("A")]])
def test_term_rejects_invalid_product_with_local_error(product):
    with pytest.raises(TypeError, match=r"Term\.product"):
        Term(1, product)


def test_linear_combination_normalizes_supported_iterables_and_is_hashable():
    A, B, _, _ = atoms()
    term_a = Term(1, A)
    term_b = Term(-1, B)

    from_tuple = LinearCombination((term_a, term_b))
    from_list = LinearCombination([term_a, term_b])
    from_generator = LinearCombination(term for term in (term_a, term_b))

    assert from_tuple == from_list == from_generator
    assert isinstance(from_list.terms, tuple)
    assert from_generator.terms == (term_a, term_b)
    assert hash(from_tuple) == hash(from_list) == hash(from_generator)
    assert {from_tuple, from_list, from_generator} == {from_tuple}
    assert {from_tuple: "combination"}[from_list] == "combination"


@pytest.mark.parametrize("value", [object(), None, "term"])
def test_linear_combination_rejects_non_term_elements(value):
    with pytest.raises(TypeError, match=r"LinearCombination\.terms"):
        LinearCombination((value,))


def test_linear_combination_rejects_duck_typed_term():
    class Fake:
        coefficient = 1
        product = "bad"

    with pytest.raises(TypeError, match=r"LinearCombination\.terms"):
        LinearCombination((Fake(),))


def test_linear_combination_rejects_non_iterable_container():
    with pytest.raises(TypeError, match=r"LinearCombination\.terms"):
        LinearCombination(None)


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
    assert {result: "empty combination"}[LinearCombination(())] == "empty combination"


def test_empty_product_represents_empty_composition():
    result = Product(())

    assert result.factors == ()
    assert {result: "empty product"}[Product(())] == "empty product"


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
