import ast
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest
import sympy as sp

import symbolic_operator_calculus as soc
from symbolic_operator_calculus import (
    I,
    R11,
    ExactBlock,
    FormalRegularizer,
    KernelAnnotatedExpression,
    KernelRepresentationRequiredError,
    ModCompactRelation,
    OperatorAtom,
    Product,
    WienerHopfModel,
    a11_exact_operator,
    a11_formal_regularizer,
    apply_atom,
    mvp_atomic_rules,
)
from semantic_helpers import regularizer_rules


def test_formal_regularizer_is_immutable_non_operator_metadata():
    association = a11_formal_regularizer()

    assert isinstance(association, FormalRegularizer)
    assert not isinstance(association, OperatorAtom)
    with pytest.raises(FrozenInstanceError):
        association.operator = I


def test_formal_regularizer_module_does_not_import_sympy():
    relations_path = (
        Path(__file__).parents[1] / "src/symbolic_operator_calculus/relations.py"
    )
    module = ast.parse(relations_path.read_text(encoding="utf-8"))
    imported_roots = {
        alias.name.split(".")[0]
        for node in ast.walk(module)
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    }

    assert "sympy" not in imported_roots


def test_a11_formal_regularizer_has_exact_target_and_existing_operator():
    association = a11_formal_regularizer()

    assert association.target == ExactBlock("A", 1, 1)
    assert association.operator is R11
    assert isinstance(association.operator, OperatorAtom)
    assert association.operator.is_formal_regularizer


def test_a11_formal_regularizer_api_is_deterministic():
    first = a11_formal_regularizer()
    second = a11_formal_regularizer()

    assert first == second == FormalRegularizer(ExactBlock("A", 1, 1), R11)


def test_a11_exact_operator_remains_the_exact_four_term_expression():
    block = a11_exact_operator()

    assert len(block.terms) == 4
    assert tuple(term.product.factors for term in block.terms) == (
        (soc.Vtilde_alpha1, soc.I),
        (soc.Vtilde_alpha1, soc.S_Rplus),
        (soc.G1, soc.I),
        (soc.G1, soc.S_Rplus),
    )


def test_r11_without_kernel_representation_is_rejected():
    x, y = sp.symbols("x y")
    f = sp.Function("f")

    with pytest.raises(KernelRepresentationRequiredError):
        apply_atom(
            R11,
            f(x),
            x,
            mvp_atomic_rules(),
            integration_variable=y,
        )


def test_explicit_r11_kernel_representation_remains_formal_and_annotated():
    x, y = sp.symbols("x y")
    f = sp.Function("f")

    result = apply_atom(
        R11,
        f(x),
        x,
        regularizer_rules(),
        integration_variable=y,
    )

    assert isinstance(result, KernelAnnotatedExpression)
    assert result.expression == sp.Integral(
        sp.Function("R11")(x, y) * f(y),
        (y, 0, sp.oo),
    )
    assert result.has(sp.Function("R11")(x, y))
    assert not result.has(sp.Inverse)


def test_products_around_a11_terms_do_not_cancel_to_identity():
    block = a11_exact_operator()
    right_products = tuple(term.product * R11 for term in block.terms)
    left_products = tuple(R11 * term.product for term in block.terms)

    assert all(isinstance(product, Product) for product in right_products)
    assert all(isinstance(product, Product) for product in left_products)
    assert all(product.factors[-1] is R11 for product in right_products)
    assert all(product.factors[0] is R11 for product in left_products)
    assert all(product != I for product in right_products + left_products)


def test_full_a11_expression_is_not_coerced_into_a_product_with_r11():
    with pytest.raises(TypeError):
        _ = a11_exact_operator() * R11
    with pytest.raises(TypeError):
        _ = R11 * a11_exact_operator()


def test_regularizer_metadata_is_neither_exact_nor_mod_compact_relation():
    association = a11_formal_regularizer()

    assert not isinstance(association, ModCompactRelation)
    assert not isinstance(association, WienerHopfModel)
    assert not hasattr(association, "left_relation")
    assert not hasattr(association, "right_relation")
    assert not hasattr(association, "as_expr")


def test_no_general_inverse_or_a22_regularizer_api_was_added():
    association = a11_formal_regularizer()

    assert not hasattr(association, "inverse")
    assert not hasattr(R11, "inverse")
    assert not hasattr(a11_exact_operator(), "inverse")
    assert not hasattr(soc, "a22_formal_regularizer")
