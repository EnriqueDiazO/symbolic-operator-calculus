import ast
from dataclasses import FrozenInstanceError
from fractions import Fraction
from pathlib import Path

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    G1,
    G2,
    R11,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    ExactBlock,
    FirstSchurReduction,
    ModCompactSchurRelation,
    OperatorAtom,
    Product,
    a11_exact_operator,
    a11_formal_regularizer,
    a12_wh_model,
    a21_wh_model,
    a22_exact_operator,
    a22_first_schur_correction,
    a22_first_schur_model,
    a22_first_schur_mod_compact_relation,
    a22_first_schur_reduction,
    apply_linear_combination_ordered,
    main_expression,
    mvp_atomic_rules,
)


HALF = Fraction(1, 2)
CORRECTION_FACTORS = (
    (Vtilde_alpha2, Wminus_21, R11, Vtilde_alpha1, Wplus_12),
    (Vtilde_alpha2, Wminus_21, R11, G1, Wplus_12),
    (G2, Wminus_21, R11, Vtilde_alpha1, Wplus_12),
    (G2, Wminus_21, R11, G1, Wplus_12),
)


def coefficients(expression):
    return tuple(term.coefficient for term in expression.terms)


def factors(expression):
    return tuple(term.product.factors for term in expression.terms)


def test_first_schur_reduction_is_immutable_non_ast_metadata():
    reduction = a22_first_schur_reduction()

    assert isinstance(reduction, FirstSchurReduction)
    assert not isinstance(reduction, OperatorAtom)
    assert not isinstance(reduction, Product)
    with pytest.raises(FrozenInstanceError):
        reduction.diagonal = ExactBlock("A", 1, 1)


def test_schur_metadata_module_does_not_import_sympy():
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


def test_exact_first_reduction_preserves_all_blocks_regularizer_and_minus_sign():
    reduction = a22_first_schur_reduction()

    assert reduction.diagonal == ExactBlock("A", 2, 2)
    assert reduction.left == ExactBlock("A", 2, 1)
    assert reduction.regularizer == a11_formal_regularizer()
    assert reduction.right == ExactBlock("A", 1, 2)
    assert reduction.offdiagonal_product_coefficient == -1
    assert not isinstance(reduction.left, OperatorAtom)
    assert not isinstance(reduction.right, OperatorAtom)


def test_offdiagonal_models_keep_normalization_and_derive_positive_correction():
    outer_schur_sign = a22_first_schur_reduction().offdiagonal_product_coefficient
    leading_a21_sign = a21_wh_model().expression.terms[0].coefficient

    assert coefficients(a12_wh_model().expression) == (1, -1)
    assert coefficients(a21_wh_model().expression) == (-1, 1)
    assert outer_schur_sign == -1
    assert leading_a21_sign == -1
    assert outer_schur_sign * leading_a21_sign == 1


def test_correction_reuses_existing_four_term_expression():
    correction = a22_first_schur_correction()

    assert correction == main_expression()
    assert coefficients(correction) == (1, -1, -1, 1)
    assert factors(correction) == CORRECTION_FACTORS


def test_complete_model_places_exact_a22_before_correction():
    diagonal = a22_exact_operator()
    correction = a22_first_schur_correction()
    model = a22_first_schur_model()

    assert model.terms[:4] == diagonal.terms
    assert model.terms[4:] == correction.terms
    assert coefficients(diagonal) == (HALF, HALF, HALF, -HALF)
    assert coefficients(model) == (HALF, HALF, HALF, -HALF, 1, -1, -1, 1)
    assert all(
        isinstance(coefficient, Fraction)
        for coefficient in coefficients(diagonal)
    )


def test_mod_compact_schur_relation_keeps_exact_and_model_distinct():
    relation = a22_first_schur_mod_compact_relation()

    assert isinstance(relation, ModCompactSchurRelation)
    assert relation.exact == a22_first_schur_reduction()
    assert relation.model == a22_first_schur_model()
    assert relation.exact != relation.model
    assert not hasattr(relation, "rewrite")
    assert not hasattr(relation, "as_expr")


def test_schur_construction_adds_no_regularizer_cancellation_or_orientation():
    reduction = a22_first_schur_reduction()
    right_products = tuple(
        term.product * R11 for term in a11_exact_operator().terms
    )
    left_products = tuple(
        R11 * term.product for term in a11_exact_operator().terms
    )

    assert not hasattr(reduction, "left_relation")
    assert not hasattr(reduction, "right_relation")
    assert not hasattr(reduction.regularizer, "left_relation")
    assert not hasattr(reduction.regularizer, "right_relation")
    assert all(product.factors[-1] is R11 for product in right_products)
    assert all(product.factors[0] is R11 for product in left_products)


def test_complete_model_application_reuses_existing_engine_and_term_actions():
    x = sp.Symbol("x")
    f = sp.Function("f")
    rules = mvp_atomic_rules()
    model = a22_first_schur_model()

    complete = apply_linear_combination_ordered(model, f(x), x, rules)
    diagonal = apply_linear_combination_ordered(a22_exact_operator(), f(x), x, rules)
    correction = apply_linear_combination_ordered(main_expression(), f(x), x, rules)

    assert len(complete.terms) == 8
    assert complete.terms[:4] == diagonal.terms
    assert complete.terms[4:] == correction.terms
    assert tuple(term.coefficient for term in complete.terms) == (
        HALF,
        HALF,
        HALF,
        -HALF,
        1,
        -1,
        -1,
        1,
    )
    assert all(isinstance(term.expression, sp.Expr) for term in complete.terms)
