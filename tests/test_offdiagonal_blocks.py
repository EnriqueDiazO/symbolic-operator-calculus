from fractions import Fraction

import sympy as sp

from symbolic_operator_calculus import (
    G1,
    G2,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    ModCompactRelation,
    a12_mod_compact_relation,
    a12_wh_model,
    a21_mod_compact_relation,
    a21_wh_model,
    apply_linear_combination_ordered,
    mvp_atomic_rules,
)


def coefficients(model):
    return tuple(term.coefficient for term in model.expression.terms)


def products(model):
    return tuple(term.product.factors for term in model.expression.terms)


def applied_expression(model):
    x = sp.Symbol("x")
    f = sp.Function("f")
    return apply_linear_combination_ordered(
        model.expression,
        f(x),
        x,
        mvp_atomic_rules(),
    )


def test_a12_wh_model_has_indices_one_two():
    model = a12_wh_model()

    assert (model.row, model.column) == (1, 2)


def test_a12_wh_model_has_exactly_two_terms():
    model = a12_wh_model()

    assert len(model.expression.terms) == 2


def test_a12_wh_model_coefficients_are_plus_one_minus_one():
    assert coefficients(a12_wh_model()) == (1, -1)


def test_a12_wh_model_products_are_normalized_ordered_terms():
    assert products(a12_wh_model()) == (
        (Vtilde_alpha1, Wplus_12),
        (G1, Wplus_12),
    )


def test_a12_wh_model_contains_no_half_coefficient():
    assert Fraction(1, 2) not in coefficients(a12_wh_model())


def test_a12_wh_model_contains_no_factor_two_coefficient():
    assert 2 not in coefficients(a12_wh_model())


def test_a12_wh_model_application_produces_two_formal_integrals():
    applied = applied_expression(a12_wh_model())

    assert len(applied.as_expr().atoms(sp.Integral)) == 2
    assert all(term.expression.has(sp.Function("Lplus_12")) for term in applied.terms)


def test_a21_wh_model_has_indices_two_one():
    model = a21_wh_model()

    assert (model.row, model.column) == (2, 1)


def test_a21_wh_model_has_exactly_two_terms():
    model = a21_wh_model()

    assert len(model.expression.terms) == 2


def test_a21_wh_model_coefficients_are_minus_one_plus_one():
    assert coefficients(a21_wh_model()) == (-1, 1)


def test_a21_wh_model_products_are_normalized_ordered_terms():
    assert products(a21_wh_model()) == (
        (Vtilde_alpha2, Wminus_21),
        (G2, Wminus_21),
    )


def test_a21_wh_model_contains_no_half_coefficient():
    assert Fraction(1, 2) not in coefficients(a21_wh_model())


def test_a21_wh_model_contains_no_minus_two_coefficient():
    assert -2 not in coefficients(a21_wh_model())


def test_a21_wh_model_application_produces_two_formal_integrals():
    applied = applied_expression(a21_wh_model())

    assert len(applied.as_expr().atoms(sp.Integral)) == 2
    assert all(term.expression.has(sp.Function("Lminus_21")) for term in applied.terms)


def test_a12_normalization_has_no_coefficient_with_absolute_value_two():
    assert all(abs(coefficient) != 2 for coefficient in coefficients(a12_wh_model()))


def test_a21_normalization_has_no_coefficient_with_absolute_value_two():
    assert all(abs(coefficient) != 2 for coefficient in coefficients(a21_wh_model()))


def test_a12_wh_model_keeps_positive_global_sign():
    assert coefficients(a12_wh_model())[0] == 1


def test_a21_wh_model_keeps_exterior_negative_in_coefficients():
    assert coefficients(a21_wh_model()) == (-1, 1)


def test_offdiagonal_relations_are_mod_compact_relations():
    assert isinstance(a12_mod_compact_relation(), ModCompactRelation)
    assert isinstance(a21_mod_compact_relation(), ModCompactRelation)


def test_relation_exact_side_does_not_identify_with_model():
    relation = a21_mod_compact_relation()

    assert relation.exact is not relation.model
    assert relation.exact != relation.model
