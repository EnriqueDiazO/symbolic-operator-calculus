from dataclasses import FrozenInstanceError

import pytest

from symbolic_operator_calculus import (
    A11_inverse,
    A12,
    A21,
    A22,
    A22_first,
    G1,
    G2,
    N2,
    N2_first,
    P_minus,
    P_plus,
    R1,
    T1_minus,
    T2_minus,
    U1_inverse,
    U2_inverse,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    Z1_inverse,
    Z2_inverse,
    CertificationStatus,
    ExactIdentity,
    FormalIdentity,
    LinearCombination,
    ModCompactEquivalence,
    NormalizedFirstSchurPivotDerivation,
    OperatorAtom,
    Product,
    a12_mod_compact_relation,
    a21_mod_compact_relation,
    a22_first_schur_exact_expression,
    build_normalized_first_schur_pivot_derivation,
    compose_ordered,
    factor_normalized_mod_compact_correction,
    normalized_a22_first_schur_expression,
    normalized_a22_with_regularizer_expression,
    normalized_first_schur_exact_expression,
    normalized_first_schur_mod_compact_expression,
    render_normalized_first_schur_pivot_exact_latex,
    render_normalized_first_schur_pivot_latex,
    render_normalized_first_schur_pivot_mod_compact_latex,
    substitute_operator_subproduct,
    t1_minus_expression,
    t2_minus_expression,
)


EXACT_CORRECTION_FACTORS = (
    Z2_inverse,
    A21,
    T1_minus,
    R1,
    Z1_inverse,
    A12,
    T2_minus,
)

MOD_COMPACT_CORRECTION_FACTORS = (
    (
        Z2_inverse,
        Vtilde_alpha2,
        Wminus_21,
        T1_minus,
        R1,
        Z1_inverse,
        Vtilde_alpha1,
        Wplus_12,
        T2_minus,
    ),
    (
        Z2_inverse,
        Vtilde_alpha2,
        Wminus_21,
        T1_minus,
        R1,
        Z1_inverse,
        G1,
        Wplus_12,
        T2_minus,
    ),
    (
        Z2_inverse,
        G2,
        Wminus_21,
        T1_minus,
        R1,
        Z1_inverse,
        Vtilde_alpha1,
        Wplus_12,
        T2_minus,
    ),
    (
        Z2_inverse,
        G2,
        Wminus_21,
        T1_minus,
        R1,
        Z1_inverse,
        G1,
        Wplus_12,
        T2_minus,
    ),
)


def factors(expression):
    return tuple(term.product.factors for term in expression.terms)


def coefficients(expression):
    return tuple(term.coefficient for term in expression.terms)


def test_all_required_normalized_pivot_atoms_are_noncommutative():
    atoms = (
        Z1_inverse,
        Z2_inverse,
        A12,
        A21,
        A22,
        A11_inverse,
        A22_first,
        T1_minus,
        T2_minus,
        U1_inverse,
        U2_inverse,
        P_plus,
        P_minus,
        R1,
        N2,
        N2_first,
        Vtilde_alpha1,
        Vtilde_alpha2,
        G1,
        G2,
        Wplus_12,
        Wminus_21,
    )

    assert all(isinstance(atom, OperatorAtom) for atom in atoms)
    assert all(not atom.is_commutative for atom in atoms)
    assert R1.is_formal_regularizer


def test_auxiliary_shift_formulas_are_stored_without_expanding_the_pivots():
    assert coefficients(t1_minus_expression()) == (1, 1)
    assert factors(t1_minus_expression()) == (
        (U1_inverse, P_plus),
        (P_minus,),
    )
    assert coefficients(t2_minus_expression()) == (1, 1)
    assert factors(t2_minus_expression()) == (
        (U2_inverse, P_plus),
        (P_minus,),
    )


def test_exact_reduction_and_normalization_preserve_all_factor_order():
    reduced = a22_first_schur_exact_expression()
    normalized = normalized_a22_first_schur_expression()

    assert coefficients(reduced) == (1, -1)
    assert factors(reduced) == (
        (A22,),
        (A21, A11_inverse, A12),
    )
    assert coefficients(normalized) == (1, -1)
    assert factors(normalized) == (
        (Z2_inverse, A22, T2_minus),
        (Z2_inverse, A21, A11_inverse, A12, T2_minus),
    )


def test_formal_inverse_substitution_changes_only_the_requested_atom():
    substituted = normalized_a22_with_regularizer_expression()

    assert coefficients(substituted) == (1, -1)
    assert factors(substituted) == (
        (Z2_inverse, A22, T2_minus),
        EXACT_CORRECTION_FACTORS,
    )
    assert all(
        A11_inverse not in term.product.factors for term in substituted.terms
    )
    assert all(
        term.product.factors.count(R1) <= 1 for term in substituted.terms
    )


def test_diagonal_recognition_produces_the_required_exact_identity():
    exact = normalized_first_schur_exact_expression()

    assert coefficients(exact) == (1, -1)
    assert factors(exact) == ((N2,), EXACT_CORRECTION_FACTORS)


def test_a21_then_a12_mod_compact_substitutions_give_expected_signs_and_order():
    expanded = normalized_first_schur_mod_compact_expression()

    assert coefficients(expanded) == (1, 1, -1, -1, 1)
    assert factors(expanded) == ((N2,), *MOD_COMPACT_CORRECTION_FACTORS)
    assert all(
        A21 not in term.product.factors and A12 not in term.product.factors
        for term in expanded.terms
    )


def test_factorized_mod_compact_form_reconstructs_controlled_expansion():
    factorization = factor_normalized_mod_compact_correction()
    reconstructed = compose_ordered(
        factorization.prefix,
        factorization.left_difference,
        factorization.left_wiener_hopf,
        factorization.bridge,
        factorization.right_difference,
        factorization.right_wiener_hopf,
        factorization.suffix,
    )

    assert reconstructed == factorization.expanded_correction
    assert factors(factorization.left_difference) == (
        (Vtilde_alpha2,),
        (G2,),
    )
    assert coefficients(factorization.left_difference) == (1, -1)
    assert factors(factorization.right_difference) == (
        (Vtilde_alpha1,),
        (G1,),
    )
    assert coefficients(factorization.right_difference) == (1, -1)


def test_exact_formal_and_mod_compact_relations_remain_distinct():
    trace = build_normalized_first_schur_pivot_derivation()

    assert isinstance(trace.reduced_definition, ExactIdentity)
    assert isinstance(trace.inverse_substitution, FormalIdentity)
    assert isinstance(trace.exact_result, ExactIdentity)
    assert isinstance(trace.mod_compact_result, ModCompactEquivalence)
    assert trace.mod_compact_result.certification_status is (
        CertificationStatus.UNCERTIFIED
    )
    assert trace.exact_result != trace.mod_compact_result
    assert trace.offdiagonal_relations == (
        a21_mod_compact_relation(),
        a12_mod_compact_relation(),
    )


def test_r1_remains_one_atomic_factor_in_every_expanded_correction_term():
    exact = normalized_first_schur_exact_expression()
    mod_compact = normalized_first_schur_mod_compact_expression()

    assert exact.terms[1].product.factors.count(R1) == 1
    assert all(
        term.product.factors.count(R1) == 1
        for term in mod_compact.terms[1:]
    )
    assert not hasattr(R1, "symbol")
    assert not hasattr(R1, "inverse")


def test_normalized_pivot_latex_matches_required_factorized_forms():
    trace = build_normalized_first_schur_pivot_derivation()

    assert render_normalized_first_schur_pivot_exact_latex(trace) == (
        r"N_{2}^{(1)} = N_{2} - Z_{2}^{-1}\,A_{2,1}\,T_{1,-}\,R_{1}\,"
        r"Z_{1}^{-1}\,A_{1,2}\,T_{2,-}"
    )
    assert render_normalized_first_schur_pivot_mod_compact_latex(trace) == (
        r"N_{2}^{(1)} \simeq N_{2} + Z_{2}^{-1}\,"
        r"\left(\widetilde V_{\alpha_2} - G_{2}\right)\,W^-_{2,1}\,"
        r"T_{1,-}\,R_{1}\,Z_{1}^{-1}\,"
        r"\left(\widetilde V_{\alpha_1} - G_{1}\right)\,W^+_{1,2}\,T_{2,-}"
    )


def test_controlled_expansion_contains_four_correction_terms_and_atomic_r1():
    trace = build_normalized_first_schur_pivot_derivation()
    latex = render_normalized_first_schur_pivot_mod_compact_latex(
        trace,
        expand_differences=True,
    )

    assert latex.count(r"Z_{2}^{-1}") == 4
    assert latex.count(r"R_{1}") == 4
    assert r"R_{1}^{-1}" not in latex
    assert r"\left(\widetilde V" not in latex


def test_derivation_and_rendering_are_deterministic_and_immutable():
    first = build_normalized_first_schur_pivot_derivation()
    second = build_normalized_first_schur_pivot_derivation()

    assert isinstance(first, NormalizedFirstSchurPivotDerivation)
    assert first == second
    assert render_normalized_first_schur_pivot_latex(first) == (
        render_normalized_first_schur_pivot_latex(second)
    )
    with pytest.raises(FrozenInstanceError):
        first.exact_result = second.exact_result


def test_ordered_substitution_is_explicit_distributive_and_noncommutative():
    left = OperatorAtom("left")
    target = OperatorAtom("target")
    right = OperatorAtom("right")
    replacement_left = OperatorAtom("replacement_left")
    replacement_right = OperatorAtom("replacement_right")
    expression = left * target * right
    replacement = replacement_left - replacement_right

    substituted = substitute_operator_subproduct(
        expression,
        target,
        replacement,
    )

    assert coefficients(substituted) == (1, -1)
    assert factors(substituted) == (
        (left, replacement_left, right),
        (left, replacement_right, right),
    )
    assert substituted != substitute_operator_subproduct(
        right * target * left,
        target,
        replacement,
    )


def test_unrequested_products_with_linear_combinations_remain_unsupported():
    difference = Vtilde_alpha1 - G1

    assert isinstance(difference, LinearCombination)
    with pytest.raises(TypeError):
        Z1_inverse * difference
    with pytest.raises(TypeError):
        difference * Z1_inverse


def test_empty_ordered_composition_is_the_structural_identity():
    identity = compose_ordered()

    assert identity == LinearCombination.from_factor(Product(()))
