import pytest

from symbolic_operator_calculus import (
    Ghat1,
    Ghat2,
    P_minus,
    P_plus,
    R1,
    T1_minus,
    T2_minus,
    U1,
    U1_inverse,
    U2,
    U2_inverse,
    Wminus_21,
    Wplus_12,
    Z1_inverse,
    Z2_inverse,
    LinearCombination,
    build_complete_correction_expansion_trace,
    render_complete_correction_expansion_latex,
)


def factors(expression):
    return tuple(term.product.factors for term in expression.terms)


def coefficients(expression):
    return tuple(term.coefficient for term in expression.terms)


def test_controlled_expansion_counts_are_exactly_1_4_4_16():
    trace = build_complete_correction_expansion_trace()

    assert trace.counts == (1, 4, 4, 16)
    assert len(trace.c0.expanded.terms) == 4
    assert tuple(group.identifier for group in trace.c2) == (
        "C2-A01",
        "C2-A02",
        "C2-A03",
        "C2-A04",
    )


def test_c1_expands_only_differences_and_keeps_auxiliaries_atomic():
    trace = build_complete_correction_expansion_trace()

    assert coefficients(trace.c1) == (1, -1, -1, 1)
    assert factors(trace.c1) == (
        (U2, Wminus_21, T1_minus, R1, U1, Wplus_12, T2_minus),
        (U2, Wminus_21, T1_minus, R1, Ghat1, Wplus_12, T2_minus),
        (Ghat2, Wminus_21, T1_minus, R1, U1, Wplus_12, T2_minus),
        (Ghat2, Wminus_21, T1_minus, R1, Ghat1, Wplus_12, T2_minus),
    )


def test_c2_expands_only_auxiliaries_and_retains_grouped_differences():
    trace = build_complete_correction_expansion_trace()
    left_pieces = tuple(group.left_auxiliary_piece.factors for group in trace.c2)
    right_pieces = tuple(group.right_auxiliary_piece.factors for group in trace.c2)

    assert left_pieces == (
        (U1_inverse, P_plus),
        (U1_inverse, P_plus),
        (P_minus,),
        (P_minus,),
    )
    assert right_pieces == (
        (U2_inverse, P_plus),
        (P_minus,),
        (U2_inverse, P_plus),
        (P_minus,),
    )
    assert all(len(group.factorization.expanded.terms) == 4 for group in trace.c2)
    assert all(
        factors(group.factorization.left_difference) == ((U2,), (Ghat2,))
        for group in trace.c2
    )
    assert all(
        factors(group.factorization.right_difference) == ((U1,), (Ghat1,))
        for group in trace.c2
    )


def test_c3_is_the_deterministic_concatenation_of_four_group_expansions():
    first = build_complete_correction_expansion_trace()
    second = build_complete_correction_expansion_trace()
    reconstructed = LinearCombination(
        tuple(
            term
            for group in first.c2
            for term in group.factorization.expanded.terms
        )
    )

    assert first.c3 == reconstructed
    assert first.c3 == second.c3
    assert coefficients(first.c3) == (1, -1, -1, 1) * 4
    assert len(set(factors(first.c3))) == 16


def test_c3_preserves_every_noncommutative_factor_position():
    trace = build_complete_correction_expansion_trace()

    assert trace.c3.terms[0].product.factors == (
        U2,
        Wminus_21,
        U1_inverse,
        P_plus,
        R1,
        U1,
        Wplus_12,
        U2_inverse,
        P_plus,
    )
    assert trace.c3.terms[-1].product.factors == (
        Ghat2,
        Wminus_21,
        P_minus,
        R1,
        Ghat1,
        Wplus_12,
        P_minus,
    )
    assert all(term.product.factors.count(R1) == 1 for term in trace.c3.terms)
    assert all(
        Z1_inverse not in term.product.factors
        and Z2_inverse not in term.product.factors
        and T1_minus not in term.product.factors
        and T2_minus not in term.product.factors
        for term in trace.c3.terms
    )


def test_all_four_expansion_levels_have_stable_latex_round_trips():
    first = build_complete_correction_expansion_trace()
    second = build_complete_correction_expansion_trace()

    for level in ("C0", "C1", "C2", "C3"):
        rendered = render_complete_correction_expansion_latex(first, level)
        assert rendered == render_complete_correction_expansion_latex(second, level)
        assert r"R_{1}" in rendered
    c3 = render_complete_correction_expansion_latex(first, "C3")
    assert c3.count("C2-T") == 16
    assert c3.count(r"R_{1}") == 16
    with pytest.raises(ValueError):
        render_complete_correction_expansion_latex(first, "all")
