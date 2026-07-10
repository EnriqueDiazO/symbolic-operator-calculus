"""Normalized m=2 off-diagonal Wiener-Hopf block models."""

from __future__ import annotations

from fractions import Fraction

from .operators import (
    G1,
    G2,
    I,
    LinearCombination,
    R11,
    S_Rplus,
    Term,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    main_expression,
)
from .relations import (
    ExactBlock,
    FirstSchurReduction,
    FormalRegularizer,
    ModCompactRelation,
    ModCompactSchurRelation,
    WienerHopfModel,
)


def pplus_operator() -> LinearCombination:
    """Return the exact structural projection ``(I + S_Rplus) / 2``."""

    half = Fraction(1, 2)
    return LinearCombination((Term(half, I), Term(half, S_Rplus)))


def pminus_operator() -> LinearCombination:
    """Return the exact structural projection ``(I - S_Rplus) / 2``."""

    half = Fraction(1, 2)
    return LinearCombination((Term(half, I), Term(-half, S_Rplus)))


def a11_exact_operator() -> LinearCombination:
    """Return the exact expanded operator ``Vtilde_alpha1 Pplus + G1 Pminus``."""

    pplus = pplus_operator()
    pminus = pminus_operator()
    return LinearCombination(
        tuple(Term(term.coefficient, Vtilde_alpha1 * term.product) for term in pplus.terms)
        + tuple(Term(term.coefficient, G1 * term.product) for term in pminus.terms)
    )


def a22_exact_operator() -> LinearCombination:
    """Return the exact expanded operator ``Vtilde_alpha2 Pplus + G2 Pminus``."""

    pplus = pplus_operator()
    pminus = pminus_operator()
    return LinearCombination(
        tuple(Term(term.coefficient, Vtilde_alpha2 * term.product) for term in pplus.terms)
        + tuple(Term(term.coefficient, G2 * term.product) for term in pminus.terms)
    )


def a11_formal_regularizer() -> FormalRegularizer:
    """Return the formal regularizer metadata associated with exact ``A11``."""

    return FormalRegularizer(
        target=ExactBlock("A", 1, 1),
        operator=R11,
    )


def a22_first_schur_reduction() -> FirstSchurReduction:
    """Return the exact metadata ``A22 - A21 R11 A12``."""

    return FirstSchurReduction(
        diagonal=ExactBlock("A", 2, 2),
        left=ExactBlock("A", 2, 1),
        regularizer=a11_formal_regularizer(),
        right=ExactBlock("A", 1, 2),
    )


def a22_first_schur_correction() -> LinearCombination:
    """Return the positive WH correction with its sign derived structurally."""

    schur_sign = a22_first_schur_reduction().offdiagonal_product_coefficient
    leading_a21_sign = a21_wh_model().expression.terms[0].coefficient
    correction_sign = schur_sign * leading_a21_sign
    return correction_sign * main_expression()


def a22_first_schur_model() -> LinearCombination:
    """Return ``A22`` followed by its positive four-term Schur correction."""

    diagonal = a22_exact_operator()
    correction = a22_first_schur_correction()
    return LinearCombination(diagonal.terms + correction.terms)


def a22_first_schur_mod_compact_relation() -> ModCompactSchurRelation:
    """Relate exact ``A22^(1)`` to its operator model modulo compacta."""

    return ModCompactSchurRelation(
        exact=a22_first_schur_reduction(),
        model=a22_first_schur_model(),
    )


def a12_wh_model() -> WienerHopfModel:
    """Return the normalized Wiener-Hopf model for ``A_{1,2}``.

    The original off-diagonal symbol has
    ``b_{1,2} = 2 chi_+ exp(-d lambda)`` and appears in
    ``1/2 (Vtilde_alpha1 - G1) W(b_{1,2})``. The factor ``2`` is already
    absorbed by the exterior ``1/2``, so the normalized model is
    ``(Vtilde_alpha1 - G1) Wplus_12``. Do not add another ``1/2`` or ``2``.
    """

    return WienerHopfModel(
        family="A",
        row=1,
        column=2,
        expression=LinearCombination(
            (
                Term(1, Vtilde_alpha1 * Wplus_12),
                Term(-1, G1 * Wplus_12),
            )
        ),
        normalized=True,
    )


def a21_wh_model() -> WienerHopfModel:
    """Return the normalized Wiener-Hopf model for ``A_{2,1}``.

    The original off-diagonal symbol has
    ``b_{2,1} = -2 chi_- exp(d lambda)`` and appears in
    ``1/2 (Vtilde_alpha2 - G2) W(b_{2,1})``. The factor ``-2`` is already
    absorbed as the exterior minus, so the normalized model is
    ``-(Vtilde_alpha2 - G2) Wminus_21``. Do not add another ``1/2``, ``-2``,
    or extra negative sign.
    """

    return WienerHopfModel(
        family="A",
        row=2,
        column=1,
        expression=LinearCombination(
            (
                Term(-1, Vtilde_alpha2 * Wminus_21),
                Term(1, G2 * Wminus_21),
            )
        ),
        normalized=True,
    )


def a12_mod_compact_relation() -> ModCompactRelation:
    """Return ``A_{1,2} ~= A_{1,2}^{WH}`` as an explicit relation object."""

    return ModCompactRelation(
        exact=ExactBlock("A", 1, 2),
        model=a12_wh_model(),
    )


def a21_mod_compact_relation() -> ModCompactRelation:
    """Return ``A_{2,1} ~= A_{2,1}^{WH}`` as an explicit relation object."""

    return ModCompactRelation(
        exact=ExactBlock("A", 2, 1),
        model=a21_wh_model(),
    )
