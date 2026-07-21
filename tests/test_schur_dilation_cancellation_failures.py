from __future__ import annotations

from dataclasses import dataclass, replace

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    AlgebraMembershipStatus,
    AssumptionContext,
    CancellationLayerStatus,
    ComplexDomain,
    DerivationRelationKind,
    ExactBlock,
    MellinExpression,
    MellinPseudodifferentialOperator,
    MellinSymbol,
    MellinSymbolDependency,
    MellinSymbolDomain,
    OperatorAtom,
    Product,
    RegularizerMellinRepresentation,
    RegularizerOperator,
    RegularizerRepresentationStatus,
    SchurBlockRepresentation,
    SchurCancellationStatus,
    SingularSet,
    WeightedDilationOperator,
    WeightedLpSpace,
    WienerHopfOperator,
    cancel_relative_dilations_in_schur_correction,
    mellin_frequency,
    output_spatial_variable,
)


@dataclass(frozen=True)
class _FailureCase:
    gamma: sp.Symbol
    assumptions: AssumptionContext
    space: WeightedLpSpace
    left_wh: WienerHopfOperator
    right_wh: WienerHopfOperator
    dilation: WeightedDilationOperator
    inverse_dilation: WeightedDilationOperator
    mellin_operator: MellinPseudodifferentialOperator
    regularizer: RegularizerOperator
    left_atom: OperatorAtom
    right_atom: OperatorAtom


def _case() -> _FailureCase:
    gamma = sp.Symbol("gamma")
    lam = sp.Symbol("lambda", real=True)
    radial = sp.Symbol("x", positive=True)
    assumptions = AssumptionContext((gamma > 0,))
    space = _space(assumptions)

    frequency = mellin_frequency(lam)
    spatial = output_spatial_variable(radial)
    domain = MellinSymbolDomain(
        frequency=frequency,
        complex_domain=ComplexDomain.real_line(lam),
        spatial_variables=(spatial,),
        assumption_context=assumptions,
        singular_set=SingularSet(),
        description="failure-test Mellin domain",
        evidence=("typed failure-test domain",),
    )
    symbol = MellinSymbol(
        MellinExpression(
            radial + lam,
            domain,
            variables=(frequency, spatial),
            evidence=("failure-test Mellin symbol",),
            description="r11",
        ),
        MellinSymbolDependency.SPACE_FREQUENCY,
        description="r11",
    )
    mellin_operator = MellinPseudodifferentialOperator(
        atom=OperatorAtom("Op_r11"),
        symbol=symbol,
        domain=space,
        codomain=space,
        symbol_class="tilde-E",
        source="failure-test Mellin PDO",
        radial_scaling_stable=True,
        stability_evidence=("radial-scaling proof",),
        evidence=("bounded Mellin PDO",),
    )
    left_wh = _wh(
        atom="W_b21",
        symbol=lam + 2,
        lam=lam,
        space=space,
    )
    right_wh = _wh(
        atom="W_b12",
        symbol=lam - 2,
        lam=lam,
        space=space,
    )
    dilation = WeightedDilationOperator(
        atom=OperatorAtom("U_gamma"),
        gamma=gamma,
        domain=space,
        codomain=space,
        assumptions=assumptions,
        source="failure-test weighted dilation",
        evidence=("weighted norm proof",),
    )
    inverse_dilation = dilation.inverse(OperatorAtom("U_gamma_inverse"))
    regularizer = RegularizerOperator(
        target=ExactBlock("A", 1, 1),
        operator=OperatorAtom("R11", kind="formal_regularizer"),
    )
    return _FailureCase(
        gamma,
        assumptions,
        space,
        left_wh,
        right_wh,
        dilation,
        inverse_dilation,
        mellin_operator,
        regularizer,
        OperatorAtom("A21"),
        OperatorAtom("A12"),
    )


def _space(
    assumptions: AssumptionContext,
    *,
    p: sp.Expr = sp.Integer(2),
    weight_exponent: sp.Expr = sp.Integer(0),
) -> WeightedLpSpace:
    return WeightedLpSpace(
        label="X_delta",
        p=p,
        weight_exponent=weight_exponent,
        assumptions=assumptions,
        source="failure-test multiplicative weight norm",
    )


def _wh(
    *,
    atom: str,
    symbol: sp.Expr,
    lam: sp.Symbol,
    space: WeightedLpSpace,
    localized: bool = False,
    localization_controlled: bool = True,
) -> WienerHopfOperator:
    return WienerHopfOperator(
        atom=OperatorAtom(atom),
        symbol=symbol,
        frequency_variable=lam,
        domain=space,
        codomain=space,
        symbol_class="weighted Fourier multipliers",
        source="failure-test Wiener--Hopf factor",
        frequency_scaling_stable=True,
        stability_evidence=("positive frequency scaling",),
        localized=localized,
        localization_controlled=localization_controlled,
        evidence=("bounded Wiener--Hopf factor",),
    )


def _blocks(
    case: _FailureCase,
    *,
    left_factors: tuple[object, ...] | None = None,
    right_factors: tuple[object, ...] | None = None,
    evidence: tuple[object, ...] = ("block representation evidence",),
    localization_controlled: bool = True,
    left_compact_ideal: object = "compact operators",
    right_compact_ideal: object = "compact operators",
) -> tuple[SchurBlockRepresentation, SchurBlockRepresentation]:
    left = SchurBlockRepresentation(
        exact_block=ExactBlock("A", 2, 1),
        atom=case.left_atom,
        factors=(
            (case.left_wh, case.dilation)
            if left_factors is None
            else left_factors
        ),
        relation_kind=DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
        source="explicit failure-test A21 representation",
        hypotheses=("bounded left model",),
        evidence=evidence,
        localization_controlled=localization_controlled,
        compact_ideal=left_compact_ideal,
    )
    right = SchurBlockRepresentation(
        exact_block=ExactBlock("A", 1, 2),
        atom=case.right_atom,
        factors=(
            (case.inverse_dilation, case.right_wh)
            if right_factors is None
            else right_factors
        ),
        relation_kind=DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
        source="explicit failure-test A12 representation",
        hypotheses=("bounded right model",),
        evidence=evidence,
        localization_controlled=localization_controlled,
        compact_ideal=right_compact_ideal,
    )
    return left, right


def _representation(
    case: _FailureCase,
    *,
    regularizer: RegularizerOperator | None = None,
    mellin_operator: MellinPseudodifferentialOperator | None = None,
) -> RegularizerMellinRepresentation:
    return RegularizerMellinRepresentation(
        regularizer=regularizer or case.regularizer,
        mellin_operator=mellin_operator or case.mellin_operator,
        status=RegularizerRepresentationStatus.CERTIFIED_EXACT,
        hypotheses=(),
        source="failure-test regularizer representation",
        evidence=("regularizer representation evidence",),
        space=(mellin_operator or case.mellin_operator).domain,
    )


def _assert_blocked(
    result: object,
    *,
    original: Product,
    reason_fragment: str,
) -> None:
    assert result.status is SchurCancellationStatus.BLOCKED_WITH_PRECISE_OBLIGATION
    assert result.relation is None
    assert not result.succeeded
    assert result.expression == original
    assert result.original_expression == original
    assert result.expression.factors == original.factors
    assert result.semantic_relation is None
    assert result.scaled_regularizer is None
    assert reason_fragment in " ".join(result.blocked_reasons)
    assert result.classification.dilation_cancellation is (
        CancellationLayerStatus.NOT_APPLIED
    )
    assert result.classification.mixed_algebra_membership is (
        CancellationLayerStatus.UNPROVED_MIXED_ALGEBRA_MEMBERSHIP
    )
    assert result.algebra_membership.status is AlgebraMembershipStatus.UNPROVED
    obligation_keys = {item.key for item in result.proof_obligations}
    assert {
        "regularizer_mellin_representation",
        "wh_mellin_wh_sandwich_membership",
        "fredholm_algebra_membership",
        "schur_correction_symbol",
    } <= obligation_keys
    assert any(key.startswith("blocked_") for key in obligation_keys)
    assert result.latex


def test_raw_atoms_are_never_expanded_or_simplified_implicitly():
    case = _case()
    original = Product((case.left_atom, case.regularizer.operator, case.right_atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=case.left_atom,
        regularizer=case.regularizer,
        right_block=case.right_atom,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="raw AST atoms are never expanded implicitly",
    )
    assert result.model_expression is None
    assert result.evidence == ()
    assert result.classification.initial_block_representations is (
        CancellationLayerStatus.NOT_APPLIED
    )
    assert result.classification.regularizer_mellin_representation is (
        CancellationLayerStatus.CONDITIONAL_ON_REGULARIZER_REPRESENTATION
    )


def test_abstract_regularizer_does_not_gain_an_implicit_mellin_representation():
    case = _case()
    left, right = _blocks(case)
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=case.regularizer,
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="regularizer remains abstract",
    )
    assert result.model_expression is not None
    assert case.regularizer.operator in result.model_expression.factors
    assert len(result.evidence) == 2
    assert result.classification.initial_block_representations is (
        CancellationLayerStatus.EQUALITY_MODULO_COMPACTS
    )
    assert result.classification.regularizer_mellin_representation is (
        CancellationLayerStatus.CONDITIONAL_ON_REGULARIZER_REPRESENTATION
    )


def test_noninverse_dilations_leave_the_original_schur_word_unchanged():
    case = _case()
    wrong_inverse = WeightedDilationOperator(
        atom=OperatorAtom("U_wrong_inverse"),
        gamma=sp.Integer(3),
        domain=case.space,
        codomain=case.space,
        assumptions=case.assumptions,
        source="deliberately noninverse test dilation",
        evidence=("bounded positive dilation",),
    )
    left, right = _blocks(
        case,
        right_factors=(wrong_inverse, case.right_wh),
    )
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="dilation scales are not inverse",
    )
    assert len(result.evidence) == 3
    assert result.classification.initial_block_representations is (
        CancellationLayerStatus.EQUALITY_MODULO_COMPACTS
    )
    assert result.classification.regularizer_mellin_representation is (
        CancellationLayerStatus.EXACT_IDENTITY
    )


@pytest.mark.parametrize(
    ("p", "weight_exponent"),
    (
        pytest.param(sp.Integer(3), sp.Integer(0), id="different-p"),
        pytest.param(
            sp.Integer(2),
            sp.Rational(1, 4),
            id="different-weight",
        ),
    ),
)
def test_different_p_or_weight_blocks_cancellation(
    p: sp.Expr,
    weight_exponent: sp.Expr,
):
    case = _case()
    other_space = _space(
        case.assumptions,
        p=p,
        weight_exponent=weight_exponent,
    )
    other_dilation = WeightedDilationOperator(
        atom=OperatorAtom("U_other_space"),
        gamma=case.gamma,
        domain=other_space,
        codomain=other_space,
        assumptions=case.assumptions,
        source="dilation on an incompatible weighted space",
        evidence=("other-space norm calculation",),
    )
    other_inverse = other_dilation.inverse(OperatorAtom("U_other_space_inverse"))
    other_right_wh = replace(
        case.right_wh,
        atom=OperatorAtom("W_b12_other_space"),
        domain=other_space,
        codomain=other_space,
    )
    left, right = _blocks(
        case,
        right_factors=(other_inverse, other_right_wh),
    )
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="identical domain, codomain, p, weight",
    )


def test_different_normalization_blocks_cancellation():
    case = _case()
    wrong_inverse = replace(
        case.inverse_dilation,
        atom=OperatorAtom("U_wrong_normalization"),
        normalization_power=case.inverse_dilation.normalization_power + 1,
    )
    left, right = _blocks(
        case,
        right_factors=(wrong_inverse, case.right_wh),
    )
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="different normalizations",
    )


def test_different_action_convention_blocks_cancellation():
    case = _case()
    wrong_inverse = replace(
        case.inverse_dilation,
        atom=OperatorAtom("U_wrong_convention"),
    )
    # The public constructor accepts only the repository convention. Corrupting a
    # frozen copy emulates incompatible metadata loaded across a serialization
    # boundary and verifies that the cancellation rule still fails closed.
    object.__setattr__(wrong_inverse, "convention", "f(x) -> f(x/gamma)")
    left, right = _blocks(
        case,
        right_factors=(wrong_inverse, case.right_wh),
    )
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="different action conventions",
    )


def test_extra_intermediate_factor_is_not_crossed_or_removed():
    case = _case()
    cutoff = OperatorAtom("chi_infinity")
    left, right = _blocks(
        case,
        left_factors=(case.left_wh, cutoff, case.dilation),
    )
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="left block must have exactly W(b21), V_gamma order",
    )
    assert result.model_expression is not None
    assert cutoff in result.model_expression.factors


def test_nonstable_mellin_symbol_class_is_not_scaled():
    case = _case()
    unstable_pdo = replace(
        case.mellin_operator,
        atom=OperatorAtom("Op_r11_unstable"),
        radial_scaling_stable=False,
        stability_evidence=(),
    )
    left, right = _blocks(case)
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case, mellin_operator=unstable_pdo),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="lacks radial-scaling stability evidence",
    )


@pytest.mark.parametrize("side", ("left", "right"))
def test_incorrect_factor_order_is_not_repaired(side: str):
    case = _case()
    left_factors: tuple[object, ...] | None = None
    right_factors: tuple[object, ...] | None = None
    expected = "left block first factor"
    if side == "left":
        left_factors = (case.dilation, case.left_wh)
    else:
        right_factors = (case.right_wh, case.inverse_dilation)
        expected = "right block first factor"
    left, right = _blocks(
        case,
        left_factors=left_factors,
        right_factors=right_factors,
    )
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(result, original=original, reason_fragment=expected)


def test_block_representations_without_evidence_are_rejected():
    case = _case()
    left, right = _blocks(case, evidence=())
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="require explicit evidence",
    )


def test_uncontrolled_localization_is_retained_and_blocks_cancellation():
    case = _case()
    localized_left = replace(
        case.left_wh,
        atom=OperatorAtom("chi_W_b21_chi"),
        localized=True,
        localization_controlled=False,
    )
    left, right = _blocks(
        case,
        left_factors=(localized_left, case.dilation),
        localization_controlled=False,
    )
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="chi_infinity is not controlled under dilation",
    )
    assert result.model_expression is not None
    assert localized_left.atom in result.model_expression.factors


def test_regularizer_for_the_wrong_target_block_is_not_substituted():
    case = _case()
    wrong_regularizer = RegularizerOperator(
        target=ExactBlock("A", 2, 2),
        operator=case.regularizer.operator,
    )
    left, right = _blocks(case)
    original = Product((left.atom, wrong_regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case, regularizer=wrong_regularizer),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="regularizer must target exact block A11",
    )


def test_incompatible_compact_ideals_block_the_schur_propagation():
    case = _case()
    left, right = _blocks(
        case,
        left_compact_ideal="ideal_L",
        right_compact_ideal="ideal_R",
    )
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="incompatible compact ideals",
    )
    assert result.evidence[:2] == (
        left.semantic_relation,
        right.semantic_relation,
    )


def test_localized_wiener_hopf_is_blocked_even_with_a_control_boolean():
    case = _case()
    localized_left = replace(
        case.left_wh,
        atom=OperatorAtom("chi_W_b21_chi"),
        localized=True,
        localization_controlled=True,
    )
    left, right = _blocks(
        case,
        left_factors=(localized_left, case.dilation),
        localization_controlled=True,
    )
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=_representation(case),
        right_block=right,
        assumptions=case.assumptions,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="cutoff rescaling",
    )


def test_incompatible_symbol_and_caller_contexts_return_a_blocked_result():
    case = _case()
    beta = sp.Symbol("beta", real=True)
    original_symbol = case.mellin_operator.symbol
    frequency = original_symbol.frequency
    spatial = original_symbol.spatial_variables[0]
    symbol_context = case.assumptions.with_assumption(beta < sp.pi)
    domain = MellinSymbolDomain(
        frequency=frequency,
        complex_domain=ComplexDomain.real_line(
            frequency.symbol,
            assumption_context=symbol_context,
        ),
        spatial_variables=(spatial,),
        assumption_context=symbol_context,
        singular_set=SingularSet(assumption_context=symbol_context),
        description="deliberately conflicting symbol context",
        evidence=("typed conflicting test domain",),
    )
    symbol = MellinSymbol(
        MellinExpression(
            original_symbol.as_expr(),
            domain,
            variables=(frequency, spatial),
            evidence=("conflicting-context test symbol",),
        ),
        MellinSymbolDependency.SPACE_FREQUENCY,
    )
    pdo = replace(case.mellin_operator, symbol=symbol)
    representation = _representation(case, mellin_operator=pdo)
    left, right = _blocks(case)
    caller_context = case.assumptions.with_assumption(beta > sp.pi)
    original = Product((left.atom, case.regularizer.operator, right.atom))

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=representation,
        right_block=right,
        assumptions=caller_context,
    )

    _assert_blocked(
        result,
        original=original,
        reason_fragment="explicitly incompatible",
    )
