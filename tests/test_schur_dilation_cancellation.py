from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from symbolic_operator_calculus import (
    AlgebraMembershipStatus,
    AssumptionContext,
    CancellationLayerStatus,
    ComplexDomain,
    DerivationRelationKind,
    ExactBlock,
    ExactIdentity,
    FormalIdentity,
    MellinExpression,
    MellinPseudodifferentialOperator,
    MellinSymbol,
    MellinSymbolDependency,
    MellinSymbolDomain,
    ModCompactEquivalence,
    OperatorAtom,
    Product,
    RegularizerMellinRepresentation,
    RegularizerOperator,
    RegularizerRepresentationStatus,
    SchurBlockRepresentation,
    SchurCancellationEvidence,
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
class _CancellationCase:
    gamma: sp.Symbol
    lam: sp.Symbol
    radial: sp.Symbol
    assumptions: AssumptionContext
    space: WeightedLpSpace
    left_wh: WienerHopfOperator
    right_wh: WienerHopfOperator
    dilation: WeightedDilationOperator
    inverse_dilation: WeightedDilationOperator
    mellin_operator: MellinPseudodifferentialOperator
    regularizer: RegularizerOperator


def _case() -> _CancellationCase:
    gamma = sp.Symbol("gamma")
    lam = sp.Symbol("lambda", real=True)
    radial = sp.Symbol("x", positive=True)
    assumptions = AssumptionContext((gamma > 0,))
    space = WeightedLpSpace(
        label="X_delta",
        p=sp.Integer(2),
        weight_exponent=sp.Integer(0),
        assumptions=assumptions,
        source="test realization of the article's multiplicative weight norm",
    )

    frequency = mellin_frequency(lam)
    spatial = output_spatial_variable(radial)
    domain = MellinSymbolDomain(
        frequency=frequency,
        complex_domain=ComplexDomain.real_line(lam),
        spatial_variables=(spatial,),
        assumption_context=assumptions,
        singular_set=SingularSet(),
        description="test domain for the radial regularizer symbol",
        evidence=("typed test domain",),
    )
    expression = MellinExpression(
        expression=radial + lam,
        domain=domain,
        variables=(frequency, spatial),
        evidence=("test symbol with space-frequency dependency",),
        description="r11",
    )
    symbol = MellinSymbol(
        expression,
        MellinSymbolDependency.SPACE_FREQUENCY,
        description="r11",
    )
    mellin_operator = MellinPseudodifferentialOperator(
        atom=OperatorAtom("Op_r11"),
        symbol=symbol,
        domain=space,
        codomain=space,
        symbol_class="tilde-E",
        source="test Mellin-PDO representation",
        radial_scaling_stable=True,
        stability_evidence=("direct radial-scaling proof",),
        evidence=("bounded test Mellin PDO",),
    )
    left_wh = WienerHopfOperator(
        atom=OperatorAtom("W_b21"),
        symbol=lam + 2,
        frequency_variable=lam,
        domain=space,
        codomain=space,
        symbol_class="weighted Fourier multipliers",
        source="test left Wiener--Hopf factor",
        frequency_scaling_stable=True,
        stability_evidence=("positive Fourier-frequency scaling",),
        evidence=("bounded left Wiener--Hopf factor",),
    )
    right_wh = WienerHopfOperator(
        atom=OperatorAtom("W_b12"),
        symbol=lam - 2,
        frequency_variable=lam,
        domain=space,
        codomain=space,
        symbol_class="weighted Fourier multipliers",
        source="test right Wiener--Hopf factor",
        frequency_scaling_stable=True,
        stability_evidence=("positive Fourier-frequency scaling",),
        evidence=("bounded right Wiener--Hopf factor",),
    )
    dilation = WeightedDilationOperator(
        atom=OperatorAtom("U_gamma"),
        gamma=gamma,
        domain=space,
        codomain=space,
        assumptions=assumptions,
        source="test normalized weighted dilation",
        evidence=("direct weighted norm calculation",),
    )
    inverse_dilation = dilation.inverse(OperatorAtom("U_gamma_inverse"))
    regularizer = RegularizerOperator(
        target=ExactBlock("A", 1, 1),
        operator=OperatorAtom("R11", kind="formal_regularizer"),
    )
    return _CancellationCase(
        gamma,
        lam,
        radial,
        assumptions,
        space,
        left_wh,
        right_wh,
        dilation,
        inverse_dilation,
        mellin_operator,
        regularizer,
    )


def _block_representations(
    case: _CancellationCase,
    *,
    relation: DerivationRelationKind,
) -> tuple[SchurBlockRepresentation, SchurBlockRepresentation]:
    left = SchurBlockRepresentation(
        exact_block=ExactBlock("A", 2, 1),
        atom=OperatorAtom("A21"),
        factors=(case.left_wh, case.dilation),
        relation_kind=relation,
        source="explicit test representation A21 = W(b21) U_gamma",
        hypotheses=("all left factors are bounded",),
        evidence=("left block factorization evidence",),
        compact_ideal="compact operators on X_delta",
    )
    right = SchurBlockRepresentation(
        exact_block=ExactBlock("A", 1, 2),
        atom=OperatorAtom("A12"),
        factors=(case.inverse_dilation, case.right_wh),
        relation_kind=relation,
        source="explicit test representation A12 = U_gamma_inverse W(b12)",
        hypotheses=("all right factors are bounded",),
        evidence=("right block factorization evidence",),
        compact_ideal="compact operators on X_delta",
    )
    return left, right


def _regularizer_representation(
    case: _CancellationCase,
    status: RegularizerRepresentationStatus,
) -> RegularizerMellinRepresentation:
    assumed = status is RegularizerRepresentationStatus.ASSUMED
    return RegularizerMellinRepresentation(
        regularizer=case.regularizer,
        mellin_operator=case.mellin_operator,
        status=status,
        hypotheses=(
            ("assume R11 is represented by the displayed Mellin PDO",)
            if assumed
            else ()
        ),
        source="explicit caller-supplied regularizer representation",
        evidence=() if assumed else ("regularizer representation evidence",),
        space=case.space,
        compact_ideal=(
            "compact operators on X_delta"
            if status is RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT
            else None
        ),
    )


def _obligation_keys(result: object) -> tuple[str, ...]:
    return tuple(item.key for item in result.proof_obligations)


def test_exact_cancellation_preserves_all_public_result_fields_and_four_layers():
    case = _case()
    left, right = _block_representations(
        case,
        relation=DerivationRelationKind.EXACT_EQUALITY,
    )
    representation = _regularizer_representation(
        case,
        RegularizerRepresentationStatus.CERTIFIED_EXACT,
    )

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=representation,
        right_block=right,
        assumptions=case.assumptions,
    )

    assert result.expression == Product(
        (
            case.left_wh.atom,
            OperatorAtom("Op_r11_radially_scaled"),
            case.right_wh.atom,
        )
    )
    assert result.original_expression == Product(
        (left.atom, case.regularizer.operator, right.atom)
    )
    assert result.model_expression == Product(
        (
            case.left_wh.atom,
            case.dilation.atom,
            case.mellin_operator.atom,
            case.inverse_dilation.atom,
            case.right_wh.atom,
        )
    )
    assert result.relation is DerivationRelationKind.EXACT_EQUALITY
    assert result.status is SchurCancellationStatus.EXACT_IDENTITY
    assert result.assumptions.contains_all(case.assumptions)
    assert result.assumptions.contains_exact(1 / case.gamma > 0)
    assert isinstance(result.semantic_relation, ExactIdentity)
    assert result.semantic_relation.left == result.original_expression
    assert result.semantic_relation.right == result.expression
    assert result.succeeded
    assert result.blocked_reasons == ()

    assert result.scaled_regularizer is not None
    assert result.scaled_regularizer.symbol.as_expr() == (
        case.gamma * case.radial + case.lam
    )
    assert case.radial / case.gamma not in sp.preorder_traversal(
        result.scaled_regularizer.symbol.as_expr()
    )
    assert r"\gamma x + \lambda" in result.latex
    assert result.latex.startswith(r"W\!\left(")
    assert r"\operatorname{Op}_M" in result.latex

    assert len(result.evidence) == 1
    evidence = result.evidence[0]
    assert isinstance(evidence, SchurCancellationEvidence)
    assert isinstance(evidence.exact_covariance, ExactIdentity)
    assert evidence.exact_covariance.left == Product(
        (
            case.dilation.atom,
            case.mellin_operator.atom,
            case.inverse_dilation.atom,
        )
    )
    assert evidence.compact_ideal_argument == ()
    assert _obligation_keys(result) == (
        "regularizer_mellin_representation",
        "wh_mellin_wh_sandwich_membership",
        "fredholm_algebra_membership",
        "schur_correction_symbol",
    )

    classification = result.classification
    assert classification.dilation_cancellation is (
        CancellationLayerStatus.EXACT_IDENTITY
    )
    assert classification.initial_block_representations is (
        CancellationLayerStatus.EXACT_IDENTITY
    )
    assert classification.regularizer_mellin_representation is (
        CancellationLayerStatus.EXACT_IDENTITY
    )
    assert classification.mixed_algebra_membership is (
        CancellationLayerStatus.UNPROVED_MIXED_ALGEBRA_MEMBERSHIP
    )
    assert result.algebra_membership.status is AlgebraMembershipStatus.UNPROVED
    assert result.algebra_membership.expression == result.expression


def test_mod_compact_inputs_produce_only_a_mod_compact_output_relation():
    case = _case()
    left, right = _block_representations(
        case,
        relation=DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
    )
    representation = _regularizer_representation(
        case,
        RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT,
    )

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=representation,
        right_block=right,
        assumptions=case.assumptions,
    )

    assert result.status is SchurCancellationStatus.EQUALITY_MODULO_COMPACTS
    assert result.relation is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
    assert isinstance(result.semantic_relation, ModCompactEquivalence)
    assert result.semantic_relation.left == result.original_expression
    assert result.semantic_relation.right == result.expression
    assert result.semantic_relation.space == case.space
    assert result.semantic_relation.compact_ideal == (
        "compact operators on X_delta"
    )
    assert result.scaled_regularizer is not None
    assert result.scaled_regularizer.symbol.as_expr() == (
        case.gamma * case.radial + case.lam
    )

    evidence = result.evidence[0]
    assert isinstance(evidence, SchurCancellationEvidence)
    assert evidence.compact_ideal_argument == (
        "AKB is compact for bounded A,B and compact K",
        "finite sums of compact operators are compact",
    )
    assert result.classification.dilation_cancellation is (
        CancellationLayerStatus.EXACT_IDENTITY
    )
    assert result.classification.initial_block_representations is (
        CancellationLayerStatus.EQUALITY_MODULO_COMPACTS
    )
    assert result.classification.regularizer_mellin_representation is (
        CancellationLayerStatus.EQUALITY_MODULO_COMPACTS
    )
    assert result.classification.mixed_algebra_membership is (
        CancellationLayerStatus.UNPROVED_MIXED_ALGEBRA_MEMBERSHIP
    )
    assert result.algebra_membership.status is AlgebraMembershipStatus.UNPROVED
    assert "wh_mellin_wh_sandwich_membership" in _obligation_keys(result)
    assert result.original_expression == Product(
        (left.atom, case.regularizer.operator, right.atom)
    )


def test_assumed_regularizer_keeps_the_cancellation_explicitly_conditional():
    case = _case()
    left, right = _block_representations(
        case,
        relation=DerivationRelationKind.EXACT_EQUALITY,
    )
    representation = _regularizer_representation(
        case,
        RegularizerRepresentationStatus.ASSUMED,
    )

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=representation,
        right_block=right,
        assumptions=case.assumptions,
    )

    assert result.status is (
        SchurCancellationStatus.CONDITIONAL_ON_REGULARIZER_REPRESENTATION
    )
    assert result.relation is DerivationRelationKind.EXACT_EQUALITY
    assert isinstance(result.semantic_relation, FormalIdentity)
    assert not isinstance(result.semantic_relation, ExactIdentity)
    assert not isinstance(result.semantic_relation, ModCompactEquivalence)
    assert isinstance(representation.semantic_relation, FormalIdentity)
    assert representation.hypotheses == (
        "assume R11 is represented by the displayed Mellin PDO",
    )
    assert representation.hypotheses[0] in result.semantic_relation.justification[3]
    assert result.classification.dilation_cancellation is (
        CancellationLayerStatus.EXACT_IDENTITY
    )
    assert result.classification.initial_block_representations is (
        CancellationLayerStatus.EXACT_IDENTITY
    )
    assert result.classification.regularizer_mellin_representation is (
        CancellationLayerStatus.CONDITIONAL_ON_REGULARIZER_REPRESENTATION
    )
    assert result.classification.mixed_algebra_membership is (
        CancellationLayerStatus.UNPROVED_MIXED_ALGEBRA_MEMBERSHIP
    )
    assert result.algebra_membership.status is AlgebraMembershipStatus.UNPROVED
    regularizer_obligation = result.proof_obligations[0]
    assert regularizer_obligation.key == "regularizer_mellin_representation"
    assert "assumed R11" in regularizer_obligation.statement
    assert result.original_expression == Product(
        (left.atom, case.regularizer.operator, right.atom)
    )
    assert result.expression != result.original_expression
    assert result.scaled_regularizer is not None
    assert result.scaled_regularizer.symbol.as_expr() == (
        case.gamma * case.radial + case.lam
    )


def test_assumed_regularizer_never_promotes_a_mod_compact_candidate_relation():
    case = _case()
    left, right = _block_representations(
        case,
        relation=DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
    )
    representation = _regularizer_representation(
        case,
        RegularizerRepresentationStatus.ASSUMED,
    )

    result = cancel_relative_dilations_in_schur_correction(
        left_block=left,
        regularizer=representation,
        right_block=right,
        assumptions=case.assumptions,
    )

    assert result.status is (
        SchurCancellationStatus.CONDITIONAL_ON_REGULARIZER_REPRESENTATION
    )
    assert result.relation is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
    assert isinstance(result.semantic_relation, FormalIdentity)
    assert not isinstance(result.semantic_relation, ExactIdentity)
    assert not isinstance(result.semantic_relation, ModCompactEquivalence)
    assert result.semantic_relation.left == result.original_expression
    assert result.semantic_relation.right == result.expression
    assert result.classification.initial_block_representations is (
        CancellationLayerStatus.EQUALITY_MODULO_COMPACTS
    )
    assert result.classification.regularizer_mellin_representation is (
        CancellationLayerStatus.CONDITIONAL_ON_REGULARIZER_REPRESENTATION
    )
    assert result.classification.mixed_algebra_membership is (
        CancellationLayerStatus.UNPROVED_MIXED_ALGEBRA_MEMBERSHIP
    )
    assert result.algebra_membership.status is AlgebraMembershipStatus.UNPROVED
