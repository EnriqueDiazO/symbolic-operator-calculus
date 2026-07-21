from __future__ import annotations

from dataclasses import dataclass

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    AssumptionContext,
    AuxiliaryOperator,
    CauchyPolarity,
    CauchyProjectionOperator,
    CoefficientMultiplicationOperator,
    ComplexDomain,
    CoreFactorizationStatus,
    CutoffEquivalenceClaim,
    CutoffEquivalenceSide,
    CutoffEquivalenceStatus,
    DerivationRelationKind,
    ExactBlock,
    ExactIdentity,
    FormalIdentity,
    FullSchurFactorizationStatus,
    FullSchurWordError,
    FullWordBlock,
    LinearCombination,
    LocalizedOperator,
    MellinConvolutionOperator,
    MellinExpression,
    MellinPseudodifferentialOperator,
    MellinSymbol,
    MellinSymbolDependency,
    MellinSymbolDomain,
    ModCompactEquivalence,
    MultiplicationOperator,
    OperatorAtom,
    Product,
    RawSchurCorrection,
    RegularizerInterfaceRole,
    RegularizerMellinRepresentation,
    RegularizerOperator,
    RegularizerRepresentationStatus,
    SchurBlockRepresentation,
    SchurPivot,
    SignedSchurContribution,
    SingularSet,
    SpatialCutoffOperator,
    Term,
    TransportedMellinCore,
    TransportedShiftOperator,
    WeightedDilationOperator,
    WeightedLpSpace,
    WienerHopfOperator,
    conjugate_cutoff_by_dilation,
    factor_full_first_schur_word,
    mellin_frequency,
    output_spatial_variable,
)


COMPACT_IDEAL = "compact operators on X_delta"


@dataclass(frozen=True)
class _FullWordCase:
    gamma: sp.Symbol
    lam: sp.Symbol
    radial: sp.Symbol
    assumptions: AssumptionContext
    space: WeightedLpSpace
    dilation: WeightedDilationOperator
    inverse_dilation: WeightedDilationOperator
    left_wh: WienerHopfOperator
    right_wh: WienerHopfOperator
    pdo: MellinPseudodifferentialOperator
    regularizer: RegularizerOperator
    cutoff: SpatialCutoffOperator
    left_localized: LocalizedOperator
    right_localized: LocalizedOperator
    pplus: CauchyProjectionOperator
    pminus: CauchyProjectionOperator


def _mellin_domain(
    *,
    lam: sp.Symbol,
    radial: sp.Symbol,
    assumptions: AssumptionContext,
    spatial: bool,
) -> MellinSymbolDomain:
    frequency = mellin_frequency(lam)
    return MellinSymbolDomain(
        frequency=frequency,
        complex_domain=ComplexDomain.real_line(lam),
        spatial_variables=(output_spatial_variable(radial),) if spatial else (),
        assumption_context=assumptions,
        singular_set=SingularSet(),
        description="typed test Mellin-symbol domain",
        evidence=("test domain evidence",),
    )


def _case() -> _FullWordCase:
    gamma = sp.Symbol("gamma")
    lam = sp.Symbol("lambda", real=True)
    radial = sp.Symbol("x", positive=True)
    assumptions = AssumptionContext((gamma > 0,))
    space = WeightedLpSpace(
        label="X_delta",
        p=sp.Integer(2),
        weight_exponent=sp.Integer(0),
        assumptions=assumptions,
        source="article multiplicative-weight convention",
    )
    dilation = WeightedDilationOperator(
        atom=OperatorAtom("V_gamma"),
        gamma=gamma,
        domain=space,
        codomain=space,
        assumptions=assumptions,
        source="normalized dilation",
        evidence=("direct norm calculation",),
    )
    inverse = dilation.inverse(OperatorAtom("V_gamma_inverse"))
    left_wh = WienerHopfOperator(
        atom=OperatorAtom("W_left"),
        symbol=lam + 2,
        frequency_variable=lam,
        domain=space,
        codomain=space,
        symbol_class="weighted Fourier multipliers",
        source="synthetic left exterior",
        frequency_scaling_stable=True,
        stability_evidence=("positive frequency scaling",),
        evidence=("bounded left WH operator",),
    )
    right_wh = WienerHopfOperator(
        atom=OperatorAtom("W_right"),
        symbol=lam - 2,
        frequency_variable=lam,
        domain=space,
        codomain=space,
        symbol_class="weighted Fourier multipliers",
        source="synthetic right exterior",
        frequency_scaling_stable=True,
        stability_evidence=("positive frequency scaling",),
        evidence=("bounded right WH operator",),
    )
    pdo_expression = MellinExpression(
        expression=radial + lam,
        domain=_mellin_domain(
            lam=lam,
            radial=radial,
            assumptions=assumptions,
            spatial=True,
        ),
        variables=(mellin_frequency(lam), output_spatial_variable(radial)),
        evidence=("space-frequency symbol evidence",),
        description="r11",
    )
    pdo = MellinPseudodifferentialOperator(
        atom=OperatorAtom("Op_r11"),
        symbol=MellinSymbol(
            pdo_expression,
            MellinSymbolDependency.SPACE_FREQUENCY,
            description="r11",
        ),
        domain=space,
        codomain=space,
        symbol_class="tilde-E",
        source="typed regularizer symbol",
        radial_scaling_stable=True,
        stability_evidence=("exact radial scaling",),
        evidence=("bounded Mellin PDO",),
    )
    regularizer = RegularizerOperator(
        target=ExactBlock("A", 1, 1),
        operator=OperatorAtom("R1", kind="formal_regularizer"),
    )
    cutoff = SpatialCutoffOperator(
        atom=OperatorAtom("M_chi_infinity"),
        cutoff_symbol=sp.Function("chi_infinity")(radial),
        radial_variable=radial,
        domain=space,
        codomain=space,
        support_region="cusp neighborhood x >= 1",
        equals_one_near="infinity, concretely x >= 2",
        equals_zero_away_from="the cusp, concretely 0 <= x <= 1",
        coordinate_system="transported half-line radial coordinate x",
        radial_scale=sp.Integer(1),
        source="article section 03, lines 747-762",
        evidence=("verified compact localization source",),
    )
    left_localized = LocalizedOperator(
        atom=OperatorAtom("W_left_localized"),
        left_cutoff=cutoff,
        core=left_wh,
        right_cutoff=cutoff,
        source="two-cutoff left Wiener--Hopf model",
        evidence=("article model retains both cuts",),
    )
    right_localized = LocalizedOperator(
        atom=OperatorAtom("W_right_localized"),
        left_cutoff=cutoff,
        core=right_wh,
        right_cutoff=cutoff,
        source="two-cutoff right Wiener--Hopf model",
        evidence=("article model retains both cuts",),
    )

    frequency_expression = MellinExpression(
        expression=lam,
        domain=_mellin_domain(
            lam=lam,
            radial=radial,
            assumptions=assumptions,
            spatial=False,
        ),
        variables=(mellin_frequency(lam),),
        evidence=("frequency-only Cauchy symbol",),
        description="p plus",
    )
    pplus_mellin = MellinConvolutionOperator(
        atom=OperatorAtom("Op_p_plus"),
        symbol=MellinSymbol(
            frequency_expression,
            MellinSymbolDependency.FREQUENCY_ONLY,
            description="p plus",
        ),
        domain=space,
        codomain=space,
        symbol_class="V(R)",
        source="half-line Cauchy Mellin representation",
        evidence=("article section 03, lines 1154-1197",),
    )
    pminus_mellin = MellinConvolutionOperator(
        atom=OperatorAtom("Op_p_minus"),
        symbol=pplus_mellin.symbol,
        domain=space,
        codomain=space,
        symbol_class="V(R)",
        source="half-line Cauchy Mellin representation",
        evidence=("article section 03, lines 1154-1197",),
    )
    pplus = CauchyProjectionOperator(
        atom=OperatorAtom("P_plus"),
        polarity=CauchyPolarity.PLUS,
        mellin_operator=pplus_mellin,
        source="article Cauchy factor P+",
        evidence=("exact Mellin-convolution representation",),
    )
    pminus = CauchyProjectionOperator(
        atom=OperatorAtom("P_minus"),
        polarity=CauchyPolarity.MINUS,
        mellin_operator=pminus_mellin,
        source="article Cauchy factor P-",
        evidence=("exact Mellin-convolution representation",),
    )
    return _FullWordCase(
        gamma,
        lam,
        radial,
        assumptions,
        space,
        dilation,
        inverse,
        left_wh,
        right_wh,
        pdo,
        regularizer,
        cutoff,
        left_localized,
        right_localized,
        pplus,
        pminus,
    )


def _regularizer_representation(
    case: _FullWordCase,
    status: RegularizerRepresentationStatus,
) -> RegularizerMellinRepresentation:
    assumed = status is RegularizerRepresentationStatus.ASSUMED
    return RegularizerMellinRepresentation(
        regularizer=case.regularizer,
        mellin_operator=case.pdo,
        status=status,
        hypotheses=("assume the displayed R1 representation",) if assumed else (),
        source="article R1=Phi^-1 Op_M(r11) Phi",
        evidence=() if assumed else ("article section 06, lines 161-168",),
        space=case.space,
        compact_ideal=(
            COMPACT_IDEAL
            if status is RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT
            else None
        ),
    )


def _transported_core(
    case: _FullWordCase,
    status: RegularizerRepresentationStatus = (
        RegularizerRepresentationStatus.CERTIFIED_EXACT
    ),
) -> TransportedMellinCore:
    return TransportedMellinCore(
        representation=_regularizer_representation(case, status),
        phi_inverse=OperatorAtom("Phi_delta_inverse", kind="space_isomorphism"),
        phi=OperatorAtom("Phi_delta", kind="space_isomorphism"),
        source="exact transported Mellin interface",
        evidence=("article section 06, lines 161-168",),
    )


def _blocks(
    case: _FullWordCase,
    *,
    relation: DerivationRelationKind = DerivationRelationKind.EXACT_EQUALITY,
    left_exterior: object | None = None,
    right_exterior: object | None = None,
    left_claims: tuple[CutoffEquivalenceClaim, ...] = (),
    right_claims: tuple[CutoffEquivalenceClaim, ...] = (),
) -> tuple[FullWordBlock, FullWordBlock]:
    left = SchurBlockRepresentation(
        exact_block=ExactBlock("A", 2, 1),
        atom=OperatorAtom("B21_primitive"),
        factors=(left_exterior or case.left_wh,),
        relation_kind=relation,
        source="explicit left primitive representation",
        hypotheses=("bounded left primitive",),
        evidence=("left factorization evidence",),
        compact_ideal=COMPACT_IDEAL,
    )
    right = SchurBlockRepresentation(
        exact_block=ExactBlock("A", 1, 2),
        atom=OperatorAtom("B12_primitive"),
        factors=(right_exterior or case.right_wh,),
        relation_kind=relation,
        source="explicit right primitive representation",
        hypotheses=("bounded right primitive",),
        evidence=("right factorization evidence",),
        compact_ideal=COMPACT_IDEAL,
    )
    return FullWordBlock(left, cutoff_claims=left_claims), FullWordBlock(
        right,
        cutoff_claims=right_claims,
    )


def _factor(
    case: _FullWordCase,
    left: FullWordBlock,
    right: FullWordBlock,
    *,
    left_auxiliary: object | None = None,
    right_auxiliary: object | None = None,
    core: object | None = None,
):
    return factor_full_first_schur_word(
        left_block=left,
        left_auxiliary=left_auxiliary or case.dilation,
        mellin_core=core or _transported_core(case),
        right_auxiliary=right_auxiliary or case.inverse_dilation,
        right_block=right,
        assumptions=case.assumptions,
    )


def _cutoff_claim(
    case: _FullWordCase,
    localized: LocalizedOperator,
    *,
    status: CutoffEquivalenceStatus,
    side: CutoffEquivalenceSide = CutoffEquivalenceSide.BOTH,
) -> CutoffEquivalenceClaim:
    scaled = conjugate_cutoff_by_dilation(
        case.dilation,
        case.cutoff,
        case.inverse_dilation,
    ).scaled_cutoff
    return CutoffEquivalenceClaim(
        original_cutoff=case.cutoff,
        replacement_cutoff=scaled,
        localized_operator=localized,
        status=status,
        ideal=COMPACT_IDEAL,
        side=side,
        assumptions=("side-specific compact cutoff difference",),
        source="caller-supplied cutoff replacement lemma",
        evidence=("certified cutoff evidence",)
        if status is CutoffEquivalenceStatus.CERTIFIED
        else (),
    )


def test_cutoff_covariance_rescales_symbol_and_preserves_geometric_metadata():
    case = _case()

    trace = conjugate_cutoff_by_dilation(
        case.dilation,
        case.cutoff,
        case.inverse_dilation,
    )

    assert trace.scaled_cutoff.cutoff_symbol == sp.Function("chi_infinity")(
        case.gamma * case.radial
    )
    assert trace.scaled_cutoff.radial_scale == case.gamma
    assert trace.scaled_cutoff.coordinate_system == case.cutoff.coordinate_system
    assert "1/gamma" in trace.scaled_cutoff.equals_one_near
    assert trace.source_cutoff == case.cutoff
    assert trace.exact_identity.left == Product(
        (case.dilation.atom, case.cutoff.atom, case.inverse_dilation.atom)
    )
    assert trace.exact_identity.right == Product((trace.scaled_cutoff.atom,))
    assert trace.scaled_cutoff != case.cutoff


def test_cutoff_is_never_a_boolean_or_untyped_multiplier():
    case = _case()
    with pytest.raises(TypeError, match="cutoff_symbol"):
        SpatialCutoffOperator(
            atom=OperatorAtom("bad_cutoff"),
            cutoff_symbol=True,  # type: ignore[arg-type]
            radial_variable=case.radial,
            domain=case.space,
            codomain=case.space,
            support_region="region",
            equals_one_near="infinity",
            equals_zero_away_from="zero",
            coordinate_system="x",
            radial_scale=sp.Integer(1),
            source="test",
        )
    with pytest.raises(TypeError, match="SpatialCutoffOperator"):
        conjugate_cutoff_by_dilation(
            case.dilation,
            case.cutoff.multiplication,  # type: ignore[arg-type]
            case.inverse_dilation,
        )


def test_cutoff_claim_strength_is_explicit_and_never_inferred():
    case = _case()
    unproved = _cutoff_claim(
        case,
        case.left_localized,
        status=CutoffEquivalenceStatus.UNPROVED,
    )
    assumed = _cutoff_claim(
        case,
        case.left_localized,
        status=CutoffEquivalenceStatus.ASSUMED,
    )
    certified = _cutoff_claim(
        case,
        case.left_localized,
        status=CutoffEquivalenceStatus.CERTIFIED,
    )

    assert unproved.semantic_relation is None
    assert isinstance(assumed.semantic_relation, FormalIdentity)
    assert not isinstance(assumed.semantic_relation, ModCompactEquivalence)
    assert isinstance(certified.semantic_relation, ModCompactEquivalence)
    assert certified.semantic_relation.compact_ideal == COMPACT_IDEAL
    assert case.left_localized.exact_definition.right == Product(
        (case.cutoff.atom, case.left_wh.atom, case.cutoff.atom)
    )


def test_typed_coefficient_shift_projection_and_auxiliary_preserve_order():
    case = _case()
    multiplier = MultiplicationOperator(
        atom=OperatorAtom("M_rho"),
        symbol=(case.gamma * case.radial + 1) / (case.radial + 1),
        radial_variable=case.radial,
        domain=case.space,
        codomain=case.space,
        source="transport coefficient rho",
        evidence=("endpoint values",),
    )
    coefficient = CoefficientMultiplicationOperator(multiplier, "rho_1")
    shift = TransportedShiftOperator(
        atom=OperatorAtom("Vtilde_1"),
        coefficient=coefficient,
        dilation=case.dilation,
        source="Vtilde_1=M_rho V_gamma",
        evidence=("article exact factorization",),
    )
    auxiliary_expression = LinearCombination(
        (
            Term(1, Product((case.inverse_dilation.atom, case.pplus.atom))),
            Term(1, case.pminus.atom),
        )
    )
    auxiliary = AuxiliaryOperator(
        atom=OperatorAtom("T1_minus"),
        expression=auxiliary_expression,
        components=(case.inverse_dilation, case.pplus, case.pminus),
        domain=case.space,
        codomain=case.space,
        source="T1_minus=U1_inverse P_plus+P_minus",
        evidence=("article section 05, lines 522-529",),
    )

    assert shift.exact_factorization.right == Product(
        (coefficient.atom, case.dilation.atom)
    )
    assert auxiliary.exact_definition.right == auxiliary_expression
    assert auxiliary.expression.terms[0].product.factors == (
        case.inverse_dilation.atom,
        case.pplus.atom,
    )
    assert auxiliary.expression.terms[1].product.factors == (case.pminus.atom,)


def test_synthetic_exact_full_word_factors_and_delegates_cancellation():
    case = _case()
    left, right = _blocks(case)

    result = _factor(case, left, right)

    assert result.status is (
        FullSchurFactorizationStatus.FULL_WORD_EXACT_FACTORIZATION
    )
    assert result.original_word == Product(
        (
            left.atom,
            case.dilation.atom,
            case.regularizer.operator,
            case.inverse_dilation.atom,
            right.atom,
        )
    )
    assert result.expanded_original_word == Product(
        (
            case.left_wh.atom,
            case.dilation.atom,
            OperatorAtom("Phi_delta_inverse", kind="space_isomorphism"),
            case.pdo.atom,
            OperatorAtom("Phi_delta", kind="space_isomorphism"),
            case.inverse_dilation.atom,
            case.right_wh.atom,
        )
    )
    assert result.factorized_word == Product(
        (
            case.left_wh.atom,
            OperatorAtom("Op_r11_radially_scaled"),
            case.right_wh.atom,
        )
    )
    assert isinstance(result.relation, ExactIdentity)
    assert result.cancellation_result is not None
    assert result.cancellation_result.succeeded
    assert result.left_core.status is CoreFactorizationStatus.EXACT
    assert result.right_core.status is CoreFactorizationStatus.EXACT
    assert result.succeeded


def test_synthetic_mod_compact_word_preserves_common_ideal_and_evidence():
    case = _case()
    left, right = _blocks(
        case,
        relation=DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
    )

    result = _factor(
        case,
        left,
        right,
        core=_transported_core(
            case,
            RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT,
        ),
    )

    assert result.status is (
        FullSchurFactorizationStatus.FULL_WORD_MOD_COMPACTS_FACTORIZATION
    )
    assert isinstance(result.relation, ModCompactEquivalence)
    assert result.relation.compact_ideal == COMPACT_IDEAL
    assert result.left_core.status is CoreFactorizationStatus.MODULO_COMPACTS
    assert result.right_core.status is CoreFactorizationStatus.MODULO_COMPACTS
    assert result.assumptions.contains_exact(case.gamma > 0)
    assert result.evidence
    assert result.cancellation_result is not None
    assert result.cancellation_result.succeeded


@pytest.mark.parametrize(
    ("claim_status", "relation_type"),
    [
        (CutoffEquivalenceStatus.ASSUMED, FormalIdentity),
        (CutoffEquivalenceStatus.CERTIFIED, ModCompactEquivalence),
    ],
)
def test_localized_factorization_requires_and_preserves_cutoff_claim_strength(
    claim_status: CutoffEquivalenceStatus,
    relation_type: type[FormalIdentity] | type[ModCompactEquivalence],
):
    case = _case()
    left_claim = _cutoff_claim(
        case,
        case.left_localized,
        status=claim_status,
    )
    right_claim = _cutoff_claim(
        case,
        case.right_localized,
        status=claim_status,
    )
    left, right = _blocks(
        case,
        left_exterior=case.left_localized,
        right_exterior=case.right_localized,
        left_claims=(left_claim,),
        right_claims=(right_claim,),
    )

    result = _factor(case, left, right)

    assert result.status is (
        FullSchurFactorizationStatus.FACTORIZATION_AFTER_CUTOFF_EQUIVALENCE
    )
    assert isinstance(result.relation, relation_type)
    assert result.factorized_word is not None
    assert result.factorized_word.factors[0] == case.left_localized.atom
    assert result.factorized_word.factors[-1] == case.right_localized.atom
    assert result.cancellation_result is None
    assert left_claim.semantic_relation in result.evidence
    assert right_claim.semantic_relation in result.evidence


def test_localized_word_without_claim_is_blocked_by_cutoff_control():
    case = _case()
    left, right = _blocks(
        case,
        left_exterior=case.left_localized,
        right_exterior=case.right_localized,
    )

    result = _factor(case, left, right)

    assert result.status is FullSchurFactorizationStatus.BLOCKED_BY_CUTOFF_CONTROL
    assert result.factorized_word is None
    assert result.relation is None
    assert result.blocking_factors == (
        case.left_localized,
        case.right_localized,
    )


def test_one_sided_claim_does_not_control_a_two_sided_localization():
    case = _case()
    left_claim = _cutoff_claim(
        case,
        case.left_localized,
        status=CutoffEquivalenceStatus.CERTIFIED,
        side=CutoffEquivalenceSide.LEFT,
    )
    right_claim = _cutoff_claim(
        case,
        case.right_localized,
        status=CutoffEquivalenceStatus.CERTIFIED,
        side=CutoffEquivalenceSide.LEFT,
    )
    left, right = _blocks(
        case,
        left_exterior=case.left_localized,
        right_exterior=case.right_localized,
        left_claims=(left_claim,),
        right_claims=(right_claim,),
    )

    result = _factor(case, left, right)

    assert result.status is FullSchurFactorizationStatus.BLOCKED_BY_CUTOFF_CONTROL


def test_projection_between_exterior_and_dilation_blocks_operator_order():
    case = _case()
    left_representation = SchurBlockRepresentation(
        exact_block=ExactBlock("A", 2, 1),
        atom=OperatorAtom("left_with_projection"),
        factors=(case.left_wh, case.pplus),
        relation_kind=DerivationRelationKind.EXACT_EQUALITY,
        source="projection remains between WH and dilation",
        hypotheses=(),
        evidence=("ordered source AST",),
    )
    right = _blocks(case)[1]

    result = _factor(case, FullWordBlock(left_representation), right)

    assert result.status is FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER
    assert case.pplus in result.left_core.blocking_factors
    assert result.original_word.factors[1] == case.dilation.atom
    assert result.factorized_word is None


def test_projection_between_inverse_dilation_and_right_exterior_blocks_order():
    case = _case()
    right_representation = SchurBlockRepresentation(
        exact_block=ExactBlock("A", 1, 2),
        atom=OperatorAtom("right_with_projection"),
        factors=(case.pminus, case.right_wh),
        relation_kind=DerivationRelationKind.EXACT_EQUALITY,
        source="projection remains before the right WH exterior",
        hypotheses=(),
        evidence=("ordered source AST",),
    )
    left = _blocks(case)[0]

    result = _factor(case, left, FullWordBlock(right_representation))

    assert result.status is FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER
    assert case.pminus in result.right_core.blocking_factors
    assert result.factorized_word is None


def test_domain_mismatch_blocks_left_core_without_retyping_spaces():
    case = _case()
    other_space = WeightedLpSpace(
        label="X_other",
        p=sp.Integer(3),
        weight_exponent=sp.Integer(0),
        assumptions=case.assumptions,
        source="different test domain",
    )
    other_wh = WienerHopfOperator(
        atom=OperatorAtom("W_other_space"),
        symbol=case.lam,
        frequency_variable=case.lam,
        domain=other_space,
        codomain=other_space,
        symbol_class="weighted Fourier multipliers",
        source="different-domain exterior",
        frequency_scaling_stable=True,
        stability_evidence=("scale stable",),
        evidence=("bounded",),
    )
    left, right = _blocks(case, left_exterior=other_wh)

    result = _factor(case, left, right)

    assert result.status is FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER
    assert result.left_core.status is CoreFactorizationStatus.BLOCKED
    assert other_wh in result.left_core.blocking_factors


def test_domain_mismatch_blocks_right_core_without_forced_symmetry():
    case = _case()
    other_space = WeightedLpSpace(
        label="X_other_right",
        p=sp.Integer(3),
        weight_exponent=sp.Integer(0),
        assumptions=case.assumptions,
        source="different right test domain",
    )
    other_wh = WienerHopfOperator(
        atom=OperatorAtom("W_other_right_space"),
        symbol=case.lam,
        frequency_variable=case.lam,
        domain=other_space,
        codomain=other_space,
        symbol_class="weighted Fourier multipliers",
        source="different-domain right exterior",
        frequency_scaling_stable=True,
        stability_evidence=("scale stable",),
        evidence=("bounded",),
    )
    left, right = _blocks(case, right_exterior=other_wh)

    result = _factor(case, left, right)

    assert result.status is FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER
    assert result.right_core.status is CoreFactorizationStatus.BLOCKED
    assert other_wh in result.right_core.blocking_factors


def test_wrong_right_dilation_orientation_is_not_cancelled():
    case = _case()
    left, right = _blocks(case)

    result = _factor(case, left, right, right_auxiliary=case.dilation)

    assert result.status is FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER
    assert result.cancellation_result is None
    assert result.factorized_word is None
    assert case.dilation in result.blocking_factors


def test_wrong_left_dilation_orientation_is_not_cancelled():
    case = _case()
    left, right = _blocks(case)

    result = _factor(case, left, right, left_auxiliary=case.inverse_dilation)

    assert result.status is FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER
    assert result.cancellation_result is None
    assert result.factorized_word is None
    assert case.inverse_dilation in result.blocking_factors


def test_assumed_regularizer_is_blocked_at_the_regularizer_interface():
    case = _case()
    left, right = _blocks(case)
    assumed = _transported_core(
        case,
        RegularizerRepresentationStatus.ASSUMED,
    )

    result = _factor(case, left, right, core=assumed)

    assert result.status is (
        FullSchurFactorizationStatus.BLOCKED_BY_REGULARIZER_INTERFACE
    )
    assert result.blocking_factors == (assumed,)
    assert result.original_word.factors[2] == case.regularizer.operator


def test_untransported_regularizer_representation_is_blocked_not_collapsed():
    case = _case()
    left, right = _blocks(case)
    representation = _regularizer_representation(
        case,
        RegularizerRepresentationStatus.CERTIFIED_EXACT,
    )

    result = _factor(case, left, right, core=representation)

    assert result.status is (
        FullSchurFactorizationStatus.BLOCKED_BY_REGULARIZER_INTERFACE
    )
    assert result.original_word.factors[2] == case.regularizer.operator
    assert case.pdo.atom not in result.original_word.factors
    assert result.factorized_word is None


def test_literal_article_labels_are_rejected_for_typed_duplicate_interfaces():
    case = _case()
    t1_minus = OperatorAtom("T1_minus")
    z1_inverse = OperatorAtom("Z1_inverse")
    left_representation = SchurBlockRepresentation(
        exact_block=ExactBlock("A", 2, 1),
        atom=OperatorAtom("B21_minus"),
        factors=(
            OperatorAtom("Vtilde2_minus_G2"),
            case.left_localized,
            t1_minus,
        ),
        relation_kind=DerivationRelationKind.EXACT_EQUALITY,
        source="article section 06, lines 382-390",
        hypotheses=(),
        evidence=("literal article definition",),
    )
    right_representation = SchurBlockRepresentation(
        exact_block=ExactBlock("A", 1, 2),
        atom=OperatorAtom("B12_plus"),
        factors=(
            z1_inverse,
            OperatorAtom("Vtilde1_minus_G1"),
            case.right_localized,
        ),
        relation_kind=DerivationRelationKind.EXACT_EQUALITY,
        source="article section 06, lines 391-399",
        hypotheses=(),
        evidence=("literal article definition",),
    )
    left = FullWordBlock(
        left_representation,
        embedded_interfaces=(RegularizerInterfaceRole.LEFT_AUXILIARY,),
    )
    right = FullWordBlock(
        right_representation,
        embedded_interfaces=(RegularizerInterfaceRole.RIGHT_AUXILIARY,),
    )

    result = _factor(
        case,
        left,
        right,
        left_auxiliary=t1_minus,
        right_auxiliary=z1_inverse,
    )

    assert result.status is (
        FullSchurFactorizationStatus.BLOCKED_BY_REGULARIZER_INTERFACE
    )
    assert result.original_word == Product(
        (
            left.atom,
            t1_minus,
            case.regularizer.operator,
            z1_inverse,
            right.atom,
        )
    )
    assert result.expanded_original_word is not None
    assert result.expanded_original_word.factors.count(t1_minus) == 2
    assert result.expanded_original_word.factors.count(z1_inverse) == 2
    assert result.factorized_word is None
    assert result.relation is None


def test_interface_detection_uses_typed_roles_not_coincidental_names():
    case = _case()
    coincidental = OperatorAtom("T1_minus")
    left_representation = SchurBlockRepresentation(
        exact_block=ExactBlock("A", 2, 1),
        atom=OperatorAtom("looks_like_article_block"),
        factors=(coincidental,),
        relation_kind=DerivationRelationKind.EXACT_EQUALITY,
        source="synthetic coincidental label",
        hypotheses=(),
        evidence=("typed role deliberately absent",),
    )
    right = _blocks(case)[1]

    result = _factor(case, FullWordBlock(left_representation), right)

    assert result.status is FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER
    assert result.status is not (
        FullSchurFactorizationStatus.BLOCKED_BY_REGULARIZER_INTERFACE
    )


def test_latex_retains_status_and_every_named_obligation():
    case = _case()
    left, right = _blocks(
        case,
        left_exterior=case.left_localized,
        right_exterior=case.right_localized,
    )

    result = _factor(case, left, right)
    keys = tuple(item.key for item in result.proof_obligations)

    assert keys == (
        "left_core_factorization",
        "right_core_factorization",
        "cutoff_replacement_mod_compacts",
        "regularizer_mellin_representation",
        "wh_mellin_wh_sandwich_membership",
        "fredholm_algebra_membership",
        "schur_correction_symbol",
    )
    assert r"BLOCKED\_BY\_CUTOFF\_CONTROL" in result.latex
    for key in keys:
        assert key.replace("_", r"\_") in result.latex


def test_raw_correction_signed_contribution_and_pivot_keep_sign_layers_separate():
    raw = RawSchurCorrection(
        OperatorAtom("A21"),
        OperatorAtom("A11_regularizer"),
        OperatorAtom("A12"),
    )
    contribution = SignedSchurContribution(raw, -1)
    pivot = SchurPivot(OperatorAtom("A22"), contribution)

    assert raw.intrinsic_sign == 1
    assert raw.word == Product(
        (OperatorAtom("A21"), OperatorAtom("A11_regularizer"), OperatorAtom("A12"))
    )
    assert contribution.term.coefficient == -1
    assert pivot.expression.terms[0].coefficient == 1
    assert pivot.expression.terms[0].product == Product((OperatorAtom("A22"),))
    assert pivot.expression.terms[1] == contribution.term


@pytest.mark.parametrize("bad_sign", [-1, 0, 2, True])
def test_raw_correction_rejects_embedded_sign_and_double_negation(bad_sign: object):
    with pytest.raises(FullSchurWordError, match="intrinsically positive"):
        RawSchurCorrection(
            OperatorAtom("A21"),
            OperatorAtom("R11"),
            OperatorAtom("A12"),
            intrinsic_sign=bad_sign,  # type: ignore[arg-type]
        )


@pytest.mark.parametrize("bad_sign", [-2, 0, 2, True])
def test_signed_contribution_accepts_only_one_external_sign(bad_sign: object):
    raw = RawSchurCorrection(
        OperatorAtom("A21"),
        OperatorAtom("R11"),
        OperatorAtom("A12"),
    )
    with pytest.raises(FullSchurWordError, match="exactly"):
        SignedSchurContribution(raw, bad_sign)  # type: ignore[arg-type]
