from __future__ import annotations

from dataclasses import dataclass, replace

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    AssumptionContext,
    CompactIdeal,
    CompactIdealMismatch,
    ComplexDomain,
    CoreRegularizerClaimStatus,
    CoreRegularizerProperty,
    DiagonalLinkIdentity,
    DiagonalLinkStatus,
    DiagonalOperator,
    DomainOrOrderError,
    ExactIdentity,
    ExactIdentityScope,
    FactoredRegularizerStatus,
    FormalIdentity,
    IDENTITY_PRODUCT,
    InvertibleAuxiliaryOperator,
    InvertibleMultiplicationOperator,
    LinearCombination,
    MellinCoreOperator,
    MellinExpression,
    MellinPseudodifferentialOperator,
    MellinSymbol,
    MellinSymbolDependency,
    MellinSymbolDomain,
    MissingRegularizerEvidence,
    ModCompactEquivalence,
    MultiplicationOperator,
    OperatorAtom,
    Product,
    ProductCertificationStatus,
    RegularizerMellinRepresentation,
    RegularizerOperator,
    RegularizerRepresentationStatus,
    SchurInsertionStatus,
    SchurSymbolStatus,
    SingularSet,
    Term,
    TransportedMellinCore,
    UnprovedClaimStatus,
    WeightedLpSpace,
    certify_factored_two_sided_diagonal_regularizer,
    insert_factored_regularizer_in_schur,
    make_mellin_core_regularizer_certificate,
    mellin_frequency,
    output_spatial_variable,
    propagate_compact_ideal,
)
from symbolic_operator_calculus.schur_full_word import AuxiliaryOperator


@dataclass(frozen=True)
class _BoundedOperator:
    atom: OperatorAtom
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    bounded: bool = True


@dataclass(frozen=True)
class _Case:
    assumptions: AssumptionContext
    space: WeightedLpSpace
    other_space: WeightedLpSpace
    diagonal: DiagonalOperator
    auxiliary: InvertibleAuxiliaryOperator
    core: MellinCoreOperator
    transported_regularizer: TransportedMellinCore
    multiplier: InvertibleMultiplicationOperator
    ideal: CompactIdeal
    left_exterior: _BoundedOperator
    right_exterior: _BoundedOperator


def _mellin_pdo(
    *,
    name: str,
    expression: sp.Expr,
    radial: sp.Symbol,
    frequency: sp.Symbol,
    assumptions: AssumptionContext,
    space: WeightedLpSpace,
) -> MellinPseudodifferentialOperator:
    typed_frequency = mellin_frequency(frequency)
    typed_radial = output_spatial_variable(radial)
    domain = MellinSymbolDomain(
        frequency=typed_frequency,
        complex_domain=ComplexDomain.real_line(frequency),
        spatial_variables=(typed_radial,),
        assumption_context=assumptions,
        singular_set=SingularSet(),
        description=f"domain for {name}",
        evidence=("typed Mellin domain",),
    )
    mellin_expression = MellinExpression(
        expression=expression,
        domain=domain,
        variables=(typed_frequency, typed_radial),
        evidence=(f"formula for {name}",),
        description=name,
    )
    symbol = MellinSymbol(
        mellin_expression,
        MellinSymbolDependency.SPACE_FREQUENCY,
        description=name,
    )
    return MellinPseudodifferentialOperator(
        atom=OperatorAtom(f"Op_{name}"),
        symbol=symbol,
        domain=space,
        codomain=space,
        symbol_class="tilde-E",
        source=f"test {name} Mellin PDO",
        radial_scaling_stable=True,
        stability_evidence=("radial scaling theorem",),
        evidence=("bounded Mellin PDO",),
    )


def _case() -> _Case:
    gamma = sp.Symbol("gamma", positive=True)
    kappa = sp.Symbol("kappa", real=True)
    radial = sp.Symbol("x", positive=True)
    frequency = sp.Symbol("lambda", real=True)
    assumptions = AssumptionContext((gamma > 0, kappa > 0, kappa < 1))
    space = WeightedLpSpace(
        label="X_delta",
        p=sp.Integer(2),
        weight_exponent=sp.Integer(0),
        assumptions=assumptions,
        source="article scalar weighted space",
    )
    other_space = WeightedLpSpace(
        label="Y_delta",
        p=sp.Integer(2),
        weight_exponent=sp.Integer(0),
        assumptions=assumptions,
        source="intentionally different weighted-space label",
    )
    diagonal = DiagonalOperator(
        atom=OperatorAtom("A11"),
        domain=space,
        codomain=space,
        source="article section 06:27-42",
        evidence=("A11=Vtilde_1 P_plus+G1 P_minus",),
    )
    u_inverse = MultiplicationOperator(
        atom=OperatorAtom("U1_inverse"),
        symbol=sp.Integer(1),
        radial_variable=radial,
        domain=space,
        codomain=space,
        source="typed test stand-in for exact U1 inverse",
        evidence=("U1 theorem",),
    )
    p_plus = MultiplicationOperator(
        atom=OperatorAtom("P_plus"),
        symbol=sp.Integer(1),
        radial_variable=radial,
        domain=space,
        codomain=space,
        source="typed P+ component",
        evidence=("Cauchy projection theorem",),
    )
    p_minus = MultiplicationOperator(
        atom=OperatorAtom("P_minus"),
        symbol=sp.Integer(1),
        radial_variable=radial,
        domain=space,
        codomain=space,
        source="typed P- component",
        evidence=("Cauchy projection theorem",),
    )
    t_minus = AuxiliaryOperator(
        atom=OperatorAtom("T1_minus"),
        expression=LinearCombination(
            (
                Term(1, Product((u_inverse.atom, p_plus.atom))),
                Term(1, p_minus.atom),
            )
        ),
        components=(u_inverse, p_plus, p_minus),
        domain=space,
        codomain=space,
        source="article section 05:515-529",
        evidence=("T1_minus=U1_inverse P_plus+P_minus",),
    )
    t_inverse_component = MultiplicationOperator(
        atom=OperatorAtom("Op_theta1_inverse"),
        symbol=sp.Integer(1),
        radial_variable=radial,
        domain=space,
        codomain=space,
        source="Mellin inverse component",
        evidence=("article section 05:770-778",),
    )
    t_minus_inverse = AuxiliaryOperator(
        atom=OperatorAtom("T1_minus_inverse"),
        expression=LinearCombination(
            (
                Term(1, t_inverse_component.atom),
                Term(1, Product((p_plus.atom, p_minus.atom))),
            )
        ),
        components=(t_inverse_component, p_plus, p_minus),
        domain=space,
        codomain=space,
        source="article section 05:770-778",
        evidence=("T1_minus inverse formula",),
    )
    auxiliary = InvertibleAuxiliaryOperator(
        operator=t_minus,
        inverse=t_minus_inverse,
        source="article exact auxiliary invertibility theorem",
        evidence=("section 05:532-782",),
    )
    op_n = _mellin_pdo(
        name="n1",
        expression=radial + frequency,
        radial=radial,
        frequency=frequency,
        assumptions=assumptions,
        space=space,
    )
    phi = OperatorAtom("Phi", kind="space_isomorphism")
    phi_inverse = OperatorAtom("Phi_inverse", kind="space_isomorphism")
    core = MellinCoreOperator(
        atom=OperatorAtom("H1"),
        mellin_operator=op_n,
        phi=phi,
        phi_inverse=phi_inverse,
        source="article sections 05:792-900 and 06:109-132",
        evidence=("Phi N1 Phi_inverse=Op(n1)",),
        alias_of="article operator N1",
    )
    op_r = _mellin_pdo(
        name="r11",
        expression=radial - frequency,
        radial=radial,
        frequency=frequency,
        assumptions=assumptions,
        space=space,
    )
    abstract_regularizer = RegularizerOperator(
        target=core.atom,
        operator=OperatorAtom("R1", kind="formal_regularizer"),
    )
    representation = RegularizerMellinRepresentation(
        regularizer=abstract_regularizer,
        mellin_operator=op_r,
        status=RegularizerRepresentationStatus.CERTIFIED_EXACT,
        hypotheses=assumptions.assumptions,
        source="article section 06:161-168",
        evidence=("R1=Phi_inverse Op(r11) Phi",),
        space=space,
    )
    transported_regularizer = TransportedMellinCore(
        representation=representation,
        phi_inverse=phi_inverse,
        phi=phi,
        source="article exact transport of R1",
        evidence=("section 06:161-168",),
    )
    zeta = (gamma * radial + 1) / (radial + 1)
    z = MultiplicationOperator(
        atom=OperatorAtom("Z1"),
        symbol=zeta,
        radial_variable=radial,
        domain=space,
        codomain=space,
        source="article section 05:339-408",
        evidence=("bounded multiplier zeta1",),
    )
    z_inverse = MultiplicationOperator(
        atom=OperatorAtom("Z1_inverse"),
        symbol=1 / zeta,
        radial_variable=radial,
        domain=space,
        codomain=space,
        source="reciprocal of article zeta1",
        evidence=("bounded reciprocal",),
    )
    multiplier = InvertibleMultiplicationOperator(
        operator=z,
        inverse=z_inverse,
        source="article exact Z1 inverse",
        evidence=("pointwise reciprocal identities",),
        nonvanishing_evidence=(
            "continuous extension and positive nonzero endpoint values",
        ),
    )
    ideal = CompactIdeal(
        label="K(X_delta)",
        space=space,
        bilateral=True,
        source="article section 06:10-12",
        evidence=("compact operators form a bilateral ideal",),
    )
    return _Case(
        assumptions=assumptions,
        space=space,
        other_space=other_space,
        diagonal=diagonal,
        auxiliary=auxiliary,
        core=core,
        transported_regularizer=transported_regularizer,
        multiplier=multiplier,
        ideal=ideal,
        left_exterior=_BoundedOperator(
            OperatorAtom("L21"),
            space,
            space,
        ),
        right_exterior=_BoundedOperator(
            OperatorAtom("R12"),
            space,
            space,
        ),
    )


def _core_certificate(
    case: _Case,
    property: CoreRegularizerProperty = (
        CoreRegularizerProperty.TWO_SIDED_REGULARIZER
    ),
    claim_status: CoreRegularizerClaimStatus = (
        CoreRegularizerClaimStatus.CERTIFIED
    ),
    ideal: CompactIdeal | None = None,
):
    selected_ideal = case.ideal if ideal is None else ideal
    if property is CoreRegularizerProperty.EXACT_INVERSE:
        selected_ideal = None
    if claim_status is CoreRegularizerClaimStatus.UNVERIFIED:
        selected_ideal = None
    return make_mellin_core_regularizer_certificate(
        mellin_operator=case.core,
        mellin_regularizer=case.transported_regularizer,
        property=property,
        claim_status=claim_status,
        compact_ideal=selected_ideal,
        assumptions=case.assumptions,
        source="article section 06:134-171",
        evidence=("KKL 2014 Theorem 5.8 and article transport",),
    )


def _link(
    case: _Case,
    status: DiagonalLinkStatus = DiagonalLinkStatus.EXACT,
    ideal: CompactIdeal | None = None,
) -> DiagonalLinkIdentity:
    left = Product(
        (
            case.multiplier.inverse_atom,
            case.diagonal.atom,
            case.auxiliary.atom,
        )
    )
    right = Product((case.core.atom,))
    selected_ideal = ideal
    if status is DiagonalLinkStatus.EXACT:
        relation = ExactIdentity(
            left,
            right,
            evidence=("article exact factorization",),
            scope=ExactIdentityScope.OPERATORIAL,
        )
    elif status is DiagonalLinkStatus.MODULO_COMPACTS:
        selected_ideal = case.ideal if ideal is None else ideal
        relation = ModCompactEquivalence(
            left,
            right,
            space=case.space,
            compact_ideal=selected_ideal,
            evidence=("hypothetical modulo-compact link",),
        )
    elif status is DiagonalLinkStatus.CONDITIONAL:
        selected_ideal = case.ideal if ideal is None else ideal
        relation = FormalIdentity(
            left,
            right,
            justification=("conditional link hypothesis",),
        )
    else:
        relation = None
    return DiagonalLinkIdentity(
        diagonal_operator=case.diagonal,
        left_interface=case.auxiliary,
        mellin_operator=case.core,
        right_interface=case.multiplier,
        right_interface_inverse=case.multiplier.inverse,
        left_expression=left,
        right_expression=right,
        relation=relation,
        compact_ideal=selected_ideal,
        domain=case.space,
        codomain=case.space,
        assumptions=case.assumptions,
        evidence=("source-backed oriented link",) if relation is not None else (),
        source="article sections 05:429-507 and 05:792-900",
        status=status,
    )


def _certify(
    case: _Case,
    *,
    core=None,
    link=None,
):
    return certify_factored_two_sided_diagonal_regularizer(
        diagonal_operator=case.diagonal,
        left_interface=case.auxiliary,
        mellin_operator=case.core,
        mellin_regularizer=core or _core_certificate(case),
        right_interface=case.multiplier,
        right_interface_inverse=case.multiplier.inverse,
        link_identity=link or _link(case),
        assumptions=case.assumptions,
    )


def test_exact_link_retains_the_article_orientation_and_typed_roles():
    case = _case()
    link = _link(case)

    assert link.status is DiagonalLinkStatus.EXACT
    assert link.left_expression == Product(
        (case.multiplier.inverse_atom, case.diagonal.atom, case.auxiliary.atom)
    )
    assert link.right_expression == Product((case.core.atom,))
    assert isinstance(link.relation, ExactIdentity)
    assert link.compact_ideal is None


def test_mod_compact_and_conditional_links_remain_semantically_distinct():
    case = _case()
    mod = _link(case, DiagonalLinkStatus.MODULO_COMPACTS)
    conditional = _link(case, DiagonalLinkStatus.CONDITIONAL)

    assert isinstance(mod.relation, ModCompactEquivalence)
    assert mod.compact_ideal == case.ideal
    assert isinstance(conditional.relation, FormalIdentity)
    assert conditional.status is DiagonalLinkStatus.CONDITIONAL


def test_unverified_link_carries_no_semantic_relation_or_evidence():
    link = _link(_case(), DiagonalLinkStatus.UNVERIFIED)

    assert link.relation is None
    assert link.evidence == ()


def test_link_rejects_wrong_noncommutative_order():
    case = _case()
    link = _link(case)

    with pytest.raises(DomainOrOrderError, match="link order"):
        replace(
            link,
            left_expression=Product(
                (
                    case.diagonal.atom,
                    case.multiplier.inverse_atom,
                    case.auxiliary.atom,
                )
            ),
        )


def test_link_rejects_duplicated_interface_atoms():
    case = _case()
    duplicate_core = replace(case.core, atom=case.auxiliary.atom)
    right = Product((duplicate_core.atom,))
    relation = ExactIdentity(
        Product(
            (
                case.multiplier.inverse_atom,
                case.diagonal.atom,
                case.auxiliary.atom,
            )
        ),
        right,
        evidence=("invalid duplicate fixture",),
        scope=ExactIdentityScope.OPERATORIAL,
    )

    with pytest.raises(DomainOrOrderError, match="duplicated interfaces"):
        replace(
            _link(case),
            mellin_operator=duplicate_core,
            right_expression=right,
            relation=relation,
        )


def test_link_rejects_an_incompatible_declared_domain():
    case = _case()

    with pytest.raises(DomainOrOrderError, match="different weighted space"):
        replace(_link(case), domain=case.other_space, codomain=case.other_space)


def test_link_rejects_a_semantic_relation_with_another_ideal():
    case = _case()
    other_ideal = CompactIdeal(
        "K_other(X_delta)",
        case.space,
        True,
        "other source",
        ("other theorem",),
    )
    valid = _link(case, DiagonalLinkStatus.MODULO_COMPACTS)
    wrong_relation = ModCompactEquivalence(
        valid.left_expression,
        valid.right_expression,
        space=case.space,
        compact_ideal=other_ideal,
        evidence=("wrong ideal",),
    )

    with pytest.raises(CompactIdealMismatch, match="another ideal"):
        replace(valid, relation=wrong_relation)


def test_link_rejects_missing_evidence():
    with pytest.raises(MissingRegularizerEvidence, match="link requires"):
        replace(_link(_case()), evidence=())


def test_z1_is_a_typed_multiplication_exact_inverse_on_both_sides():
    multiplier = _case().multiplier

    assert isinstance(multiplier.operator, MultiplicationOperator)
    assert isinstance(multiplier.inverse, MultiplicationOperator)
    assert multiplier.operator_times_inverse.left == Product(
        (multiplier.atom, multiplier.inverse_atom)
    )
    assert multiplier.inverse_times_operator.left == Product(
        (multiplier.inverse_atom, multiplier.atom)
    )
    assert multiplier.operator_times_inverse.right == IDENTITY_PRODUCT
    assert multiplier.inverse_times_operator.right == IDENTITY_PRODUCT
    assert "article" in multiplier.source


def test_z1_rejects_a_nonreciprocal_inverse_symbol():
    case = _case()
    wrong = replace(case.multiplier.inverse, symbol=sp.Integer(2))

    with pytest.raises(DomainOrOrderError, match="not exact reciprocals"):
        replace(case.multiplier, inverse=wrong)


def test_z1_rejects_a_symbol_that_can_be_certified_as_zero():
    case = _case()
    zero = replace(case.multiplier.operator, symbol=sp.Integer(0))
    inverse = replace(case.multiplier.inverse, symbol=sp.Integer(1))

    with pytest.raises(DomainOrOrderError, match="certified zero"):
        replace(case.multiplier, operator=zero, inverse=inverse)


def test_z1_rejects_a_positive_half_line_zero_even_with_a_formal_reciprocal():
    case = _case()
    radial = case.multiplier.operator.radial_variable
    vanishing = replace(case.multiplier.operator, symbol=radial - 1)
    formal_reciprocal = replace(
        case.multiplier.inverse,
        symbol=1 / (radial - 1),
    )

    with pytest.raises(DomainOrOrderError, match="positive half-line"):
        replace(
            case.multiplier,
            operator=vanishing,
            inverse=formal_reciprocal,
        )


def test_z1_rejects_missing_nonvanishing_evidence():
    with pytest.raises(MissingRegularizerEvidence, match="nonvanishing"):
        replace(_case().multiplier, nonvanishing_evidence=())


def test_t1_minus_preserves_its_noncommutative_sum_and_exact_inverse():
    auxiliary = _case().auxiliary

    assert isinstance(auxiliary.operator.expression, LinearCombination)
    assert len(auxiliary.operator.expression.terms) == 2
    assert auxiliary.operator.expression.terms[0].product.factors == (
        OperatorAtom("U1_inverse"),
        OperatorAtom("P_plus"),
    )
    assert auxiliary.operator_times_inverse.right == IDENTITY_PRODUCT
    assert auxiliary.inverse_times_operator.right == IDENTITY_PRODUCT


def test_t1_minus_rejects_a_bare_shift_as_its_inverse():
    case = _case()

    with pytest.raises(TypeError, match="bare shift"):
        InvertibleAuxiliaryOperator(
            operator=case.auxiliary.operator,
            inverse=case.auxiliary.inverse.components[0],  # type: ignore[arg-type]
            source="unsupported",
            evidence=("none",),
        )


def test_t1_minus_rejects_unbacked_invertibility():
    with pytest.raises(MissingRegularizerEvidence, match="theorem evidence"):
        replace(_case().auxiliary, evidence=())


@pytest.mark.parametrize(
    ("property", "has_left", "has_right"),
    (
        (CoreRegularizerProperty.LEFT_REGULARIZER, True, False),
        (CoreRegularizerProperty.RIGHT_REGULARIZER, False, True),
        (CoreRegularizerProperty.TWO_SIDED_REGULARIZER, True, True),
        (CoreRegularizerProperty.CALKIN_INVERSE, True, True),
        (CoreRegularizerProperty.EXACT_INVERSE, True, True),
    ),
)
def test_core_regularizer_sidedness_is_never_promoted(
    property: CoreRegularizerProperty,
    has_left: bool,
    has_right: bool,
):
    certificate = _core_certificate(_case(), property)

    assert (certificate.left_relation is not None) is has_left
    assert (certificate.right_relation is not None) is has_right
    if property is CoreRegularizerProperty.EXACT_INVERSE:
        assert isinstance(certificate.left_relation, ExactIdentity)
        assert certificate.compact_ideal is None
    elif certificate.left_relation is not None:
        assert isinstance(certificate.left_relation, ModCompactEquivalence)


def test_core_regularizer_retains_symbol_class_transport_ideal_and_source():
    case = _case()
    certificate = _core_certificate(case)

    assert certificate.mellin_operator.symbol == case.core.symbol
    assert certificate.mellin_operator.symbol_class == "tilde-E"
    assert isinstance(case.core.exact_transport, ExactIdentity)
    assert isinstance(
        certificate.mellin_regularizer.semantic_relation,
        ExactIdentity,
    )
    assert certificate.compact_ideal == case.ideal
    assert "section 06" in certificate.source


def test_core_regularizer_rejects_missing_certification_evidence():
    case = _case()

    with pytest.raises(MissingRegularizerEvidence, match="requires evidence"):
        make_mellin_core_regularizer_certificate(
            mellin_operator=case.core,
            mellin_regularizer=case.transported_regularizer,
            property=CoreRegularizerProperty.TWO_SIDED_REGULARIZER,
            claim_status=CoreRegularizerClaimStatus.CERTIFIED,
            compact_ideal=case.ideal,
            assumptions=case.assumptions,
            source="missing evidence fixture",
            evidence=(),
        )


def test_bilateral_certificate_builds_the_exact_candidate_and_both_products():
    case = _case()
    result = _certify(case)

    assert result.status is (
        FactoredRegularizerStatus.CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER
    )
    assert result.candidate_regularizer.ast_product == Product(
        (
            case.auxiliary.atom,
            case.transported_regularizer.atom,
            case.multiplier.inverse_atom,
        )
    )
    assert result.left_product == Product(
        (case.diagonal.atom, *result.candidate_regularizer.ast_product.factors)
    )
    assert result.right_product == Product(
        (*result.candidate_regularizer.ast_product.factors, case.diagonal.atom)
    )
    assert isinstance(result.left_relation, ModCompactEquivalence)
    assert isinstance(result.right_relation, ModCompactEquivalence)
    assert result.left_relation.compact_ideal == case.ideal
    assert result.right_relation.compact_ideal == case.ideal
    assert result.left_status is (
        ProductCertificationStatus.LEFT_REGULARIZER_CERTIFIED
    )
    assert result.right_status is (
        ProductCertificationStatus.RIGHT_REGULARIZER_CERTIFIED
    )


def test_left_and_right_traces_are_independent_and_use_opposite_core_products():
    result = _certify(_case())

    assert tuple(step.step for step in result.left_reduction_trace) == (
        "left-1-link",
        "left-2-core-regularizer",
        "left-3-exact-cancellation",
    )
    assert tuple(step.step for step in result.right_reduction_trace) == (
        "right-1-link",
        "right-2-core-regularizer",
        "right-3-exact-cancellation",
    )
    assert "H R" in result.left_reduction_trace[1].rule
    assert "R H" in result.right_reduction_trace[1].rule
    assert result.left_reduction_trace is not result.right_reduction_trace


def test_mod_compact_link_and_core_keep_one_common_bilateral_ideal():
    case = _case()
    result = _certify(
        case,
        link=_link(case, DiagonalLinkStatus.MODULO_COMPACTS),
    )

    assert result.status is (
        FactoredRegularizerStatus.CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER
    )
    assert isinstance(result.left_reduction_trace[0].relation, ModCompactEquivalence)
    assert isinstance(result.right_reduction_trace[0].relation, ModCompactEquivalence)
    assert result.compact_ideal == case.ideal


def test_left_only_core_produces_no_right_relation_or_trace():
    case = _case()
    result = _certify(
        case,
        core=_core_certificate(
            case,
            CoreRegularizerProperty.LEFT_REGULARIZER,
        ),
    )

    assert result.status is FactoredRegularizerStatus.CERTIFIED_LEFT_REGULARIZER_ONLY
    assert result.left_relation is not None
    assert result.right_relation is None
    assert result.right_reduction_trace == ()
    assert "regularizer_then_operator_relation" in {
        item.key for item in result.proof_obligations
    }


def test_right_only_core_produces_no_left_relation_or_trace():
    case = _case()
    result = _certify(
        case,
        core=_core_certificate(
            case,
            CoreRegularizerProperty.RIGHT_REGULARIZER,
        ),
    )

    assert result.status is FactoredRegularizerStatus.CERTIFIED_RIGHT_REGULARIZER_ONLY
    assert result.left_relation is None
    assert result.left_reduction_trace == ()
    assert result.right_relation is not None
    assert "operator_then_regularizer_relation" in {
        item.key for item in result.proof_obligations
    }


def test_conditional_link_produces_only_a_conditional_two_sided_result():
    case = _case()
    result = _certify(
        case,
        link=_link(case, DiagonalLinkStatus.CONDITIONAL),
    )

    assert result.status is (
        FactoredRegularizerStatus.CONDITIONAL_TWO_SIDED_FACTORED_REGULARIZER
    )
    assert isinstance(result.left_relation, FormalIdentity)
    assert isinstance(result.right_relation, FormalIdentity)
    assert result.left_status is ProductCertificationStatus.CONDITIONAL
    assert result.right_status is ProductCertificationStatus.CONDITIONAL


def test_conditional_core_is_not_promoted_to_certified():
    case = _case()
    core = _core_certificate(
        case,
        CoreRegularizerProperty.TWO_SIDED_REGULARIZER,
        CoreRegularizerClaimStatus.CONDITIONAL,
    )
    result = _certify(case, core=core)

    assert result.status is (
        FactoredRegularizerStatus.CONDITIONAL_TWO_SIDED_FACTORED_REGULARIZER
    )
    assert isinstance(result.left_relation, FormalIdentity)
    assert isinstance(result.right_relation, FormalIdentity)


def test_unverified_link_blocks_without_reconstructing_it_from_names():
    case = _case()
    result = _certify(
        case,
        link=_link(case, DiagonalLinkStatus.UNVERIFIED),
    )

    assert result.status is FactoredRegularizerStatus.BLOCKED_BY_LINK_IDENTITY
    assert result.left_relation is None
    assert result.right_relation is None
    assert result.left_reduction_trace == ()
    assert result.right_reduction_trace == ()


def test_unverified_core_regularizer_blocks_both_products():
    case = _case()
    core = _core_certificate(
        case,
        CoreRegularizerProperty.TWO_SIDED_REGULARIZER,
        CoreRegularizerClaimStatus.UNVERIFIED,
    )
    result = _certify(case, core=core)

    assert result.status is (
        FactoredRegularizerStatus.BLOCKED_BY_MELLIN_CORE_REGULARIZER
    )
    assert result.left_relation is None
    assert result.right_relation is None


def test_api_blocks_when_link_references_a_different_typed_diagonal():
    case = _case()
    other_diagonal = replace(case.diagonal, atom=OperatorAtom("other_A11"))
    other_left = Product(
        (
            case.multiplier.inverse_atom,
            other_diagonal.atom,
            case.auxiliary.atom,
        )
    )
    other_link = replace(
        _link(case),
        diagonal_operator=other_diagonal,
        left_expression=other_left,
        relation=ExactIdentity(
            other_left,
            Product((case.core.atom,)),
            evidence=("other exact link",),
            scope=ExactIdentityScope.OPERATORIAL,
        ),
    )
    result = _certify(case, link=other_link)

    assert result.status is FactoredRegularizerStatus.BLOCKED_BY_DOMAIN_OR_ORDER


def test_api_raises_compact_ideal_mismatch_before_combining_residues():
    case = _case()
    other_ideal = CompactIdeal(
        "K_other(X_delta)",
        case.space,
        True,
        "different compact ideal",
        ("different evidence",),
    )

    with pytest.raises(CompactIdealMismatch, match="different labels"):
        _certify(
            case,
            link=_link(
                case,
                DiagonalLinkStatus.MODULO_COMPACTS,
                other_ideal,
            ),
        )


def test_nonbilateral_ideal_is_rejected_before_regularizer_construction():
    case = _case()
    one_sided = CompactIdeal(
        "one-sided ideal",
        case.space,
        False,
        "invalid ideal fixture",
        ("no bilateral theorem",),
    )

    with pytest.raises(CompactIdealMismatch, match="not bilateral"):
        _core_certificate(case, ideal=one_sided)


def test_exact_core_inverse_yields_exact_complete_relations_without_an_ideal():
    case = _case()
    result = _certify(
        case,
        core=_core_certificate(
            case,
            CoreRegularizerProperty.EXACT_INVERSE,
        ),
    )

    assert result.status is (
        FactoredRegularizerStatus.CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER
    )
    assert isinstance(result.left_relation, ExactIdentity)
    assert isinstance(result.right_relation, ExactIdentity)
    assert result.compact_ideal is None


def test_main_status_vocabulary_is_exactly_the_nine_admissible_values():
    assert {item.value for item in FactoredRegularizerStatus} == {
        "CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER",
        "CERTIFIED_LEFT_REGULARIZER_ONLY",
        "CERTIFIED_RIGHT_REGULARIZER_ONLY",
        "CONDITIONAL_TWO_SIDED_FACTORED_REGULARIZER",
        "BLOCKED_BY_LINK_IDENTITY",
        "BLOCKED_BY_MELLIN_CORE_REGULARIZER",
        "BLOCKED_BY_DOMAIN_OR_ORDER",
        "BLOCKED_BY_COMPACT_IDEAL_MISMATCH",
        "BLOCKED_BY_SOURCE_GAP",
    }


def test_certificate_retains_original_expressions_latex_evidence_and_obligations():
    result = _certify(_case())

    assert result.original_expression == (result.left_product, result.right_product)
    assert r"\mathcal T_{1,-}\mathcal R_1Z_1^{-1}" in result.latex
    assert result.status.value in result.latex
    assert result.evidence
    assert {
        "full_regularizer_algebra_membership",
        "single_mellin_pdo_representation",
        "cutoff_replacement_mod_compacts",
        "wh_mellin_wh_sandwich_membership",
        "fredholm_algebra_membership",
        "schur_correction_symbol",
    }.issubset({item.key for item in result.proof_obligations})


def test_compact_ideal_propagation_records_lxr_minus_lyr():
    case = _case()
    core = _core_certificate(case)
    assert isinstance(core.left_relation, ModCompactEquivalence)

    propagated = propagate_compact_ideal(
        core.left_relation,
        left_factors=(case.multiplier,),
        right_factors=(case.multiplier.inverse,),
        compact_ideal=case.ideal,
        evidence=("bounded ideal theorem",),
    )

    assert propagated.left == Product(
        (
            case.multiplier.atom,
            case.core.atom,
            core.atom,
            case.multiplier.inverse_atom,
        )
    )
    assert propagated.right == Product(
        (
            case.multiplier.atom,
            IDENTITY_PRODUCT.factors[0],
            case.multiplier.inverse_atom,
        )
    )


def test_compact_ideal_propagation_rejects_unbounded_exterior_factors():
    case = _case()
    core = _core_certificate(case)
    assert isinstance(core.left_relation, ModCompactEquivalence)
    unbounded = replace(case.left_exterior, bounded=False)

    with pytest.raises(DomainOrOrderError, match="explicitly bounded"):
        propagate_compact_ideal(
            core.left_relation,
            left_factors=(unbounded,),
            right_factors=(case.multiplier.inverse,),
            compact_ideal=case.ideal,
            evidence=("invalid",),
        )


def test_compact_ideal_propagation_rejects_another_space():
    case = _case()
    core = _core_certificate(case)
    assert isinstance(core.left_relation, ModCompactEquivalence)
    foreign = _BoundedOperator(
        OperatorAtom("foreign"),
        case.other_space,
        case.other_space,
    )

    with pytest.raises(CompactIdealMismatch, match="different weighted space"):
        propagate_compact_ideal(
            core.left_relation,
            left_factors=(foreign,),
            right_factors=(),
            compact_ideal=case.ideal,
            evidence=("invalid",),
        )


def test_schur_insertion_keeps_the_certified_regularizer_opaque():
    case = _case()
    certificate = _certify(case)
    insertion = insert_factored_regularizer_in_schur(
        left_exterior=case.left_exterior,
        central_regularizer=certificate,
        right_exterior=case.right_exterior,
    )

    assert insertion.schur_insertion is SchurInsertionStatus.STRUCTURALLY_VALID
    assert insertion.expression == Product(
        (
            case.left_exterior.atom,
            certificate.candidate_regularizer.atom,
            case.right_exterior.atom,
        )
    )
    assert insertion.expanded_expression is None
    assert insertion.common_algebra_membership is UnprovedClaimStatus.UNPROVED
    assert (
        insertion.single_mellin_pdo_representation
        is UnprovedClaimStatus.UNPROVED
    )
    assert insertion.schur_symbol is SchurSymbolStatus.NOT_CERTIFIED
    assert insertion.central_regularizer.status is (
        FactoredRegularizerStatus.CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER
    )


@pytest.mark.parametrize(
    "property",
    (
        CoreRegularizerProperty.LEFT_REGULARIZER,
        CoreRegularizerProperty.RIGHT_REGULARIZER,
    ),
)
def test_schur_insertion_rejects_unilateral_regularizers(
    property: CoreRegularizerProperty,
):
    case = _case()
    certificate = _certify(case, core=_core_certificate(case, property))

    with pytest.raises(ValueError, match="certified two-sided"):
        insert_factored_regularizer_in_schur(
            left_exterior=case.left_exterior,
            central_regularizer=certificate,
            right_exterior=case.right_exterior,
        )
