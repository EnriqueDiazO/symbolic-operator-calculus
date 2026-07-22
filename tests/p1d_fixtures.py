"""Reusable, source-shaped P1-D test data.

The objects use simple scalar symbols, but every operator interface and
semantic relation has the same typed shape as the article construction.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from symbolic_operator_calculus import (
    AssumptionContext,
    Branch,
    BranchDifferenceOperator,
    CoefficientMultiplicationOperator,
    CompactIdeal,
    ComplexDomain,
    CoreRegularizerClaimStatus,
    CoreRegularizerProperty,
    DiagonalLinkIdentity,
    DiagonalLinkStatus,
    DiagonalOperator,
    ExactIdentity,
    ExactIdentityScope,
    InvertibleAuxiliaryOperator,
    InvertibleMultiplicationOperator,
    LinearCombination,
    MellinCoreOperator,
    MellinExpression,
    MellinPseudodifferentialOperator,
    MellinSymbol,
    MellinSymbolDependency,
    MellinSymbolDomain,
    MultiplicationOperator,
    OperatorAtom,
    Product,
    RegularizerMellinRepresentation,
    RegularizerOperator,
    RegularizerRepresentationStatus,
    SchurExteriorFactor,
    SchurSourceBlock,
    SingularSet,
    SourceBlockSide,
    SpatialCutoffOperator,
    Term,
    TransportedMellinCore,
    TransportedShiftOperator,
    WeightedDilationOperator,
    WeightedLpSpace,
    WienerHopfOperator,
    certify_factored_two_sided_diagonal_regularizer,
    make_mellin_core_regularizer_certificate,
    mellin_frequency,
    output_spatial_variable,
)
from symbolic_operator_calculus.schur_full_word import (
    AuxiliaryOperator,
    LocalizedOperator,
)


@dataclass(frozen=True)
class P1DCase:
    assumptions: AssumptionContext
    space: WeightedLpSpace
    other_space: WeightedLpSpace
    ideal: CompactIdeal
    gamma1: sp.Symbol
    gamma2: sp.Symbol
    radial: sp.Symbol
    frequency: sp.Symbol
    central: object
    u1: WeightedDilationOperator
    u1_inverse: WeightedDilationOperator
    u2: WeightedDilationOperator
    left: SchurExteriorFactor
    right: SchurExteriorFactor
    left_source: SchurSourceBlock
    right_source: SchurSourceBlock


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
    expression_model = MellinExpression(
        expression=expression,
        domain=domain,
        variables=(typed_frequency, typed_radial),
        evidence=(f"formula for {name}",),
        description=name,
    )
    return MellinPseudodifferentialOperator(
        atom=OperatorAtom(f"Op_{name}"),
        symbol=MellinSymbol(
            expression_model,
            MellinSymbolDependency.SPACE_FREQUENCY,
            description=name,
        ),
        domain=space,
        codomain=space,
        symbol_class="tilde-E",
        source=f"source-shaped {name} Mellin PDO",
        radial_scaling_stable=True,
        stability_evidence=("radial scaling theorem",),
        evidence=("bounded Mellin PDO",),
    )


def _central_regularizer(
    *,
    assumptions: AssumptionContext,
    space: WeightedLpSpace,
    ideal: CompactIdeal,
    radial: sp.Symbol,
    frequency: sp.Symbol,
    u1_inverse: WeightedDilationOperator,
):
    p_plus = MultiplicationOperator(
        atom=OperatorAtom("P_plus"),
        symbol=sp.Integer(1),
        radial_variable=radial,
        domain=space,
        codomain=space,
        source="typed Cauchy projection stand-in",
        evidence=("P-plus source",),
    )
    p_minus = MultiplicationOperator(
        atom=OperatorAtom("P_minus"),
        symbol=sp.Integer(1),
        radial_variable=radial,
        domain=space,
        codomain=space,
        source="typed Cauchy projection stand-in",
        evidence=("P-minus source",),
    )
    t_minus_raw = AuxiliaryOperator(
        atom=OperatorAtom("T1_minus"),
        expression=LinearCombination(
            (
                Term(1, Product((u1_inverse.atom, p_plus.atom))),
                Term(1, p_minus.atom),
            )
        ),
        components=(u1_inverse, p_plus, p_minus),
        domain=space,
        codomain=space,
        source="article section 05:515-529",
        evidence=("T1_minus=U1_inverse P_plus+P_minus",),
    )
    inverse_component = MultiplicationOperator(
        atom=OperatorAtom("Op_theta1_inverse"),
        symbol=sp.Integer(1),
        radial_variable=radial,
        domain=space,
        codomain=space,
        source="auxiliary inverse component",
        evidence=("article auxiliary inverse theorem",),
    )
    t_inverse_raw = AuxiliaryOperator(
        atom=OperatorAtom("T1_minus_inverse"),
        expression=LinearCombination(
            (
                Term(1, inverse_component.atom),
                Term(1, Product((p_plus.atom, p_minus.atom))),
            )
        ),
        components=(inverse_component, p_plus, p_minus),
        domain=space,
        codomain=space,
        source="article section 05:770-778",
        evidence=("T1-minus inverse source",),
    )
    auxiliary = InvertibleAuxiliaryOperator(
        operator=t_minus_raw,
        inverse=t_inverse_raw,
        source="article exact auxiliary invertibility",
        evidence=("section 05:532-782",),
    )
    diagonal = DiagonalOperator(
        atom=OperatorAtom("A11"),
        domain=space,
        codomain=space,
        source="article section 06:27-42",
        evidence=("A11 exact definition",),
    )
    op_n = _mellin_pdo(
        name="n1",
        expression=radial + frequency,
        radial=radial,
        frequency=frequency,
        assumptions=assumptions,
        space=space,
    )
    phi = OperatorAtom("Phi_delta", kind="space_isomorphism")
    phi_inverse = OperatorAtom("Phi_delta_inverse", kind="space_isomorphism")
    core = MellinCoreOperator(
        atom=OperatorAtom("N1"),
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
    regularizer = RegularizerOperator(
        target=core.atom,
        operator=OperatorAtom("R1", kind="formal_regularizer"),
    )
    representation = RegularizerMellinRepresentation(
        regularizer=regularizer,
        mellin_operator=op_r,
        status=RegularizerRepresentationStatus.CERTIFIED_EXACT,
        hypotheses=assumptions.assumptions,
        source="article section 06:161-168",
        evidence=("R1=Phi_inverse Op(r11) Phi",),
        space=space,
    )
    transported = TransportedMellinCore(
        representation=representation,
        phi_inverse=phi_inverse,
        phi=phi,
        source="article exact transport of R1",
        evidence=("section 06:161-168",),
    )
    core_regularizer = make_mellin_core_regularizer_certificate(
        mellin_operator=core,
        mellin_regularizer=transported,
        property=CoreRegularizerProperty.TWO_SIDED_REGULARIZER,
        claim_status=CoreRegularizerClaimStatus.CERTIFIED,
        compact_ideal=ideal,
        assumptions=assumptions,
        source="article section 06:134-171",
        evidence=("KKL 2014 regularizer theorem and article transport",),
    )
    zeta = (sp.Symbol("gamma1", positive=True) * radial + 1) / (radial + 1)
    z1 = MultiplicationOperator(
        atom=OperatorAtom("Z1"),
        symbol=zeta,
        radial_variable=radial,
        domain=space,
        codomain=space,
        source="article section 05:339-408",
        evidence=("bounded nonvanishing zeta1",),
    )
    z1_inverse = MultiplicationOperator(
        atom=OperatorAtom("Z1_inverse"),
        symbol=1 / zeta,
        radial_variable=radial,
        domain=space,
        codomain=space,
        source="reciprocal of article zeta1",
        evidence=("bounded reciprocal",),
    )
    multiplier = InvertibleMultiplicationOperator(
        operator=z1,
        inverse=z1_inverse,
        source="article exact Z1 inverse",
        evidence=("pointwise reciprocal identity",),
        nonvanishing_evidence=("positive nonzero endpoint values",),
    )
    link_left = Product((multiplier.inverse_atom, diagonal.atom, auxiliary.atom))
    link_right = Product((core.atom,))
    link = DiagonalLinkIdentity(
        diagonal_operator=diagonal,
        left_interface=auxiliary,
        mellin_operator=core,
        right_interface=multiplier,
        right_interface_inverse=multiplier.inverse,
        left_expression=link_left,
        right_expression=link_right,
        relation=ExactIdentity(
            link_left,
            link_right,
            evidence=("article exact link",),
            scope=ExactIdentityScope.OPERATORIAL,
        ),
        compact_ideal=None,
        domain=space,
        codomain=space,
        assumptions=assumptions,
        evidence=("source-backed oriented link",),
        source="article sections 05:429-507 and 05:792-900",
        status=DiagonalLinkStatus.EXACT,
    )
    return certify_factored_two_sided_diagonal_regularizer(
        diagonal_operator=diagonal,
        left_interface=auxiliary,
        mellin_operator=core,
        mellin_regularizer=core_regularizer,
        right_interface=multiplier,
        right_interface_inverse=multiplier.inverse,
        link_identity=link,
        assumptions=assumptions,
    )


def _exterior(
    *,
    index: int,
    branch_from: int,
    branch_to: int,
    dilation: WeightedDilationOperator,
    assumptions: AssumptionContext,
    space: WeightedLpSpace,
    radial: sp.Symbol,
    frequency: sp.Symbol,
) -> SchurExteriorFactor:
    transport_multiplier = CoefficientMultiplicationOperator(
        MultiplicationOperator(
            atom=OperatorAtom(f"M_rho{index}"),
            symbol=(dilation.gamma * radial + 1) / (radial + 1),
            radial_variable=radial,
            domain=space,
            codomain=space,
            source=f"article rho_{index}",
            evidence=("bounded transport coefficient",),
        ),
        coefficient_role=f"transport coefficient rho_{index}",
    )
    shift = TransportedShiftOperator(
        atom=OperatorAtom(f"Vtilde_alpha{index}"),
        coefficient=transport_multiplier,
        dilation=dilation,
        source=f"article transported shift on branch {index}",
        evidence=("Vtilde=M_rho U",),
    )
    coefficient = CoefficientMultiplicationOperator(
        MultiplicationOperator(
            atom=OperatorAtom(f"G{index}"),
            symbol=sp.Integer(index),
            radial_variable=radial,
            domain=space,
            codomain=space,
            source=f"article coefficient G_{index}",
            evidence=("bounded coefficient",),
        ),
        coefficient_role=f"diagonal coefficient G_{index}",
    )
    difference = BranchDifferenceOperator(
        atom=OperatorAtom(f"D{index}"),
        transported_shift=shift,
        coefficient=coefficient,
        source=f"article Vtilde_alpha{index}-G{index}",
        assumptions=assumptions,
        evidence=("normalized exterior difference",),
    )
    polarity = "minus" if index == 2 else "plus"
    wh = WienerHopfOperator(
        atom=OperatorAtom(f"W_core_{polarity}"),
        symbol=sp.exp(-sp.Abs(frequency)) + index * frequency,
        frequency_variable=frequency,
        domain=space,
        codomain=space,
        symbol_class="article exponential half-line symbol",
        source=f"article b_{polarity}",
        frequency_scaling_stable=True,
        stability_evidence=("half-line indicator is scale invariant",),
        localized=False,
        localization_controlled=True,
        evidence=("exact Wiener--Hopf symbol",),
    )
    chi = sp.Function("chi_infty")
    left_cutoff = SpatialCutoffOperator(
        atom=OperatorAtom(f"M_chi_{polarity}_left"),
        cutoff_symbol=chi(radial),
        radial_variable=radial,
        domain=space,
        codomain=space,
        support_region="transition on [1,2]",
        equals_one_near="[2,infinity)",
        equals_zero_away_from="[0,1]",
        coordinate_system="half-line radial coordinate",
        radial_scale=sp.Integer(1),
        source="article section 03:747-764",
        evidence=("chi_infty definition",),
    )
    right_cutoff = SpatialCutoffOperator(
        atom=OperatorAtom(f"M_chi_{polarity}_right"),
        cutoff_symbol=chi(radial),
        radial_variable=radial,
        domain=space,
        codomain=space,
        support_region="transition on [1,2]",
        equals_one_near="[2,infinity)",
        equals_zero_away_from="[0,1]",
        coordinate_system="half-line radial coordinate",
        radial_scale=sp.Integer(1),
        source="article section 03:747-764",
        evidence=("chi_infty definition",),
    )
    localized = LocalizedOperator(
        atom=OperatorAtom(f"W{21 if index == 2 else 12}_{polarity}"),
        left_cutoff=left_cutoff,
        core=wh,
        right_cutoff=right_cutoff,
        source="article normalized doubly localized Wiener--Hopf block",
        evidence=("section 06:60-95",),
    )
    return SchurExteriorFactor(
        atom=OperatorAtom("L21" if index == 2 else "R12"),
        difference=difference,
        localized_wiener_hopf=localized,
        branch_from=Branch(branch_from),
        branch_to=Branch(branch_to),
        source=f"article exterior branch {branch_from}->{branch_to}",
        assumptions=assumptions,
        evidence=("section 06 normalized block",),
    )


def build_p1d_case() -> P1DCase:
    gamma1 = sp.Symbol("gamma1", positive=True)
    gamma2 = sp.Symbol("gamma2", positive=True)
    radial = sp.Symbol("x", positive=True)
    frequency = sp.Symbol("lambda", real=True)
    kappa = sp.Symbol("kappa", positive=True)
    assumptions = AssumptionContext((gamma1 > 0, gamma2 > 0, kappa < 1))
    space = WeightedLpSpace(
        label="X=L^p(R_+,w_delta)",
        p=sp.Integer(2),
        weight_exponent=sp.Integer(0),
        assumptions=assumptions,
        source="article scalar weighted space",
    )
    other_space = WeightedLpSpace(
        label="Y",
        p=sp.Integer(2),
        weight_exponent=sp.Integer(0),
        assumptions=assumptions,
        source="different test space",
    )
    ideal = CompactIdeal(
        label="K(X)",
        space=space,
        bilateral=True,
        source="article section 06:10-12",
        evidence=("compact operators form a bilateral ideal",),
    )
    u1 = WeightedDilationOperator(
        atom=OperatorAtom("U1"),
        gamma=gamma1,
        domain=space,
        codomain=space,
        assumptions=assumptions,
        source="normalized branch-one dilation",
        evidence=("article normalized dilation",),
    )
    u1_inverse = u1.inverse(OperatorAtom("U1_inverse"))
    u2 = WeightedDilationOperator(
        atom=OperatorAtom("U2"),
        gamma=gamma2,
        domain=space,
        codomain=space,
        assumptions=assumptions,
        source="normalized branch-two dilation",
        evidence=("article normalized dilation",),
    )
    central = _central_regularizer(
        assumptions=assumptions,
        space=space,
        ideal=ideal,
        radial=radial,
        frequency=frequency,
        u1_inverse=u1_inverse,
    )
    left = _exterior(
        index=2,
        branch_from=1,
        branch_to=2,
        dilation=u2,
        assumptions=assumptions,
        space=space,
        radial=radial,
        frequency=frequency,
    )
    right = _exterior(
        index=1,
        branch_from=2,
        branch_to=1,
        dilation=u1,
        assumptions=assumptions,
        space=space,
        radial=radial,
        frequency=frequency,
    )
    left_source = SchurSourceBlock(
        atom=OperatorAtom("B21_minus"),
        side=SourceBlockSide.LEFT,
        exterior=left,
        interface=central.candidate_regularizer.left_interface,
        source="article section 06:382-390",
        evidence=("B21_minus=L21 T1_minus",),
    )
    right_source = SchurSourceBlock(
        atom=OperatorAtom("B12_plus"),
        side=SourceBlockSide.RIGHT,
        exterior=right,
        interface=central.candidate_regularizer.right_interface_inverse,
        source="article section 06:391-399",
        evidence=("B12_plus=Z1_inverse R12",),
    )
    return P1DCase(
        assumptions=assumptions,
        space=space,
        other_space=other_space,
        ideal=ideal,
        gamma1=gamma1,
        gamma2=gamma2,
        radial=radial,
        frequency=frequency,
        central=central,
        u1=u1,
        u1_inverse=u1_inverse,
        u2=u2,
        left=left,
        right=right,
        left_source=left_source,
        right_source=right_source,
    )

