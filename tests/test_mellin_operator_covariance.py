from dataclasses import replace

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    AssumptionContext,
    ComplexDomain,
    DerivationRelationKind,
    ExactIdentityScope,
    MellinConvolutionOperator,
    MellinExpression,
    MellinOperatorError,
    MellinPseudodifferentialOperator,
    MellinSymbol,
    MellinSymbolDependency,
    MellinSymbolDomain,
    MellinVariableRole,
    MultiplicationOperator,
    OperatorAtom,
    Product,
    RadialScalingEvidence,
    RuleCertificationStatus,
    WeightedDilationOperator,
    WeightedLpSpace,
    WienerHopfOperator,
    conjugate_mellin_convolution_by_dilation,
    conjugate_mellin_pdo_by_dilation,
    conjugate_multiplication_by_dilation,
    conjugate_wiener_hopf_by_dilation,
    mellin_frequency,
    output_spatial_variable,
    scale_radial_symbol,
)


def _space(
    *,
    p: sp.Expr = sp.Integer(2),
    delta: sp.Expr = sp.Rational(1, 4),
    label: str = "X_delta",
) -> WeightedLpSpace:
    return WeightedLpSpace(
        label=label,
        p=p,
        weight_exponent=delta,
        assumptions=AssumptionContext(),
        source="article weighted half-line convention",
    )


def _dilation_pair(
    space: WeightedLpSpace,
) -> tuple[
    sp.Symbol,
    AssumptionContext,
    WeightedDilationOperator,
    WeightedDilationOperator,
]:
    gamma = sp.Symbol("gamma", real=True)
    context = AssumptionContext((gamma > 0,))
    dilation = WeightedDilationOperator(
        atom=OperatorAtom("D_gamma", kind="weighted_dilation"),
        gamma=gamma,
        domain=space,
        codomain=space,
        assumptions=context,
        source="normalized positive pullback dilation",
        evidence=("weighted norm calculation",),
    )
    inverse = dilation.inverse(
        OperatorAtom("D_gamma_inverse", kind="weighted_dilation")
    )
    return gamma, context, dilation, inverse


def _mellin_domain(
    lam: sp.Symbol,
    *,
    spatial: sp.Symbol | None = None,
) -> MellinSymbolDomain:
    return MellinSymbolDomain(
        frequency=mellin_frequency(lam),
        complex_domain=ComplexDomain.real_line(lam),
        spatial_variables=(
            (output_spatial_variable(spatial),) if spatial is not None else ()
        ),
        description="typed Mellin covariance test domain",
        evidence=("explicit variable roles",),
    )


def _frequency_symbol(lam: sp.Symbol) -> MellinSymbol:
    domain = _mellin_domain(lam)
    expression = MellinExpression(
        lam + 1,
        domain,
        variables=(domain.frequency,),
        evidence=("frequency symbol evidence",),
        description="radial-independent Mellin multiplier",
    )
    return MellinSymbol(
        expression,
        MellinSymbolDependency.FREQUENCY_ONLY,
        "Mellin convolution symbol",
    )


def _space_frequency_symbol(lam: sp.Symbol, x: sp.Symbol) -> MellinSymbol:
    domain = _mellin_domain(lam, spatial=x)
    expression = MellinExpression(
        x * (lam**2 + 1),
        domain,
        variables=(domain.frequency, domain.spatial_variables[0]),
        evidence=("space-frequency symbol evidence",),
        description="radial Mellin PDO symbol",
    )
    return MellinSymbol(
        expression,
        MellinSymbolDependency.SPACE_FREQUENCY,
        "radial Mellin PDO symbol",
    )


def _mellin_pdo(
    space: WeightedLpSpace,
    lam: sp.Symbol,
    x: sp.Symbol,
    *,
    radial_scaling_stable: bool = True,
) -> MellinPseudodifferentialOperator:
    return MellinPseudodifferentialOperator(
        atom=OperatorAtom("Op_a", kind="mellin_pdo"),
        symbol=_space_frequency_symbol(lam, x),
        domain=space,
        codomain=space,
        symbol_class="tilde-E(R_+,V(R))",
        source="Mellin quantization definition",
        radial_scaling_stable=radial_scaling_stable,
        stability_evidence=("radial-scale stability theorem",)
        if radial_scaling_stable
        else (),
        evidence=("original PDO evidence",),
    )


def _wiener_hopf(
    space: WeightedLpSpace,
    lam: sp.Symbol,
    *,
    frequency_scaling_stable: bool = True,
    localized: bool = False,
    localization_controlled: bool = True,
) -> WienerHopfOperator:
    return WienerHopfOperator(
        atom=OperatorAtom("W_b", kind="wiener_hopf"),
        symbol=sp.exp(-2 * lam),
        frequency_variable=lam,
        domain=space,
        codomain=space,
        symbol_class="M_p(R)",
        source="Fourier multiplier definition",
        frequency_scaling_stable=frequency_scaling_stable,
        stability_evidence=("positive frequency scaling theorem",)
        if frequency_scaling_stable
        else (),
        localized=localized,
        localization_controlled=localization_controlled,
        evidence=("original Wiener--Hopf evidence",),
    )


def test_multiplication_covariance_is_exact_and_preserves_roles_and_evidence() -> None:
    space = _space()
    gamma, _, dilation, inverse = _dilation_pair(space)
    x = sp.Symbol("x", positive=True)
    multiplication = MultiplicationOperator(
        atom=OperatorAtom("M_a", kind="multiplication"),
        symbol=x**2 + 1,
        radial_variable=x,
        domain=space,
        codomain=space,
        source="pointwise multiplier",
        evidence=("original multiplication evidence",),
    )

    trace = conjugate_multiplication_by_dilation(
        dilation,
        multiplication,
        inverse,
    )
    transformed = trace.transformed_operator

    assert isinstance(transformed, MultiplicationOperator)
    assert transformed.symbol == gamma**2 * x**2 + 1
    assert transformed.radial_variable is x
    assert transformed.domain == transformed.codomain == space
    assert transformed.evidence[:-1] == multiplication.evidence
    assert trace.source_product == Product(
        (dilation.atom, multiplication.atom, inverse.atom)
    )
    assert trace.result_product == Product((transformed.atom,))
    assert trace.exact_identity.scope is ExactIdentityScope.OPERATORIAL
    assert trace.exact_identity.evidence.operator_family == "multiplication"
    assert trace.exact_identity.evidence.orientation_check == "a_gamma(x)=a(gamma*x)"
    assert trace.rule.relation_kind is DerivationRelationKind.EXACT_EQUALITY
    assert trace.rule.certification_status is RuleCertificationStatus.CERTIFIED_EXACT


def test_wiener_hopf_covariance_uses_b_of_lambda_over_gamma() -> None:
    space = _space()
    gamma, _, dilation, inverse = _dilation_pair(space)
    lam = sp.Symbol("lambda", real=True)
    operator = _wiener_hopf(space, lam)

    trace = conjugate_wiener_hopf_by_dilation(dilation, operator, inverse)
    transformed = trace.transformed_operator

    assert isinstance(transformed, WienerHopfOperator)
    assert transformed.symbol == sp.exp(-2 * lam / gamma)
    assert transformed.frequency_variable is lam
    assert transformed.symbol_class == operator.symbol_class
    assert transformed.stability_evidence == operator.stability_evidence
    assert transformed.fourier_convention == operator.fourier_convention
    assert transformed.evidence[:-1] == operator.evidence
    assert trace.exact_identity.evidence.orientation_check == (
        "b_gamma(lambda)=b(lambda/gamma)"
    )
    assert trace.rule.certification_status is RuleCertificationStatus.CERTIFIED_EXACT


def test_mellin_convolution_is_exactly_invariant_under_dilation() -> None:
    space = _space()
    _, _, dilation, inverse = _dilation_pair(space)
    lam = sp.Symbol("lambda", real=True)
    operator = MellinConvolutionOperator(
        atom=OperatorAtom("M_c", kind="mellin_convolution"),
        symbol=_frequency_symbol(lam),
        domain=space,
        codomain=space,
        symbol_class="V(R)",
        source="radial-independent Mellin multiplier",
        evidence=("original convolution evidence",),
    )

    trace = conjugate_mellin_convolution_by_dilation(
        dilation,
        operator,
        inverse,
    )

    assert trace.transformed_operator is operator
    assert trace.result_product == Product((operator.atom,))
    assert trace.original_operator.symbol.dependency is (
        MellinSymbolDependency.FREQUENCY_ONLY
    )
    assert trace.original_operator.symbol.evidence == (
        "frequency symbol evidence",
    )
    assert trace.exact_identity.evidence.orientation_check == (
        "c(lambda) is independent of the radial variable"
    )
    assert trace.rule.certification_status is RuleCertificationStatus.CERTIFIED_EXACT


def test_mellin_pdo_covariance_scales_x_forward_and_preserves_metadata() -> None:
    space = _space()
    gamma, context, dilation, inverse = _dilation_pair(space)
    lam = sp.Symbol("lambda", real=True)
    x = sp.Symbol("x", positive=True)
    operator = _mellin_pdo(space, lam, x)

    trace = conjugate_mellin_pdo_by_dilation(
        dilation,
        operator,
        inverse,
        context,
    )
    transformed = trace.transformed_operator

    assert isinstance(transformed, MellinPseudodifferentialOperator)
    assert sp.simplify(
        transformed.symbol.as_expr() - gamma * x * (lam**2 + 1)
    ) == 0
    assert transformed.symbol.frequency.role is MellinVariableRole.FREQUENCY
    assert transformed.symbol.spatial_variables[0].role is (
        MellinVariableRole.OUTPUT_SPATIAL
    )
    gamma_parameter = next(
        item for item in transformed.symbol.parameters if item.symbol == gamma
    )
    assert gamma_parameter.role is MellinVariableRole.PARAMETER
    assert transformed.symbol.assumptions.contains_exact(gamma > 0)
    assert transformed.symbol.evidence[0] == "space-frequency symbol evidence"
    scaling_evidence = tuple(
        item
        for item in transformed.symbol.evidence
        if isinstance(item, RadialScalingEvidence)
    )
    assert len(scaling_evidence) == 1
    assert scaling_evidence[0].radial_variable is x
    assert scaling_evidence[0].original_expression == operator.symbol.as_expr()
    assert scaling_evidence[0].scaled_expression == transformed.symbol.as_expr()
    assert scaling_evidence[0].symbol_class == operator.symbol_class
    assert scaling_evidence[0].stability_evidence == operator.stability_evidence
    assert scaling_evidence[0] in transformed.symbol.domain.evidence
    assert transformed.symbol_class == operator.symbol_class
    assert transformed.stability_evidence == operator.stability_evidence
    assert transformed.evidence[:-1] == operator.evidence
    assert trace.exact_identity.evidence.orientation_check == (
        "a_gamma(x,lambda)=a(gamma*x,lambda)"
    )
    assert trace.rule.certification_status is RuleCertificationStatus.CERTIFIED_EXACT


def test_inverse_orientation_scales_the_mellin_symbol_by_x_over_gamma() -> None:
    space = _space()
    gamma, _, dilation, inverse = _dilation_pair(space)
    lam = sp.Symbol("lambda", real=True)
    x = sp.Symbol("x", positive=True)
    operator = _mellin_pdo(space, lam, x)

    trace = conjugate_mellin_pdo_by_dilation(
        inverse,
        operator,
        dilation,
        inverse.assumptions,
    )
    transformed = trace.transformed_operator

    assert isinstance(transformed, MellinPseudodifferentialOperator)
    assert sp.simplify(
        transformed.symbol.as_expr() - x * (lam**2 + 1) / gamma
    ) == 0
    assert sp.simplify(trace.dilation.gamma - 1 / gamma) == 0
    assert trace.exact_identity.evidence.orientation_check == (
        "a_gamma(x,lambda)=a(gamma*x,lambda)"
    )


def test_public_radial_scaling_requires_symbol_class_stability_evidence() -> None:
    gamma = sp.Symbol("gamma", positive=True)
    lam = sp.Symbol("lambda", real=True)
    x = sp.Symbol("x", positive=True)
    symbol = _space_frequency_symbol(lam, x)

    with pytest.raises(MellinOperatorError, match="stability evidence"):
        scale_radial_symbol(
            symbol,
            gamma,
            AssumptionContext(),
            symbol_class="tilde-E",
            stability_evidence=(),
        )


def test_mellin_pdo_rejects_unhashable_stability_evidence_at_construction() -> None:
    space = _space()
    lam = sp.Symbol("lambda", real=True)
    x = sp.Symbol("x", positive=True)

    with pytest.raises(TypeError, match="stability_evidence items must be hashable"):
        replace(
            _mellin_pdo(space, lam, x),
            stability_evidence=({"proof": "external"},),
        )


def test_covariance_rejects_incompatible_fourier_or_mellin_conventions() -> None:
    space = _space()
    _, context, dilation, inverse = _dilation_pair(space)
    lam = sp.Symbol("lambda", real=True)
    x = sp.Symbol("x", positive=True)
    wrong_fourier = replace(
        _wiener_hopf(space, lam),
        fourier_convention="incompatible positive-exponent Fourier transform",
    )
    wrong_quantization = replace(
        _mellin_pdo(space, lam, x),
        quantization="incompatible right Mellin quantization",
    )

    with pytest.raises(MellinOperatorError, match="repository Fourier convention"):
        conjugate_wiener_hopf_by_dilation(dilation, wrong_fourier, inverse)
    with pytest.raises(MellinOperatorError, match="repository Mellin quantization"):
        conjugate_mellin_pdo_by_dilation(
            dilation,
            wrong_quantization,
            inverse,
            context,
        )


def test_covariance_rejects_explicitly_conflicting_dilation_contexts() -> None:
    space = _space()
    gamma, _, dilation, _ = _dilation_pair(space)
    conflicting = AssumptionContext((gamma < 0, 1 / gamma > 0))
    inverse = WeightedDilationOperator(
        atom=OperatorAtom("D_conflicting_inverse", kind="weighted_dilation"),
        gamma=1 / gamma,
        domain=space,
        codomain=space,
        assumptions=conflicting,
        source="deliberately contradictory inverse metadata",
    )
    x = sp.Symbol("x", positive=True)
    operator = MultiplicationOperator(
        atom=OperatorAtom("M_a", kind="multiplication"),
        symbol=x + 1,
        radial_variable=x,
        domain=space,
        codomain=space,
        source="bounded multiplier",
    )

    with pytest.raises(MellinOperatorError, match="explicitly incompatible"):
        conjugate_multiplication_by_dilation(dilation, operator, inverse)


def test_exact_covariance_retains_the_actual_dilation_hypotheses() -> None:
    space = _space()
    gamma, _, dilation, inverse = _dilation_pair(space)
    x = sp.Symbol("x", positive=True)
    operator = MultiplicationOperator(
        atom=OperatorAtom("M_a", kind="multiplication"),
        symbol=x + 1,
        radial_variable=x,
        domain=space,
        codomain=space,
        source="bounded multiplier",
    )

    trace = conjugate_multiplication_by_dilation(dilation, operator, inverse)

    assert (gamma > 0) in trace.exact_identity.hypotheses
    assert (1 / gamma > 0) in trace.exact_identity.hypotheses
    assert trace.assumptions.contains_exact(gamma > 0)
    assert trace.assumptions.contains_exact(1 / gamma > 0)


def test_constant_wiener_hopf_symbol_is_valid_and_scale_invariant() -> None:
    space = _space()
    _, _, dilation, inverse = _dilation_pair(space)
    lam = sp.Symbol("lambda", real=True)
    operator = WienerHopfOperator(
        atom=OperatorAtom("W_constant", kind="wiener_hopf"),
        symbol=sp.Integer(3),
        frequency_variable=lam,
        domain=space,
        codomain=space,
        symbol_class="constant Fourier multipliers",
        source="constant multiplier",
        frequency_scaling_stable=True,
        stability_evidence=("constants are scale invariant",),
    )

    trace = conjugate_wiener_hopf_by_dilation(dilation, operator, inverse)

    assert trace.transformed_operator.symbol == 3


def test_mellin_pdo_covariance_rejects_a_nonstable_symbol_class() -> None:
    space = _space()
    _, context, dilation, inverse = _dilation_pair(space)
    lam = sp.Symbol("lambda", real=True)
    x = sp.Symbol("x", positive=True)
    operator = _mellin_pdo(
        space,
        lam,
        x,
        radial_scaling_stable=False,
    )

    with pytest.raises(MellinOperatorError, match="not radially scale-stable"):
        conjugate_mellin_pdo_by_dilation(
            dilation,
            operator,
            inverse,
            context,
        )


def test_wiener_hopf_covariance_rejects_a_nonstable_symbol_class() -> None:
    space = _space()
    _, _, dilation, inverse = _dilation_pair(space)
    lam = sp.Symbol("lambda", real=True)
    operator = _wiener_hopf(
        space,
        lam,
        frequency_scaling_stable=False,
    )

    with pytest.raises(MellinOperatorError, match="not scale-stable"):
        conjugate_wiener_hopf_by_dilation(dilation, operator, inverse)


@pytest.mark.parametrize("localization_controlled", (False, True))
def test_localized_wiener_hopf_is_rejected_until_cutoff_rescaling_is_typed(
    localization_controlled: bool,
) -> None:
    space = _space()
    _, _, dilation, inverse = _dilation_pair(space)
    lam = sp.Symbol("lambda", real=True)
    operator = _wiener_hopf(
        space,
        lam,
        localized=True,
        localization_controlled=localization_controlled,
    )

    with pytest.raises(MellinOperatorError, match="cutoff rescaling"):
        conjugate_wiener_hopf_by_dilation(dilation, operator, inverse)


def test_all_exact_covariance_helpers_reject_unbounded_operators() -> None:
    space = _space()
    _, context, dilation, inverse = _dilation_pair(space)
    lam = sp.Symbol("lambda", real=True)
    x = sp.Symbol("x", positive=True)
    multiplication = MultiplicationOperator(
        atom=OperatorAtom("M_unbounded", kind="multiplication"),
        symbol=x + 1,
        radial_variable=x,
        domain=space,
        codomain=space,
        source="negative boundedness test",
        bounded=False,
    )
    wiener_hopf = replace(_wiener_hopf(space, lam), bounded=False)
    convolution = MellinConvolutionOperator(
        atom=OperatorAtom("M_c_unbounded", kind="mellin_convolution"),
        symbol=_frequency_symbol(lam),
        domain=space,
        codomain=space,
        symbol_class="V(R)",
        source="negative boundedness test",
        bounded=False,
    )
    pdo = replace(_mellin_pdo(space, lam, x), bounded=False)
    cases = (
        (conjugate_multiplication_by_dilation, (dilation, multiplication, inverse)),
        (conjugate_wiener_hopf_by_dilation, (dilation, wiener_hopf, inverse)),
        (conjugate_mellin_convolution_by_dilation, (dilation, convolution, inverse)),
        (conjugate_mellin_pdo_by_dilation, (dilation, pdo, inverse, context)),
    )

    for covariance, arguments in cases:
        with pytest.raises(MellinOperatorError, match="bounded operator"):
            covariance(*arguments)


@pytest.mark.parametrize(
    "operator_space",
    (
        _space(delta=sp.Rational(1, 5), label="different_weight"),
        _space(p=sp.Integer(3), label="different_p"),
    ),
)
def test_covariance_rejects_operator_space_or_weight_mismatch(
    operator_space: WeightedLpSpace,
) -> None:
    dilation_space = _space()
    _, _, dilation, inverse = _dilation_pair(dilation_space)
    x = sp.Symbol("x", positive=True)
    operator = MultiplicationOperator(
        atom=OperatorAtom("M_mismatch", kind="multiplication"),
        symbol=x + 1,
        radial_variable=x,
        domain=operator_space,
        codomain=operator_space,
        source="mismatched operator",
    )

    with pytest.raises(MellinOperatorError, match="do not match"):
        conjugate_multiplication_by_dilation(dilation, operator, inverse)


def test_covariance_rejects_inverse_with_a_different_normalization() -> None:
    space = _space()
    gamma, context, dilation, _ = _dilation_pair(space)
    wrong_inverse = WeightedDilationOperator(
        atom=OperatorAtom("wrong_inverse", kind="weighted_dilation"),
        gamma=1 / gamma,
        domain=space,
        codomain=space,
        assumptions=AssumptionContext((*context.assumptions, 1 / gamma > 0)),
        source="wrong normalization",
        normalization_power=space.normalization_exponent + 1,
    )
    x = sp.Symbol("x", positive=True)
    operator = MultiplicationOperator(
        atom=OperatorAtom("M_a", kind="multiplication"),
        symbol=x + 1,
        radial_variable=x,
        domain=space,
        codomain=space,
        source="pointwise multiplier",
    )

    with pytest.raises(MellinOperatorError, match="different normalizations"):
        conjugate_multiplication_by_dilation(dilation, operator, wrong_inverse)
