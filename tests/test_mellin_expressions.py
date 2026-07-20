import pytest
import sympy as sp

from symbolic_operator_calculus import (
    AssumptionContext,
    ComplexDomain,
    ConditionalVerificationStatus,
    MellinBranchConditionRequiredError,
    MellinExpression,
    MellinExpressionError,
    MellinSingularMetadataRequiredError,
    MellinSymbolDomain,
    SingularSet,
    UndeclaredMellinVariableError,
    detect_singularities,
    mellin_frequency,
    mellin_parameter,
    output_spatial_variable,
)


def _domain_for(
    expression: sp.Expr,
    *,
    lam: sp.Symbol,
    parameters: tuple[sp.Symbol, ...] = (),
    spatial: tuple[sp.Symbol, ...] = (),
    context: AssumptionContext = AssumptionContext(),
    exclusions: tuple[sp.Basic, ...] = (),
) -> MellinSymbolDomain:
    frequency = mellin_frequency(lam)
    parameter_variables = tuple(mellin_parameter(item) for item in parameters)
    spatial_variables = tuple(output_spatial_variable(item) for item in spatial)
    complex_domain = ComplexDomain.real_line(
        lam,
        assumption_context=context,
        exclusions=exclusions,
    )
    singular_set = detect_singularities(
        expression,
        lam,
        assumptions=context,
        domain=complex_domain,
    )
    return MellinSymbolDomain(
        frequency=frequency,
        complex_domain=complex_domain,
        spatial_variables=spatial_variables,
        parameters=parameter_variables,
        assumption_context=context,
        singular_set=singular_set,
    )


def test_expression_requires_existing_sympy_expression_and_explicit_domain() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain_for(lam, lam=lam)

    with pytest.raises(TypeError, match="SymPy scalar expression"):
        MellinExpression("lambda", domain)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="MellinSymbolDomain"):
        MellinExpression(lam, object())  # type: ignore[arg-type]


def test_every_free_symbol_must_be_declared() -> None:
    lam, unknown = sp.symbols("lambda unknown")
    domain = _domain_for(lam + unknown, lam=lam)

    with pytest.raises(UndeclaredMellinVariableError, match="unknown"):
        MellinExpression(
            lam + unknown,
            domain,
            variables=(domain.frequency,),
        )


def test_omitted_parameter_is_rejected() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    domain = _domain_for(lam + kappa, lam=lam, parameters=(kappa,))

    with pytest.raises(UndeclaredMellinVariableError, match="kappa"):
        MellinExpression(
            lam + kappa,
            domain,
            variables=(domain.frequency,),
        )


def test_parameter_cannot_be_declared_as_an_expression_variable() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    domain = _domain_for(lam + kappa, lam=lam, parameters=(kappa,))

    with pytest.raises(MellinExpressionError, match="parameters"):
        MellinExpression(
            lam + kappa,
            domain,
            variables=(domain.frequency, domain.parameters[0]),
        )


def test_frequency_and_parameter_declarations_are_preserved() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    context = AssumptionContext((kappa > 0, kappa < 1))
    domain = _domain_for(
        lam + sp.I * kappa,
        lam=lam,
        parameters=(kappa,),
        context=context,
    )
    expression = MellinExpression(
        lam + sp.I * kappa,
        domain,
        variables=(domain.frequency,),
        parameters=domain.parameters,
    )

    assert expression.variables == (domain.frequency,)
    assert expression.parameters == domain.parameters
    assert expression.domain.assumption_context.contains_all(context)


def test_limited_detected_singularities_must_be_declared() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    context = AssumptionContext((kappa > 0, kappa < 1))
    expression = sp.coth(sp.pi * (lam + sp.I * kappa))
    domain = MellinSymbolDomain(
        frequency=mellin_frequency(lam),
        complex_domain=ComplexDomain.real_line(
            lam, assumption_context=context
        ),
        parameters=(mellin_parameter(kappa),),
        assumption_context=context,
        singular_set=SingularSet(assumption_context=context),
    )

    with pytest.raises(MellinSingularMetadataRequiredError, match="pole"):
        MellinExpression(
            expression,
            domain,
            variables=(domain.frequency,),
            parameters=domain.parameters,
        )


def test_branch_sensitive_parameter_requires_explicit_positive_base() -> None:
    lam, gamma = sp.symbols("lambda gamma")
    expression = sp.Pow(gamma, sp.I * lam, evaluate=False)
    domain = _domain_for(expression, lam=lam, parameters=(gamma,))

    with pytest.raises(MellinBranchConditionRequiredError, match="gamma > 0"):
        MellinExpression(
            expression,
            domain,
            variables=(domain.frequency,),
            parameters=domain.parameters,
        )


def test_explicit_negative_base_is_rejected() -> None:
    lam, gamma = sp.symbols("lambda gamma")
    expression = sp.Pow(gamma, sp.I * lam, evaluate=False)
    context = AssumptionContext((gamma < 0,))
    domain = _domain_for(
        expression,
        lam=lam,
        parameters=(gamma,),
        context=context,
    )

    with pytest.raises(MellinBranchConditionRequiredError, match="nonpositive"):
        MellinExpression(
            expression,
            domain,
            variables=(domain.frequency,),
            parameters=domain.parameters,
        )


def test_positive_branch_condition_and_metadata_are_retained() -> None:
    lam, gamma = sp.symbols("lambda gamma")
    expression = sp.Pow(gamma, sp.I * lam, evaluate=False)
    context = AssumptionContext((gamma > 0,))
    domain = _domain_for(
        expression,
        lam=lam,
        parameters=(gamma,),
        context=context,
    )
    formal = MellinExpression(
        expression,
        domain,
        variables=(domain.frequency,),
        parameters=domain.parameters,
    )

    assert formal.domain.assumption_context.contains_exact(gamma > 0)
    assert formal.domain.singular_set.avoidances


def test_evidence_is_not_automatic_internal_verification() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain_for(lam, lam=lam)
    declared = MellinExpression(
        lam,
        domain,
        variables=(domain.frequency,),
        evidence=("external lemma",),
    )

    assert declared.verification_status is ConditionalVerificationStatus.DECLARED
    assert declared.evidence == ("external lemma",)

    with pytest.raises(MellinExpressionError, match="requires"):
        MellinExpression(
            lam,
            domain,
            variables=(domain.frequency,),
            verification_status=ConditionalVerificationStatus.EVIDENCE_SUPPLIED,
        )


def test_standalone_expression_cannot_claim_internal_identity_check() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain_for(lam, lam=lam)

    with pytest.raises(MellinExpressionError, match="standalone"):
        MellinExpression(
            lam,
            domain,
            variables=(domain.frequency,),
            verification_status=(
                ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
            ),
            evidence=("not a check record",),
        )


def test_binary_operations_intersect_domains_and_use_weaker_status() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    first_context = AssumptionContext((kappa > 0,))
    second_context = AssumptionContext((kappa < 1,))
    first_domain = _domain_for(
        lam,
        lam=lam,
        parameters=(kappa,),
        context=first_context,
        exclusions=(0,),
    )
    second_domain = _domain_for(
        lam + kappa,
        lam=lam,
        parameters=(kappa,),
        context=second_context,
        exclusions=(sp.I,),
    )
    first = MellinExpression(
        lam,
        first_domain,
        variables=(first_domain.frequency,),
        verification_status=ConditionalVerificationStatus.EVIDENCE_SUPPLIED,
        evidence=("paper",),
    )
    second = MellinExpression(
        lam + kappa,
        second_domain,
        variables=(second_domain.frequency,),
        parameters=second_domain.parameters,
    )

    for result in (first + second, first * second):
        assert result.verification_status is ConditionalVerificationStatus.DECLARED
        assert result.evidence == ("paper",)
        assert result.domain.assumption_context.contains_all(first_context)
        assert result.domain.assumption_context.contains_all(second_context)
        assert set(result.domain.complex_domain.exclusions) == {0, sp.I}


def test_negation_preserves_all_semantic_metadata() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain_for(lam, lam=lam, exclusions=(0,))
    formal = MellinExpression(
        lam,
        domain,
        variables=(domain.frequency,),
        evidence=("source",),
    )

    negated = -formal

    assert negated.as_expr() == -lam
    assert negated.domain == formal.domain
    assert negated.evidence == formal.evidence


def test_parameter_substitution_updates_formula_domain_and_roles() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    context = AssumptionContext((kappa > 0, kappa < 1))
    domain = _domain_for(
        lam + kappa,
        lam=lam,
        parameters=(kappa,),
        context=context,
    )
    formal = MellinExpression(
        lam + kappa,
        domain,
        variables=(domain.frequency,),
        parameters=domain.parameters,
    )

    substituted = formal.substitute({kappa: sp.Rational(1, 2)})

    assert substituted.as_expr() == lam + sp.Rational(1, 2)
    assert substituted.parameters == ()
    assert substituted.domain.parameters == ()


def test_frequency_rename_preserves_frequency_role() -> None:
    lam, mu = sp.symbols("lambda mu")
    domain = _domain_for(lam, lam=lam)
    formal = MellinExpression(lam, domain, variables=(domain.frequency,))

    renamed = formal.substitute({lam: mu})

    assert renamed.as_expr() == mu
    assert renamed.domain.frequency.symbol == mu
    assert renamed.variables == (renamed.domain.frequency,)


def test_substitution_cannot_change_a_role() -> None:
    lam, x = sp.symbols("lambda x")
    domain = _domain_for(lam + x, lam=lam, spatial=(x,))
    formal = MellinExpression(
        lam + x,
        domain,
        variables=(domain.frequency, domain.spatial_variables[0]),
    )

    with pytest.raises(ValueError, match="role"):
        formal.substitute({x: lam})


def test_differentiation_adds_pole_order_without_losing_old_metadata() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    context = AssumptionContext((kappa > 0, kappa < 1))
    formula = sp.coth(sp.pi * (lam + sp.I * kappa))
    domain = _domain_for(
        formula,
        lam=lam,
        parameters=(kappa,),
        context=context,
    )
    formal = MellinExpression(
        formula,
        domain,
        variables=(domain.frequency,),
        parameters=domain.parameters,
    )

    derivative = formal.differentiate(lam)

    assert sp.simplify(derivative.as_expr() + sp.pi / sp.sinh(
        sp.pi * (lam + sp.I * kappa)
    ) ** 2) == 0
    assert len(derivative.domain.singular_set.singularities) >= len(
        formal.domain.singular_set.singularities
    )
    assert derivative.domain.assumption_context == formal.domain.assumption_context


def test_formal_conjugation_does_not_drop_metadata() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain_for(lam + sp.I, lam=lam, exclusions=(0,))
    formal = MellinExpression(
        lam + sp.I,
        domain,
        variables=(domain.frequency,),
        evidence=("source",),
    )

    conjugated = formal.conjugate()

    assert conjugated.as_expr().func is sp.conjugate
    assert conjugated.domain.complex_domain.exclusions == (0,)
    assert conjugated.evidence == ("source",)


def test_bound_variables_require_explicit_dummy_metadata() -> None:
    lam = sp.Symbol("lambda")
    n = sp.Symbol("n", integer=True)
    formula = sp.Sum(n * lam, (n, 1, 3))
    domain = _domain_for(formula, lam=lam)

    with pytest.raises(UndeclaredMellinVariableError, match="dummy_variables"):
        MellinExpression(formula, domain, variables=(domain.frequency,))

    formal = MellinExpression(
        formula,
        domain,
        variables=(domain.frequency,),
        dummy_variables=(n,),
    )
    assert formal.dummy_variables == (n,)


def test_projection_does_not_reconstruct_strong_metadata() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain_for(lam, lam=lam)
    formal = MellinExpression(lam, domain, variables=(domain.frequency,))

    projected = formal.as_expr()

    assert isinstance(projected, sp.Expr)
    assert not isinstance(projected, MellinExpression)
    with pytest.raises(TypeError):
        MellinExpression(projected)  # type: ignore[call-arg]
