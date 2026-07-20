import pytest
import sympy as sp

from symbolic_operator_calculus import (
    ApproximateEquality,
    AssumptionContext,
    ComplexDomain,
    ConditionalIdentity,
    ExactIdentityRequiredError,
    MellinDependencyError,
    MellinExpression,
    MellinSymbol,
    MellinSymbolDependency,
    MellinSymbolDomain,
    SingularSet,
    mellin_frequency,
    mellin_parameter,
    output_spatial_variable,
    relative_multiplicative_variable,
    require_exact_identity,
)


def _domain(
    lam: sp.Symbol,
    *,
    spatial: sp.Symbol | None = None,
    relative: sp.Symbol | None = None,
    parameters: tuple[sp.Symbol, ...] = (),
    description: str | None = None,
) -> MellinSymbolDomain:
    return MellinSymbolDomain(
        frequency=mellin_frequency(lam),
        complex_domain=ComplexDomain.real_line(lam),
        spatial_variables=(
            (output_spatial_variable(spatial),) if spatial is not None else ()
        ),
        relative_variable=(
            relative_multiplicative_variable(relative)
            if relative is not None
            else None
        ),
        parameters=tuple(mellin_parameter(item) for item in parameters),
        singular_set=SingularSet(),
        description=description,
    )


def test_frequency_dependency_uses_role_not_name() -> None:
    x = sp.Symbol("x")
    domain = _domain(x)
    formal = MellinExpression(x, domain, variables=(domain.frequency,))
    symbol = MellinSymbol(formal, MellinSymbolDependency.FREQUENCY_ONLY)

    assert symbol.frequency.symbol == x
    assert symbol.dependency is MellinSymbolDependency.FREQUENCY_ONLY


def test_space_frequency_dependency_accepts_lambda_named_spatial_variable() -> None:
    x = sp.Symbol("x")
    lambda_named = sp.Symbol("lambda")
    domain = _domain(x, spatial=lambda_named)
    formal = MellinExpression(
        x + lambda_named,
        domain,
        variables=(domain.frequency, domain.spatial_variables[0]),
    )

    symbol = MellinSymbol(formal, MellinSymbolDependency.SPACE_FREQUENCY)

    assert symbol.spatial_variables[0].symbol == lambda_named


def test_relative_frequency_dependency_is_distinct() -> None:
    lam, r = sp.symbols("lambda r")
    domain = _domain(lam, relative=r)
    formal = MellinExpression(
        lam * r,
        domain,
        variables=(domain.frequency, domain.relative_variable),  # type: ignore[arg-type]
    )
    symbol = MellinSymbol(formal, MellinSymbolDependency.RELATIVE_FREQUENCY)

    assert symbol.relative_variable == domain.relative_variable


def test_parametric_constant_and_strict_constant_are_distinct() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    domain = _domain(lam, parameters=(kappa,))
    parametric = MellinSymbol(
        MellinExpression(kappa, domain, parameters=domain.parameters),
        MellinSymbolDependency.PARAMETRIC_CONSTANT,
    )
    constant = MellinSymbol.constant(sp.Integer(2), domain=domain)

    assert parametric.dependency is MellinSymbolDependency.PARAMETRIC_CONSTANT
    assert constant.dependency is MellinSymbolDependency.CONSTANT
    assert constant.expression.variables == ()
    assert constant.frequency not in constant.expression.variables


def test_strict_constant_rejects_a_hidden_frequency() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain(lam)

    with pytest.raises(MellinDependencyError, match="free symbols"):
        MellinSymbol.constant(lam, domain=domain)


def test_symbol_requires_typed_expression_and_explicit_dependency() -> None:
    lam = sp.Symbol("lambda")

    with pytest.raises(TypeError, match="MellinExpression"):
        MellinSymbol(lam, MellinSymbolDependency.FREQUENCY_ONLY)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        MellinSymbol(lam)  # type: ignore[call-arg]


def test_incorrect_dependency_is_rejected() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain(lam)
    formal = MellinExpression(lam, domain, variables=(domain.frequency,))

    with pytest.raises(MellinDependencyError, match="conflicts"):
        MellinSymbol(formal, MellinSymbolDependency.CONSTANT)


def test_unsupported_spatial_only_dependency_is_rejected() -> None:
    lam, x = sp.symbols("lambda x")
    domain = _domain(lam, spatial=x)
    formal = MellinExpression(
        x,
        domain,
        variables=(domain.spatial_variables[0],),
    )

    with pytest.raises(MellinDependencyError, match="supported"):
        MellinSymbol(formal, MellinSymbolDependency.SPACE_FREQUENCY)


def test_structural_and_sympy_equality_are_explicitly_different() -> None:
    lam = sp.Symbol("lambda")
    first_domain = _domain(lam, description="first")
    second_domain = _domain(lam, description="second")
    first = MellinSymbol(
        MellinExpression(lam, first_domain, variables=(first_domain.frequency,)),
        MellinSymbolDependency.FREQUENCY_ONLY,
    )
    second = MellinSymbol(
        MellinExpression(lam, second_domain, variables=(second_domain.frequency,)),
        MellinSymbolDependency.FREQUENCY_ONLY,
    )

    assert not first.structurally_equals(second)
    assert first.sympy_equals(second) is True


def test_identity_between_symbols_remains_conditional_and_not_exact() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain(lam)
    first = MellinSymbol(
        MellinExpression(lam + 1, domain, variables=(domain.frequency,)),
        MellinSymbolDependency.FREQUENCY_ONLY,
    )
    second = MellinSymbol(
        MellinExpression(1 + lam, domain, variables=(domain.frequency,)),
        MellinSymbolDependency.FREQUENCY_ONLY,
    )

    identity = first.conditional_identity(second)

    assert isinstance(identity, ConditionalIdentity)
    with pytest.raises(ExactIdentityRequiredError):
        require_exact_identity(identity)


def test_approximate_relation_remains_a_separate_type() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain(lam)
    first = MellinSymbol(
        MellinExpression(lam, domain, variables=(domain.frequency,)),
        MellinSymbolDependency.FREQUENCY_ONLY,
    )
    second = MellinSymbol(
        MellinExpression(lam + 1, domain, variables=(domain.frequency,)),
        MellinSymbolDependency.FREQUENCY_ONLY,
    )

    relation = first.approximate_equality(
        second,
        tolerance=sp.Rational(1, 10),
        norm="pointwise residual",
    )

    assert isinstance(relation, ApproximateEquality)
    assert relation.residual == -1


def test_symbol_has_no_operator_fredholm_or_compact_api() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain(lam)
    symbol = MellinSymbol(
        MellinExpression(lam, domain, variables=(domain.frequency,)),
        MellinSymbolDependency.FREQUENCY_ONLY,
    )

    for attribute in ("operator", "fredholm", "compact"):
        with pytest.raises(AttributeError):
            getattr(symbol, attribute)


def test_projection_cannot_reconstruct_a_symbol_without_metadata() -> None:
    lam = sp.Symbol("lambda")
    domain = _domain(lam)
    symbol = MellinSymbol(
        MellinExpression(lam, domain, variables=(domain.frequency,)),
        MellinSymbolDependency.FREQUENCY_ONLY,
    )

    projected = symbol.as_expr()

    assert not isinstance(projected, MellinSymbol)
    with pytest.raises(TypeError):
        MellinSymbol(projected)  # type: ignore[call-arg]


def test_symbol_assumptions_are_the_p0b_context() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    context = AssumptionContext((kappa > 0,))
    domain = MellinSymbolDomain(
        frequency=mellin_frequency(lam),
        complex_domain=ComplexDomain.real_line(lam, assumption_context=context),
        parameters=(mellin_parameter(kappa),),
        assumption_context=context,
    )
    symbol = MellinSymbol(
        MellinExpression(
            lam + kappa,
            domain,
            variables=(domain.frequency,),
            parameters=domain.parameters,
        ),
        MellinSymbolDependency.FREQUENCY_ONLY,
    )

    assert symbol.assumptions is domain.assumption_context
