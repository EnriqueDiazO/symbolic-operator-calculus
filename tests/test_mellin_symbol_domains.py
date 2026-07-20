import pytest
import sympy as sp

from symbolic_operator_calculus import (
    AssumptionContext,
    ComplexDomain,
    IncompatibleDomainError,
    MembershipStatus,
    MellinDomainError,
    MellinDomainRoleError,
    MellinFrequencyMismatchError,
    MellinSymbolDomain,
    MellinVariableRole,
    Singularity,
    SingularityKind,
    SingularSet,
    input_spatial_variable,
    mellin_frequency,
    mellin_parameter,
    output_spatial_variable,
    relative_multiplicative_variable,
)


def _domain() -> tuple[MellinSymbolDomain, dict[str, sp.Symbol]]:
    lam, x, y, r, p, delta, kappa = sp.symbols(
        "lambda x y r p delta kappa"
    )
    frequency = mellin_frequency(lam)
    context = AssumptionContext(
        (
            p > 1,
            kappa > 0,
            kappa < 1,
            sp.Eq(kappa, delta + 1 / p),
        )
    )
    complex_domain = ComplexDomain.real_line(
        lam,
        assumption_context=context,
        exclusions=(sp.I * (sp.Symbol("n", integer=True) - kappa),),
    )
    singular_set = SingularSet(
        singularities=(
            Singularity(
                location=sp.I * (sp.Symbol("n", integer=True) - kappa),
                kind=SingularityKind.POLE,
                variable=lam,
                conditions=context,
            ),
        ),
        assumption_context=context,
    )
    return (
        MellinSymbolDomain(
            frequency=frequency,
            complex_domain=complex_domain,
            spatial_variables=(
                output_spatial_variable(x),
                input_spatial_variable(y),
            ),
            relative_variable=relative_multiplicative_variable(r),
            parameters=tuple(
                mellin_parameter(item) for item in (p, delta, kappa)
            ),
            assumption_context=context,
            singular_set=singular_set,
            description="standard scalar Mellin domain",
        ),
        {
            "lambda": lam,
            "x": x,
            "y": y,
            "r": r,
            "p": p,
            "delta": delta,
            "kappa": kappa,
        },
    )


def test_domain_preserves_roles_assumptions_exclusions_and_singularities() -> None:
    domain, symbols = _domain()

    assert domain.frequency.symbol == symbols["lambda"]
    assert domain.complex_domain.variable == symbols["lambda"]
    assert domain.assumption_context.contains_exact(symbols["kappa"] > 0)
    assert domain.complex_domain.exclusions
    assert len(domain.singular_set.singularities) == 1
    assert domain.require_role(
        domain.frequency, MellinVariableRole.FREQUENCY
    ) is domain.frequency


def test_domain_rejects_complex_domain_for_another_frequency() -> None:
    lam, mu = sp.symbols("lambda mu")

    with pytest.raises(MellinFrequencyMismatchError):
        MellinSymbolDomain(
            frequency=mellin_frequency(lam),
            complex_domain=ComplexDomain.real_line(mu),
        )


def test_domain_rejects_second_frequency_in_another_slot() -> None:
    lam, mu = sp.symbols("lambda mu")

    with pytest.raises(MellinDomainRoleError):
        MellinSymbolDomain(
            frequency=mellin_frequency(lam),
            complex_domain=ComplexDomain.real_line(lam),
            spatial_variables=(mellin_frequency(mu),),
        )


def test_domain_rejects_parameter_or_relative_role_confusion() -> None:
    lam, x, r = sp.symbols("lambda x r")
    common = {
        "frequency": mellin_frequency(lam),
        "complex_domain": ComplexDomain.real_line(lam),
    }

    with pytest.raises(MellinDomainRoleError):
        MellinSymbolDomain(**common, parameters=(output_spatial_variable(x),))
    with pytest.raises(MellinDomainRoleError):
        MellinSymbolDomain(**common, relative_variable=mellin_parameter(r))


def test_domain_rejects_one_symbol_with_two_roles() -> None:
    lam = sp.Symbol("lambda")

    with pytest.raises(MellinDomainRoleError, match="multiple roles"):
        MellinSymbolDomain(
            frequency=mellin_frequency(lam),
            complex_domain=ComplexDomain.real_line(lam),
            parameters=(mellin_parameter(lam),),
        )


def test_domain_rejects_singularities_for_another_frequency() -> None:
    lam, mu = sp.symbols("lambda mu")
    singular_set = SingularSet(
        singularities=(
            Singularity(0, SingularityKind.POLE, variable=mu),
        )
    )

    with pytest.raises(MellinFrequencyMismatchError):
        MellinSymbolDomain(
            frequency=mellin_frequency(lam),
            complex_domain=ComplexDomain.real_line(lam),
            singular_set=singular_set,
        )


def test_with_assumptions_never_discards_existing_conditions() -> None:
    domain, symbols = _domain()
    narrowed = domain.with_assumptions(symbols["delta"] > -1)

    assert narrowed.assumption_context.contains_all(domain.assumption_context)
    assert narrowed.assumption_context.contains_exact(symbols["delta"] > -1)
    assert narrowed.complex_domain.exclusions == domain.complex_domain.exclusions
    assert narrowed.singular_set.singularities == domain.singular_set.singularities


def test_contradictory_assumptions_are_rejected() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    domain = MellinSymbolDomain(
        frequency=mellin_frequency(lam),
        complex_domain=ComplexDomain.real_line(lam),
        parameters=(mellin_parameter(kappa),),
        assumption_context=AssumptionContext((kappa > 0,)),
    )

    with pytest.raises(IncompatibleDomainError):
        domain.with_assumptions(kappa < 0)


def test_with_singularities_is_a_union() -> None:
    lam = sp.Symbol("lambda")
    domain = MellinSymbolDomain(
        frequency=mellin_frequency(lam),
        complex_domain=ComplexDomain.real_line(lam),
        singular_set=SingularSet(
            (Singularity(0, SingularityKind.POLE, lam),)
        ),
    )
    extra = SingularSet((Singularity(sp.I, SingularityKind.POLE, lam),))

    combined = domain.with_singularities(extra)

    assert {item.location for item in combined.singular_set.singularities} == {
        0,
        sp.I,
    }


def test_intersection_retains_both_domains_metadata() -> None:
    lam, kappa = sp.symbols("lambda kappa")
    frequency = mellin_frequency(lam)
    first = MellinSymbolDomain(
        frequency=frequency,
        complex_domain=ComplexDomain.real_line(lam).with_exclusions(0),
        parameters=(mellin_parameter(kappa),),
        assumption_context=AssumptionContext((kappa > 0,)),
    )
    second = MellinSymbolDomain(
        frequency=frequency,
        complex_domain=ComplexDomain.real_line(lam).with_exclusions(sp.I),
        parameters=(mellin_parameter(kappa),),
        assumption_context=AssumptionContext((kappa < 1,)),
    )

    intersection = first.intersect(second)

    assert intersection.assumption_context.contains_all(first.assumption_context)
    assert intersection.assumption_context.contains_all(second.assumption_context)
    assert set(intersection.complex_domain.exclusions) == {0, sp.I}
    assert first.is_compatible_with(second) is MembershipStatus.YES


def test_incompatible_frequency_reports_no_without_promotion() -> None:
    lam, mu = sp.symbols("lambda mu")
    first = MellinSymbolDomain(
        mellin_frequency(lam), ComplexDomain.real_line(lam)
    )
    second = MellinSymbolDomain(
        mellin_frequency(mu), ComplexDomain.real_line(mu)
    )

    assert first.is_compatible_with(second) is MembershipStatus.NO


def test_undetermined_assumptions_are_not_promoted_to_compatible() -> None:
    lam, p, delta, kappa = sp.symbols("lambda p delta kappa")
    frequency = mellin_frequency(lam)
    first = MellinSymbolDomain(
        frequency,
        ComplexDomain.real_line(lam),
        parameters=tuple(mellin_parameter(item) for item in (p, delta, kappa)),
        assumption_context=AssumptionContext((sp.Eq(kappa, delta + 1 / p),)),
    )
    second = MellinSymbolDomain(frequency, ComplexDomain.real_line(lam))

    assert first.is_compatible_with(second) is MembershipStatus.UNDETERMINED


def test_require_role_rejects_wrong_or_foreign_declarations() -> None:
    domain, symbols = _domain()

    with pytest.raises(MellinDomainRoleError):
        domain.require_role(domain.frequency, MellinVariableRole.PARAMETER)
    with pytest.raises(MellinDomainRoleError):
        domain.require_role(
            mellin_parameter(sp.Symbol("foreign")), MellinVariableRole.PARAMETER
        )


def test_empty_descriptions_and_nonhashable_evidence_are_rejected() -> None:
    lam = sp.Symbol("lambda")
    with pytest.raises(MellinDomainError):
        MellinSymbolDomain(
            mellin_frequency(lam),
            ComplexDomain.real_line(lam),
            description=" ",
        )
    with pytest.raises(TypeError, match="hashable"):
        MellinSymbolDomain(
            mellin_frequency(lam),
            ComplexDomain.real_line(lam),
            evidence=[],
        )
