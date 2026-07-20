from dataclasses import FrozenInstanceError

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    AssumptionContext,
    CertificationStatus,
    ComplexDomain,
    ConsistencyStatus,
    DomainRegionKind,
    MembershipStatus,
)
from symbolic_operator_calculus.domains import (
    DomainSubstitutionError,
    DomainVariableMismatchError,
    IncompatibleDomainError,
    InvalidAssumptionError,
)


def test_empty_assumption_context_is_immutable_and_consistent():
    context = AssumptionContext()

    assert context.is_empty
    assert context.assumptions == ()
    assert context.as_boolean_expression() is sp.true
    assert context.consistency_status is ConsistencyStatus.CONSISTENT
    with pytest.raises(FrozenInstanceError):
        context.assumptions = (sp.true,)


def test_assumption_context_deduplicates_exact_relations_in_presentation_order():
    p, kappa = sp.symbols("p kappa")
    context = AssumptionContext((p > 1, kappa > 0, p > 1, kappa < 1))

    assert context.assumptions == (p > 1, kappa > 0, kappa < 1)
    assert str(context) == (
        "AssumptionContext(p > 1; kappa > 0; kappa < 1; consistent)"
    )


@pytest.mark.parametrize("invalid", (sp.Symbol("x"), sp.Integer(1), "x", object()))
def test_assumption_context_rejects_non_propositions(invalid):
    with pytest.raises(InvalidAssumptionError):
        AssumptionContext((invalid,))


def test_context_combination_retains_every_distinct_hypothesis():
    p, kappa = sp.symbols("p kappa")
    left = AssumptionContext((p > 1, kappa > 0))
    right = AssumptionContext((kappa > 0, kappa < 1))

    combined = left.combine(right)

    assert combined.assumptions == (p > 1, kappa > 0, kappa < 1)
    assert left.assumptions == (p > 1, kappa > 0)
    assert right.assumptions == (kappa > 0, kappa < 1)


def test_assumption_context_hash_is_stable_for_equal_contexts():
    x = sp.Symbol("x")
    first = AssumptionContext((x > 0, x < 1))
    second = AssumptionContext((x > 0, x < 1, x > 0))

    assert first == second
    assert hash(first) == hash(second)
    assert {first, second} == {first}


@pytest.mark.parametrize(
    "assumptions",
    (
        (sp.false,),
        (sp.true, sp.false),
        lambda x: (x > 0, x <= 0),
        lambda x: (x < 2, x >= 2),
        lambda x: (sp.Eq(x, 3), sp.Ne(x, 3)),
        lambda x: (x > 1, x < 0),
    ),
)
def test_elementary_contradictions_are_detected(assumptions):
    x = sp.Symbol("x")
    values = assumptions(x) if callable(assumptions) else assumptions

    assert (
        AssumptionContext(values).consistency_status
        is ConsistencyStatus.INCONSISTENT
    )


def test_parameter_dependent_consistency_is_undetermined_not_consistent():
    p, kappa, delta = sp.symbols("p kappa delta")
    context = AssumptionContext(
        (p > 1, kappa > 0, kappa < 1, sp.Eq(kappa, delta + 1 / p))
    )

    assert context.consistency_status is ConsistencyStatus.UNDETERMINED
    assert context.consistency_status is not ConsistencyStatus.CONSISTENT


def test_context_never_needs_symbolic_bool_conversion(monkeypatch):
    x = sp.Symbol("x")

    def forbidden_bool(self):
        raise AssertionError("symbolic bool conversion was attempted")

    monkeypatch.setattr(type(x > 0), "__bool__", forbidden_bool)
    context = AssumptionContext((x > 0, x < 1))

    assert context.contains_exact(x > 0)
    assert context.consistency_status is ConsistencyStatus.CONSISTENT


def test_reference_hypotheses_are_preserved_exactly():
    p, kappa, delta = sp.symbols("p kappa delta")
    assumptions = (
        p > 1,
        kappa > 0,
        kappa < 1,
        sp.Eq(kappa, delta + 1 / p),
    )

    context = AssumptionContext(assumptions)

    assert context.assumptions == assumptions
    assert all(context.contains_exact(item) for item in assumptions)
    assert context.as_boolean_expression() == sp.And(*assumptions)


def test_substitution_rechecks_hypotheses_instead_of_dropping_them():
    kappa = sp.Symbol("kappa")
    context = AssumptionContext((kappa > 0, kappa < 1))

    substituted = context.substitute({kappa: 0})

    assert substituted.assumptions == (sp.false, sp.true)
    assert substituted.consistency_status is ConsistencyStatus.INCONSISTENT


def test_real_line_construction_has_explicit_region_and_metadata():
    lam = sp.Symbol("lambda", real=True)
    domain = ComplexDomain.real_line(
        lam,
        description="real spectral line",
        evidence="external domain declaration",
    )

    assert domain.variable == lam
    assert tuple(region.kind for region in domain.regions) == (
        DomainRegionKind.REAL_LINE,
    )
    assert domain.description == "real spectral line"
    assert domain.evidence == ("external domain declaration",)
    assert domain.certification_status is CertificationStatus.EVIDENCE_SUPPLIED


def test_vertical_and_horizontal_strip_construction_is_explicit():
    z = sp.Symbol("z")
    beta = sp.Symbol("beta", real=True)
    vertical = ComplexDomain.vertical_strip(z, -sp.pi, sp.pi)
    horizontal = ComplexDomain.horizontal_strip(z, 0, 1).with_assumptions(
        beta > -sp.pi,
        beta < sp.pi,
    )

    assert vertical.regions[0].kind is DomainRegionKind.VERTICAL_STRIP
    assert vertical.regions[0].lower == -sp.pi
    assert vertical.regions[0].upper == sp.pi
    assert horizontal.regions[0].kind is DomainRegionKind.HORIZONTAL_STRIP
    assert horizontal.assumption_context.assumptions == (
        beta > -sp.pi,
        beta < sp.pi,
    )


def test_domain_intersection_preserves_assumptions_exclusions_and_evidence():
    z = sp.Symbol("z")
    kappa = sp.Symbol("kappa")
    left = ComplexDomain.vertical_strip(
        z,
        -2,
        2,
        assumption_context=AssumptionContext((kappa > 0,)),
        exclusions=(sp.Integer(0),),
        evidence=("left evidence",),
    )
    right = ComplexDomain.horizontal_strip(
        z,
        -1,
        1,
        assumption_context=AssumptionContext((kappa < 1,)),
        exclusions=(sp.I,),
        evidence=("right evidence",),
    )

    intersection = left.intersect(right)

    assert intersection.assumption_context.assumptions == (kappa > 0, kappa < 1)
    assert intersection.exclusions == (0, sp.I)
    assert tuple(region.kind for region in intersection.regions) == (
        DomainRegionKind.VERTICAL_STRIP,
        DomainRegionKind.HORIZONTAL_STRIP,
    )
    assert intersection.evidence == ("left evidence", "right evidence")
    assert intersection.is_subset_of(left) is MembershipStatus.YES
    assert left.is_subset_of(intersection) is MembershipStatus.NO


@pytest.mark.parametrize(
    "builder",
    (
        lambda z: ComplexDomain.real_line(z).intersect(
            ComplexDomain.upper_half_plane(z)
        ),
        lambda z: ComplexDomain.upper_half_plane(z).intersect(
            ComplexDomain.lower_half_plane(z)
        ),
        lambda z: ComplexDomain.vertical_strip(z, -2, -1).intersect(
            ComplexDomain.vertical_strip(z, 0, 1)
        ),
        lambda z: ComplexDomain.real_line(z).intersect(
            ComplexDomain.horizontal_strip(z, 0, 1)
        ),
    ),
)
def test_elementarily_incompatible_domains_are_rejected(builder):
    with pytest.raises(IncompatibleDomainError):
        builder(sp.Symbol("z"))


def test_domain_membership_returns_yes_no_and_undetermined():
    z = sp.Symbol("z")
    unknown = sp.Symbol("w")
    domain = ComplexDomain.upper_half_plane(z).with_exclusions(1 + sp.I)

    assert domain.contains_point(2 + 3 * sp.I) is MembershipStatus.YES
    assert domain.contains_point(2 - sp.I) is MembershipStatus.NO
    assert domain.contains_point(1 + sp.I) is MembershipStatus.NO
    assert domain.contains_point(unknown) is MembershipStatus.UNDETERMINED


def test_assumption_defined_domain_uses_only_explicit_condition():
    z = sp.Symbol("z")
    domain = ComplexDomain.defined_by(z, sp.Abs(z) < 2)

    assert domain.contains_point(1) is MembershipStatus.YES
    assert domain.contains_point(3) is MembershipStatus.NO
    assert domain.contains_point(sp.Symbol("w")) is MembershipStatus.UNDETERMINED


def test_undetermined_domain_assumptions_prevent_positive_membership_claim():
    z, a, b = sp.symbols("z a b")
    domain = ComplexDomain.real_line(
        z,
        assumption_context=AssumptionContext((sp.Eq(a, b + 1),)),
    )

    assert domain.contains_point(1) is MembershipStatus.UNDETERMINED


@pytest.mark.parametrize("invalid", (1, sp.I, sp.sin(sp.Symbol("z")), "z"))
def test_domain_rejects_invalid_principal_variable(invalid):
    with pytest.raises(TypeError):
        ComplexDomain.real_line(invalid)


def test_domain_operations_never_enlarge_original_domain():
    z = sp.Symbol("z")
    original = ComplexDomain.vertical_strip(z, -2, 2)
    narrowed = original.with_exclusions(0).with_assumptions(sp.Symbol("k") > 0)

    assert narrowed.is_subset_of(original) is MembershipStatus.YES
    assert original.is_subset_of(narrowed) is MembershipStatus.NO
    assert original.exclusions == ()
    assert original.assumption_context.is_empty


def test_domain_rejects_mismatched_variables_and_capture_prone_substitution():
    z, w, kappa = sp.symbols("z w kappa")
    with pytest.raises(DomainVariableMismatchError):
        ComplexDomain.real_line(z).intersect(ComplexDomain.real_line(w))
    with pytest.raises(DomainSubstitutionError):
        ComplexDomain.real_line(z).with_assumptions(kappa > 0).substitute(
            {kappa: z + 1}
        )
    with pytest.raises(DomainSubstitutionError):
        ComplexDomain.real_line(z).substitute({z: z + 1})


def test_domain_hash_and_string_representation_are_deterministic():
    z, kappa = sp.symbols("z kappa")
    first = ComplexDomain.real_line(z).with_assumptions(kappa > 0, kappa < 1)
    second = ComplexDomain.real_line(
        z,
        assumption_context=AssumptionContext((kappa > 0, kappa < 1)),
    )

    assert first == second
    assert hash(first) == hash(second)
    assert str(first) == (
        "ComplexDomain(variable=z; region=real_line; assumptions=2; "
        "exclusions=none)"
    )
