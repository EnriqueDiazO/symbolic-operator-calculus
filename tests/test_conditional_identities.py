from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    AssumptionContext,
    CertificationStatus,
    ConditionalIdentity,
    ConditionalVerificationStatus,
    ComplexDomain,
    ExactIdentityRequiredError,
    ExactIdentityScope,
    FormalIdentity,
    IdentityApplicabilityError,
    IdentityScopeError,
    IdentityVerificationIndeterminateError,
    IncompatibleDomainError,
    MembershipStatus,
    Singularity,
    SingularityAvoidance,
    SingularityKind,
    SingularityOrigin,
    SingularSet,
    SingularMetadataRequiredError,
    build_coth_conditional_identities,
    detect_singularities,
    render_conditional_identity_latex,
    render_semantic_relation_latex,
    require_applicable_identity,
    require_exact_identity,
)


def test_limited_coth_rule_records_integer_lattice_poles():
    z = sp.Symbol("z")
    singular_set = detect_singularities(sp.coth(sp.pi * z), z)
    pole = singular_set.singularities[0]
    n = sp.Symbol("n", integer=True)

    assert len(singular_set.singularities) == 1
    assert pole.kind is SingularityKind.POLE
    assert pole.order == 1
    assert pole.origin is SingularityOrigin.INTERNAL_RULE
    assert pole.location == sp.ImageSet(
        sp.Lambda(n, sp.I * n), sp.S.Integers
    )
    assert singular_set.contains_point(
        2 * sp.I, variable=z
    ) is MembershipStatus.YES


def test_limited_reciprocal_sinh_rule_records_same_simple_poles():
    z = sp.Symbol("z")
    singular_set = detect_singularities(1 / sp.sinh(sp.pi * z), z)

    assert len(singular_set.singularities) == 1
    assert singular_set.singularities[0].kind is SingularityKind.POLE
    assert singular_set.singularities[0].order == 1
    assert singular_set.contains_point(0, variable=z) is MembershipStatus.YES


def test_shifted_coth_rule_specializes_poles_to_real_lambda_zero():
    lam, kappa = sp.symbols("lambda kappa", real=True)
    singular_set = detect_singularities(
        sp.coth(sp.pi * (lam + sp.I * kappa)), lam
    )
    pole = singular_set.singularities[0]

    assert pole.location == 0
    assert pole.variable == lam
    assert pole.conditions.contains_exact(sp.Contains(kappa, sp.S.Integers))
    assert singular_set.contains_point(
        0,
        variable=lam,
        available_context=AssumptionContext((sp.Eq(kappa, 0),)),
    ) is MembershipStatus.YES


def test_open_unit_kappa_interval_records_conditioned_pole_avoidance():
    lam, kappa = sp.symbols("lambda kappa", real=True)
    context = AssumptionContext((kappa > 0, kappa < 1))
    domain = ComplexDomain.real_line(lam, assumption_context=context)
    singular_set = detect_singularities(
        sp.coth(sp.pi * (lam + sp.I * kappa)),
        lam,
        assumptions=context,
    )

    assert len(singular_set.singularities) == 1
    assert len(singular_set.avoidances) == 1
    assert singular_set.avoidances[0].conditions == context
    assert singular_set.avoidance_status(domain) is MembershipStatus.YES
    assert singular_set.singularities[0] in singular_set.singularities


@pytest.mark.parametrize("singular_value", (0, 1))
def test_integer_endpoint_kappa_makes_lambda_zero_singular(singular_value):
    lam, kappa = sp.symbols("lambda kappa", real=True)
    singular_set = detect_singularities(
        sp.coth(sp.pi * (lam + sp.I * kappa)), lam
    )
    context = AssumptionContext((sp.Eq(kappa, singular_value),))
    domain = ComplexDomain.real_line(lam, assumption_context=context)

    assert singular_set.contains_point(
        0,
        variable=lam,
        available_context=context,
    ) is MembershipStatus.YES
    assert singular_set.avoidance_status(domain) is MembershipStatus.NO


def test_shifted_pole_status_remains_undetermined_without_kappa_information():
    lam, kappa = sp.symbols("lambda kappa", real=True)
    singular_set = detect_singularities(
        1 / sp.sinh(sp.pi * (lam + sp.I * kappa)), lam
    )

    assert singular_set.contains_point(
        0,
        variable=lam,
        available_context=AssumptionContext(),
    ) is MembershipStatus.UNDETERMINED
    assert singular_set.avoidance_status(
        ComplexDomain.real_line(lam)
    ) is MembershipStatus.UNDETERMINED


def test_noninteger_power_records_principal_branch_point_and_cut():
    z, alpha = sp.symbols("z alpha")
    singular_set = detect_singularities(z**alpha, z)

    assert tuple(item.kind for item in singular_set.singularities) == (
        SingularityKind.BRANCH_POINT,
        SingularityKind.BRANCH_CUT,
    )
    assert singular_set.singularities[0].location == 0
    assert singular_set.singularities[1].location == sp.Interval(-sp.oo, 0)
    assert all(
        item.origin is SingularityOrigin.INTERNAL_RULE
        for item in singular_set.singularities
    )


def test_logarithm_records_sympy_principal_branch_without_choosing_another():
    z = sp.Symbol("z")
    singular_set = detect_singularities(sp.log(z), z)

    assert {item.kind for item in singular_set.singularities} == {
        SingularityKind.BRANCH_POINT,
        SingularityKind.BRANCH_CUT,
    }
    assert all("principal SymPy" in item.description for item in singular_set.singularities)


def test_positive_gamma_power_treatment_requires_explicit_context():
    gamma = sp.Symbol("gamma")
    lam = sp.Symbol("lambda", real=True)
    expression = gamma ** (sp.I * lam)
    no_context = detect_singularities(expression, lam)
    explicit = AssumptionContext((gamma > 0,))
    with_context = detect_singularities(
        expression,
        lam,
        assumptions=explicit,
    )

    assert no_context.positive_base_treatment_status(
        gamma, AssumptionContext()
    ) is MembershipStatus.UNDETERMINED
    assert with_context.positive_base_treatment_status(
        gamma, explicit
    ) is MembershipStatus.YES
    assert with_context.assumption_context == explicit
    assert len(with_context.avoidances) == 2
    assert len(with_context.singularities) == 2


def test_singular_set_union_and_substitution_never_remove_records():
    z, w = sp.symbols("z w")
    pole = Singularity(0, SingularityKind.POLE, z, order=1)
    branch = Singularity(0, SingularityKind.BRANCH_POINT, z)
    singular_set = SingularSet((pole,)).union(SingularSet((branch,)))

    substituted = singular_set.substitute({z: w})

    assert len(singular_set.singularities) == 2
    assert len(substituted.singularities) == 2
    assert {item.variable for item in substituted.singularities} == {w}
    assert not hasattr(singular_set, "remove")


def test_declared_detected_and_external_origins_remain_distinct():
    z = sp.Symbol("z")
    declared = Singularity(
        0,
        SingularityKind.UNKNOWN,
        z,
        origin=SingularityOrigin.USER_DECLARED,
    )
    external = Singularity(
        1,
        SingularityKind.REMOVABLE,
        z,
        origin=SingularityOrigin.EXTERNAL_EVIDENCE,
        evidence="reference lemma",
    )
    detected = detect_singularities(sp.coth(sp.pi * z), z).singularities[0]

    assert declared.origin is SingularityOrigin.USER_DECLARED
    assert detected.origin is SingularityOrigin.INTERNAL_RULE
    assert external.origin is SingularityOrigin.EXTERNAL_EVIDENCE
    assert external.evidence == "reference lemma"
    assert external.certification_status is CertificationStatus.EVIDENCE_SUPPLIED
    assert detected.certification_status is CertificationStatus.UNCERTIFIED


def test_singularities_and_avoidances_are_immutable_and_hashable():
    z = sp.Symbol("z")
    pole = Singularity(0, SingularityKind.POLE, z, order=1)
    avoidance = SingularityAvoidance(
        0,
        SingularityKind.POLE,
        z,
        AssumptionContext((sp.Symbol("kappa") > 0,)),
    )
    singular_set = SingularSet((pole,), avoidances=(avoidance,))

    assert hash(pole)
    assert hash(avoidance)
    assert hash(singular_set)
    with pytest.raises(FrozenInstanceError):
        pole.location = 1


def test_detection_source_contains_no_numerical_sampling_strategy():
    source = (
        Path(__file__).parents[1]
        / "src"
        / "symbolic_operator_calculus"
        / "singularities.py"
    ).read_text(encoding="utf-8")

    assert "lambdify" not in source
    assert "random" not in source
    assert "numpy" not in source


@pytest.fixture
def coth_family():
    lam, kappa = sp.symbols("lambda kappa", real=True)
    return lam, kappa, build_coth_conditional_identities(lam, kappa)


def test_reference_coth_family_has_all_scalar_conditional_identities(coth_family):
    lam, kappa, family = coth_family
    argument = sp.pi * (lam + sp.I * kappa)
    denominator = sp.sinh(argument) ** 2

    assert family.symbol == sp.coth(argument)
    assert family.p_minus == (1 - family.symbol) / 2
    assert family.p_plus == (1 + family.symbol) / 2
    assert sp.simplify(
        family.product.lhs - family.p_minus * family.p_plus
    ) == 0
    assert family.product.rhs == -1 / (4 * denominator)
    assert family.symbol_derivative.rhs == -sp.pi / denominator
    assert family.p_minus_derivative.rhs == sp.pi / (2 * denominator)
    assert family.p_plus_derivative.rhs == -sp.pi / (2 * denominator)
    for identity in (
        family.product,
        family.symbol_derivative,
        family.p_minus_derivative,
        family.p_plus_derivative,
    ):
        assert identity.assumption_context.assumptions == (kappa > 0, kappa < 1)
        assert identity.domain.variable == lam
        assert identity.scope is ExactIdentityScope.SCALAR_SYMBOLIC
        assert identity.verification_status is (
            ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
        )
        assert identity.singular_set.singularities


def test_conditional_identity_is_applicable_with_all_conditions(coth_family):
    _, _, family = coth_family
    identity = family.product

    assert require_applicable_identity(
        identity,
        available_context=identity.assumption_context,
        available_domain=identity.domain,
    ) is identity


@pytest.mark.parametrize("retained_index", (0, 1))
def test_applicability_rejects_each_missing_kappa_bound(coth_family, retained_index):
    lam, _, family = coth_family
    retained = family.product.assumption_context.assumptions[retained_index]
    incomplete = AssumptionContext((retained,))
    domain = ComplexDomain.real_line(lam, assumption_context=incomplete)

    with pytest.raises(IdentityApplicabilityError, match="omits"):
        require_applicable_identity(
            family.product,
            available_context=incomplete,
            available_domain=domain,
        )


@pytest.mark.parametrize("singular_value", (0, 1))
def test_applicability_rejects_singular_kappa_endpoints(coth_family, singular_value):
    lam, kappa, family = coth_family
    contradictory = AssumptionContext(
        (kappa > 0, kappa < 1, sp.Eq(kappa, singular_value))
    )

    assert contradictory.consistency_status.value == "inconsistent"
    with pytest.raises(IncompatibleDomainError):
        ComplexDomain.real_line(lam, assumption_context=contradictory)
    with pytest.raises(IdentityApplicabilityError):
        family.product.substitute({kappa: singular_value})


def test_applicability_rejects_a_domain_broader_than_real_lambda(coth_family):
    lam, _, family = coth_family
    broad = ComplexDomain.complex_plane(
        lam,
        assumption_context=family.product.assumption_context,
    )

    with pytest.raises(IdentityApplicabilityError, match="broader"):
        require_applicable_identity(
            family.product,
            available_context=family.product.assumption_context,
            available_domain=broad,
        )


def test_conditional_identity_cannot_be_used_as_exact_or_operatorial(coth_family):
    _, _, family = coth_family

    with pytest.raises(ExactIdentityRequiredError):
        require_exact_identity(family.product)
    with pytest.raises(IdentityScopeError):
        require_applicable_identity(
            family.product,
            available_context=family.product.assumption_context,
            available_domain=family.product.domain,
            required_scope=ExactIdentityScope.OPERATORIAL,
        )
    with pytest.raises(IdentityScopeError):
        ConditionalIdentity(
            lhs=family.product.lhs,
            rhs=family.product.rhs,
            assumption_context=family.product.assumption_context,
            domain=family.product.domain,
            singular_set=family.product.singular_set,
            scope=ExactIdentityScope.OPERATORIAL,
        )


def test_differentiation_preserves_context_domain_and_all_singularities(coth_family):
    lam, _, family = coth_family
    differentiated = family.product.differentiate(lam)

    assert differentiated.assumption_context == family.product.assumption_context
    assert differentiated.domain == family.product.domain
    assert set(family.product.singular_set.singularities).issubset(
        differentiated.singular_set.singularities
    )
    assert differentiated.scope is ExactIdentityScope.SCALAR_SYMBOLIC
    assert differentiated.verification_status is (
        ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
    )
    with pytest.raises(IdentityScopeError):
        family.product.differentiate(sp.Symbol("x"))


def test_safe_substitution_updates_every_component_and_retains_conditions(coth_family):
    _, kappa, family = coth_family
    specialized = family.product.substitute({kappa: sp.Rational(1, 2)})

    assert specialized.assumption_context.assumptions == (sp.true,)
    assert specialized.domain.assumption_context.assumptions == (sp.true,)
    assert len(specialized.singular_set.singularities) == len(
        family.product.singular_set.singularities
    )
    assert kappa not in specialized.lhs.free_symbols
    assert specialized.verification_status is (
        ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
    )


def test_identity_equality_hashing_and_text_are_deterministic(coth_family):
    lam, kappa, family = coth_family
    rebuilt = build_coth_conditional_identities(lam, kappa).product

    assert family.product == rebuilt
    assert hash(family.product) == hash(rebuilt)
    text = str(family.product)
    assert "ConditionalIdentity" in text
    assert "kappa > 0" in text
    assert "scalar_symbolic" in text
    assert "symbolically_checked_under_assumptions" in text


def test_conditional_latex_keeps_conditions_domain_singularities_and_status(coth_family):
    _, _, family = coth_family
    latex = render_conditional_identity_latex(family.product)

    assert r"\text{conditional scalar identity}" in latex
    assert r"\kappa > 0" in latex
    assert r"\kappa < 1" in latex
    assert r"\lambda\in\mathbb{R}" in latex
    assert r"\text{pole}" in latex
    assert "symbolically checked under assumptions" in latex


def test_semantic_renderer_visually_distinguishes_formal_and_conditional(coth_family):
    _, _, family = coth_family
    conditional = render_semantic_relation_latex(family.product)
    formal = render_semantic_relation_latex(FormalIdentity(sp.Symbol("a"), sp.Symbol("b")))

    assert "conditional scalar identity" in conditional
    assert r"\text{formal}" in formal
    assert conditional != formal
    assert "Exact definition" not in conditional


def test_symbolic_check_refuses_an_indeterminate_residual():
    x = sp.Symbol("x", real=True)
    f = sp.Function("f")
    g = sp.Function("g")
    context = AssumptionContext((x > 0,))
    identity = ConditionalIdentity(
        lhs=f(x),
        rhs=g(x),
        assumption_context=context,
        domain=ComplexDomain.real_line(x, assumption_context=context),
        singular_set=SingularSet(assumption_context=context),
        scope=ExactIdentityScope.SCALAR_SYMBOLIC,
    )

    with pytest.raises(IdentityVerificationIndeterminateError):
        identity.symbolically_check()


def test_evidence_does_not_automatically_promote_a_declared_identity(coth_family):
    _, _, family = coth_family
    source = family.product
    declared = ConditionalIdentity(
        lhs=source.lhs,
        rhs=source.rhs,
        assumption_context=source.assumption_context,
        domain=source.domain,
        singular_set=source.singular_set,
        scope=source.scope,
        evidence="external paper",
    )
    supplied = ConditionalIdentity(
        lhs=source.lhs,
        rhs=source.rhs,
        assumption_context=source.assumption_context,
        domain=source.domain,
        singular_set=source.singular_set,
        scope=source.scope,
        verification_status=ConditionalVerificationStatus.EVIDENCE_SUPPLIED,
        evidence="external paper",
    )

    assert declared.verification_status is ConditionalVerificationStatus.DECLARED
    assert supplied.verification_status is (
        ConditionalVerificationStatus.EVIDENCE_SUPPLIED
    )
    assert supplied.verification_status is not (
        ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
    )


def test_limited_branch_detection_is_required_for_conditional_identity():
    z, alpha = sp.symbols("z alpha")
    context = AssumptionContext((alpha > 0,))
    domain = ComplexDomain.complex_plane(z, assumption_context=context)

    with pytest.raises(SingularMetadataRequiredError, match="branch"):
        ConditionalIdentity(
            lhs=z**alpha,
            rhs=sp.exp(alpha * sp.log(z)),
            assumption_context=context,
            domain=domain,
            singular_set=SingularSet(assumption_context=context),
            scope=ExactIdentityScope.SCALAR_SYMBOLIC,
        )


def test_as_expr_is_one_way_metadata_loss_not_a_reconstruction(coth_family):
    _, _, family = coth_family
    projected = family.product.as_expr()

    assert isinstance(projected, sp.Equality)
    assert not isinstance(projected, ConditionalIdentity)
    assert not hasattr(projected, "assumption_context")
    with pytest.raises(TypeError):
        require_applicable_identity(
            projected,
            available_context=family.product.assumption_context,
            available_domain=family.product.domain,
        )


def test_combination_rejects_contradictory_contexts(coth_family):
    lam, kappa, family = coth_family
    opposite_context = AssumptionContext((kappa <= 0,))
    opposite = ConditionalIdentity(
        lhs=sp.Integer(1),
        rhs=sp.Integer(1),
        assumption_context=opposite_context,
        domain=ComplexDomain.real_line(lam, assumption_context=opposite_context),
        singular_set=SingularSet(assumption_context=opposite_context),
        scope=ExactIdentityScope.SCALAR_SYMBOLIC,
    )

    combined = family.product.assumption_context.combine(opposite_context)
    assert combined.consistency_status.value == "inconsistent"
    with pytest.raises(IdentityApplicabilityError, match="contradictory"):
        family.product.multiply(opposite)


def test_applicability_rejects_domain_for_another_principal_variable(coth_family):
    _, _, family = coth_family
    other = sp.Symbol("mu", real=True)
    wrong_domain = ComplexDomain.real_line(
        other,
        assumption_context=family.product.assumption_context,
    )

    with pytest.raises(IdentityApplicabilityError, match="broader"):
        require_applicable_identity(
            family.product,
            available_context=family.product.assumption_context,
            available_domain=wrong_domain,
        )


def test_undetermined_consistency_is_never_accepted_as_consistent():
    x, a, b = sp.symbols("x a b", real=True)
    context = AssumptionContext((sp.Eq(a, b + 1),))
    identity = ConditionalIdentity(
        lhs=sp.Integer(1),
        rhs=sp.Integer(1),
        assumption_context=context,
        domain=ComplexDomain.real_line(x, assumption_context=context),
        singular_set=SingularSet(assumption_context=context),
        scope=ExactIdentityScope.SCALAR_SYMBOLIC,
    )

    assert context.consistency_status.value == "undetermined"
    with pytest.raises(IdentityVerificationIndeterminateError):
        identity.symbolically_check()
    with pytest.raises(IdentityApplicabilityError, match="undetermined"):
        require_applicable_identity(
            identity,
            available_context=context,
            available_domain=identity.domain,
        )


def test_multiplication_unions_metadata_and_never_strengthens_status(coth_family):
    _, _, family = coth_family
    left = family.symbol_derivative
    right = family.p_plus_derivative

    product = left.multiply(right)

    assert product.assumption_context == left.assumption_context
    assert product.domain == left.domain
    assert set(left.singular_set.singularities).issubset(
        product.singular_set.singularities
    )
    assert set(right.singular_set.singularities).issubset(
        product.singular_set.singularities
    )
    assert product.verification_status is ConditionalVerificationStatus.DECLARED
    assert product.evidence == (*left.evidence, *right.evidence)


def test_evidence_supplied_status_requires_actual_evidence(coth_family):
    _, _, family = coth_family
    source = family.product

    with pytest.raises(ValueError, match="requires an evidence"):
        ConditionalIdentity(
            lhs=source.lhs,
            rhs=source.rhs,
            assumption_context=source.assumption_context,
            domain=source.domain,
            singular_set=source.singular_set,
            scope=source.scope,
            verification_status=ConditionalVerificationStatus.EVIDENCE_SUPPLIED,
        )


def test_gamma_power_applicability_requires_explicit_positive_base_context():
    gamma, lam = sp.symbols("gamma lambda", real=True)
    lhs = gamma ** (sp.I * lam)
    rhs = sp.exp(sp.I * lam * sp.log(gamma))
    empty = AssumptionContext()
    bare_domain = ComplexDomain.real_line(lam)
    bare = ConditionalIdentity(
        lhs=lhs,
        rhs=rhs,
        assumption_context=empty,
        domain=bare_domain,
        singular_set=detect_singularities(lhs - rhs, lam),
        scope=ExactIdentityScope.SCALAR_SYMBOLIC,
    )

    with pytest.raises(IdentityApplicabilityError, match="Singularities|singularities"):
        require_applicable_identity(
            bare,
            available_context=empty,
            available_domain=bare_domain,
        )

    positive = AssumptionContext((gamma > 0,))
    positive_domain = ComplexDomain.real_line(lam, assumption_context=positive)
    conditioned = ConditionalIdentity(
        lhs=lhs,
        rhs=rhs,
        assumption_context=positive,
        domain=positive_domain,
        singular_set=detect_singularities(lhs - rhs, lam, assumptions=positive),
        scope=ExactIdentityScope.SCALAR_SYMBOLIC,
    )

    assert require_applicable_identity(
        conditioned,
        available_context=positive,
        available_domain=positive_domain,
    ) is conditioned
    assert conditioned.verification_status is ConditionalVerificationStatus.DECLARED


def test_composite_logarithm_uses_preimage_branch_sets():
    z = sp.Symbol("z")
    singular_set = detect_singularities(sp.log(z + 1), z)
    point, cut = singular_set.singularities

    assert isinstance(point.location, sp.ConditionSet)
    assert isinstance(cut.location, sp.ConditionSet)
    assert point.variable == z
    assert cut.variable == z
    ExactIdentityRequiredError,
    ExactIdentityScope,
    FormalIdentity,
    IdentityApplicabilityError,
    IdentityScopeError,
    IdentityVerificationIndeterminateError,
