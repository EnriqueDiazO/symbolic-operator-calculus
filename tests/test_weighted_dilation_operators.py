import pytest
import sympy as sp

from symbolic_operator_calculus import (
    AssumptionContext,
    DilationConvention,
    MellinOperatorError,
    OperatorAtom,
    WeightNormConvention,
    WeightedDilationOperator,
    WeightedLpSpace,
)


def _space(
    *,
    p: sp.Expr = sp.Integer(2),
    delta: sp.Expr = sp.Rational(1, 4),
    assumptions: AssumptionContext | None = None,
    label: str = "X_delta",
) -> WeightedLpSpace:
    return WeightedLpSpace(
        label=label,
        p=p,
        weight_exponent=delta,
        assumptions=assumptions or AssumptionContext(),
        source="article weighted half-line convention",
    )


def _dilation(
    gamma: sp.Expr,
    space: WeightedLpSpace,
    *,
    assumptions: AssumptionContext | None = None,
    normalization_power: sp.Expr | None = None,
    name: str = "D_gamma",
) -> WeightedDilationOperator:
    return WeightedDilationOperator(
        atom=OperatorAtom(name, kind="weighted_dilation"),
        gamma=gamma,
        domain=space,
        codomain=space,
        assumptions=assumptions or AssumptionContext(),
        source="direct weighted norm calculation",
        normalization_power=normalization_power,
        evidence=("sections/03_transport_and_kernels.tex:1368-1469",),
    )


def test_identity_dilation_has_identity_scale_action_and_norm() -> None:
    space = _space()
    identity = _dilation(sp.Integer(1), space, name="D_1")
    x = sp.Symbol("x", positive=True)
    f = sp.Function("f")

    action = identity.normalization_factor * f(identity.gamma * x)

    assert action == f(x)
    assert sp.simplify(identity.normalization_factor - 1) == 0
    assert sp.simplify(identity.norm_factor - 1) == 0
    assert identity.is_isometric_normalization is True
    assert identity.convention is DilationConvention.ARGUMENT_MULTIPLICATION
    assert sp.simplify(identity.inverse().gamma - 1) == 0


def test_inverse_preserves_space_normalization_convention_and_evidence() -> None:
    space = _space(p=sp.Integer(3), delta=sp.Rational(1, 6))
    dilation = _dilation(sp.Integer(4), space)
    inverse = dilation.inverse(
        OperatorAtom("D_gamma_inverse", kind="weighted_dilation")
    )

    assert sp.simplify(dilation.gamma * inverse.gamma - 1) == 0
    assert dilation.domain == inverse.codomain == space
    assert dilation.codomain == inverse.domain == space
    assert inverse.normalization_power == dilation.normalization_power
    assert inverse.convention is dilation.convention
    assert inverse.evidence == dilation.evidence
    assert inverse.bounded is dilation.bounded
    assert sp.simplify(
        dilation.normalization_factor * inverse.normalization_factor - 1
    ) == 0
    assert sp.simplify(inverse.inverse().gamma - dilation.gamma) == 0


def test_composition_multiplies_scales_and_normalizing_factors() -> None:
    space = _space(p=sp.Integer(3), delta=sp.Rational(1, 6))
    first = _dilation(sp.Integer(2), space, name="D_2")
    second = _dilation(sp.Integer(3), space, name="D_3")
    composite = _dilation(sp.Integer(6), space, name="D_6")
    x = sp.Symbol("x", positive=True)
    f = sp.Function("f")

    sequential_action = (
        first.normalization_factor
        * second.normalization_factor
        * f(second.gamma * first.gamma * x)
    )
    composite_action = composite.normalization_factor * f(composite.gamma * x)

    assert sp.simplify(first.gamma * second.gamma - composite.gamma) == 0
    assert sp.simplify(
        first.normalization_factor
        * second.normalization_factor
        - composite.normalization_factor
    ) == 0
    assert sequential_action == composite_action


def test_kappa_and_operator_norm_are_derived_from_the_actual_weighted_norm() -> None:
    p = sp.Integer(3)
    delta = sp.Rational(1, 6)
    gamma = sp.Integer(4)
    space = _space(p=p, delta=delta)
    kappa = delta + 1 / p
    normalized = _dilation(gamma, space)
    unnormalized = _dilation(
        gamma,
        space,
        normalization_power=sp.Integer(0),
        name="V_gamma",
    )

    assert kappa == sp.Rational(1, 2)
    assert space.normalization_exponent == kappa
    assert normalized.normalization_power == kappa
    assert sp.simplify(normalized.normalization_factor - gamma**kappa) == 0
    assert sp.simplify(normalized.norm_factor - 1) == 0
    assert normalized.is_isometric_normalization is True
    assert sp.simplify(unnormalized.norm_factor - gamma ** (-kappa)) == 0
    assert unnormalized.is_isometric_normalization is False


def test_direct_integral_confirms_the_normalized_dilation_is_an_isometry() -> None:
    p = sp.Integer(2)
    delta = sp.Rational(1, 4)
    gamma = sp.Integer(4)
    kappa = delta + 1 / p
    x = sp.Symbol("x", positive=True)
    function = sp.exp(-x)

    original_norm_power = sp.integrate(
        (x**delta * function) ** p,
        (x, 0, sp.oo),
    )
    dilated_norm_power = sp.integrate(
        (x**delta * gamma**kappa * function.xreplace({x: gamma * x})) ** p,
        (x, 0, sp.oo),
    )
    unnormalized_norm_power = sp.integrate(
        (x**delta * function.xreplace({x: gamma * x})) ** p,
        (x, 0, sp.oo),
    )

    assert sp.simplify(dilated_norm_power / original_norm_power - 1) == 0
    assert sp.simplify(
        unnormalized_norm_power / original_norm_power - gamma ** (-p * kappa)
    ) == 0


@pytest.mark.parametrize("gamma", (sp.Integer(0), sp.Integer(-1), sp.Rational(-1, 3)))
def test_nonpositive_gamma_is_rejected(gamma: sp.Expr) -> None:
    with pytest.raises(MellinOperatorError, match="strictly positive"):
        _dilation(gamma, _space())


def test_symbolic_gamma_requires_an_explicit_positive_hypothesis() -> None:
    gamma = sp.Symbol("gamma", real=True)
    space = _space()

    with pytest.raises(MellinOperatorError, match="gamma>0 must be provable"):
        _dilation(gamma, space)

    context = AssumptionContext((gamma > 0,))
    dilation = _dilation(gamma, space, assumptions=context)

    assert dilation.gamma == gamma
    assert dilation.assumptions.contains_exact(gamma > 0)


@pytest.mark.parametrize(
    ("domain", "codomain"),
    (
        (_space(delta=sp.Rational(1, 4)), _space(delta=sp.Rational(1, 5))),
        (_space(p=sp.Integer(2)), _space(p=sp.Integer(3))),
    ),
)
def test_dilation_rejects_incompatible_spaces_and_weights(
    domain: WeightedLpSpace,
    codomain: WeightedLpSpace,
) -> None:
    with pytest.raises(MellinOperatorError, match="preserve one weighted space"):
        WeightedDilationOperator(
            atom=OperatorAtom("D", kind="weighted_dilation"),
            gamma=sp.Integer(2),
            domain=domain,
            codomain=codomain,
            assumptions=AssumptionContext(),
            source="incompatible test",
        )


def test_space_records_multiplicative_weight_convention_not_density_weight() -> None:
    space = _space(p=sp.Integer(5), delta=sp.Rational(-1, 10))

    assert space.norm_convention is WeightNormConvention.MULTIPLICATIVE_POWER_WEIGHT
    assert space.weight_exponent == sp.Rational(-1, 10)
    assert space.normalization_exponent == sp.Rational(1, 10)
    assert space.measure == "dx on (0,infinity)"


def test_symbolic_p_requires_the_p_greater_than_one_assumption() -> None:
    p = sp.Symbol("p", real=True)

    with pytest.raises(MellinOperatorError, match="p>1"):
        _space(p=p)

    context = AssumptionContext((p > 1,))
    space = _space(p=p, assumptions=context)

    assert space.p == p
    assert space.assumptions.contains_exact(p > 1)
