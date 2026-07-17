from dataclasses import FrozenInstanceError

import pytest
import sympy as sp

import symbolic_operator_calculus as soc
import symbolic_operator_calculus.relative_wiener_hopf as rwh
from symbolic_operator_calculus.substitution import collect_bound_symbols


@pytest.fixture(scope="module")
def variables():
    gamma1, gamma2, decay = sp.symbols(
        "gamma1 gamma2 decay",
        positive=True,
    )
    x, y, s, z, u = sp.symbols("x y s z u")
    return gamma1, gamma2, decay, x, y, s, z, u


@pytest.fixture(scope="module")
def trace12(variables):
    gamma1, gamma2, decay, *_ = variables
    return rwh.build_relative_wiener_hopf_trace(
        1,
        2,
        gamma1,
        gamma2,
        0,
        decay,
    )


@pytest.fixture(scope="module")
def trace21(variables):
    gamma1, gamma2, decay, *_ = variables
    return rwh.build_relative_wiener_hopf_trace(
        2,
        1,
        gamma2,
        gamma1,
        decay,
        0,
    )


def _kernel_at(factor, output, input_):
    model = factor.factorization
    kernel = {
        "original": model.original_kernel,
        "left": model.left_convolution_kernel,
        "right": model.right_convolution_kernel,
    }[factor.placement]
    return kernel.xreplace({model.time_variable: output - input_})


def _computed_actions(trace, operand, output, input_):
    return tuple(
        rwh.apply_ordered_relative_product(
            product,
            operand,
            output_variable=output,
            input_variable=input_,
        )
        for product in (
            trace.identity.original,
            trace.identity.left,
            trace.identity.right,
        )
    )


def test_new_action_api_is_publicly_exported():
    for name in (
        "RelativeProductAction",
        "RelativeProductActionVerification",
        "apply_relative_factor",
        "apply_ordered_relative_product",
        "verify_relative_product_actions",
    ):
        assert getattr(soc, name) is getattr(rwh, name)


def test_direct_branch_dilation_action(trace12, variables):
    gamma1, _, _, x, *_ = variables
    f = sp.Function("f")
    factor = trace12.original_operator.factors[0]

    result = rwh.apply_relative_factor(factor, f(x), output_variable=x)
    expression = sp.exp(-x) + x**2
    expression_result = rwh.apply_relative_factor(
        factor,
        expression,
        output_variable=x,
    )

    assert result == f(gamma1 * x)
    assert expression_result == sp.exp(-gamma1 * x) + gamma1**2 * x**2


def test_inverse_branch_dilation_action(trace12, variables):
    _, gamma2, _, x, *_ = variables
    f = sp.Function("f")
    factor = trace12.original_operator.factors[2]

    result = rwh.apply_relative_factor(factor, f(x), output_variable=x)

    assert result == f(x / gamma2)


def test_relative_dilation_action(trace12, variables):
    gamma1, gamma2, _, x, *_ = variables
    f = sp.Function("f")
    factor = trace12.left_operator.factors[0]

    result = rwh.apply_relative_factor(factor, f(x), output_variable=x)

    assert result == f(gamma1 * x / gamma2)


@pytest.mark.parametrize(
    ("product_name", "factor_position", "placement"),
    [
        ("original_operator", 1, "original"),
        ("left_operator", 1, "left"),
        ("right_operator", 0, "right"),
    ],
)
def test_each_wiener_hopf_placement_builds_its_own_kernel_integral(
    trace12,
    variables,
    product_name,
    factor_position,
    placement,
):
    *_, x, _, _, _, u = variables
    f = sp.Function("f")
    factor = getattr(trace12, product_name).factors[factor_position]

    result = rwh.apply_relative_factor(
        factor,
        f(x),
        output_variable=x,
        integration_variables=(u,),
    )

    assert factor.placement == placement
    assert result == sp.Integral(
        _kernel_at(factor, x, u) * f(u),
        (u, 0, sp.oo),
    )
    assert result.limits == ((u, 0, sp.oo),)


def test_original_product_is_applied_right_to_left(trace12, variables):
    gamma1, gamma2, _, x, y, *_, u = variables
    f = sp.Function("f")

    action = rwh.apply_ordered_relative_product(
        trace12.original_operator,
        f(y),
        output_variable=x,
        input_variable=y,
        integration_variables=(u,),
    )
    wh_factor = trace12.original_operator.factors[1]

    assert action.direct_expression == sp.Integral(
        _kernel_at(wh_factor, gamma1 * x, u) * f(u / gamma2),
        (u, 0, sp.oo),
    )
    assert action.direct_expression != trace12.action.integral


def test_original_normalization_derives_the_jacobian(trace12, variables):
    gamma1, gamma2, _, x, y, *_ = variables
    f = sp.Function("f")
    action = rwh.apply_ordered_relative_product(
        trace12.original_operator,
        f(y),
        output_variable=x,
        input_variable=y,
    )
    l1 = trace12.factorization
    expected_kernel = gamma2 * l1.original_kernel.xreplace(
        {l1.time_variable: gamma1 * x - gamma2 * y}
    )

    assert sp.simplify(action.kernel - expected_kernel) == 0
    assert action.normalized_expression == sp.Integral(
        action.kernel * f(y),
        (y, 0, sp.oo),
    )


def test_left_product_direct_and_normalized_actions(trace12, variables):
    gamma1, gamma2, _, x, y, *_, u = variables
    f = sp.Function("f")
    action = rwh.apply_ordered_relative_product(
        trace12.left_operator,
        f(y),
        output_variable=x,
        input_variable=y,
        integration_variables=(u,),
    )
    wh_factor = trace12.left_operator.factors[1]

    assert action.direct_expression == sp.Integral(
        _kernel_at(wh_factor, gamma1 * x / gamma2, u) * f(u),
        (u, 0, sp.oo),
    )
    assert action.normalized_expression == sp.Integral(
        action.kernel * f(y),
        (y, 0, sp.oo),
    )


def test_right_product_direct_and_normalized_actions(trace12, variables):
    gamma1, gamma2, _, x, y, *_, u = variables
    f = sp.Function("f")
    action = rwh.apply_ordered_relative_product(
        trace12.right_operator,
        f(y),
        output_variable=x,
        input_variable=y,
        integration_variables=(u,),
    )
    wh_factor = trace12.right_operator.factors[0]

    assert action.direct_expression == sp.Integral(
        _kernel_at(wh_factor, x, u) * f(gamma1 * u / gamma2),
        (u, 0, sp.oo),
    )
    assert action.normalized_expression == sp.Integral(
        action.kernel * f(y),
        (y, 0, sp.oo),
    )
    assert action.direct_expression != action.normalized_expression


@pytest.mark.parametrize("trace_fixture", ["trace12", "trace21"])
def test_all_three_computed_kernels_are_equal_and_match_l1(
    request,
    trace_fixture,
    variables,
):
    trace = request.getfixturevalue(trace_fixture)
    *_, x, y, _, _, _ = variables
    f = sp.Function("f")
    actions = _computed_actions(trace, f(y), x, y)
    l1 = trace.factorization
    expected = l1.integral_kernel.xreplace(
        {l1.output_variable: x, l1.input_variable: y}
    )

    assert sp.simplify(actions[0].kernel - actions[1].kernel) == 0
    assert sp.simplify(actions[0].kernel - actions[2].kernel) == 0
    assert all(sp.simplify(action.kernel - expected) == 0 for action in actions)


@pytest.mark.parametrize(
    "operand_factory",
    [
        lambda variable: sp.Function("f")(variable),
        lambda variable: sp.exp(-variable),
    ],
)
def test_symbolic_function_and_concrete_expression_keep_integral_structure(
    trace12,
    variables,
    operand_factory,
):
    *_, x, y, _, _, _ = variables
    operand = operand_factory(y)

    for action in _computed_actions(trace12, operand, x, y):
        assert isinstance(action.direct_expression, sp.Integral)
        assert isinstance(action.normalized_expression, sp.Integral)
        assert action.normalized_expression.function == action.kernel * operand
        assert action.normalized_expression.limits == ((y, 0, sp.oo),)


def test_alternative_variables_are_used_without_hard_coded_x_or_y(
    trace12,
    variables,
):
    *_, s, z, _ = variables
    f = sp.Function("f")

    for action in _computed_actions(trace12, f(z), s, z):
        names = {symbol.name for symbol in action.normalized_expression.atoms(sp.Symbol)}
        assert "x" not in names
        assert "y" not in names
        assert action.output_variable == s
        assert action.input_variable == z


def test_free_u_in_operand_is_not_captured_by_direct_integral(trace12, variables):
    *_, x, y, _, _, u = variables
    f = sp.Function("f")
    operand = f(y) + u

    action = rwh.apply_ordered_relative_product(
        trace12.original_operator,
        operand,
        output_variable=x,
        input_variable=y,
    )

    assert action.direct_integration_variable != u
    assert u in action.direct_expression.free_symbols
    assert action.direct_integration_variable in collect_bound_symbols(
        action.direct_expression
    )


def test_output_symbol_already_free_in_operand_is_not_shifted_as_input(
    trace12,
    variables,
):
    gamma2 = variables[1]
    x, y = variables[3:5]
    f = sp.Function("f")
    operand = f(y) + x

    action = rwh.apply_ordered_relative_product(
        trace12.original_operator,
        operand,
        output_variable=x,
        input_variable=y,
    )
    direct_variable = action.direct_integration_variable

    assert action.direct_expression.function.has(f(direct_variable / gamma2) + x)


def test_bound_variable_collision_is_alpha_safe(trace12, variables):
    *_, x, y, _, _, u = variables
    f = sp.Function("f")
    operand = sp.Integral(f(y) + u, (u, 0, 1))

    action = rwh.apply_ordered_relative_product(
        trace12.left_operator,
        operand,
        output_variable=x,
        input_variable=y,
        integration_variables=(u,),
    )

    bound = tuple(action.direct_expression.variables)
    assert len(set(bound)) == len(bound)
    assert action.direct_integration_variable != u


def test_two_wiener_hopf_applications_use_distinct_bound_variables(
    trace12,
    variables,
):
    *_, x, _, _, _, _ = variables
    f = sp.Function("f")
    factor = trace12.left_operator.factors[1]

    first = rwh.apply_relative_factor(factor, f(x), output_variable=x)
    second = rwh.apply_relative_factor(factor, first, output_variable=x)

    bound = collect_bound_symbols(second)
    assert len(bound) == 2


def test_repeated_integration_variables_are_rejected(trace12, variables):
    *_, x, y, _, _, u = variables
    f = sp.Function("f")
    with pytest.raises(ValueError):
        rwh.apply_ordered_relative_product(
            trace12.original_operator,
            f(y),
            output_variable=x,
            input_variable=y,
            integration_variables=(u, u),
        )


def test_invalid_action_inputs_fail_early(trace12, variables):
    *_, x, y, _, _, _ = variables
    f = sp.Function("f")
    with pytest.raises(TypeError):
        rwh.apply_ordered_relative_product(
            object(),
            f(y),
            output_variable=x,
            input_variable=y,
        )
    with pytest.raises(TypeError):
        rwh.apply_ordered_relative_product(
            trace12.original_operator,
            f(y),
            output_variable="x",
            input_variable=y,
        )
    with pytest.raises(TypeError):
        rwh.apply_relative_factor(object(), f(x), output_variable=x)
    with pytest.raises(TypeError):
        rwh.apply_relative_factor(
            trace12.branch_k,
            f(x),
            output_variable=x,
        )
    with pytest.raises(ValueError):
        rwh.OrderedRelativeOperatorProduct(())


@pytest.mark.parametrize("trace_fixture", ["trace12", "trace21"])
def test_verification_evidence_is_derived_for_both_directions(
    request,
    trace_fixture,
):
    trace = request.getfixturevalue(trace_fixture)

    evidence = rwh.verify_relative_product_actions(trace.identity)

    assert isinstance(evidence, rwh.RelativeProductActionVerification)
    assert evidence.actions_verified is True
    assert evidence.identity == trace.identity
    assert evidence.original.product == trace.identity.original
    assert evidence.left.product == trace.identity.left
    assert evidence.right.product == trace.identity.right
    assert trace.identity.exact is True
    with pytest.raises(FrozenInstanceError):
        evidence.actions_verified = False


def test_computed_actions_preserve_the_existing_common_action(trace12):
    common = trace12.action
    evidence = rwh.verify_relative_product_actions(trace12.identity)

    for action in (evidence.original, evidence.left, evidence.right):
        expected = common.kernel.xreplace(
            {
                trace12.factorization.output_variable: action.output_variable,
                trace12.factorization.input_variable: action.input_variable,
            }
        )
        assert sp.simplify(action.kernel - expected) == 0


def test_product_actions_and_verification_evidence_are_hashable(trace12):
    evidence = rwh.verify_relative_product_actions(trace12.identity)

    assert isinstance(hash(evidence), int)
    assert all(
        isinstance(hash(action), int)
        for action in (evidence.original, evidence.left, evidence.right)
    )
