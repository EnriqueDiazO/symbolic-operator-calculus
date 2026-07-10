from dataclasses import FrozenInstanceError
from fractions import Fraction

import pytest
import sympy as sp

import symbolic_operator_calculus.schur as schur_module
from symbolic_operator_calculus import (
    AppliedLinearCombination,
    CompactSchurActionError,
    FirstSchurCompactModelAction,
    FirstSchurReduction,
    LinearCombination,
    OperatorAtom,
    PrincipalValue,
    Product,
    a22_exact_operator,
    a22_first_schur_model,
    a22_first_schur_mod_compact_relation,
    apply_a22_first_schur_model_compact,
    apply_combined_kernel_c22,
    apply_linear_combination_ordered,
    collect_bound_symbols,
    explicit_wiener_hopf_rules,
    extract_applied_kernels,
    mvp_atomic_rules,
    positive_decay_symbol,
)


HALF = Fraction(1, 2)


def nested_integrand(kernel):
    outer_factors = sp.Mul.make_args(kernel.function)
    inner = next(factor for factor in outer_factors if isinstance(factor, sp.Integral))
    outer = sp.Mul(*(factor for factor in outer_factors if factor is not inner))
    return outer * inner.function


def compact_kernel_integrand(action, input_function):
    input_variable = action.correction.variables[0]
    input_atom = input_function(input_variable)
    factors = sp.Mul.make_args(action.correction.function)
    c22 = next(factor for factor in factors if factor != input_atom)
    outer_factors = sp.Mul.make_args(c22.function)
    inner = next(factor for factor in outer_factors if isinstance(factor, sp.Integral))
    outer = sp.Mul(*(factor for factor in outer_factors if factor is not inner))
    return outer * inner.function


def test_compact_model_action_is_immutable_non_ast_action_metadata():
    x = sp.Symbol("x")
    action = apply_a22_first_schur_model_compact(sp.Function("f")(x), x)

    assert isinstance(action, FirstSchurCompactModelAction)
    assert not isinstance(action, OperatorAtom)
    assert not isinstance(action, Product)
    assert not isinstance(action, LinearCombination)
    with pytest.raises(FrozenInstanceError):
        action.correction = sp.Integral(1, (x, 0, 1))


def test_compact_action_is_explicitly_associated_with_mod_compact_model():
    x = sp.Symbol("x")
    action = apply_a22_first_schur_model_compact(sp.Function("f")(x), x)

    assert action.relation == a22_first_schur_mod_compact_relation()
    assert isinstance(action.relation.exact, FirstSchurReduction)
    assert action.relation.exact != action.relation.model
    assert not hasattr(action.relation, "rewrite")
    assert not hasattr(action.relation, "as_expr")


def test_compact_action_delegates_to_diagonal_and_combined_kernel_sources(monkeypatch):
    x, y, u, v = sp.symbols("x y u v")
    f = sp.Function("f")
    calls = []
    original_diagonal = schur_module.a22_exact_operator
    original_correction = schur_module.apply_combined_kernel_c22

    def tracked_diagonal():
        calls.append("diagonal")
        return original_diagonal()

    def tracked_correction(*args, **kwargs):
        calls.append("correction")
        return original_correction(*args, **kwargs)

    monkeypatch.setattr(schur_module, "a22_exact_operator", tracked_diagonal)
    monkeypatch.setattr(
        schur_module,
        "apply_combined_kernel_c22",
        tracked_correction,
    )

    result = apply_a22_first_schur_model_compact(
        f(x),
        x,
        input_variable=y,
        outer_variable=u,
        middle_variable=v,
    )

    assert calls == ["diagonal", "correction"]
    assert result.diagonal == apply_linear_combination_ordered(
        a22_exact_operator(),
        f(x),
        x,
        mvp_atomic_rules(),
        integration_variable=y,
    )
    assert result.correction == apply_combined_kernel_c22(
        x,
        f,
        y,
        outer_variable=u,
        middle_variable=v,
        rules=mvp_atomic_rules(),
    )


def test_formal_compact_action_preserves_order_principal_value_and_formal_kernels():
    x = sp.Symbol("x")
    action = apply_a22_first_schur_model_compact(sp.Function("f")(x), x)

    assert isinstance(action.diagonal, AppliedLinearCombination)
    assert len(action.diagonal.terms) == 4
    assert tuple(term.coefficient for term in action.diagonal.terms) == (
        HALF,
        HALF,
        HALF,
        -HALF,
    )
    assert action.diagonal.as_expr().has(PrincipalValue)
    assert isinstance(action.correction, sp.Integral)
    assert action.correction.has(
        sp.Function("Lplus_12"),
        sp.Function("Lminus_21"),
        sp.Function("R11"),
    )


def test_explicit_compact_action_changes_only_wiener_hopf_kernel_realization():
    x = sp.Symbol("x")
    f = sp.Function("f")
    rules = explicit_wiener_hopf_rules(decay=positive_decay_symbol())
    formal = apply_a22_first_schur_model_compact(f(x), x)
    explicit = apply_a22_first_schur_model_compact(f(x), x, rules=rules)

    assert explicit.diagonal == formal.diagonal
    assert not explicit.correction.has(
        sp.Function("Lplus_12"),
        sp.Function("Lminus_21"),
    )
    assert explicit.correction.has(
        sp.Function("chi_infinity"),
        sp.Function("R11"),
    )


def test_compact_action_uses_three_fresh_distinct_internal_variables():
    x = sp.Symbol("x")
    action = apply_a22_first_schur_model_compact(sp.Function("f")(x), x)
    integrals = action.correction.atoms(sp.Integral)
    variables = {
        integral.variables[0]
        for integral in integrals
    }

    assert variables == set(sp.symbols("y u v"))
    assert variables.isdisjoint(action.correction.free_symbols)


def test_compact_action_avoids_capture_when_external_variable_uses_dummy_name():
    y = sp.Symbol("y")
    action = apply_a22_first_schur_model_compact(sp.Function("f")(y), y)
    expression = action.as_expr()
    bound = collect_bound_symbols(expression)

    assert y in expression.free_symbols
    assert y not in bound
    assert bound == set(sp.symbols("u v w"))


def test_explicit_internal_variables_are_respected_and_collisions_fail():
    x, y, u, v = sp.symbols("x y u v")
    f = sp.Function("f")
    action = apply_a22_first_schur_model_compact(
        f(x),
        x,
        input_variable=y,
        outer_variable=u,
        middle_variable=v,
    )

    assert action.correction.variables == [y]
    assert collect_bound_symbols(action.correction) == {y, u, v}
    with pytest.raises(CompactSchurActionError):
        apply_a22_first_schur_model_compact(f(x), x, input_variable=x)
    with pytest.raises(ValueError):
        apply_a22_first_schur_model_compact(
            f(x), x, input_variable=y, outer_variable=y
        )
    with pytest.raises(ValueError):
        apply_a22_first_schur_model_compact(
            f(x),
            x,
            input_variable=y,
            outer_variable=u,
            middle_variable=u,
        )


@pytest.mark.parametrize(
    "operand",
    (
        sp.Symbol("x") + 1,
        sp.Function("f")(sp.Symbol("x")) + sp.Function("g")(sp.Symbol("x")),
        sp.Function("f")(sp.Symbol("x") ** 2),
        sp.sin(sp.Symbol("x")),
    ),
)
def test_compact_action_rejects_operands_outside_function_scope(operand):
    with pytest.raises(CompactSchurActionError):
        apply_a22_first_schur_model_compact(operand, sp.Symbol("x"))


def test_sympy_projection_is_local_sum_without_losing_separated_structure():
    x = sp.Symbol("x")
    action = apply_a22_first_schur_model_compact(sp.Function("f")(x), x)

    assert action.as_expr() == action.diagonal.as_expr() + action.correction
    assert len(action.diagonal.terms) == 4
    assert isinstance(action.correction, sp.Integral)
    assert not isinstance(action.as_expr(), LinearCombination)


@pytest.mark.parametrize("explicit", (False, True))
def test_compact_correction_matches_four_expanded_model_kernels_end_to_end(explicit):
    x, y, u, v = sp.symbols("x y u v")
    f = sp.Function("f")
    rules = (
        explicit_wiener_hopf_rules(decay=positive_decay_symbol())
        if explicit
        else mvp_atomic_rules()
    )
    compact = apply_a22_first_schur_model_compact(
        f(x),
        x,
        input_variable=y,
        outer_variable=u,
        middle_variable=v,
        rules=rules,
    )
    expanded = apply_linear_combination_ordered(
        a22_first_schur_model(),
        f(x),
        x,
        rules,
        integration_variable=y,
    )
    expanded_correction = AppliedLinearCombination(expanded.terms[4:])
    extracted = extract_applied_kernels(expanded_correction, x, f, y)
    factorized_terms = sp.Add.make_args(
        sp.expand(compact_kernel_integrand(compact, f))
    )
    expected_terms = tuple(
        term.coefficient * nested_integrand(term.kernel)
        for term in extracted.terms
    )

    assert compact.diagonal.terms == expanded.terms[:4]
    assert tuple(term.coefficient for term in extracted.terms) == (1, -1, -1, 1)
    assert len(factorized_terms) == len(expected_terms) == 4
    assert all(
        any(sp.simplify(candidate - expected) == 0 for candidate in factorized_terms)
        for expected in expected_terms
    )


def test_no_exact_first_schur_action_or_rewrite_was_added():
    relation = a22_first_schur_mod_compact_relation()

    assert not hasattr(relation.exact, "apply")
    assert not hasattr(relation, "rewrite")
    assert not hasattr(relation, "as_expr")
