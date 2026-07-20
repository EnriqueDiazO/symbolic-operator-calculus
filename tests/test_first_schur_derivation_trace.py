import ast
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest
import sympy as sp

import symbolic_operator_calculus.derivations as derivation_module
from symbolic_operator_calculus import (
    R11,
    DerivationTraceError,
    FirstSchurCorrectionSignTrace,
    FirstSchurDerivationTrace,
    LinearCombination,
    ModCompactSchurRelation,
    OperatorAtom,
    Product,
    Term,
    a12_mod_compact_relation,
    a21_mod_compact_relation,
    a22_first_schur_correction,
    a22_first_schur_model,
    a22_first_schur_mod_compact_relation,
    a22_first_schur_reduction,
    apply_a22_first_schur_model_compact,
    build_first_schur_derivation_trace,
    combined_kernel_c22,
    explicit_wiener_hopf_rules,
    factor_first_schur_correction,
    kminus_kernel,
    kplus_kernel,
    m12_kernel,
    m21_kernel,
    positive_decay_symbol,
)
from semantic_helpers import explicit_r11_kernel_representation


def trace_variables():
    return sp.symbols("x y u v")


def build_trace(*, rules=None):
    x, y, u, v = trace_variables()
    return build_first_schur_derivation_trace(
        sp.Function("f")(x),
        x,
        input_variable=y,
        outer_variable=u,
        middle_variable=v,
        rules=rules,
        regularizer_kernel=explicit_r11_kernel_representation(),
    )


def reconstruct(factorization):
    return LinearCombination(
        tuple(
            Term(
                left.coefficient * right.coefficient,
                Product(
                    left.product.factors
                    + (factorization.regularizer,)
                    + right.product.factors
                ),
            )
            for left in factorization.left.terms
            for right in factorization.right.terms
        )
    )


def test_first_schur_derivation_trace_is_immutable_non_ast_metadata():
    trace = build_trace()

    assert isinstance(trace, FirstSchurDerivationTrace)
    assert not isinstance(trace, OperatorAtom)
    assert not isinstance(trace, Product)
    assert not isinstance(trace, LinearCombination)
    assert not isinstance(trace, ModCompactSchurRelation)
    with pytest.raises(FrozenInstanceError):
        trace.m21 = sp.Integer(0)


def test_builder_uses_existing_public_derivation_apis(monkeypatch):
    called = set()
    names = (
        "a21_mod_compact_relation",
        "a12_mod_compact_relation",
        "a22_first_schur_reduction",
        "a22_first_schur_mod_compact_relation",
        "a22_first_schur_correction",
        "factor_first_schur_correction",
        "m21_kernel",
        "m12_kernel",
        "combined_kernel_c22",
        "apply_a22_first_schur_model_compact",
    )

    for name in names:
        original = getattr(derivation_module, name)

        def tracked(*args, _name=name, _original=original, **kwargs):
            called.add(_name)
            return _original(*args, **kwargs)

        monkeypatch.setattr(derivation_module, name, tracked)

    trace = build_trace()

    assert called == set(names)
    assert isinstance(trace, FirstSchurDerivationTrace)


def test_offdiagonal_relations_are_public_a21_then_a12_objects():
    trace = build_trace()

    assert trace.offdiagonal_relations == (
        a21_mod_compact_relation(),
        a12_mod_compact_relation(),
    )
    assert trace.offdiagonal_relations[0].exact.row == 2
    assert trace.offdiagonal_relations[1].exact.row == 1


def test_sign_trace_derives_all_values_from_current_relations():
    trace = build_trace()
    reduction = a22_first_schur_reduction()
    a21_model = a21_mod_compact_relation().model

    assert isinstance(trace.sign_trace, FirstSchurCorrectionSignTrace)
    assert trace.sign_trace.schur_sign == reduction.offdiagonal_product_coefficient
    assert (
        trace.sign_trace.leading_a21_sign
        == a21_model.expression.terms[0].coefficient
    )
    assert trace.sign_trace.correction_sign == (
        trace.sign_trace.schur_sign * trace.sign_trace.leading_a21_sign
    ) == 1
    with pytest.raises(DerivationTraceError):
        FirstSchurCorrectionSignTrace(-1, -1, -1)


def test_reduced_relation_correction_and_factorization_reuse_public_objects():
    trace = build_trace()

    assert trace.reduced_relation == a22_first_schur_mod_compact_relation()
    assert trace.reduced_relation.exact == a22_first_schur_reduction()
    assert trace.reduced_relation.model == a22_first_schur_model()
    assert trace.correction == a22_first_schur_correction()
    assert trace.factorization == factor_first_schur_correction(trace.correction)
    assert trace.factorization.regularizer is R11
    assert reconstruct(trace.factorization) == trace.correction


def test_kernels_and_compact_action_match_existing_public_apis():
    x, y, u, v = trace_variables()
    f = sp.Function("f")
    trace = build_trace()

    assert trace.m21 == m21_kernel(x, u)
    assert trace.m12 == m12_kernel(v, y)
    assert trace.combined_kernel == combined_kernel_c22(
        x,
        y,
        outer_variable=u,
        middle_variable=v,
        regularizer_kernel=explicit_r11_kernel_representation(),
    )
    assert trace.compact_action == apply_a22_first_schur_model_compact(
        f(x),
        x,
        input_variable=y,
        outer_variable=u,
        middle_variable=v,
        regularizer_kernel=explicit_r11_kernel_representation(),
    )
    assert trace.compact_action.relation == trace.reduced_relation


def test_trace_preserves_all_explicit_deterministic_variables():
    trace = build_trace()
    variables = trace_variables()

    assert (
        trace.output_variable,
        trace.input_variable,
        trace.outer_variable,
        trace.middle_variable,
    ) == variables
    assert len(set(variables)) == 4
    assert trace.combined_kernel.variables == [trace.outer_variable]
    assert trace.compact_action.correction.variables == [trace.input_variable]


@pytest.mark.parametrize(
    ("input_name", "outer_name", "middle_name"),
    (("x", "u", "v"), ("y", "x", "v"), ("y", "u", "x"), ("y", "u", "u")),
)
def test_trace_rejects_colliding_explicit_variables(
    input_name,
    outer_name,
    middle_name,
):
    symbols = {name: sp.Symbol(name) for name in "xyuv"}
    x = symbols["x"]

    with pytest.raises(DerivationTraceError):
        build_first_schur_derivation_trace(
            sp.Function("f")(x),
            x,
            input_variable=symbols[input_name],
            outer_variable=symbols[outer_name],
            middle_variable=symbols[middle_name],
        )


def test_formal_trace_contains_formal_wiener_hopf_and_regularizer_kernels():
    trace = build_trace()

    assert trace.m21.has(sp.Function("Lminus_21"))
    assert trace.m12.has(sp.Function("Lplus_12"))
    assert trace.combined_kernel.has(
        sp.Function("Lminus_21"),
        sp.Function("Lplus_12"),
        sp.Function("R11"),
    )


def test_explicit_trace_uses_normalized_fourier_kernels_and_formal_r11():
    x, y, u, v = trace_variables()
    d = positive_decay_symbol()
    rules = explicit_wiener_hopf_rules(decay=d)
    trace = build_trace(rules=rules)

    assert isinstance(trace, FirstSchurDerivationTrace)
    assert not trace.combined_kernel.has(
        sp.Function("Lminus_21"),
        sp.Function("Lplus_12"),
    )
    assert trace.m21.has(
        kminus_kernel(sp.Symbol("gamma2") * x - u, d)
    )
    assert trace.m12.has(
        kplus_kernel(sp.Symbol("gamma1") * v - y, d)
    )
    assert trace.combined_kernel.has(
        sp.Function("chi_infinity"),
        sp.Function("R11"),
    )


def test_derivation_layer_has_no_renderers_or_manual_kernel_formulas():
    trace = build_trace()
    source_path = (
        Path(__file__).parents[1]
        / "src/symbolic_operator_calculus/derivations.py"
    )
    module = ast.parse(source_path.read_text(encoding="utf-8"))
    function_names = {
        node.name
        for node in ast.walk(module)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    assigned_names = {
        target.id
        for node in ast.walk(module)
        if isinstance(node, (ast.Assign, ast.AnnAssign))
        for target in (
            node.targets if isinstance(node, ast.Assign) else (node.target,)
        )
        if isinstance(target, ast.Name)
    }
    forbidden_names = {
        "latex",
        "to_latex",
        "render_latex",
        "markdown",
        "render_markdown",
        "explain",
        "pretty_print",
        "narrative",
    }
    manual_formula_names = {
        "rho1",
        "rho2",
        "G1_scalar",
        "G2_scalar",
        "gamma1",
        "gamma2",
        "Lplus_12_kernel",
        "Lminus_21_kernel",
        "R11_kernel",
    }

    assert function_names.isdisjoint(forbidden_names)
    assert assigned_names.isdisjoint(manual_formula_names)
    assert all(not hasattr(trace, name) for name in forbidden_names)
