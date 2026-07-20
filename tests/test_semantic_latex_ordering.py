import ast
from pathlib import Path

import sympy as sp
from semantic_helpers import explicit_r11_kernel_representation

from symbolic_operator_calculus import (
    G1,
    G2,
    KernelCombination,
    build_first_schur_derivation_trace,
    explicit_wiener_hopf_rules,
    m12_kernel,
    m12_kernel_combination,
    m21_kernel,
    m21_kernel_combination,
    positive_decay_symbol,
    render_first_schur_derivation_latex,
    render_operator_atom_latex,
    render_scalar_latex,
)


def _trace(*, explicit=False):
    x, y, u, v = sp.symbols("x y u v")
    rules = (
        explicit_wiener_hopf_rules(decay=positive_decay_symbol())
        if explicit
        else None
    )
    return build_first_schur_derivation_trace(
        sp.Function("f")(x),
        x,
        input_variable=y,
        outer_variable=u,
        middle_variable=v,
        rules=rules,
        regularizer_kernel=explicit_r11_kernel_representation(),
    )


def _step(trace, key):
    return next(
        step.latex
        for step in render_first_schur_derivation_latex(trace).steps
        if step.key == key
    )


def test_operator_atoms_and_applied_scalar_functions_have_distinct_notation():
    x = sp.Symbol("x")

    assert render_operator_atom_latex(G1) == r"G_{1}\,I"
    assert render_operator_atom_latex(G2) == r"G_{2}\,I"
    assert render_scalar_latex(sp.Function("G1")(x)) == r"G_{1}\!\left(x\right)"
    assert render_scalar_latex(sp.Function("G2")(x)) == r"G_{2}\!\left(x\right)"


def test_offdiagonal_and_sign_steps_use_multiplication_operators():
    trace = _trace()
    offdiagonal = _step(trace, "offdiagonal_models")
    sign = _step(trace, "correction_sign")

    assert (
        r"A_{2,1} &\simeq - \left(\widetilde V_{\alpha_2} - G_{2}\,I\right)\,W^-_{2,1}"
        in offdiagonal
    )
    assert (
        r"A_{1,2} &\simeq \left(\widetilde V_{\alpha_1} - G_{1}\,I\right)\,W^+_{1,2}"
        in offdiagonal
    )
    assert r"B_{2,1} &:= \left(\widetilde V_{\alpha_2} - G_{2}\,I\right)" in sign
    assert r"B_{1,2} &:= \left(\widetilde V_{\alpha_1} - G_{1}\,I\right)" in sign
    assert offdiagonal.count(r"\simeq") == 2


def test_operator_correction_wraps_every_g_multiplication_atom():
    correction = _step(_trace(), "correction_expansion")

    assert correction.count(r"\left(G_{1}\,I\right)") == 2
    assert correction.count(r"\left(G_{2}\,I\right)") == 2


def test_ordered_lateral_kernel_apis_are_the_scalar_source_of_truth():
    x, y, u, v = sp.symbols("x y u v")
    m21 = m21_kernel_combination(x, u)
    m12 = m12_kernel_combination(v, y)

    assert isinstance(m21, KernelCombination)
    assert isinstance(m12, KernelCombination)
    assert tuple(term.coefficient for term in m21.terms) == (1, -1)
    assert tuple(term.coefficient for term in m12.terms) == (1, -1)
    assert m21_kernel(x, u) == m21.as_expr()
    assert m12_kernel(v, y) == m12.as_expr()


def test_trace_keeps_ordered_combinations_and_their_scalar_projections():
    trace = _trace()

    assert isinstance(trace.m21_combination, KernelCombination)
    assert isinstance(trace.m12_combination, KernelCombination)
    assert trace.m21 == trace.m21_combination.as_expr()
    assert trace.m12 == trace.m12_combination.as_expr()


def test_formal_lateral_kernels_preserve_term_and_factor_provenance_order():
    trace = _trace()
    left = _step(trace, "left_kernel")
    right = _step(trace, "right_kernel")

    rho2 = left.index(r"\rho_2\!\left(x\right)")
    shifted_lminus = left.index(r"L^-_{2,1}\!\left(\gamma_2 x,u\right)")
    minus_g2 = left.index(r"- G_{2}\!\left(x\right)")
    rho1 = right.index(r"\rho_1\!\left(v\right)")
    shifted_lplus = right.index(r"L^+_{1,2}\!\left(\gamma_1 v,y\right)")
    minus_g1 = right.index(r"- G_{1}\!\left(v\right)")

    assert rho2 < shifted_lminus < minus_g2
    assert rho1 < shifted_lplus < minus_g1
    assert r"G_{1}\,I\!\left" not in right
    assert r"G_{2}\,I\!\left" not in left


def test_explicit_mode_preserves_semantic_order_without_operator_i_in_kernels():
    trace = _trace(explicit=True)
    rendered = render_first_schur_derivation_latex(trace)
    left = _step(trace, "left_kernel")
    right = _step(trace, "right_kernel")

    assert left.index(r"\rho_2\!\left(x\right)") < left.index(r"\chi_\infty")
    assert left.index(r"\chi_\infty") < left.index(r"- G_{2}\!\left(x\right)")
    assert right.index(r"\rho_1\!\left(v\right)") < right.index(r"\chi_\infty")
    assert right.index(r"\chi_\infty") < right.index(r"- G_{1}\!\left(v\right)")
    assert r"G_{1}\,I\!\left" not in left + right
    assert r"G_{2}\,I\!\left" not in left + right
    assert r"G_{1}\,I" in rendered.as_latex()
    assert r"G_{2}\,I" in rendered.as_latex()


def test_renderer_contains_no_string_or_sympy_order_reconstruction():
    path = (
        Path(__file__).parents[1]
        / "src"
        / "symbolic_operator_calculus"
        / "rendering.py"
    )
    tree = ast.parse(path.read_text(encoding="utf-8"))
    calls = [node.func for node in ast.walk(tree) if isinstance(node, ast.Call)]

    assert all(
        not (isinstance(call, ast.Attribute) and call.attr == "replace")
        for call in calls
    )
    assert all(
        not (isinstance(call, ast.Name) and call.id in {"sorted", "default_sort_key", "srepr"})
        for call in calls
    )
    assert all(
        not (
            isinstance(call, ast.Attribute)
            and isinstance(call.value, ast.Name)
            and call.value.id == "sp"
            and call.attr in {"Add", "Mul"}
        )
        for call in calls
    )
    assert all(
        not (isinstance(node, ast.Import) and any(alias.name == "re" for alias in node.names))
        for node in ast.walk(tree)
    )
