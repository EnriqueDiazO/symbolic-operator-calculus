import ast
from pathlib import Path

import sympy as sp

from symbolic_operator_calculus import (
    build_first_schur_derivation_trace,
    explicit_wiener_hopf_rules,
    positive_decay_symbol,
    render_first_schur_derivation_latex,
    render_product_latex,
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
    )


def _step(trace, key):
    rendered = render_first_schur_derivation_latex(trace)
    return next(step.latex for step in rendered.steps if step.key == key)


def test_correction_expansion_uses_one_ast_term_per_semantic_line():
    trace = _trace()
    latex = _step(trace, "correction_expansion")
    products = tuple(
        render_product_latex(term.product) for term in trace.correction.terms
    )

    assert latex.startswith(r"\begin{aligned}")
    assert latex.endswith(r"\end{aligned}")
    assert tuple(term.coefficient for term in trace.correction.terms) == (1, -1, -1, 1)
    assert all(latex.count(product) == 1 for product in products)
    assert [latex.index(product) for product in products] == sorted(
        latex.index(product) for product in products
    )
    assert latex.count(r" \\ ") == 3
    assert latex.count(r"&{}- ") == 2
    assert latex.count(r"&{}+ ") == 1
    assert r"\left(G_{1}\,I\right)" in latex
    assert r"\left(G_{2}\,I\right)" in latex


def test_combined_kernel_is_compact_first_then_structurally_expanded():
    trace = _trace()
    latex = _step(trace, "combined_kernel")

    compact_m21 = latex.index(r"M_{2,1}")
    compact_r11 = latex.index(r"R_{1,1}")
    compact_m12 = latex.index(r"M_{1,2}")
    expanded_rho2 = latex.index(r"\rho_2\!\left(x\right)")
    expanded_g2 = latex.index(r"G_{2}\!\left(x\right)")
    expanded_rho1 = latex.index(r"\rho_1\!\left(v\right)")
    expanded_g1 = latex.index(r"G_{1}\!\left(v\right)")

    assert latex.startswith(r"\begin{aligned}")
    assert compact_m21 < compact_r11 < compact_m12 < expanded_rho2
    assert expanded_rho2 < expanded_g2 < expanded_rho1 < expanded_g1
    assert r"\times\Biggl[" in latex
    assert r"L^-_{2,1}" in latex and r"L^+_{1,2}" in latex
    assert latex.count(r"C_{2,2}^{(1)}") == 2


def test_explicit_combined_kernel_uses_same_layout_and_normalized_factors():
    latex = _step(_trace(explicit=True), "combined_kernel")

    assert r"L^-_{2,1}" not in latex and r"L^+_{1,2}" not in latex
    assert r"\chi_\infty" in latex
    assert r"2 \pi d" in latex
    assert r"R_{1,1}" in latex
    assert latex.index(r"\rho_2") < latex.index(r"G_{2}")
    assert latex.index(r"\rho_1") < latex.index(r"G_{1}")


def test_compact_action_references_c22_once_and_splits_real_diagonal_terms():
    trace = _trace()
    latex = _step(trace, "compact_model_action")
    diagonal_expressions = tuple(
        render_scalar_latex(term.expression)
        for term in trace.compact_action.diagonal.terms
    )

    assert latex.startswith(r"\begin{gathered}\begin{aligned}")
    assert latex.index(r"A_{2,2}^{(1)}") < latex.index(
        r"\left(A_{2,2}^{(1),\mathrm{model}} f\right)"
    )
    assert r"\left(A_{2,2}f\right)" in latex
    assert latex.count(r"C_{2,2}^{(1)}") == 1
    assert r"\int_0^\infty C_{2,2}^{(1)}" in latex
    assert render_scalar_latex(trace.compact_action.correction) not in latex
    assert len(trace.compact_action.diagonal.terms) == 4
    assert all(latex.count(expression) == 1 for expression in diagonal_expressions)
    assert latex.count(r" \\ ") == 5
    assert r"\\[0.8em]" in latex
    assert r"\operatorname{p.v.}" in latex
    assert r"G_{2}\!\left(x\right)" in latex
    assert r"G_{2}\,I\!\left(x\right)" not in latex


def test_layout_is_structural_and_exporter_has_no_mathematical_layout_policy():
    root = Path(__file__).parents[1] / "src" / "symbolic_operator_calculus"
    rendering_tree = ast.parse((root / "rendering.py").read_text(encoding="utf-8"))
    exporting_source = (root / "exporting.py").read_text(encoding="utf-8")
    calls = [
        node.func for node in ast.walk(rendering_tree) if isinstance(node, ast.Call)
    ]

    assert all(
        not (isinstance(call, ast.Attribute) and call.attr == "replace")
        for call in calls
    )
    assert all(
        not (
            isinstance(call, ast.Name)
            and call.id in {"sorted", "default_sort_key", "wrap", "fill"}
        )
        for call in calls
    )
    rendering_source = (root / "rendering.py").read_text(encoding="utf-8")
    for forbidden in (r"\resizebox", r"\scalebox", r"\tiny", r"\scriptsize"):
        assert forbidden not in rendering_source
    assert "if step.key" not in exporting_source
    assert ".kernels" not in exporting_source
    assert ".blocks" not in exporting_source
