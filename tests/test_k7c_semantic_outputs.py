import json
from pathlib import Path

import sympy as sp

from symbolic_operator_calculus import (
    a22_first_schur_correction,
    apply_combined_kernel_c22,
    apply_linear_combination_ordered,
    build_first_schur_derivation_trace,
    combined_kernel_c22,
    mvp_atomic_rules,
    render_combined_kernel_action_latex,
    render_first_schur_derivation_latex,
    render_scalar_latex,
)


INTERNAL_PRESENTATION_NAMES = (
    "Lminus_21",
    "Lplus_12",
    "Wminus_21",
    "Wplus_12",
    "Vtilde_alpha1",
    "Vtilde_alpha2",
    "chi_infinity",
    "rho1",
    "rho2",
    "gamma1",
    "gamma2",
    "R11",
)


def _objects():
    x, y, u, v = sp.symbols("x y u v")
    f = sp.Function("f")
    rules = mvp_atomic_rules()
    trace = build_first_schur_derivation_trace(
        f(x),
        x,
        input_variable=y,
        outer_variable=u,
        middle_variable=v,
        rules=rules,
    )
    rendered = render_first_schur_derivation_latex(trace)
    return x, y, u, v, f, rules, trace, rendered


def test_all_rendered_derivation_output_uses_the_central_semantic_notation():
    *_, rendered = _objects()
    latex = rendered.as_latex()

    assert all(name not in latex for name in INTERNAL_PRESENTATION_NAMES)
    assert r"L^-_{2,1}" in latex
    assert r"L^+_{1,2}" in latex
    assert r"W^-_{2,1}" in latex
    assert r"W^+_{1,2}" in latex
    assert r"R_{1,1}" in latex
    assert r"G_{1}\,I" in latex and r"G_{2}\,I" in latex
    assert r"G_{1}\,I\!\left" not in latex
    assert r"G_{2}\,I\!\left" not in latex


def test_four_applied_terms_are_rendered_semantically_without_duplicate_signs():
    x, _, _, _, f, rules, *_ = _objects()
    applied = apply_linear_combination_ordered(
        a22_first_schur_correction(),
        f(x),
        x,
        rules,
    )
    visible = tuple(render_scalar_latex(term.expression) for term in applied.terms)

    assert len(applied.terms) == 4
    assert tuple(term.coefficient for term in applied.terms) == (1, -1, -1, 1)
    assert r"L^-_{2,1}" in visible[0] and r"L^+_{1,2}" in visible[0]
    assert r"G_{1}\!\left(v\right)" in visible[1]
    assert r"G_{2}\!\left(x\right)" in visible[2]
    assert r"G_{2}\!\left(x\right)" in visible[3]
    assert r"G_{1}\!\left(v\right)" in visible[3]
    assert all(not expression.lstrip().startswith("-") for expression in visible)
    assert all(
        "Lminus_21" not in expression and "Lplus_12" not in expression
        for expression in visible
    )


def test_combined_kernel_and_action_come_from_public_apis_and_render_semantically():
    x, y, u, v, f, rules, trace, rendered = _objects()
    kernel = combined_kernel_c22(
        x,
        y,
        outer_variable=u,
        middle_variable=v,
        rules=rules,
    )
    action = apply_combined_kernel_c22(
        x,
        f,
        y,
        outer_variable=u,
        middle_variable=v,
        rules=rules,
    )
    combined_step = next(
        step.latex for step in rendered.steps if step.key == "combined_kernel"
    )
    action_latex = render_combined_kernel_action_latex(action, x, f, y)

    assert kernel == trace.combined_kernel
    assert r"M_{2,1}" in combined_step
    assert r"R_{1,1}" in combined_step
    assert r"M_{1,2}" in combined_step
    assert r"\rho_2\!\left(x\right)" in combined_step
    assert r"\rho_1\!\left(v\right)" in combined_step
    assert render_scalar_latex(action) in action_latex
    assert r"C_{2,2}^{(1)}" in action_latex
    assert render_scalar_latex(f(y)) in action_latex
    assert r"L^-_{2,1}" in action_latex and r"L^+_{1,2}" in action_latex
    assert all(name not in action_latex for name in INTERNAL_PRESENTATION_NAMES)


def test_notebook_structurally_consumes_public_pipeline_and_semantic_renderers():
    path = Path(__file__).parents[1] / "notebooks" / "test.ipynb"
    notebook = json.loads(path.read_text(encoding="utf-8"))
    code_cells = [
        cell for cell in notebook["cells"] if cell["cell_type"] == "code"
    ]
    source = "\n".join("".join(cell["source"]) for cell in code_cells)

    assert path.exists()
    for public_api in (
        "a22_first_schur_correction",
        "apply_linear_combination_ordered",
        "combined_kernel_c22",
        "apply_combined_kernel_c22",
        "a21_mod_compact_relation",
        "a12_mod_compact_relation",
        "m21_kernel_combination",
        "m12_kernel_combination",
        "explicit_wiener_hopf_rules",
        "render_first_schur_derivation_latex",
        "render_scalar_latex",
        "export_first_schur_derivation_tex",
    ):
        assert f"soc.{public_api}" in source
    assert "display(term.expression)" not in source
    assert "display(combined_kernel)" not in source
    assert "display(kernel_action)" not in source
    assert "sp.Integral(" not in source
    assert "soc.Product(" not in source
    assert "soc.Term(" not in source
    assert "display(Math(" in source
    assert all(cell["outputs"] == [] for cell in code_cells)
