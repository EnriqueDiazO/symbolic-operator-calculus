import ast
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
from pathlib import Path

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    G1,
    G2,
    I,
    R11,
    S_Rplus,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    FirstSchurCorrectionSignTrace,
    LatexRenderingError,
    Product,
    RenderedDerivationStep,
    RenderedFirstSchurDerivation,
    Term,
    a22_exact_operator,
    build_first_schur_derivation_trace,
    explicit_wiener_hopf_rules,
    positive_decay_symbol,
    render_first_schur_derivation_latex,
    render_kernel_combination_latex,
    render_linear_combination_latex,
    render_operator_atom_latex,
    render_product_latex,
    render_scalar_latex,
    render_term_latex,
)


EXPECTED_KEYS = (
    "exact_definition",
    "offdiagonal_models",
    "correction_sign",
    "reduced_model",
    "correction_expansion",
    "left_kernel",
    "right_kernel",
    "combined_kernel",
    "compact_model_action",
)


@pytest.fixture(scope="module")
def formal_trace():
    x, y, u, v = sp.symbols("x y u v")
    return build_first_schur_derivation_trace(
        sp.Function("f")(x),
        x,
        input_variable=y,
        outer_variable=u,
        middle_variable=v,
    )


@pytest.fixture(scope="module")
def explicit_trace():
    x, y, u, v = sp.symbols("x y u v")
    return build_first_schur_derivation_trace(
        sp.Function("f")(x),
        x,
        input_variable=y,
        outer_variable=u,
        middle_variable=v,
        rules=explicit_wiener_hopf_rules(decay=positive_decay_symbol()),
    )


def step_by_key(rendered, key):
    return next(step for step in rendered.steps if step.key == key)


def test_rendered_step_and_derivation_are_immutable_structures():
    step = RenderedDerivationStep("key", "Title", "x = y")
    rendered = RenderedFirstSchurDerivation((step,))

    assert isinstance(step, RenderedDerivationStep)
    assert isinstance(rendered, RenderedFirstSchurDerivation)
    with pytest.raises(FrozenInstanceError):
        step.latex = "z"
    with pytest.raises(FrozenInstanceError):
        rendered.steps = ()
    with pytest.raises(LatexRenderingError):
        RenderedDerivationStep("", "Title", "x")


def test_renderer_accepts_only_trace_and_preserves_stable_step_order(formal_trace):
    rendered = render_first_schur_derivation_latex(formal_trace)

    assert tuple(step.key for step in rendered.steps) == EXPECTED_KEYS
    assert rendered.as_latex() == "\n\n".join(
        step.latex for step in rendered.steps
    )
    with pytest.raises(TypeError):
        render_first_schur_derivation_latex(formal_trace.reduced_relation.exact)
    with pytest.raises(TypeError):
        render_first_schur_derivation_latex(formal_trace.correction)


@pytest.mark.parametrize(
    ("atom", "expected"),
    (
        (I, r"I"),
        (R11, r"R_{1,1}"),
        (S_Rplus, r"S_{\mathbb R_+}"),
        (G1, r"G_{1}\,I"),
        (G2, r"G_{2}\,I"),
        (Vtilde_alpha1, r"\widetilde V_{\alpha_1}"),
        (Vtilde_alpha2, r"\widetilde V_{\alpha_2}"),
        (Wplus_12, r"W^+_{1,2}"),
        (Wminus_21, r"W^-_{2,1}"),
    ),
)
def test_operator_atoms_use_central_explicit_notation(atom, expected):
    assert render_operator_atom_latex(atom) == expected


def test_product_and_linear_combination_preserve_ast_order_and_exact_halves():
    product = Vtilde_alpha2 * Wminus_21 * R11 * Vtilde_alpha1 * Wplus_12
    combination = a22_exact_operator()
    negative_half = Term(Fraction(-1, 2), Product((G2, S_Rplus)))

    assert render_product_latex(product) == (
        r"\widetilde V_{\alpha_2}\,W^-_{2,1}\,R_{1,1}\,"
        r"\widetilde V_{\alpha_1}\,W^+_{1,2}"
    )
    rendered = render_linear_combination_latex(combination)
    assert rendered.index(r"\widetilde V_{\alpha_2}\,I") < rendered.index(
        r"\widetilde V_{\alpha_2}\,S_{\mathbb R_+}"
    ) < rendered.index(r"\left(G_{2}\,I\right)\,I") < rendered.index(
        r"\left(G_{2}\,I\right)\,S_{\mathbb R_+}"
    )
    assert r"\frac{1}{2}" in rendered
    assert "0.5" not in rendered
    assert render_term_latex(negative_half).startswith(r"- \frac{1}{2}")


def test_scalar_printer_controls_known_names_and_principal_value():
    x, y = sp.symbols("x y")
    gamma1, gamma2 = sp.symbols("gamma1 gamma2")
    expression = (
        sp.Function("Lplus_12")(gamma1 * x, y)
        + sp.Function("Lminus_21")(gamma2 * x, y)
        + sp.Function("chi_infinity")(x)
        + sp.Function("rho1")(x)
        + sp.Function("rho2")(x)
        + sp.Function("R11")(x, y)
    )
    principal_value = sp.Function("f")(x) + sp.Function("G1")(x)
    from symbolic_operator_calculus import PrincipalValue

    principal_value += PrincipalValue(
        sp.Integral(sp.Function("f")(y) / (y - x), (y, 0, sp.oo))
    )
    latex = render_scalar_latex(expression)
    pv_latex = render_scalar_latex(principal_value)

    assert r"L^+_{1,2}" in latex
    assert r"L^-_{2,1}" in latex
    assert r"\chi_\infty" in latex
    assert r"\rho_1" in latex and r"\rho_2" in latex
    assert r"\gamma_1" in latex and r"\gamma_2" in latex
    assert r"R_{1,1}" in latex
    assert r"\operatorname{p.v.}" in pv_latex
    assert "PrincipalValue" not in pv_latex


def test_exact_and_mod_compact_steps_use_distinct_relation_symbols(formal_trace):
    rendered = render_first_schur_derivation_latex(formal_trace)
    exact = step_by_key(rendered, "exact_definition").latex
    offdiagonal = step_by_key(rendered, "offdiagonal_models").latex
    reduced = step_by_key(rendered, "reduced_model").latex
    compact = step_by_key(rendered, "compact_model_action").latex

    assert " = " in exact and r"\simeq" not in exact
    assert offdiagonal.count(r"\simeq") == 2
    assert r"\simeq" in reduced
    assert r"\simeq" in compact
    assert r"A_{2,2}^{(1)} = A_{2,2}^{(1),\mathrm{model}}" not in compact


def test_each_step_renders_the_corresponding_trace_field(formal_trace):
    rendered = render_first_schur_derivation_latex(formal_trace)

    assert r"A_{2,2}^{(1)}" in step_by_key(rendered, "exact_definition").latex
    assert step_by_key(rendered, "offdiagonal_models").latex.count(
        r"\simeq"
    ) == len(formal_trace.offdiagonal_relations)
    sign_step = step_by_key(rendered, "correction_sign").latex
    assert sign_step.count("-") >= 2
    assert "+B_{2,1}" in sign_step
    correction_step = step_by_key(rendered, "correction_expansion").latex
    assert all(
        render_product_latex(term.product) in correction_step
        for term in formal_trace.correction.terms
    )
    assert render_kernel_combination_latex(formal_trace.m21_combination) in step_by_key(
        rendered, "left_kernel"
    ).latex
    assert render_kernel_combination_latex(formal_trace.m12_combination) in step_by_key(
        rendered, "right_kernel"
    ).latex
    combined_step = step_by_key(rendered, "combined_kernel").latex
    assert r"M_{2,1}" in combined_step and r"M_{1,2}" in combined_step
    compact_step = step_by_key(rendered, "compact_model_action").latex
    assert all(
        render_scalar_latex(term.expression) in compact_step
        for term in formal_trace.compact_action.diagonal.terms
    )
    assert render_scalar_latex(formal_trace.compact_action.correction) not in compact_step


def test_sign_step_changes_when_trace_sign_metadata_changes(formal_trace):
    changed = replace(
        formal_trace,
        sign_trace=FirstSchurCorrectionSignTrace(1, 1, 1),
    )

    original_latex = step_by_key(
        render_first_schur_derivation_latex(formal_trace),
        "correction_sign",
    ).latex
    changed_latex = step_by_key(
        render_first_schur_derivation_latex(changed),
        "correction_sign",
    ).latex

    assert original_latex != changed_latex
    assert r"+\left[+B_{2,1}\right]" in changed_latex


def test_formal_output_uses_mathematical_not_internal_names(formal_trace):
    latex = render_first_schur_derivation_latex(formal_trace).as_latex()
    forbidden = (
        "Lplus_12",
        "Lminus_21",
        "Vtilde_alpha1",
        "Vtilde_alpha2",
        "Wplus_12",
        "Wminus_21",
        "S_Rplus",
        "chi_infinity",
        "rho1",
        "rho2",
        "gamma1",
        "gamma2",
        "PrincipalValue",
    )

    assert r"L^+_{1,2}" in latex and r"L^-_{2,1}" in latex
    assert r"\widetilde V_{\alpha_1}" in latex
    assert r"\widetilde V_{\alpha_2}" in latex
    assert r"W^+_{1,2}" in latex and r"W^-_{2,1}" in latex
    assert r"R_{1,1}" in latex
    assert r"\operatorname{p.v.}" in latex
    assert all(name not in latex for name in forbidden)


def test_explicit_trace_uses_same_renderer_without_formal_l_names(
    formal_trace,
    explicit_trace,
):
    formal_rendered = render_first_schur_derivation_latex(formal_trace)
    explicit_rendered = render_first_schur_derivation_latex(explicit_trace)
    explicit_latex = explicit_rendered.as_latex()

    assert type(formal_rendered) is type(explicit_rendered)
    assert tuple(step.key for step in explicit_rendered.steps) == EXPECTED_KEYS
    assert r"L^+_{1,2}" not in explicit_latex
    assert r"L^-_{2,1}" not in explicit_latex
    assert r"\chi_\infty" in explicit_latex
    assert r"R_{1,1}" in explicit_latex


def test_renderer_is_deterministic_and_does_not_recalculate_trace(formal_trace):
    first = render_first_schur_derivation_latex(formal_trace)
    second = render_first_schur_derivation_latex(formal_trace)

    assert first == second
    assert first.as_latex() == second.as_latex()


def test_rendering_layer_has_no_global_replacement_regex_ast_add_or_file_output():
    root = Path(__file__).parents[1] / "src/symbolic_operator_calculus"
    rendering_path = root / "rendering.py"
    module = ast.parse(rendering_path.read_text(encoding="utf-8"))
    imported_modules = {
        alias.name
        for node in ast.walk(module)
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    }
    imported_from_modules = {
        node.module
        for node in ast.walk(module)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    }
    called_attributes = {
        node.func.attr
        for node in ast.walk(module)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    }
    called_names = {
        node.func.id
        for node in ast.walk(module)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }

    assert "re" not in imported_modules
    assert imported_from_modules.isdisjoint(
        {
            "symbolic_operator_calculus.blocks",
            "symbolic_operator_calculus.schur",
            "blocks",
            "schur",
        }
    )
    assert "replace" not in called_attributes
    assert "Add" not in called_attributes
    assert called_names.isdisjoint({"open", "write", "markdown", "render_markdown"})
    for path in root.glob("*.py"):
        if path.name in {"__init__.py", "rendering.py", "exporting.py"}:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        assert all(
            not (
                isinstance(node, ast.ImportFrom)
                and node.module is not None
                and node.module.endswith("rendering")
            )
            for node in ast.walk(tree)
        )
