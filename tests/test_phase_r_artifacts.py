import json
from pathlib import Path

import sympy as sp

from symbolic_operator_calculus import (
    build_first_pivot_left_core_closure,
    render_first_pivot_left_cores_latex,
    render_left_dilation_covariance_latex,
)


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
NOTEBOOK = ROOT / "notebooks" / "left_cauchy_dilation_closure.ipynb"


def test_covariance_and_scaling_reports_record_exact_formulas_and_boundaries():
    covariance = (DOCS / "LEFT_DILATION_MELLIN_COVARIANCE_DERIVATION.md").read_text()
    scaling = (DOCS / "RADIAL_SCALING_INVARIANCE_OF_MELLIN_SYMBOLS.md").read_text()

    assert r"D_\gamma^{-1}\operatorname{Op}(a)" in covariance
    assert r"d_{\gamma^{-1}}(\lambda)=\gamma^{-i\lambda}" in covariance
    assert r"D_\gamma^{-1}\operatorname{Op}(a)D_\gamma" in covariance
    assert r"\operatorname{cm}^{\infty}_{t/\gamma}(a)" in scaling
    assert "Uniform derivative tails" in scaling
    assert "ANALYTICALLY_PROVED" in scaling


def test_cauchy_audit_has_source_provenance_and_no_projection_falsehood():
    audit = (DOCS / "CAUCHY_FACTOR_MELLIN_SEMIPRODUCT_AUDIT.md").read_text()

    assert "CERTIFIED_MOD_COMPACT" in audit
    assert "Theorem 3.3" in audit
    assert "printed page | 942" in audit
    assert "3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec" in audit
    assert r"P^-P^+=P^+P^-\ne0" in audit
    assert "P^+P^-=0" not in audit
    assert "Theorem 3.5" in audit and "is not used" in audit


def test_consolidated_reports_have_four_relations_and_no_fredholm_claim():
    markdown = (DOCS / "FIRST_PIVOT_LEFT_CORE_CLOSURE.md").read_text()
    latex = (DOCS / "FIRST_PIVOT_LEFT_CORE_CLOSURE.tex").read_text()

    assert markdown.count(r"\simeq") >= 8
    assert "right-factorized" in markdown
    assert "standard `E-tilde`" in markdown
    assert "P-01 and P-05 remain\nblocked" in markdown
    assert "documentclass" not in latex
    assert "No conclusion about" in latex
    assert "is Fredholm" not in latex


def test_public_phase_r_renderers_are_deterministic():
    lam, gamma, x = sp.symbols("lambda gamma x")
    first = build_first_pivot_left_core_closure(lam, gamma, x)
    second = build_first_pivot_left_core_closure(lam, gamma, x)

    assert render_first_pivot_left_cores_latex(first) == (
        render_first_pivot_left_cores_latex(second)
    )
    assert render_left_dilation_covariance_latex(first.covariance) == (
        render_left_dilation_covariance_latex(second.covariance)
    )


def test_phase_r_notebook_has_thirteen_executed_sections_without_errors():
    notebook = json.loads(NOTEBOOK.read_text())
    markdown = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "markdown"
    )
    code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]

    assert len(code_cells) == 13
    assert [cell["execution_count"] for cell in code_cells] == list(range(1, 14))
    assert all(
        output.get("output_type") != "error"
        for cell in code_cells
        for output in cell["outputs"]
    )
    assert all(f"# {number}." in markdown for number in range(1, 14))
    outputs = json.dumps(
        [output for cell in code_cells for output in cell["outputs"]],
        ensure_ascii=False,
    )
    assert "gamma**(-I*lambda)" in outputs
    assert "ANALYTICALLY_PROVED" in outputs
    assert "decision: NONE" in outputs
