import json
from pathlib import Path

from symbolic_operator_calculus import (
    build_complete_correction_expansion_trace,
    render_complete_correction_expansion_latex,
    render_correction_interface_matrix_markdown,
    render_correction_term_catalog_latex,
    render_correction_term_catalog_markdown,
)


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
NOTEBOOK = ROOT / "notebooks" / "first_pivot_complete_correction.ipynb"


def test_catalog_artifacts_are_generated_from_the_public_api():
    trace = build_complete_correction_expansion_trace()
    markdown = (DOCS / "FIRST_PIVOT_CORRECTION_TERM_CATALOG.md").read_text()
    latex = (DOCS / "FIRST_PIVOT_CORRECTION_TERM_CATALOG.tex").read_text()

    assert render_correction_term_catalog_markdown(trace) in markdown
    assert render_correction_term_catalog_latex(trace) in latex
    assert markdown.count("C2-T") == 16
    assert latex.count("C2-T") == 16
    assert "metadatos declarativos, no pruebas" in markdown


def test_interface_reference_and_route_reports_remain_conservative():
    trace = build_complete_correction_expansion_trace()
    interfaces = (DOCS / "FIRST_PIVOT_CORRECTION_INTERFACES.md").read_text()
    references = (DOCS / "FIRST_PIVOT_CORRECTION_REFERENCE_MAP.md").read_text()
    decision = (DOCS / "FIRST_PIVOT_CORRECTION_ROUTE_DECISION.md").read_text()

    assert render_correction_interface_matrix_markdown(trace) in interfaces
    assert "unproved_adaptation" in references
    assert "source verification pending" in references
    assert "Karlovich2017HasemanSO" in references
    assert "**C: insufficient evidence**" in decision
    assert "A: supported by verified closure rules" in decision
    assert "no está sustentada" in decision
    assert "B: supported by verified algebra membership" in decision


def test_complete_latex_contains_all_controlled_api_expansions():
    trace = build_complete_correction_expansion_trace()
    latex = (DOCS / "FIRST_PIVOT_COMPLETE_CORRECTION.tex").read_text()

    for level in ("C0", "C1", "C2", "C3"):
        assert render_complete_correction_expansion_latex(trace, level) in latex
    assert latex.count(r"\mathrm{C2-T") == 16
    assert r"N_{2}^{(1)}\simeq N_{2}+\mathcal C_{2}" in latex
    assert r"\boxed{\text{C: insufficient evidence}}" in latex
    assert "no certificada" in latex


def test_notebook_is_fully_executed_and_contains_no_error_output():
    notebook = json.loads(NOTEBOOK.read_text())
    code_cells = [
        cell for cell in notebook["cells"] if cell["cell_type"] == "code"
    ]
    markdown = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "markdown"
    )

    assert len(code_cells) == 16
    assert [cell["execution_count"] for cell in code_cells] == list(range(1, 17))
    assert all(
        output.get("output_type") != "error"
        for cell in code_cells
        for output in cell["outputs"]
    )
    assert "# 1. Propósito y límites" in markdown
    assert "## 16. Frontera de capacidad" in markdown
    assert any(
        "C2-T16" in json.dumps(output, ensure_ascii=False)
        for cell in code_cells
        for output in cell["outputs"]
    )
