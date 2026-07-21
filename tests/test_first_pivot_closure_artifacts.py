import json
from pathlib import Path

from symbolic_operator_calculus import (
    build_first_pivot_closure_analysis,
    render_closure_evidence_manifest_yaml,
    render_closure_graph_markdown,
    render_closure_interface_matrix_markdown,
    render_closure_obligations_yaml,
    render_closure_words_latex,
    render_minimal_closure_decision_yaml,
)


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
NOTEBOOK = ROOT / "notebooks" / "first_pivot_integrated_closure.ipynb"


def test_machine_readable_artifacts_are_exact_public_api_outputs():
    analysis = build_first_pivot_closure_analysis()

    assert (DOCS / "FIRST_PIVOT_EVIDENCE_MANIFEST.yaml").read_text() == (
        render_closure_evidence_manifest_yaml(analysis.evidence)
    )
    assert (DOCS / "FIRST_PIVOT_CLOSURE_OBLIGATIONS.yaml").read_text() == (
        render_closure_obligations_yaml(analysis.obligations)
    )


def test_human_reports_embed_generated_matrix_graph_and_decision():
    analysis = build_first_pivot_closure_analysis()
    matrix = (DOCS / "FIRST_PIVOT_RESULT_INTERFACE_MATRIX.md").read_text()
    graph = (DOCS / "FIRST_PIVOT_CLOSURE_GRAPH.md").read_text()
    decision = (DOCS / "FIRST_PIVOT_MINIMAL_LEMMA_DECISION.md").read_text()

    assert render_closure_interface_matrix_markdown(analysis.interface_matrix) in matrix
    assert render_closure_graph_markdown(analysis.obligations) in graph
    assert render_minimal_closure_decision_yaml(analysis.decision) in decision
    assert "decision: NONE" in decision
    assert "confidence: high" in decision


def test_route_and_identification_reports_retain_conservative_outcomes():
    w_audit = (DOCS / "FIRST_PIVOT_W12_PLUS_IDENTIFICATION_AUDIT.md").read_text()
    mellin = (DOCS / "FIRST_PIVOT_MELLIN_ROUTE_ANALYSIS.md").read_text()
    cusp = (DOCS / "FIRST_PIVOT_CUSP_ROUTE_ANALYSIS.md").read_text()

    assert "conclusion: NOT_IDENTIFIED" in w_audit
    assert "relation_to_R_y: NONE_PROVED" in w_audit
    assert "final_symbol_generated: false" in mellin
    assert "conclusion: NOT_SUPPORTED" in cusp
    assert "fredholm_criterion_applied: false" in cusp


def test_latex_export_contains_api_words_and_explicit_status_warning():
    analysis = build_first_pivot_closure_analysis()
    latex = (DOCS / "FIRST_PIVOT_MINIMAL_CLOSURE_LEMMA.tex").read_text()

    assert render_closure_words_latex(analysis.words) in latex
    assert r"\boxed{\text{decision: NONE}" in latex
    assert "does not certify Mellin-PDO membership" in latex
    assert "documentclass" not in latex


def test_notebook_has_fourteen_executed_sections_and_no_error_outputs():
    notebook = json.loads(NOTEBOOK.read_text())
    markdown = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "markdown"
    )
    code_cells = [
        cell for cell in notebook["cells"] if cell["cell_type"] == "code"
    ]

    assert len(code_cells) == 14
    assert [cell["execution_count"] for cell in code_cells] == list(range(1, 15))
    assert all(
        output.get("output_type") != "error"
        for cell in code_cells
        for output in cell["outputs"]
    )
    assert all(f"# {index}." in markdown for index in range(1, 15))
    assert any(
        "decision: NONE" in json.dumps(output, ensure_ascii=False)
        for cell in code_cells
        for output in cell["outputs"]
    )
    assert any(
        "'K++': 4" in json.dumps(output, ensure_ascii=False)
        for cell in code_cells
        for output in cell["outputs"]
    )


def test_notebook_consumes_phase_p_and_phase_o_public_apis():
    notebook = json.loads(NOTEBOOK.read_text())
    code_source = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )

    assert "build_first_pivot_closure_analysis" in code_source
    assert "build_complete_correction_expansion_trace" in code_source
    assert "render_closure_words_markdown" in code_source
    assert "render_closure_interface_matrix_markdown" in code_source
    assert "render_closure_graph_markdown" in code_source
