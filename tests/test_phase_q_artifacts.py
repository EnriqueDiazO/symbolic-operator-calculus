import json
from pathlib import Path

from symbolic_operator_calculus import (
    build_first_pivot_closure_analysis,
    render_closure_graph_markdown,
    render_closure_interface_matrix_markdown,
    render_closure_obligations_yaml,
    render_minimal_closure_decision_yaml,
)


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
NOTEBOOK = ROOT / "notebooks" / "right_dilation_mellin_composition.ipynb"


def test_phase_q_classification_contains_nonmembership_and_damping_proofs():
    text = (DOCS / "OSCILLATORY_DILATION_SYMBOL_CLASSIFICATION.md").read_text()

    assert r"d_\gamma\notin V(\mathbb R)" in text
    assert r"2T|\log\gamma|" in text
    assert r"\operatorname{Var}(d_\gamma h)" in text
    assert r"\operatorname{Var}(h)" in text
    assert "not declared to be an\nalgebra" in text
    assert "false\ninclusion" in text


def test_coefficient_audit_retains_mod_compact_strength_and_provenance():
    text = (DOCS / "REGULARIZER_COEFFICIENT_SEMIPRODUCT_AUDIT.md").read_text()

    assert "CERTIFIED_MOD_COMPACT" in text
    assert "Theorem 3.3" in text
    assert "printed 942; PDF 8" in text
    assert "3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec" in text
    assert "not exact equality" in text


def test_phase_p_artifacts_embed_the_updated_public_outputs():
    analysis = build_first_pivot_closure_analysis()

    assert (DOCS / "FIRST_PIVOT_CLOSURE_OBLIGATIONS.yaml").read_text() == (
        render_closure_obligations_yaml(analysis.obligations)
    )
    assert render_closure_graph_markdown(analysis.obligations) in (
        DOCS / "FIRST_PIVOT_CLOSURE_GRAPH.md"
    ).read_text()
    assert render_closure_interface_matrix_markdown(analysis.interface_matrix) in (
        DOCS / "FIRST_PIVOT_RESULT_INTERFACE_MATRIX.md"
    ).read_text()
    assert render_minimal_closure_decision_yaml(analysis.decision) in (
        DOCS / "FIRST_PIVOT_MINIMAL_LEMMA_DECISION.md"
    ).read_text()


def test_phase_q_notebook_has_thirteen_executed_sections_without_errors():
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
    output_text = json.dumps(notebook["cells"], ensure_ascii=False)
    assert "gamma**(I*lambda)" in output_text
    assert "CERTIFIED_MOD_COMPACT" in output_text
    assert "decision: NONE" in output_text
