import json
from pathlib import Path


NOTEBOOK_PATH = (
    Path(__file__).parents[1]
    / "notebooks"
    / "first_schur_pivot_normalization.ipynb"
)


def load_notebook():
    return json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))


def test_notebook_has_all_fourteen_required_sections_in_order():
    notebook = load_notebook()
    markdown = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "markdown"
    )
    section_titles = (
        "1. Propósito y alcance",
        "2. Convenciones de igualdad",
        "3. Definición de átomos no conmutativos",
        "4. Construcción de $A_{2,2}^{(1)}$",
        "5. Normalización por $Z_2^{-1}$ y $T_{2,-}$",
        "6. Sustitución de $A_{1,1}^{(-1)}$",
        "7. Reconocimiento del término $N_2$",
        "8. Derivación factorizada de $N_2^{(1)}$",
        "9. Sustitución módulo compactos de $A_{2,1}$ y $A_{1,2}$",
        "10. Expansión controlada de las diferencias lineales",
        "11. Tabla de clasificación de factores",
        "12. Obligaciones analíticas pendientes",
        "13. Frontera de capacidad",
        "14. Exportación LaTeX",
    )

    positions = tuple(markdown.index(title) for title in section_titles)
    assert positions == tuple(sorted(positions))


def test_notebook_consumes_public_normalized_pivot_api():
    notebook = load_notebook()
    code = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )

    assert "import symbolic_operator_calculus as soc" in code
    assert "soc.build_normalized_first_schur_pivot_derivation()" in code
    assert "soc.render_normalized_first_schur_pivot_exact_latex(trace)" in code
    assert "soc.render_normalized_first_schur_pivot_mod_compact_latex(" in code
    assert "soc.render_normalized_factor_classification_markdown(trace)" in code
    assert "soc.render_normalized_proof_obligations_markdown(trace)" in code
    assert "sympy" not in code.lower()


def test_notebook_is_fully_executed_without_error_outputs():
    notebook = load_notebook()
    code_cells = [
        cell for cell in notebook["cells"] if cell["cell_type"] == "code"
    ]
    outputs = [output for cell in code_cells for output in cell["outputs"]]

    assert code_cells
    assert all(isinstance(cell["execution_count"], int) for cell in code_cells)
    assert all(output["output_type"] != "error" for output in outputs)


def test_notebook_stores_both_required_main_results():
    notebook = load_notebook()
    outputs = json.dumps(
        [
            output
            for cell in notebook["cells"]
            if cell["cell_type"] == "code"
            for output in cell["outputs"]
        ],
        ensure_ascii=False,
    )

    assert (
        r"N_{2}^{(1)} = N_{2} - Z_{2}^{-1}\\,A_{2,1}\\,T_{1,-}"
        in outputs
    )
    assert (
        r"N_{2}^{(1)} \\simeq N_{2} + Z_{2}^{-1}\\,\\left("
        r"\\widetilde V_{\\alpha_2} - G_{2}\\right)"
        in outputs
    )
    assert "MellinPDORegularizer" in outputs
    assert "Justificar la sustitución Wiener--Hopf" in outputs
