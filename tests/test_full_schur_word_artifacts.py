from __future__ import annotations

import json
from pathlib import Path

import symbolic_operator_calculus as soc


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "full_schur_word_factorization.ipynb"


def test_full_word_api_is_public_and_has_no_duplicate_exports():
    expected = {
        "AuxiliaryOperator",
        "CauchyProjectionOperator",
        "CoefficientMultiplicationOperator",
        "CutoffEquivalenceClaim",
        "FullSchurFactorizationStatus",
        "FullSchurWordFactorizationResult",
        "FullWordBlock",
        "LocalizedOperator",
        "RawSchurCorrection",
        "SchurPivot",
        "SignedSchurContribution",
        "SpatialCutoffOperator",
        "TransportedMellinCore",
        "TransportedShiftOperator",
        "conjugate_cutoff_by_dilation",
        "factor_full_first_schur_word",
    }

    assert expected.issubset(set(soc.__all__))
    assert len(soc.__all__) == len(set(soc.__all__))
    for name in expected:
        assert getattr(soc, name) is not None


def test_first_commit_documents_exist_and_keep_the_conditional_verdicts():
    expected = {
        "PHASE_P1B_FULL_SCHUR_WORD_SOURCE_AUDIT.md": (
            "BLOCKED_BY_REGULARIZER_INTERFACE",
            "T_{1,-}^2",
        ),
        "PHASE_P1B_LEFT_CORE_FACTORIZATION.md": ("BLOCKED", "P^+"),
        "PHASE_P1B_RIGHT_CORE_FACTORIZATION.md": ("BLOCKED", "Z_1^{-1}Z_1^{-1}"),
        "PHASE_P1B_CUTOFF_CONTROL.md": (
            "BLOCKED_BY_CUTOFF_CONTROL",
            "status=ASSUMED|CERTIFIED|UNPROVED",
        ),
        "PHASE_P1B_SCHUR_SIGN_CONVENTION.md": (
            "RawSchurCorrection",
            "SignedSchurContribution",
        ),
    }

    for filename, needles in expected.items():
        text = (ROOT / "docs" / filename).read_text(encoding="utf-8")
        for needle in needles:
            assert needle in text


def test_final_math_documents_keep_blocked_article_boundary_and_obligations():
    math = (ROOT / "docs" / "PHASE_P1B_FULL_SCHUR_WORD_MATH.md").read_text()
    review = (
        ROOT / "docs" / "PHASE_P1B_FULL_SCHUR_WORD_ADVERSARIAL_REVIEW.md"
    ).read_text()
    latex = (ROOT / "docs" / "PHASE_P1B_FULL_SCHUR_WORD_LATEX.md").read_text()
    obligations = (
        "left_core_factorization",
        "right_core_factorization",
        "cutoff_replacement_mod_compacts",
        "regularizer_mellin_representation",
        "wh_mellin_wh_sandwich_membership",
        "fredholm_algebra_membership",
        "schur_correction_symbol",
    )

    assert "BLOCKED_BY_REGULARIZER_INTERFACE" in math
    assert "synthetic tests validate the inference engine" in math.lower()
    assert all(key in math for key in obligations)
    assert review.count("`PASS`") >= 12
    assert "| `FAIL` |" not in review
    assert "The pivot is not declared solved" in review
    assert latex.startswith("# conditional draft — not ready for insertion")
    assert "not a proved\narticle lemma" in latex
    assert all(key in latex for key in obligations)


def test_full_word_notebook_has_eleven_executed_sections_and_no_errors():
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    markdown = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "markdown"
    )
    code_cells = [
        cell for cell in notebook["cells"] if cell["cell_type"] == "code"
    ]

    assert len(code_cells) == 11
    assert [cell["execution_count"] for cell in code_cells] == list(range(1, 12))
    assert all(
        output.get("output_type") != "error"
        for cell in code_cells
        for output in cell["outputs"]
    )
    assert all(f"# {number}." in markdown for number in range(1, 12))


def test_full_word_notebook_uses_public_api_and_stores_both_outcomes():
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    code = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )
    outputs = json.dumps(
        [
            output
            for cell in notebook["cells"]
            if cell["cell_type"] == "code"
            for output in cell["outputs"]
        ]
    )

    assert "import symbolic_operator_calculus as soc" in code
    assert "soc.factor_full_first_schur_word(" in code
    assert "soc.SpatialCutoffOperator(" in code
    assert "soc.conjugate_cutoff_by_dilation(" in code
    assert "soc.RegularizerInterfaceRole.LEFT_AUXILIARY" in code
    assert "FULL_WORD_EXACT_FACTORIZATION" in outputs
    assert "BLOCKED_BY_REGULARIZER_INTERFACE" in outputs
    assert "Phi_delta_inverse Op_r11 Phi_delta" in outputs
    assert "chi_infinity(gamma*x)" in outputs
    assert "schur_correction_symbol" in outputs
