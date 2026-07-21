from pathlib import Path

from symbolic_operator_calculus import (
    build_normalized_first_schur_pivot_derivation,
    export_normalized_first_schur_pivot_tex,
    normalized_first_schur_pivot_tex_fragment,
    render_normalized_first_schur_pivot_latex,
)


def rendered_derivation():
    trace = build_normalized_first_schur_pivot_derivation()
    return render_normalized_first_schur_pivot_latex(trace)


def test_fragment_contains_required_material_without_a_document_preamble():
    fragment = normalized_first_schur_pivot_tex_fragment(rendered_derivation())

    assert r"\documentclass" not in fragment
    assert r"\begin{document}" not in fragment
    assert r"\end{document}" not in fragment
    assert "newcommand" not in fragment
    assert r"A_{2,2}^{(1)} = A_{2,2}" in fragment
    assert r"N_{2}^{(1)} = Z_{2}^{-1}" in fragment
    assert r"A_{1,1}^{(-1)}" in fragment
    assert (
        r"N_{2}^{(1)} = N_{2} - Z_{2}^{-1}\,A_{2,1}\,T_{1,-}"
        in fragment
    )
    assert (
        r"N_{2}^{(1)} \simeq N_{2} + Z_{2}^{-1}\,\left("
        r"\widetilde V_{\alpha_2} - G_{2}\right)"
        in fragment
    )
    assert fragment.count(r"\item ") == 5


def test_fragment_is_deterministic_and_export_writes_it_verbatim(tmp_path):
    rendered = rendered_derivation()
    first_fragment = normalized_first_schur_pivot_tex_fragment(rendered)
    second_fragment = normalized_first_schur_pivot_tex_fragment(rendered)
    destination = tmp_path / "first_schur_pivot.tex"

    returned = export_normalized_first_schur_pivot_tex(rendered, destination)

    assert first_fragment == second_fragment
    assert returned == destination
    assert destination.read_text(encoding="utf-8") == first_fragment


def test_repository_fragment_is_exactly_the_public_api_output():
    destination = (
        Path(__file__).parents[1]
        / "docs"
        / "FIRST_SCHUR_PIVOT_DERIVATION.tex"
    )

    assert destination.read_text(encoding="utf-8") == (
        normalized_first_schur_pivot_tex_fragment(rendered_derivation())
    )
