import ast
from pathlib import Path

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    build_first_schur_derivation_trace,
    export_first_schur_derivation_tex,
    render_first_schur_derivation_latex,
)


@pytest.fixture(scope="module")
def rendered_derivation():
    x, y, u, v = sp.symbols("x y u v")
    trace = build_first_schur_derivation_trace(
        sp.Function("f")(x),
        x,
        input_variable=y,
        outer_variable=u,
        middle_variable=v,
    )
    return render_first_schur_derivation_latex(trace)


def test_export_writes_a_minimal_document_with_all_steps_in_order(
    tmp_path,
    rendered_derivation,
):
    destination = tmp_path / "first-schur.tex"

    returned = export_first_schur_derivation_tex(
        rendered_derivation,
        destination,
    )
    document = destination.read_text(encoding="utf-8")

    assert returned == destination
    assert document.startswith(r"\documentclass{article}")
    assert r"\usepackage{amsmath}" in document
    assert r"\usepackage{amssymb}" in document
    assert r"\usepackage[margin=1in]{geometry}" in document
    assert document.endswith("\\end{document}\n")
    assert document.count(r"\section*") == 9
    positions = [document.index(step.latex) for step in rendered_derivation.steps]
    assert positions == sorted(positions)
    for step in rendered_derivation.steps:
        assert rf"\section*{{{step.title}}}" in document
    assert "Lminus_21" not in document
    assert "Lplus_12" not in document
    assert r"L^-_{2,1}" in document
    assert r"L^+_{1,2}" in document


def test_export_is_deterministic_and_rejects_unrendered_input(
    tmp_path,
    rendered_derivation,
):
    first = tmp_path / "first.tex"
    second = tmp_path / "second.tex"

    export_first_schur_derivation_tex(rendered_derivation, first)
    export_first_schur_derivation_tex(rendered_derivation, second)

    assert first.read_bytes() == second.read_bytes()
    with pytest.raises(TypeError):
        export_first_schur_derivation_tex(object(), tmp_path / "bad.tex")


def test_exporter_only_consumes_rendered_structures_and_writes_text():
    path = (
        Path(__file__).parents[1]
        / "src"
        / "symbolic_operator_calculus"
        / "exporting.py"
    )
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imports = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    }
    calls = [node.func for node in ast.walk(tree) if isinstance(node, ast.Call)]

    assert imports.isdisjoint({"kernels", "blocks", "derivations"})
    assert all(
        not (isinstance(call, ast.Name) and call.id in {"open", "run", "Popen"})
        for call in calls
    )
    assert all(
        not (
            isinstance(call, ast.Attribute)
            and call.attr in {"render", "render_first_schur_derivation_latex"}
        )
        for call in calls
    )
