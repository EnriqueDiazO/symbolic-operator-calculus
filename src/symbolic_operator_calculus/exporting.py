"""File export for already-rendered symbolic derivations."""

from __future__ import annotations

from os import PathLike
from pathlib import Path

from .rendering import (
    RenderedFirstSchurDerivation,
    RenderedNormalizedFirstSchurPivot,
)


def export_first_schur_derivation_tex(
    rendered_derivation: RenderedFirstSchurDerivation,
    destination: str | PathLike[str],
) -> Path:
    """Write one rendered first-Schur derivation as a minimal LaTeX document."""

    if not isinstance(rendered_derivation, RenderedFirstSchurDerivation):
        raise TypeError(
            "rendered_derivation must be a RenderedFirstSchurDerivation."
        )
    path = Path(destination)
    document = _latex_document(rendered_derivation)
    path.write_text(document, encoding="utf-8")
    return path


def _latex_document(rendered: RenderedFirstSchurDerivation) -> str:
    sections = "\n\n".join(
        rf"\section*{{{step.title}}}" "\n" rf"\[{step.latex}\]"
        for step in rendered.steps
    )
    return (
        "\\documentclass{article}\n"
        "\\usepackage{amsmath}\n"
        "\\usepackage{amssymb}\n"
        "\\usepackage[margin=1in]{geometry}\n"
        "\\begin{document}\n\n"
        f"{sections}\n\n"
        "\\end{document}\n"
    )


def normalized_first_schur_pivot_tex_fragment(
    rendered_derivation: RenderedNormalizedFirstSchurPivot,
) -> str:
    """Build a reusable no-preamble LaTeX fragment from rendered trace data."""

    if not isinstance(rendered_derivation, RenderedNormalizedFirstSchurPivot):
        raise TypeError(
            "rendered_derivation must be a RenderedNormalizedFirstSchurPivot."
        )
    by_key = {step.key: step for step in rendered_derivation.steps}
    required_keys = (
        "reduced_definition",
        "pivot_definition",
        "normalized_expansion",
        "inverse_rule",
        "inverse_substitution",
        "diagonal_recognition",
        "exact_result",
        "mod_compact_result",
        "controlled_expansion",
    )
    missing = tuple(key for key in required_keys if key not in by_key)
    if missing:
        raise ValueError(f"rendered derivation is missing steps: {missing!r}.")

    def display(key: str) -> str:
        return "\\[\n" + by_key[key].latex + "\n\\]"

    obligations = "\n".join(
        rf"\item {statement}"
        for statement in rendered_derivation.proof_obligations
    )
    return (
        "% Reusable fragment generated from symbolic_operator_calculus.\n"
        "% Requires the paper's standard amsmath support; defines no private macros.\n\n"
        "\\subsection*{Definición del primer pivote de Schur normalizado}\n"
        f"{display('reduced_definition')}\n"
        f"{display('pivot_definition')}\n\n"
        "\\subsection*{Normalización ordenada}\n"
        f"{display('normalized_expansion')}\n\n"
        "\\subsection*{Sustitución formal del regularizador}\n"
        f"{display('inverse_rule')}\n"
        f"{display('inverse_substitution')}\n"
        f"{display('diagonal_recognition')}\n\n"
        "\\subsection*{Identidad formal exacta}\n"
        f"{display('exact_result')}\n\n"
        "\\subsection*{Equivalencia módulo compactos}\n"
        f"{display('mod_compact_result')}\n\n"
        "\\subsection*{Expansión controlada}\n"
        f"{display('controlled_expansion')}\n\n"
        "\\subsection*{Obligaciones analíticas pendientes}\n"
        "\\begin{enumerate}\n"
        f"{obligations}\n"
        "\\end{enumerate}\n"
    )


def export_normalized_first_schur_pivot_tex(
    rendered_derivation: RenderedNormalizedFirstSchurPivot,
    destination: str | PathLike[str],
) -> Path:
    """Write the normalized first-pivot derivation as a reusable fragment."""

    path = Path(destination)
    path.write_text(
        normalized_first_schur_pivot_tex_fragment(rendered_derivation),
        encoding="utf-8",
    )
    return path
