"""File export for already-rendered symbolic derivations."""

from __future__ import annotations

from os import PathLike
from pathlib import Path

from .rendering import RenderedFirstSchurDerivation


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
