"""Text and LaTeX projections for typed scalar Mellin objects."""

from __future__ import annotations

import sympy as sp

from ..conditional import ConditionalIdentity, ConditionalVerificationStatus
from ..rendering import (
    render_assumption_context_latex,
    render_complex_domain_latex,
    render_conditional_identity_latex,
    render_singular_set_latex,
)
from .expressions import MellinExpression
from .symbols import MellinSymbol


def _status_label(status: ConditionalVerificationStatus) -> str:
    labels = {
        ConditionalVerificationStatus.DECLARED: "Formal Mellin expression",
        ConditionalVerificationStatus.EVIDENCE_SUPPLIED: "Evidence supplied",
        ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS: (
            "Symbolically checked under assumptions"
        ),
    }
    return labels[status]


def render_mellin_expression_text(expression: MellinExpression) -> str:
    """Render all semantic fields without suggesting operator semantics."""

    if not isinstance(expression, MellinExpression):
        raise TypeError("expression must be a MellinExpression.")
    roles = ", ".join(
        f"{item.symbol}:{item.role.value}"
        for item in (*expression.variables, *expression.parameters)
    ) or "none"
    assumptions = "; ".join(
        sp.sstr(item) for item in expression.domain.assumption_context.assumptions
    ) or "none"
    evidence = ", ".join(str(item) for item in expression.evidence) or "none"
    return "\n".join(
        (
            "Formal Mellin expression",
            f"Expression: {sp.sstr(expression.as_expr())}",
            f"Roles: {roles}",
            f"Assumptions: {assumptions}",
            f"Domain: {expression.domain.complex_domain}",
            "Singularities: "
            f"{len(expression.domain.singular_set.singularities)}",
            f"Status: {_status_label(expression.verification_status)}",
            f"Evidence: {evidence}",
        )
    )


def render_mellin_symbol_text(symbol: MellinSymbol) -> str:
    """Render a scalar symbol and its dependency classification."""

    if not isinstance(symbol, MellinSymbol):
        raise TypeError("symbol must be a MellinSymbol.")
    return "\n".join(
        (
            "Mellin scalar symbol",
            f"Dependency: {symbol.dependency.value}",
            render_mellin_expression_text(symbol.expression),
        )
    )


def render_mellin_expression_latex(expression: MellinExpression) -> str:
    """Return a presentation projection; strings carry no identity semantics."""

    if not isinstance(expression, MellinExpression):
        raise TypeError("expression must be a MellinExpression.")
    roles = r",\;".join(
        rf"{sp.latex(item.symbol)}:\mathrm{{{item.role.value}}}"
        for item in (*expression.variables, *expression.parameters)
    ) or r"\mathrm{none}"
    evidence = (
        r",\;".join(rf"\text{{{str(item)}}}" for item in expression.evidence)
        or r"\mathrm{none}"
    )
    return (
        r"\begin{aligned}"
        r"&\text{Formal Mellin expression}\\"
        rf"&{sp.latex(expression.as_expr())}\\"
        rf"&\text{{roles: }}{roles}\\"
        rf"&\text{{assumptions: }}"
        rf"{render_assumption_context_latex(expression.domain.assumption_context)}\\"
        rf"&\text{{domain: }}"
        rf"{render_complex_domain_latex(expression.domain.complex_domain)}\\"
        rf"&\text{{singularities: }}"
        rf"{render_singular_set_latex(expression.domain.singular_set)}\\"
        rf"&\text{{status: {_status_label(expression.verification_status)}}}\\"
        rf"&\text{{evidence: }}{evidence}"
        r"\end{aligned}"
    )


def render_mellin_symbol_latex(symbol: MellinSymbol) -> str:
    """Render a scalar Mellin symbol without operator terminology."""

    if not isinstance(symbol, MellinSymbol):
        raise TypeError("symbol must be a MellinSymbol.")
    return (
        r"\begin{aligned}"
        r"&\text{Mellin scalar symbol}\\"
        rf"&\text{{dependency: {symbol.dependency.value}}}\\"
        rf"&{render_mellin_expression_latex(symbol.expression)}"
        r"\end{aligned}"
    )


def render_conditional_scalar_identity_text(
    identity: ConditionalIdentity,
) -> str:
    """Render a P0-B conditional identity with its limited verification label."""

    if not isinstance(identity, ConditionalIdentity):
        raise TypeError("identity must be a ConditionalIdentity.")
    assumptions = "; ".join(
        sp.sstr(item) for item in identity.assumption_context.assumptions
    ) or "none"
    return "\n".join(
        (
            "Conditional scalar identity",
            f"Relation: {sp.sstr(identity.lhs)} = {sp.sstr(identity.rhs)}",
            f"Assumptions: {assumptions}",
            f"Domain: {identity.domain}",
            f"Singularities: {len(identity.singular_set.singularities)}",
            f"Status: {_status_label(identity.verification_status)}",
        )
    )


def render_conditional_scalar_identity_latex(
    identity: ConditionalIdentity,
) -> str:
    """Delegate formula rendering while explicitly retaining scalar wording."""

    if not isinstance(identity, ConditionalIdentity):
        raise TypeError("identity must be a ConditionalIdentity.")
    return (
        r"\begin{aligned}"
        r"&\text{Conditional scalar identity}\\"
        rf"&{render_conditional_identity_latex(identity)}"
        r"\end{aligned}"
    )
