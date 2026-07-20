"""Typed infrastructure for scalar Mellin symbols; no operators are defined."""

from .domains import (
    MellinDomainError,
    MellinDomainRoleError,
    MellinFrequencyMismatchError,
    MellinSymbolDomain,
)
from .examples import (
    DiagonalMellinSymbolFamily,
    build_diagonal_mellin_symbols,
    build_dilation_multiplier,
    build_exponential_sinh_symbol,
    build_relative_frequency_symbol,
    build_space_frequency_symbol,
)
from .expressions import (
    MellinBranchConditionRequiredError,
    MellinExpression,
    MellinExpressionError,
    MellinSingularMetadataRequiredError,
    UndeclaredMellinVariableError,
)
from .rendering import (
    render_conditional_scalar_identity_latex,
    render_conditional_scalar_identity_text,
    render_mellin_expression_latex,
    render_mellin_expression_text,
    render_mellin_symbol_latex,
    render_mellin_symbol_text,
)
from .symbols import (
    MellinDependencyError,
    MellinSymbol,
    MellinSymbolDependency,
    MellinSymbolError,
    classify_mellin_dependency,
)
from .variables import (
    MellinVariable,
    MellinVariableError,
    MellinVariableRole,
    MellinVariableRoleError,
    input_spatial_variable,
    mellin_frequency,
    mellin_parameter,
    output_spatial_variable,
    relative_multiplicative_variable,
)

__all__ = [
    "MellinDomainError",
    "MellinDomainRoleError",
    "MellinBranchConditionRequiredError",
    "MellinDependencyError",
    "MellinExpression",
    "MellinExpressionError",
    "MellinFrequencyMismatchError",
    "MellinSingularMetadataRequiredError",
    "MellinSymbol",
    "MellinSymbolDependency",
    "MellinSymbolDomain",
    "MellinSymbolError",
    "MellinVariable",
    "MellinVariableError",
    "MellinVariableRole",
    "MellinVariableRoleError",
    "UndeclaredMellinVariableError",
    "DiagonalMellinSymbolFamily",
    "build_diagonal_mellin_symbols",
    "build_dilation_multiplier",
    "build_exponential_sinh_symbol",
    "build_relative_frequency_symbol",
    "build_space_frequency_symbol",
    "classify_mellin_dependency",
    "input_spatial_variable",
    "mellin_frequency",
    "mellin_parameter",
    "output_spatial_variable",
    "relative_multiplicative_variable",
    "render_conditional_scalar_identity_latex",
    "render_conditional_scalar_identity_text",
    "render_mellin_expression_latex",
    "render_mellin_expression_text",
    "render_mellin_symbol_latex",
    "render_mellin_symbol_text",
]
