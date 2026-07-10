"""Compact action of the first reduced Schur model for the m=2 case."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import sympy as sp

from .actions import (
    AppliedLinearCombination,
    AtomicAction,
    apply_linear_combination_ordered,
    mvp_atomic_rules,
)
from .blocks import (
    a22_exact_operator,
    a22_first_schur_mod_compact_relation,
)
from .kernels import apply_combined_kernel_c22
from .operators import OperatorAtom
from .relations import ModCompactSchurRelation
from .substitution import collect_bound_symbols, fresh_symbol


class CompactSchurActionError(ValueError):
    """Raised when the compact m=2 model action cannot be constructed."""


@dataclass(frozen=True)
class FirstSchurCompactModelAction:
    """Ordered diagonal action plus one compact Schur-correction integral."""

    relation: ModCompactSchurRelation
    diagonal: AppliedLinearCombination
    correction: sp.Integral

    def __post_init__(self) -> None:
        if not isinstance(self.relation, ModCompactSchurRelation):
            raise TypeError("relation must be a ModCompactSchurRelation.")
        if not isinstance(self.diagonal, AppliedLinearCombination):
            raise TypeError("diagonal must be an AppliedLinearCombination.")
        if not isinstance(self.correction, sp.Integral):
            raise TypeError("correction must be a SymPy Integral.")

    def as_expr(self) -> sp.Expr:
        """Project the separated model action to its final scalar expression."""

        return self.diagonal.as_expr() + self.correction


def apply_a22_first_schur_model_compact(
    operand: sp.Expr,
    variable: sp.Symbol,
    *,
    input_variable: sp.Symbol | None = None,
    outer_variable: sp.Symbol | None = None,
    middle_variable: sp.Symbol | None = None,
    rules: Mapping[OperatorAtom, AtomicAction] | None = None,
) -> FirstSchurCompactModelAction:
    """Apply the modulo-compact first Schur model to an operand ``f(x)``.

    The compact correction infrastructure requires an undefined scalar
    function applied directly to the external variable. General scalar
    operands are intentionally outside the scope of this m=2 API.
    """

    input_function = _input_function_from_operand(operand, variable)
    input_variable = _resolve_input_variable(
        operand,
        variable,
        input_variable,
    )
    resolved_rules = mvp_atomic_rules() if rules is None else rules
    relation = a22_first_schur_mod_compact_relation()
    diagonal = apply_linear_combination_ordered(
        a22_exact_operator(),
        operand,
        variable,
        resolved_rules,
        integration_variable=input_variable,
    )
    correction = apply_combined_kernel_c22(
        variable,
        input_function,
        input_variable,
        outer_variable=outer_variable,
        middle_variable=middle_variable,
        rules=resolved_rules,
    )
    return FirstSchurCompactModelAction(
        relation=relation,
        diagonal=diagonal,
        correction=correction,
    )


def _input_function_from_operand(
    operand: sp.Expr,
    variable: sp.Symbol,
) -> sp.FunctionClass:
    if not isinstance(variable, sp.Symbol):
        raise CompactSchurActionError("variable must be a SymPy Symbol.")
    if (
        not isinstance(operand, sp.core.function.AppliedUndef)
        or len(operand.args) != 1
        or operand.args[0] != variable
        or not isinstance(operand.func, sp.FunctionClass)
    ):
        raise CompactSchurActionError(
            "operand must be an undefined scalar function applied as f(variable)."
        )
    return operand.func


def _resolve_input_variable(
    operand: sp.Expr,
    variable: sp.Symbol,
    input_variable: sp.Symbol | None,
) -> sp.Symbol:
    avoid = (
        operand.free_symbols
        | operand.atoms(sp.Symbol)
        | collect_bound_symbols(operand)
        | {variable}
    )
    if input_variable is None:
        return fresh_symbol(avoid)
    if not isinstance(input_variable, sp.Symbol):
        raise CompactSchurActionError(
            "input_variable must be a SymPy Symbol."
        )
    if input_variable in avoid:
        raise CompactSchurActionError(
            "input_variable collides with the operand or external variable."
        )
    return input_variable
