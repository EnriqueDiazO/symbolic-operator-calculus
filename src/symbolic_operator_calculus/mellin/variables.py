"""Explicit variable roles for scalar Mellin notation."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import final

import sympy as sp

from ..domains import AssumptionContext


class MellinVariableRole(str, Enum):
    """A declared mathematical role; names never determine this value."""

    OUTPUT_SPATIAL = "output_spatial"
    INPUT_SPATIAL = "input_spatial"
    RELATIVE_MULTIPLICATIVE = "relative_multiplicative"
    FREQUENCY = "frequency"
    PARAMETER = "parameter"


class MellinVariableError(ValueError):
    """Raised when a Mellin variable declaration is malformed."""


class MellinVariableRoleError(MellinVariableError):
    """Raised when a variable is used with an incompatible declared role."""


@final
@dataclass(frozen=True)
class MellinVariable:
    """A SymPy symbol paired with one explicit role and optional hypotheses."""

    symbol: sp.Symbol
    role: MellinVariableRole
    assumption_context: AssumptionContext = AssumptionContext()
    description: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.symbol, sp.Symbol):
            raise TypeError("symbol must be an existing SymPy Symbol.")
        if not isinstance(self.role, MellinVariableRole):
            raise TypeError("role must be a MellinVariableRole.")
        if not isinstance(self.assumption_context, AssumptionContext):
            raise TypeError("assumption_context must be an AssumptionContext.")
        if self.description is not None and not self.description.strip():
            raise MellinVariableError("description must be nonempty when supplied.")

    @property
    def assumptions(self) -> tuple[sp.Basic, ...]:
        """Return recorded hypotheses without altering SymPy assumptions."""

        return self.assumption_context.assumptions

    def with_assumptions(self, *assumptions: object) -> MellinVariable:
        """Return the same declaration narrowed by more explicit hypotheses."""

        additional = (
            assumptions[0]
            if len(assumptions) == 1
            and isinstance(assumptions[0], AssumptionContext)
            else AssumptionContext(tuple(assumptions))
        )
        return MellinVariable(
            symbol=self.symbol,
            role=self.role,
            assumption_context=self.assumption_context.combine(additional),
            description=self.description,
        )


def _declare(
    symbol: sp.Symbol,
    role: MellinVariableRole,
    *,
    assumption_context: AssumptionContext | None = None,
    description: str | None = None,
) -> MellinVariable:
    if not isinstance(symbol, sp.Symbol):
        raise TypeError("Mellin variable factories require an existing SymPy Symbol.")
    return MellinVariable(
        symbol=symbol,
        role=role,
        assumption_context=(
            AssumptionContext()
            if assumption_context is None
            else assumption_context
        ),
        description=description,
    )


def mellin_frequency(
    symbol: sp.Symbol,
    *,
    assumption_context: AssumptionContext | None = None,
    description: str | None = None,
) -> MellinVariable:
    return _declare(
        symbol,
        MellinVariableRole.FREQUENCY,
        assumption_context=assumption_context,
        description=description,
    )


def output_spatial_variable(
    symbol: sp.Symbol,
    *,
    assumption_context: AssumptionContext | None = None,
    description: str | None = None,
) -> MellinVariable:
    return _declare(
        symbol,
        MellinVariableRole.OUTPUT_SPATIAL,
        assumption_context=assumption_context,
        description=description,
    )


def input_spatial_variable(
    symbol: sp.Symbol,
    *,
    assumption_context: AssumptionContext | None = None,
    description: str | None = None,
) -> MellinVariable:
    return _declare(
        symbol,
        MellinVariableRole.INPUT_SPATIAL,
        assumption_context=assumption_context,
        description=description,
    )


def relative_multiplicative_variable(
    symbol: sp.Symbol,
    *,
    assumption_context: AssumptionContext | None = None,
    description: str | None = None,
) -> MellinVariable:
    return _declare(
        symbol,
        MellinVariableRole.RELATIVE_MULTIPLICATIVE,
        assumption_context=assumption_context,
        description=description,
    )


def mellin_parameter(
    symbol: sp.Symbol,
    *,
    assumption_context: AssumptionContext | None = None,
    description: str | None = None,
) -> MellinVariable:
    return _declare(
        symbol,
        MellinVariableRole.PARAMETER,
        assumption_context=assumption_context,
        description=description,
    )
