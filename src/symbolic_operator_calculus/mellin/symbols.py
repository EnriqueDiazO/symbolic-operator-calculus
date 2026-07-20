"""Scalar Mellin symbols, deliberately separate from expressions and operators."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import final

import sympy as sp

from ..conditional import ConditionalIdentity, ConditionalVerificationStatus
from ..semantics import ApproximateEquality, ExactIdentityScope
from ..singularities import detect_singularities
from .domains import MellinSymbolDomain
from .expressions import MellinExpression
from .variables import MellinVariableRole


class MellinSymbolDependency(str, Enum):
    """Supported scalar dependencies, classified from explicit variable roles."""

    FREQUENCY_ONLY = "frequency_only"
    SPACE_FREQUENCY = "space_frequency"
    RELATIVE_FREQUENCY = "relative_frequency"
    PARAMETRIC_CONSTANT = "parametric_constant"
    CONSTANT = "constant"


class MellinSymbolError(ValueError):
    """Raised when a scalar Mellin symbol is malformed or misclassified."""


class MellinDependencyError(MellinSymbolError):
    """Raised when free-variable roles do not match a supported dependency."""


def classify_mellin_dependency(
    expression: MellinExpression,
) -> MellinSymbolDependency:
    """Classify only from free symbols and their explicit roles, never names."""

    if not isinstance(expression, MellinExpression):
        raise TypeError("expression must be a MellinExpression.")
    free = expression.free_symbols
    used_variables = tuple(
        item for item in expression.variables if item.symbol in free
    )
    used_parameters = tuple(
        item for item in expression.parameters if item.symbol in free
    )
    roles = {item.role for item in used_variables}
    frequency = MellinVariableRole.FREQUENCY in roles
    spatial = bool(
        roles
        & {
            MellinVariableRole.OUTPUT_SPATIAL,
            MellinVariableRole.INPUT_SPATIAL,
        }
    )
    relative = MellinVariableRole.RELATIVE_MULTIPLICATIVE in roles
    if not roles:
        return (
            MellinSymbolDependency.PARAMETRIC_CONSTANT
            if used_parameters
            else MellinSymbolDependency.CONSTANT
        )
    if frequency and not spatial and not relative and len(roles) == 1:
        return MellinSymbolDependency.FREQUENCY_ONLY
    if frequency and spatial and not relative:
        return MellinSymbolDependency.SPACE_FREQUENCY
    if frequency and relative and not spatial:
        return MellinSymbolDependency.RELATIVE_FREQUENCY
    raise MellinDependencyError(
        "free-variable roles do not match a supported scalar Mellin dependency."
    )


@final
@dataclass(frozen=True)
class MellinSymbol:
    """A typed scalar symbol; it carries no operator or analytic semantics."""

    expression: MellinExpression
    dependency: MellinSymbolDependency
    description: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.expression, MellinExpression):
            raise TypeError(
                "expression must be a MellinExpression with an explicit domain."
            )
        if not isinstance(self.dependency, MellinSymbolDependency):
            raise TypeError("dependency must be a MellinSymbolDependency.")
        classified = classify_mellin_dependency(self.expression)
        if self.dependency is not classified:
            raise MellinDependencyError(
                f"declared dependency {self.dependency.value} conflicts with "
                f"role-based classification {classified.value}."
            )
        if self.description is not None and not self.description.strip():
            raise MellinSymbolError("description must be nonempty when supplied.")

    @classmethod
    def constant(
        cls,
        expression: sp.Expr,
        *,
        domain: MellinSymbolDomain,
        verification_status: ConditionalVerificationStatus = (
            ConditionalVerificationStatus.DECLARED
        ),
        evidence: object | tuple[object, ...] | None = (),
        description: str | None = None,
    ) -> MellinSymbol:
        """Construct a strict scalar constant without inventing frequency use."""

        if not isinstance(expression, sp.Expr):
            raise TypeError("constant expression must be a SymPy expression.")
        if expression.free_symbols:
            raise MellinDependencyError(
                "a strict Mellin constant cannot contain free symbols."
            )
        formal = MellinExpression(
            expression=expression,
            domain=domain,
            verification_status=verification_status,
            evidence=evidence,
            description=description,
        )
        return cls(
            expression=formal,
            dependency=MellinSymbolDependency.CONSTANT,
            description=description,
        )

    @property
    def domain(self) -> MellinSymbolDomain:
        return self.expression.domain

    @property
    def frequency(self):
        return self.domain.frequency

    @property
    def variables(self):
        return self.expression.variables

    @property
    def spatial_variables(self):
        return tuple(
            item
            for item in self.variables
            if item.role
            in {
                MellinVariableRole.OUTPUT_SPATIAL,
                MellinVariableRole.INPUT_SPATIAL,
            }
        )

    @property
    def relative_variable(self):
        return next(
            (
                item
                for item in self.variables
                if item.role is MellinVariableRole.RELATIVE_MULTIPLICATIVE
            ),
            None,
        )

    @property
    def parameters(self):
        return self.expression.parameters

    @property
    def assumptions(self):
        return self.domain.assumption_context

    @property
    def singularities(self):
        return self.domain.singular_set

    @property
    def verification_status(self):
        return self.expression.verification_status

    @property
    def evidence(self):
        return self.expression.evidence

    def structurally_equals(self, other: MellinSymbol) -> bool:
        """Compare the complete typed objects, not mathematical truth."""

        return isinstance(other, MellinSymbol) and self == other

    def sympy_equals(self, other: MellinSymbol) -> bool | None:
        """Ask SymPy about a scalar residual without changing semantic types."""

        if not isinstance(other, MellinSymbol):
            raise TypeError("other must be a MellinSymbol.")
        residual = sp.simplify(self.as_expr() - other.as_expr())
        return residual.is_zero

    def conditional_identity(
        self,
        other: MellinSymbol,
        *,
        description: str | None = None,
    ) -> ConditionalIdentity:
        """Build a P0-B conditional scalar identity, never an ExactIdentity."""

        if not isinstance(other, MellinSymbol):
            raise TypeError("other must be a MellinSymbol.")
        domain = self.domain.intersect(other.domain)
        detected = detect_singularities(
            self.as_expr() - other.as_expr(),
            domain.frequency.symbol,
            assumptions=domain.assumption_context,
            domain=domain.complex_domain,
        )
        domain = domain.with_singularities(detected)
        return ConditionalIdentity(
            lhs=self.as_expr(),
            rhs=other.as_expr(),
            assumption_context=domain.assumption_context,
            domain=domain.complex_domain,
            singular_set=domain.singular_set,
            scope=ExactIdentityScope.SCALAR_SYMBOLIC,
            verification_status=ConditionalVerificationStatus.DECLARED,
            evidence=(*self.evidence, *other.evidence),
            description=description,
        )

    def approximate_equality(
        self,
        other: MellinSymbol,
        *,
        tolerance: sp.Expr,
        norm: object,
    ) -> ApproximateEquality:
        """Create the separate P0-A numerical relation type."""

        if not isinstance(other, MellinSymbol):
            raise TypeError("other must be a MellinSymbol.")
        return ApproximateEquality(
            left=self,
            right=other,
            tolerance=tolerance,
            norm=norm,
            residual=self.as_expr() - other.as_expr(),
        )

    def as_expr(self) -> sp.Expr:
        """Explicitly project the scalar expression, losing all metadata."""

        return self.expression.as_expr()

    def __str__(self) -> str:
        return (
            f"Mellin scalar symbol({sp.sstr(self.as_expr())}; "
            f"dependency={self.dependency.value}; "
            f"status={self.verification_status.value})"
        )
