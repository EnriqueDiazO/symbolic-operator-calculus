"""Explicit semantic relation and kernel-representation types.

The classes in this module are deliberately small.  They record what a caller
is claiming and which evidence, hypotheses, or comparison criterion was
supplied.  They do not prove compactness, convergence, boundedness, or any
other analytic property.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import final

import sympy as sp


class CertificationStatus(str, Enum):
    """Whether externally supplied evidence accompanies a semantic claim."""

    UNCERTIFIED = "uncertified"
    EVIDENCE_SUPPLIED = "evidence_supplied"


class ExactIdentityScope(str, Enum):
    """Scope in which a caller explicitly claims that an identity is exact."""

    STRUCTURAL = "structural"
    SCALAR_SYMBOLIC = "scalar_symbolic"
    WITHIN_MODEL = "within_model"
    OPERATORIAL = "operatorial"


class KernelRepresentationStatus(str, Enum):
    """Caller-assigned semantic strength of a kernel representation."""

    FORMAL = "formal"
    ASSUMED = "assumed"
    EXTERNALLY_CERTIFIED = "externally_certified"


class ExactIdentityRequiredError(TypeError):
    """Raised when a non-exact semantic relation is used as an exact one."""


class KernelRepresentationRequiredError(ValueError):
    """Raised when a formal regularizer is projected to a kernel implicitly."""


@final
@dataclass(frozen=True)
class ExactIdentity:
    """An exactness claim with explicit scope and caller-supplied evidence.

    Construction records the claim; it does not verify that the evidence is a
    proof or promote a model identity to an operatorial theorem.
    """

    left: object
    right: object
    evidence: object
    scope: ExactIdentityScope
    hypotheses: tuple[object, ...] = ()

    def __post_init__(self) -> None:
        if self.evidence is None:
            raise ValueError("ExactIdentity requires explicit evidence.")
        if not isinstance(self.scope, ExactIdentityScope):
            raise TypeError("scope must be an ExactIdentityScope.")
        if isinstance(self.hypotheses, (str, bytes)):
            raise TypeError("hypotheses must be an iterable of hypothesis objects.")
        try:
            hypotheses = tuple(self.hypotheses)
        except TypeError as exc:
            raise TypeError("hypotheses must be an iterable.") from exc
        object.__setattr__(self, "hypotheses", hypotheses)


@final
@dataclass(frozen=True)
class FormalIdentity:
    """A formal identity that makes no claim of analytic or exact validity."""

    left: object
    right: object
    justification: object | None = None


@final
@dataclass(frozen=True)
class ModCompactEquivalence:
    """A modulo-compact declaration with optional external evidence.

    Omitting ``evidence`` is valid, but the relation then remains explicitly
    uncertified.  This class never infers that ``residual`` is compact.
    """

    left: object
    right: object
    space: object | None = None
    compact_ideal: object | None = None
    residual: object | None = None
    evidence: object | None = None
    certification_status: CertificationStatus = field(init=False)

    def __post_init__(self) -> None:
        status = (
            CertificationStatus.EVIDENCE_SUPPLIED
            if self.evidence is not None
            else CertificationStatus.UNCERTIFIED
        )
        object.__setattr__(self, "certification_status", status)


@final
@dataclass(frozen=True)
class ApproximateEquality:
    """A numerical comparison with an explicit tolerance, norm, and residual."""

    left: object
    right: object
    tolerance: sp.Expr
    norm: object
    residual: object

    def __post_init__(self) -> None:
        if isinstance(self.tolerance, bool):
            raise TypeError("tolerance must be a positive finite real scalar.")
        tolerance = sp.sympify(self.tolerance)
        if (
            tolerance.is_real is not True
            or tolerance.is_positive is not True
            or tolerance.is_finite is not True
        ):
            raise ValueError("tolerance must be provably positive, finite, and real.")
        if self.norm is None:
            raise ValueError("ApproximateEquality requires an explicit norm or criterion.")
        if self.residual is None:
            raise ValueError("ApproximateEquality requires residual information.")
        object.__setattr__(self, "tolerance", tolerance)


def require_exact_identity(relation: object) -> ExactIdentity:
    """Return an exact identity or reject every other semantic relation type."""

    if not isinstance(relation, ExactIdentity):
        raise ExactIdentityRequiredError(
            "This operation requires ExactIdentity; conditional, formal, "
            "modulo-compact, and approximate relations are not interchangeable "
            "with it."
        )
    return relation


@dataclass(frozen=True)
class RegularizerOperator:
    """A regularizer operator with no implicit kernel or inverse semantics."""

    target: object
    operator: object

    def __post_init__(self) -> None:
        if self.target is None:
            raise ValueError("RegularizerOperator.target must be supplied.")
        if self.operator is None:
            raise ValueError("RegularizerOperator.operator must be supplied.")


@dataclass(frozen=True)
class IntegrationDomain:
    """Explicit integration interval; no convergence claim is implied."""

    lower: sp.Expr
    upper: sp.Expr
    description: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "lower", sp.sympify(self.lower))
        object.__setattr__(self, "upper", sp.sympify(self.upper))
        if self.description is not None and not self.description.strip():
            raise ValueError("IntegrationDomain.description must be nonempty when set.")


@dataclass(frozen=True)
class KernelRepresentation:
    """An explicitly supplied kernel representation of one regularizer.

    ``kernel_expression`` is written in terms of ``output_variable`` and
    ``input_variable``.  Its presence records a caller-supplied representation;
    it does not establish convergence or existence of a bounded operator.
    """

    operator: RegularizerOperator
    kernel_expression: sp.Expr
    output_variable: sp.Symbol
    input_variable: sp.Symbol
    integration_domain: IntegrationDomain
    semantic_status: KernelRepresentationStatus
    hypotheses: tuple[object, ...] = ()
    evidence: object | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.operator, RegularizerOperator):
            raise TypeError("operator must be a RegularizerOperator.")
        if not isinstance(self.kernel_expression, sp.Expr):
            raise TypeError("kernel_expression must be a SymPy expression.")
        if not isinstance(self.output_variable, sp.Symbol):
            raise TypeError("output_variable must be a SymPy Symbol.")
        if not isinstance(self.input_variable, sp.Symbol):
            raise TypeError("input_variable must be a SymPy Symbol.")
        if self.output_variable == self.input_variable:
            raise ValueError("kernel output and input variables must be distinct.")
        if not isinstance(self.integration_domain, IntegrationDomain):
            raise TypeError("integration_domain must be an IntegrationDomain.")
        if not isinstance(self.semantic_status, KernelRepresentationStatus):
            raise TypeError("semantic_status must be a KernelRepresentationStatus.")
        if isinstance(self.hypotheses, (str, bytes)):
            raise TypeError("hypotheses must be an iterable of hypothesis objects.")
        try:
            hypotheses = tuple(self.hypotheses)
        except TypeError as exc:
            raise TypeError("hypotheses must be an iterable.") from exc
        object.__setattr__(self, "hypotheses", hypotheses)
        if (
            self.semantic_status
            is KernelRepresentationStatus.EXTERNALLY_CERTIFIED
            and self.evidence is None
        ):
            raise ValueError(
                "An externally certified KernelRepresentation requires caller-supplied "
                "evidence; the program does not verify that evidence."
            )

    def instantiate(
        self,
        output_variable: sp.Symbol,
        input_variable: sp.Symbol,
    ) -> sp.Expr:
        """Instantiate the supplied kernel in two explicit variables."""

        if not isinstance(output_variable, sp.Symbol):
            raise TypeError("output_variable must be a SymPy Symbol.")
        if not isinstance(input_variable, sp.Symbol):
            raise TypeError("input_variable must be a SymPy Symbol.")
        if output_variable == input_variable:
            raise ValueError("kernel output and input variables must be distinct.")
        return self.kernel_expression.xreplace(
            {
                self.output_variable: output_variable,
                self.input_variable: input_variable,
            }
        )


@dataclass(frozen=True)
class KernelAnnotatedExpression:
    """A SymPy expression retaining every explicit kernel representation used."""

    expression: sp.Expr
    kernel_representations: tuple[KernelRepresentation, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.expression, sp.Expr):
            raise TypeError("expression must be a SymPy expression.")
        try:
            representations = tuple(self.kernel_representations)
        except TypeError as exc:
            raise TypeError("kernel_representations must be an iterable.") from exc
        if not representations or not all(
            isinstance(item, KernelRepresentation) for item in representations
        ):
            raise TypeError(
                "kernel_representations must contain KernelRepresentation objects."
            )
        object.__setattr__(self, "kernel_representations", representations)

    @property
    def hypotheses(self) -> tuple[object, ...]:
        """Return all hypotheses retained from the kernel representations."""

        return tuple(
            hypothesis
            for representation in self.kernel_representations
            for hypothesis in representation.hypotheses
        )

    @property
    def semantic_statuses(self) -> tuple[KernelRepresentationStatus, ...]:
        """Return the semantic status of every kernel representation used."""

        return tuple(
            representation.semantic_status
            for representation in self.kernel_representations
        )

    @property
    def evidences(self) -> tuple[object | None, ...]:
        """Return evidence objects without interpreting them as verified proofs."""

        return tuple(
            representation.evidence
            for representation in self.kernel_representations
        )

    @property
    def represented_operators(self) -> tuple[RegularizerOperator, ...]:
        """Return the regularizer associated with each retained representation."""

        return tuple(
            representation.operator
            for representation in self.kernel_representations
        )

    def as_expr(self) -> sp.Expr:
        """Explicitly project to the underlying SymPy expression."""

        return self.expression

    @property
    def free_symbols(self) -> set[sp.Symbol]:
        """Expose free symbols without discarding the semantic wrapper."""

        return self.expression.free_symbols

    @property
    def variables(self) -> list[sp.Symbol]:
        """Expose integral variables when the wrapped expression is an Integral."""

        if not isinstance(self.expression, sp.Integral):
            raise AttributeError("variables is available only for an Integral expression.")
        return self.expression.variables

    def atoms(self, *types: type[sp.Basic]) -> set[sp.Basic]:
        """Delegate SymPy atom inspection while retaining semantic metadata."""

        return self.expression.atoms(*types)

    def has(self, *patterns: object) -> bool:
        """Delegate structural containment checks to the wrapped expression."""

        return self.expression.has(*patterns)
