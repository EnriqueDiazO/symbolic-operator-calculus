"""Metadata-preserving formal scalar Mellin expressions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import final

import sympy as sp

from ..conditional import ConditionalVerificationStatus
from ..domains import ConsistencyStatus
from ..singularities import detect_singularities
from .domains import MellinDomainError, MellinSymbolDomain
from .variables import MellinVariable, MellinVariableRole


class MellinExpressionError(ValueError):
    """Base error for malformed or unsupported scalar Mellin expressions."""


class UndeclaredMellinVariableError(MellinExpressionError):
    """Raised when a free or bound symbol lacks an explicit declaration."""


class MellinBranchConditionRequiredError(MellinExpressionError):
    """Raised when a branch-sensitive power lacks a positive-base condition."""


class MellinSingularMetadataRequiredError(MellinExpressionError):
    """Raised when limited P0-B detection finds an undeclared singularity."""


def _normalize_evidence(
    evidence: object | tuple[object, ...] | None,
) -> tuple[object, ...]:
    if evidence is None:
        values: tuple[object, ...] = ()
    elif isinstance(evidence, tuple):
        values = evidence
    else:
        values = (evidence,)
    result: list[object] = []
    for item in values:
        try:
            hash(item)
        except TypeError as exc:
            raise TypeError("expression evidence objects must be hashable.") from exc
        if item not in result:
            result.append(item)
    return tuple(result)


def _as_variable_tuple(value: object, name: str) -> tuple[MellinVariable, ...]:
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be an iterable of MellinVariable objects.")
    try:
        supplied = tuple(value)  # type: ignore[arg-type]
    except TypeError as exc:
        raise TypeError(
            f"{name} must be an iterable of MellinVariable objects."
        ) from exc
    if not all(isinstance(item, MellinVariable) for item in supplied):
        raise TypeError(f"{name} must contain MellinVariable objects.")
    return tuple(dict.fromkeys(supplied))


def _bound_symbols(expression: sp.Expr) -> set[sp.Symbol]:
    result: set[sp.Symbol] = set()
    for binder in expression.atoms(sp.Integral, sp.Sum, sp.Product):
        result.update(binder.variables)
    for function in expression.atoms(sp.Lambda):
        variables = function.variables
        if isinstance(variables, tuple):
            result.update(variables)
        else:
            result.add(variables)
    return result


def _singularity_key(item: object) -> tuple[object, ...]:
    return (
        getattr(item, "location", None),
        getattr(item, "kind", None),
        getattr(item, "variable", None),
        getattr(item, "order", None),
        getattr(item, "conditions", None),
    )


def _missing_singularities(
    expression: sp.Expr,
    domain: MellinSymbolDomain,
) -> tuple[object, ...]:
    detected = detect_singularities(
        expression,
        domain.frequency.symbol,
        assumptions=domain.assumption_context,
        domain=domain.complex_domain,
    )
    supplied = {
        _singularity_key(item) for item in domain.singular_set.singularities
    }
    return tuple(
        item
        for item in detected.singularities
        if _singularity_key(item) not in supplied
    )


def _branch_status_is_safe(
    expression: sp.Expr,
    domain: MellinSymbolDomain,
) -> None:
    parameter_symbols = {item.symbol for item in domain.parameters}
    frequency = domain.frequency.symbol
    for power in expression.atoms(sp.Pow):
        base, exponent = power.as_base_exp()
        if (
            not isinstance(base, sp.Symbol)
            or base not in parameter_symbols
            or frequency not in exponent.free_symbols
            or exponent.is_integer is True
        ):
            continue
        status = domain.singular_set.positive_base_treatment_status(
            base,
            domain.assumption_context,
        )
        if status.value == "yes":
            continue
        if status.value == "no":
            raise MellinBranchConditionRequiredError(
                f"branch-sensitive base {base} is declared nonpositive."
            )
        raise MellinBranchConditionRequiredError(
            f"branch-sensitive base {base} requires the explicit assumption {base} > 0."
        )


def _weaker_status(
    left: ConditionalVerificationStatus,
    right: ConditionalVerificationStatus,
) -> ConditionalVerificationStatus:
    if left is right:
        return left
    return ConditionalVerificationStatus.DECLARED


def _detected_union(
    expression: sp.Expr,
    domain: MellinSymbolDomain,
) -> MellinSymbolDomain:
    detected = detect_singularities(
        expression,
        domain.frequency.symbol,
        assumptions=domain.assumption_context,
        domain=domain.complex_domain,
    )
    return domain.with_singularities(detected)


@final
@dataclass(frozen=True)
class MellinExpression:
    """A formal scalar expression whose semantic metadata is never implicit."""

    expression: sp.Expr
    domain: MellinSymbolDomain
    variables: tuple[MellinVariable, ...] = ()
    parameters: tuple[MellinVariable, ...] = ()
    dummy_variables: tuple[sp.Symbol, ...] = ()
    verification_status: ConditionalVerificationStatus = (
        ConditionalVerificationStatus.DECLARED
    )
    evidence: object | tuple[object, ...] | None = ()
    description: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.expression, sp.Expr):
            raise TypeError("expression must be an existing SymPy scalar expression.")
        if not isinstance(self.domain, MellinSymbolDomain):
            raise TypeError("domain must be a MellinSymbolDomain.")
        variables = _as_variable_tuple(self.variables, "variables")
        parameters = _as_variable_tuple(self.parameters, "parameters")
        if any(item.role is MellinVariableRole.PARAMETER for item in variables):
            raise MellinExpressionError(
                "parameter declarations belong in parameters, not variables."
            )
        if any(item.role is not MellinVariableRole.PARAMETER for item in parameters):
            raise MellinExpressionError(
                "parameters must contain only PARAMETER declarations."
            )
        if any(item not in self.domain.declared_variables for item in variables):
            raise MellinExpressionError(
                "each expression variable must be declared by its Mellin domain."
            )
        if any(item not in self.domain.parameters for item in parameters):
            raise MellinExpressionError(
                "each expression parameter must be declared by its Mellin domain."
            )
        declared_symbols = [item.symbol for item in (*variables, *parameters)]
        if len(set(declared_symbols)) != len(declared_symbols):
            raise MellinExpressionError(
                "one expression symbol cannot be declared with multiple roles."
            )
        try:
            dummy_variables = tuple(self.dummy_variables)
        except TypeError as exc:
            raise TypeError("dummy_variables must be an iterable of Symbols.") from exc
        if not all(isinstance(item, sp.Symbol) for item in dummy_variables):
            raise TypeError("dummy_variables must contain only SymPy Symbols.")
        dummy_variables = tuple(dict.fromkeys(dummy_variables))
        required_dummies = _bound_symbols(self.expression) | set(
            self.expression.atoms(sp.Dummy)
        )
        if not required_dummies.issubset(set(dummy_variables)):
            missing = required_dummies.difference(dummy_variables)
            raise UndeclaredMellinVariableError(
                "bound or Dummy symbols require explicit dummy_variables metadata: "
                + ", ".join(sorted(sp.sstr(item) for item in missing))
            )
        allowed = set(declared_symbols) | set(dummy_variables)
        missing_free = self.expression.free_symbols.difference(allowed)
        if missing_free:
            raise UndeclaredMellinVariableError(
                "free symbols lack Mellin role declarations: "
                + ", ".join(sorted(sp.sstr(item) for item in missing_free))
            )
        if not isinstance(
            self.verification_status, ConditionalVerificationStatus
        ):
            raise TypeError(
                "verification_status must be a ConditionalVerificationStatus."
            )
        evidence = _normalize_evidence(self.evidence)
        if (
            self.verification_status
            is ConditionalVerificationStatus.EVIDENCE_SUPPLIED
            and not evidence
        ):
            raise MellinExpressionError(
                "EVIDENCE_SUPPLIED requires an external evidence object."
            )
        if (
            self.verification_status
            is ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
        ):
            raise MellinExpressionError(
                "a standalone expression cannot be promoted to an internally checked identity."
            )
        if self.description is not None and not self.description.strip():
            raise MellinExpressionError(
                "description must be nonempty when supplied."
            )
        if (
            self.domain.assumption_context.consistency_status
            is ConsistencyStatus.INCONSISTENT
        ):
            raise MellinDomainError(
                "an expression cannot use an inconsistent Mellin domain."
            )
        _branch_status_is_safe(self.expression, self.domain)
        missing_singularities = _missing_singularities(
            self.expression,
            self.domain,
        )
        if missing_singularities:
            kinds = ", ".join(
                sorted({item.kind.value for item in missing_singularities})
            )
            raise MellinSingularMetadataRequiredError(
                "limited P0-B rules detected undeclared singular metadata: "
                f"{kinds}."
            )
        object.__setattr__(self, "variables", variables)
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "dummy_variables", dummy_variables)
        object.__setattr__(self, "evidence", evidence)

    @property
    def status(self) -> ConditionalVerificationStatus:
        return self.verification_status

    @property
    def free_symbols(self) -> set[sp.Symbol]:
        return set(self.expression.free_symbols)

    def _combine(
        self,
        other: MellinExpression,
        operation: object,
        description: str,
    ) -> MellinExpression:
        if not isinstance(other, MellinExpression):
            raise TypeError("scalar Mellin operations require another MellinExpression.")
        domain = self.domain.intersect(other.domain)
        expression = operation(self.expression, other.expression)  # type: ignore[operator]
        domain = _detected_union(expression, domain)
        return MellinExpression(
            expression=expression,
            domain=domain,
            variables=tuple(dict.fromkeys((*self.variables, *other.variables))),
            parameters=tuple(dict.fromkeys((*self.parameters, *other.parameters))),
            dummy_variables=tuple(
                dict.fromkeys((*self.dummy_variables, *other.dummy_variables))
            ),
            verification_status=_weaker_status(
                self.verification_status,
                other.verification_status,
            ),
            evidence=(*self.evidence, *other.evidence),
            description=description,
        )

    def __add__(self, other: MellinExpression) -> MellinExpression:
        return self._combine(other, lambda left, right: left + right, "Scalar sum")

    def __mul__(self, other: MellinExpression) -> MellinExpression:
        return self._combine(
            other,
            lambda left, right: left * right,
            "Scalar product",
        )

    def __neg__(self) -> MellinExpression:
        return MellinExpression(
            expression=-self.expression,
            domain=self.domain,
            variables=self.variables,
            parameters=self.parameters,
            dummy_variables=self.dummy_variables,
            verification_status=self.verification_status,
            evidence=self.evidence,
            description=(
                f"Negative of {self.description}"
                if self.description is not None
                else "Scalar negation"
            ),
        )

    def substitute(
        self,
        substitutions: dict[sp.Basic, sp.Basic],
    ) -> MellinExpression:
        """Substitute all scalar and domain metadata without inventing roles."""

        if not isinstance(substitutions, dict):
            raise TypeError("substitutions must be a dict.")
        normalized = {
            sp.sympify(source): sp.sympify(target)
            for source, target in substitutions.items()
        }
        if any(source in self.dummy_variables for source in normalized):
            raise MellinExpressionError(
                "dummy-variable substitutions are unsupported; rebuild explicitly."
            )
        domain = self.domain.substitute(normalized)
        expression = self.expression.xreplace(normalized)
        domain = _detected_union(expression, domain)
        declarations = {
            item.symbol: item for item in domain.declared_variables
        }
        used = expression.free_symbols.difference(self.dummy_variables)
        variables = tuple(
            item
            for item in domain.declared_variables
            if item.symbol in used and item.role is not MellinVariableRole.PARAMETER
        )
        parameters = tuple(
            item for item in domain.parameters if item.symbol in used
        )
        if used.difference(declarations):
            raise UndeclaredMellinVariableError(
                "substitution introduced free symbols without declared roles."
            )
        return MellinExpression(
            expression=expression,
            domain=domain,
            variables=variables,
            parameters=parameters,
            dummy_variables=self.dummy_variables,
            verification_status=self.verification_status,
            evidence=self.evidence,
            description=self.description,
        )

    def differentiate(self, variable: sp.Symbol) -> MellinExpression:
        """Differentiate a scalar formula and conservatively union singularities."""

        if not isinstance(variable, sp.Symbol):
            raise TypeError("variable must be a SymPy Symbol.")
        declared = {item.symbol for item in (*self.variables, *self.parameters)}
        if variable not in declared:
            raise UndeclaredMellinVariableError(
                "differentiation variable must have a declared Mellin role."
            )
        expression = sp.diff(self.expression, variable)
        domain = _detected_union(expression, self.domain)
        return MellinExpression(
            expression=expression,
            domain=domain,
            variables=self.variables,
            parameters=self.parameters,
            dummy_variables=self.dummy_variables,
            verification_status=self.verification_status,
            evidence=self.evidence,
            description=(
                f"Scalar derivative of {self.description}"
                if self.description is not None
                else "Scalar derivative"
            ),
        )

    def conjugate(self) -> MellinExpression:
        """Form an unevaluated scalar conjugate without simplifying branches."""

        expression = sp.conjugate(self.expression, evaluate=False)
        domain = _detected_union(expression, self.domain)
        return MellinExpression(
            expression=expression,
            domain=domain,
            variables=self.variables,
            parameters=self.parameters,
            dummy_variables=self.dummy_variables,
            verification_status=self.verification_status,
            evidence=self.evidence,
            description=(
                f"Formal conjugate of {self.description}"
                if self.description is not None
                else "Formal scalar conjugate"
            ),
        )

    def as_expr(self) -> sp.Expr:
        """Explicitly project to SymPy, deliberately losing all metadata."""

        return self.expression

    def __str__(self) -> str:
        return (
            f"Formal Mellin expression({sp.sstr(self.expression)}; "
            f"variables={len(self.variables)}; parameters={len(self.parameters)}; "
            f"status={self.verification_status.value})"
        )
