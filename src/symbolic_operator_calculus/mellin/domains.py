"""Typed scalar Mellin domains composed from the P0-B safety objects."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import final

import sympy as sp

from ..domains import (
    AssumptionContext,
    ComplexDomain,
    ConsistencyStatus,
    IncompatibleDomainError,
    MembershipStatus,
)
from ..semantics import CertificationStatus
from ..singularities import SingularSet
from .variables import MellinVariable, MellinVariableRole, MellinVariableRoleError


class MellinDomainError(ValueError):
    """Raised when typed Mellin-domain metadata is inconsistent."""


class MellinFrequencyMismatchError(MellinDomainError):
    """Raised when semantic components refer to different frequencies."""


class MellinDomainRoleError(MellinDomainError, MellinVariableRoleError):
    """Raised when a domain slot receives the wrong variable role."""


def _as_tuple(value: object, name: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be an iterable of MellinVariable objects.")
    try:
        return tuple(value)  # type: ignore[arg-type]
    except TypeError as exc:
        raise TypeError(
            f"{name} must be an iterable of MellinVariable objects."
        ) from exc


def _deduplicate_variables(values: tuple[MellinVariable, ...]) -> tuple[MellinVariable, ...]:
    return tuple(dict.fromkeys(values))


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
            raise TypeError("domain evidence objects must be hashable.") from exc
        if item not in result:
            result.append(item)
    return tuple(result)


def _merge_variable_contexts(
    context: AssumptionContext,
    variables: tuple[MellinVariable, ...],
) -> AssumptionContext:
    result = context
    for variable in variables:
        result = result.combine(variable.assumption_context)
    return result


@final
@dataclass(frozen=True)
class MellinSymbolDomain:
    """A scalar-symbol domain that conservatively composes all P0-B metadata."""

    frequency: MellinVariable
    complex_domain: ComplexDomain
    spatial_variables: tuple[MellinVariable, ...] = ()
    relative_variable: MellinVariable | None = None
    parameters: tuple[MellinVariable, ...] = ()
    assumption_context: AssumptionContext = AssumptionContext()
    singular_set: SingularSet = SingularSet()
    description: str | None = None
    evidence: object | tuple[object, ...] | None = ()
    certification_status: CertificationStatus = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.frequency, MellinVariable):
            raise TypeError("frequency must be a MellinVariable.")
        if self.frequency.role is not MellinVariableRole.FREQUENCY:
            raise MellinDomainRoleError("frequency must have role FREQUENCY.")
        if not isinstance(self.complex_domain, ComplexDomain):
            raise TypeError("complex_domain must be a ComplexDomain.")
        if self.complex_domain.variable != self.frequency.symbol:
            raise MellinFrequencyMismatchError(
                "the ComplexDomain principal variable must be the Mellin frequency."
            )
        if not isinstance(self.assumption_context, AssumptionContext):
            raise TypeError("assumption_context must be an AssumptionContext.")
        if not isinstance(self.singular_set, SingularSet):
            raise TypeError("singular_set must be a SingularSet.")

        spatial_raw = _as_tuple(self.spatial_variables, "spatial_variables")
        parameter_raw = _as_tuple(self.parameters, "parameters")
        if not all(isinstance(item, MellinVariable) for item in spatial_raw):
            raise TypeError("spatial_variables must contain MellinVariable objects.")
        if not all(isinstance(item, MellinVariable) for item in parameter_raw):
            raise TypeError("parameters must contain MellinVariable objects.")
        spatial = _deduplicate_variables(spatial_raw)  # type: ignore[arg-type]
        parameters = _deduplicate_variables(parameter_raw)  # type: ignore[arg-type]

        spatial_roles = {
            MellinVariableRole.OUTPUT_SPATIAL,
            MellinVariableRole.INPUT_SPATIAL,
        }
        if any(item.role not in spatial_roles for item in spatial):
            raise MellinDomainRoleError(
                "spatial_variables must have OUTPUT_SPATIAL or INPUT_SPATIAL roles."
            )
        if any(item.role is not MellinVariableRole.PARAMETER for item in parameters):
            raise MellinDomainRoleError("parameters must have role PARAMETER.")
        if self.relative_variable is not None:
            if not isinstance(self.relative_variable, MellinVariable):
                raise TypeError("relative_variable must be a MellinVariable or None.")
            if (
                self.relative_variable.role
                is not MellinVariableRole.RELATIVE_MULTIPLICATIVE
            ):
                raise MellinDomainRoleError(
                    "relative_variable must have role RELATIVE_MULTIPLICATIVE."
                )

        declared = (
            self.frequency,
            *spatial,
            *((self.relative_variable,) if self.relative_variable is not None else ()),
            *parameters,
        )
        roles_by_symbol: dict[object, MellinVariableRole] = {}
        for variable in declared:
            previous = roles_by_symbol.get(variable.symbol)
            if previous is not None and previous is not variable.role:
                raise MellinDomainRoleError(
                    "one SymPy symbol cannot have multiple roles inside one domain."
                )
            roles_by_symbol[variable.symbol] = variable.role
        if len(roles_by_symbol) != len(declared):
            raise MellinDomainError(
                "each symbolic variable may occupy only one domain slot."
            )
        declared_symbols = {item.symbol for item in declared}
        if any(
            item.variable not in declared_symbols
            for item in self.singular_set.singularities
        ) or any(
            item.variable not in declared_symbols
            for item in self.singular_set.avoidances
        ):
            raise MellinFrequencyMismatchError(
                "Mellin-domain singular metadata must use a declared typed variable."
            )

        context = self.assumption_context.combine(
            self.complex_domain.assumption_context
        ).combine(self.singular_set.assumption_context)
        context = _merge_variable_contexts(context, declared)
        if context.consistency_status is ConsistencyStatus.INCONSISTENT:
            raise IncompatibleDomainError(
                "a Mellin domain cannot contain contradictory assumptions."
            )
        complex_domain = self.complex_domain.with_assumptions(context)
        singular_set = SingularSet(
            singularities=self.singular_set.singularities,
            assumption_context=context,
            avoidances=self.singular_set.avoidances,
            description=self.singular_set.description,
        )
        if self.description is not None and not self.description.strip():
            raise MellinDomainError("description must be nonempty when supplied.")
        evidence = _normalize_evidence(
            (*_normalize_evidence(self.evidence), *complex_domain.evidence)
        )

        object.__setattr__(self, "spatial_variables", spatial)
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "assumption_context", context)
        object.__setattr__(self, "complex_domain", complex_domain)
        object.__setattr__(self, "singular_set", singular_set)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(
            self,
            "certification_status",
            (
                CertificationStatus.EVIDENCE_SUPPLIED
                if evidence
                else CertificationStatus.UNCERTIFIED
            ),
        )

    @property
    def declared_variables(self) -> tuple[MellinVariable, ...]:
        """Return every declaration in a deterministic semantic order."""

        return (
            self.frequency,
            *self.spatial_variables,
            *((self.relative_variable,) if self.relative_variable is not None else ()),
            *self.parameters,
        )

    def require_role(
        self,
        variable: MellinVariable,
        role: MellinVariableRole,
    ) -> MellinVariable:
        """Return a stored variable only when its explicit role matches."""

        if not isinstance(variable, MellinVariable):
            raise TypeError("variable must be a MellinVariable.")
        if not isinstance(role, MellinVariableRole):
            raise TypeError("role must be a MellinVariableRole.")
        if variable not in self.declared_variables:
            raise MellinDomainRoleError("variable is not declared in this domain.")
        if variable.role is not role:
            raise MellinDomainRoleError(
                f"expected role {role.value}, got {variable.role.value}."
            )
        return variable

    def with_assumptions(self, *assumptions: object) -> MellinSymbolDomain:
        """Narrow the domain while retaining all prior semantic metadata."""

        additional = (
            assumptions[0]
            if len(assumptions) == 1
            and isinstance(assumptions[0], AssumptionContext)
            else AssumptionContext(tuple(assumptions))
        )
        return MellinSymbolDomain(
            frequency=self.frequency,
            complex_domain=self.complex_domain,
            spatial_variables=self.spatial_variables,
            relative_variable=self.relative_variable,
            parameters=self.parameters,
            assumption_context=self.assumption_context.combine(additional),
            singular_set=self.singular_set,
            description=self.description,
            evidence=self.evidence,
        )

    def with_singularities(self, singular_set: SingularSet) -> MellinSymbolDomain:
        """Conservatively union singular metadata; records are never removed."""

        if not isinstance(singular_set, SingularSet):
            raise TypeError("singular_set must be a SingularSet.")
        return MellinSymbolDomain(
            frequency=self.frequency,
            complex_domain=self.complex_domain,
            spatial_variables=self.spatial_variables,
            relative_variable=self.relative_variable,
            parameters=self.parameters,
            assumption_context=self.assumption_context,
            singular_set=self.singular_set.union(singular_set),
            description=self.description,
            evidence=self.evidence,
        )

    def intersect(self, other: MellinSymbolDomain) -> MellinSymbolDomain:
        """Return a conservative typed intersection without widening either side."""

        if not isinstance(other, MellinSymbolDomain):
            raise TypeError("other must be a MellinSymbolDomain.")
        if self.frequency != other.frequency:
            raise MellinFrequencyMismatchError(
                "Mellin domains with different frequencies cannot be intersected."
            )
        if (
            self.relative_variable is not None
            and other.relative_variable is not None
            and self.relative_variable != other.relative_variable
        ):
            raise MellinDomainRoleError(
                "Mellin domains declare incompatible relative variables."
            )
        spatial = _deduplicate_variables(
            (*self.spatial_variables, *other.spatial_variables)
        )
        parameters = _deduplicate_variables((*self.parameters, *other.parameters))
        relative = self.relative_variable or other.relative_variable
        descriptions = tuple(
            dict.fromkeys(
                item
                for item in (self.description, other.description)
                if item is not None
            )
        )
        return MellinSymbolDomain(
            frequency=self.frequency,
            complex_domain=self.complex_domain.intersect(other.complex_domain),
            spatial_variables=spatial,
            relative_variable=relative,
            parameters=parameters,
            assumption_context=self.assumption_context.combine(
                other.assumption_context
            ),
            singular_set=self.singular_set.union(other.singular_set),
            description=" ∩ ".join(descriptions) if descriptions else None,
            evidence=(*self.evidence, *other.evidence),
        )

    def is_compatible_with(self, other: MellinSymbolDomain) -> MembershipStatus:
        """Report only compatibility established by the limited stored rules."""

        if not isinstance(other, MellinSymbolDomain):
            raise TypeError("other must be a MellinSymbolDomain.")
        if self.frequency != other.frequency:
            return MembershipStatus.NO
        try:
            intersection = self.intersect(other)
        except (IncompatibleDomainError, MellinDomainError):
            return MembershipStatus.NO
        if (
            intersection.assumption_context.consistency_status
            is ConsistencyStatus.UNDETERMINED
        ):
            return MembershipStatus.UNDETERMINED
        return MembershipStatus.YES

    def substitute(
        self,
        substitutions: dict[sp.Basic, sp.Basic],
    ) -> MellinSymbolDomain:
        """Apply role-preserving scalar substitutions to every metadata layer."""

        if not isinstance(substitutions, dict):
            raise TypeError("substitutions must be a dict.")
        normalized = {
            sp.sympify(source): sp.sympify(target)
            for source, target in substitutions.items()
        }
        declared_by_symbol = {
            variable.symbol: variable for variable in self.declared_variables
        }
        for source, target in normalized.items():
            if source not in declared_by_symbol:
                continue
            source_role = declared_by_symbol[source].role
            for target_symbol in target.free_symbols:
                target_declaration = declared_by_symbol.get(target_symbol)
                if (
                    target_declaration is not None
                    and target_declaration.role is not source_role
                ):
                    raise MellinDomainRoleError(
                        "substitution cannot change a declared variable role."
                    )
                if target_declaration is None and target_symbol != source:
                    if isinstance(target, sp.Symbol):
                        continue
                    raise MellinDomainRoleError(
                        "substitution cannot introduce an undeclared free symbol."
                    )

        substituted_complex_domain = self.complex_domain.substitute(normalized)
        substituted_singular_set = self.singular_set.substitute(normalized)

        def substitute_declaration(
            declaration: MellinVariable,
        ) -> MellinVariable | None:
            target = normalized.get(declaration.symbol, declaration.symbol)
            if not isinstance(target, sp.Symbol):
                if target.free_symbols:
                    existing = tuple(
                        variable
                        for variable in self.declared_variables
                        if variable.symbol in target.free_symbols
                        and variable.role is declaration.role
                    )
                    if len(existing) == 1 and len(target.free_symbols) == 1:
                        return existing[0]
                    raise MellinDomainRoleError(
                        "a non-symbol replacement cannot define a new variable role."
                    )
                return None
            return MellinVariable(
                symbol=target,
                role=declaration.role,
                assumption_context=declaration.assumption_context.substitute(
                    normalized
                ),
                description=declaration.description,
            )

        frequency = substitute_declaration(self.frequency)
        if frequency is None:
            raise MellinFrequencyMismatchError(
                "the Mellin frequency can only be renamed to another Symbol."
            )
        spatial = tuple(
            replacement
            for item in self.spatial_variables
            if (replacement := substitute_declaration(item)) is not None
        )
        relative = (
            None
            if self.relative_variable is None
            else substitute_declaration(self.relative_variable)
        )
        parameters = tuple(
            replacement
            for item in self.parameters
            if (replacement := substitute_declaration(item)) is not None
        )
        return MellinSymbolDomain(
            frequency=frequency,
            complex_domain=substituted_complex_domain,
            spatial_variables=spatial,
            relative_variable=relative,
            parameters=parameters,
            assumption_context=self.assumption_context.substitute(normalized),
            singular_set=substituted_singular_set,
            description=self.description,
            evidence=self.evidence,
        )

    def __str__(self) -> str:
        return (
            "MellinSymbolDomain("
            f"frequency={self.frequency.symbol}; "
            f"spatial={len(self.spatial_variables)}; "
            f"relative={self.relative_variable is not None}; "
            f"parameters={len(self.parameters)}; "
            f"assumptions={len(self.assumption_context.assumptions)}; "
            f"singularities={len(self.singular_set.singularities)}; "
            f"status={self.certification_status.value})"
        )
