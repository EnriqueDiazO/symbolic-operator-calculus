"""Declared singular sets and deliberately limited internal detection rules."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import final

import sympy as sp

from .domains import (
    AssumptionContext,
    ComplexDomain,
    ConsistencyStatus,
    MembershipStatus,
)
from .semantics import CertificationStatus


class SingularityKind(str, Enum):
    """Kinds recorded by P0-B without claiming complete complex analysis."""

    POLE = "pole"
    BRANCH_POINT = "branch_point"
    BRANCH_CUT = "branch_cut"
    REMOVABLE = "removable"
    UNKNOWN = "unknown"
    GENERAL_EXCLUSION = "general_exclusion"


class SingularityOrigin(str, Enum):
    """Where a singularity record came from."""

    USER_DECLARED = "user_declared"
    INTERNAL_RULE = "internal_rule"
    EXTERNAL_EVIDENCE = "external_evidence"


class SingularityDetectionError(ValueError):
    """Raised when a requested limited singularity analysis is malformed."""


def _deduplicate_records(records: tuple[object, ...]) -> tuple[object, ...]:
    result: list[object] = []
    seen: set[str] = set()
    for record in records:
        key = repr(record)
        if key not in seen:
            seen.add(key)
            result.append(record)
    return tuple(result)


@final
@dataclass(frozen=True)
class Singularity:
    """One declared or specifically detected singularity."""

    location: sp.Basic
    kind: SingularityKind
    variable: sp.Symbol
    order: int | None = None
    conditions: AssumptionContext = AssumptionContext()
    origin: SingularityOrigin = SingularityOrigin.USER_DECLARED
    evidence: object | None = None
    description: str | None = None
    certification_status: CertificationStatus = field(init=False)

    def __post_init__(self) -> None:
        try:
            location = sp.sympify(self.location)
        except (sp.SympifyError, TypeError) as exc:
            raise SingularityDetectionError(
                "singularity location must be a SymPy expression or set."
            ) from exc
        if not isinstance(self.kind, SingularityKind):
            raise TypeError("kind must be a SingularityKind.")
        if not isinstance(self.variable, sp.Symbol):
            raise TypeError("variable must be a SymPy Symbol.")
        if self.order is not None and (
            not isinstance(self.order, int)
            or isinstance(self.order, bool)
            or self.order <= 0
        ):
            raise ValueError("singularity order must be a positive integer.")
        if self.kind is not SingularityKind.POLE and self.order is not None:
            raise ValueError("only a pole may carry an order in P0-B.")
        if not isinstance(self.conditions, AssumptionContext):
            raise TypeError("conditions must be an AssumptionContext.")
        if not isinstance(self.origin, SingularityOrigin):
            raise TypeError("origin must be a SingularityOrigin.")
        if (
            self.origin is SingularityOrigin.EXTERNAL_EVIDENCE
            and self.evidence is None
        ):
            raise ValueError(
                "an externally supplied singularity requires an evidence object."
            )
        if self.evidence is not None:
            try:
                hash(self.evidence)
            except TypeError as exc:
                raise TypeError("singularity evidence must be hashable.") from exc
        if self.description is not None and not self.description.strip():
            raise ValueError("description must be nonempty when supplied.")
        object.__setattr__(self, "location", location)
        object.__setattr__(
            self,
            "certification_status",
            (
                CertificationStatus.EVIDENCE_SUPPLIED
                if self.evidence is not None
                else CertificationStatus.UNCERTIFIED
            ),
        )

    def substitute(
        self,
        substitutions: dict[sp.Basic, sp.Basic],
    ) -> Singularity:
        """Substitute without changing kind, origin, or evidence."""

        new_variable = substitutions.get(self.variable, self.variable)
        if not isinstance(new_variable, sp.Symbol):
            raise SingularityDetectionError(
                "a singularity variable can only be renamed to a Symbol."
            )
        return Singularity(
            location=self.location.xreplace(substitutions),
            kind=self.kind,
            variable=new_variable,
            order=self.order,
            conditions=self.conditions.substitute(substitutions),
            origin=self.origin,
            evidence=self.evidence,
            description=self.description,
        )


@final
@dataclass(frozen=True)
class SingularityAvoidance:
    """A limited rule recording when one stored singularity is avoided."""

    location: sp.Basic
    kind: SingularityKind
    variable: sp.Symbol
    conditions: AssumptionContext
    origin: SingularityOrigin = SingularityOrigin.INTERNAL_RULE
    evidence: object | None = None
    description: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "location", sp.sympify(self.location))
        if not isinstance(self.kind, SingularityKind):
            raise TypeError("kind must be a SingularityKind.")
        if not isinstance(self.variable, sp.Symbol):
            raise TypeError("variable must be a SymPy Symbol.")
        if not isinstance(self.conditions, AssumptionContext):
            raise TypeError("conditions must be an AssumptionContext.")
        if not isinstance(self.origin, SingularityOrigin):
            raise TypeError("origin must be a SingularityOrigin.")
        if (
            self.origin is SingularityOrigin.EXTERNAL_EVIDENCE
            and self.evidence is None
        ):
            raise ValueError(
                "an externally supplied avoidance requires an evidence object."
            )
        if self.evidence is not None:
            try:
                hash(self.evidence)
            except TypeError as exc:
                raise TypeError("avoidance evidence must be hashable.") from exc
        if self.description is not None and not self.description.strip():
            raise ValueError("description must be nonempty when supplied.")

    def substitute(
        self,
        substitutions: dict[sp.Basic, sp.Basic],
    ) -> SingularityAvoidance:
        new_variable = substitutions.get(self.variable, self.variable)
        if not isinstance(new_variable, sp.Symbol):
            raise SingularityDetectionError(
                "an avoidance variable can only be renamed to a Symbol."
            )
        return SingularityAvoidance(
            location=self.location.xreplace(substitutions),
            kind=self.kind,
            variable=new_variable,
            conditions=self.conditions.substitute(substitutions),
            origin=self.origin,
            evidence=self.evidence,
            description=self.description,
        )


def _contains_integer_status(
    parameter: sp.Symbol,
    context: AssumptionContext,
) -> MembershipStatus:
    integer_condition = sp.Contains(parameter, sp.S.Integers)
    if context.contains_exact(integer_condition):
        return MembershipStatus.YES
    if context.contains_exact(sp.Not(integer_condition)):
        return MembershipStatus.NO
    for assumption in context.assumptions:
        if not isinstance(assumption, sp.Equality):
            continue
        if assumption.lhs == parameter and assumption.rhs.is_number:
            value = assumption.rhs
        elif assumption.rhs == parameter and assumption.lhs.is_number:
            value = assumption.lhs
        else:
            continue
        if value.is_integer is True:
            return MembershipStatus.YES
        if value.is_integer is False:
            return MembershipStatus.NO
    if context.contains_exact(parameter > 0) and context.contains_exact(parameter < 1):
        return MembershipStatus.NO
    return MembershipStatus.UNDETERMINED


def _condition_status(
    conditions: AssumptionContext,
    available_context: AssumptionContext,
) -> MembershipStatus:
    if not conditions.assumptions:
        return MembershipStatus.YES
    statuses: list[MembershipStatus] = []
    for condition in conditions.assumptions:
        if available_context.contains_exact(condition):
            statuses.append(MembershipStatus.YES)
            continue
        if available_context.contains_exact(sp.Not(condition)):
            return MembershipStatus.NO
        if isinstance(condition, sp.Contains) and condition.args[1] == sp.S.Integers:
            element = condition.args[0]
            if isinstance(element, sp.Symbol):
                statuses.append(
                    _contains_integer_status(element, available_context)
                )
                continue
        statuses.append(MembershipStatus.UNDETERMINED)
    if MembershipStatus.NO in statuses:
        return MembershipStatus.NO
    if MembershipStatus.UNDETERMINED in statuses:
        return MembershipStatus.UNDETERMINED
    return MembershipStatus.YES


def _location_contains(location: sp.Basic, point: sp.Expr) -> MembershipStatus:
    proposition = (
        location.contains(point)
        if isinstance(location, sp.Set)
        else sp.Eq(point, location)
    )
    simplified = sp.simplify(proposition)
    if simplified is sp.true:
        return MembershipStatus.YES
    if simplified is sp.false:
        return MembershipStatus.NO
    return MembershipStatus.UNDETERMINED


def _avoidance_matches(
    singularity: Singularity,
    avoidance: SingularityAvoidance,
) -> bool:
    return (
        singularity.location == avoidance.location
        and singularity.kind is avoidance.kind
        and singularity.variable == avoidance.variable
    )


@final
@dataclass(frozen=True)
class SingularSet:
    """Immutable singularities plus explicit conditioned avoidance records."""

    singularities: tuple[Singularity, ...] = ()
    assumption_context: AssumptionContext = AssumptionContext()
    avoidances: tuple[SingularityAvoidance, ...] = ()
    description: str | None = None

    def __post_init__(self) -> None:
        try:
            singularities = tuple(self.singularities)
            avoidances = tuple(self.avoidances)
        except TypeError as exc:
            raise TypeError("singularities and avoidances must be iterables.") from exc
        if not all(isinstance(item, Singularity) for item in singularities):
            raise TypeError("singularities must contain Singularity objects.")
        if not all(isinstance(item, SingularityAvoidance) for item in avoidances):
            raise TypeError("avoidances must contain SingularityAvoidance objects.")
        if not isinstance(self.assumption_context, AssumptionContext):
            raise TypeError("assumption_context must be an AssumptionContext.")
        if self.description is not None and not self.description.strip():
            raise ValueError("description must be nonempty when supplied.")
        object.__setattr__(
            self,
            "singularities",
            _deduplicate_records(singularities),
        )
        object.__setattr__(self, "avoidances", _deduplicate_records(avoidances))

    @property
    def is_empty(self) -> bool:
        return not self.singularities

    def with_singularity(self, singularity: Singularity) -> SingularSet:
        """Return a union that cannot discard existing singularities."""

        if not isinstance(singularity, Singularity):
            raise TypeError("singularity must be a Singularity.")
        return SingularSet(
            (*self.singularities, singularity),
            self.assumption_context,
            self.avoidances,
            self.description,
        )

    def union(self, other: SingularSet) -> SingularSet:
        """Conservatively union records, contexts, and separate avoidances."""

        if not isinstance(other, SingularSet):
            raise TypeError("other must be a SingularSet.")
        descriptions = [
            item for item in (self.description, other.description) if item is not None
        ]
        return SingularSet(
            singularities=(*self.singularities, *other.singularities),
            assumption_context=self.assumption_context.combine(
                other.assumption_context
            ),
            avoidances=(*self.avoidances, *other.avoidances),
            description=" ∪ ".join(descriptions) if descriptions else None,
        )

    def contains_point(
        self,
        point: object,
        *,
        variable: sp.Symbol,
        available_context: AssumptionContext | None = None,
    ) -> MembershipStatus:
        """Check whether a stored singularity is active at one point."""

        if not isinstance(variable, sp.Symbol):
            raise TypeError("variable must be a SymPy Symbol.")
        context = (
            self.assumption_context
            if available_context is None
            else available_context
        )
        if not isinstance(context, AssumptionContext):
            raise TypeError("available_context must be an AssumptionContext.")
        if context.consistency_status is ConsistencyStatus.INCONSISTENT:
            return MembershipStatus.UNDETERMINED
        point_expression = sp.sympify(point)
        result = MembershipStatus.NO
        for singularity in self.singularities:
            if singularity.variable != variable:
                continue
            location = _location_contains(singularity.location, point_expression)
            condition = _condition_status(singularity.conditions, context)
            if location is MembershipStatus.YES and condition is MembershipStatus.YES:
                return MembershipStatus.YES
            if MembershipStatus.UNDETERMINED in (location, condition):
                result = MembershipStatus.UNDETERMINED
        return result

    def avoidance_status(
        self,
        domain: ComplexDomain,
        *,
        available_context: AssumptionContext | None = None,
    ) -> MembershipStatus:
        """Check whether all relevant stored singularities avoid a domain."""

        if not isinstance(domain, ComplexDomain):
            raise TypeError("domain must be a ComplexDomain.")
        context = domain.assumption_context
        if available_context is not None:
            if not isinstance(available_context, AssumptionContext):
                raise TypeError("available_context must be an AssumptionContext.")
            context = context.combine(available_context)
        if context.consistency_status is ConsistencyStatus.INCONSISTENT:
            return MembershipStatus.UNDETERMINED
        result = MembershipStatus.YES
        for singularity in self.singularities:
            matching_avoidance = any(
                _avoidance_matches(singularity, avoidance)
                and context.contains_all(avoidance.conditions)
                and avoidance.conditions.consistency_status
                is not ConsistencyStatus.INCONSISTENT
                for avoidance in self.avoidances
            )
            if matching_avoidance:
                continue
            if singularity.variable != domain.variable:
                result = MembershipStatus.UNDETERMINED
                continue
            condition = _condition_status(singularity.conditions, context)
            if condition is MembershipStatus.NO:
                continue
            if isinstance(singularity.location, sp.Set):
                result = MembershipStatus.UNDETERMINED
                continue
            location = domain.contains_point(singularity.location)
            if condition is MembershipStatus.YES and location is MembershipStatus.YES:
                return MembershipStatus.NO
            if MembershipStatus.UNDETERMINED in (condition, location):
                result = MembershipStatus.UNDETERMINED
        return result

    def positive_base_treatment_status(
        self,
        base: sp.Symbol,
        available_context: AssumptionContext,
    ) -> MembershipStatus:
        """Recognize only an explicit ``base > 0`` branch condition."""

        if not isinstance(base, sp.Symbol):
            raise TypeError("base must be a SymPy Symbol.")
        if not isinstance(available_context, AssumptionContext):
            raise TypeError("available_context must be an AssumptionContext.")
        if available_context.contains_exact(base > 0):
            return MembershipStatus.YES
        if available_context.contains_exact(base <= 0):
            return MembershipStatus.NO
        return MembershipStatus.UNDETERMINED

    def substitute(self, substitutions: dict[sp.Basic, sp.Basic]) -> SingularSet:
        """Update every record; no singularity is removed by substitution."""

        if not isinstance(substitutions, dict):
            raise TypeError("substitutions must be a dict.")
        normalized = {
            sp.sympify(source): sp.sympify(target)
            for source, target in substitutions.items()
        }
        return SingularSet(
            singularities=tuple(
                singularity.substitute(normalized)
                for singularity in self.singularities
            ),
            assumption_context=self.assumption_context.substitute(normalized),
            avoidances=tuple(
                avoidance.substitute(normalized) for avoidance in self.avoidances
            ),
            description=self.description,
        )


def _integer_lattice() -> sp.ImageSet:
    integer = sp.Symbol("n", integer=True)
    return sp.ImageSet(sp.Lambda(integer, sp.I * integer), sp.S.Integers)


def _scaled_hyperbolic_argument(
    argument: sp.Expr,
    variable: sp.Symbol,
) -> tuple[str, sp.Symbol | None] | None:
    scaled = sp.expand(argument / sp.pi)
    if scaled == variable:
        return "direct", None
    shift = sp.expand(scaled - variable)
    parameter = sp.simplify(shift / sp.I)
    if isinstance(parameter, sp.Symbol) and sp.expand(
        scaled - (variable + sp.I * parameter)
    ) == 0:
        return "imaginary_shift", parameter
    return None


def _hyperbolic_pole_records(
    argument: sp.Expr,
    variable: sp.Symbol,
    order: int,
    assumptions: AssumptionContext,
    function_name: str,
) -> tuple[tuple[Singularity, ...], tuple[SingularityAvoidance, ...]]:
    match = _scaled_hyperbolic_argument(argument, variable)
    if match is None:
        return (), ()
    match_kind, parameter = match
    if match_kind == "direct":
        return (
            Singularity(
                location=_integer_lattice(),
                kind=SingularityKind.POLE,
                variable=variable,
                order=order,
                origin=SingularityOrigin.INTERNAL_RULE,
                description=(
                    f"limited {function_name}(pi*z) rule: z = I*n, n in Integers"
                ),
            ),
        ), ()
    if parameter is None:
        raise SingularityDetectionError("imaginary-shift rule lost its parameter.")
    singularity = Singularity(
        location=sp.Integer(0),
        kind=SingularityKind.POLE,
        variable=variable,
        order=order,
        conditions=AssumptionContext(
            (sp.Contains(parameter, sp.S.Integers),)
        ),
        origin=SingularityOrigin.INTERNAL_RULE,
        description=(
            f"limited {function_name}(pi*(lambda+I*kappa)) rule: "
            "lambda = 0 when kappa is an integer"
        ),
    )
    avoidances: tuple[SingularityAvoidance, ...] = ()
    if assumptions.contains_exact(parameter > 0) and assumptions.contains_exact(
        parameter < 1
    ):
        avoidances = (
            SingularityAvoidance(
                location=0,
                kind=SingularityKind.POLE,
                variable=variable,
                conditions=AssumptionContext((parameter > 0, parameter < 1)),
                description=(
                    "the limited integer rule certifies no integer parameter "
                    "strictly between 0 and 1"
                ),
            ),
        )
    return (singularity,), avoidances


def _branch_records(
    expression: sp.Expr,
    variable: sp.Symbol,
    assumptions: AssumptionContext,
) -> tuple[tuple[Singularity, ...], tuple[SingularityAvoidance, ...]]:
    singularities: list[Singularity] = []
    avoidances: list[SingularityAvoidance] = []
    seen_bases: set[str] = set()
    def branch_locations(
        base: sp.Expr,
    ) -> tuple[sp.Symbol, sp.Basic, sp.Basic]:
        if isinstance(base, sp.Symbol):
            return base, sp.Integer(0), sp.Interval(-sp.oo, 0)
        branch_point = sp.ConditionSet(
            variable,
            sp.Eq(base, 0),
            sp.S.Complexes,
        )
        branch_cut = sp.ConditionSet(
            variable,
            sp.And(sp.Eq(sp.im(base), 0), sp.re(base) <= 0),
            sp.S.Complexes,
        )
        return variable, branch_point, branch_cut

    for logarithm in expression.atoms(sp.log):
        base = logarithm.args[0]
        if variable not in base.free_symbols:
            continue
        branch_variable, branch_point, branch_cut = branch_locations(base)
        singularities.extend(
            (
                Singularity(
                    location=branch_point,
                    kind=SingularityKind.BRANCH_POINT,
                    variable=branch_variable,
                    origin=SingularityOrigin.INTERNAL_RULE,
                    description="principal SymPy logarithm branch point",
                ),
                Singularity(
                    location=branch_cut,
                    kind=SingularityKind.BRANCH_CUT,
                    variable=branch_variable,
                    origin=SingularityOrigin.INTERNAL_RULE,
                    description="principal SymPy logarithm branch cut",
                ),
            )
        )
        seen_bases.add(sp.srepr(base))
    for power in expression.atoms(sp.Pow):
        base, exponent = power.as_base_exp()
        if exponent.is_integer is True or base.is_number:
            continue
        if variable not in power.free_symbols:
            continue
        base_key = sp.srepr(base)
        if base_key in seen_bases:
            continue
        relevant_variable, branch_point, branch_cut = branch_locations(base)
        singularities.extend(
            (
                Singularity(
                    location=branch_point,
                    kind=SingularityKind.BRANCH_POINT,
                    variable=relevant_variable,
                    origin=SingularityOrigin.INTERNAL_RULE,
                    description="noninteger power is sensitive to the principal branch",
                ),
                Singularity(
                    location=branch_cut,
                    kind=SingularityKind.BRANCH_CUT,
                    variable=relevant_variable,
                    origin=SingularityOrigin.INTERNAL_RULE,
                    description="noninteger power uses SymPy's principal branch",
                ),
            )
        )
        if isinstance(base, sp.Symbol) and assumptions.contains_exact(base > 0):
            avoidances.extend(
                (
                    SingularityAvoidance(
                        location=0,
                        kind=SingularityKind.BRANCH_POINT,
                        variable=base,
                        conditions=AssumptionContext((base > 0,)),
                        description="positive real base avoids the branch point",
                    ),
                    SingularityAvoidance(
                        location=branch_cut,
                        kind=SingularityKind.BRANCH_CUT,
                        variable=base,
                        conditions=AssumptionContext((base > 0,)),
                        description="positive real base uses the real logarithm",
                    ),
                )
            )
    return tuple(singularities), tuple(avoidances)


def detect_singularities(
    expression: sp.Expr,
    variable: sp.Symbol,
    *,
    assumptions: AssumptionContext | None = None,
) -> SingularSet:
    """Apply only the explicit P0-B rules; this is not a universal detector."""

    if not isinstance(expression, sp.Expr):
        raise TypeError("expression must be a SymPy expression.")
    if not isinstance(variable, sp.Symbol):
        raise TypeError("variable must be a SymPy Symbol.")
    context = AssumptionContext() if assumptions is None else assumptions
    if not isinstance(context, AssumptionContext):
        raise TypeError("assumptions must be an AssumptionContext.")
    singularities: list[Singularity] = []
    avoidances: list[SingularityAvoidance] = []
    for coth_term in expression.atoms(sp.coth):
        found, avoided = _hyperbolic_pole_records(
            coth_term.args[0], variable, 1, context, "coth"
        )
        singularities.extend(found)
        avoidances.extend(avoided)
    for power in expression.atoms(sp.Pow):
        base, exponent = power.as_base_exp()
        if base.func is not sp.sinh or not exponent.is_Integer or exponent >= 0:
            continue
        found, avoided = _hyperbolic_pole_records(
            base.args[0], variable, int(-exponent), context, "1/sinh"
        )
        singularities.extend(found)
        avoidances.extend(avoided)
    branch_singularities, branch_avoidances = _branch_records(
        expression, variable, context
    )
    singularities.extend(branch_singularities)
    avoidances.extend(branch_avoidances)
    return SingularSet(
        singularities=tuple(singularities),
        assumption_context=context,
        avoidances=tuple(avoidances),
        description="limited structural singularity analysis",
    )
