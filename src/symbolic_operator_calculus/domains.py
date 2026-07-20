"""Conservative assumptions and declared complex domains.

This module records semantic conditions.  It deliberately implements only
small, auditable consistency and membership rules; ``UNDETERMINED`` never
means that a condition has been proved consistent or a point has been proved
to belong to a domain.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import final

import sympy as sp
from sympy.core.relational import (
    Equality,
    GreaterThan,
    LessThan,
    Relational,
    StrictGreaterThan,
    StrictLessThan,
    Unequality,
)
from sympy.logic.boolalg import Boolean, BooleanAtom, BooleanFunction

from .semantics import CertificationStatus


class ConsistencyStatus(str, Enum):
    """Result of the deliberately limited assumption-consistency checks."""

    CONSISTENT = "consistent"
    INCONSISTENT = "inconsistent"
    UNDETERMINED = "undetermined"


class InvalidAssumptionError(TypeError):
    """Raised when a value is not an explicit SymPy proposition."""


class MembershipStatus(str, Enum):
    """Trivalent result for declared-domain membership and containment."""

    YES = "yes"
    NO = "no"
    UNDETERMINED = "undetermined"


class DomainRegionKind(str, Enum):
    """Supported elementary regions; arbitrary geometry is out of scope."""

    COMPLEX_PLANE = "complex_plane"
    REAL_LINE = "real_line"
    UPPER_HALF_PLANE = "upper_half_plane"
    LOWER_HALF_PLANE = "lower_half_plane"
    VERTICAL_STRIP = "vertical_strip"
    HORIZONTAL_STRIP = "horizontal_strip"
    ASSUMPTION_DEFINED = "assumption_defined"


class InvalidDomainError(ValueError):
    """Raised when a declared domain is malformed."""


class IncompatibleDomainError(InvalidDomainError):
    """Raised for an elementarily certified empty domain intersection."""


class DomainVariableMismatchError(InvalidDomainError):
    """Raised when domains use different principal variables."""


class DomainSubstitutionError(InvalidDomainError):
    """Raised when substitution would capture or invalidate a domain variable."""


def _structural_key(value: sp.Basic) -> str:
    """Return a deterministic key without asking for symbolic truth."""

    return sp.srepr(value)


def _as_proposition(value: object) -> sp.Basic:
    try:
        proposition = sp.sympify(value)
    except (sp.SympifyError, TypeError) as exc:
        raise InvalidAssumptionError(
            "assumptions must be SymPy Boolean propositions."
        ) from exc
    valid_boolean = isinstance(
        proposition, (Relational, BooleanFunction, BooleanAtom)
    ) or (
        isinstance(proposition, Boolean)
        and not isinstance(proposition, sp.Expr)
    )
    if not valid_boolean:
        raise InvalidAssumptionError(
            "assumptions must be explicit relations or Boolean expressions; "
            "plain symbolic expressions are not propositions."
        )
    return proposition


def _has_literal_contradiction(assumptions: tuple[sp.Basic, ...]) -> bool:
    keys = {_structural_key(item) for item in assumptions}
    if _structural_key(sp.false) in keys:
        return True
    for assumption in assumptions:
        negated = sp.Not(assumption)
        if _structural_key(negated) in keys:
            return True
    return False


@dataclass
class _NumericBounds:
    lower: sp.Expr | None = None
    lower_strict: bool = False
    upper: sp.Expr | None = None
    upper_strict: bool = False
    equality: sp.Expr | None = None
    excluded: set[sp.Expr] = field(default_factory=set)


def _symbol_and_numeric_value(
    relation: Relational,
) -> tuple[sp.Symbol, sp.Expr, type[Relational]] | None:
    lhs, rhs = relation.lhs, relation.rhs
    if isinstance(lhs, sp.Symbol) and rhs.is_number and rhs.is_real is True:
        return lhs, rhs, type(relation)
    if isinstance(rhs, sp.Symbol) and lhs.is_number and lhs.is_real is True:
        reversed_type: dict[type[Relational], type[Relational]] = {
            StrictGreaterThan: StrictLessThan,
            GreaterThan: LessThan,
            StrictLessThan: StrictGreaterThan,
            LessThan: GreaterThan,
            Equality: Equality,
            Unequality: Unequality,
        }
        return rhs, lhs, reversed_type[type(relation)]
    return None


def _number_order(left: sp.Expr, right: sp.Expr) -> int | None:
    """Compare certified real numbers, returning -1, 0, 1, or unknown."""

    difference = sp.simplify(left - right)
    if difference.is_zero is True:
        return 0
    if difference.is_negative is True:
        return -1
    if difference.is_positive is True:
        return 1
    return None


def _stronger_lower(
    current: sp.Expr | None,
    current_strict: bool,
    candidate: sp.Expr,
    candidate_strict: bool,
) -> tuple[sp.Expr, bool] | None:
    if current is None:
        return candidate, candidate_strict
    comparison = _number_order(candidate, current)
    if comparison is None:
        return None
    if comparison > 0:
        return candidate, candidate_strict
    if comparison < 0:
        return current, current_strict
    return current, current_strict or candidate_strict


def _stronger_upper(
    current: sp.Expr | None,
    current_strict: bool,
    candidate: sp.Expr,
    candidate_strict: bool,
) -> tuple[sp.Expr, bool] | None:
    if current is None:
        return candidate, candidate_strict
    comparison = _number_order(candidate, current)
    if comparison is None:
        return None
    if comparison < 0:
        return candidate, candidate_strict
    if comparison > 0:
        return current, current_strict
    return current, current_strict or candidate_strict


def _certify_numeric_bound_system(
    assumptions: tuple[sp.Basic, ...],
) -> ConsistencyStatus:
    bounds: dict[sp.Symbol, _NumericBounds] = {}
    for assumption in assumptions:
        if assumption is sp.true:
            continue
        if not isinstance(assumption, Relational):
            return ConsistencyStatus.UNDETERMINED
        parsed = _symbol_and_numeric_value(assumption)
        if parsed is None:
            return ConsistencyStatus.UNDETERMINED
        symbol, value, relation_type = parsed
        state = bounds.setdefault(symbol, _NumericBounds())
        if relation_type is StrictGreaterThan:
            result = _stronger_lower(state.lower, state.lower_strict, value, True)
            if result is None:
                return ConsistencyStatus.UNDETERMINED
            state.lower, state.lower_strict = result
        elif relation_type is GreaterThan:
            result = _stronger_lower(state.lower, state.lower_strict, value, False)
            if result is None:
                return ConsistencyStatus.UNDETERMINED
            state.lower, state.lower_strict = result
        elif relation_type is StrictLessThan:
            result = _stronger_upper(state.upper, state.upper_strict, value, True)
            if result is None:
                return ConsistencyStatus.UNDETERMINED
            state.upper, state.upper_strict = result
        elif relation_type is LessThan:
            result = _stronger_upper(state.upper, state.upper_strict, value, False)
            if result is None:
                return ConsistencyStatus.UNDETERMINED
            state.upper, state.upper_strict = result
        elif relation_type is Equality:
            if state.equality is not None:
                comparison = _number_order(state.equality, value)
                if comparison is None:
                    return ConsistencyStatus.UNDETERMINED
                if comparison != 0:
                    return ConsistencyStatus.INCONSISTENT
            state.equality = value
        elif relation_type is Unequality:
            state.excluded.add(value)

    for state in bounds.values():
        if state.lower is not None and state.upper is not None:
            comparison = _number_order(state.lower, state.upper)
            if comparison is None:
                return ConsistencyStatus.UNDETERMINED
            if comparison > 0 or (
                comparison == 0 and (state.lower_strict or state.upper_strict)
            ):
                return ConsistencyStatus.INCONSISTENT
        if state.equality is None:
            continue
        if state.equality in state.excluded:
            return ConsistencyStatus.INCONSISTENT
        if state.lower is not None:
            comparison = _number_order(state.equality, state.lower)
            if comparison is None:
                return ConsistencyStatus.UNDETERMINED
            if comparison < 0 or (comparison == 0 and state.lower_strict):
                return ConsistencyStatus.INCONSISTENT
        if state.upper is not None:
            comparison = _number_order(state.equality, state.upper)
            if comparison is None:
                return ConsistencyStatus.UNDETERMINED
            if comparison > 0 or (comparison == 0 and state.upper_strict):
                return ConsistencyStatus.INCONSISTENT
    return ConsistencyStatus.CONSISTENT


def _consistency_status(
    assumptions: tuple[sp.Basic, ...],
) -> ConsistencyStatus:
    if not assumptions:
        return ConsistencyStatus.CONSISTENT
    if _has_literal_contradiction(assumptions):
        return ConsistencyStatus.INCONSISTENT
    return _certify_numeric_bound_system(assumptions)


@final
@dataclass(frozen=True)
class AssumptionContext:
    """An immutable, presentation-ordered collection of explicit hypotheses.

    Exact duplicates are removed by structural SymPy representation while the
    first presentation order is retained.  Consistency checking is restricted
    to literal contradictions and elementary real numeric bounds.
    """

    assumptions: tuple[sp.Basic, ...] = ()
    consistency_status: ConsistencyStatus = field(init=False)

    def __post_init__(self) -> None:
        if isinstance(self.assumptions, (str, bytes)):
            raise InvalidAssumptionError(
                "assumptions must be an iterable of SymPy propositions."
            )
        try:
            supplied = tuple(self.assumptions)
        except TypeError as exc:
            raise InvalidAssumptionError(
                "assumptions must be an iterable of SymPy propositions."
            ) from exc
        normalized: list[sp.Basic] = []
        seen: set[str] = set()
        for value in supplied:
            proposition = _as_proposition(value)
            if proposition is sp.true:
                continue
            key = _structural_key(proposition)
            if key not in seen:
                seen.add(key)
                normalized.append(proposition)
        canonical = tuple(normalized)
        object.__setattr__(self, "assumptions", canonical)
        object.__setattr__(
            self,
            "consistency_status",
            _consistency_status(canonical),
        )

    @property
    def is_empty(self) -> bool:
        """Whether no explicit hypotheses were supplied."""

        return not self.assumptions

    def with_assumption(self, assumption: object) -> AssumptionContext:
        """Return a new context retaining every existing hypothesis."""

        return AssumptionContext((*self.assumptions, _as_proposition(assumption)))

    def combine(self, other: AssumptionContext) -> AssumptionContext:
        """Return the conservative ordered union of two contexts."""

        if not isinstance(other, AssumptionContext):
            raise TypeError("other must be an AssumptionContext.")
        return AssumptionContext((*self.assumptions, *other.assumptions))

    def contains_exact(self, assumption: object) -> bool:
        """Check structural presence without evaluating symbolic truth."""

        key = _structural_key(_as_proposition(assumption))
        return any(_structural_key(item) == key for item in self.assumptions)

    def contains_all(self, other: AssumptionContext) -> bool:
        """Check that all required assumptions occur structurally."""

        if not isinstance(other, AssumptionContext):
            raise TypeError("other must be an AssumptionContext.")
        return all(self.contains_exact(item) for item in other.assumptions)

    def as_boolean_expression(self) -> sp.Basic:
        """Explicitly project the context to a SymPy Boolean expression."""

        if not self.assumptions:
            return sp.true
        return sp.And(*self.assumptions)

    def substitute(self, substitutions: dict[sp.Basic, sp.Basic]) -> AssumptionContext:
        """Apply explicit structural substitutions and recheck consistency."""

        if not isinstance(substitutions, dict):
            raise TypeError("substitutions must be a dict.")
        normalized = {
            sp.sympify(source): sp.sympify(target)
            for source, target in substitutions.items()
        }
        return AssumptionContext(
            tuple(assumption.xreplace(normalized) for assumption in self.assumptions)
        )

    def __str__(self) -> str:
        if not self.assumptions:
            return "AssumptionContext(empty; consistent)"
        rendered = "; ".join(sp.sstr(item) for item in self.assumptions)
        return f"AssumptionContext({rendered}; {self.consistency_status.value})"


def _deduplicate_basic(values: tuple[sp.Basic, ...]) -> tuple[sp.Basic, ...]:
    result: list[sp.Basic] = []
    seen: set[str] = set()
    for value in values:
        key = _structural_key(value)
        if key not in seen:
            seen.add(key)
            result.append(value)
    return tuple(result)


def _relational_membership(proposition: sp.Basic) -> MembershipStatus:
    if proposition is sp.true:
        return MembershipStatus.YES
    if proposition is sp.false:
        return MembershipStatus.NO
    simplified = sp.simplify(proposition)
    if simplified is sp.true:
        return MembershipStatus.YES
    if simplified is sp.false:
        return MembershipStatus.NO
    return MembershipStatus.UNDETERMINED


@final
@dataclass(frozen=True)
class DomainRegion:
    """One elementary geometric or explicitly conditioned region."""

    kind: DomainRegionKind
    lower: sp.Expr | None = None
    upper: sp.Expr | None = None
    condition: sp.Basic | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.kind, DomainRegionKind):
            raise TypeError("kind must be a DomainRegionKind.")
        lower = None if self.lower is None else sp.sympify(self.lower)
        upper = None if self.upper is None else sp.sympify(self.upper)
        condition = (
            None if self.condition is None else _as_proposition(self.condition)
        )
        strip_kinds = {
            DomainRegionKind.VERTICAL_STRIP,
            DomainRegionKind.HORIZONTAL_STRIP,
        }
        if self.kind in strip_kinds:
            if lower is None or upper is None:
                raise InvalidDomainError("a strip requires lower and upper bounds.")
            if (
                lower.is_extended_real is not True
                or upper.is_extended_real is not True
            ):
                raise InvalidDomainError(
                    "strip bounds must be provably extended-real values."
                )
            comparison = _number_order(lower, upper)
            if comparison is not None and comparison >= 0:
                raise IncompatibleDomainError(
                    "a strict strip requires a certified lower bound below its upper bound."
                )
        elif lower is not None or upper is not None:
            raise InvalidDomainError("bounds are supported only for strip regions.")
        if self.kind is DomainRegionKind.ASSUMPTION_DEFINED:
            if condition is None:
                raise InvalidDomainError(
                    "an assumption-defined region requires a Boolean condition."
                )
        elif condition is not None:
            raise InvalidDomainError(
                "a condition is supported only for an assumption-defined region."
            )
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)
        object.__setattr__(self, "condition", condition)

    def substitute(
        self,
        substitutions: dict[sp.Basic, sp.Basic],
    ) -> DomainRegion:
        """Return this region after an explicit structural substitution."""

        return DomainRegion(
            self.kind,
            None if self.lower is None else self.lower.xreplace(substitutions),
            None if self.upper is None else self.upper.xreplace(substitutions),
            (
                None
                if self.condition is None
                else self.condition.xreplace(substitutions)
            ),
        )


def _strip_intervals_are_incompatible(regions: tuple[DomainRegion, ...]) -> bool:
    for strip_kind in (
        DomainRegionKind.VERTICAL_STRIP,
        DomainRegionKind.HORIZONTAL_STRIP,
    ):
        relevant = [region for region in regions if region.kind is strip_kind]
        for index, left in enumerate(relevant):
            for right in relevant[index + 1 :]:
                left_before_right = _number_order(left.upper, right.lower)
                right_before_left = _number_order(right.upper, left.lower)
                if left_before_right is not None and left_before_right <= 0:
                    return True
                if right_before_left is not None and right_before_left <= 0:
                    return True
    return False


def _zero_belongs_to_strict_interval(region: DomainRegion) -> MembershipStatus:
    lower = _relational_membership(sp.StrictLessThan(region.lower, 0))
    upper = _relational_membership(sp.StrictLessThan(0, region.upper))
    if MembershipStatus.NO in (lower, upper):
        return MembershipStatus.NO
    if MembershipStatus.UNDETERMINED in (lower, upper):
        return MembershipStatus.UNDETERMINED
    return MembershipStatus.YES


def _validate_region_intersection(regions: tuple[DomainRegion, ...]) -> None:
    kinds = {region.kind for region in regions}
    if {
        DomainRegionKind.UPPER_HALF_PLANE,
        DomainRegionKind.LOWER_HALF_PLANE,
    }.issubset(kinds):
        raise IncompatibleDomainError(
            "upper and lower open half-planes have empty intersection."
        )
    if DomainRegionKind.REAL_LINE in kinds and (
        DomainRegionKind.UPPER_HALF_PLANE in kinds
        or DomainRegionKind.LOWER_HALF_PLANE in kinds
    ):
        raise IncompatibleDomainError(
            "the real line does not intersect either open imaginary half-plane."
        )
    if _strip_intervals_are_incompatible(regions):
        raise IncompatibleDomainError("the declared strict strips are disjoint.")
    horizontal = [
        region
        for region in regions
        if region.kind is DomainRegionKind.HORIZONTAL_STRIP
    ]
    if DomainRegionKind.REAL_LINE in kinds and any(
        _zero_belongs_to_strict_interval(region) is MembershipStatus.NO
        for region in horizontal
    ):
        raise IncompatibleDomainError(
            "the horizontal strip excludes the real line."
        )


def _region_membership(
    region: DomainRegion,
    variable: sp.Symbol,
    point: sp.Expr,
) -> MembershipStatus:
    if region.kind is DomainRegionKind.COMPLEX_PLANE:
        return MembershipStatus.YES
    if region.kind is DomainRegionKind.REAL_LINE:
        return _relational_membership(sp.Eq(sp.im(point), 0))
    if region.kind is DomainRegionKind.UPPER_HALF_PLANE:
        return _relational_membership(sp.StrictGreaterThan(sp.im(point), 0))
    if region.kind is DomainRegionKind.LOWER_HALF_PLANE:
        return _relational_membership(sp.StrictLessThan(sp.im(point), 0))
    coordinate = (
        sp.re(point)
        if region.kind is DomainRegionKind.VERTICAL_STRIP
        else sp.im(point)
    )
    if region.kind in {
        DomainRegionKind.VERTICAL_STRIP,
        DomainRegionKind.HORIZONTAL_STRIP,
    }:
        lower = _relational_membership(
            sp.StrictLessThan(region.lower, coordinate)
        )
        upper = _relational_membership(
            sp.StrictLessThan(coordinate, region.upper)
        )
        if MembershipStatus.NO in (lower, upper):
            return MembershipStatus.NO
        if MembershipStatus.UNDETERMINED in (lower, upper):
            return MembershipStatus.UNDETERMINED
        return MembershipStatus.YES
    if region.condition is None:
        raise InvalidDomainError("assumption-defined region lost its condition.")
    return _relational_membership(region.condition.xreplace({variable: point}))


def _exclusion_membership(exclusion: sp.Basic, point: sp.Expr) -> MembershipStatus:
    if isinstance(exclusion, sp.Set):
        return _relational_membership(exclusion.contains(point))
    return _relational_membership(sp.Eq(point, exclusion))


@final
@dataclass(frozen=True)
class ComplexDomain:
    """A declared complex domain with hypotheses and explicit exclusions.

    Construction records a region of validity.  It does not establish
    convergence, boundedness, or existence of an associated operator.
    """

    variable: sp.Symbol
    regions: tuple[DomainRegion, ...] = (
        DomainRegion(DomainRegionKind.COMPLEX_PLANE),
    )
    assumption_context: AssumptionContext = AssumptionContext()
    exclusions: tuple[sp.Basic, ...] = ()
    description: str | None = None
    evidence: object | tuple[object, ...] | None = ()
    certification_status: CertificationStatus = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.variable, sp.Symbol):
            raise TypeError("variable must be a SymPy Symbol.")
        if not isinstance(self.assumption_context, AssumptionContext):
            raise TypeError("assumption_context must be an AssumptionContext.")
        if (
            self.assumption_context.consistency_status
            is ConsistencyStatus.INCONSISTENT
        ):
            raise IncompatibleDomainError(
                "a domain cannot be constructed with inconsistent assumptions."
            )
        try:
            supplied_regions = tuple(self.regions)
        except TypeError as exc:
            raise TypeError("regions must be an iterable of DomainRegion objects.") from exc
        if not supplied_regions or not all(
            isinstance(region, DomainRegion) for region in supplied_regions
        ):
            raise InvalidDomainError(
                "regions must contain at least one DomainRegion."
            )
        regions = tuple(
            region
            for region in supplied_regions
            if region.kind is not DomainRegionKind.COMPLEX_PLANE
        )
        if not regions:
            regions = (DomainRegion(DomainRegionKind.COMPLEX_PLANE),)
        regions = tuple(
            dict.fromkeys(regions)
        )
        _validate_region_intersection(regions)
        try:
            raw_exclusions = tuple(self.exclusions)
        except TypeError as exc:
            raise TypeError("exclusions must be an iterable.") from exc
        normalized_exclusions: list[sp.Basic] = []
        for exclusion in raw_exclusions:
            try:
                normalized = sp.sympify(exclusion)
            except (sp.SympifyError, TypeError) as exc:
                raise InvalidDomainError(
                    "each exclusion must be a SymPy point or set."
                ) from exc
            if isinstance(normalized, (Relational, BooleanFunction, BooleanAtom)):
                raise InvalidDomainError(
                    "Boolean exclusions must be represented as SymPy sets."
                )
            normalized_exclusions.append(normalized)
        exclusions = _deduplicate_basic(tuple(normalized_exclusions))
        if self.description is not None and not self.description.strip():
            raise InvalidDomainError("description must be nonempty when supplied.")
        if self.evidence is None:
            evidence: tuple[object, ...] = ()
        elif isinstance(self.evidence, tuple):
            evidence = self.evidence
        else:
            evidence = (self.evidence,)
        for item in evidence:
            try:
                hash(item)
            except TypeError as exc:
                raise TypeError("domain evidence objects must be hashable.") from exc
        object.__setattr__(self, "regions", regions)
        object.__setattr__(self, "exclusions", exclusions)
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

    @classmethod
    def complex_plane(
        cls,
        variable: sp.Symbol,
        **metadata: object,
    ) -> ComplexDomain:
        return cls(variable=variable, **metadata)

    @classmethod
    def real_line(
        cls,
        variable: sp.Symbol,
        **metadata: object,
    ) -> ComplexDomain:
        return cls(
            variable=variable,
            regions=(DomainRegion(DomainRegionKind.REAL_LINE),),
            **metadata,
        )

    @classmethod
    def upper_half_plane(
        cls,
        variable: sp.Symbol,
        **metadata: object,
    ) -> ComplexDomain:
        return cls(
            variable=variable,
            regions=(DomainRegion(DomainRegionKind.UPPER_HALF_PLANE),),
            **metadata,
        )

    @classmethod
    def lower_half_plane(
        cls,
        variable: sp.Symbol,
        **metadata: object,
    ) -> ComplexDomain:
        return cls(
            variable=variable,
            regions=(DomainRegion(DomainRegionKind.LOWER_HALF_PLANE),),
            **metadata,
        )

    @classmethod
    def vertical_strip(
        cls,
        variable: sp.Symbol,
        lower: object,
        upper: object,
        **metadata: object,
    ) -> ComplexDomain:
        return cls(
            variable=variable,
            regions=(
                DomainRegion(DomainRegionKind.VERTICAL_STRIP, lower, upper),
            ),
            **metadata,
        )

    @classmethod
    def horizontal_strip(
        cls,
        variable: sp.Symbol,
        lower: object,
        upper: object,
        **metadata: object,
    ) -> ComplexDomain:
        return cls(
            variable=variable,
            regions=(
                DomainRegion(DomainRegionKind.HORIZONTAL_STRIP, lower, upper),
            ),
            **metadata,
        )

    @classmethod
    def defined_by(
        cls,
        variable: sp.Symbol,
        condition: object,
        **metadata: object,
    ) -> ComplexDomain:
        return cls(
            variable=variable,
            regions=(
                DomainRegion(
                    DomainRegionKind.ASSUMPTION_DEFINED,
                    condition=condition,
                ),
            ),
            **metadata,
        )

    def with_exclusions(self, *exclusions: object) -> ComplexDomain:
        """Return a narrowed domain retaining every prior exclusion."""

        return ComplexDomain(
            self.variable,
            self.regions,
            self.assumption_context,
            (*self.exclusions, *exclusions),
            self.description,
            self.evidence,
        )

    def with_assumptions(
        self,
        *assumptions: object,
    ) -> ComplexDomain:
        """Return a domain narrowed by an additional explicit context."""

        if len(assumptions) == 1 and isinstance(
            assumptions[0], AssumptionContext
        ):
            additional = assumptions[0]
        else:
            additional = AssumptionContext(tuple(assumptions))
        return ComplexDomain(
            self.variable,
            self.regions,
            self.assumption_context.combine(additional),
            self.exclusions,
            self.description,
            self.evidence,
        )

    def intersect(self, other: ComplexDomain) -> ComplexDomain:
        """Return the conservative intersection of two declared domains."""

        if not isinstance(other, ComplexDomain):
            raise TypeError("other must be a ComplexDomain.")
        if self.variable != other.variable:
            raise DomainVariableMismatchError(
                "domains with different principal variables cannot be intersected."
            )
        description_parts = list(
            dict.fromkeys(
                part
                for part in (self.description, other.description)
                if part is not None
            )
        )
        return ComplexDomain(
            variable=self.variable,
            regions=(*self.regions, *other.regions),
            assumption_context=self.assumption_context.combine(
                other.assumption_context
            ),
            exclusions=(*self.exclusions, *other.exclusions),
            description=(
                " ∩ ".join(description_parts) if description_parts else None
            ),
            evidence=(*self.evidence, *other.evidence),
        )

    def contains_point(self, point: object) -> MembershipStatus:
        """Return only membership that the limited internal rules establish."""

        point_expression = sp.sympify(point)
        complex_membership = _relational_membership(
            sp.S.Complexes.contains(point_expression)
        )
        if complex_membership is MembershipStatus.NO:
            return MembershipStatus.NO
        if (
            self.assumption_context.consistency_status
            is ConsistencyStatus.UNDETERMINED
        ):
            return MembershipStatus.UNDETERMINED
        result = complex_membership
        for region in self.regions:
            membership = _region_membership(region, self.variable, point_expression)
            if membership is MembershipStatus.NO:
                return MembershipStatus.NO
            if membership is MembershipStatus.UNDETERMINED:
                result = MembershipStatus.UNDETERMINED
        for exclusion in self.exclusions:
            excluded = _exclusion_membership(exclusion, point_expression)
            if excluded is MembershipStatus.YES:
                return MembershipStatus.NO
            if excluded is MembershipStatus.UNDETERMINED:
                result = MembershipStatus.UNDETERMINED
        return result

    def is_subset_of(self, other: ComplexDomain) -> MembershipStatus:
        """Conservatively compare domains by explicit stored restrictions."""

        if not isinstance(other, ComplexDomain):
            raise TypeError("other must be a ComplexDomain.")
        if self.variable != other.variable:
            return MembershipStatus.NO
        if not self.assumption_context.contains_all(other.assumption_context):
            return MembershipStatus.NO
        self_exclusions = {_structural_key(item) for item in self.exclusions}
        if any(
            _structural_key(item) not in self_exclusions
            for item in other.exclusions
        ):
            return MembershipStatus.NO
        required_regions = tuple(
            region
            for region in other.regions
            if region.kind is not DomainRegionKind.COMPLEX_PLANE
        )
        if not required_regions:
            return MembershipStatus.YES
        available_regions = set(self.regions)
        if all(region in available_regions for region in required_regions):
            return MembershipStatus.YES
        if all(
            region.kind is DomainRegionKind.COMPLEX_PLANE for region in self.regions
        ):
            return MembershipStatus.NO
        return MembershipStatus.UNDETERMINED

    def substitute(self, substitutions: dict[sp.Basic, sp.Basic]) -> ComplexDomain:
        """Substitute parameters without capturing the principal variable."""

        if not isinstance(substitutions, dict):
            raise TypeError("substitutions must be a dict.")
        normalized = {
            sp.sympify(source): sp.sympify(target)
            for source, target in substitutions.items()
        }
        for source, target in normalized.items():
            if source != self.variable and self.variable in target.free_symbols:
                raise DomainSubstitutionError(
                    "a parameter substitution cannot capture the domain variable."
                )
        new_variable = normalized.get(self.variable, self.variable)
        if not isinstance(new_variable, sp.Symbol):
            raise DomainSubstitutionError(
                "the principal domain variable can only be renamed to a Symbol."
            )
        return ComplexDomain(
            variable=new_variable,
            regions=tuple(
                region.substitute(normalized) for region in self.regions
            ),
            assumption_context=self.assumption_context.substitute(normalized),
            exclusions=tuple(
                exclusion.xreplace(normalized) for exclusion in self.exclusions
            ),
            description=self.description,
            evidence=self.evidence,
        )

    def __str__(self) -> str:
        region_text = " ∩ ".join(region.kind.value for region in self.regions)
        exclusions = ", ".join(sp.sstr(item) for item in self.exclusions) or "none"
        return (
            f"ComplexDomain(variable={self.variable}; region={region_text}; "
            f"assumptions={len(self.assumption_context.assumptions)}; "
            f"exclusions={exclusions})"
        )
