"""Conditional scalar identities with explicit applicability barriers."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import final

import sympy as sp

from .domains import (
    AssumptionContext,
    ComplexDomain,
    ConsistencyStatus,
    IncompatibleDomainError,
    MembershipStatus,
)
from .semantics import ExactIdentityScope
from .singularities import SingularSet, detect_singularities


class ConditionalVerificationStatus(str, Enum):
    """Strength of the limited verification recorded for an identity."""

    DECLARED = "declared"
    SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS = (
        "symbolically_checked_under_assumptions"
    )
    EVIDENCE_SUPPLIED = "evidence_supplied"


class ConditionalIdentityError(ValueError):
    """Base error for malformed or unsupported conditional identities."""


class IdentityApplicabilityError(ConditionalIdentityError):
    """Raised when declared conditions do not justify using an identity."""


class IdentityVerificationError(ConditionalIdentityError):
    """Raised when the limited scalar check cannot verify an identity."""


class IdentityVerificationIndeterminateError(IdentityVerificationError):
    """Raised when the limited scalar check cannot decide a zero residual."""


class IdentityScopeError(ConditionalIdentityError):
    """Raised for unsupported or stronger identity scopes."""


class SingularMetadataRequiredError(ConditionalIdentityError):
    """Raised when limited detection finds an undeclared singular feature."""


@final
@dataclass(frozen=True)
class SymbolicCheckRecord:
    """Auditable result of one concrete scalar algebra check."""

    method: str
    residual: sp.Expr
    note: str = (
        "A zero scalar residual is not an operatorial or analytic theorem."
    )

    def __post_init__(self) -> None:
        if self.method != "sympy_simplify_difference":
            raise ValueError("unsupported symbolic-check method.")
        if not isinstance(self.residual, sp.Expr):
            raise TypeError("residual must be a SymPy expression.")
        if self.residual.is_zero is not True:
            raise ValueError("a symbolic check record requires a certified zero residual.")


def _expression_with_both_sides(lhs: sp.Expr, rhs: sp.Expr) -> sp.Expr:
    return sp.Add(lhs, sp.Mul(-1, rhs, evaluate=False), evaluate=False)


def _missing_detected_singularities(
    lhs: sp.Expr,
    rhs: sp.Expr,
    domain: ComplexDomain,
    assumptions: AssumptionContext,
    supplied: SingularSet,
) -> tuple[object, ...]:
    detected = detect_singularities(
        _expression_with_both_sides(lhs, rhs),
        domain.variable,
        assumptions=assumptions,
        domain=domain,
    )
    supplied_keys = {
        (
            item.location,
            item.kind,
            item.variable,
            item.order,
            item.conditions,
        )
        for item in supplied.singularities
    }
    return tuple(
        item
        for item in detected.singularities
        if (
            item.location,
            item.kind,
            item.variable,
            item.order,
            item.conditions,
        )
        not in supplied_keys
    )


def _normalize_evidence(
    evidence: object | tuple[object, ...] | None,
) -> tuple[object, ...]:
    if evidence is None:
        return ()
    if isinstance(evidence, tuple):
        values = evidence
    else:
        values = (evidence,)
    for item in values:
        try:
            hash(item)
        except TypeError as exc:
            raise TypeError("identity evidence objects must be hashable.") from exc
    return values


@final
@dataclass(frozen=True)
class ConditionalIdentity:
    """A scalar identity usable only under its stored semantic conditions."""

    lhs: sp.Expr
    rhs: sp.Expr
    assumption_context: AssumptionContext
    domain: ComplexDomain
    singular_set: SingularSet
    scope: ExactIdentityScope
    verification_status: ConditionalVerificationStatus = (
        ConditionalVerificationStatus.DECLARED
    )
    evidence: object | tuple[object, ...] | None = ()
    description: str | None = None
    verification_method: str | None = None

    def __post_init__(self) -> None:
        try:
            lhs = sp.sympify(self.lhs)
            rhs = sp.sympify(self.rhs)
        except (sp.SympifyError, TypeError) as exc:
            raise TypeError("lhs and rhs must be SymPy scalar expressions.") from exc
        if not isinstance(lhs, sp.Expr) or not isinstance(rhs, sp.Expr):
            raise TypeError("lhs and rhs must be SymPy scalar expressions.")
        if not isinstance(self.assumption_context, AssumptionContext):
            raise TypeError("assumption_context must be an AssumptionContext.")
        if not isinstance(self.domain, ComplexDomain):
            raise TypeError("domain must be a ComplexDomain.")
        if not isinstance(self.singular_set, SingularSet):
            raise TypeError("singular_set must be a SingularSet.")
        if self.scope is not ExactIdentityScope.SCALAR_SYMBOLIC:
            raise IdentityScopeError(
                "P0-B ConditionalIdentity supports only scalar symbolic scope; "
                "operatorial identities are out of scope."
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
            raise ValueError("EVIDENCE_SUPPLIED requires an evidence object.")
        if (
            self.verification_status
            is ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
        ):
            if self.verification_method != "sympy_simplify_difference":
                raise ValueError(
                    "a symbolic check requires its concrete verification method."
                )
            if len(evidence) != 1 or not isinstance(
                evidence[0], SymbolicCheckRecord
            ):
                raise ValueError(
                    "a symbolic check requires exactly one SymbolicCheckRecord."
                )
        if self.description is not None and not self.description.strip():
            raise ValueError("description must be nonempty when supplied.")
        missing = _missing_detected_singularities(
            lhs,
            rhs,
            self.domain,
            self.assumption_context,
            self.singular_set,
        )
        if missing:
            kinds = ", ".join(sorted({item.kind.value for item in missing}))
            raise SingularMetadataRequiredError(
                "limited internal rules detected undeclared singular metadata: "
                f"{kinds}."
            )
        object.__setattr__(self, "lhs", lhs)
        object.__setattr__(self, "rhs", rhs)
        object.__setattr__(self, "evidence", evidence)

    def symbolically_check(self) -> ConditionalIdentity:
        """Run the single supported scalar check and retain every condition."""

        if (
            self.assumption_context.consistency_status
            is not ConsistencyStatus.CONSISTENT
        ):
            raise IdentityVerificationIndeterminateError(
                "the assumption context is not internally certified consistent."
            )
        residual = sp.simplify(self.lhs - self.rhs)
        if residual.is_zero is not True:
            if residual.is_zero is False:
                raise IdentityVerificationError(
                    f"the scalar residual is nonzero: {sp.sstr(residual)}"
                )
            raise IdentityVerificationIndeterminateError(
                "SymPy could not decide whether the scalar residual is zero."
            )
        record = SymbolicCheckRecord(
            method="sympy_simplify_difference",
            residual=residual,
        )
        return ConditionalIdentity(
            lhs=self.lhs,
            rhs=self.rhs,
            assumption_context=self.assumption_context,
            domain=self.domain,
            singular_set=self.singular_set,
            scope=self.scope,
            verification_status=(
                ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
            ),
            evidence=(record,),
            description=self.description,
            verification_method=record.method,
        )

    def substitute(
        self,
        substitutions: dict[sp.Basic, sp.Basic],
    ) -> ConditionalIdentity:
        """Substitute every semantic component or reject an invalid context."""

        if not isinstance(substitutions, dict):
            raise TypeError("substitutions must be a dict.")
        normalized = {
            sp.sympify(source): sp.sympify(target)
            for source, target in substitutions.items()
        }
        context = self.assumption_context.substitute(normalized)
        if context.consistency_status is ConsistencyStatus.INCONSISTENT:
            raise IdentityApplicabilityError(
                "substitution explicitly violates a required hypothesis."
            )
        try:
            domain = self.domain.substitute(normalized)
        except IncompatibleDomainError as exc:
            raise IdentityApplicabilityError(
                "substitution makes the declared domain incompatible."
            ) from exc
        result = ConditionalIdentity(
            lhs=self.lhs.xreplace(normalized),
            rhs=self.rhs.xreplace(normalized),
            assumption_context=context,
            domain=domain,
            singular_set=self.singular_set.substitute(normalized),
            scope=self.scope,
            verification_status=ConditionalVerificationStatus.DECLARED,
            evidence=self.evidence,
            description=self.description,
        )
        if (
            self.verification_status
            is ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
        ):
            return result.symbolically_check()
        return result

    def differentiate(self, variable: sp.Symbol) -> ConditionalIdentity:
        """Differentiate a scalar identity while retaining or adding metadata."""

        if not isinstance(variable, sp.Symbol):
            raise TypeError("variable must be a SymPy Symbol.")
        if variable != self.domain.variable:
            raise IdentityScopeError(
                "differentiation is supported only in the principal domain variable."
            )
        if (
            self.verification_status
            is not ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
        ):
            raise IdentityVerificationError(
                "differentiation requires an internally checked scalar identity."
            )
        require_applicable_identity(
            self,
            available_context=self.assumption_context,
            available_domain=self.domain,
            required_verification_status=(
                ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
            ),
        )
        lhs = sp.diff(self.lhs, variable)
        rhs = sp.diff(self.rhs, variable)
        detected = detect_singularities(
            _expression_with_both_sides(lhs, rhs),
            variable,
            assumptions=self.assumption_context,
            domain=self.domain,
        )
        result = ConditionalIdentity(
            lhs=lhs,
            rhs=rhs,
            assumption_context=self.assumption_context,
            domain=self.domain,
            singular_set=self.singular_set.union(detected),
            scope=self.scope,
            verification_status=ConditionalVerificationStatus.DECLARED,
            evidence=self.evidence,
            description=(
                f"Derivative of {self.description}"
                if self.description is not None
                else "Derivative of a conditional scalar identity"
            ),
        )
        if (
            self.verification_status
            is ConditionalVerificationStatus.SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS
        ):
            return result.symbolically_check()
        return result

    def multiply(self, other: ConditionalIdentity) -> ConditionalIdentity:
        """Multiply scalar identities with conservative metadata combination."""

        if not isinstance(other, ConditionalIdentity):
            raise TypeError("other must be a ConditionalIdentity.")
        context = self.assumption_context.combine(other.assumption_context)
        if context.consistency_status is ConsistencyStatus.INCONSISTENT:
            raise IdentityApplicabilityError(
                "conditional identities have contradictory contexts."
            )
        try:
            domain = self.domain.intersect(other.domain)
        except IncompatibleDomainError as exc:
            raise IdentityApplicabilityError(
                "conditional identities have incompatible domains."
            ) from exc
        lhs = self.lhs * other.lhs
        rhs = self.rhs * other.rhs
        detected = detect_singularities(
            _expression_with_both_sides(lhs, rhs),
            domain.variable,
            assumptions=context,
            domain=domain,
        )
        return ConditionalIdentity(
            lhs=lhs,
            rhs=rhs,
            assumption_context=context,
            domain=domain,
            singular_set=(
                self.singular_set.union(other.singular_set).union(detected)
            ),
            scope=ExactIdentityScope.SCALAR_SYMBOLIC,
            verification_status=ConditionalVerificationStatus.DECLARED,
            evidence=(*self.evidence, *other.evidence),
            description="Product of two conditional scalar identities",
        )

    def as_expr(self) -> sp.Equality:
        """Explicitly project to SymPy, deliberately losing semantic metadata."""

        return sp.Eq(self.lhs, self.rhs, evaluate=False)

    def __str__(self) -> str:
        assumptions = "; ".join(
            sp.sstr(item) for item in self.assumption_context.assumptions
        ) or "none"
        return (
            f"ConditionalIdentity({sp.sstr(self.lhs)} = {sp.sstr(self.rhs)}; "
            f"assumptions={assumptions}; domain={self.domain}; "
            f"singularities={len(self.singular_set.singularities)}; "
            f"scope={self.scope.value}; status={self.verification_status.value})"
        )


def require_applicable_identity(
    identity: ConditionalIdentity,
    *,
    available_context: AssumptionContext,
    available_domain: ComplexDomain,
    point: object | None = None,
    required_scope: ExactIdentityScope = ExactIdentityScope.SCALAR_SYMBOLIC,
    required_verification_status: ConditionalVerificationStatus = (
        ConditionalVerificationStatus.DECLARED
    ),
) -> ConditionalIdentity:
    """Return an identity only after every declared applicability check passes."""

    if not isinstance(identity, ConditionalIdentity):
        raise TypeError("identity must be a ConditionalIdentity.")
    if not isinstance(available_context, AssumptionContext):
        raise TypeError("available_context must be an AssumptionContext.")
    if not isinstance(available_domain, ComplexDomain):
        raise TypeError("available_domain must be a ComplexDomain.")
    if not isinstance(required_scope, ExactIdentityScope):
        raise TypeError("required_scope must be an ExactIdentityScope.")
    if not isinstance(
        required_verification_status,
        ConditionalVerificationStatus,
    ):
        raise TypeError(
            "required_verification_status must be a "
            "ConditionalVerificationStatus."
        )
    if required_scope is not ExactIdentityScope.SCALAR_SYMBOLIC:
        raise IdentityScopeError(
            "a scalar ConditionalIdentity cannot satisfy a stronger or different scope."
        )
    if (
        required_verification_status is not ConditionalVerificationStatus.DECLARED
        and identity.verification_status is not required_verification_status
    ):
        raise IdentityApplicabilityError(
            "the identity does not have the required verification status."
        )
    context = available_domain.assumption_context.combine(available_context)
    if not context.contains_all(identity.assumption_context):
        raise IdentityApplicabilityError(
            "the available context omits one or more required hypotheses."
        )
    if context.consistency_status is ConsistencyStatus.INCONSISTENT:
        raise IdentityApplicabilityError("the available context is inconsistent.")
    if context.consistency_status is ConsistencyStatus.UNDETERMINED:
        raise IdentityApplicabilityError(
            "the available context has undetermined consistency."
        )
    containment = available_domain.is_subset_of(identity.domain)
    if containment is not MembershipStatus.YES:
        raise IdentityApplicabilityError(
            "the available domain is broader than, or not provably contained in, "
            "the identity domain."
        )
    if point is not None:
        if available_domain.contains_point(point) is not MembershipStatus.YES:
            raise IdentityApplicabilityError(
                "the requested point is not provably in the available domain."
            )
        singular = identity.singular_set.contains_point(
            point,
            variable=identity.domain.variable,
            available_context=context,
        )
        if singular is not MembershipStatus.NO:
            raise IdentityApplicabilityError(
                "the requested point is singular or singularity status is undetermined."
            )
    else:
        avoided = identity.singular_set.avoidance_status(
            available_domain,
            available_context=context,
        )
        if avoided is not MembershipStatus.YES:
            raise IdentityApplicabilityError(
                "singularities are not provably avoided on the available domain."
            )
    return identity


@final
@dataclass(frozen=True)
class CothConditionalIdentityFamily:
    """The required conditioned scalar ``coth`` reference identities."""

    symbol: sp.Expr
    p_minus: sp.Expr
    p_plus: sp.Expr
    product: ConditionalIdentity
    symbol_derivative: ConditionalIdentity
    p_minus_derivative: ConditionalIdentity
    p_plus_derivative: ConditionalIdentity


def _checked_identity(
    lhs: sp.Expr,
    rhs: sp.Expr,
    *,
    context: AssumptionContext,
    domain: ComplexDomain,
    description: str,
) -> ConditionalIdentity:
    singular_set = detect_singularities(
        _expression_with_both_sides(lhs, rhs),
        domain.variable,
        assumptions=context,
        domain=domain,
    )
    return ConditionalIdentity(
        lhs=lhs,
        rhs=rhs,
        assumption_context=context,
        domain=domain,
        singular_set=singular_set,
        scope=ExactIdentityScope.SCALAR_SYMBOLIC,
        description=description,
    ).symbolically_check()


def build_coth_conditional_identities(
    lam: sp.Symbol,
    kappa: sp.Symbol,
) -> CothConditionalIdentityFamily:
    """Build the four required scalar identities under ``0 < kappa < 1``."""

    if not isinstance(lam, sp.Symbol) or not isinstance(kappa, sp.Symbol):
        raise TypeError("lam and kappa must be SymPy Symbols.")
    if lam == kappa:
        raise ValueError("lam and kappa must be distinct symbols.")
    context = AssumptionContext((kappa > 0, kappa < 1))
    domain = ComplexDomain.real_line(
        lam,
        assumption_context=context,
        description="real lambda line under 0 < kappa < 1",
    )
    argument = sp.pi * (lam + sp.I * kappa)
    symbol = sp.coth(argument)
    p_minus = (1 - symbol) / 2
    p_plus = (1 + symbol) / 2
    denominator = sp.sinh(argument) ** 2
    return CothConditionalIdentityFamily(
        symbol=symbol,
        p_minus=p_minus,
        p_plus=p_plus,
        product=_checked_identity(
            p_minus * p_plus,
            -1 / (4 * denominator),
            context=context,
            domain=domain,
            description="Conditional scalar product identity for p_minus and p_plus",
        ),
        symbol_derivative=_checked_identity(
            sp.diff(symbol, lam),
            -sp.pi / denominator,
            context=context,
            domain=domain,
            description="Conditional scalar derivative identity for coth",
        ),
        p_minus_derivative=_checked_identity(
            sp.diff(p_minus, lam),
            sp.pi / (2 * denominator),
            context=context,
            domain=domain,
            description="Conditional scalar derivative identity for p_minus",
        ),
        p_plus_derivative=_checked_identity(
            sp.diff(p_plus, lam),
            -sp.pi / (2 * denominator),
            context=context,
            domain=domain,
            description="Conditional scalar derivative identity for p_plus",
        ),
    )
