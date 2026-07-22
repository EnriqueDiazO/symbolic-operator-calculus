"""Fail-closed certification of a factored diagonal regularizer.

The module proves only the two ordered products supplied by a typed link
identity and an explicitly sided Mellin-core regularizer.  It never identifies
the Mellin regularizer with the full word, commutes factors, infers membership
in a common algebra, or turns the word into one Mellin PDO.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import sympy as sp

from .domains import AssumptionContext, ConsistencyStatus
from .mellin.operators import (
    MellinPseudodifferentialOperator,
    MultiplicationOperator,
    WeightedLpSpace,
)
from .operators import OperatorAtom, Product
from .schur_full_word import AuxiliaryOperator, TransportedMellinCore
from .semantics import (
    AnalyticProofObligation,
    ExactIdentity,
    ExactIdentityScope,
    FormalIdentity,
    ModCompactEquivalence,
)


IDENTITY_OPERATOR = OperatorAtom("I", kind="identity")
IDENTITY_PRODUCT = Product((IDENTITY_OPERATOR,))


class DiagonalRegularizerError(ValueError):
    """Base error for malformed or unsupported diagonal certificates."""


class DomainOrOrderError(DiagonalRegularizerError):
    """Raised when spaces, roles, or noncommutative order disagree."""


class CompactIdealMismatch(DiagonalRegularizerError):
    """Raised before differently typed compact ideals can be combined."""


class MissingRegularizerEvidence(DiagonalRegularizerError):
    """Raised when a certified or conditional relation lacks its evidence."""


class DiagonalLinkStatus(str, Enum):
    EXACT = "EXACT"
    MODULO_COMPACTS = "MODULO_COMPACTS"
    CONDITIONAL = "CONDITIONAL"
    UNVERIFIED = "UNVERIFIED"


class CoreRegularizerProperty(str, Enum):
    """Sidedness of the explicitly recorded core relations.

    ``LEFT_REGULARIZER`` names the relation with the core on the left,
    ``H R ~= I``. ``RIGHT_REGULARIZER`` names ``R H ~= I``.  The product
    fields themselves remain authoritative and remove nomenclature ambiguity.
    """

    LEFT_REGULARIZER = "LEFT_REGULARIZER"
    RIGHT_REGULARIZER = "RIGHT_REGULARIZER"
    TWO_SIDED_REGULARIZER = "TWO_SIDED_REGULARIZER"
    EXACT_INVERSE = "EXACT_INVERSE"
    CALKIN_INVERSE = "CALKIN_INVERSE"


class CoreRegularizerClaimStatus(str, Enum):
    CERTIFIED = "CERTIFIED"
    CONDITIONAL = "CONDITIONAL"
    UNVERIFIED = "UNVERIFIED"


class FactoredRegularizerStatus(str, Enum):
    CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER = (
        "CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER"
    )
    CERTIFIED_LEFT_REGULARIZER_ONLY = "CERTIFIED_LEFT_REGULARIZER_ONLY"
    CERTIFIED_RIGHT_REGULARIZER_ONLY = "CERTIFIED_RIGHT_REGULARIZER_ONLY"
    CONDITIONAL_TWO_SIDED_FACTORED_REGULARIZER = (
        "CONDITIONAL_TWO_SIDED_FACTORED_REGULARIZER"
    )
    BLOCKED_BY_LINK_IDENTITY = "BLOCKED_BY_LINK_IDENTITY"
    BLOCKED_BY_MELLIN_CORE_REGULARIZER = "BLOCKED_BY_MELLIN_CORE_REGULARIZER"
    BLOCKED_BY_DOMAIN_OR_ORDER = "BLOCKED_BY_DOMAIN_OR_ORDER"
    BLOCKED_BY_COMPACT_IDEAL_MISMATCH = "BLOCKED_BY_COMPACT_IDEAL_MISMATCH"
    BLOCKED_BY_SOURCE_GAP = "BLOCKED_BY_SOURCE_GAP"


class ProductCertificationStatus(str, Enum):
    LEFT_REGULARIZER_CERTIFIED = "LEFT_REGULARIZER_CERTIFIED"
    RIGHT_REGULARIZER_CERTIFIED = "RIGHT_REGULARIZER_CERTIFIED"
    CONDITIONAL = "CONDITIONAL"
    BLOCKED = "BLOCKED"


class SchurInsertionStatus(str, Enum):
    STRUCTURALLY_VALID = "STRUCTURALLY_VALID"


class UnprovedClaimStatus(str, Enum):
    UNPROVED = "UNPROVED"


class SchurSymbolStatus(str, Enum):
    NOT_CERTIFIED = "NOT_CERTIFIED"


def _nonempty(value: object, field_name: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{field_name} must be a non-empty str.")
    if not value.strip():
        raise ValueError(f"{field_name} must be a non-empty str.")
    return value


def _items(value: object, field_name: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{field_name} must be an iterable, not text.")
    try:
        return tuple(value)  # type: ignore[arg-type]
    except TypeError as exc:
        raise TypeError(f"{field_name} must be an iterable.") from exc


def _atom(value: object, field_name: str) -> OperatorAtom:
    if isinstance(value, OperatorAtom):
        return value
    atom = getattr(value, "atom", None)
    if not isinstance(atom, OperatorAtom):
        raise TypeError(f"{field_name} must expose an OperatorAtom as .atom.")
    return atom


def _space(value: object, field_name: str) -> WeightedLpSpace:
    domain = getattr(value, "domain", None)
    codomain = getattr(value, "codomain", None)
    if not isinstance(domain, WeightedLpSpace) or not isinstance(
        codomain, WeightedLpSpace
    ):
        raise TypeError(
            f"{field_name} must expose WeightedLpSpace domain and codomain."
        )
    if domain != codomain:
        raise DomainOrOrderError(f"{field_name} must preserve one weighted space.")
    return domain


def _require_bounded(value: object, field_name: str) -> None:
    if getattr(value, "bounded", None) is not True:
        raise DomainOrOrderError(
            f"{field_name} must be explicitly bounded before ideal propagation."
        )


def _require_same_space(
    expected: WeightedLpSpace,
    *values: tuple[object, str],
) -> None:
    for value, field_name in values:
        if _space(value, field_name) != expected:
            raise DomainOrOrderError(
                f"{field_name} acts on a different weighted space."
            )


def _provably_has_positive_zero(
    expression: sp.Expr,
    radial_variable: sp.Symbol,
) -> bool:
    """Detect only zeros SymPy certifies inside the open positive half-line."""

    try:
        zeros = sp.solveset(
            expression,
            radial_variable,
            domain=sp.Interval.open(0, sp.oo),
        )
    except (NotImplementedError, ValueError):
        return False
    return zeros.is_empty is False


def _relation_endpoints(
    relation: ExactIdentity | FormalIdentity | ModCompactEquivalence,
) -> tuple[object, object]:
    return relation.left, relation.right


def _relation_matches(
    relation: ExactIdentity | FormalIdentity | ModCompactEquivalence,
    left: object,
    right: object,
) -> bool:
    return _relation_endpoints(relation) == (left, right)


def _prepend_append(
    expression: object,
    left_factors: tuple[object, ...],
    right_factors: tuple[object, ...],
) -> Product:
    if isinstance(expression, OperatorAtom):
        middle = (expression,)
    elif isinstance(expression, Product):
        middle = expression.factors
    else:
        raise TypeError("ideal propagation endpoints must be operator products.")
    return Product(
        tuple(_atom(item, "left factor") for item in left_factors)
        + middle
        + tuple(_atom(item, "right factor") for item in right_factors)
    )


@dataclass(frozen=True)
class CompactIdeal:
    """One explicitly bilateral compact ideal on one scalar operator space."""

    label: str
    space: WeightedLpSpace
    bilateral: bool
    source: str
    evidence: tuple[object, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "label", _nonempty(self.label, "label"))
        if not isinstance(self.space, WeightedLpSpace):
            raise TypeError("space must be a WeightedLpSpace.")
        if not isinstance(self.bilateral, bool):
            raise TypeError("bilateral must be a bool.")
        source = _nonempty(self.source, "source")
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise MissingRegularizerEvidence(
                "a compact ideal requires source evidence."
            )
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)

    @property
    def identity_key(self) -> tuple[str, WeightedLpSpace]:
        """Return the fields that identify the ideal independently of prose."""

        return self.label, self.space

    def require_same(self, other: CompactIdeal) -> None:
        if not isinstance(other, CompactIdeal):
            raise CompactIdealMismatch("both ideals must be CompactIdeal objects.")
        if self.identity_key != other.identity_key:
            raise CompactIdealMismatch(
                "compact ideals have different labels or weighted spaces."
            )


@dataclass(frozen=True)
class DiagonalOperator:
    """A bounded diagonal operator retaining its authoritative AST atom."""

    atom: OperatorAtom
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    source: str
    evidence: tuple[object, ...]
    bounded: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        _space(self, "diagonal operator")
        if not isinstance(self.bounded, bool):
            raise TypeError("bounded must be a bool.")
        source = _nonempty(self.source, "source")
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise MissingRegularizerEvidence(
                "a diagonal operator requires definition evidence."
            )
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)


@dataclass(frozen=True)
class MellinCoreOperator:
    """An article-level core with an exact transported Mellin realization."""

    atom: OperatorAtom
    mellin_operator: MellinPseudodifferentialOperator
    phi: OperatorAtom
    phi_inverse: OperatorAtom
    source: str
    evidence: tuple[object, ...]
    alias_of: str
    exact_transport: ExactIdentity = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.mellin_operator, MellinPseudodifferentialOperator):
            raise TypeError("mellin_operator must be a Mellin PDO.")
        if not isinstance(self.phi, OperatorAtom) or not isinstance(
            self.phi_inverse, OperatorAtom
        ):
            raise TypeError("phi and phi_inverse must be OperatorAtom objects.")
        _space(self.mellin_operator, "Mellin core realization")
        source = _nonempty(self.source, "source")
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise MissingRegularizerEvidence(
                "an exact Mellin-core transport requires evidence."
            )
        alias = _nonempty(self.alias_of, "alias_of")
        relation = ExactIdentity(
            Product((self.phi, self.atom, self.phi_inverse)),
            Product((self.mellin_operator.atom,)),
            evidence=(source, evidence, alias),
            scope=ExactIdentityScope.OPERATORIAL,
        )
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "alias_of", alias)
        object.__setattr__(self, "exact_transport", relation)

    @property
    def domain(self) -> WeightedLpSpace:
        return self.mellin_operator.domain

    @property
    def codomain(self) -> WeightedLpSpace:
        return self.mellin_operator.codomain

    @property
    def bounded(self) -> bool:
        return self.mellin_operator.bounded

    @property
    def symbol(self) -> object:
        return self.mellin_operator.symbol

    @property
    def symbol_class(self) -> str:
        return self.mellin_operator.symbol_class


@dataclass(frozen=True)
class InvertibleAuxiliaryOperator:
    """A typed auxiliary and its source-backed exact inverse."""

    operator: AuxiliaryOperator
    inverse: AuxiliaryOperator
    source: str
    evidence: tuple[object, ...]
    operator_times_inverse: ExactIdentity = field(init=False)
    inverse_times_operator: ExactIdentity = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.operator, AuxiliaryOperator) or not isinstance(
            self.inverse, AuxiliaryOperator
        ):
            raise TypeError(
                "operator and inverse must be AuxiliaryOperator objects; "
                "a bare shift is not an auxiliary inverse."
            )
        space = _space(self.operator, "auxiliary")
        _require_same_space(space, (self.inverse, "auxiliary inverse"))
        _require_bounded(self.operator, "auxiliary")
        _require_bounded(self.inverse, "auxiliary inverse")
        if self.operator.atom == self.inverse.atom:
            raise DomainOrOrderError(
                "an auxiliary and its inverse need distinct authoritative atoms."
            )
        source = _nonempty(self.source, "source")
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise MissingRegularizerEvidence(
                "exact auxiliary invertibility requires theorem evidence."
            )
        forward = ExactIdentity(
            Product((self.operator.atom, self.inverse.atom)),
            IDENTITY_PRODUCT,
            evidence=(source, evidence),
            scope=ExactIdentityScope.OPERATORIAL,
        )
        backward = ExactIdentity(
            Product((self.inverse.atom, self.operator.atom)),
            IDENTITY_PRODUCT,
            evidence=(source, evidence),
            scope=ExactIdentityScope.OPERATORIAL,
        )
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "operator_times_inverse", forward)
        object.__setattr__(self, "inverse_times_operator", backward)

    @property
    def atom(self) -> OperatorAtom:
        return self.operator.atom

    @property
    def inverse_atom(self) -> OperatorAtom:
        return self.inverse.atom

    @property
    def domain(self) -> WeightedLpSpace:
        return self.operator.domain

    @property
    def codomain(self) -> WeightedLpSpace:
        return self.operator.codomain

    @property
    def bounded(self) -> bool:
        return self.operator.bounded and self.inverse.bounded


@dataclass(frozen=True)
class InvertibleMultiplicationOperator:
    """A nonvanishing multiplier paired with its exact multiplier inverse."""

    operator: MultiplicationOperator
    inverse: MultiplicationOperator
    source: str
    evidence: tuple[object, ...]
    nonvanishing_evidence: tuple[object, ...]
    operator_times_inverse: ExactIdentity = field(init=False)
    inverse_times_operator: ExactIdentity = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.operator, MultiplicationOperator) or not isinstance(
            self.inverse, MultiplicationOperator
        ):
            raise TypeError(
                "operator and inverse must be MultiplicationOperator objects."
            )
        space = _space(self.operator, "multiplier")
        _require_same_space(space, (self.inverse, "inverse multiplier"))
        _require_bounded(self.operator, "multiplier")
        _require_bounded(self.inverse, "inverse multiplier")
        if self.operator.radial_variable != self.inverse.radial_variable:
            raise DomainOrOrderError(
                "a multiplier and its inverse need the same radial variable."
            )
        if self.operator.atom == self.inverse.atom:
            raise DomainOrOrderError(
                "a multiplier and its inverse need distinct authoritative atoms."
            )
        source = _nonempty(self.source, "source")
        evidence = _items(self.evidence, "evidence")
        nonvanishing = _items(
            self.nonvanishing_evidence,
            "nonvanishing_evidence",
        )
        if not evidence or not nonvanishing:
            raise MissingRegularizerEvidence(
                "an exact multiplier inverse requires identity and "
                "nonvanishing evidence."
            )
        if self.operator.symbol.is_zero is True or _provably_has_positive_zero(
            self.operator.symbol,
            self.operator.radial_variable,
        ):
            raise DomainOrOrderError(
                "the multiplier symbol has a certified zero on the positive half-line."
            )
        inverse_check = sp.simplify(
            self.operator.symbol * self.inverse.symbol - sp.Integer(1)
        )
        if inverse_check != 0:
            raise DomainOrOrderError(
                "the supplied multiplication symbols are not exact reciprocals."
            )
        forward = ExactIdentity(
            Product((self.operator.atom, self.inverse.atom)),
            IDENTITY_PRODUCT,
            evidence=(source, evidence, nonvanishing),
            scope=ExactIdentityScope.OPERATORIAL,
        )
        backward = ExactIdentity(
            Product((self.inverse.atom, self.operator.atom)),
            IDENTITY_PRODUCT,
            evidence=(source, evidence, nonvanishing),
            scope=ExactIdentityScope.OPERATORIAL,
        )
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "nonvanishing_evidence", nonvanishing)
        object.__setattr__(self, "operator_times_inverse", forward)
        object.__setattr__(self, "inverse_times_operator", backward)

    @property
    def atom(self) -> OperatorAtom:
        return self.operator.atom

    @property
    def inverse_atom(self) -> OperatorAtom:
        return self.inverse.atom

    @property
    def domain(self) -> WeightedLpSpace:
        return self.operator.domain

    @property
    def codomain(self) -> WeightedLpSpace:
        return self.operator.codomain

    @property
    def bounded(self) -> bool:
        return self.operator.bounded and self.inverse.bounded


@dataclass(frozen=True)
class DiagonalLinkIdentity:
    """One oriented link supplied as typed structure, never reconstructed."""

    diagonal_operator: DiagonalOperator
    left_interface: InvertibleAuxiliaryOperator
    mellin_operator: MellinCoreOperator
    right_interface: InvertibleMultiplicationOperator
    right_interface_inverse: MultiplicationOperator
    left_expression: Product
    right_expression: Product
    relation: ExactIdentity | FormalIdentity | ModCompactEquivalence | None
    compact_ideal: CompactIdeal | None
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    assumptions: AssumptionContext
    evidence: tuple[object, ...]
    source: str
    status: DiagonalLinkStatus

    def __post_init__(self) -> None:
        if not isinstance(self.diagonal_operator, DiagonalOperator):
            raise TypeError("diagonal_operator must be a DiagonalOperator.")
        if not isinstance(self.left_interface, InvertibleAuxiliaryOperator):
            raise TypeError("left_interface must be an invertible auxiliary.")
        if not isinstance(self.mellin_operator, MellinCoreOperator):
            raise TypeError("mellin_operator must be a MellinCoreOperator.")
        if not isinstance(self.right_interface, InvertibleMultiplicationOperator):
            raise TypeError("right_interface must be an invertible multiplier.")
        if not isinstance(
            self.right_interface_inverse, MultiplicationOperator
        ) or self.right_interface_inverse != self.right_interface.inverse:
            raise DomainOrOrderError(
                "right_interface_inverse must be the typed inverse of right_interface."
            )
        if not isinstance(self.left_expression, Product) or not isinstance(
            self.right_expression, Product
        ):
            raise TypeError("link endpoints must be Product objects.")
        expected_left = Product(
            (
                self.right_interface.inverse_atom,
                self.diagonal_operator.atom,
                self.left_interface.atom,
            )
        )
        expected_right = Product((self.mellin_operator.atom,))
        if self.left_expression != expected_left or self.right_expression != expected_right:
            raise DomainOrOrderError(
                "link order must be Z_inverse A T = H with the supplied typed roles."
            )
        role_atoms = (
            self.diagonal_operator.atom,
            self.left_interface.atom,
            self.mellin_operator.atom,
            self.right_interface.atom,
            self.right_interface.inverse_atom,
        )
        if len(set(role_atoms)) != len(role_atoms):
            raise DomainOrOrderError("link roles contain duplicated interfaces.")
        if not isinstance(self.domain, WeightedLpSpace) or not isinstance(
            self.codomain, WeightedLpSpace
        ):
            raise TypeError("link domain/codomain must be WeightedLpSpace objects.")
        if self.domain != self.codomain:
            raise DomainOrOrderError("the link must preserve one weighted space.")
        _require_same_space(
            self.domain,
            (self.diagonal_operator, "diagonal operator"),
            (self.left_interface, "left interface"),
            (self.mellin_operator, "Mellin core"),
            (self.right_interface, "right interface"),
            (self.right_interface_inverse, "right inverse"),
        )
        for value, label in (
            (self.diagonal_operator, "diagonal operator"),
            (self.left_interface, "left interface"),
            (self.mellin_operator, "Mellin core"),
            (self.right_interface, "right interface"),
            (self.right_interface_inverse, "right inverse"),
        ):
            _require_bounded(value, label)
        if not isinstance(self.assumptions, AssumptionContext):
            raise TypeError("assumptions must be an AssumptionContext.")
        if self.assumptions.consistency_status is ConsistencyStatus.INCONSISTENT:
            raise DomainOrOrderError("link assumptions are inconsistent.")
        if not isinstance(self.status, DiagonalLinkStatus):
            raise TypeError("status must be a DiagonalLinkStatus.")
        source = _nonempty(self.source, "source")
        evidence = _items(self.evidence, "evidence")
        if self.status is not DiagonalLinkStatus.UNVERIFIED and not evidence:
            raise MissingRegularizerEvidence(
                "a non-unverified link requires explicit evidence."
            )
        if self.status is DiagonalLinkStatus.EXACT:
            if not isinstance(self.relation, ExactIdentity):
                raise DomainOrOrderError("an exact link requires ExactIdentity.")
            if self.compact_ideal is not None:
                raise CompactIdealMismatch("an exact link must not name an ideal.")
        elif self.status is DiagonalLinkStatus.MODULO_COMPACTS:
            if not isinstance(self.relation, ModCompactEquivalence):
                raise DomainOrOrderError(
                    "a modulo-compact link requires ModCompactEquivalence."
                )
            if not isinstance(self.compact_ideal, CompactIdeal):
                raise CompactIdealMismatch(
                    "a modulo-compact link requires a typed compact ideal."
                )
            if not self.compact_ideal.bilateral:
                raise CompactIdealMismatch("the link ideal is not bilateral.")
            if self.compact_ideal.space != self.domain:
                raise CompactIdealMismatch("the link ideal uses another space.")
            if self.relation.compact_ideal != self.compact_ideal:
                raise CompactIdealMismatch(
                    "the semantic link relation names another ideal."
                )
        elif self.status is DiagonalLinkStatus.CONDITIONAL:
            if not isinstance(self.relation, FormalIdentity):
                raise DomainOrOrderError(
                    "a conditional link requires a FormalIdentity."
                )
            if self.assumptions.is_empty:
                raise MissingRegularizerEvidence(
                    "a conditional link requires explicit assumptions."
                )
            if self.compact_ideal is not None:
                if not self.compact_ideal.bilateral:
                    raise CompactIdealMismatch(
                        "a conditional link ideal must be bilateral."
                    )
                if self.compact_ideal.space != self.domain:
                    raise CompactIdealMismatch(
                        "the conditional link ideal uses another space."
                    )
        elif self.relation is not None:
            raise DomainOrOrderError("an unverified link cannot carry a relation.")
        if self.relation is not None and not _relation_matches(
            self.relation,
            self.left_expression,
            self.right_expression,
        ):
            raise DomainOrOrderError(
                "link semantic endpoints do not match the declared expressions."
            )
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)


@dataclass(frozen=True)
class MellinCoreRegularizerCertificate:
    """Oriented regularizer evidence for a transported Mellin core."""

    mellin_operator: MellinCoreOperator
    mellin_regularizer: TransportedMellinCore
    property: CoreRegularizerProperty
    claim_status: CoreRegularizerClaimStatus
    operator_then_regularizer_relation: (
        ExactIdentity | FormalIdentity | ModCompactEquivalence | None
    )
    regularizer_then_operator_relation: (
        ExactIdentity | FormalIdentity | ModCompactEquivalence | None
    )
    compact_ideal: CompactIdeal | None
    assumptions: AssumptionContext
    source: str
    evidence: tuple[object, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.mellin_operator, MellinCoreOperator):
            raise TypeError("mellin_operator must be a MellinCoreOperator.")
        if not isinstance(self.mellin_regularizer, TransportedMellinCore):
            raise TypeError("mellin_regularizer must be a TransportedMellinCore.")
        if not isinstance(self.property, CoreRegularizerProperty):
            raise TypeError("property must be a CoreRegularizerProperty.")
        if not isinstance(self.claim_status, CoreRegularizerClaimStatus):
            raise TypeError("claim_status must be a CoreRegularizerClaimStatus.")
        space = _space(self.mellin_operator, "Mellin core")
        _require_same_space(
            space,
            (self.mellin_regularizer, "Mellin regularizer"),
        )
        _require_bounded(self.mellin_operator, "Mellin core")
        representation = self.mellin_regularizer.representation
        if representation.mellin_operator.bounded is not True:
            raise DomainOrOrderError("Mellin regularizer must be bounded.")
        target = representation.regularizer.target
        if target != self.mellin_operator.atom and target != self.mellin_operator:
            raise DomainOrOrderError(
                "the regularizer target is not the supplied Mellin core."
            )
        if not isinstance(self.assumptions, AssumptionContext):
            raise TypeError("assumptions must be an AssumptionContext.")
        if self.assumptions.consistency_status is ConsistencyStatus.INCONSISTENT:
            raise DomainOrOrderError("core assumptions are inconsistent.")
        source = _nonempty(self.source, "source")
        evidence = _items(self.evidence, "evidence")
        if self.claim_status is not CoreRegularizerClaimStatus.UNVERIFIED and not evidence:
            raise MissingRegularizerEvidence(
                "a certified or conditional regularizer requires evidence."
            )
        if (
            self.claim_status is CoreRegularizerClaimStatus.CONDITIONAL
            and self.assumptions.is_empty
        ):
            raise MissingRegularizerEvidence(
                "a conditional regularizer requires explicit assumptions."
            )
        needs_left = self.property in {
            CoreRegularizerProperty.LEFT_REGULARIZER,
            CoreRegularizerProperty.TWO_SIDED_REGULARIZER,
            CoreRegularizerProperty.EXACT_INVERSE,
            CoreRegularizerProperty.CALKIN_INVERSE,
        }
        needs_right = self.property in {
            CoreRegularizerProperty.RIGHT_REGULARIZER,
            CoreRegularizerProperty.TWO_SIDED_REGULARIZER,
            CoreRegularizerProperty.EXACT_INVERSE,
            CoreRegularizerProperty.CALKIN_INVERSE,
        }
        left_expression = Product(
            (self.mellin_operator.atom, self.mellin_regularizer.atom)
        )
        right_expression = Product(
            (self.mellin_regularizer.atom, self.mellin_operator.atom)
        )
        if self.claim_status is CoreRegularizerClaimStatus.UNVERIFIED:
            if (
                self.operator_then_regularizer_relation is not None
                or self.regularizer_then_operator_relation is not None
            ):
                raise DomainOrOrderError(
                    "an unverified core regularizer cannot carry relations."
                )
        else:
            if needs_left != (
                self.operator_then_regularizer_relation is not None
            ) or needs_right != (
                self.regularizer_then_operator_relation is not None
            ):
                raise DomainOrOrderError(
                    "core relation availability does not match declared sidedness."
                )
        for relation, expected_left in (
            (self.operator_then_regularizer_relation, left_expression),
            (self.regularizer_then_operator_relation, right_expression),
        ):
            if relation is None:
                continue
            if not _relation_matches(relation, expected_left, IDENTITY_PRODUCT):
                raise DomainOrOrderError(
                    "a core relation has incorrect order or endpoints."
                )
            if self.claim_status is CoreRegularizerClaimStatus.CONDITIONAL:
                if not isinstance(relation, FormalIdentity):
                    raise DomainOrOrderError(
                        "conditional core relations must be FormalIdentity objects."
                    )
            elif self.property is CoreRegularizerProperty.EXACT_INVERSE:
                if not isinstance(relation, ExactIdentity):
                    raise DomainOrOrderError(
                        "EXACT_INVERSE requires exact relations on both sides."
                    )
            elif not isinstance(relation, ModCompactEquivalence):
                raise DomainOrOrderError(
                    "certified non-exact regularizers require modulo-compact relations."
                )
        if self.property is CoreRegularizerProperty.EXACT_INVERSE:
            if self.compact_ideal is not None:
                raise CompactIdealMismatch("an exact inverse must not name an ideal.")
        elif self.claim_status is not CoreRegularizerClaimStatus.UNVERIFIED:
            if not isinstance(self.compact_ideal, CompactIdeal):
                raise CompactIdealMismatch(
                    "a non-exact regularizer requires a typed compact ideal."
                )
            if not self.compact_ideal.bilateral:
                raise CompactIdealMismatch("the regularizer ideal is not bilateral.")
            if self.compact_ideal.space != space:
                raise CompactIdealMismatch(
                    "the core regularizer ideal uses another weighted space."
                )
            for relation in (
                self.operator_then_regularizer_relation,
                self.regularizer_then_operator_relation,
            ):
                if isinstance(relation, ModCompactEquivalence) and (
                    relation.compact_ideal != self.compact_ideal
                ):
                    raise CompactIdealMismatch(
                        "a core relation names a different compact ideal."
                    )
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)

    @property
    def atom(self) -> OperatorAtom:
        return self.mellin_regularizer.atom

    @property
    def domain(self) -> WeightedLpSpace:
        return self.mellin_regularizer.domain

    @property
    def codomain(self) -> WeightedLpSpace:
        return self.mellin_regularizer.codomain

    @property
    def bounded(self) -> bool:
        return self.mellin_regularizer.representation.mellin_operator.bounded

    @property
    def left_relation(
        self,
    ) -> ExactIdentity | FormalIdentity | ModCompactEquivalence | None:
        return self.operator_then_regularizer_relation

    @property
    def right_relation(
        self,
    ) -> ExactIdentity | FormalIdentity | ModCompactEquivalence | None:
        return self.regularizer_then_operator_relation


def make_mellin_core_regularizer_certificate(
    *,
    mellin_operator: MellinCoreOperator,
    mellin_regularizer: TransportedMellinCore,
    property: CoreRegularizerProperty,
    claim_status: CoreRegularizerClaimStatus,
    compact_ideal: CompactIdeal | None,
    assumptions: AssumptionContext,
    source: str,
    evidence: tuple[object, ...],
) -> MellinCoreRegularizerCertificate:
    """Build exactly the oriented relations declared by ``property``."""

    if not isinstance(property, CoreRegularizerProperty):
        raise TypeError("property must be a CoreRegularizerProperty.")
    if not isinstance(claim_status, CoreRegularizerClaimStatus):
        raise TypeError("claim_status must be a CoreRegularizerClaimStatus.")
    left_needed = property in {
        CoreRegularizerProperty.LEFT_REGULARIZER,
        CoreRegularizerProperty.TWO_SIDED_REGULARIZER,
        CoreRegularizerProperty.EXACT_INVERSE,
        CoreRegularizerProperty.CALKIN_INVERSE,
    }
    right_needed = property in {
        CoreRegularizerProperty.RIGHT_REGULARIZER,
        CoreRegularizerProperty.TWO_SIDED_REGULARIZER,
        CoreRegularizerProperty.EXACT_INVERSE,
        CoreRegularizerProperty.CALKIN_INVERSE,
    }
    left_expression = Product((mellin_operator.atom, mellin_regularizer.atom))
    right_expression = Product((mellin_regularizer.atom, mellin_operator.atom))

    def relation(expression: Product) -> (
        ExactIdentity | FormalIdentity | ModCompactEquivalence | None
    ):
        if claim_status is CoreRegularizerClaimStatus.UNVERIFIED:
            return None
        if claim_status is CoreRegularizerClaimStatus.CONDITIONAL:
            return FormalIdentity(
                expression,
                IDENTITY_PRODUCT,
                justification=(source, assumptions, evidence, compact_ideal),
            )
        if property is CoreRegularizerProperty.EXACT_INVERSE:
            return ExactIdentity(
                expression,
                IDENTITY_PRODUCT,
                evidence=(source, evidence),
                scope=ExactIdentityScope.OPERATORIAL,
                hypotheses=assumptions.assumptions,
            )
        return ModCompactEquivalence(
            expression,
            IDENTITY_PRODUCT,
            space=mellin_operator.domain,
            compact_ideal=compact_ideal,
            evidence=(source, evidence, assumptions),
        )

    return MellinCoreRegularizerCertificate(
        mellin_operator=mellin_operator,
        mellin_regularizer=mellin_regularizer,
        property=property,
        claim_status=claim_status,
        operator_then_regularizer_relation=(
            relation(left_expression) if left_needed else None
        ),
        regularizer_then_operator_relation=(
            relation(right_expression) if right_needed else None
        ),
        compact_ideal=compact_ideal,
        assumptions=assumptions,
        source=source,
        evidence=evidence,
    )


def propagate_compact_ideal(
    relation: ModCompactEquivalence,
    *,
    left_factors: tuple[object, ...],
    right_factors: tuple[object, ...],
    compact_ideal: CompactIdeal,
    evidence: tuple[object, ...],
) -> ModCompactEquivalence:
    """Implement ``X-Y in K => LXR-LYR in K`` with typed guards."""

    if not isinstance(relation, ModCompactEquivalence):
        raise TypeError("relation must be a ModCompactEquivalence.")
    if not isinstance(compact_ideal, CompactIdeal):
        raise TypeError("compact_ideal must be a CompactIdeal.")
    if not compact_ideal.bilateral:
        raise CompactIdealMismatch("compact-ideal propagation needs a bilateral ideal.")
    if relation.compact_ideal != compact_ideal:
        raise CompactIdealMismatch(
            "the relation and propagation use different compact ideals."
        )
    if relation.space != compact_ideal.space:
        raise CompactIdealMismatch(
            "the relation and compact ideal use different spaces."
        )
    left = _items(left_factors, "left_factors")
    right = _items(right_factors, "right_factors")
    supplied_evidence = _items(evidence, "evidence")
    if not supplied_evidence:
        raise MissingRegularizerEvidence(
            "compact-ideal propagation requires explicit evidence."
        )
    for index, factor in enumerate((*left, *right)):
        _require_bounded(factor, f"exterior factor {index}")
        if _space(factor, f"exterior factor {index}") != compact_ideal.space:
            raise CompactIdealMismatch(
                "an exterior factor acts on a different weighted space."
            )
    return ModCompactEquivalence(
        _prepend_append(relation.left, left, right),
        _prepend_append(relation.right, left, right),
        space=compact_ideal.space,
        compact_ideal=compact_ideal,
        evidence=(
            relation,
            supplied_evidence,
            "bilateral ideal property: L(X-Y)R is compact",
        ),
    )


@dataclass(frozen=True)
class FactoredDiagonalRegularizer:
    """The complete ordered word, retained as one structured Schur object."""

    atom: OperatorAtom
    left_interface: InvertibleAuxiliaryOperator
    mellin_regularizer: MellinCoreRegularizerCertificate
    right_interface_inverse: MultiplicationOperator
    ast_product: Product = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.left_interface, InvertibleAuxiliaryOperator):
            raise TypeError("left_interface must be an invertible auxiliary.")
        if not isinstance(
            self.mellin_regularizer, MellinCoreRegularizerCertificate
        ):
            raise TypeError("mellin_regularizer must be a core certificate.")
        if not isinstance(self.right_interface_inverse, MultiplicationOperator):
            raise TypeError("right_interface_inverse must be a multiplier.")
        space = _space(self.left_interface, "left interface")
        _require_same_space(
            space,
            (self.mellin_regularizer, "Mellin regularizer"),
            (self.right_interface_inverse, "right inverse"),
        )
        object.__setattr__(
            self,
            "ast_product",
            Product(
                (
                    self.left_interface.atom,
                    self.mellin_regularizer.atom,
                    self.right_interface_inverse.atom,
                )
            ),
        )

    @property
    def domain(self) -> WeightedLpSpace:
        return self.left_interface.domain

    @property
    def codomain(self) -> WeightedLpSpace:
        return self.left_interface.codomain

    @property
    def bounded(self) -> bool:
        return (
            self.left_interface.bounded
            and self.mellin_regularizer.bounded
            and self.right_interface_inverse.bounded
        )


@dataclass(frozen=True)
class RegularizerReductionStep:
    step: str
    input_expression: Product
    rule: str
    output_expression: Product
    relation: ExactIdentity | FormalIdentity | ModCompactEquivalence
    evidence: tuple[object, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "step", _nonempty(self.step, "step"))
        object.__setattr__(self, "rule", _nonempty(self.rule, "rule"))
        if not isinstance(self.input_expression, Product) or not isinstance(
            self.output_expression, Product
        ):
            raise TypeError("trace expressions must be Product objects.")
        if not isinstance(
            self.relation,
            (ExactIdentity, FormalIdentity, ModCompactEquivalence),
        ):
            raise TypeError("trace relation has an unsupported semantic type.")
        if not _relation_matches(
            self.relation,
            self.input_expression,
            self.output_expression,
        ):
            raise DomainOrOrderError(
                "trace relation endpoints must match its input and output."
            )
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise MissingRegularizerEvidence("every trace step requires evidence.")
        object.__setattr__(self, "evidence", evidence)


@dataclass(frozen=True)
class FactoredTwoSidedRegularizerCertificate:
    operator: DiagonalOperator
    candidate_regularizer: FactoredDiagonalRegularizer
    left_product: Product
    right_product: Product
    left_reduction_trace: tuple[RegularizerReductionStep, ...]
    right_reduction_trace: tuple[RegularizerReductionStep, ...]
    left_relation: ExactIdentity | FormalIdentity | ModCompactEquivalence | None
    right_relation: ExactIdentity | FormalIdentity | ModCompactEquivalence | None
    compact_ideal: CompactIdeal | None
    status: FactoredRegularizerStatus
    left_status: ProductCertificationStatus
    right_status: ProductCertificationStatus
    assumptions: AssumptionContext
    evidence: tuple[object, ...]
    proof_obligations: tuple[AnalyticProofObligation, ...]
    blocking_factors: tuple[object, ...]
    latex: str

    def __post_init__(self) -> None:
        if not isinstance(self.operator, DiagonalOperator):
            raise TypeError("operator must be a DiagonalOperator.")
        if not isinstance(self.candidate_regularizer, FactoredDiagonalRegularizer):
            raise TypeError(
                "candidate_regularizer must be a FactoredDiagonalRegularizer."
            )
        expected_left = Product(
            (self.operator.atom, *self.candidate_regularizer.ast_product.factors)
        )
        expected_right = Product(
            (*self.candidate_regularizer.ast_product.factors, self.operator.atom)
        )
        if self.left_product != expected_left or self.right_product != expected_right:
            raise DomainOrOrderError("certificate products lost the candidate order.")
        if not isinstance(self.status, FactoredRegularizerStatus):
            raise TypeError("status must be a FactoredRegularizerStatus.")
        if not isinstance(self.left_status, ProductCertificationStatus) or not isinstance(
            self.right_status, ProductCertificationStatus
        ):
            raise TypeError("product statuses have invalid types.")
        for trace in (self.left_reduction_trace, self.right_reduction_trace):
            if not all(isinstance(item, RegularizerReductionStep) for item in trace):
                raise TypeError("reduction traces contain invalid steps.")
        if self.left_relation is not None and not _relation_matches(
            self.left_relation, self.left_product, IDENTITY_PRODUCT
        ):
            raise DomainOrOrderError("left conclusion has incorrect endpoints.")
        if self.right_relation is not None and not _relation_matches(
            self.right_relation, self.right_product, IDENTITY_PRODUCT
        ):
            raise DomainOrOrderError("right conclusion has incorrect endpoints.")
        if not isinstance(self.assumptions, AssumptionContext):
            raise TypeError("assumptions must be an AssumptionContext.")
        if not all(
            isinstance(item, AnalyticProofObligation)
            for item in self.proof_obligations
        ):
            raise TypeError("proof_obligations contain an invalid item.")
        object.__setattr__(
            self,
            "left_reduction_trace",
            _items(self.left_reduction_trace, "left_reduction_trace"),
        )
        object.__setattr__(
            self,
            "right_reduction_trace",
            _items(self.right_reduction_trace, "right_reduction_trace"),
        )
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))
        object.__setattr__(
            self,
            "blocking_factors",
            _items(self.blocking_factors, "blocking_factors"),
        )
        object.__setattr__(self, "latex", _nonempty(self.latex, "latex"))

    @property
    def original_expression(self) -> tuple[Product, Product]:
        """Return both independently certified original product ASTs."""

        return self.left_product, self.right_product


def _open_schur_obligations() -> tuple[AnalyticProofObligation, ...]:
    return (
        AnalyticProofObligation(
            "full_regularizer_algebra_membership",
            "Prove common-algebra membership for the complete factored regularizer.",
        ),
        AnalyticProofObligation(
            "single_mellin_pdo_representation",
            "Prove a single Mellin-PDO representation before making that claim.",
        ),
        AnalyticProofObligation(
            "cutoff_replacement_mod_compacts",
            "Prove every cutoff replacement used in the surrounding Schur word.",
        ),
        AnalyticProofObligation(
            "wh_mellin_wh_sandwich_membership",
            "Prove membership of the ordered WH--Mellin--WH sandwich.",
        ),
        AnalyticProofObligation(
            "fredholm_algebra_membership",
            "Prove membership in an algebra carrying the Fredholm calculus.",
        ),
        AnalyticProofObligation(
            "schur_correction_symbol",
            "Construct a Schur-correction symbol only after membership is proved.",
        ),
    )


def _final_relation(
    original: Product,
    *,
    exact: bool,
    conditional: bool,
    compact_ideal: CompactIdeal | None,
    assumptions: AssumptionContext,
    evidence: tuple[object, ...],
) -> ExactIdentity | FormalIdentity | ModCompactEquivalence:
    if conditional:
        return FormalIdentity(
            original,
            IDENTITY_PRODUCT,
            justification=(evidence, assumptions, compact_ideal),
        )
    if exact:
        return ExactIdentity(
            original,
            IDENTITY_PRODUCT,
            evidence=evidence,
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=assumptions.assumptions,
        )
    if compact_ideal is None:
        raise CompactIdealMismatch(
            "a modulo-compact conclusion requires a common compact ideal."
        )
    return ModCompactEquivalence(
        original,
        IDENTITY_PRODUCT,
        space=compact_ideal.space,
        compact_ideal=compact_ideal,
        evidence=evidence,
    )


def _intermediate_relation(
    left: Product,
    right: Product,
    *,
    link: DiagonalLinkIdentity,
    evidence: tuple[object, ...],
) -> ExactIdentity | FormalIdentity | ModCompactEquivalence:
    if link.status is DiagonalLinkStatus.EXACT:
        return ExactIdentity(
            left,
            right,
            evidence=evidence,
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=link.assumptions.assumptions,
        )
    if link.status is DiagonalLinkStatus.MODULO_COMPACTS:
        return ModCompactEquivalence(
            left,
            right,
            space=link.domain,
            compact_ideal=link.compact_ideal,
            evidence=evidence,
        )
    return FormalIdentity(left, right, justification=evidence)


def _propagate_core_relation(
    relation: ExactIdentity | FormalIdentity | ModCompactEquivalence,
    *,
    left_factor: object,
    right_factor: object,
    compact_ideal: CompactIdeal | None,
    assumptions: AssumptionContext,
    evidence: tuple[object, ...],
) -> ExactIdentity | FormalIdentity | ModCompactEquivalence:
    output_left = _prepend_append(relation.left, (left_factor,), (right_factor,))
    output_right = _prepend_append(relation.right, (left_factor,), (right_factor,))
    if isinstance(relation, ModCompactEquivalence):
        if compact_ideal is None:
            raise CompactIdealMismatch("core propagation lost its compact ideal.")
        return propagate_compact_ideal(
            relation,
            left_factors=(left_factor,),
            right_factors=(right_factor,),
            compact_ideal=compact_ideal,
            evidence=evidence,
        )
    if isinstance(relation, ExactIdentity):
        return ExactIdentity(
            output_left,
            output_right,
            evidence=(relation, evidence),
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=assumptions.assumptions,
        )
    return FormalIdentity(
        output_left,
        output_right,
        justification=(relation, evidence, assumptions),
    )


def _derive_left_product(
    diagonal: DiagonalOperator,
    candidate: FactoredDiagonalRegularizer,
    left_interface: InvertibleAuxiliaryOperator,
    core: MellinCoreOperator,
    regularizer: MellinCoreRegularizerCertificate,
    right_interface: InvertibleMultiplicationOperator,
    link: DiagonalLinkIdentity,
    assumptions: AssumptionContext,
    compact_ideal: CompactIdeal | None,
) -> tuple[
    tuple[RegularizerReductionStep, ...],
    ExactIdentity | FormalIdentity | ModCompactEquivalence | None,
]:
    core_relation = regularizer.left_relation
    if core_relation is None:
        return (), None
    original = Product((diagonal.atom, *candidate.ast_product.factors))
    linked = Product(
        (
            right_interface.atom,
            core.atom,
            regularizer.atom,
            right_interface.inverse_atom,
        )
    )
    link_relation = _intermediate_relation(
        original,
        linked,
        link=link,
        evidence=(
            link.relation,
            right_interface.inverse_times_operator,
            "left-multiply the link by Z and associate; no commutation",
        ),
    )
    step_one = RegularizerReductionStep(
        "left-1-link",
        original,
        "apply Z_inverse A T = H and cancel Z Z_inverse exactly",
        linked,
        link_relation,
        (link, right_interface.inverse_times_operator),
    )
    regularized = Product(
        (
            right_interface.atom,
            IDENTITY_OPERATOR,
            right_interface.inverse_atom,
        )
    )
    core_propagation = _propagate_core_relation(
        core_relation,
        left_factor=right_interface,
        right_factor=right_interface.inverse,
        compact_ideal=compact_ideal,
        assumptions=assumptions,
        evidence=(
            regularizer,
            "bounded bilateral-ideal propagation by Z and Z_inverse",
        ),
    )
    step_two = RegularizerReductionStep(
        "left-2-core-regularizer",
        linked,
        "apply H R ~= I inside bounded Z and Z_inverse",
        regularized,
        core_propagation,
        (core_relation, compact_ideal),
    )
    cancellation = ExactIdentity(
        regularized,
        IDENTITY_PRODUCT,
        evidence=(
            right_interface.operator_times_inverse,
            "remove the structural identity and cancel Z Z_inverse",
        ),
        scope=ExactIdentityScope.OPERATORIAL,
    )
    step_three = RegularizerReductionStep(
        "left-3-exact-cancellation",
        regularized,
        "cancel Z Z_inverse exactly",
        IDENTITY_PRODUCT,
        cancellation,
        (right_interface.operator_times_inverse,),
    )
    conditional = (
        link.status is DiagonalLinkStatus.CONDITIONAL
        or regularizer.claim_status is CoreRegularizerClaimStatus.CONDITIONAL
    )
    exact = (
        link.status is DiagonalLinkStatus.EXACT
        and isinstance(core_relation, ExactIdentity)
    )
    conclusion = _final_relation(
        original,
        exact=exact,
        conditional=conditional,
        compact_ideal=compact_ideal,
        assumptions=assumptions,
        evidence=(step_one, step_two, step_three),
    )
    return (step_one, step_two, step_three), conclusion


def _derive_right_product(
    diagonal: DiagonalOperator,
    candidate: FactoredDiagonalRegularizer,
    left_interface: InvertibleAuxiliaryOperator,
    core: MellinCoreOperator,
    regularizer: MellinCoreRegularizerCertificate,
    link: DiagonalLinkIdentity,
    assumptions: AssumptionContext,
    compact_ideal: CompactIdeal | None,
) -> tuple[
    tuple[RegularizerReductionStep, ...],
    ExactIdentity | FormalIdentity | ModCompactEquivalence | None,
]:
    core_relation = regularizer.right_relation
    if core_relation is None:
        return (), None
    original = Product((*candidate.ast_product.factors, diagonal.atom))
    linked = Product(
        (
            left_interface.atom,
            regularizer.atom,
            core.atom,
            left_interface.inverse_atom,
        )
    )
    link_relation = _intermediate_relation(
        original,
        linked,
        link=link,
        evidence=(
            link.relation,
            left_interface.operator_times_inverse,
            "right-multiply the link by T_inverse and associate; no commutation",
        ),
    )
    step_one = RegularizerReductionStep(
        "right-1-link",
        original,
        "derive Z_inverse A = H T_inverse from the oriented link",
        linked,
        link_relation,
        (link, left_interface.operator_times_inverse),
    )
    regularized = Product(
        (
            left_interface.atom,
            IDENTITY_OPERATOR,
            left_interface.inverse_atom,
        )
    )
    core_propagation = _propagate_core_relation(
        core_relation,
        left_factor=left_interface,
        right_factor=left_interface.inverse,
        compact_ideal=compact_ideal,
        assumptions=assumptions,
        evidence=(
            regularizer,
            "bounded bilateral-ideal propagation by T and T_inverse",
        ),
    )
    step_two = RegularizerReductionStep(
        "right-2-core-regularizer",
        linked,
        "apply R H ~= I inside bounded T and T_inverse",
        regularized,
        core_propagation,
        (core_relation, compact_ideal),
    )
    cancellation = ExactIdentity(
        regularized,
        IDENTITY_PRODUCT,
        evidence=(
            left_interface.operator_times_inverse,
            "remove the structural identity and cancel T T_inverse",
        ),
        scope=ExactIdentityScope.OPERATORIAL,
    )
    step_three = RegularizerReductionStep(
        "right-3-exact-cancellation",
        regularized,
        "cancel T T_inverse exactly",
        IDENTITY_PRODUCT,
        cancellation,
        (left_interface.operator_times_inverse,),
    )
    conditional = (
        link.status is DiagonalLinkStatus.CONDITIONAL
        or regularizer.claim_status is CoreRegularizerClaimStatus.CONDITIONAL
    )
    exact = (
        link.status is DiagonalLinkStatus.EXACT
        and isinstance(core_relation, ExactIdentity)
    )
    conclusion = _final_relation(
        original,
        exact=exact,
        conditional=conditional,
        compact_ideal=compact_ideal,
        assumptions=assumptions,
        evidence=(step_one, step_two, step_three),
    )
    return (step_one, step_two, step_three), conclusion


def _blocked_certificate(
    *,
    diagonal_operator: DiagonalOperator,
    candidate: FactoredDiagonalRegularizer,
    status: FactoredRegularizerStatus,
    assumptions: AssumptionContext,
    evidence: tuple[object, ...],
    blocking_factors: tuple[object, ...],
    obligations: tuple[AnalyticProofObligation, ...],
    compact_ideal: CompactIdeal | None,
) -> FactoredTwoSidedRegularizerCertificate:
    left = Product((diagonal_operator.atom, *candidate.ast_product.factors))
    right = Product((*candidate.ast_product.factors, diagonal_operator.atom))
    return FactoredTwoSidedRegularizerCertificate(
        operator=diagonal_operator,
        candidate_regularizer=candidate,
        left_product=left,
        right_product=right,
        left_reduction_trace=(),
        right_reduction_trace=(),
        left_relation=None,
        right_relation=None,
        compact_ideal=compact_ideal,
        status=status,
        left_status=ProductCertificationStatus.BLOCKED,
        right_status=ProductCertificationStatus.BLOCKED,
        assumptions=assumptions,
        evidence=evidence,
        proof_obligations=obligations,
        blocking_factors=blocking_factors,
        latex=(
            r"B_1=\mathcal T_{1,-}\mathcal R_1Z_1^{-1}"
            rf"\quad\text{{{status.value}}}"
        ),
    )


def certify_factored_two_sided_diagonal_regularizer(
    *,
    diagonal_operator: DiagonalOperator,
    left_interface: InvertibleAuxiliaryOperator,
    mellin_operator: MellinCoreOperator,
    mellin_regularizer: MellinCoreRegularizerCertificate,
    right_interface: InvertibleMultiplicationOperator,
    right_interface_inverse: MultiplicationOperator,
    link_identity: DiagonalLinkIdentity,
    assumptions: AssumptionContext,
) -> FactoredTwoSidedRegularizerCertificate:
    """Certify the two complete ordered products independently."""

    if not isinstance(diagonal_operator, DiagonalOperator):
        raise TypeError("diagonal_operator must be a DiagonalOperator.")
    if not isinstance(left_interface, InvertibleAuxiliaryOperator):
        raise TypeError("left_interface must be an invertible auxiliary.")
    if not isinstance(mellin_operator, MellinCoreOperator):
        raise TypeError("mellin_operator must be a MellinCoreOperator.")
    if not isinstance(mellin_regularizer, MellinCoreRegularizerCertificate):
        raise TypeError("mellin_regularizer must be a core certificate.")
    if not isinstance(right_interface, InvertibleMultiplicationOperator):
        raise TypeError("right_interface must be an invertible multiplier.")
    if not isinstance(right_interface_inverse, MultiplicationOperator):
        raise TypeError("right_interface_inverse must be a multiplier.")
    if not isinstance(link_identity, DiagonalLinkIdentity):
        raise TypeError("link_identity must be a DiagonalLinkIdentity.")
    if not isinstance(assumptions, AssumptionContext):
        raise TypeError("assumptions must be an AssumptionContext.")
    if assumptions.consistency_status is ConsistencyStatus.INCONSISTENT:
        raise DomainOrOrderError("certification assumptions are inconsistent.")
    candidate = FactoredDiagonalRegularizer(
        atom=OperatorAtom(
            f"{diagonal_operator.atom.name}_factored_regularizer",
            kind="factored_regularizer",
        ),
        left_interface=left_interface,
        mellin_regularizer=mellin_regularizer,
        right_interface_inverse=right_interface_inverse,
    )
    obligations = _open_schur_obligations()
    typed_roles_match = (
        link_identity.diagonal_operator == diagonal_operator
        and link_identity.left_interface == left_interface
        and link_identity.mellin_operator == mellin_operator
        and link_identity.right_interface == right_interface
        and link_identity.right_interface_inverse == right_interface_inverse
        and mellin_regularizer.mellin_operator == mellin_operator
        and right_interface.inverse == right_interface_inverse
    )
    if not typed_roles_match:
        return _blocked_certificate(
            diagonal_operator=diagonal_operator,
            candidate=candidate,
            status=FactoredRegularizerStatus.BLOCKED_BY_DOMAIN_OR_ORDER,
            assumptions=assumptions,
            evidence=(link_identity, mellin_regularizer),
            blocking_factors=(
                diagonal_operator,
                left_interface,
                mellin_operator,
                mellin_regularizer,
                right_interface,
                right_interface_inverse,
            ),
            obligations=(
                AnalyticProofObligation(
                    "typed_link_roles",
                    "Supply one link referencing exactly the API's typed operators.",
                ),
                *obligations,
            ),
            compact_ideal=mellin_regularizer.compact_ideal,
        )
    space = _space(diagonal_operator, "diagonal operator")
    try:
        _require_same_space(
            space,
            (left_interface, "left interface"),
            (mellin_operator, "Mellin core"),
            (mellin_regularizer, "Mellin regularizer"),
            (right_interface, "right interface"),
            (right_interface_inverse, "right inverse"),
        )
    except DomainOrOrderError:
        return _blocked_certificate(
            diagonal_operator=diagonal_operator,
            candidate=candidate,
            status=FactoredRegularizerStatus.BLOCKED_BY_DOMAIN_OR_ORDER,
            assumptions=assumptions,
            evidence=(link_identity, mellin_regularizer),
            blocking_factors=(diagonal_operator, candidate),
            obligations=(
                AnalyticProofObligation(
                    "common_weighted_space",
                    "Make every factor preserve the same weighted Lp space.",
                ),
                *obligations,
            ),
            compact_ideal=mellin_regularizer.compact_ideal,
        )
    core_ideal = mellin_regularizer.compact_ideal
    link_ideal = link_identity.compact_ideal
    if core_ideal is not None and link_ideal is not None:
        core_ideal.require_same(link_ideal)
    compact_ideal = core_ideal or link_ideal
    if compact_ideal is not None:
        if not compact_ideal.bilateral or compact_ideal.space != space:
            raise CompactIdealMismatch(
                "certification requires one bilateral compact ideal on its space."
            )
    if link_identity.status is DiagonalLinkStatus.UNVERIFIED:
        return _blocked_certificate(
            diagonal_operator=diagonal_operator,
            candidate=candidate,
            status=FactoredRegularizerStatus.BLOCKED_BY_LINK_IDENTITY,
            assumptions=assumptions,
            evidence=(link_identity,),
            blocking_factors=(link_identity,),
            obligations=(
                AnalyticProofObligation(
                    "diagonal_link_identity",
                    "Prove the oriented Z_inverse A T = H link.",
                ),
                *obligations,
            ),
            compact_ideal=compact_ideal,
        )
    if mellin_regularizer.claim_status is CoreRegularizerClaimStatus.UNVERIFIED:
        return _blocked_certificate(
            diagonal_operator=diagonal_operator,
            candidate=candidate,
            status=FactoredRegularizerStatus.BLOCKED_BY_MELLIN_CORE_REGULARIZER,
            assumptions=assumptions,
            evidence=(mellin_regularizer,),
            blocking_factors=(mellin_regularizer,),
            obligations=(
                AnalyticProofObligation(
                    "mellin_core_regularizer",
                    "Prove at least one explicitly oriented core regularizer relation.",
                ),
                *obligations,
            ),
            compact_ideal=compact_ideal,
        )
    left_trace, left_relation = _derive_left_product(
        diagonal_operator,
        candidate,
        left_interface,
        mellin_operator,
        mellin_regularizer,
        right_interface,
        link_identity,
        assumptions,
        compact_ideal,
    )
    right_trace, right_relation = _derive_right_product(
        diagonal_operator,
        candidate,
        left_interface,
        mellin_operator,
        mellin_regularizer,
        link_identity,
        assumptions,
        compact_ideal,
    )
    conditional = (
        link_identity.status is DiagonalLinkStatus.CONDITIONAL
        or mellin_regularizer.claim_status
        is CoreRegularizerClaimStatus.CONDITIONAL
    )
    if left_relation is not None and right_relation is not None:
        status = (
            FactoredRegularizerStatus.CONDITIONAL_TWO_SIDED_FACTORED_REGULARIZER
            if conditional
            else FactoredRegularizerStatus.CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER
        )
    elif conditional:
        status = FactoredRegularizerStatus.BLOCKED_BY_MELLIN_CORE_REGULARIZER
    elif left_relation is not None:
        status = FactoredRegularizerStatus.CERTIFIED_LEFT_REGULARIZER_ONLY
    elif right_relation is not None:
        status = FactoredRegularizerStatus.CERTIFIED_RIGHT_REGULARIZER_ONLY
    else:
        status = FactoredRegularizerStatus.BLOCKED_BY_MELLIN_CORE_REGULARIZER
    left_status = (
        ProductCertificationStatus.CONDITIONAL
        if conditional and left_relation is not None
        else ProductCertificationStatus.LEFT_REGULARIZER_CERTIFIED
        if left_relation is not None
        else ProductCertificationStatus.BLOCKED
    )
    right_status = (
        ProductCertificationStatus.CONDITIONAL
        if conditional and right_relation is not None
        else ProductCertificationStatus.RIGHT_REGULARIZER_CERTIFIED
        if right_relation is not None
        else ProductCertificationStatus.BLOCKED
    )
    extra_obligations: tuple[AnalyticProofObligation, ...] = ()
    if left_relation is None:
        extra_obligations += (
            AnalyticProofObligation(
                "operator_then_regularizer_relation",
                "Prove H R ~= I for the product with A on the left.",
            ),
        )
    if right_relation is None:
        extra_obligations += (
            AnalyticProofObligation(
                "regularizer_then_operator_relation",
                "Prove R H ~= I for the product with A on the right.",
            ),
        )
    left_product = Product((diagonal_operator.atom, *candidate.ast_product.factors))
    right_product = Product((*candidate.ast_product.factors, diagonal_operator.atom))
    return FactoredTwoSidedRegularizerCertificate(
        operator=diagonal_operator,
        candidate_regularizer=candidate,
        left_product=left_product,
        right_product=right_product,
        left_reduction_trace=left_trace,
        right_reduction_trace=right_trace,
        left_relation=left_relation,
        right_relation=right_relation,
        compact_ideal=compact_ideal,
        status=status,
        left_status=left_status,
        right_status=right_status,
        assumptions=assumptions,
        evidence=(
            link_identity,
            mellin_regularizer,
            left_trace,
            right_trace,
        ),
        proof_obligations=(*extra_obligations, *obligations),
        blocking_factors=(
            ()
            if left_relation is not None and right_relation is not None
            else (mellin_regularizer,)
        ),
        latex=(
            r"B_1=\mathcal T_{1,-}\mathcal R_1Z_1^{-1},\quad "
            r"\mathcal A_{1,1}B_1\simeq I,\quad "
            r"B_1\mathcal A_{1,1}\simeq I"
            rf"\quad\text{{{status.value}}}"
        ),
    )


@dataclass(frozen=True)
class SchurRegularizerInsertion:
    """A structural Schur sandwich with an opaque certified center."""

    left_exterior: object
    central_regularizer: FactoredTwoSidedRegularizerCertificate
    right_exterior: object
    expression: Product
    schur_insertion: SchurInsertionStatus
    common_algebra_membership: UnprovedClaimStatus
    single_mellin_pdo_representation: UnprovedClaimStatus
    schur_symbol: SchurSymbolStatus
    expanded_expression: None = None

    def __post_init__(self) -> None:
        if self.central_regularizer.status is not (
            FactoredRegularizerStatus.CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER
        ):
            raise DiagonalRegularizerError(
                "Schur insertion requires a certified two-sided regularizer."
            )
        space = self.central_regularizer.operator.domain
        _require_same_space(
            space,
            (self.left_exterior, "left Schur exterior"),
            (self.right_exterior, "right Schur exterior"),
        )
        _require_bounded(self.left_exterior, "left Schur exterior")
        _require_bounded(self.right_exterior, "right Schur exterior")
        expected = Product(
            (
                _atom(self.left_exterior, "left Schur exterior"),
                self.central_regularizer.candidate_regularizer.atom,
                _atom(self.right_exterior, "right Schur exterior"),
            )
        )
        if self.expression != expected:
            raise DomainOrOrderError(
                "Schur insertion must retain the regularizer as one central atom."
            )
        if self.expanded_expression is not None:
            raise DomainOrOrderError("automatic Schur expansion is forbidden.")


def insert_factored_regularizer_in_schur(
    *,
    left_exterior: object,
    central_regularizer: FactoredTwoSidedRegularizerCertificate,
    right_exterior: object,
) -> SchurRegularizerInsertion:
    """Insert a certified center without expansion, symbol, or membership."""

    if not isinstance(
        central_regularizer, FactoredTwoSidedRegularizerCertificate
    ):
        raise TypeError("central_regularizer must be a factored certificate.")
    expression = Product(
        (
            _atom(left_exterior, "left Schur exterior"),
            central_regularizer.candidate_regularizer.atom,
            _atom(right_exterior, "right Schur exterior"),
        )
    )
    return SchurRegularizerInsertion(
        left_exterior=left_exterior,
        central_regularizer=central_regularizer,
        right_exterior=right_exterior,
        expression=expression,
        schur_insertion=SchurInsertionStatus.STRUCTURALLY_VALID,
        common_algebra_membership=UnprovedClaimStatus.UNPROVED,
        single_mellin_pdo_representation=UnprovedClaimStatus.UNPROVED,
        schur_symbol=SchurSymbolStatus.NOT_CERTIFIED,
    )


__all__ = [
    "CompactIdeal",
    "CompactIdealMismatch",
    "CoreRegularizerClaimStatus",
    "CoreRegularizerProperty",
    "DiagonalLinkIdentity",
    "DiagonalLinkStatus",
    "DiagonalOperator",
    "DiagonalRegularizerError",
    "DomainOrOrderError",
    "FactoredDiagonalRegularizer",
    "FactoredRegularizerStatus",
    "FactoredTwoSidedRegularizerCertificate",
    "IDENTITY_OPERATOR",
    "IDENTITY_PRODUCT",
    "InvertibleAuxiliaryOperator",
    "InvertibleMultiplicationOperator",
    "MellinCoreOperator",
    "MellinCoreRegularizerCertificate",
    "MissingRegularizerEvidence",
    "ProductCertificationStatus",
    "RegularizerReductionStep",
    "SchurInsertionStatus",
    "SchurRegularizerInsertion",
    "SchurSymbolStatus",
    "UnprovedClaimStatus",
    "certify_factored_two_sided_diagonal_regularizer",
    "insert_factored_regularizer_in_schur",
    "make_mellin_core_regularizer_certificate",
    "propagate_compact_ideal",
]
