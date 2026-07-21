"""Fail-closed analysis of the complete first-Schur operator word.

The authoritative expression is always the five-factor AST supplied to
``factor_full_first_schur_word``.  Expanded models are evidence attached to
that expression; they never replace it implicitly.  In particular, this
module does not commute localizers, Cauchy factors, multipliers, or
Wiener--Hopf operators in order to manufacture a dilation pattern.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import sympy as sp

from .domains import AssumptionContext, ConsistencyStatus, IncompatibleDomainError
from .mellin.operators import (
    DilationCovarianceTrace,
    MellinConvolutionOperator,
    MellinOperatorError,
    MellinPseudodifferentialOperator,
    MultiplicationOperator,
    RegularizerMellinRepresentation,
    RegularizerRepresentationStatus,
    WeightedDilationOperator,
    WeightedLpSpace,
    WienerHopfOperator,
    conjugate_mellin_pdo_by_dilation,
    conjugate_multiplication_by_dilation,
)
from .operators import LinearCombination, OperatorAtom, Product, Term
from .relations import ExactBlock
from .schur_dilation import (
    SchurBlockRepresentation,
    SchurDilationCancellationResult,
    SchurDilationError,
    cancel_relative_dilations_in_schur_correction,
)
from .semantics import (
    AnalyticProofObligation,
    DerivationRelationKind,
    ExactIdentity,
    ExactIdentityScope,
    FormalIdentity,
    ModCompactEquivalence,
)


class FullSchurWordError(ValueError):
    """Raised when full-word metadata are malformed or inconsistent."""


class CauchyPolarity(str, Enum):
    PLUS = "PLUS"
    MINUS = "MINUS"


class CutoffEquivalenceStatus(str, Enum):
    ASSUMED = "ASSUMED"
    CERTIFIED = "CERTIFIED"
    UNPROVED = "UNPROVED"


class CutoffEquivalenceSide(str, Enum):
    LEFT = "LEFT"
    RIGHT = "RIGHT"
    BOTH = "BOTH"


class RegularizerInterfaceRole(str, Enum):
    LEFT_AUXILIARY = "LEFT_AUXILIARY"
    RIGHT_AUXILIARY = "RIGHT_AUXILIARY"


class CoreFactorizationStatus(str, Enum):
    EXACT = "EXACT"
    MODULO_COMPACTS = "MODULO_COMPACTS"
    CONDITIONAL_ON_CUTOFF = "CONDITIONAL_ON_CUTOFF"
    BLOCKED = "BLOCKED"


class FullSchurFactorizationStatus(str, Enum):
    FULL_WORD_EXACT_FACTORIZATION = "FULL_WORD_EXACT_FACTORIZATION"
    FULL_WORD_MOD_COMPACTS_FACTORIZATION = (
        "FULL_WORD_MOD_COMPACTS_FACTORIZATION"
    )
    FACTORIZATION_AFTER_CUTOFF_EQUIVALENCE = (
        "FACTORIZATION_AFTER_CUTOFF_EQUIVALENCE"
    )
    BLOCKED_BY_CUTOFF_CONTROL = "BLOCKED_BY_CUTOFF_CONTROL"
    BLOCKED_BY_OPERATOR_ORDER = "BLOCKED_BY_OPERATOR_ORDER"
    BLOCKED_BY_REGULARIZER_INTERFACE = "BLOCKED_BY_REGULARIZER_INTERFACE"


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


def _space_pair(value: object, field_name: str) -> tuple[WeightedLpSpace, WeightedLpSpace]:
    domain = getattr(value, "domain", None)
    codomain = getattr(value, "codomain", None)
    if not isinstance(domain, WeightedLpSpace) or not isinstance(
        codomain, WeightedLpSpace
    ):
        raise TypeError(
            f"{field_name} must expose typed WeightedLpSpace domain/codomain."
        )
    return domain, codomain


def _same_space(*values: object) -> WeightedLpSpace | None:
    expected: WeightedLpSpace | None = None
    for value in values:
        pair = _space_pair(value, "operator")
        if pair[0] != pair[1]:
            return None
        if expected is None:
            expected = pair[0]
        elif pair[0] != expected:
            return None
    return expected


@dataclass(frozen=True)
class SpatialCutoffOperator:
    """A multiplier cutoff whose geometric and radial data are explicit."""

    atom: OperatorAtom
    cutoff_symbol: sp.Expr
    radial_variable: sp.Symbol
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    support_region: str
    equals_one_near: str
    equals_zero_away_from: str
    coordinate_system: str
    radial_scale: sp.Expr
    source: str
    evidence: tuple[object, ...] = ()
    bounded: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.cutoff_symbol, sp.Expr):
            raise TypeError("cutoff_symbol must be a SymPy expression.")
        if not isinstance(self.radial_variable, sp.Symbol):
            raise TypeError("radial_variable must be a SymPy Symbol.")
        _space_pair(self, "cutoff")
        if self.domain != self.codomain:
            raise FullSchurWordError("a cutoff must preserve one weighted space.")
        scale = sp.sympify(self.radial_scale)
        if isinstance(self.radial_scale, bool) or scale.is_nonpositive is True:
            raise FullSchurWordError("radial_scale must not be nonpositive.")
        for name in (
            "support_region",
            "equals_one_near",
            "equals_zero_away_from",
            "coordinate_system",
            "source",
        ):
            object.__setattr__(self, name, _nonempty(getattr(self, name), name))
        if not isinstance(self.bounded, bool):
            raise TypeError("bounded must be a bool.")
        object.__setattr__(self, "radial_scale", scale)
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))

    @property
    def multiplication(self) -> MultiplicationOperator:
        return MultiplicationOperator(
            atom=self.atom,
            symbol=self.cutoff_symbol,
            radial_variable=self.radial_variable,
            domain=self.domain,
            codomain=self.codomain,
            source=self.source,
            evidence=self.evidence,
            bounded=self.bounded,
        )


@dataclass(frozen=True)
class CutoffDilationCovarianceTrace:
    source_cutoff: SpatialCutoffOperator
    scaled_cutoff: SpatialCutoffOperator
    multiplication_trace: DilationCovarianceTrace
    exact_identity: ExactIdentity = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.source_cutoff, SpatialCutoffOperator) or not isinstance(
            self.scaled_cutoff, SpatialCutoffOperator
        ):
            raise TypeError("cutoff covariance requires typed cutoff operators.")
        if not isinstance(self.multiplication_trace, DilationCovarianceTrace):
            raise TypeError("multiplication_trace must be a DilationCovarianceTrace.")
        identity = ExactIdentity(
            self.multiplication_trace.source_product,
            Product((self.scaled_cutoff.atom,)),
            evidence=(
                self.multiplication_trace.exact_identity,
                "the cutoff geometry was rescaled, not erased",
            ),
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=self.multiplication_trace.exact_identity.hypotheses,
        )
        object.__setattr__(self, "exact_identity", identity)


def conjugate_cutoff_by_dilation(
    dilation: WeightedDilationOperator,
    cutoff: SpatialCutoffOperator,
    inverse: WeightedDilationOperator,
) -> CutoffDilationCovarianceTrace:
    """Return the exact rule ``V_gamma M_chi V_gamma^-1=M_chi(gamma dot)``."""

    if not isinstance(cutoff, SpatialCutoffOperator):
        raise TypeError("cutoff must be a SpatialCutoffOperator.")
    trace = conjugate_multiplication_by_dilation(
        dilation,
        cutoff.multiplication,
        inverse,
    )
    transformed = trace.transformed_operator
    if not isinstance(transformed, MultiplicationOperator):
        raise FullSchurWordError("multiplication covariance returned a wrong type.")
    scaled = SpatialCutoffOperator(
        atom=transformed.atom,
        cutoff_symbol=transformed.symbol,
        radial_variable=cutoff.radial_variable,
        domain=cutoff.domain,
        codomain=cutoff.codomain,
        support_region=f"{cutoff.support_region}; radial coordinate divided by gamma",
        equals_one_near=f"scaled({cutoff.equals_one_near}, 1/gamma)",
        equals_zero_away_from=f"scaled({cutoff.equals_zero_away_from}, 1/gamma)",
        coordinate_system=cutoff.coordinate_system,
        radial_scale=sp.simplify(cutoff.radial_scale * dilation.gamma),
        source=f"{cutoff.source}; exact pointwise dilation covariance",
        evidence=(*cutoff.evidence, trace.exact_identity),
        bounded=cutoff.bounded,
    )
    return CutoffDilationCovarianceTrace(cutoff, scaled, trace)


@dataclass(frozen=True)
class LocalizedOperator:
    """An ordered two-cutoff localization ``M_left K M_right``."""

    atom: OperatorAtom
    left_cutoff: SpatialCutoffOperator
    core: object
    right_cutoff: SpatialCutoffOperator
    source: str
    evidence: tuple[object, ...]
    bounded: bool = True
    ast_product: Product = field(init=False)
    domain: WeightedLpSpace = field(init=False)
    codomain: WeightedLpSpace = field(init=False)
    exact_definition: ExactIdentity = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.left_cutoff, SpatialCutoffOperator) or not isinstance(
            self.right_cutoff, SpatialCutoffOperator
        ):
            raise TypeError("localized operators require two typed cutoffs.")
        core_atom = _atom(self.core, "core")
        space = _same_space(self.left_cutoff, self.core, self.right_cutoff)
        if space is None:
            raise FullSchurWordError(
                "localized operator factors must preserve the same weighted space."
            )
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise FullSchurWordError("a localized operator requires evidence.")
        if not isinstance(self.bounded, bool):
            raise TypeError("bounded must be a bool.")
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "domain", space)
        object.__setattr__(self, "codomain", space)
        object.__setattr__(
            self,
            "ast_product",
            Product((self.left_cutoff.atom, core_atom, self.right_cutoff.atom)),
        )
        object.__setattr__(
            self,
            "exact_definition",
            ExactIdentity(
                self.atom,
                Product((self.left_cutoff.atom, core_atom, self.right_cutoff.atom)),
                evidence=(self.source, evidence),
                scope=ExactIdentityScope.OPERATORIAL,
            ),
        )


@dataclass(frozen=True)
class CauchyProjectionOperator:
    atom: OperatorAtom
    polarity: CauchyPolarity
    mellin_operator: MellinConvolutionOperator
    source: str
    evidence: tuple[object, ...]
    bounded: bool = True
    domain: WeightedLpSpace = field(init=False)
    codomain: WeightedLpSpace = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.polarity, CauchyPolarity):
            raise TypeError("polarity must be a CauchyPolarity.")
        if not isinstance(self.mellin_operator, MellinConvolutionOperator):
            raise TypeError(
                "mellin_operator must be a MellinConvolutionOperator."
            )
        domain, codomain = _space_pair(self.mellin_operator, "mellin_operator")
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise FullSchurWordError("a Cauchy factor requires representation evidence.")
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "domain", domain)
        object.__setattr__(self, "codomain", codomain)


@dataclass(frozen=True)
class CoefficientMultiplicationOperator:
    multiplication: MultiplicationOperator
    coefficient_role: str

    def __post_init__(self) -> None:
        if not isinstance(self.multiplication, MultiplicationOperator):
            raise TypeError("multiplication must be a MultiplicationOperator.")
        object.__setattr__(
            self,
            "coefficient_role",
            _nonempty(self.coefficient_role, "coefficient_role"),
        )

    @property
    def atom(self) -> OperatorAtom:
        return self.multiplication.atom

    @property
    def domain(self) -> WeightedLpSpace:
        return self.multiplication.domain

    @property
    def codomain(self) -> WeightedLpSpace:
        return self.multiplication.codomain

    @property
    def bounded(self) -> bool:
        return self.multiplication.bounded


@dataclass(frozen=True)
class TransportedShiftOperator:
    """Typed exact product of a transport multiplier and a dilation."""

    atom: OperatorAtom
    coefficient: CoefficientMultiplicationOperator
    dilation: WeightedDilationOperator
    source: str
    evidence: tuple[object, ...]
    exact_factorization: ExactIdentity = field(init=False)
    domain: WeightedLpSpace = field(init=False)
    codomain: WeightedLpSpace = field(init=False)
    bounded: bool = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.coefficient, CoefficientMultiplicationOperator):
            raise TypeError("coefficient must be typed coefficient multiplication.")
        if not isinstance(self.dilation, WeightedDilationOperator):
            raise TypeError("dilation must be a WeightedDilationOperator.")
        space = _same_space(self.coefficient, self.dilation)
        if space is None:
            raise FullSchurWordError("transported-shift factors use different spaces.")
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise FullSchurWordError("transported shift requires exact evidence.")
        source = _nonempty(self.source, "source")
        identity = ExactIdentity(
            self.atom,
            Product((self.coefficient.atom, self.dilation.atom)),
            evidence=(source, evidence),
            scope=ExactIdentityScope.OPERATORIAL,
        )
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "exact_factorization", identity)
        object.__setattr__(self, "domain", space)
        object.__setattr__(self, "codomain", space)
        object.__setattr__(
            self,
            "bounded",
            self.coefficient.bounded and self.dilation.bounded,
        )


@dataclass(frozen=True)
class AuxiliaryOperator:
    """A typed auxiliary such as ``T_-=U^-1 P^+ +P^-``."""

    atom: OperatorAtom
    expression: LinearCombination
    components: tuple[object, ...]
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    source: str
    evidence: tuple[object, ...]
    bounded: bool = True
    exact_definition: ExactIdentity = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.expression, LinearCombination):
            raise TypeError("expression must be a LinearCombination.")
        components = _items(self.components, "components")
        if not components:
            raise FullSchurWordError("an auxiliary needs explicit components.")
        for component in components:
            _atom(component, "auxiliary component")
        component_space = _same_space(*components)
        if component_space is None or component_space != self.domain:
            raise FullSchurWordError(
                "auxiliary components must preserve its declared weighted space."
            )
        if self.domain != self.codomain:
            raise FullSchurWordError("an auxiliary must preserve one weighted space.")
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise FullSchurWordError("an auxiliary definition requires evidence.")
        source = _nonempty(self.source, "source")
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(
            self,
            "exact_definition",
            ExactIdentity(
                self.atom,
                self.expression,
                evidence=(source, evidence),
                scope=ExactIdentityScope.OPERATORIAL,
            ),
        )


@dataclass(frozen=True)
class CutoffEquivalenceClaim:
    original_cutoff: SpatialCutoffOperator
    replacement_cutoff: SpatialCutoffOperator
    localized_operator: LocalizedOperator
    status: CutoffEquivalenceStatus
    ideal: object
    side: CutoffEquivalenceSide
    assumptions: tuple[object, ...]
    source: str
    evidence: tuple[object, ...] = ()
    semantic_relation: FormalIdentity | ModCompactEquivalence | None = field(
        init=False
    )

    def __post_init__(self) -> None:
        if not isinstance(self.original_cutoff, SpatialCutoffOperator) or not isinstance(
            self.replacement_cutoff, SpatialCutoffOperator
        ):
            raise TypeError("cutoff claims require two typed cutoffs.")
        if not isinstance(self.localized_operator, LocalizedOperator):
            raise TypeError("localized_operator must be a LocalizedOperator.")
        if not isinstance(self.status, CutoffEquivalenceStatus):
            raise TypeError("status must be a CutoffEquivalenceStatus.")
        if not isinstance(self.side, CutoffEquivalenceSide):
            raise TypeError("side must be a CutoffEquivalenceSide.")
        if self.ideal is None:
            raise FullSchurWordError("a cutoff claim must name its ideal.")
        assumptions = _items(self.assumptions, "assumptions")
        evidence = _items(self.evidence, "evidence")
        source = _nonempty(self.source, "source")
        if self.status is CutoffEquivalenceStatus.CERTIFIED and not evidence:
            raise FullSchurWordError("a certified cutoff claim requires evidence.")
        if self.status is CutoffEquivalenceStatus.ASSUMED and not assumptions:
            raise FullSchurWordError("an assumed cutoff claim requires assumptions.")
        if self.original_cutoff.domain != self.replacement_cutoff.domain:
            raise FullSchurWordError("cutoff claim spaces do not agree.")
        if self.side is CutoffEquivalenceSide.LEFT and (
            self.original_cutoff != self.localized_operator.left_cutoff
        ):
            raise FullSchurWordError(
                "a left cutoff claim must reference the localized left cutoff."
            )
        if self.side is CutoffEquivalenceSide.RIGHT and (
            self.original_cutoff != self.localized_operator.right_cutoff
        ):
            raise FullSchurWordError(
                "a right cutoff claim must reference the localized right cutoff."
            )
        if self.side is CutoffEquivalenceSide.BOTH and not (
            self.original_cutoff == self.localized_operator.left_cutoff
            and self.original_cutoff == self.localized_operator.right_cutoff
        ):
            raise FullSchurWordError(
                "a BOTH claim requires the same original cutoff on both sides."
            )
        core_atom = _atom(self.localized_operator.core, "localized core")
        left = Product((self.original_cutoff.atom, core_atom))
        right = Product((self.replacement_cutoff.atom, core_atom))
        if self.side is CutoffEquivalenceSide.RIGHT:
            left = Product((core_atom, self.original_cutoff.atom))
            right = Product((core_atom, self.replacement_cutoff.atom))
        elif self.side is CutoffEquivalenceSide.BOTH:
            left = Product(
                (
                    self.original_cutoff.atom,
                    core_atom,
                    self.original_cutoff.atom,
                )
            )
            right = Product(
                (
                    self.replacement_cutoff.atom,
                    core_atom,
                    self.replacement_cutoff.atom,
                )
            )
        relation: FormalIdentity | ModCompactEquivalence | None
        if self.status is CutoffEquivalenceStatus.CERTIFIED:
            relation = ModCompactEquivalence(
                left,
                right,
                space=self.original_cutoff.domain,
                compact_ideal=self.ideal,
                evidence=(source, evidence, assumptions),
            )
        elif self.status is CutoffEquivalenceStatus.ASSUMED:
            relation = FormalIdentity(
                left,
                right,
                justification=(source, assumptions, evidence, self.ideal),
            )
        else:
            relation = None
        object.__setattr__(self, "assumptions", assumptions)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "semantic_relation", relation)


@dataclass(frozen=True)
class FullWordBlock:
    """An existing Schur block plus full-word-only interface provenance."""

    representation: SchurBlockRepresentation
    embedded_interfaces: tuple[RegularizerInterfaceRole, ...] = ()
    cutoff_claims: tuple[CutoffEquivalenceClaim, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.representation, SchurBlockRepresentation):
            raise TypeError("representation must be a SchurBlockRepresentation.")
        interfaces = _items(self.embedded_interfaces, "embedded_interfaces")
        claims = _items(self.cutoff_claims, "cutoff_claims")
        if not all(isinstance(item, RegularizerInterfaceRole) for item in interfaces):
            raise TypeError("embedded_interfaces contains an invalid role.")
        if not all(isinstance(item, CutoffEquivalenceClaim) for item in claims):
            raise TypeError("cutoff_claims must contain CutoffEquivalenceClaim objects.")
        object.__setattr__(self, "embedded_interfaces", interfaces)
        object.__setattr__(self, "cutoff_claims", claims)

    @property
    def atom(self) -> OperatorAtom:
        return self.representation.atom


@dataclass(frozen=True)
class TransportedMellinCore:
    """Keep ``Phi^-1 Op_M(r) Phi`` distinct from the full regularizer."""

    representation: RegularizerMellinRepresentation
    phi_inverse: OperatorAtom
    phi: OperatorAtom
    source: str
    evidence: tuple[object, ...]
    transported_word: Product = field(init=False)
    semantic_relation: ExactIdentity | FormalIdentity | ModCompactEquivalence = field(
        init=False
    )
    domain: WeightedLpSpace = field(init=False)
    codomain: WeightedLpSpace = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.representation, RegularizerMellinRepresentation):
            raise TypeError("representation must be a RegularizerMellinRepresentation.")
        if not isinstance(self.phi_inverse, OperatorAtom) or not isinstance(
            self.phi, OperatorAtom
        ):
            raise TypeError("phi_inverse and phi must be OperatorAtom objects.")
        evidence = _items(self.evidence, "evidence")
        source = _nonempty(self.source, "source")
        if not evidence:
            raise FullSchurWordError("transported Mellin core requires evidence.")
        word = Product(
            (
                self.phi_inverse,
                self.representation.mellin_operator.atom,
                self.phi,
            )
        )
        left = self.representation.regularizer.operator
        status = self.representation.status
        relation: ExactIdentity | FormalIdentity | ModCompactEquivalence
        if status is RegularizerRepresentationStatus.CERTIFIED_EXACT:
            relation = ExactIdentity(
                left,
                word,
                evidence=(source, evidence, self.representation.semantic_relation),
                scope=ExactIdentityScope.OPERATORIAL,
                hypotheses=self.representation.hypotheses,
            )
        elif status is RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT:
            relation = ModCompactEquivalence(
                left,
                word,
                space=self.representation.mellin_operator.domain,
                compact_ideal=self.representation.compact_ideal
                or "compact operators",
                evidence=(source, evidence, self.representation.semantic_relation),
            )
        else:
            relation = FormalIdentity(
                left,
                word,
                justification=(source, evidence, self.representation.hypotheses),
            )
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "transported_word", word)
        object.__setattr__(self, "semantic_relation", relation)
        object.__setattr__(self, "domain", self.representation.mellin_operator.domain)
        object.__setattr__(self, "codomain", self.representation.mellin_operator.codomain)

    @property
    def atom(self) -> OperatorAtom:
        operator = self.representation.regularizer.operator
        if not isinstance(operator, OperatorAtom):
            raise FullSchurWordError("regularizer core lost its authoritative atom.")
        return operator


@dataclass(frozen=True)
class SchurExteriorCoreResult:
    original_word: Product
    expanded_word: Product | None
    status: CoreFactorizationStatus
    exterior: object | None
    dilation: WeightedDilationOperator | None
    relation: ExactIdentity | FormalIdentity | ModCompactEquivalence | None
    evidence: tuple[object, ...]
    proof_obligations: tuple[AnalyticProofObligation, ...]
    blocking_factors: tuple[object, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.original_word, Product):
            raise TypeError("original_word must be a Product.")
        if self.expanded_word is not None and not isinstance(
            self.expanded_word, Product
        ):
            raise TypeError("expanded_word must be a Product or None.")
        if not isinstance(self.status, CoreFactorizationStatus):
            raise TypeError("status must be a CoreFactorizationStatus.")
        if self.status is CoreFactorizationStatus.BLOCKED:
            if self.relation is not None or not self.blocking_factors:
                raise FullSchurWordError(
                    "a blocked exterior core needs factors and no relation."
                )
        elif self.exterior is None or self.dilation is None or self.relation is None:
            raise FullSchurWordError(
                "a successful exterior core needs exterior, dilation, and relation."
            )


def _base_obligations() -> tuple[AnalyticProofObligation, ...]:
    return (
        AnalyticProofObligation(
            "left_core_factorization",
            "Prove the order-preserving left exterior/dilation factorization.",
        ),
        AnalyticProofObligation(
            "right_core_factorization",
            "Prove the order-preserving inverse-dilation/right exterior factorization.",
        ),
        AnalyticProofObligation(
            "cutoff_replacement_mod_compacts",
            "Prove every side-specific scaled-cutoff replacement used by the word.",
        ),
        AnalyticProofObligation(
            "regularizer_mellin_representation",
            "Retain and verify R1=Phi^-1 Op_M(r11) Phi on the stated spaces.",
        ),
        AnalyticProofObligation(
            "wh_mellin_wh_sandwich_membership",
            "Prove membership of the ordered exterior--Mellin--exterior sandwich.",
        ),
        AnalyticProofObligation(
            "fredholm_algebra_membership",
            "Prove membership in an algebra carrying the required Fredholm calculus.",
        ),
        AnalyticProofObligation(
            "schur_correction_symbol",
            "Construct the correction symbol only after its algebra membership.",
        ),
    )


def _relation_kind_status(
    relation: DerivationRelationKind,
) -> CoreFactorizationStatus:
    if relation is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE:
        return CoreFactorizationStatus.MODULO_COMPACTS
    return CoreFactorizationStatus.EXACT


def _analyze_left_core(
    block: FullWordBlock,
    auxiliary: object,
) -> SchurExteriorCoreResult:
    representation = block.representation
    original = Product((block.atom, _atom(auxiliary, "left_auxiliary")))
    factors = (*representation.factors, auxiliary)
    expanded = Product(tuple(_atom(item, "left core factor") for item in factors))
    obligations = _base_obligations()
    if len(factors) != 2 or not isinstance(
        factors[0], (WienerHopfOperator, LocalizedOperator)
    ) or not isinstance(factors[1], WeightedDilationOperator):
        return SchurExteriorCoreResult(
            original,
            expanded,
            CoreFactorizationStatus.BLOCKED,
            None,
            None,
            None,
            representation.evidence,
            obligations,
            tuple(factors),
        )
    exterior, dilation = factors
    if _same_space(exterior, dilation) is None:
        return SchurExteriorCoreResult(
            original,
            expanded,
            CoreFactorizationStatus.BLOCKED,
            None,
            None,
            None,
            representation.evidence,
            obligations,
            (exterior, dilation),
        )
    model = Product((_atom(exterior, "left exterior"), dilation.atom))
    relation: ExactIdentity | ModCompactEquivalence
    if representation.relation_kind is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE:
        relation = ModCompactEquivalence(
            original,
            model,
            space=dilation.domain,
            compact_ideal=representation.compact_ideal,
            evidence=representation.semantic_relation,
        )
    else:
        relation = ExactIdentity(
            original,
            model,
            evidence=representation.semantic_relation,
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=representation.hypotheses,
        )
    return SchurExteriorCoreResult(
        original,
        expanded,
        _relation_kind_status(representation.relation_kind),
        exterior,
        dilation,
        relation,
        (representation.semantic_relation,),
        obligations,
        (),
    )


def _analyze_right_core(
    auxiliary: object,
    block: FullWordBlock,
) -> SchurExteriorCoreResult:
    representation = block.representation
    original = Product((_atom(auxiliary, "right_auxiliary"), block.atom))
    factors = (auxiliary, *representation.factors)
    expanded = Product(tuple(_atom(item, "right core factor") for item in factors))
    obligations = _base_obligations()
    if len(factors) != 2 or not isinstance(
        factors[0], WeightedDilationOperator
    ) or not isinstance(factors[1], (WienerHopfOperator, LocalizedOperator)):
        return SchurExteriorCoreResult(
            original,
            expanded,
            CoreFactorizationStatus.BLOCKED,
            None,
            None,
            None,
            representation.evidence,
            obligations,
            tuple(factors),
        )
    dilation, exterior = factors
    if _same_space(dilation, exterior) is None:
        return SchurExteriorCoreResult(
            original,
            expanded,
            CoreFactorizationStatus.BLOCKED,
            None,
            None,
            None,
            representation.evidence,
            obligations,
            (dilation, exterior),
        )
    model = Product((dilation.atom, _atom(exterior, "right exterior")))
    relation: ExactIdentity | ModCompactEquivalence
    if representation.relation_kind is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE:
        relation = ModCompactEquivalence(
            original,
            model,
            space=dilation.domain,
            compact_ideal=representation.compact_ideal,
            evidence=representation.semantic_relation,
        )
    else:
        relation = ExactIdentity(
            original,
            model,
            evidence=representation.semantic_relation,
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=representation.hypotheses,
        )
    return SchurExteriorCoreResult(
        original,
        expanded,
        _relation_kind_status(representation.relation_kind),
        exterior,
        dilation,
        relation,
        (representation.semantic_relation,),
        obligations,
        (),
    )


@dataclass(frozen=True)
class FullSchurWordFactorizationResult:
    original_word: Product
    left_core: SchurExteriorCoreResult
    mellin_core: object
    right_core: SchurExteriorCoreResult
    factorized_word: Product | None
    relation: ExactIdentity | FormalIdentity | ModCompactEquivalence | None
    status: FullSchurFactorizationStatus
    evidence: tuple[object, ...]
    assumptions: AssumptionContext
    proof_obligations: tuple[AnalyticProofObligation, ...]
    blocking_factors: tuple[object, ...]
    latex: str
    expanded_original_word: Product | None = None
    cancellation_result: SchurDilationCancellationResult | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.original_word, Product):
            raise TypeError("original_word must be a Product.")
        if not isinstance(self.left_core, SchurExteriorCoreResult) or not isinstance(
            self.right_core, SchurExteriorCoreResult
        ):
            raise TypeError("left_core and right_core must be typed core results.")
        if not isinstance(self.status, FullSchurFactorizationStatus):
            raise TypeError("status must be a FullSchurFactorizationStatus.")
        if not isinstance(self.assumptions, AssumptionContext):
            raise TypeError("assumptions must be an AssumptionContext.")
        if not all(
            isinstance(item, AnalyticProofObligation)
            for item in self.proof_obligations
        ):
            raise TypeError("proof_obligations must be typed.")
        blocked = self.status in {
            FullSchurFactorizationStatus.BLOCKED_BY_CUTOFF_CONTROL,
            FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER,
            FullSchurFactorizationStatus.BLOCKED_BY_REGULARIZER_INTERFACE,
        }
        if blocked:
            if self.factorized_word is not None or self.relation is not None:
                raise FullSchurWordError("a blocked result cannot claim a factorization.")
            if not self.blocking_factors:
                raise FullSchurWordError("a blocked result needs precise factors.")
        else:
            if self.factorized_word is None or self.relation is None:
                raise FullSchurWordError("a successful result requires a relation.")
            if self.relation.left != self.original_word:
                raise FullSchurWordError("the relation lost the authoritative word.")
            if self.relation.right != self.factorized_word:
                raise FullSchurWordError("the relation lost the factorized word.")

    @property
    def succeeded(self) -> bool:
        return self.status not in {
            FullSchurFactorizationStatus.BLOCKED_BY_CUTOFF_CONTROL,
            FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER,
            FullSchurFactorizationStatus.BLOCKED_BY_REGULARIZER_INTERFACE,
        }


def _render_latex(
    word: Product,
    status: FullSchurFactorizationStatus,
    obligations: tuple[AnalyticProofObligation, ...],
    factorized: Product | None = None,
) -> str:
    def render(product: Product) -> str:
        names = (factor.name.replace("_", r"\_") for factor in product.factors)
        return r"\,".join(rf"\operatorname{{{name}}}" for name in names)

    original = render(word)
    obligation_text = ", ".join(
        item.key.replace("_", r"\_") for item in obligations
    )
    rendered_status = status.value.replace("_", r"\_")
    if factorized is None:
        return (
            original
            + rf"\quad\text{{[{rendered_status}; open: {obligation_text}]}}"
        )
    symbol = "=" if status is FullSchurFactorizationStatus.FULL_WORD_EXACT_FACTORIZATION else r"\simeq"
    return (
        original
        + rf"\;{symbol}\;"
        + render(factorized)
        + rf"\quad\text{{[{rendered_status}; open: {obligation_text}]}}"
    )


def _blocked_core(word: Product, factors: tuple[object, ...]) -> SchurExteriorCoreResult:
    return SchurExteriorCoreResult(
        word,
        None,
        CoreFactorizationStatus.BLOCKED,
        None,
        None,
        None,
        (),
        _base_obligations(),
        factors,
    )


def _blocked_result(
    *,
    original: Product,
    left_core: SchurExteriorCoreResult,
    mellin_core: object,
    right_core: SchurExteriorCoreResult,
    status: FullSchurFactorizationStatus,
    assumptions: AssumptionContext,
    factors: tuple[object, ...],
    evidence: tuple[object, ...] = (),
    expanded: Product | None = None,
) -> FullSchurWordFactorizationResult:
    obligations = _base_obligations()
    return FullSchurWordFactorizationResult(
        original,
        left_core,
        mellin_core,
        right_core,
        None,
        None,
        status,
        evidence,
        assumptions,
        obligations,
        factors,
        _render_latex(original, status, obligations),
        expanded,
        None,
    )


def _claims_for_localization(
    blocks: tuple[FullWordBlock, FullWordBlock],
    localized: tuple[LocalizedOperator, ...],
) -> tuple[CutoffEquivalenceClaim, ...] | None:
    claims = tuple(claim for block in blocks for claim in block.cutoff_claims)
    if not localized:
        return ()
    if not claims or any(
        claim.status is CutoffEquivalenceStatus.UNPROVED for claim in claims
    ):
        return None
    covered = {claim.localized_operator.atom for claim in claims}
    if any(item.atom not in covered for item in localized):
        return None
    for item in localized:
        sides = {
            claim.side
            for claim in claims
            if claim.localized_operator.atom == item.atom
        }
        if CutoffEquivalenceSide.BOTH not in sides and not {
            CutoffEquivalenceSide.LEFT,
            CutoffEquivalenceSide.RIGHT,
        }.issubset(sides):
            return None
    return claims


def _mellin_core_atom(value: object) -> OperatorAtom:
    if isinstance(value, TransportedMellinCore):
        return value.atom
    if isinstance(value, RegularizerMellinRepresentation):
        operator = value.regularizer.operator
        if isinstance(operator, OperatorAtom):
            return operator
    operator = getattr(value, "operator", None)
    if isinstance(operator, OperatorAtom):
        return operator
    return _atom(value, "mellin_core")


def _common_compact_ideal(
    left: FullWordBlock,
    core: TransportedMellinCore,
    right: FullWordBlock,
    claims: tuple[CutoffEquivalenceClaim, ...],
) -> object | None:
    ideals: list[object] = []
    for block in (left, right):
        if (
            block.representation.relation_kind
            is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
        ):
            ideals.append(block.representation.compact_ideal)
    if (
        core.representation.status
        is RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT
    ):
        ideals.append(core.representation.compact_ideal or "compact operators")
    ideals.extend(
        claim.ideal
        for claim in claims
        if claim.status is CutoffEquivalenceStatus.CERTIFIED
    )
    if not ideals:
        return None
    first = ideals[0]
    if any((item == first) is not True for item in ideals[1:]):
        raise FullSchurWordError("modulo-compact inputs use different ideals.")
    return first


def factor_full_first_schur_word(
    *,
    left_block: FullWordBlock,
    left_auxiliary: object,
    mellin_core: object,
    right_auxiliary: object,
    right_block: FullWordBlock,
    assumptions: AssumptionContext,
) -> FullSchurWordFactorizationResult:
    """Factor a typed five-factor word without changing its real order.

    A successful exact synthetic pattern expands contiguously as
    ``W_left, V_gamma, Phi^-1, Op_M(r), Phi, V_gamma^-1, W_right``.
    Localized patterns additionally require explicit cutoff claims.  Raw atoms
    and incomplete regularizer interfaces are rejected rather than guessed.
    """

    if not isinstance(assumptions, AssumptionContext):
        raise TypeError("assumptions must be an AssumptionContext.")
    if assumptions.consistency_status is ConsistencyStatus.INCONSISTENT:
        raise FullSchurWordError("full-word analysis rejects inconsistent assumptions.")
    if not isinstance(left_block, FullWordBlock) or not isinstance(
        right_block, FullWordBlock
    ):
        raise TypeError("left_block and right_block must be FullWordBlock objects.")
    left_aux_atom = _atom(left_auxiliary, "left_auxiliary")
    right_aux_atom = _atom(right_auxiliary, "right_auxiliary")
    core_atom = _mellin_core_atom(mellin_core)
    original = Product(
        (
            left_block.atom,
            left_aux_atom,
            core_atom,
            right_aux_atom,
            right_block.atom,
        )
    )
    left_word = Product((left_block.atom, left_aux_atom))
    right_word = Product((right_aux_atom, right_block.atom))
    provisional_left = _blocked_core(left_word, (left_block, left_auxiliary))
    provisional_right = _blocked_core(right_word, (right_auxiliary, right_block))

    duplicate_factors: list[object] = []
    if (
        RegularizerInterfaceRole.LEFT_AUXILIARY
        in left_block.embedded_interfaces
    ):
        duplicate_factors.extend((left_block, left_auxiliary))
    if (
        RegularizerInterfaceRole.RIGHT_AUXILIARY
        in right_block.embedded_interfaces
    ):
        duplicate_factors.extend((right_auxiliary, right_block))
    if duplicate_factors:
        expanded = Product(
            tuple(
                _atom(item, "expanded article factor")
                for item in (
                    *left_block.representation.factors,
                    left_auxiliary,
                    core_atom,
                    right_auxiliary,
                    *right_block.representation.factors,
                )
            )
        )
        return _blocked_result(
            original=original,
            left_core=provisional_left,
            mellin_core=mellin_core,
            right_core=provisional_right,
            status=FullSchurFactorizationStatus.BLOCKED_BY_REGULARIZER_INTERFACE,
            assumptions=assumptions,
            factors=tuple(duplicate_factors),
            evidence=(
                "typed embedded-interface metadata shows duplicated T1_minus or Z1_inverse",
            ),
            expanded=expanded,
        )

    if not isinstance(mellin_core, TransportedMellinCore) or (
        mellin_core.representation.status
        is RegularizerRepresentationStatus.ASSUMED
    ):
        return _blocked_result(
            original=original,
            left_core=provisional_left,
            mellin_core=mellin_core,
            right_core=provisional_right,
            status=FullSchurFactorizationStatus.BLOCKED_BY_REGULARIZER_INTERFACE,
            assumptions=assumptions,
            factors=(mellin_core,),
        )

    left_core = _analyze_left_core(left_block, left_auxiliary)
    right_core = _analyze_right_core(right_auxiliary, right_block)
    expanded = Product(
        (
            *(left_core.expanded_word.factors if left_core.expanded_word else ()),
            *mellin_core.transported_word.factors,
            *(right_core.expanded_word.factors if right_core.expanded_word else ()),
        )
    )
    if (
        left_core.status is CoreFactorizationStatus.BLOCKED
        or right_core.status is CoreFactorizationStatus.BLOCKED
    ):
        factors = (*left_core.blocking_factors, *right_core.blocking_factors)
        return _blocked_result(
            original=original,
            left_core=left_core,
            mellin_core=mellin_core,
            right_core=right_core,
            status=FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER,
            assumptions=assumptions,
            factors=factors,
            evidence=(left_core.evidence, right_core.evidence),
            expanded=expanded,
        )

    assert left_core.dilation is not None
    assert right_core.dilation is not None
    assert left_core.exterior is not None
    assert right_core.exterior is not None
    localized = tuple(
        item
        for item in (left_core.exterior, right_core.exterior)
        if isinstance(item, LocalizedOperator)
    )
    claims = _claims_for_localization((left_block, right_block), localized)
    if claims is None:
        return _blocked_result(
            original=original,
            left_core=left_core,
            mellin_core=mellin_core,
            right_core=right_core,
            status=FullSchurFactorizationStatus.BLOCKED_BY_CUTOFF_CONTROL,
            assumptions=assumptions,
            factors=localized,
            evidence=(left_core.relation, right_core.relation),
            expanded=expanded,
        )

    pdo = mellin_core.representation.mellin_operator
    if _same_space(
        left_core.exterior,
        left_core.dilation,
        pdo,
        right_core.dilation,
        right_core.exterior,
    ) is None:
        return _blocked_result(
            original=original,
            left_core=left_core,
            mellin_core=mellin_core,
            right_core=right_core,
            status=FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER,
            assumptions=assumptions,
            factors=(
                left_core.exterior,
                left_core.dilation,
                pdo,
                right_core.dilation,
                right_core.exterior,
            ),
            expanded=expanded,
        )
    try:
        covariance = conjugate_mellin_pdo_by_dilation(
            left_core.dilation,
            pdo,
            right_core.dilation,
            assumptions,
        )
    except (
        IncompatibleDomainError,
        MellinOperatorError,
        SchurDilationError,
    ) as exc:
        return _blocked_result(
            original=original,
            left_core=left_core,
            mellin_core=mellin_core,
            right_core=right_core,
            status=FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER,
            assumptions=assumptions,
            factors=(left_core.dilation, pdo, right_core.dilation, str(exc)),
            expanded=expanded,
        )
    scaled = covariance.transformed_operator
    if not isinstance(scaled, MellinPseudodifferentialOperator):
        raise FullSchurWordError("Mellin covariance returned a wrong type.")
    factorized = Product(
        (
            _atom(left_core.exterior, "left exterior"),
            scaled.atom,
            _atom(right_core.exterior, "right exterior"),
        )
    )

    cancellation: SchurDilationCancellationResult | None = None
    if isinstance(left_core.exterior, WienerHopfOperator) and isinstance(
        right_core.exterior, WienerHopfOperator
    ):
        combined_left = SchurBlockRepresentation(
            exact_block=ExactBlock("A", 2, 1),
            atom=OperatorAtom("left_full_word_core", kind="derived_block"),
            factors=(left_core.exterior, left_core.dilation),
            relation_kind=left_block.representation.relation_kind,
            source=left_block.representation.source,
            hypotheses=left_block.representation.hypotheses,
            evidence=left_block.representation.evidence,
            compact_ideal=left_block.representation.compact_ideal,
        )
        combined_right = SchurBlockRepresentation(
            exact_block=ExactBlock("A", 1, 2),
            atom=OperatorAtom("right_full_word_core", kind="derived_block"),
            factors=(right_core.dilation, right_core.exterior),
            relation_kind=right_block.representation.relation_kind,
            source=right_block.representation.source,
            hypotheses=right_block.representation.hypotheses,
            evidence=right_block.representation.evidence,
            compact_ideal=right_block.representation.compact_ideal,
        )
        cancellation = cancel_relative_dilations_in_schur_correction(
            left_block=combined_left,
            regularizer=mellin_core.representation,
            right_block=combined_right,
            assumptions=assumptions,
        )
        if not cancellation.succeeded:
            return _blocked_result(
                original=original,
                left_core=left_core,
                mellin_core=mellin_core,
                right_core=right_core,
                status=FullSchurFactorizationStatus.BLOCKED_BY_OPERATOR_ORDER,
                assumptions=assumptions,
                factors=cancellation.blocked_reasons,
                evidence=(cancellation,),
                expanded=expanded,
            )
        factorized = cancellation.expression

    claims_tuple = tuple(claims)
    any_assumed_cutoff = any(
        claim.status is CutoffEquivalenceStatus.ASSUMED for claim in claims_tuple
    )
    any_mod = any(
        block.representation.relation_kind
        is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
        for block in (left_block, right_block)
    ) or (
        mellin_core.representation.status
        is RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT
    )
    compact_ideal = _common_compact_ideal(
        left_block,
        mellin_core,
        right_block,
        claims_tuple,
    )
    evidence = (
        left_core.relation,
        mellin_core.semantic_relation,
        right_core.relation,
        covariance.exact_identity,
        *(claim.semantic_relation for claim in claims_tuple),
        *((cancellation,) if cancellation is not None else ()),
    )
    hypotheses = (
        *assumptions.assumptions,
        *left_block.representation.hypotheses,
        *mellin_core.representation.hypotheses,
        *right_block.representation.hypotheses,
        *(item for claim in claims_tuple for item in claim.assumptions),
    )
    relation: ExactIdentity | FormalIdentity | ModCompactEquivalence
    if claims_tuple:
        status = FullSchurFactorizationStatus.FACTORIZATION_AFTER_CUTOFF_EQUIVALENCE
        if any_assumed_cutoff:
            relation = FormalIdentity(
                original,
                factorized,
                justification=(evidence, hypotheses, "cutoff replacement assumed"),
            )
        else:
            if compact_ideal is None:
                compact_ideal = claims_tuple[0].ideal
            relation = ModCompactEquivalence(
                original,
                factorized,
                space=pdo.domain,
                compact_ideal=compact_ideal,
                evidence=evidence,
            )
    elif any_mod:
        status = FullSchurFactorizationStatus.FULL_WORD_MOD_COMPACTS_FACTORIZATION
        if compact_ideal is None:
            raise FullSchurWordError("a modulo-compact result lost its ideal.")
        relation = ModCompactEquivalence(
            original,
            factorized,
            space=pdo.domain,
            compact_ideal=compact_ideal,
            evidence=evidence,
        )
    else:
        status = FullSchurFactorizationStatus.FULL_WORD_EXACT_FACTORIZATION
        relation = ExactIdentity(
            original,
            factorized,
            evidence=evidence,
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=hypotheses,
        )
    obligations = _base_obligations()
    return FullSchurWordFactorizationResult(
        original,
        left_core,
        mellin_core,
        right_core,
        factorized,
        relation,
        status,
        evidence,
        assumptions,
        obligations,
        (),
        _render_latex(original, status, obligations, factorized),
        expanded,
        cancellation,
    )


@dataclass(frozen=True)
class RawSchurCorrection:
    left_block: OperatorAtom
    regularizer: OperatorAtom
    right_block: OperatorAtom
    intrinsic_sign: int = 1
    word: Product = field(init=False)

    def __post_init__(self) -> None:
        for name in ("left_block", "regularizer", "right_block"):
            if not isinstance(getattr(self, name), OperatorAtom):
                raise TypeError(f"{name} must be an OperatorAtom.")
        if isinstance(self.intrinsic_sign, bool) or self.intrinsic_sign != 1:
            raise FullSchurWordError(
                "RawSchurCorrection is intrinsically positive; store sign outside it."
            )
        object.__setattr__(
            self,
            "word",
            Product((self.left_block, self.regularizer, self.right_block)),
        )


@dataclass(frozen=True)
class SignedSchurContribution:
    raw: RawSchurCorrection
    sign: int
    term: Term = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.raw, RawSchurCorrection):
            raise TypeError("raw must be a RawSchurCorrection.")
        if isinstance(self.sign, bool) or self.sign not in {-1, 1}:
            raise FullSchurWordError("sign must be exactly +1 or -1.")
        if self.raw.intrinsic_sign != 1:
            raise FullSchurWordError("double-negated raw corrections are forbidden.")
        object.__setattr__(self, "term", Term(self.sign, self.raw.word))


@dataclass(frozen=True)
class SchurPivot:
    diagonal: OperatorAtom
    contribution: SignedSchurContribution
    expression: LinearCombination = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.diagonal, OperatorAtom):
            raise TypeError("diagonal must be an OperatorAtom.")
        if not isinstance(self.contribution, SignedSchurContribution):
            raise TypeError("contribution must be a SignedSchurContribution.")
        object.__setattr__(
            self,
            "expression",
            LinearCombination(
                (
                    Term(1, self.diagonal),
                    self.contribution.term,
                )
            ),
        )


__all__ = [
    "AuxiliaryOperator",
    "CauchyPolarity",
    "CauchyProjectionOperator",
    "CoefficientMultiplicationOperator",
    "CoreFactorizationStatus",
    "CutoffDilationCovarianceTrace",
    "CutoffEquivalenceClaim",
    "CutoffEquivalenceSide",
    "CutoffEquivalenceStatus",
    "FullSchurFactorizationStatus",
    "FullSchurWordError",
    "FullSchurWordFactorizationResult",
    "FullWordBlock",
    "LocalizedOperator",
    "RawSchurCorrection",
    "RegularizerInterfaceRole",
    "SchurExteriorCoreResult",
    "SchurPivot",
    "SignedSchurContribution",
    "SpatialCutoffOperator",
    "TransportedMellinCore",
    "TransportedShiftOperator",
    "conjugate_cutoff_by_dilation",
    "factor_full_first_schur_word",
]
