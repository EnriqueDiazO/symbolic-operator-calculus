"""Certified assembly and fail-closed algebra audit for the first Schur correction.

The exact assembly implemented here is deliberately independent of membership
in a Fredholm algebra.  The P1-C regularizer stays opaque in the grouped word,
while a separate evidence ledger records its source-backed three-factor
expansion.  No operation in this module commutes factors, distributes a linear
combination, removes a cutoff, or manufactures a symbol map.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import sympy as sp

from .diagonal_regularizer import (
    CompactIdeal,
    FactoredRegularizerStatus,
    FactoredTwoSidedRegularizerCertificate,
    IDENTITY_PRODUCT,
    InvertibleAuxiliaryOperator,
)
from .domains import AssumptionContext, ConsistencyStatus
from .mellin.operators import (
    MultiplicationOperator,
    WeightedDilationOperator,
    WeightedLpSpace,
    WienerHopfOperator,
)
from .operators import LinearCombination, OperatorAtom, Product, Term
from .schur_full_word import (
    CoefficientMultiplicationOperator,
    CutoffDilationCovarianceTrace,
    LocalizedOperator,
    SpatialCutoffOperator,
    TransportedShiftOperator,
)
from .semantics import (
    AnalyticProofObligation,
    ExactIdentity,
    ExactIdentityScope,
    FormalIdentity,
    ModCompactEquivalence,
)


class SchurCorrectionError(ValueError):
    """Base error for malformed or overclaimed P1-D data."""


class SchurDomainOrOrderError(SchurCorrectionError):
    """Raised when a typed branch, space, role, or product order is invalid."""


class CutoffReplacementNotCertified(SchurCorrectionError):
    """Raised when a requested cutoff replacement has no certified theorem."""


class SchurSymbolNotAvailable(SchurCorrectionError):
    """Raised when a blocked symbol result is requested as certified."""


class AlgebraMembershipNotCertified(SchurSymbolNotAvailable):
    """The complete correction is not certified in a suitable common algebra."""


class CutoffTreatmentNotCertified(SchurSymbolNotAvailable):
    """At least one cutoff operation required by the symbol is unproved."""


class RelativeDilationNotResolved(SchurSymbolNotAvailable):
    """A required relative-dilation reduction remains blocked or conditional."""


class NoncommutativeProductRuleMissing(SchurSymbolNotAvailable):
    """No certified ordered-product rule is available for the symbol map."""


class SchurAssemblyStatus(str, Enum):
    CERTIFIED_EXACT_SCHUR_ASSEMBLY = "CERTIFIED_EXACT_SCHUR_ASSEMBLY"
    CERTIFIED_MOD_COMPACTS_SCHUR_ASSEMBLY = (
        "CERTIFIED_MOD_COMPACTS_SCHUR_ASSEMBLY"
    )
    CONDITIONAL_SCHUR_ASSEMBLY = "CONDITIONAL_SCHUR_ASSEMBLY"
    BLOCKED_BY_SOURCE_ORDER_MISMATCH = "BLOCKED_BY_SOURCE_ORDER_MISMATCH"
    BLOCKED_BY_DOMAIN_OR_CODOMAIN = "BLOCKED_BY_DOMAIN_OR_CODOMAIN"
    BLOCKED_BY_SIGN_CONVENTION = "BLOCKED_BY_SIGN_CONVENTION"
    BLOCKED_BY_MISSING_FACTOR = "BLOCKED_BY_MISSING_FACTOR"


class AlgebraMembershipStatus(str, Enum):
    CERTIFIED_COMMON_ALGEBRA_MEMBERSHIP = (
        "CERTIFIED_COMMON_ALGEBRA_MEMBERSHIP"
    )
    CERTIFIED_MEMBERSHIP_IN_EXPLICIT_EXTENSION = (
        "CERTIFIED_MEMBERSHIP_IN_EXPLICIT_EXTENSION"
    )
    CONDITIONAL_ALGEBRA_MEMBERSHIP = "CONDITIONAL_ALGEBRA_MEMBERSHIP"
    BLOCKED_BY_CUTOFF_REPLACEMENT = "BLOCKED_BY_CUTOFF_REPLACEMENT"
    BLOCKED_BY_RELATIVE_DILATION = "BLOCKED_BY_RELATIVE_DILATION"
    BLOCKED_BY_OFF_DIAGONAL_WH_MEMBERSHIP = (
        "BLOCKED_BY_OFF_DIAGONAL_WH_MEMBERSHIP"
    )
    BLOCKED_BY_CENTRAL_REGULARIZER_MEMBERSHIP = (
        "BLOCKED_BY_CENTRAL_REGULARIZER_MEMBERSHIP"
    )
    BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP = (
        "BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP"
    )
    NOT_ATTEMPTED_BECAUSE_ASSEMBLY_FAILED = (
        "NOT_ATTEMPTED_BECAUSE_ASSEMBLY_FAILED"
    )


class SchurCorrectionSymbolStatus(str, Enum):
    CERTIFIED_SCHUR_CORRECTION_SYMBOL = "CERTIFIED_SCHUR_CORRECTION_SYMBOL"
    CONDITIONAL_SCHUR_CORRECTION_SYMBOL = "CONDITIONAL_SCHUR_CORRECTION_SYMBOL"
    BLOCKED_BY_ALGEBRA_MEMBERSHIP = "BLOCKED_BY_ALGEBRA_MEMBERSHIP"
    BLOCKED_BY_CUTOFFS = "BLOCKED_BY_CUTOFFS"
    BLOCKED_BY_NONCOMMUTATIVE_PRODUCT = "BLOCKED_BY_NONCOMMUTATIVE_PRODUCT"
    NOT_ATTEMPTED = "NOT_ATTEMPTED"


class CutoffTreatmentStatus(str, Enum):
    EXACT_COVARIANCE_ONLY = "EXACT_COVARIANCE_ONLY"
    MOD_COMPACTS_REPLACEMENT_CERTIFIED = (
        "MOD_COMPACTS_REPLACEMENT_CERTIFIED"
    )
    CONDITIONAL_REPLACEMENT = "CONDITIONAL_REPLACEMENT"
    REPLACEMENT_BLOCKED = "REPLACEMENT_BLOCKED"
    NO_REPLACEMENT_NEEDED = "NO_REPLACEMENT_NEEDED"


class RelativeDilationStatus(str, Enum):
    EXACT_ADJACENT_CANCELLATION = "EXACT_ADJACENT_CANCELLATION"
    EXACT_COVARIANCE_REDUCTION = "EXACT_COVARIANCE_REDUCTION"
    CONDITIONAL_REDUCTION = "CONDITIONAL_REDUCTION"
    BLOCKED_BY_INTERVENING_FACTORS = "BLOCKED_BY_INTERVENING_FACTORS"
    NO_RELATIVE_DILATION_PAIR = "NO_RELATIVE_DILATION_PAIR"


class FactorCoverageClassification(str, Enum):
    EXACT_THEOREM = "EXACT_THEOREM"
    DIRECT_ALGEBRAIC_CONSEQUENCE = "DIRECT_ALGEBRAIC_CONSEQUENCE"
    MATCH_AFTER_NOTATION_TRANSLATION = "MATCH_AFTER_NOTATION_TRANSLATION"
    PLAUSIBLE_ADAPTATION = "PLAUSIBLE_ADAPTATION"
    NOT_COVERED = "NOT_COVERED"
    CONTRADICTED_BY_HYPOTHESES = "CONTRADICTED_BY_HYPOTHESES"
    SOURCE_GAP = "SOURCE_GAP"


class ExteriorFactorStatus(str, Enum):
    EXACT_SOURCE_FACTOR = "EXACT_SOURCE_FACTOR"
    CONDITIONAL_SOURCE_FACTOR = "CONDITIONAL_SOURCE_FACTOR"
    UNVERIFIED_SOURCE_FACTOR = "UNVERIFIED_SOURCE_FACTOR"


class SourceBlockSide(str, Enum):
    LEFT = "LEFT"
    RIGHT = "RIGHT"


class CutoffSide(str, Enum):
    LEFT = "LEFT"
    RIGHT = "RIGHT"


class PivotFredholmStatus(str, Enum):
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


def _operator_space(value: object, field_name: str) -> WeightedLpSpace:
    domain = getattr(value, "domain", None)
    codomain = getattr(value, "codomain", None)
    if not isinstance(domain, WeightedLpSpace) or not isinstance(
        codomain, WeightedLpSpace
    ):
        raise TypeError(
            f"{field_name} must expose WeightedLpSpace domain and codomain."
        )
    if domain != codomain:
        raise SchurDomainOrOrderError(
            f"{field_name} must preserve one scalar weighted space."
        )
    return domain


def _require_bounded(value: object, field_name: str) -> None:
    if getattr(value, "bounded", None) is not True:
        raise SchurDomainOrOrderError(f"{field_name} must be explicitly bounded.")


def _require_assumptions(value: object) -> AssumptionContext:
    if not isinstance(value, AssumptionContext):
        raise TypeError("assumptions must be an AssumptionContext.")
    if value.consistency_status is ConsistencyStatus.INCONSISTENT:
        raise SchurDomainOrOrderError("assumptions must be consistent.")
    return value


def _as_product(expression: OperatorAtom | Product, field_name: str) -> Product:
    if isinstance(expression, OperatorAtom):
        return Product((expression,))
    if isinstance(expression, Product):
        return expression
    raise TypeError(f"{field_name} must be an OperatorAtom or Product.")


@dataclass(frozen=True, order=True)
class Branch:
    """A branch identifier whose orientation is never inferred from a name."""

    index: int

    def __post_init__(self) -> None:
        if not isinstance(self.index, int) or isinstance(self.index, bool):
            raise TypeError("branch index must be an integer.")
        if self.index <= 0:
            raise ValueError("branch index must be positive.")


@dataclass(frozen=True)
class BranchedWeightedSpace:
    branch: Branch
    space: WeightedLpSpace

    def __post_init__(self) -> None:
        if not isinstance(self.branch, Branch):
            raise TypeError("branch must be a Branch.")
        if not isinstance(self.space, WeightedLpSpace):
            raise TypeError("space must be a WeightedLpSpace.")


@dataclass(frozen=True)
class BranchDifferenceOperator:
    """Typed exact difference ``transported shift - coefficient``."""

    atom: OperatorAtom
    transported_shift: TransportedShiftOperator
    coefficient: CoefficientMultiplicationOperator
    source: str
    assumptions: AssumptionContext
    evidence: tuple[object, ...]
    status: ExteriorFactorStatus = ExteriorFactorStatus.EXACT_SOURCE_FACTOR
    expanded_expression: LinearCombination = field(init=False)
    exact_definition: ExactIdentity | FormalIdentity | None = field(init=False)
    domain: WeightedLpSpace = field(init=False)
    codomain: WeightedLpSpace = field(init=False)
    bounded: bool = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.transported_shift, TransportedShiftOperator):
            raise TypeError("transported_shift must be a TransportedShiftOperator.")
        if not isinstance(self.coefficient, CoefficientMultiplicationOperator):
            raise TypeError("coefficient must be a CoefficientMultiplicationOperator.")
        space = _operator_space(self.transported_shift, "transported_shift")
        if _operator_space(self.coefficient, "coefficient") != space:
            raise SchurDomainOrOrderError(
                "the two difference terms act on different weighted spaces."
            )
        assumptions = _require_assumptions(self.assumptions)
        source = _nonempty(self.source, "source")
        evidence = _items(self.evidence, "evidence")
        if not isinstance(self.status, ExteriorFactorStatus):
            raise TypeError("status must be an ExteriorFactorStatus.")
        if self.status is ExteriorFactorStatus.EXACT_SOURCE_FACTOR and not evidence:
            raise SchurCorrectionError("an exact branch difference requires evidence.")
        expression = LinearCombination(
            (
                Term(1, self.transported_shift.atom),
                Term(-1, self.coefficient.atom),
            )
        )
        relation: ExactIdentity | FormalIdentity | None
        if self.status is ExteriorFactorStatus.EXACT_SOURCE_FACTOR:
            relation = ExactIdentity(
                self.atom,
                expression,
                evidence=(source, evidence),
                scope=ExactIdentityScope.OPERATORIAL,
                hypotheses=assumptions.assumptions,
            )
        elif self.status is ExteriorFactorStatus.CONDITIONAL_SOURCE_FACTOR:
            if assumptions.is_empty:
                raise SchurCorrectionError(
                    "a conditional branch difference requires assumptions."
                )
            relation = FormalIdentity(
                self.atom,
                expression,
                justification=(source, assumptions, evidence),
            )
        else:
            relation = None
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "expanded_expression", expression)
        object.__setattr__(self, "exact_definition", relation)
        object.__setattr__(self, "domain", space)
        object.__setattr__(self, "codomain", space)
        object.__setattr__(
            self,
            "bounded",
            self.transported_shift.bounded and self.coefficient.bounded,
        )


@dataclass(frozen=True)
class SchurExteriorFactor:
    """One branch-oriented exterior ``(Vtilde-G) W_localized``."""

    atom: OperatorAtom
    difference: BranchDifferenceOperator
    localized_wiener_hopf: LocalizedOperator
    branch_from: Branch
    branch_to: Branch
    source: str
    assumptions: AssumptionContext
    evidence: tuple[object, ...]
    status: ExteriorFactorStatus = ExteriorFactorStatus.EXACT_SOURCE_FACTOR
    source_expression: Product = field(init=False)
    expanded_expression: Product = field(init=False)
    domain: WeightedLpSpace = field(init=False)
    codomain: WeightedLpSpace = field(init=False)
    weight: sp.Expr = field(init=False)
    localizers: tuple[SpatialCutoffOperator, SpatialCutoffOperator] = field(init=False)
    relative_dilations: tuple[WeightedDilationOperator, ...] = field(init=False)
    wiener_hopf_symbols: tuple[sp.Expr, ...] = field(init=False)
    multipliers: tuple[CoefficientMultiplicationOperator, ...] = field(init=False)
    bounded: bool = field(init=False)
    exact_definition: ExactIdentity | FormalIdentity | None = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.difference, BranchDifferenceOperator):
            raise TypeError("difference must be a BranchDifferenceOperator.")
        if not isinstance(self.localized_wiener_hopf, LocalizedOperator):
            raise TypeError("localized_wiener_hopf must be a LocalizedOperator.")
        if not isinstance(self.branch_from, Branch) or not isinstance(
            self.branch_to, Branch
        ):
            raise TypeError("branch_from and branch_to must be Branch objects.")
        if self.branch_from == self.branch_to:
            raise SchurDomainOrOrderError("a Schur exterior must change branch.")
        space = _operator_space(self.difference, "difference")
        if _operator_space(self.localized_wiener_hopf, "localized_wiener_hopf") != space:
            raise SchurDomainOrOrderError(
                "the difference and Wiener--Hopf block use different spaces."
            )
        assumptions = _require_assumptions(self.assumptions)
        source = _nonempty(self.source, "source")
        evidence = _items(self.evidence, "evidence")
        if not isinstance(self.status, ExteriorFactorStatus):
            raise TypeError("status must be an ExteriorFactorStatus.")
        expression = Product((self.difference.atom, self.localized_wiener_hopf.atom))
        relation: ExactIdentity | FormalIdentity | None
        if self.status is ExteriorFactorStatus.EXACT_SOURCE_FACTOR:
            if not evidence or not isinstance(
                self.difference.exact_definition, ExactIdentity
            ):
                raise SchurCorrectionError(
                    "an exact exterior requires exact difference evidence."
                )
            relation = ExactIdentity(
                self.atom,
                expression,
                evidence=(source, evidence, self.difference.exact_definition),
                scope=ExactIdentityScope.OPERATORIAL,
                hypotheses=assumptions.assumptions,
            )
        elif self.status is ExteriorFactorStatus.CONDITIONAL_SOURCE_FACTOR:
            if assumptions.is_empty:
                raise SchurCorrectionError(
                    "a conditional exterior requires explicit assumptions."
                )
            relation = FormalIdentity(
                self.atom,
                expression,
                justification=(source, assumptions, evidence),
            )
        else:
            relation = None
        core = self.localized_wiener_hopf.core
        wh_symbols = (core.symbol,) if isinstance(core, WienerHopfOperator) else ()
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "source_expression", expression)
        object.__setattr__(self, "expanded_expression", expression)
        object.__setattr__(self, "domain", space)
        object.__setattr__(self, "codomain", space)
        object.__setattr__(self, "weight", space.weight_exponent)
        object.__setattr__(
            self,
            "localizers",
            (
                self.localized_wiener_hopf.left_cutoff,
                self.localized_wiener_hopf.right_cutoff,
            ),
        )
        object.__setattr__(
            self,
            "relative_dilations",
            (self.difference.transported_shift.dilation,),
        )
        object.__setattr__(self, "wiener_hopf_symbols", wh_symbols)
        object.__setattr__(
            self,
            "multipliers",
            (
                self.difference.transported_shift.coefficient,
                self.difference.coefficient,
            ),
        )
        object.__setattr__(
            self,
            "bounded",
            self.difference.bounded and self.localized_wiener_hopf.bounded,
        )
        object.__setattr__(self, "exact_definition", relation)

    @property
    def branch_domain(self) -> BranchedWeightedSpace:
        return BranchedWeightedSpace(self.branch_from, self.domain)

    @property
    def branch_codomain(self) -> BranchedWeightedSpace:
        return BranchedWeightedSpace(self.branch_to, self.codomain)


@dataclass(frozen=True)
class SchurSourceBlock:
    """Article factor ``B^- = L T`` or ``B^+ = Z^-1 R``."""

    atom: OperatorAtom
    side: SourceBlockSide
    exterior: SchurExteriorFactor
    interface: object
    source: str
    evidence: tuple[object, ...]
    status: ExteriorFactorStatus = ExteriorFactorStatus.EXACT_SOURCE_FACTOR
    expanded_expression: Product = field(init=False)
    exact_definition: ExactIdentity | FormalIdentity | None = field(init=False)
    domain: WeightedLpSpace = field(init=False)
    codomain: WeightedLpSpace = field(init=False)
    branch_from: Branch = field(init=False)
    branch_to: Branch = field(init=False)
    bounded: bool = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.side, SourceBlockSide):
            raise TypeError("side must be a SourceBlockSide.")
        if not isinstance(self.exterior, SchurExteriorFactor):
            raise TypeError("exterior must be a SchurExteriorFactor.")
        interface_atom = _atom(self.interface, "interface")
        if _operator_space(self.interface, "interface") != self.exterior.domain:
            raise SchurDomainOrOrderError(
                "source-block interface acts on a different weighted space."
            )
        _require_bounded(self.interface, "source-block interface")
        source = _nonempty(self.source, "source")
        evidence = _items(self.evidence, "evidence")
        if not isinstance(self.status, ExteriorFactorStatus):
            raise TypeError("status must be an ExteriorFactorStatus.")
        if self.side is SourceBlockSide.LEFT:
            if not isinstance(self.interface, InvertibleAuxiliaryOperator):
                raise TypeError("the left source block requires the T interface.")
            expression = Product((self.exterior.atom, interface_atom))
        else:
            if not isinstance(self.interface, MultiplicationOperator):
                raise TypeError("the right source block requires the Z inverse multiplier.")
            expression = Product((interface_atom, self.exterior.atom))
        relation: ExactIdentity | FormalIdentity | None
        if self.status is ExteriorFactorStatus.EXACT_SOURCE_FACTOR:
            if not evidence:
                raise SchurCorrectionError("an exact source block requires evidence.")
            relation = ExactIdentity(
                self.atom,
                expression,
                evidence=(source, evidence, self.exterior.exact_definition),
                scope=ExactIdentityScope.OPERATORIAL,
                hypotheses=self.exterior.assumptions.assumptions,
            )
        elif self.status is ExteriorFactorStatus.CONDITIONAL_SOURCE_FACTOR:
            relation = FormalIdentity(
                self.atom,
                expression,
                justification=(source, evidence, self.exterior.assumptions),
            )
        else:
            relation = None
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "expanded_expression", expression)
        object.__setattr__(self, "exact_definition", relation)
        object.__setattr__(self, "domain", self.exterior.domain)
        object.__setattr__(self, "codomain", self.exterior.codomain)
        object.__setattr__(self, "branch_from", self.exterior.branch_from)
        object.__setattr__(self, "branch_to", self.exterior.branch_to)
        object.__setattr__(
            self,
            "bounded",
            self.exterior.bounded and getattr(self.interface, "bounded", False),
        )


@dataclass(frozen=True)
class AssemblyTraceStep:
    key: str
    input_expression: Product
    output_expression: Product
    rule: str
    relation: ExactIdentity | FormalIdentity | ModCompactEquivalence
    evidence: tuple[object, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "key", _nonempty(self.key, "key"))
        object.__setattr__(self, "rule", _nonempty(self.rule, "rule"))
        if not isinstance(self.input_expression, Product) or not isinstance(
            self.output_expression, Product
        ):
            raise TypeError("trace endpoints must be Product objects.")
        if not isinstance(
            self.relation, (ExactIdentity, FormalIdentity, ModCompactEquivalence)
        ):
            raise TypeError("trace relation has an unsupported type.")
        if (self.relation.left, self.relation.right) != (
            self.input_expression,
            self.output_expression,
        ):
            raise SchurDomainOrOrderError(
                "trace relation endpoints disagree with the trace."
            )
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise SchurCorrectionError("every assembly trace needs evidence.")
        object.__setattr__(self, "evidence", evidence)


@dataclass(frozen=True)
class SchurCorrectionWordLedger:
    raw_schur_word: Product | None
    source_factorized_word: Product | None
    grouped_regularizer_word: Product | None
    expanded_word: Product | None

    def __post_init__(self) -> None:
        for name in (
            "raw_schur_word",
            "source_factorized_word",
            "grouped_regularizer_word",
            "expanded_word",
        ):
            value = getattr(self, name)
            if value is not None and not isinstance(value, Product):
                raise TypeError(f"{name} must be a Product or None.")


def _schur_obligations() -> tuple[AnalyticProofObligation, ...]:
    statements = (
        (
            "full_regularizer_algebra_membership",
            "Prove common-algebra membership for the complete P1-C regularizer.",
        ),
        (
            "cutoff_replacement_mod_compacts",
            "Prove each cutoff replacement if a later derivation requests one.",
        ),
        (
            "off_diagonal_wh_membership",
            "Retain source-backed membership of both localized Wiener--Hopf blocks.",
        ),
        (
            "relative_dilation_membership",
            "Prove membership of nontrivial relative dilations in a common algebra.",
        ),
        (
            "wh_mellin_wh_sandwich_membership",
            "Prove closure for the exact ordered WH--Mellin--WH sandwich.",
        ),
        (
            "fredholm_algebra_membership",
            "Prove membership in one inverse-closed algebra with a Fredholm theorem.",
        ),
        (
            "schur_correction_symbol",
            "Construct a symbol only after every algebraic gate is certified.",
        ),
        (
            "first_schur_pivot_fredholmness",
            "Apply a valid pivot Fredholm criterion; assembly alone is insufficient.",
        ),
    )
    return tuple(AnalyticProofObligation(key, statement) for key, statement in statements)


@dataclass(frozen=True)
class SchurCorrectionAssemblyCertificate:
    left_factor: SchurExteriorFactor | None
    central_regularizer: FactoredTwoSidedRegularizerCertificate | None
    right_factor: SchurExteriorFactor | None
    source_blocks: tuple[SchurSourceBlock, ...]
    correction_atom: OperatorAtom
    grouped_expression: Product | None
    expanded_expression: Product | None
    source_expression: Product | None
    grouped_to_expanded_trace: tuple[AssemblyTraceStep, ...]
    domain: BranchedWeightedSpace | None
    codomain: BranchedWeightedSpace | None
    compact_ideal: CompactIdeal | None
    assembly_relation: ExactIdentity | FormalIdentity | ModCompactEquivalence | None
    source_relation: ExactIdentity | FormalIdentity | ModCompactEquivalence | None
    correction_definition: ExactIdentity | FormalIdentity | None
    assembly_status: SchurAssemblyStatus
    assumptions: AssumptionContext
    evidence: tuple[object, ...]
    proof_obligations: tuple[AnalyticProofObligation, ...]
    blocking_factors: tuple[object, ...]
    word_ledger: SchurCorrectionWordLedger
    latex: str

    def __post_init__(self) -> None:
        if not isinstance(self.correction_atom, OperatorAtom):
            raise TypeError("correction_atom must be an OperatorAtom.")
        if not isinstance(self.assembly_status, SchurAssemblyStatus):
            raise TypeError("assembly_status must be a SchurAssemblyStatus.")
        _require_assumptions(self.assumptions)
        for name in ("grouped_expression", "expanded_expression", "source_expression"):
            value = getattr(self, name)
            if value is not None and not isinstance(value, Product):
                raise TypeError(f"{name} must be a Product or None.")
        if not all(
            isinstance(item, SchurSourceBlock) for item in self.source_blocks
        ):
            raise TypeError("source_blocks contains an invalid item.")
        if not all(
            isinstance(item, AssemblyTraceStep)
            for item in self.grouped_to_expanded_trace
        ):
            raise TypeError("grouped_to_expanded_trace contains an invalid item.")
        if not all(
            isinstance(item, AnalyticProofObligation)
            for item in self.proof_obligations
        ):
            raise TypeError("proof_obligations contains an invalid item.")
        if not isinstance(self.word_ledger, SchurCorrectionWordLedger):
            raise TypeError("word_ledger must be a SchurCorrectionWordLedger.")
        successful = self.assembly_status in {
            SchurAssemblyStatus.CERTIFIED_EXACT_SCHUR_ASSEMBLY,
            SchurAssemblyStatus.CERTIFIED_MOD_COMPACTS_SCHUR_ASSEMBLY,
            SchurAssemblyStatus.CONDITIONAL_SCHUR_ASSEMBLY,
        }
        if successful and (
            self.grouped_expression is None
            or self.expanded_expression is None
            or self.source_expression is None
            or self.assembly_relation is None
            or self.source_relation is None
            or self.correction_definition is None
            or self.domain is None
            or self.codomain is None
        ):
            raise SchurCorrectionError(
                "a successful assembly requires complete expressions and relations."
            )
        if not successful and (
            self.assembly_relation is not None
            or self.source_relation is not None
            or self.correction_definition is not None
        ):
            raise SchurCorrectionError(
                "a blocked assembly cannot carry a derived semantic relation."
            )
        object.__setattr__(self, "source_blocks", tuple(self.source_blocks))
        object.__setattr__(
            self,
            "grouped_to_expanded_trace",
            tuple(self.grouped_to_expanded_trace),
        )
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))
        object.__setattr__(
            self, "blocking_factors", _items(self.blocking_factors, "blocking_factors")
        )
        object.__setattr__(self, "latex", _nonempty(self.latex, "latex"))


def _blocked_assembly(
    *,
    status: SchurAssemblyStatus,
    assumptions: AssumptionContext,
    left_factor: SchurExteriorFactor | None,
    central_regularizer: FactoredTwoSidedRegularizerCertificate | None,
    right_factor: SchurExteriorFactor | None,
    source_blocks: tuple[SchurSourceBlock, ...],
    blockers: tuple[object, ...],
    grouped: Product | None = None,
    expanded: Product | None = None,
    source: Product | None = None,
) -> SchurCorrectionAssemblyCertificate:
    return SchurCorrectionAssemblyCertificate(
        left_factor=left_factor,
        central_regularizer=central_regularizer,
        right_factor=right_factor,
        source_blocks=source_blocks,
        correction_atom=OperatorAtom("C22_1", kind="schur_correction"),
        grouped_expression=grouped,
        expanded_expression=expanded,
        source_expression=source,
        grouped_to_expanded_trace=(),
        domain=None,
        codomain=None,
        compact_ideal=(
            central_regularizer.compact_ideal
            if central_regularizer is not None
            else None
        ),
        assembly_relation=None,
        source_relation=None,
        correction_definition=None,
        assembly_status=status,
        assumptions=assumptions,
        evidence=(),
        proof_obligations=_schur_obligations(),
        blocking_factors=blockers,
        word_ledger=SchurCorrectionWordLedger(None, source, grouped, expanded),
        latex=rf"\mathcal C_{{2,2}}^{{(1)}}\quad\text{{{status.value}}}",
    )


def assemble_first_schur_correction(
    *,
    left_factor: SchurExteriorFactor | None,
    diagonal_regularizer: FactoredTwoSidedRegularizerCertificate | None,
    right_factor: SchurExteriorFactor | None,
    source_blocks: tuple[SchurSourceBlock, ...],
    assumptions: AssumptionContext,
) -> SchurCorrectionAssemblyCertificate:
    """Assemble the article correction by exact structural substitution only."""

    assumptions = _require_assumptions(assumptions)
    blocks = _items(source_blocks, "source_blocks")
    typed_blocks = tuple(item for item in blocks if isinstance(item, SchurSourceBlock))
    if (
        left_factor is None
        or diagonal_regularizer is None
        or right_factor is None
        or len(blocks) != 2
        or len(typed_blocks) != 2
    ):
        return _blocked_assembly(
            status=SchurAssemblyStatus.BLOCKED_BY_MISSING_FACTOR,
            assumptions=assumptions,
            left_factor=left_factor,
            central_regularizer=diagonal_regularizer,
            right_factor=right_factor,
            source_blocks=typed_blocks,
            blockers=tuple(
                item
                for item in (left_factor, diagonal_regularizer, right_factor, *blocks)
                if item is None or not isinstance(
                    item,
                    (
                        SchurExteriorFactor,
                        FactoredTwoSidedRegularizerCertificate,
                        SchurSourceBlock,
                    ),
                )
            ),
        )
    if not isinstance(left_factor, SchurExteriorFactor) or not isinstance(
        right_factor, SchurExteriorFactor
    ):
        raise TypeError("left_factor and right_factor must be SchurExteriorFactor objects.")
    if not isinstance(
        diagonal_regularizer, FactoredTwoSidedRegularizerCertificate
    ):
        raise TypeError("diagonal_regularizer must be a P1-C certificate.")
    left_source, right_source = typed_blocks
    candidate = diagonal_regularizer.candidate_regularizer
    grouped = Product((left_factor.atom, candidate.atom, right_factor.atom))
    expanded = Product(
        (
            *left_factor.expanded_expression.factors,
            *candidate.ast_product.factors,
            *right_factor.expanded_expression.factors,
        )
    )
    source = Product(
        (
            left_source.atom,
            candidate.mellin_regularizer.atom,
            right_source.atom,
        )
    )
    source_expanded = Product(
        (
            *left_source.exterior.expanded_expression.factors,
            _atom(left_source.interface, "left source interface"),
            candidate.mellin_regularizer.atom,
            _atom(right_source.interface, "right source interface"),
            *right_source.exterior.expanded_expression.factors,
        )
    )
    expected_branch = Branch(1)
    domain_ok = (
        right_factor.branch_from == Branch(2)
        and right_factor.branch_to == expected_branch
        and left_factor.branch_from == expected_branch
        and left_factor.branch_to == Branch(2)
        and left_factor.domain == right_factor.domain
        and candidate.domain == left_factor.domain
        and diagonal_regularizer.operator.domain == left_factor.domain
        and diagonal_regularizer.compact_ideal is not None
        and diagonal_regularizer.compact_ideal.space == left_factor.domain
    )
    if not domain_ok:
        return _blocked_assembly(
            status=SchurAssemblyStatus.BLOCKED_BY_DOMAIN_OR_CODOMAIN,
            assumptions=assumptions,
            left_factor=left_factor,
            central_regularizer=diagonal_regularizer,
            right_factor=right_factor,
            source_blocks=typed_blocks,
            blockers=(left_factor, diagonal_regularizer, right_factor),
            grouped=grouped,
            expanded=expanded,
            source=source,
        )
    source_ok = (
        left_source.side is SourceBlockSide.LEFT
        and right_source.side is SourceBlockSide.RIGHT
        and left_source.exterior == left_factor
        and right_source.exterior == right_factor
        and left_source.interface == candidate.left_interface
        and right_source.interface == candidate.right_interface_inverse
        and source_expanded == expanded
    )
    if not source_ok:
        return _blocked_assembly(
            status=SchurAssemblyStatus.BLOCKED_BY_SOURCE_ORDER_MISMATCH,
            assumptions=assumptions,
            left_factor=left_factor,
            central_regularizer=diagonal_regularizer,
            right_factor=right_factor,
            source_blocks=typed_blocks,
            blockers=typed_blocks,
            grouped=grouped,
            expanded=expanded,
            source=source,
        )
    if diagonal_regularizer.status is not (
        FactoredRegularizerStatus.CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER
    ):
        status = (
            SchurAssemblyStatus.CONDITIONAL_SCHUR_ASSEMBLY
            if diagonal_regularizer.status
            is FactoredRegularizerStatus.CONDITIONAL_TWO_SIDED_FACTORED_REGULARIZER
            else SchurAssemblyStatus.BLOCKED_BY_MISSING_FACTOR
        )
        if status is not SchurAssemblyStatus.CONDITIONAL_SCHUR_ASSEMBLY:
            return _blocked_assembly(
                status=status,
                assumptions=assumptions,
                left_factor=left_factor,
                central_regularizer=diagonal_regularizer,
                right_factor=right_factor,
                source_blocks=typed_blocks,
                blockers=(diagonal_regularizer,),
                grouped=grouped,
                expanded=expanded,
                source=source,
            )
    statuses = {
        left_factor.status,
        right_factor.status,
        left_source.status,
        right_source.status,
    }
    if ExteriorFactorStatus.UNVERIFIED_SOURCE_FACTOR in statuses:
        return _blocked_assembly(
            status=SchurAssemblyStatus.BLOCKED_BY_MISSING_FACTOR,
            assumptions=assumptions,
            left_factor=left_factor,
            central_regularizer=diagonal_regularizer,
            right_factor=right_factor,
            source_blocks=typed_blocks,
            blockers=(left_factor, right_factor, *typed_blocks),
            grouped=grouped,
            expanded=expanded,
            source=source,
        )
    conditional = (
        ExteriorFactorStatus.CONDITIONAL_SOURCE_FACTOR in statuses
        or diagonal_regularizer.status
        is FactoredRegularizerStatus.CONDITIONAL_TWO_SIDED_FACTORED_REGULARIZER
    )
    assembly_status = (
        SchurAssemblyStatus.CONDITIONAL_SCHUR_ASSEMBLY
        if conditional
        else SchurAssemblyStatus.CERTIFIED_EXACT_SCHUR_ASSEMBLY
    )
    if conditional:
        grouped_relation: ExactIdentity | FormalIdentity = FormalIdentity(
            grouped,
            expanded,
            justification=(assumptions, left_factor, diagonal_regularizer, right_factor),
        )
        source_relation: ExactIdentity | FormalIdentity = FormalIdentity(
            source,
            expanded,
            justification=(assumptions, left_source, right_source),
        )
        correction_definition: ExactIdentity | FormalIdentity = FormalIdentity(
            OperatorAtom("C22_1", kind="schur_correction"),
            grouped,
            justification=(assumptions, "conditional model-correction definition"),
        )
    else:
        grouped_relation = ExactIdentity(
            grouped,
            expanded,
            evidence=(
                left_factor.exact_definition,
                candidate.ast_product,
                right_factor.exact_definition,
                "structural substitution without commutation or distribution",
            ),
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=assumptions.assumptions,
        )
        source_relation = ExactIdentity(
            source,
            expanded,
            evidence=(
                left_source.exact_definition,
                right_source.exact_definition,
                "article B-minus/R1/B-plus definitions",
            ),
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=assumptions.assumptions,
        )
        correction_definition = ExactIdentity(
            OperatorAtom("C22_1", kind="schur_correction"),
            grouped,
            evidence=(
                "article exact model-correction definition",
                grouped_relation,
                source_relation,
            ),
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=assumptions.assumptions,
        )
    trace = (
        AssemblyTraceStep(
            key="expand_grouped_regularizer_word",
            input_expression=grouped,
            output_expression=expanded,
            rule="substitute exact exterior and P1-C candidate definitions in place",
            relation=grouped_relation,
            evidence=(left_factor, diagonal_regularizer, right_factor),
        ),
    )
    correction_atom = OperatorAtom("C22_1", kind="schur_correction")
    return SchurCorrectionAssemblyCertificate(
        left_factor=left_factor,
        central_regularizer=diagonal_regularizer,
        right_factor=right_factor,
        source_blocks=typed_blocks,
        correction_atom=correction_atom,
        grouped_expression=grouped,
        expanded_expression=expanded,
        source_expression=source,
        grouped_to_expanded_trace=trace,
        domain=BranchedWeightedSpace(Branch(2), right_factor.domain),
        codomain=BranchedWeightedSpace(Branch(2), left_factor.codomain),
        compact_ideal=diagonal_regularizer.compact_ideal,
        assembly_relation=grouped_relation,
        source_relation=source_relation,
        correction_definition=correction_definition,
        assembly_status=assembly_status,
        assumptions=assumptions,
        evidence=(
            left_factor,
            diagonal_regularizer,
            right_factor,
            left_source,
            right_source,
        ),
        proof_obligations=_schur_obligations(),
        blocking_factors=(),
        word_ledger=SchurCorrectionWordLedger(None, source, grouped, expanded),
        latex=(
            r"\mathcal C_{2,2}^{(1)}="
            r"\mathcal L_{2,1}"
            r"(\mathcal T_{1,-}\mathcal R_1Z_1^{-1})"
            r"\mathcal R_{1,2}"
            rf"\quad\text{{{assembly_status.value}}}"
        ),
    )


@dataclass(frozen=True)
class CutoffAuditRecord:
    cutoff: SpatialCutoffOperator
    cutoff_function: sp.Expr
    argument: sp.Expr
    scale: sp.Expr
    branch: Branch
    side: CutoffSide
    adjacent_operator: object
    support_claim: str
    source: str
    status: CutoffTreatmentStatus
    replacement_cutoff: SpatialCutoffOperator | None = None
    covariance_trace: CutoffDilationCovarianceTrace | None = None
    semantic_relation: ExactIdentity | FormalIdentity | ModCompactEquivalence | None = None
    assumptions: tuple[object, ...] = ()
    evidence: tuple[object, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.cutoff, SpatialCutoffOperator):
            raise TypeError("cutoff must be a SpatialCutoffOperator.")
        if self.cutoff_function != self.cutoff.cutoff_symbol:
            raise SchurCorrectionError("cutoff_function disagrees with the cutoff.")
        if not isinstance(self.branch, Branch):
            raise TypeError("branch must be a Branch.")
        if not isinstance(self.side, CutoffSide):
            raise TypeError("side must be a CutoffSide.")
        if not isinstance(self.status, CutoffTreatmentStatus):
            raise TypeError("status must be a CutoffTreatmentStatus.")
        if sp.simplify(self.scale - self.cutoff.radial_scale) != 0:
            raise SchurCorrectionError("cutoff scale disagrees with the typed cutoff.")
        object.__setattr__(self, "support_claim", _nonempty(self.support_claim, "support_claim"))
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(self, "assumptions", _items(self.assumptions, "assumptions"))
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))
        if self.status is CutoffTreatmentStatus.EXACT_COVARIANCE_ONLY and not isinstance(
            self.covariance_trace, CutoffDilationCovarianceTrace
        ):
            raise CutoffReplacementNotCertified(
                "EXACT_COVARIANCE_ONLY requires an exact covariance trace."
            )
        if self.status is CutoffTreatmentStatus.MOD_COMPACTS_REPLACEMENT_CERTIFIED:
            if not isinstance(self.semantic_relation, ModCompactEquivalence):
                raise CutoffReplacementNotCertified(
                    "a certified replacement requires a modulo-compact relation."
                )
            if not self.evidence or self.replacement_cutoff is None:
                raise CutoffReplacementNotCertified(
                    "a certified replacement requires evidence and a replacement."
                )
        if self.status is CutoffTreatmentStatus.CONDITIONAL_REPLACEMENT and (
            not isinstance(self.semantic_relation, FormalIdentity)
            or not self.assumptions
        ):
            raise CutoffReplacementNotCertified(
                "a conditional replacement requires assumptions and a formal relation."
            )
        if self.status in {
            CutoffTreatmentStatus.NO_REPLACEMENT_NEEDED,
            CutoffTreatmentStatus.REPLACEMENT_BLOCKED,
        } and self.semantic_relation is not None:
            raise CutoffReplacementNotCertified(
                "no-replacement and blocked records cannot carry a replacement relation."
            )


def audit_schur_cutoffs(
    assembly: SchurCorrectionAssemblyCertificate,
) -> tuple[CutoffAuditRecord, ...]:
    """Enumerate the four retained localization multipliers in source order."""

    if assembly.assembly_status not in {
        SchurAssemblyStatus.CERTIFIED_EXACT_SCHUR_ASSEMBLY,
        SchurAssemblyStatus.CERTIFIED_MOD_COMPACTS_SCHUR_ASSEMBLY,
        SchurAssemblyStatus.CONDITIONAL_SCHUR_ASSEMBLY,
    }:
        return ()
    if assembly.left_factor is None or assembly.right_factor is None:
        return ()
    records: list[CutoffAuditRecord] = []
    for exterior in (assembly.left_factor, assembly.right_factor):
        localized = exterior.localized_wiener_hopf
        for side, cutoff, branch in (
            (CutoffSide.LEFT, localized.left_cutoff, exterior.branch_to),
            (CutoffSide.RIGHT, localized.right_cutoff, exterior.branch_from),
        ):
            records.append(
                CutoffAuditRecord(
                    cutoff=cutoff,
                    cutoff_function=cutoff.cutoff_symbol,
                    argument=cutoff.radial_variable,
                    scale=cutoff.radial_scale,
                    branch=branch,
                    side=side,
                    adjacent_operator=localized.core,
                    support_claim=(
                        f"zero: {cutoff.equals_zero_away_from}; "
                        f"one: {cutoff.equals_one_near}; "
                        f"support: {cutoff.support_region}"
                    ),
                    source=cutoff.source,
                    status=CutoffTreatmentStatus.NO_REPLACEMENT_NEEDED,
                    evidence=(
                        localized.exact_definition,
                        "cutoff retained in the exact correction word",
                    ),
                )
            )
    return tuple(records)


def require_cutoff_replacement(
    record: CutoffAuditRecord,
) -> ModCompactEquivalence:
    """Return only a separately certified modulo-compact replacement."""

    if not isinstance(record, CutoffAuditRecord):
        raise TypeError("record must be a CutoffAuditRecord.")
    if (
        record.status
        is not CutoffTreatmentStatus.MOD_COMPACTS_REPLACEMENT_CERTIFIED
        or not isinstance(record.semantic_relation, ModCompactEquivalence)
    ):
        raise CutoffReplacementNotCertified(
            "exact covariance or support geometry is not a replacement theorem."
        )
    return record.semantic_relation


@dataclass(frozen=True)
class RelativeDilationAudit:
    first: WeightedDilationOperator
    second: WeightedDilationOperator
    are_inverse: bool
    adjacent: bool
    intervening_factors: tuple[object, ...]
    exact_covariance: ExactIdentity | None
    status: RelativeDilationStatus
    semantic_relation: ExactIdentity | FormalIdentity | None
    source: str
    assumptions: tuple[object, ...]
    evidence: tuple[object, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.first, WeightedDilationOperator) or not isinstance(
            self.second, WeightedDilationOperator
        ):
            raise TypeError("first and second must be WeightedDilationOperator objects.")
        if not isinstance(self.are_inverse, bool) or not isinstance(self.adjacent, bool):
            raise TypeError("are_inverse and adjacent must be bool values.")
        if not isinstance(self.status, RelativeDilationStatus):
            raise TypeError("status must be a RelativeDilationStatus.")
        intervening = _items(self.intervening_factors, "intervening_factors")
        if self.adjacent != (not intervening):
            raise SchurCorrectionError(
                "adjacency must agree with the intervening-factor ledger."
            )
        if self.status is RelativeDilationStatus.EXACT_ADJACENT_CANCELLATION:
            if not self.are_inverse or not self.adjacent or not isinstance(
                self.semantic_relation, ExactIdentity
            ):
                raise SchurCorrectionError(
                    "exact cancellation requires an adjacent typed inverse pair."
                )
        if self.status is RelativeDilationStatus.EXACT_COVARIANCE_REDUCTION:
            if not isinstance(self.exact_covariance, ExactIdentity) or not isinstance(
                self.semantic_relation, ExactIdentity
            ):
                raise SchurCorrectionError(
                    "exact covariance reduction requires a certified exact trace."
                )
        if self.status is RelativeDilationStatus.BLOCKED_BY_INTERVENING_FACTORS and (
            not self.are_inverse or not intervening or self.semantic_relation is not None
        ):
            raise SchurCorrectionError(
                "an intervening-factor blocker needs an inverse pair and no relation."
            )
        object.__setattr__(self, "intervening_factors", intervening)
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(self, "assumptions", _items(self.assumptions, "assumptions"))
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))


def _inverse_dilations(
    first: WeightedDilationOperator,
    second: WeightedDilationOperator,
) -> bool:
    return (
        first.domain == second.domain
        and first.convention is second.convention
        and sp.simplify(first.gamma * second.gamma - 1) == 0
        and sp.simplify(first.normalization_power - second.normalization_power) == 0
    )


def audit_relative_dilation_pair(
    first: WeightedDilationOperator,
    second: WeightedDilationOperator,
    *,
    intervening_factors: tuple[object, ...] = (),
    covariance_relation: ExactIdentity | None = None,
    assumptions: tuple[object, ...] = (),
    source: str,
    evidence: tuple[object, ...] = (),
) -> RelativeDilationAudit:
    """Audit one ordered pair without searching, commuting, or distributing."""

    if not isinstance(first, WeightedDilationOperator) or not isinstance(
        second, WeightedDilationOperator
    ):
        raise TypeError("first and second must be WeightedDilationOperator objects.")
    if first.codomain != second.domain:
        raise SchurDomainOrOrderError("the ordered dilation pair has incompatible spaces.")
    between = _items(intervening_factors, "intervening_factors")
    assumptions = _items(assumptions, "assumptions")
    evidence = _items(evidence, "evidence")
    inverse = _inverse_dilations(first, second)
    adjacent = not between
    relation: ExactIdentity | FormalIdentity | None = None
    if covariance_relation is not None:
        if not isinstance(covariance_relation, ExactIdentity):
            raise TypeError("covariance_relation must be an ExactIdentity.")
        status = RelativeDilationStatus.EXACT_COVARIANCE_REDUCTION
        relation = covariance_relation
    elif inverse and adjacent:
        status = RelativeDilationStatus.EXACT_ADJACENT_CANCELLATION
        relation = ExactIdentity(
            Product((first.atom, second.atom)),
            IDENTITY_PRODUCT,
            evidence=(source, evidence, first.evidence, second.evidence),
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=assumptions,
        )
    elif inverse:
        status = RelativeDilationStatus.BLOCKED_BY_INTERVENING_FACTORS
    else:
        status = RelativeDilationStatus.NO_RELATIVE_DILATION_PAIR
    return RelativeDilationAudit(
        first=first,
        second=second,
        are_inverse=inverse,
        adjacent=adjacent,
        intervening_factors=between,
        exact_covariance=covariance_relation,
        status=status,
        semantic_relation=relation,
        source=source,
        assumptions=assumptions,
        evidence=evidence,
    )


@dataclass(frozen=True)
class AlgebraFactorEvidence:
    factor: object
    factor_label: str
    candidate_algebra: str | None
    classification: FactorCoverageClassification
    source: str
    hypotheses: tuple[object, ...]
    evidence: tuple[object, ...]
    compact_ideal: CompactIdeal | None = None

    def __post_init__(self) -> None:
        _atom(self.factor, "factor")
        object.__setattr__(self, "factor_label", _nonempty(self.factor_label, "factor_label"))
        if self.candidate_algebra is not None:
            object.__setattr__(
                self,
                "candidate_algebra",
                _nonempty(self.candidate_algebra, "candidate_algebra"),
            )
        if not isinstance(self.classification, FactorCoverageClassification):
            raise TypeError("classification must be a FactorCoverageClassification.")
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(self, "hypotheses", _items(self.hypotheses, "hypotheses"))
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))
        if self.classification in {
            FactorCoverageClassification.EXACT_THEOREM,
            FactorCoverageClassification.DIRECT_ALGEBRAIC_CONSEQUENCE,
            FactorCoverageClassification.MATCH_AFTER_NOTATION_TRANSLATION,
        } and (self.candidate_algebra is None or not self.evidence):
            raise SchurCorrectionError(
                "positive factor coverage requires an algebra and evidence."
            )


_CERTIFIED_COVERAGE = {
    FactorCoverageClassification.EXACT_THEOREM,
    FactorCoverageClassification.DIRECT_ALGEBRAIC_CONSEQUENCE,
    FactorCoverageClassification.MATCH_AFTER_NOTATION_TRANSLATION,
}


@dataclass(frozen=True)
class CommonAlgebraCertificate:
    assembly: SchurCorrectionAssemblyCertificate
    factor_evidence: tuple[AlgebraFactorEvidence, ...]
    common_algebra: str | None
    compact_ideal: CompactIdeal | None
    closure_theorem: str | None
    inverse_closedness_theorem: str | None
    symbol_homomorphism: str | None
    explicit_extension: bool
    status: AlgebraMembershipStatus
    assumptions: AssumptionContext
    evidence: tuple[object, ...]
    missing_generator: str | None = None
    missing_closure_property: str | None = None
    missing_cutoff_theorem: str | None = None
    missing_dilation_membership: str | None = None
    missing_inverse_closedness: str | None = None
    missing_symbol_homomorphism: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.assembly, SchurCorrectionAssemblyCertificate):
            raise TypeError("assembly must be a SchurCorrectionAssemblyCertificate.")
        factors = _items(self.factor_evidence, "factor_evidence")
        if not all(isinstance(item, AlgebraFactorEvidence) for item in factors):
            raise TypeError("factor_evidence contains an invalid item.")
        if not isinstance(self.status, AlgebraMembershipStatus):
            raise TypeError("status must be an AlgebraMembershipStatus.")
        if not isinstance(self.explicit_extension, bool):
            raise TypeError("explicit_extension must be a bool.")
        _require_assumptions(self.assumptions)
        certified = self.status in {
            AlgebraMembershipStatus.CERTIFIED_COMMON_ALGEBRA_MEMBERSHIP,
            AlgebraMembershipStatus.CERTIFIED_MEMBERSHIP_IN_EXPLICIT_EXTENSION,
        }
        if certified:
            algebra = _nonempty(self.common_algebra, "common_algebra")
            if not isinstance(self.compact_ideal, CompactIdeal):
                raise SchurCorrectionError(
                    "certified common membership requires one typed compact ideal."
                )
            for field_name in (
                "closure_theorem",
                "inverse_closedness_theorem",
                "symbol_homomorphism",
            ):
                _nonempty(getattr(self, field_name), field_name)
            if not factors:
                raise SchurCorrectionError("certified membership needs factor evidence.")
            expected_atoms = (
                self.assembly.expanded_expression.factors
                if self.assembly.expanded_expression is not None
                else ()
            )
            covered_atoms = tuple(_atom(item.factor, "factor evidence") for item in factors)
            if covered_atoms != expected_atoms:
                raise SchurDomainOrOrderError(
                    "factor evidence must follow and cover the exact expanded word."
                )
            if any(
                item.classification not in _CERTIFIED_COVERAGE
                or item.candidate_algebra != algebra
                or item.compact_ideal != self.compact_ideal
                for item in factors
            ):
                raise SchurCorrectionError(
                    "certified membership cannot mix algebras, ideals, or source gaps."
                )
            if (
                self.status
                is AlgebraMembershipStatus.CERTIFIED_MEMBERSHIP_IN_EXPLICIT_EXTENSION
            ) != self.explicit_extension:
                raise SchurCorrectionError(
                    "explicit-extension status disagrees with its typed flag."
                )
            if any(
                value is not None
                for value in (
                    self.missing_generator,
                    self.missing_closure_property,
                    self.missing_cutoff_theorem,
                    self.missing_dilation_membership,
                    self.missing_inverse_closedness,
                    self.missing_symbol_homomorphism,
                )
            ):
                raise SchurCorrectionError(
                    "a certified algebra result cannot retain missing-theorem fields."
                )
        elif self.status is AlgebraMembershipStatus.CONDITIONAL_ALGEBRA_MEMBERSHIP:
            if self.assumptions.is_empty:
                raise SchurCorrectionError(
                    "conditional algebra membership requires assumptions."
                )
        elif self.status is AlgebraMembershipStatus.BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP:
            if not any(
                (
                    self.missing_generator,
                    self.missing_closure_property,
                    self.missing_cutoff_theorem,
                    self.missing_dilation_membership,
                    self.missing_inverse_closedness,
                    self.missing_symbol_homomorphism,
                )
            ):
                raise SchurCorrectionError(
                    "the common-source-gap status must identify a missing theorem field."
                )
        object.__setattr__(self, "factor_evidence", factors)
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))

    @property
    def is_certified(self) -> bool:
        return self.status in {
            AlgebraMembershipStatus.CERTIFIED_COMMON_ALGEBRA_MEMBERSHIP,
            AlgebraMembershipStatus.CERTIFIED_MEMBERSHIP_IN_EXPLICIT_EXTENSION,
        }


def article_factor_membership_evidence(
    assembly: SchurCorrectionAssemblyCertificate,
) -> tuple[AlgebraFactorEvidence, ...]:
    """Return the audited factor matrix for the article instance."""

    if assembly.expanded_expression is None or assembly.central_regularizer is None:
        return ()
    factors = assembly.expanded_expression.factors
    if len(factors) != 7:
        return ()
    ideal = assembly.compact_ideal
    data = (
        (
            "Vtilde_alpha2-G2",
            "project exact factor; no common dilation-stable algebra theorem",
            None,
            FactorCoverageClassification.SOURCE_GAP,
            (),
        ),
        (
            "W21_minus",
            "Karlovich 2025 half-line Wiener--Hopf algebra identification",
            "Karlovich-2025 half-line algebra",
            FactorCoverageClassification.EXACT_THEOREM,
            ("both localization cutoffs retained",),
        ),
        (
            "T1_minus",
            "article exact invertibility; common-algebra membership absent",
            None,
            FactorCoverageClassification.SOURCE_GAP,
            (),
        ),
        (
            "R1",
            "KKL 2014 radial Mellin-PDO calculus",
            "KKL-2014 radial Mellin PDO algebra",
            FactorCoverageClassification.MATCH_AFTER_NOTATION_TRANSLATION,
            ("transport by Phi_delta retained",),
        ),
        (
            "Z1_inverse",
            "bounded invertible coefficient multiplier",
            "Karlovich-2025 half-line algebra",
            FactorCoverageClassification.DIRECT_ALGEBRAIC_CONSEQUENCE,
            ("nonvanishing multiplier evidence from P1-C",),
        ),
        (
            "Vtilde_alpha1-G1",
            "project exact factor; no common dilation-stable algebra theorem",
            None,
            FactorCoverageClassification.SOURCE_GAP,
            (),
        ),
        (
            "W12_plus",
            "Karlovich 2025 half-line Wiener--Hopf algebra identification",
            "Karlovich-2025 half-line algebra",
            FactorCoverageClassification.EXACT_THEOREM,
            ("both localization cutoffs retained",),
        ),
    )
    return tuple(
        AlgebraFactorEvidence(
            factor=factor,
            factor_label=label,
            candidate_algebra=algebra,
            classification=classification,
            source=source,
            hypotheses=assembly.assumptions.assumptions,
            evidence=evidence,
            compact_ideal=ideal,
        )
        for factor, (label, source, algebra, classification, evidence) in zip(
            factors, data, strict=True
        )
    )


def audit_common_algebra_membership(
    *,
    assembly: SchurCorrectionAssemblyCertificate,
    factor_evidence: tuple[AlgebraFactorEvidence, ...],
    assumptions: AssumptionContext,
    common_algebra: str | None = None,
    compact_ideal: CompactIdeal | None = None,
    closure_theorem: str | None = None,
    inverse_closedness_theorem: str | None = None,
    symbol_homomorphism: str | None = None,
    explicit_extension: bool = False,
    cutoff_audits: tuple[CutoffAuditRecord, ...] = (),
    dilation_audits: tuple[RelativeDilationAudit, ...] = (),
    requires_dilation_reduction: bool = False,
    evidence: tuple[object, ...] = (),
) -> CommonAlgebraCertificate:
    """Certify a common algebra only from one complete, sourced evidence tuple."""

    assumptions = _require_assumptions(assumptions)
    factors = _items(factor_evidence, "factor_evidence")
    cutoffs = _items(cutoff_audits, "cutoff_audits")
    dilations = _items(dilation_audits, "dilation_audits")
    supplied_evidence = _items(evidence, "evidence")
    successful_assembly = assembly.assembly_status in {
        SchurAssemblyStatus.CERTIFIED_EXACT_SCHUR_ASSEMBLY,
        SchurAssemblyStatus.CERTIFIED_MOD_COMPACTS_SCHUR_ASSEMBLY,
        SchurAssemblyStatus.CONDITIONAL_SCHUR_ASSEMBLY,
    }
    status = AlgebraMembershipStatus.BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP
    missing = {
        "missing_generator": None,
        "missing_closure_property": None,
        "missing_cutoff_theorem": None,
        "missing_dilation_membership": None,
        "missing_inverse_closedness": None,
        "missing_symbol_homomorphism": None,
    }
    if not successful_assembly:
        status = AlgebraMembershipStatus.NOT_ATTEMPTED_BECAUSE_ASSEMBLY_FAILED
    elif any(
        isinstance(item, CutoffAuditRecord)
        and item.status
        in {
            CutoffTreatmentStatus.CONDITIONAL_REPLACEMENT,
            CutoffTreatmentStatus.REPLACEMENT_BLOCKED,
        }
        for item in cutoffs
    ):
        status = AlgebraMembershipStatus.BLOCKED_BY_CUTOFF_REPLACEMENT
        missing["missing_cutoff_theorem"] = (
            "the requested scaled-cutoff replacement theorem"
        )
    elif requires_dilation_reduction and any(
        isinstance(item, RelativeDilationAudit)
        and item.status
        in {
            RelativeDilationStatus.CONDITIONAL_REDUCTION,
            RelativeDilationStatus.BLOCKED_BY_INTERVENING_FACTORS,
        }
        for item in dilations
    ):
        status = AlgebraMembershipStatus.BLOCKED_BY_RELATIVE_DILATION
        missing["missing_dilation_membership"] = (
            "a valid reduction or algebra membership for the nonadjacent dilation"
        )
    else:
        exact_word = (
            assembly.expanded_expression.factors
            if assembly.expanded_expression is not None
            else ()
        )
        evidence_atoms = tuple(
            _atom(item.factor, "factor evidence")
            for item in factors
            if isinstance(item, AlgebraFactorEvidence)
        )
        common = common_algebra.strip() if isinstance(common_algebra, str) else None
        certifiable = (
            bool(common)
            and evidence_atoms == exact_word
            and all(
                isinstance(item, AlgebraFactorEvidence)
                and item.classification in _CERTIFIED_COVERAGE
                and item.candidate_algebra == common
                and item.compact_ideal == compact_ideal
                for item in factors
            )
            and isinstance(compact_ideal, CompactIdeal)
            and bool(closure_theorem)
            and bool(inverse_closedness_theorem)
            and bool(symbol_homomorphism)
            and bool(supplied_evidence)
        )
        if certifiable:
            status = (
                AlgebraMembershipStatus.CERTIFIED_MEMBERSHIP_IN_EXPLICIT_EXTENSION
                if explicit_extension
                else AlgebraMembershipStatus.CERTIFIED_COMMON_ALGEBRA_MEMBERSHIP
            )
        elif (
            common
            and factors
            and assumptions.is_empty is False
            and all(
                isinstance(item, AlgebraFactorEvidence)
                and item.classification
                in _CERTIFIED_COVERAGE
                | {FactorCoverageClassification.PLAUSIBLE_ADAPTATION}
                for item in factors
            )
        ):
            status = AlgebraMembershipStatus.CONDITIONAL_ALGEBRA_MEMBERSHIP
        else:
            wh_positions = {1, 6}
            if any(
                index in wh_positions
                and isinstance(item, AlgebraFactorEvidence)
                and item.classification
                in {
                    FactorCoverageClassification.NOT_COVERED,
                    FactorCoverageClassification.CONTRADICTED_BY_HYPOTHESES,
                }
                for index, item in enumerate(factors)
            ):
                status = AlgebraMembershipStatus.BLOCKED_BY_OFF_DIAGONAL_WH_MEMBERSHIP
            elif len(factors) > 3 and isinstance(factors[3], AlgebraFactorEvidence) and (
                factors[3].classification
                in {
                    FactorCoverageClassification.NOT_COVERED,
                    FactorCoverageClassification.CONTRADICTED_BY_HYPOTHESES,
                }
            ):
                status = (
                    AlgebraMembershipStatus.BLOCKED_BY_CENTRAL_REGULARIZER_MEMBERSHIP
                )
            else:
                status = AlgebraMembershipStatus.BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP
                missing.update(
                    {
                        "missing_generator": (
                            "one algebra containing the transported radial Mellin PDO "
                            "and nontrivial dilation interfaces"
                        ),
                        "missing_closure_property": (
                            "closure under the exact seven-factor ordered product"
                        ),
                        "missing_cutoff_theorem": (
                            "one theorem coordinating all retained double cutoffs"
                        ),
                        "missing_dilation_membership": (
                            "membership of V_gamma and T1_minus in the common algebra"
                        ),
                        "missing_inverse_closedness": (
                            "inverse closedness for the actual mixed extension"
                        ),
                        "missing_symbol_homomorphism": (
                            "a multiplicative Fredholm symbol on that mixed algebra"
                        ),
                    }
                )
    return CommonAlgebraCertificate(
        assembly=assembly,
        factor_evidence=tuple(
            item for item in factors if isinstance(item, AlgebraFactorEvidence)
        ),
        common_algebra=common_algebra,
        compact_ideal=compact_ideal,
        closure_theorem=closure_theorem,
        inverse_closedness_theorem=inverse_closedness_theorem,
        symbol_homomorphism=symbol_homomorphism,
        explicit_extension=explicit_extension,
        status=status,
        assumptions=assumptions,
        evidence=supplied_evidence,
        **missing,
    )


@dataclass(frozen=True)
class FactorSymbol:
    factor: OperatorAtom
    symbol: object
    algebra: str
    source: str
    evidence: tuple[object, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.factor, OperatorAtom):
            raise TypeError("factor must be an OperatorAtom.")
        if self.symbol is None:
            raise ValueError("symbol must be supplied.")
        object.__setattr__(self, "algebra", _nonempty(self.algebra, "algebra"))
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise SchurCorrectionError("a factor symbol requires evidence.")
        object.__setattr__(self, "evidence", evidence)


@dataclass(frozen=True)
class NoncommutativeProductRuleCertificate:
    algebra: str
    ordered: bool
    multiplicative: bool
    source: str
    evidence: tuple[object, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "algebra", _nonempty(self.algebra, "algebra"))
        if not isinstance(self.ordered, bool) or not isinstance(
            self.multiplicative, bool
        ):
            raise TypeError("ordered and multiplicative must be bool values.")
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise SchurCorrectionError("a product rule requires evidence.")
        object.__setattr__(self, "evidence", evidence)


@dataclass(frozen=True)
class CertifiedSchurCorrectionSymbol:
    correction: OperatorAtom
    algebra: str
    ordered_factor_symbols: tuple[FactorSymbol, ...]
    product_rule: NoncommutativeProductRuleCertificate
    expression: tuple[object, ...] = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.correction, OperatorAtom):
            raise TypeError("correction must be an OperatorAtom.")
        object.__setattr__(self, "algebra", _nonempty(self.algebra, "algebra"))
        symbols = _items(self.ordered_factor_symbols, "ordered_factor_symbols")
        if not symbols or not all(isinstance(item, FactorSymbol) for item in symbols):
            raise TypeError("ordered_factor_symbols must contain FactorSymbol objects.")
        if not isinstance(self.product_rule, NoncommutativeProductRuleCertificate):
            raise TypeError("product_rule has an invalid type.")
        if self.product_rule.algebra != self.algebra or any(
            item.algebra != self.algebra for item in symbols
        ):
            raise SchurCorrectionError("all symbol data must use one algebra.")
        if not self.product_rule.ordered or not self.product_rule.multiplicative:
            raise NoncommutativeProductRuleMissing(
                "the supplied rule does not certify ordered multiplication."
            )
        object.__setattr__(self, "ordered_factor_symbols", symbols)
        object.__setattr__(self, "expression", tuple(item.symbol for item in symbols))


@dataclass(frozen=True)
class SchurCorrectionSymbolResult:
    assembly: SchurCorrectionAssemblyCertificate
    algebra_certificate: CommonAlgebraCertificate
    symbol_status: SchurCorrectionSymbolStatus
    symbol: CertifiedSchurCorrectionSymbol | None
    blockers: tuple[object, ...]
    evidence: tuple[object, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.symbol_status, SchurCorrectionSymbolStatus):
            raise TypeError("symbol_status must be a SchurCorrectionSymbolStatus.")
        certified = self.symbol_status is (
            SchurCorrectionSymbolStatus.CERTIFIED_SCHUR_CORRECTION_SYMBOL
        )
        if certified != isinstance(self.symbol, CertifiedSchurCorrectionSymbol):
            raise SchurCorrectionError(
                "only a certified result may carry a certified symbol."
            )
        object.__setattr__(self, "blockers", _items(self.blockers, "blockers"))
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))


def derive_schur_correction_symbol(
    *,
    assembly: SchurCorrectionAssemblyCertificate,
    algebra_certificate: CommonAlgebraCertificate,
    factor_symbols: tuple[FactorSymbol, ...] = (),
    product_rule: NoncommutativeProductRuleCertificate | None = None,
    cutoff_audits: tuple[CutoffAuditRecord, ...] = (),
    dilation_audits: tuple[RelativeDilationAudit, ...] = (),
) -> SchurCorrectionSymbolResult:
    """Apply the strict symbol gate without producing a fallback formula."""

    if assembly.assembly_status not in {
        SchurAssemblyStatus.CERTIFIED_EXACT_SCHUR_ASSEMBLY,
        SchurAssemblyStatus.CERTIFIED_MOD_COMPACTS_SCHUR_ASSEMBLY,
        SchurAssemblyStatus.CONDITIONAL_SCHUR_ASSEMBLY,
    }:
        return SchurCorrectionSymbolResult(
            assembly,
            algebra_certificate,
            SchurCorrectionSymbolStatus.NOT_ATTEMPTED,
            None,
            (assembly,),
            (),
        )
    if not algebra_certificate.is_certified and algebra_certificate.status is not (
        AlgebraMembershipStatus.CONDITIONAL_ALGEBRA_MEMBERSHIP
    ):
        return SchurCorrectionSymbolResult(
            assembly,
            algebra_certificate,
            SchurCorrectionSymbolStatus.BLOCKED_BY_ALGEBRA_MEMBERSHIP,
            None,
            (algebra_certificate,),
            (),
        )
    cutoffs = _items(cutoff_audits, "cutoff_audits")
    if len(cutoffs) != 4 or any(
        not isinstance(item, CutoffAuditRecord)
        or item.status
        not in {
            CutoffTreatmentStatus.NO_REPLACEMENT_NEEDED,
            CutoffTreatmentStatus.MOD_COMPACTS_REPLACEMENT_CERTIFIED,
        }
        for item in cutoffs
    ):
        return SchurCorrectionSymbolResult(
            assembly,
            algebra_certificate,
            SchurCorrectionSymbolStatus.BLOCKED_BY_CUTOFFS,
            None,
            cutoffs,
            (),
        )
    dilations = _items(dilation_audits, "dilation_audits")
    if not dilations or any(
        not isinstance(item, RelativeDilationAudit)
        or item.status
        in {
            RelativeDilationStatus.CONDITIONAL_REDUCTION,
            RelativeDilationStatus.BLOCKED_BY_INTERVENING_FACTORS,
        }
        for item in dilations
    ):
        return SchurCorrectionSymbolResult(
            assembly,
            algebra_certificate,
            SchurCorrectionSymbolStatus.BLOCKED_BY_NONCOMMUTATIVE_PRODUCT,
            None,
            dilations,
            (),
        )
    if product_rule is None or not product_rule.ordered or not product_rule.multiplicative:
        return SchurCorrectionSymbolResult(
            assembly,
            algebra_certificate,
            SchurCorrectionSymbolStatus.BLOCKED_BY_NONCOMMUTATIVE_PRODUCT,
            None,
            ((product_rule,) if product_rule is not None else ()),
            (),
        )
    symbols = _items(factor_symbols, "factor_symbols")
    expected = assembly.expanded_expression.factors if assembly.expanded_expression else ()
    algebra = algebra_certificate.common_algebra
    if (
        not isinstance(algebra, str)
        or tuple(
            item.factor for item in symbols if isinstance(item, FactorSymbol)
        )
        != expected
        or not all(
            isinstance(item, FactorSymbol) and item.algebra == algebra
            for item in symbols
        )
        or product_rule.algebra != algebra
    ):
        return SchurCorrectionSymbolResult(
            assembly,
            algebra_certificate,
            SchurCorrectionSymbolStatus.BLOCKED_BY_NONCOMMUTATIVE_PRODUCT,
            None,
            (symbols, product_rule),
            (),
        )
    certified_symbol = CertifiedSchurCorrectionSymbol(
        correction=assembly.correction_atom,
        algebra=algebra,
        ordered_factor_symbols=tuple(symbols),
        product_rule=product_rule,
    )
    status = (
        SchurCorrectionSymbolStatus.CERTIFIED_SCHUR_CORRECTION_SYMBOL
        if algebra_certificate.is_certified
        and assembly.assembly_status
        is not SchurAssemblyStatus.CONDITIONAL_SCHUR_ASSEMBLY
        else SchurCorrectionSymbolStatus.CONDITIONAL_SCHUR_CORRECTION_SYMBOL
    )
    if status is SchurCorrectionSymbolStatus.CONDITIONAL_SCHUR_CORRECTION_SYMBOL:
        return SchurCorrectionSymbolResult(
            assembly,
            algebra_certificate,
            status,
            None,
            ("conditional inputs cannot produce a certified symbol object",),
            (symbols, product_rule),
        )
    return SchurCorrectionSymbolResult(
        assembly,
        algebra_certificate,
        status,
        certified_symbol,
        (),
        (symbols, product_rule, cutoffs, dilations),
    )


def require_certified_schur_symbol(
    result: SchurCorrectionSymbolResult,
) -> CertifiedSchurCorrectionSymbol:
    if not isinstance(result, SchurCorrectionSymbolResult):
        raise TypeError("result must be a SchurCorrectionSymbolResult.")
    if result.symbol is not None:
        return result.symbol
    if result.symbol_status is SchurCorrectionSymbolStatus.BLOCKED_BY_ALGEBRA_MEMBERSHIP:
        raise AlgebraMembershipNotCertified("common algebra membership is not certified.")
    if result.symbol_status is SchurCorrectionSymbolStatus.BLOCKED_BY_CUTOFFS:
        raise CutoffTreatmentNotCertified("cutoff treatment is not certified.")
    if result.symbol_status is SchurCorrectionSymbolStatus.BLOCKED_BY_NONCOMMUTATIVE_PRODUCT:
        if any(isinstance(item, RelativeDilationAudit) for item in result.blockers):
            raise RelativeDilationNotResolved("a relative dilation remains unresolved.")
        raise NoncommutativeProductRuleMissing(
            "the ordered noncommutative product rule is incomplete."
        )
    raise SchurSymbolNotAvailable("the Schur correction symbol is not available.")


@dataclass(frozen=True)
class BranchedDiagonalOperator:
    atom: OperatorAtom
    branch: Branch
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    source: str
    evidence: tuple[object, ...]
    bounded: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.branch, Branch):
            raise TypeError("branch must be a Branch.")
        _operator_space(self, "diagonal block")
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise SchurCorrectionError("a diagonal block requires evidence.")
        object.__setattr__(self, "evidence", evidence)
        if not isinstance(self.bounded, bool):
            raise TypeError("bounded must be a bool.")


@dataclass(frozen=True)
class OffDiagonalBlockModel:
    block_atom: OperatorAtom
    branch_from: Branch
    branch_to: Branch
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    normalized_exterior: SchurExteriorFactor
    model_sign: int
    compact_ideal: CompactIdeal
    source: str
    evidence: tuple[object, ...]
    semantic_relation: ModCompactEquivalence = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.block_atom, OperatorAtom):
            raise TypeError("block_atom must be an OperatorAtom.")
        if not isinstance(self.branch_from, Branch) or not isinstance(
            self.branch_to, Branch
        ):
            raise TypeError("branch_from and branch_to must be Branch objects.")
        if not isinstance(self.normalized_exterior, SchurExteriorFactor):
            raise TypeError("normalized_exterior must be a SchurExteriorFactor.")
        if isinstance(self.model_sign, bool) or self.model_sign not in {-1, 1}:
            raise SchurCorrectionError("model_sign must be exactly -1 or +1.")
        if (
            self.branch_from != self.normalized_exterior.branch_from
            or self.branch_to != self.normalized_exterior.branch_to
            or self.domain != self.normalized_exterior.domain
            or self.codomain != self.normalized_exterior.codomain
        ):
            raise SchurDomainOrOrderError(
                "off-diagonal model orientation disagrees with its exterior."
            )
        if not isinstance(self.compact_ideal, CompactIdeal):
            raise TypeError("compact_ideal must be a CompactIdeal.")
        if self.compact_ideal.space != self.domain or not self.compact_ideal.bilateral:
            raise SchurDomainOrOrderError(
                "off-diagonal relation requires the same bilateral compact ideal."
            )
        source = _nonempty(self.source, "source")
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise SchurCorrectionError("off-diagonal compactness requires evidence.")
        model = LinearCombination(
            (Term(self.model_sign, self.normalized_exterior.atom),)
        )
        relation = ModCompactEquivalence(
            self.block_atom,
            model,
            space=self.domain,
            compact_ideal=self.compact_ideal,
            evidence=(source, evidence),
        )
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "semantic_relation", relation)


@dataclass(frozen=True)
class SchurSignConvention:
    raw_schur_sign: int
    left_model_sign: int
    right_model_sign: int
    correction_sign_in_model_pivot: int
    source: str
    evidence: tuple[object, ...]

    def __post_init__(self) -> None:
        signs = (
            self.raw_schur_sign,
            self.left_model_sign,
            self.right_model_sign,
            self.correction_sign_in_model_pivot,
        )
        if any(isinstance(value, bool) or value not in {-1, 1} for value in signs):
            raise SchurCorrectionError("all Schur signs must be exactly -1 or +1.")
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        evidence = _items(self.evidence, "evidence")
        if not evidence:
            raise SchurCorrectionError("the sign convention requires source evidence.")
        object.__setattr__(self, "evidence", evidence)

    @property
    def is_consistent(self) -> bool:
        return (
            self.raw_schur_sign * self.left_model_sign * self.right_model_sign
            == self.correction_sign_in_model_pivot
        )


@dataclass(frozen=True)
class FirstSchurPivotCertificate:
    diagonal_block: BranchedDiagonalOperator
    correction: SchurCorrectionAssemblyCertificate
    sign: SchurSignConvention
    pivot_expression: LinearCombination | None
    source_expression: LinearCombination | None
    raw_schur_word: Product | None
    domain: BranchedWeightedSpace
    codomain: BranchedWeightedSpace
    assembly_status: SchurAssemblyStatus
    pivot_relation: ModCompactEquivalence | None
    fredholm_status: PivotFredholmStatus
    symbol_status: SchurCorrectionSymbolStatus
    proof_obligations: tuple[AnalyticProofObligation, ...]
    latex: str

    def __post_init__(self) -> None:
        if not isinstance(self.diagonal_block, BranchedDiagonalOperator):
            raise TypeError("diagonal_block must be a BranchedDiagonalOperator.")
        if not isinstance(self.correction, SchurCorrectionAssemblyCertificate):
            raise TypeError("correction must be a SchurCorrectionAssemblyCertificate.")
        if not isinstance(self.sign, SchurSignConvention):
            raise TypeError("sign must be a SchurSignConvention.")
        if not isinstance(self.assembly_status, SchurAssemblyStatus):
            raise TypeError("assembly_status must be a SchurAssemblyStatus.")
        if not isinstance(self.fredholm_status, PivotFredholmStatus):
            raise TypeError("fredholm_status must be a PivotFredholmStatus.")
        if not isinstance(self.symbol_status, SchurCorrectionSymbolStatus):
            raise TypeError("symbol_status must be a SchurCorrectionSymbolStatus.")
        if self.assembly_status is SchurAssemblyStatus.BLOCKED_BY_SIGN_CONVENTION:
            if self.pivot_expression is not None or self.pivot_relation is not None:
                raise SchurCorrectionError(
                    "a sign-blocked pivot cannot carry a model expression or relation."
                )
        elif not isinstance(self.pivot_relation, ModCompactEquivalence):
            raise SchurCorrectionError(
                "a successful pivot needs the source/model modulo-compact relation."
            )
        object.__setattr__(self, "latex", _nonempty(self.latex, "latex"))


def certify_first_schur_pivot(
    *,
    diagonal_block: BranchedDiagonalOperator,
    correction: SchurCorrectionAssemblyCertificate,
    left_block: OffDiagonalBlockModel,
    right_block: OffDiagonalBlockModel,
    sign: SchurSignConvention,
    symbol_result: SchurCorrectionSymbolResult | None = None,
) -> FirstSchurPivotCertificate:
    """Translate the raw Schur minus into the positive normalized correction."""

    if not isinstance(diagonal_block, BranchedDiagonalOperator):
        raise TypeError("diagonal_block must be a BranchedDiagonalOperator.")
    if not isinstance(left_block, OffDiagonalBlockModel) or not isinstance(
        right_block, OffDiagonalBlockModel
    ):
        raise TypeError("left_block and right_block must be OffDiagonalBlockModel objects.")
    if correction.central_regularizer is None:
        raise SchurCorrectionError("the correction has no central regularizer.")
    ideal = correction.compact_ideal
    if ideal is None:
        raise SchurCorrectionError("the correction has no typed compact ideal.")
    if (
        diagonal_block.branch != Branch(2)
        or diagonal_block.domain != ideal.space
        or left_block.compact_ideal != ideal
        or right_block.compact_ideal != ideal
        or left_block.normalized_exterior != correction.left_factor
        or right_block.normalized_exterior != correction.right_factor
    ):
        raise SchurDomainOrOrderError(
            "pivot blocks, correction, branches, or compact ideals disagree."
        )
    symbol_status = (
        symbol_result.symbol_status
        if symbol_result is not None
        else SchurCorrectionSymbolStatus.NOT_ATTEMPTED
    )
    raw_word = Product(
        (
            left_block.block_atom,
            correction.central_regularizer.candidate_regularizer.atom,
            right_block.block_atom,
        )
    )
    raw_pivot = LinearCombination(
        (Term(1, diagonal_block.atom), Term(sign.raw_schur_sign, raw_word))
    )
    sign_ok = (
        sign.is_consistent
        and sign.raw_schur_sign == -1
        and sign.left_model_sign == left_block.model_sign == -1
        and sign.right_model_sign == right_block.model_sign == 1
        and sign.correction_sign_in_model_pivot == 1
    )
    domain = BranchedWeightedSpace(diagonal_block.branch, diagonal_block.domain)
    if not sign_ok:
        return FirstSchurPivotCertificate(
            diagonal_block=diagonal_block,
            correction=correction,
            sign=sign,
            pivot_expression=None,
            source_expression=raw_pivot,
            raw_schur_word=raw_word,
            domain=domain,
            codomain=domain,
            assembly_status=SchurAssemblyStatus.BLOCKED_BY_SIGN_CONVENTION,
            pivot_relation=None,
            fredholm_status=PivotFredholmStatus.NOT_CERTIFIED,
            symbol_status=symbol_status,
            proof_obligations=correction.proof_obligations,
            latex=r"\mathcal A_{2,2}^{(1)}\quad\text{BLOCKED_BY_SIGN_CONVENTION}",
        )
    model_pivot = LinearCombination(
        (
            Term(1, diagonal_block.atom),
            Term(sign.correction_sign_in_model_pivot, correction.correction_atom),
        )
    )
    relation = ModCompactEquivalence(
        raw_pivot,
        model_pivot,
        space=diagonal_block.domain,
        compact_ideal=ideal,
        evidence=(
            left_block.semantic_relation,
            right_block.semantic_relation,
            correction.assembly_relation,
            sign.evidence,
            "bilateral compact-ideal propagation through bounded factors",
        ),
    )
    return FirstSchurPivotCertificate(
        diagonal_block=diagonal_block,
        correction=correction,
        sign=sign,
        pivot_expression=model_pivot,
        source_expression=raw_pivot,
        raw_schur_word=raw_word,
        domain=domain,
        codomain=domain,
        assembly_status=correction.assembly_status,
        pivot_relation=relation,
        fredholm_status=PivotFredholmStatus.NOT_CERTIFIED,
        symbol_status=symbol_status,
        proof_obligations=correction.proof_obligations,
        latex=(
            r"\mathcal A_{2,2}^{(1)}="
            r"\mathcal A_{2,2}-\mathcal A_{2,1}"
            r"\mathcal A_{1,1}^{(-1)}\mathcal A_{1,2}"
            r"\simeq\mathcal A_{2,2}+\mathcal C_{2,2}^{(1)}"
        ),
    )


__all__ = [
    "AlgebraFactorEvidence",
    "AlgebraMembershipNotCertified",
    "AlgebraMembershipStatus",
    "AssemblyTraceStep",
    "Branch",
    "BranchDifferenceOperator",
    "BranchedDiagonalOperator",
    "BranchedWeightedSpace",
    "CertifiedSchurCorrectionSymbol",
    "CommonAlgebraCertificate",
    "CutoffAuditRecord",
    "CutoffReplacementNotCertified",
    "CutoffSide",
    "CutoffTreatmentNotCertified",
    "CutoffTreatmentStatus",
    "ExteriorFactorStatus",
    "FactorCoverageClassification",
    "FactorSymbol",
    "FirstSchurPivotCertificate",
    "NoncommutativeProductRuleMissing",
    "NoncommutativeProductRuleCertificate",
    "OffDiagonalBlockModel",
    "PivotFredholmStatus",
    "RelativeDilationAudit",
    "RelativeDilationNotResolved",
    "RelativeDilationStatus",
    "SchurAssemblyStatus",
    "SchurCorrectionAssemblyCertificate",
    "SchurCorrectionError",
    "SchurCorrectionSymbolResult",
    "SchurCorrectionSymbolStatus",
    "SchurCorrectionWordLedger",
    "SchurDomainOrOrderError",
    "SchurExteriorFactor",
    "SchurSignConvention",
    "SchurSourceBlock",
    "SchurSymbolNotAvailable",
    "SourceBlockSide",
    "article_factor_membership_evidence",
    "assemble_first_schur_correction",
    "audit_common_algebra_membership",
    "audit_relative_dilation_pair",
    "audit_schur_cutoffs",
    "certify_first_schur_pivot",
    "derive_schur_correction_symbol",
    "require_certified_schur_symbol",
    "require_cutoff_replacement",
]
