"""Safe relative-dilation cancellation for a narrowly typed Schur model.

The public rule in this module applies only after callers explicitly supply
the two off-diagonal block representations and a separate Mellin
representation of the abstract regularizer.  It intentionally does not match
or simplify the repository's current first-pivot builders, whose actual
factor order contains additional Cauchy, coefficient, and localization
operators.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import sympy as sp

from .domains import AssumptionContext, ConsistencyStatus, IncompatibleDomainError
from .mellin.domains import MellinDomainError
from .mellin.operators import (
    AlgebraMembershipClaim,
    AlgebraMembershipStatus,
    MellinOperatorError,
    MellinPseudodifferentialOperator,
    RegularizerMellinRepresentation,
    RegularizerRepresentationStatus,
    WeightedDilationOperator,
    WeightedLpSpace,
    WienerHopfOperator,
    conjugate_mellin_pdo_by_dilation,
)
from .operators import OperatorAtom, Product
from .relations import ExactBlock
from .semantics import (
    AnalyticProofObligation,
    DerivationRelationKind,
    ExactIdentity,
    ExactIdentityScope,
    FormalIdentity,
    ModCompactEquivalence,
    RegularizerOperator,
)


class SchurDilationError(ValueError):
    """Raised when a Schur representation object is internally malformed."""


class CancellationLayerStatus(str, Enum):
    """Orthogonal logical states required by the phase audit."""

    EXACT_IDENTITY = "EXACT_IDENTITY"
    EQUALITY_MODULO_COMPACTS = "EQUALITY_MODULO_COMPACTS"
    CONDITIONAL_ON_REGULARIZER_REPRESENTATION = (
        "CONDITIONAL_ON_REGULARIZER_REPRESENTATION"
    )
    UNPROVED_MIXED_ALGEBRA_MEMBERSHIP = (
        "UNPROVED_MIXED_ALGEBRA_MEMBERSHIP"
    )
    NOT_APPLIED = "NOT_APPLIED"


class SchurCancellationStatus(str, Enum):
    """Top-level status; the four layers remain separately inspectable."""

    EXACT_IDENTITY = "EXACT_IDENTITY"
    EQUALITY_MODULO_COMPACTS = "EQUALITY_MODULO_COMPACTS"
    CONDITIONAL_ON_REGULARIZER_REPRESENTATION = (
        "CONDITIONAL_ON_REGULARIZER_REPRESENTATION"
    )
    BLOCKED_WITH_PRECISE_OBLIGATION = "BLOCKED_WITH_PRECISE_OBLIGATION"


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


def _factor_atom(factor: object) -> OperatorAtom:
    if isinstance(factor, OperatorAtom):
        return factor
    atom = getattr(factor, "atom", None)
    if not isinstance(atom, OperatorAtom):
        raise TypeError(
            "each block-model factor must be an OperatorAtom or expose one as .atom."
        )
    return atom


@dataclass(frozen=True)
class BlockRepresentationEvidence:
    source: str
    supplied_evidence: tuple[object, ...]


@dataclass(frozen=True)
class SchurBlockRepresentation:
    """One exact or modulo-compact model of an off-diagonal block."""

    exact_block: ExactBlock
    atom: OperatorAtom
    factors: tuple[object, ...]
    relation_kind: DerivationRelationKind
    source: str
    hypotheses: tuple[object, ...]
    evidence: tuple[object, ...]
    bounded: bool = True
    localization_controlled: bool = True
    compact_ideal: object = "compact operators"
    model: Product = field(init=False)
    semantic_relation: ExactIdentity | ModCompactEquivalence = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.exact_block, ExactBlock):
            raise TypeError("exact_block must be an ExactBlock.")
        if self.exact_block.row == self.exact_block.column:
            raise SchurDilationError("a Schur off-diagonal representation needs row != column.")
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        factors = _items(self.factors, "factors")
        if not factors:
            raise SchurDilationError("a block representation needs a nonempty model.")
        model = Product(tuple(_factor_atom(factor) for factor in factors))
        if self.relation_kind not in {
            DerivationRelationKind.EXACT_EQUALITY,
            DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
        }:
            raise SchurDilationError(
                "a block representation must be exact or modulo compact operators."
            )
        source = _nonempty(self.source, "source")
        hypotheses = _items(self.hypotheses, "hypotheses")
        evidence = _items(self.evidence, "evidence")
        if not isinstance(self.bounded, bool) or not isinstance(
            self.localization_controlled, bool
        ):
            raise TypeError("bounded and localization_controlled must be bool values.")
        relation: ExactIdentity | ModCompactEquivalence
        if self.relation_kind is DerivationRelationKind.EXACT_EQUALITY:
            relation = ExactIdentity(
                self.atom,
                model,
                evidence=BlockRepresentationEvidence(source, evidence),
                scope=ExactIdentityScope.OPERATORIAL,
                hypotheses=hypotheses,
            )
        else:
            relation = ModCompactEquivalence(
                self.atom,
                model,
                compact_ideal=self.compact_ideal,
                evidence=evidence or None,
            )
        object.__setattr__(self, "factors", factors)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "hypotheses", hypotheses)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "model", model)
        object.__setattr__(self, "semantic_relation", relation)


@dataclass(frozen=True)
class SchurSandwichExpression:
    """Typed WH--Mellin-PDO--WH output tied to its AST product."""

    left_wiener_hopf: WienerHopfOperator
    mellin_operator: MellinPseudodifferentialOperator
    right_wiener_hopf: WienerHopfOperator
    ast_product: Product = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.left_wiener_hopf, WienerHopfOperator) or not isinstance(
            self.right_wiener_hopf, WienerHopfOperator
        ):
            raise TypeError("the exterior factors must be WienerHopfOperator objects.")
        if not isinstance(self.mellin_operator, MellinPseudodifferentialOperator):
            raise TypeError("the middle factor must be a MellinPseudodifferentialOperator.")
        object.__setattr__(
            self,
            "ast_product",
            Product(
                (
                    self.left_wiener_hopf.atom,
                    self.mellin_operator.atom,
                    self.right_wiener_hopf.atom,
                )
            ),
        )


@dataclass(frozen=True)
class SchurCancellationClassification:
    dilation_cancellation: CancellationLayerStatus
    initial_block_representations: CancellationLayerStatus
    regularizer_mellin_representation: CancellationLayerStatus
    mixed_algebra_membership: CancellationLayerStatus

    def __post_init__(self) -> None:
        if not all(
            isinstance(item, CancellationLayerStatus)
            for item in (
                self.dilation_cancellation,
                self.initial_block_representations,
                self.regularizer_mellin_representation,
                self.mixed_algebra_membership,
            )
        ):
            raise TypeError("all classification fields must be CancellationLayerStatus.")


@dataclass(frozen=True)
class SchurCancellationEvidence:
    """Evidence bundle that keeps every input relation and the exact core."""

    left_block_relation: object
    regularizer_relation: object
    right_block_relation: object
    exact_covariance: ExactIdentity
    compact_ideal_argument: tuple[str, ...]


@dataclass(frozen=True)
class SchurDilationCancellationResult:
    """Complete result, including original expression and open obligations."""

    expression: Product
    relation: DerivationRelationKind | None
    status: SchurCancellationStatus
    assumptions: AssumptionContext
    evidence: tuple[object, ...]
    proof_obligations: tuple[AnalyticProofObligation, ...]
    latex: str
    original_expression: Product
    model_expression: Product | None
    semantic_relation: ExactIdentity | FormalIdentity | ModCompactEquivalence | None
    classification: SchurCancellationClassification
    algebra_membership: AlgebraMembershipClaim
    scaled_regularizer: MellinPseudodifferentialOperator | None
    blocked_reasons: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.expression, Product) or not isinstance(
            self.original_expression, Product
        ):
            raise TypeError("expression and original_expression must be Product objects.")
        if self.model_expression is not None and not isinstance(
            self.model_expression, Product
        ):
            raise TypeError("model_expression must be a Product or None.")
        if self.relation is not None and not isinstance(
            self.relation, DerivationRelationKind
        ):
            raise TypeError("relation must be a DerivationRelationKind or None.")
        if self.semantic_relation is not None and not isinstance(
            self.semantic_relation,
            (ExactIdentity, FormalIdentity, ModCompactEquivalence),
        ):
            raise TypeError(
                "semantic_relation must be exact, formal, modulo compact, or None."
            )
        if not isinstance(self.status, SchurCancellationStatus):
            raise TypeError("status must be a SchurCancellationStatus.")
        if not isinstance(self.assumptions, AssumptionContext):
            raise TypeError("assumptions must be an AssumptionContext.")
        if not isinstance(self.classification, SchurCancellationClassification):
            raise TypeError("classification must be a SchurCancellationClassification.")
        if not isinstance(self.algebra_membership, AlgebraMembershipClaim):
            raise TypeError("algebra_membership must be an AlgebraMembershipClaim.")
        if not all(
            isinstance(item, AnalyticProofObligation)
            for item in self.proof_obligations
        ):
            raise TypeError("proof_obligations must contain AnalyticProofObligation objects.")
        if self.status is SchurCancellationStatus.BLOCKED_WITH_PRECISE_OBLIGATION:
            if self.relation is not None or self.semantic_relation is not None:
                raise SchurDilationError(
                    "a blocked result cannot carry a derived semantic relation."
                )
            if self.scaled_regularizer is not None or self.expression != self.original_expression:
                raise SchurDilationError(
                    "a blocked result must retain the original expression unchanged."
                )
            if not self.blocked_reasons:
                raise SchurDilationError("a blocked result requires a precise reason.")
        else:
            if self.relation is None or self.semantic_relation is None:
                raise SchurDilationError(
                    "a nonblocked result requires an explicit semantic relation."
                )
            if self.scaled_regularizer is None:
                raise SchurDilationError(
                    "a nonblocked result must retain its scaled regularizer."
                )
            if (
                self.semantic_relation.left != self.original_expression
                or self.semantic_relation.right != self.expression
            ):
                raise SchurDilationError(
                    "the semantic relation endpoints must match the result ASTs."
                )
        if self.status is SchurCancellationStatus.EXACT_IDENTITY and not isinstance(
            self.semantic_relation, ExactIdentity
        ):
            raise SchurDilationError("EXACT_IDENTITY requires an ExactIdentity.")
        if (
            self.status is SchurCancellationStatus.EQUALITY_MODULO_COMPACTS
            and not isinstance(self.semantic_relation, ModCompactEquivalence)
        ):
            raise SchurDilationError(
                "EQUALITY_MODULO_COMPACTS requires ModCompactEquivalence."
            )
        if (
            self.status
            is SchurCancellationStatus.CONDITIONAL_ON_REGULARIZER_REPRESENTATION
            and not isinstance(self.semantic_relation, FormalIdentity)
        ):
            raise SchurDilationError(
                "a conditional regularizer result must remain a FormalIdentity."
            )

    @property
    def succeeded(self) -> bool:
        return self.status is not SchurCancellationStatus.BLOCKED_WITH_PRECISE_OBLIGATION


def _base_obligations(
    *,
    regularizer_is_conditional: bool,
) -> tuple[AnalyticProofObligation, ...]:
    regularizer_statement = (
        "Prove the assumed R11 Mellin-PDO representation in the article's "
        "weighted operator space."
        if regularizer_is_conditional
        else "Identify and verify the article-level regularizer Mellin representation "
        "used in the Schur word."
    )
    return (
        AnalyticProofObligation(
            "regularizer_mellin_representation",
            regularizer_statement,
        ),
        AnalyticProofObligation(
            "wh_mellin_wh_sandwich_membership",
            "Prove that the final WH--Mellin-PDO--WH sandwich belongs to a closed "
            "operator algebra.",
        ),
        AnalyticProofObligation(
            "fredholm_algebra_membership",
            "Prove membership in an algebra carrying the required Fredholm calculus.",
        ),
        AnalyticProofObligation(
            "schur_correction_symbol",
            "Construct the symbol of the Schur correction only after algebra membership.",
        ),
    )


def _regularizer_atom(
    regularizer: RegularizerMellinRepresentation | RegularizerOperator,
) -> OperatorAtom:
    abstract = (
        regularizer.regularizer
        if isinstance(regularizer, RegularizerMellinRepresentation)
        else regularizer
    )
    if not isinstance(abstract, RegularizerOperator):
        raise TypeError(
            "regularizer must be a RegularizerOperator or RegularizerMellinRepresentation."
        )
    if not isinstance(abstract.operator, OperatorAtom):
        raise TypeError("the abstract regularizer must retain an OperatorAtom.")
    return abstract.operator


def _input_atom(value: SchurBlockRepresentation | OperatorAtom, name: str) -> OperatorAtom:
    if isinstance(value, SchurBlockRepresentation):
        return value.atom
    if isinstance(value, OperatorAtom):
        return value
    raise TypeError(f"{name} must be a SchurBlockRepresentation or OperatorAtom.")


def _blocked_result(
    *,
    original: Product,
    assumptions: AssumptionContext,
    reasons: tuple[str, ...],
    left_block: SchurBlockRepresentation | OperatorAtom,
    regularizer: RegularizerMellinRepresentation | RegularizerOperator,
    right_block: SchurBlockRepresentation | OperatorAtom,
    model: Product | None = None,
) -> SchurDilationCancellationResult:
    regularizer_is_conditional = not isinstance(
        regularizer, RegularizerMellinRepresentation
    ) or regularizer.status is RegularizerRepresentationStatus.ASSUMED
    obligations = (
        *_base_obligations(
            regularizer_is_conditional=regularizer_is_conditional
        ),
        *(
            AnalyticProofObligation(f"blocked_{index}", reason)
            for index, reason in enumerate(reasons, start=1)
        ),
    )
    membership = AlgebraMembershipClaim(
        expression=model or original,
        algebra="Fredholm algebra for mixed Wiener--Hopf--Mellin products",
        status=AlgebraMembershipStatus.UNPROVED,
        assumptions=assumptions.assumptions,
        source="no membership theorem is applied by the cancellation API",
    )
    block_relations = tuple(
        item.semantic_relation
        for item in (left_block, right_block)
        if isinstance(item, SchurBlockRepresentation)
    )
    if len(block_relations) != 2:
        block_layer = CancellationLayerStatus.NOT_APPLIED
    elif any(
        item.relation_kind is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
        for item in (left_block, right_block)
        if isinstance(item, SchurBlockRepresentation)
    ):
        block_layer = CancellationLayerStatus.EQUALITY_MODULO_COMPACTS
    else:
        block_layer = CancellationLayerStatus.EXACT_IDENTITY
    if not isinstance(regularizer, RegularizerMellinRepresentation) or (
        regularizer.status is RegularizerRepresentationStatus.ASSUMED
    ):
        regularizer_layer = (
            CancellationLayerStatus.CONDITIONAL_ON_REGULARIZER_REPRESENTATION
        )
    elif regularizer.status is RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT:
        regularizer_layer = CancellationLayerStatus.EQUALITY_MODULO_COMPACTS
    else:
        regularizer_layer = CancellationLayerStatus.EXACT_IDENTITY
    input_evidence = (
        *block_relations,
        *(
            (regularizer.semantic_relation,)
            if isinstance(regularizer, RegularizerMellinRepresentation)
            else ()
        ),
    )
    classification = SchurCancellationClassification(
        CancellationLayerStatus.NOT_APPLIED,
        block_layer,
        regularizer_layer,
        CancellationLayerStatus.UNPROVED_MIXED_ALGEBRA_MEMBERSHIP,
    )
    rendered_atoms = tuple(
        factor.name.replace("_", r"\_") for factor in original.factors
    )
    latex = "\\,".join(
        rf"\operatorname{{{rendered}}}" for rendered in rendered_atoms
    )
    return SchurDilationCancellationResult(
        expression=original,
        relation=None,
        status=SchurCancellationStatus.BLOCKED_WITH_PRECISE_OBLIGATION,
        assumptions=assumptions,
        evidence=input_evidence,
        proof_obligations=tuple(obligations),
        latex=latex,
        original_expression=original,
        model_expression=model,
        semantic_relation=None,
        classification=classification,
        algebra_membership=membership,
        scaled_regularizer=None,
        blocked_reasons=reasons,
    )


def _validate_block_roles(
    left: SchurBlockRepresentation,
    right: SchurBlockRepresentation,
) -> tuple[str, ...]:
    reasons: list[str] = []
    if left.exact_block != ExactBlock("A", 2, 1):
        reasons.append("left_block must represent exact block A21.")
    if right.exact_block != ExactBlock("A", 1, 2):
        reasons.append("right_block must represent exact block A12.")
    if not left.evidence or not right.evidence:
        reasons.append("both block representations require explicit evidence.")
    if not left.bounded or not right.bounded:
        reasons.append("compact-ideal propagation requires bounded block models.")
    if not left.localization_controlled or not right.localization_controlled:
        reasons.append("localization by chi_infinity is not controlled under dilation.")
    return tuple(reasons)


def _typed_pattern(
    left: SchurBlockRepresentation,
    right: SchurBlockRepresentation,
) -> tuple[
    WienerHopfOperator | None,
    WeightedDilationOperator | None,
    WeightedDilationOperator | None,
    WienerHopfOperator | None,
    tuple[str, ...],
]:
    reasons: list[str] = []
    if len(left.factors) != 2:
        reasons.append("left block must have exactly W(b21), V_gamma order.")
    if len(right.factors) != 2:
        reasons.append("right block must have exactly V_gamma_inverse, W(b12) order.")
    if reasons:
        return None, None, None, None, tuple(reasons)
    left_wh, left_dilation = left.factors
    right_dilation, right_wh = right.factors
    if not isinstance(left_wh, WienerHopfOperator):
        reasons.append("left block first factor is not a WienerHopfOperator.")
    if not isinstance(left_dilation, WeightedDilationOperator):
        reasons.append("left block second factor is not a WeightedDilationOperator.")
    if not isinstance(right_dilation, WeightedDilationOperator):
        reasons.append("right block first factor is not a WeightedDilationOperator.")
    if not isinstance(right_wh, WienerHopfOperator):
        reasons.append("right block second factor is not a WienerHopfOperator.")
    return (
        left_wh if isinstance(left_wh, WienerHopfOperator) else None,
        left_dilation if isinstance(left_dilation, WeightedDilationOperator) else None,
        right_dilation if isinstance(right_dilation, WeightedDilationOperator) else None,
        right_wh if isinstance(right_wh, WienerHopfOperator) else None,
        tuple(reasons),
    )


def _operators_are_bounded(*operators: object) -> bool:
    return all(getattr(operator, "bounded", False) is True for operator in operators)


def _common_space_reason(
    space: WeightedLpSpace,
    *operators: object,
) -> str | None:
    for operator in operators:
        if getattr(operator, "domain", None) != space or getattr(
            operator, "codomain", None
        ) != space:
            return (
                "all factors must have identical domain, codomain, p, weight, "
                "and norm convention."
            )
    return None


def _common_compact_ideal(
    left: SchurBlockRepresentation,
    regularizer: RegularizerMellinRepresentation,
    right: SchurBlockRepresentation,
) -> tuple[object | None, str | None]:
    ideals: list[object | None] = []
    for block in (left, right):
        if (
            block.relation_kind
            is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
        ):
            ideals.append(block.semantic_relation.compact_ideal)
    if regularizer.status is RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT:
        if not isinstance(regularizer.semantic_relation, ModCompactEquivalence):
            raise SchurDilationError(
                "a certified modulo-compact regularizer lost its semantic relation."
            )
        ideals.append(regularizer.semantic_relation.compact_ideal)
    if not ideals:
        return None, None
    if any(ideal is None for ideal in ideals):
        return None, "every modulo-compact input must name its compact ideal."
    common = ideals[0]
    if any((ideal == common) is not True for ideal in ideals[1:]):
        return None, (
            "modulo-compact inputs use incompatible compact ideals; a common "
            "ideal must be supplied explicitly."
        )
    return common, None


def _render_success_latex(expression: SchurSandwichExpression) -> str:
    left = sp.latex(expression.left_wiener_hopf.symbol)
    middle = sp.latex(expression.mellin_operator.symbol.as_expr())
    right = sp.latex(expression.right_wiener_hopf.symbol)
    return (
        rf"W\!\left({left}\right)\,"
        rf"\operatorname{{Op}}_M\!\left({middle}\right)\,"
        rf"W\!\left({right}\right)"
    )


def cancel_relative_dilations_in_schur_correction(
    *,
    left_block: SchurBlockRepresentation | OperatorAtom,
    regularizer: RegularizerMellinRepresentation | RegularizerOperator,
    right_block: SchurBlockRepresentation | OperatorAtom,
    assumptions: AssumptionContext,
) -> SchurDilationCancellationResult:
    """Cancel one explicitly adjacent inverse-dilation pair, or stop safely.

    The only accepted model pattern is

    ``W(b21), V_gamma, Op_M(r11), V_gamma_inverse, W(b12)``.

    No factor is commuted to create that pattern.  A blocked call returns the
    unchanged original expression plus precise proof obligations.
    """

    if not isinstance(assumptions, AssumptionContext):
        raise TypeError("assumptions must be an AssumptionContext.")
    if assumptions.consistency_status is ConsistencyStatus.INCONSISTENT:
        raise SchurDilationError("cancellation cannot use inconsistent assumptions.")
    left_atom = _input_atom(left_block, "left_block")
    right_atom = _input_atom(right_block, "right_block")
    regularizer_atom = _regularizer_atom(regularizer)
    original = Product((left_atom, regularizer_atom, right_atom))

    if not isinstance(left_block, SchurBlockRepresentation) or not isinstance(
        right_block, SchurBlockRepresentation
    ):
        return _blocked_result(
            original=original,
            assumptions=assumptions,
            left_block=left_block,
            regularizer=regularizer,
            right_block=right_block,
            reasons=(
                "explicit left and right block representations are required; "
                "raw AST atoms are never expanded implicitly.",
            ),
        )
    model_without_regularizer = Product(
        (*left_block.model.factors, regularizer_atom, *right_block.model.factors)
    )
    reasons = list(_validate_block_roles(left_block, right_block))
    left_wh, left_dilation, right_dilation, right_wh, pattern_reasons = _typed_pattern(
        left_block, right_block
    )
    reasons.extend(pattern_reasons)
    if not isinstance(regularizer, RegularizerMellinRepresentation):
        reasons.append(
            "prove R11 has an explicit Mellin-PDO representation; the regularizer "
            "remains abstract."
        )
    if reasons:
        return _blocked_result(
            original=original,
            assumptions=assumptions,
            left_block=left_block,
            regularizer=regularizer,
            right_block=right_block,
            reasons=tuple(dict.fromkeys(reasons)),
            model=model_without_regularizer,
        )

    assert left_wh is not None
    assert left_dilation is not None
    assert right_dilation is not None
    assert right_wh is not None
    assert isinstance(regularizer, RegularizerMellinRepresentation)
    pdo = regularizer.mellin_operator
    model = Product(
        (
            left_wh.atom,
            left_dilation.atom,
            pdo.atom,
            right_dilation.atom,
            right_wh.atom,
        )
    )
    if regularizer.regularizer.operator != regularizer_atom:
        reasons.append("regularizer representation does not reference the input operator.")
    if regularizer.regularizer.target != ExactBlock("A", 1, 1):
        reasons.append("the regularizer must target exact block A11.")
    if not _operators_are_bounded(
        left_wh,
        left_dilation,
        pdo,
        right_dilation,
        right_wh,
    ):
        reasons.append("every factor must be explicitly bounded.")
    if left_wh.localized or right_wh.localized:
        reasons.append(
            "localized Wiener--Hopf factors are not accepted: the current model "
            "does not represent and certify the cutoff rescaling "
            "chi(x) -> chi(gamma*x)."
        )
    common_reason = _common_space_reason(
        pdo.domain,
        left_wh,
        left_dilation,
        pdo,
        right_dilation,
        right_wh,
    )
    if common_reason is not None:
        reasons.append(common_reason)
    if not pdo.radial_scaling_stable or not pdo.stability_evidence:
        reasons.append("the Mellin symbol class lacks radial-scaling stability evidence.")
    if regularizer.status is RegularizerRepresentationStatus.ASSUMED and not (
        regularizer.hypotheses
    ):
        reasons.append("an assumed regularizer representation needs explicit hypotheses.")
    is_mod_compact = any(
        relation is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
        for relation in (left_block.relation_kind, right_block.relation_kind)
    ) or (
        regularizer.status
        is RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT
    )
    compact_ideal, compact_ideal_reason = _common_compact_ideal(
        left_block,
        regularizer,
        right_block,
    )
    if compact_ideal_reason is not None:
        reasons.append(compact_ideal_reason)
    if reasons:
        return _blocked_result(
            original=original,
            assumptions=assumptions,
            left_block=left_block,
            regularizer=regularizer,
            right_block=right_block,
            reasons=tuple(dict.fromkeys(reasons)),
            model=model,
        )

    try:
        covariance = conjugate_mellin_pdo_by_dilation(
            left_dilation,
            pdo,
            right_dilation,
            assumptions,
        )
    except (MellinOperatorError, MellinDomainError, IncompatibleDomainError) as exc:
        return _blocked_result(
            original=original,
            assumptions=assumptions,
            left_block=left_block,
            regularizer=regularizer,
            right_block=right_block,
            reasons=(str(exc),),
            model=model,
        )
    scaled = covariance.transformed_operator
    if not isinstance(scaled, MellinPseudodifferentialOperator):
        raise SchurDilationError("Mellin covariance returned the wrong operator type.")
    sandwich = SchurSandwichExpression(left_wh, scaled, right_wh)

    regularizer_conditional = (
        regularizer.status is RegularizerRepresentationStatus.ASSUMED
    )
    compact_argument = (
        "AKB is compact for bounded A,B and compact K",
        "finite sums of compact operators are compact",
    )
    evidence_bundle = SchurCancellationEvidence(
        left_block.semantic_relation,
        regularizer.semantic_relation,
        right_block.semantic_relation,
        covariance.exact_identity,
        compact_argument if is_mod_compact else (),
    )
    hypotheses = (
        *left_block.hypotheses,
        *regularizer.hypotheses,
        *right_block.hypotheses,
        *covariance.assumptions.assumptions,
    )
    semantic_relation: ExactIdentity | FormalIdentity | ModCompactEquivalence
    relation_kind: DerivationRelationKind
    if regularizer_conditional:
        relation_kind = (
            DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
            if is_mod_compact
            else DerivationRelationKind.EXACT_EQUALITY
        )
        semantic_relation = FormalIdentity(
            original,
            sandwich.ast_product,
            justification=(
                "candidate conclusion only under the explicitly assumed "
                "regularizer Mellin representation",
                relation_kind,
                evidence_bundle,
                hypotheses,
            ),
        )
    elif is_mod_compact:
        if compact_ideal is None:
            raise SchurDilationError(
                "a modulo-compact result requires the validated common ideal."
            )
        relation_kind = DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
        semantic_relation = ModCompactEquivalence(
            original,
            sandwich.ast_product,
            space=pdo.domain,
            compact_ideal=compact_ideal,
            evidence=evidence_bundle,
        )
    else:
        relation_kind = DerivationRelationKind.EXACT_EQUALITY
        semantic_relation = ExactIdentity(
            original,
            sandwich.ast_product,
            evidence=evidence_bundle,
            scope=ExactIdentityScope.OPERATORIAL,
            hypotheses=hypotheses,
        )

    block_layer = (
        CancellationLayerStatus.EQUALITY_MODULO_COMPACTS
        if any(
            item.relation_kind
            is DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
            for item in (left_block, right_block)
        )
        else CancellationLayerStatus.EXACT_IDENTITY
    )
    if regularizer.status is RegularizerRepresentationStatus.ASSUMED:
        regularizer_layer = (
            CancellationLayerStatus.CONDITIONAL_ON_REGULARIZER_REPRESENTATION
        )
    elif regularizer.status is RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT:
        regularizer_layer = CancellationLayerStatus.EQUALITY_MODULO_COMPACTS
    else:
        regularizer_layer = CancellationLayerStatus.EXACT_IDENTITY
    classification = SchurCancellationClassification(
        CancellationLayerStatus.EXACT_IDENTITY,
        block_layer,
        regularizer_layer,
        CancellationLayerStatus.UNPROVED_MIXED_ALGEBRA_MEMBERSHIP,
    )
    if regularizer_conditional:
        status = (
            SchurCancellationStatus.CONDITIONAL_ON_REGULARIZER_REPRESENTATION
        )
    elif is_mod_compact:
        status = SchurCancellationStatus.EQUALITY_MODULO_COMPACTS
    else:
        status = SchurCancellationStatus.EXACT_IDENTITY
    membership = AlgebraMembershipClaim(
        expression=sandwich.ast_product,
        algebra="Fredholm algebra for mixed Wiener--Hopf--Mellin products",
        status=AlgebraMembershipStatus.UNPROVED,
        assumptions=hypotheses,
        source="the cancellation identity does not imply mixed-algebra membership",
    )
    obligations = _base_obligations(
        regularizer_is_conditional=regularizer_conditional
    )
    return SchurDilationCancellationResult(
        expression=sandwich.ast_product,
        relation=relation_kind,
        status=status,
        assumptions=covariance.assumptions,
        evidence=(evidence_bundle,),
        proof_obligations=obligations,
        latex=_render_success_latex(sandwich),
        original_expression=original,
        model_expression=model,
        semantic_relation=semantic_relation,
        classification=classification,
        algebra_membership=membership,
        scaled_regularizer=scaled,
    )


__all__ = [
    "BlockRepresentationEvidence",
    "CancellationLayerStatus",
    "SchurBlockRepresentation",
    "SchurCancellationClassification",
    "SchurCancellationEvidence",
    "SchurCancellationStatus",
    "SchurDilationCancellationResult",
    "SchurDilationError",
    "SchurSandwichExpression",
    "cancel_relative_dilations_in_schur_correction",
]
