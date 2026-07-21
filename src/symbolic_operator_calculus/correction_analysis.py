"""Controlled expansion and structural analysis of the complete correction.

All authoritative algebraic expressions use the existing operator AST. The
types in this module annotate terms, exact contiguous motifs, and unresolved
interfaces; they do not form a second expression tree and do not prove class
membership or analytic closure.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import Enum

from .complete_correction import (
    CompleteCorrectionFactorization,
    build_complete_correction_derivation,
)
from .normalized_schur import t1_minus_expression, t2_minus_expression
from .operators import (
    Ghat1,
    Ghat2,
    P_minus,
    P_plus,
    R1,
    T1_minus,
    T2_minus,
    U1,
    U1_inverse,
    U2,
    U2_inverse,
    Wminus_21,
    Wplus_12,
    Z1_inverse,
    Z2_inverse,
    LinearCombination,
    OperatorAtom,
    Product,
    Scalar,
    Term,
    compose_ordered,
)
from .semantics import AnalyticProofObligation, DerivationRelationKind


class CorrectionAnalysisError(ValueError):
    """Raised when controlled expansion metadata is inconsistent."""


class TermOperatorClass(str, Enum):
    """Declared term-signature classes; no membership proof is implied."""

    MULTIPLICATION_OPERATOR = "MultiplicationOperator"
    NORMALIZED_COEFFICIENT = "NormalizedCoefficient"
    DILATION_OPERATOR = "DilationOperator"
    CAUCHY_FACTOR_PLUS = "CauchyFactorPlus"
    CAUCHY_FACTOR_MINUS = "CauchyFactorMinus"
    LOCALIZED_WIENER_HOPF_PLUS = "LocalizedWienerHopfPlus"
    LOCALIZED_WIENER_HOPF_MINUS = "LocalizedWienerHopfMinus"
    AUXILIARY_SHIFT_OPERATOR = "AuxiliaryShiftOperator"
    MELLIN_PDO_REGULARIZER = "MellinPDORegularizer"
    UNKNOWN_COMPOSITE = "UnknownComposite"


class MotifStatus(str, Enum):
    """Logical status of one exactly located contiguous motif."""

    EXACT_RECOGNIZED = "EXACT_RECOGNIZED"
    SUPPLIED_MOD_COMPACT = "SUPPLIED_MOD_COMPACT"
    UNRESOLVED = "UNRESOLVED"
    FORBIDDEN_WITHOUT_COMMUTATION = "FORBIDDEN_WITHOUT_COMMUTATION"


class InterfaceStatus(str, Enum):
    """Allowed status values for contiguous operator interfaces."""

    CERTIFIED_EXACT = "CERTIFIED_EXACT"
    CERTIFIED_MOD_COMPACT = "CERTIFIED_MOD_COMPACT"
    SOURCE_VERIFIED_RULE_MISSING = "SOURCE_VERIFIED_RULE_MISSING"
    SOURCE_UNVERIFIED = "SOURCE_UNVERIFIED"
    NO_RULE_FOUND = "NO_RULE_FOUND"
    NOT_APPLICABLE = "NOT_APPLICABLE"


@dataclass(frozen=True)
class FactorSignature:
    """One concrete factor and its declarative semantics at a fixed position."""

    position: int
    factor: OperatorAtom
    operator_class: TermOperatorClass
    branch: int | None
    polarity: str | None
    relation_kind: DerivationRelationKind

    def __post_init__(self) -> None:
        if not isinstance(self.position, int) or isinstance(self.position, bool):
            raise TypeError("position must be a positive integer.")
        if self.position <= 0:
            raise ValueError("position must be positive.")
        if not isinstance(self.factor, OperatorAtom):
            raise TypeError("factor must be an OperatorAtom.")
        if not isinstance(self.operator_class, TermOperatorClass):
            raise TypeError("operator_class must be a TermOperatorClass.")
        if self.branch is not None and self.branch not in (1, 2):
            raise ValueError("branch must be 1, 2, or None.")
        if self.polarity not in (None, "+", "-"):
            raise ValueError("polarity must be '+', '-', or None.")
        if not isinstance(self.relation_kind, DerivationRelationKind):
            raise TypeError("relation_kind must be a DerivationRelationKind.")


@dataclass(frozen=True)
class MotifMatch:
    """A contiguous match found without rearranging the operator word."""

    name: str
    status: MotifStatus
    start_position: int
    end_position: int
    factors: tuple[OperatorAtom, ...]
    source: str

    def __post_init__(self) -> None:
        if not self.name or not self.source:
            raise ValueError("motif name and source must be non-empty.")
        if not isinstance(self.status, MotifStatus):
            raise TypeError("status must be a MotifStatus.")
        if self.start_position <= 0 or self.end_position < self.start_position:
            raise ValueError("motif positions must define a positive interval.")
        if len(self.factors) != self.end_position - self.start_position + 1:
            raise ValueError("motif factors and positions disagree.")


@dataclass(frozen=True)
class OpenOperatorInterface:
    """One unresolved contiguous pair in an expanded term."""

    category: str
    position: int
    left: OperatorAtom
    right: OperatorAtom
    status: InterfaceStatus
    source: str

    def __post_init__(self) -> None:
        if not self.category or not self.source:
            raise ValueError("interface category and source must be non-empty.")
        if self.position <= 0:
            raise ValueError("interface position must be positive.")
        if not isinstance(self.left, OperatorAtom) or not isinstance(
            self.right,
            OperatorAtom,
        ):
            raise TypeError("interface endpoints must be OperatorAtoms.")
        if not isinstance(self.status, InterfaceStatus):
            raise TypeError("status must be an InterfaceStatus.")


@dataclass(frozen=True)
class CorrectionTermRecord:
    """Stable metadata for one authoritative AST term."""

    identifier: str
    term: Term
    signature: tuple[FactorSignature, ...]
    origin_relation: DerivationRelationKind
    origin: str
    motifs: tuple[MotifMatch, ...]
    open_interfaces: tuple[OpenOperatorInterface, ...]
    proof_obligations: tuple[AnalyticProofObligation, ...]

    def __post_init__(self) -> None:
        if not self.identifier.startswith("C2-T") or len(self.identifier) != 6:
            raise ValueError("identifier must have stable form C2-T01 through C2-T16.")
        if not isinstance(self.term, Term):
            raise TypeError("term must be a Term.")
        if not isinstance(self.origin_relation, DerivationRelationKind):
            raise TypeError("origin_relation must be a DerivationRelationKind.")
        if not self.origin:
            raise ValueError("origin must be non-empty.")
        if tuple(item.position for item in self.signature) != tuple(
            range(1, len(self.term.product.factors) + 1)
        ):
            raise CorrectionAnalysisError("signature positions are not consecutive.")
        if tuple(item.factor for item in self.signature) != self.term.product.factors:
            raise CorrectionAnalysisError("signature factors differ from the AST term.")

    @property
    def coefficient(self) -> Scalar:
        return self.term.coefficient


@dataclass(frozen=True)
class AuxiliaryGroupedTerm:
    """One C2 summand retaining both differences as grouped metadata."""

    identifier: str
    coefficient: Scalar
    factorization: CompleteCorrectionFactorization
    left_auxiliary_piece: Product
    right_auxiliary_piece: Product

    def __post_init__(self) -> None:
        if not self.identifier.startswith("C2-A"):
            raise ValueError("auxiliary group identifier must start with C2-A.")
        if self.coefficient != 1:
            raise CorrectionAnalysisError(
                "the supplied auxiliary decompositions have unit coefficients."
            )
        if not isinstance(self.factorization, CompleteCorrectionFactorization):
            raise TypeError("factorization must be a CompleteCorrectionFactorization.")
        if not isinstance(self.left_auxiliary_piece, Product) or not isinstance(
            self.right_auxiliary_piece,
            Product,
        ):
            raise TypeError("auxiliary pieces must be Products.")


@dataclass(frozen=True)
class MotifFrequency:
    name: str
    status: MotifStatus
    count: int


@dataclass(frozen=True)
class InterfaceSummary:
    """One row of the global contiguous-interface matrix."""

    interface: str
    term_count: int
    existing_rule: str
    source: str
    status: InterfaceStatus


@dataclass(frozen=True)
class CompleteCorrectionExpansionTrace:
    """The controlled C0/C1/C2/C3 representations and their annotations."""

    c0: CompleteCorrectionFactorization
    c1: LinearCombination
    c2: tuple[AuxiliaryGroupedTerm, ...]
    c3: LinearCombination
    terms: tuple[CorrectionTermRecord, ...]
    motif_frequency: tuple[MotifFrequency, ...]
    interface_matrix: tuple[InterfaceSummary, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.c0, CompleteCorrectionFactorization):
            raise TypeError("c0 must be factorization metadata.")
        if not isinstance(self.c1, LinearCombination) or len(self.c1.terms) != 4:
            raise CorrectionAnalysisError("C1 must contain exactly four AST terms.")
        if len(self.c2) != 4 or not all(
            isinstance(item, AuxiliaryGroupedTerm) for item in self.c2
        ):
            raise CorrectionAnalysisError("C2 must contain four grouped terms.")
        if not isinstance(self.c3, LinearCombination) or len(self.c3.terms) != 16:
            raise CorrectionAnalysisError("C3 must contain exactly 16 AST terms.")
        if len(self.terms) != 16 or tuple(item.term for item in self.terms) != (
            self.c3.terms
        ):
            raise CorrectionAnalysisError("term catalog must annotate C3 in order.")
        if any(
            term.term.product.factors.count(R1) != 1 for term in self.terms
        ):
            raise CorrectionAnalysisError("every C3 term must keep one atomic R1.")

    @property
    def counts(self) -> tuple[int, int, int, int]:
        return (1, len(self.c1.terms), len(self.c2), len(self.c3.terms))


def build_factor_signature(position: int, factor: OperatorAtom) -> FactorSignature:
    """Build declarative metadata for one existing AST factor.

    The returned class and relation labels are annotations only. They do not
    establish analytic membership or discharge any proof obligation.
    """
    if factor in (Z1_inverse, Z2_inverse):
        operator_class = TermOperatorClass.MULTIPLICATION_OPERATOR
    elif factor in (Ghat1, Ghat2):
        operator_class = TermOperatorClass.NORMALIZED_COEFFICIENT
    elif factor in (U1, U2, U1_inverse, U2_inverse):
        operator_class = TermOperatorClass.DILATION_OPERATOR
    elif factor is P_plus:
        operator_class = TermOperatorClass.CAUCHY_FACTOR_PLUS
    elif factor is P_minus:
        operator_class = TermOperatorClass.CAUCHY_FACTOR_MINUS
    elif factor is Wplus_12:
        operator_class = TermOperatorClass.LOCALIZED_WIENER_HOPF_PLUS
    elif factor is Wminus_21:
        operator_class = TermOperatorClass.LOCALIZED_WIENER_HOPF_MINUS
    elif factor in (T1_minus, T2_minus):
        operator_class = TermOperatorClass.AUXILIARY_SHIFT_OPERATOR
    elif factor is R1:
        operator_class = TermOperatorClass.MELLIN_PDO_REGULARIZER
    else:
        operator_class = TermOperatorClass.UNKNOWN_COMPOSITE

    branch_by_factor = {
        U1: 1,
        U1_inverse: 1,
        Ghat1: 1,
        Wplus_12: 1,
        T1_minus: 1,
        R1: 1,
        U2: 2,
        U2_inverse: 2,
        Ghat2: 2,
        Wminus_21: 2,
        T2_minus: 2,
        Z1_inverse: 1,
        Z2_inverse: 2,
    }
    polarity_by_factor = {
        P_plus: "+",
        Wplus_12: "+",
        P_minus: "-",
        Wminus_21: "-",
    }
    if factor in (Wplus_12, Wminus_21):
        relation_kind = DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
    elif factor is R1:
        relation_kind = DerivationRelationKind.FORMAL_SUBSTITUTION
    else:
        relation_kind = DerivationRelationKind.EXACT_EQUALITY
    return FactorSignature(
        position=position,
        factor=factor,
        operator_class=operator_class,
        branch=branch_by_factor.get(factor),
        polarity=polarity_by_factor.get(factor),
        relation_kind=relation_kind,
    )


def _motif(
    name: str,
    status: MotifStatus,
    start_index: int,
    factors: tuple[OperatorAtom, ...],
    source: str,
) -> MotifMatch:
    return MotifMatch(
        name=name,
        status=status,
        start_position=start_index + 1,
        end_position=start_index + len(factors),
        factors=factors,
        source=source,
    )


def detect_correction_motifs(product: Product) -> tuple[MotifMatch, ...]:
    """Detect only exact contiguous patterns in their stored order."""

    if not isinstance(product, Product):
        raise TypeError("product must be a Product.")
    factors = product.factors
    matches: list[MotifMatch] = []
    for index in range(len(factors) - 1):
        pair = factors[index : index + 2]
        if pair[0] in (U2, Ghat2) and pair[1] is Wminus_21:
            matches.append(
                _motif(
                    "normalized_left_shift",
                    MotifStatus.SUPPLIED_MOD_COMPACT,
                    index,
                    pair,
                    "Phase N WH model plus exact branch-2 normalization",
                )
            )
        if pair[0] in (U1, Ghat1) and pair[1] is Wplus_12:
            matches.append(
                _motif(
                    "normalized_right_shift",
                    MotifStatus.SUPPLIED_MOD_COMPACT,
                    index,
                    pair,
                    "Phase N WH model plus exact branch-1 normalization",
                )
            )
        if pair in ((Wminus_21, P_minus), (Wplus_12, P_minus)):
            name = (
                "left_auxiliary_piece"
                if pair[0] is Wminus_21
                else "right_auxiliary_piece"
            )
            matches.append(
                _motif(
                    name,
                    MotifStatus.EXACT_RECOGNIZED,
                    index,
                    pair,
                    "exact supplied T_{k,-} expansion",
                )
            )
        if pair in ((U1_inverse, P_plus), (U2_inverse, P_plus)):
            matches.append(
                _motif(
                    "mixed_cauchy_dilation",
                    MotifStatus.EXACT_RECOGNIZED,
                    index,
                    pair,
                    "contiguous U_k^{-1} P^+ piece of the exact T_{k,-} expansion",
                )
            )
    for index in range(len(factors) - 2):
        triple = factors[index : index + 3]
        if triple in (
            (Wminus_21, U1_inverse, P_plus),
            (Wplus_12, U2_inverse, P_plus),
        ):
            name = (
                "left_auxiliary_piece"
                if triple[0] is Wminus_21
                else "right_auxiliary_piece"
            )
            matches.append(
                _motif(
                    name,
                    MotifStatus.EXACT_RECOGNIZED,
                    index,
                    triple,
                    "exact supplied T_{k,-} expansion",
                )
            )
    regularizer_index = factors.index(R1)
    core = factors[regularizer_index - 1 : regularizer_index + 2]
    matches.append(
        _motif(
            "mellin_regularizer_core",
            MotifStatus.UNRESOLVED,
            regularizer_index - 1,
            core,
            "R1 remains atomic; adjacent semiproduct rules are not established",
        )
    )
    return tuple(
        sorted(
            matches,
            key=lambda item: (item.start_position, item.end_position, item.name),
        )
    )


def _open_interface(
    category: str,
    position: int,
    left: OperatorAtom,
    right: OperatorAtom,
    status: InterfaceStatus,
    source: str,
) -> OpenOperatorInterface:
    return OpenOperatorInterface(category, position, left, right, status, source)


def detect_open_interfaces(product: Product) -> tuple[OpenOperatorInterface, ...]:
    """Return unresolved contiguous pairs; exact U^-1 P+ pieces are omitted."""

    if not isinstance(product, Product):
        raise TypeError("product must be a Product.")
    interfaces: list[OpenOperatorInterface] = []
    for index, (left, right) in enumerate(
        zip(product.factors, product.factors[1:]),
        start=1,
    ):
        item: OpenOperatorInterface | None = None
        if left in (Ghat1, Ghat2) and right in (Wminus_21, Wplus_12):
            item = _open_interface(
                "multiplication_x_wiener_hopf",
                index,
                left,
                right,
                InterfaceStatus.SOURCE_VERIFIED_RULE_MISSING,
                "paper verifies normalized factor; no semiproduct rule implemented",
            )
        elif left in (U1, U2) and right in (Wminus_21, Wplus_12):
            item = _open_interface(
                "dilation_x_wiener_hopf",
                index,
                left,
                right,
                InterfaceStatus.SOURCE_VERIFIED_RULE_MISSING,
                "paper verifies normalized factor; no closure rule implemented",
            )
        elif left in (Wminus_21, Wplus_12) and right in (
            U1_inverse,
            U2_inverse,
        ):
            item = _open_interface(
                "wiener_hopf_x_dilation",
                index,
                left,
                right,
                InterfaceStatus.SOURCE_VERIFIED_RULE_MISSING,
                "exact auxiliary expansion; WH semiproduct rule missing",
            )
        elif left in (Wminus_21, Wplus_12) and right is P_minus:
            item = _open_interface(
                "wiener_hopf_x_cauchy",
                index,
                left,
                right,
                InterfaceStatus.NO_RULE_FOUND,
                "no applicable contiguous interface rule is implemented",
            )
        elif left in (P_plus, P_minus) and right is R1:
            item = _open_interface(
                "cauchy_x_r1",
                index,
                left,
                right,
                InterfaceStatus.SOURCE_UNVERIFIED,
                "KKL 2014 applicability to this semiproduct remains to verify",
            )
        elif left is R1 and right in (Ghat1, Ghat2):
            item = _open_interface(
                "r1_x_multiplication",
                index,
                left,
                right,
                InterfaceStatus.SOURCE_UNVERIFIED,
                "Mellin semiproduct rule and symbol class remain to verify",
            )
        elif left is R1 and right in (U1, U2):
            item = _open_interface(
                "r1_x_dilation",
                index,
                left,
                right,
                InterfaceStatus.SOURCE_UNVERIFIED,
                "Mellin-dilation semiproduct rule remains to verify",
            )
        elif left in (U1_inverse, U2_inverse) and right is P_plus:
            item = None
        if item is not None:
            interfaces.append(item)
    return tuple(interfaces)


def _term_record(index: int, term: Term) -> CorrectionTermRecord:
    identifier = f"C2-T{index:02d}"
    signature = tuple(
        build_factor_signature(position, factor)
        for position, factor in enumerate(term.product.factors, start=1)
    )
    interfaces = detect_open_interfaces(term.product)
    obligations = tuple(
        AnalyticProofObligation(
            key=f"{identifier.lower().replace('-', '_')}_{category}",
            statement=(
                f"Justify the unresolved {category} interfaces in {identifier} "
                "without commuting factors."
            ),
        )
        for category in dict.fromkeys(item.category for item in interfaces)
    )
    return CorrectionTermRecord(
        identifier=identifier,
        term=term,
        signature=signature,
        origin_relation=DerivationRelationKind.EXACT_EQUALITY,
        origin="controlled distribution of two exact differences and two exact auxiliary sums",
        motifs=detect_correction_motifs(term.product),
        open_interfaces=interfaces,
        proof_obligations=obligations,
    )


def _motif_frequencies(
    records: tuple[CorrectionTermRecord, ...],
) -> tuple[MotifFrequency, ...]:
    counter = Counter(
        (match.name, match.status)
        for record in records
        for match in record.motifs
    )
    return tuple(
        MotifFrequency(name, status, count)
        for (name, status), count in sorted(
            counter.items(),
            key=lambda item: (item[0][0], item[0][1].value),
        )
    )


def _interface_matrix(
    records: tuple[CorrectionTermRecord, ...],
    c1: LinearCombination,
) -> tuple[InterfaceSummary, ...]:
    counts = Counter(
        interface.category
        for record in records
        for interface in record.open_interfaces
    )
    wh_auxiliary_count = sum(
        sum(
            pair in ((Wminus_21, T1_minus), (Wplus_12, T2_minus))
            for pair in zip(
                term.product.factors,
                term.product.factors[1:],
            )
        )
        for term in c1.terms
    )
    rows = (
        (
            "multiplicación × Wiener--Hopf",
            counts["multiplication_x_wiener_hopf"],
            "no implemented semiproduct rule",
            "paper exact normalization; closure rule absent",
            InterfaceStatus.SOURCE_VERIFIED_RULE_MISSING,
        ),
        (
            "dilatación × Wiener--Hopf",
            counts["dilation_x_wiener_hopf"],
            "no implemented closure rule",
            "paper exact normalization; closure rule absent",
            InterfaceStatus.SOURCE_VERIFIED_RULE_MISSING,
        ),
        (
            "Wiener--Hopf × dilatación",
            counts["wiener_hopf_x_dilation"],
            "exact auxiliary expansion only",
            "T_{k,-} definition; semiproduct rule absent",
            InterfaceStatus.SOURCE_VERIFIED_RULE_MISSING,
        ),
        (
            "Wiener--Hopf × Cauchy factor",
            counts["wiener_hopf_x_cauchy"],
            "none",
            "no rule found in the implemented API",
            InterfaceStatus.NO_RULE_FOUND,
        ),
        (
            "Cauchy factor × R_1",
            counts["cauchy_x_r1"],
            "none verified",
            "KKL 2014 source verification pending",
            InterfaceStatus.SOURCE_UNVERIFIED,
        ),
        (
            "Wiener--Hopf × R_1",
            0,
            "not contiguous in C3",
            "structural inspection",
            InterfaceStatus.NOT_APPLICABLE,
        ),
        (
            "R_1 × multiplicación",
            counts["r1_x_multiplication"],
            "none verified",
            "Mellin semiproduct source verification pending",
            InterfaceStatus.SOURCE_UNVERIFIED,
        ),
        (
            "R_1 × dilatación",
            counts["r1_x_dilation"],
            "none verified",
            "Mellin-dilation source verification pending",
            InterfaceStatus.SOURCE_UNVERIFIED,
        ),
        (
            "R_1 × Wiener--Hopf",
            0,
            "not contiguous in C3",
            "structural inspection",
            InterfaceStatus.NOT_APPLICABLE,
        ),
        (
            "Wiener--Hopf × auxiliar",
            wh_auxiliary_count,
            "T_{k,-}=U_k^{-1}P^+ + P^-",
            "exact supplied auxiliary identity; interface rule absent",
            InterfaceStatus.SOURCE_VERIFIED_RULE_MISSING,
        ),
        (
            "productos con cutoffs/localización explícitos",
            0,
            "no explicit cutoff atom in C3",
            "localized status is metadata on W factors",
            InterfaceStatus.NOT_APPLICABLE,
        ),
    )
    return tuple(InterfaceSummary(*row) for row in rows)


def build_complete_correction_expansion_trace(
) -> CompleteCorrectionExpansionTrace:
    """Build C0/C1/C2/C3 and all term-level structural annotations."""

    correction = build_complete_correction_derivation()
    c0 = correction.normalized_factorization
    c1 = c0.expanded
    left_pieces = tuple(term.product for term in t1_minus_expression().terms)
    right_pieces = tuple(term.product for term in t2_minus_expression().terms)
    grouped: list[AuxiliaryGroupedTerm] = []
    full_terms: list[Term] = []
    group_index = 1
    for left_piece in left_pieces:
        for right_piece in right_pieces:
            expanded = compose_ordered(
                c0.left_difference,
                c0.left_wiener_hopf,
                left_piece,
                R1,
                c0.right_difference,
                c0.right_wiener_hopf,
                right_piece,
            )
            factorization = CompleteCorrectionFactorization(
                expanded=expanded,
                left_difference=c0.left_difference,
                left_wiener_hopf=c0.left_wiener_hopf,
                bridge=Product(left_piece.factors + (R1,)),
                right_difference=c0.right_difference,
                right_wiener_hopf=c0.right_wiener_hopf,
                suffix=right_piece,
            )
            grouped.append(
                AuxiliaryGroupedTerm(
                    identifier=f"C2-A{group_index:02d}",
                    coefficient=1,
                    factorization=factorization,
                    left_auxiliary_piece=left_piece,
                    right_auxiliary_piece=right_piece,
                )
            )
            full_terms.extend(expanded.terms)
            group_index += 1
    c3 = LinearCombination(tuple(full_terms))
    records = tuple(
        _term_record(index, term)
        for index, term in enumerate(c3.terms, start=1)
    )
    return CompleteCorrectionExpansionTrace(
        c0=c0,
        c1=c1,
        c2=tuple(grouped),
        c3=c3,
        terms=records,
        motif_frequency=_motif_frequencies(records),
        interface_matrix=_interface_matrix(records, c1),
    )
