"""Evidence-aware structural model for the first-pivot closure problem.

The authoritative operator words remain :class:`Product` instances from the
existing noncommutative AST.  The records in this module only annotate source
evidence, unresolved interfaces, candidate lemmas, and dependencies.  They do
not infer Mellin membership, compactness, algebra closure, or Fredholmness.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from .correction_analysis import FactorSignature, build_factor_signature
from .operators import (
    Ghat1,
    P_minus,
    P_plus,
    R1,
    U1,
    U1_inverse,
    Wplus_12,
    OperatorAtom,
    Product,
)


class ClosureAnalysisError(ValueError):
    """Raised when closure metadata contradicts its authoritative AST word."""


class EvidenceApplicationKind(str, Enum):
    """Permitted source-to-problem application classifications."""

    DIRECT = "DIRECT"
    SPECIALIZATION = "SPECIALIZATION"
    ANALOGY = "ANALOGY"
    UNPROVED_ADAPTATION = "UNPROVED_ADAPTATION"
    NO_SOURCE_FOUND = "NO_SOURCE_FOUND"


class ClosureInterfaceStatus(str, Enum):
    """Evidence status for one closure-relevant interface."""

    CERTIFIED_EXACT = "CERTIFIED_EXACT"
    CERTIFIED_MOD_COMPACT = "CERTIFIED_MOD_COMPACT"
    SPECIALIZATION_TO_PROVE = "SPECIALIZATION_TO_PROVE"
    ANALOGY_ONLY = "ANALOGY_ONLY"
    NO_RULE = "NO_RULE"
    BLOCKED = "BLOCKED"


class ClosureObligationStatus(str, Enum):
    """Lifecycle labels that never imply discharge by construction."""

    OPEN = "OPEN"
    BLOCKED = "BLOCKED"
    SOURCE_VERIFIED = "SOURCE_VERIFIED"
    FORMALLY_REDUCED = "FORMALLY_REDUCED"
    ANALYTICALLY_PROVED = "ANALYTICALLY_PROVED"


class MinimalLemmaChoice(str, Enum):
    H1 = "H1"
    H2 = "H2"
    H3 = "H3"
    NONE = "NONE"


class DecisionConfidence(str, Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


def _nonempty(value: object, field_name: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{field_name} must be a non-empty str.")
    if not value.strip():
        raise ValueError(f"{field_name} must be a non-empty str.")
    return value


def _text_tuple(value: object, field_name: str) -> tuple[str, ...]:
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{field_name} must be an iterable of strings.")
    try:
        result = tuple(value)  # type: ignore[arg-type]
    except TypeError as exc:
        raise TypeError(f"{field_name} must be an iterable of strings.") from exc
    if not result or not all(isinstance(item, str) and item.strip() for item in result):
        raise ValueError(f"{field_name} must contain non-empty strings.")
    return result


@dataclass(frozen=True)
class EvidenceReference:
    """One verified literature record with mandatory provenance fields."""

    bibkey: str
    result_type: str
    number: str
    status: str
    printed_pages: tuple[str, ...]
    pdf_pages: tuple[int, ...]
    source_checksum: str
    application_kind: EvidenceApplicationKind
    hypotheses_used: tuple[str, ...]
    conclusion_used: str
    limitations: tuple[str, ...]

    def __post_init__(self) -> None:
        for field_name in (
            "bibkey",
            "result_type",
            "number",
            "status",
            "source_checksum",
            "conclusion_used",
        ):
            object.__setattr__(self, field_name, _nonempty(getattr(self, field_name), field_name))
        object.__setattr__(
            self,
            "printed_pages",
            _text_tuple(self.printed_pages, "printed_pages"),
        )
        if not self.pdf_pages or not all(
            isinstance(page, int) and not isinstance(page, bool) and page > 0
            for page in self.pdf_pages
        ):
            raise ValueError("pdf_pages must contain positive integers.")
        if not isinstance(self.application_kind, EvidenceApplicationKind):
            raise TypeError("application_kind must be EvidenceApplicationKind.")
        object.__setattr__(
            self,
            "hypotheses_used",
            _text_tuple(self.hypotheses_used, "hypotheses_used"),
        )
        object.__setattr__(
            self,
            "limitations",
            _text_tuple(self.limitations, "limitations"),
        )

    @property
    def citation_key(self) -> str:
        return f"{self.bibkey}:{self.number}"


@dataclass(frozen=True)
class ClosureInterfaceRecord:
    """One contiguous interface or grouped contiguous subword."""

    name: str
    start_position: int
    end_position: int
    factors: tuple[OperatorAtom, ...]
    candidate_results: tuple[str, ...]
    required_hypotheses: tuple[str, ...]
    verified_hypotheses: tuple[str, ...]
    missing_hypotheses: tuple[str, ...]
    status: ClosureInterfaceStatus

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _nonempty(self.name, "name"))
        if self.start_position <= 0 or self.end_position < self.start_position:
            raise ValueError("interface positions must define a positive interval.")
        if len(self.factors) != self.end_position - self.start_position + 1:
            raise ValueError("interface positions and factors disagree.")
        if not all(isinstance(factor, OperatorAtom) for factor in self.factors):
            raise TypeError("interface factors must be OperatorAtom objects.")
        for field_name in (
            "candidate_results",
            "required_hypotheses",
            "verified_hypotheses",
            "missing_hypotheses",
        ):
            object.__setattr__(
                self,
                field_name,
                _text_tuple(getattr(self, field_name), field_name),
            )
        if not isinstance(self.status, ClosureInterfaceStatus):
            raise TypeError("status must be ClosureInterfaceStatus.")


@dataclass(frozen=True)
class ClosureResultInterfaceRow:
    """One global result-by-interface assessment."""

    interface: str
    candidate_result: str
    required_hypotheses: tuple[str, ...]
    verified_hypotheses: tuple[str, ...]
    missing_hypotheses: tuple[str, ...]
    status: ClosureInterfaceStatus

    def __post_init__(self) -> None:
        object.__setattr__(self, "interface", _nonempty(self.interface, "interface"))
        object.__setattr__(
            self,
            "candidate_result",
            _nonempty(self.candidate_result, "candidate_result"),
        )
        for field_name in (
            "required_hypotheses",
            "verified_hypotheses",
            "missing_hypotheses",
        ):
            object.__setattr__(
                self,
                field_name,
                _text_tuple(getattr(self, field_name), field_name),
            )
        if not isinstance(self.status, ClosureInterfaceStatus):
            raise TypeError("status must be ClosureInterfaceStatus.")


@dataclass(frozen=True)
class ClosureWordRecord:
    """Typed metadata over one existing noncommutative Product."""

    identifier: str
    latex_label: str
    product: Product
    signature: tuple[FactorSignature, ...]
    domain: str
    codomain: str
    source_branch: int
    target_branch: int
    current_relation: str
    interfaces: tuple[ClosureInterfaceRecord, ...]
    exact_relations: tuple[str, ...]
    mod_compact_relations: tuple[str, ...]
    candidate_rules: tuple[str, ...]
    obligation_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        for field_name in (
            "identifier",
            "latex_label",
            "domain",
            "codomain",
            "current_relation",
        ):
            object.__setattr__(self, field_name, _nonempty(getattr(self, field_name), field_name))
        if not isinstance(self.product, Product):
            raise TypeError("product must be a Product from the operator AST.")
        if tuple(item.factor for item in self.signature) != self.product.factors:
            raise ClosureAnalysisError("signature factors differ from the AST word.")
        if tuple(item.position for item in self.signature) != tuple(
            range(1, len(self.product.factors) + 1)
        ):
            raise ClosureAnalysisError("signature positions are not consecutive.")
        if self.product.factors.count(R1) != 1:
            raise ClosureAnalysisError("each closure word must contain one atomic R1.")
        if self.product.factors[-1] is not Wplus_12:
            raise ClosureAnalysisError("Wplus_12 must remain the rightmost factor.")
        if self.source_branch != 2 or self.target_branch != 1:
            raise ClosureAnalysisError("closure words must map branch 2 to branch 1.")
        for interface in self.interfaces:
            if self.product.factors[
                interface.start_position - 1 : interface.end_position
            ] != interface.factors:
                raise ClosureAnalysisError("an interface is not a contiguous AST slice.")
        for field_name in (
            "exact_relations",
            "mod_compact_relations",
            "candidate_rules",
            "obligation_ids",
        ):
            object.__setattr__(
                self,
                field_name,
                _text_tuple(getattr(self, field_name), field_name),
            )


@dataclass(frozen=True)
class ClosureProofObligation:
    """One node in the declarative proof-dependency graph."""

    identifier: str
    statement: str
    depends_on: tuple[str, ...]
    sources: tuple[str, ...]
    factors_affected: tuple[OperatorAtom, ...]
    status: ClosureObligationStatus
    evidence_required: str
    blocking_effect: str
    closure_criterion: str

    def __post_init__(self) -> None:
        for field_name in (
            "identifier",
            "statement",
            "evidence_required",
            "blocking_effect",
            "closure_criterion",
        ):
            object.__setattr__(self, field_name, _nonempty(getattr(self, field_name), field_name))
        if not self.identifier.startswith("P-"):
            raise ValueError("obligation identifiers must start with P-.")
        if self.depends_on:
            object.__setattr__(
                self,
                "depends_on",
                tuple(_nonempty(item, "depends_on item") for item in self.depends_on),
            )
        if self.sources:
            object.__setattr__(
                self,
                "sources",
                tuple(_nonempty(item, "sources item") for item in self.sources),
            )
        if not self.factors_affected or not all(
            isinstance(factor, OperatorAtom) for factor in self.factors_affected
        ):
            raise ValueError("factors_affected must contain OperatorAtom objects.")
        if not isinstance(self.status, ClosureObligationStatus):
            raise TypeError("status must be ClosureObligationStatus.")


@dataclass(frozen=True)
class ClosureLemmaCandidate:
    """An explicitly unproved H1, H2, or H3 candidate."""

    identifier: MinimalLemmaChoice
    statement: str
    hypotheses: tuple[str, ...]
    sources: tuple[str, ...]
    established_parts: tuple[str, ...]
    formal_parts: tuple[str, ...]
    open_parts: tuple[str, ...]
    affected_term_count: int
    risk: str
    status: ClosureObligationStatus

    def __post_init__(self) -> None:
        if not isinstance(self.identifier, MinimalLemmaChoice):
            raise TypeError("identifier must be MinimalLemmaChoice.")
        if self.identifier is MinimalLemmaChoice.NONE:
            raise ValueError("NONE is a decision, not a lemma candidate.")
        object.__setattr__(self, "statement", _nonempty(self.statement, "statement"))
        for field_name in (
            "hypotheses",
            "sources",
            "established_parts",
            "formal_parts",
            "open_parts",
        ):
            object.__setattr__(
                self,
                field_name,
                _text_tuple(getattr(self, field_name), field_name),
            )
        if self.affected_term_count != 16:
            raise ClosureAnalysisError("every candidate must be evaluated on all 16 terms.")
        object.__setattr__(self, "risk", _nonempty(self.risk, "risk"))
        if not isinstance(self.status, ClosureObligationStatus):
            raise TypeError("status must be ClosureObligationStatus.")
        if self.status is ClosureObligationStatus.ANALYTICALLY_PROVED:
            raise ClosureAnalysisError("Phase P candidates cannot be marked proved.")


@dataclass(frozen=True)
class MinimalClosureDecision:
    decision: MinimalLemmaChoice
    confidence: DecisionConfidence
    blocking_obligations: tuple[str, ...]
    evidence: tuple[str, ...]
    rationale: str
    prerequisite_statement: str

    def __post_init__(self) -> None:
        if not isinstance(self.decision, MinimalLemmaChoice):
            raise TypeError("decision must be MinimalLemmaChoice.")
        if not isinstance(self.confidence, DecisionConfidence):
            raise TypeError("confidence must be DecisionConfidence.")
        object.__setattr__(
            self,
            "blocking_obligations",
            _text_tuple(self.blocking_obligations, "blocking_obligations"),
        )
        object.__setattr__(self, "evidence", _text_tuple(self.evidence, "evidence"))
        object.__setattr__(self, "rationale", _nonempty(self.rationale, "rationale"))
        object.__setattr__(
            self,
            "prerequisite_statement",
            _nonempty(self.prerequisite_statement, "prerequisite_statement"),
        )


@dataclass(frozen=True)
class FirstPivotClosureAnalysis:
    """Declarative Phase P output with later obligation-status updates."""

    evidence: tuple[EvidenceReference, ...]
    words: tuple[ClosureWordRecord, ...]
    interface_matrix: tuple[ClosureResultInterfaceRow, ...]
    obligations: tuple[ClosureProofObligation, ...]
    candidates: tuple[ClosureLemmaCandidate, ...]
    decision: MinimalClosureDecision

    def __post_init__(self) -> None:
        if len(self.words) != 4:
            raise ClosureAnalysisError("exactly four closure words are required.")
        if tuple(item.identifier for item in self.words) != (
            "K++",
            "K+Ghat",
            "K-+",
            "K-Ghat",
        ):
            raise ClosureAnalysisError("closure-word identifiers are not stable.")
        if len(self.interface_matrix) != 9:
            raise ClosureAnalysisError("the global interface matrix must have nine rows.")
        evidence_keys = tuple(item.citation_key for item in self.evidence)
        if len(set(evidence_keys)) != len(evidence_keys):
            raise ClosureAnalysisError("evidence citation keys must be unique.")
        obligation_ids = tuple(item.identifier for item in self.obligations)
        if obligation_ids != tuple(f"P-{index:02d}" for index in range(1, 13)):
            raise ClosureAnalysisError("obligations must have stable IDs P-01 through P-12.")
        known = set(obligation_ids)
        for index, obligation in enumerate(self.obligations):
            earlier = set(obligation_ids[:index])
            if not set(obligation.depends_on) <= earlier:
                raise ClosureAnalysisError("proof graph must be acyclic and ordered.")
        if not set(self.decision.blocking_obligations) <= known:
            raise ClosureAnalysisError("decision refers to unknown obligations.")
        if tuple(item.identifier for item in self.candidates) != (
            MinimalLemmaChoice.H1,
            MinimalLemmaChoice.H2,
            MinimalLemmaChoice.H3,
        ):
            raise ClosureAnalysisError("candidates must retain H1/H2/H3 order.")
        proved = {
            item.identifier
            for item in self.obligations
            if item.status is ClosureObligationStatus.ANALYTICALLY_PROVED
        }
        if not proved <= {"P-02", "P-03"}:
            raise ClosureAnalysisError("only independently discharged nodes may be proved.")


def _evidence(
    bibkey: str,
    result_type: str,
    number: str,
    printed_pages: tuple[str, ...],
    pdf_pages: tuple[int, ...],
    checksum: str,
    application_kind: EvidenceApplicationKind,
    hypotheses: tuple[str, ...],
    conclusion: str,
    limitations: tuple[str, ...],
) -> EvidenceReference:
    return EvidenceReference(
        bibkey=bibkey,
        result_type=result_type,
        number=number,
        status="verified",
        printed_pages=printed_pages,
        pdf_pages=pdf_pages,
        source_checksum=checksum,
        application_kind=application_kind,
        hypotheses_used=hypotheses,
        conclusion_used=conclusion,
        limitations=limitations,
    )


def closure_evidence_manifest() -> tuple[EvidenceReference, ...]:
    """Return every verified record consumed in the integrated Phase P audit."""

    k17 = "9f92e37599da89f88cbdcdb85d70ef69645aecdc987013f2b54fb878915e9091"
    k14 = "3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec"
    k09 = "b37ae39fe0e49d5efd67341382b8d1a6df44b05579cffae7a18dd312ba485c99"
    reg = "7ed69d096f9e31983cb9b4701e8b486fd15736d82163144f27f9bae244ce41f6"
    cusp = "301f8435688416e6ccd88be1a44f85e12372c9a557bdc1dbea301318b846be23"
    direct = EvidenceApplicationKind.DIRECT
    analogy = EvidenceApplicationKind.ANALOGY
    adaptation = EvidenceApplicationKind.UNPROVED_ADAPTATION
    return (
        _evidence(
            "Karlovich2017HasemanSO",
            "theorem",
            "3.5",
            ("469",),
            (7, 8, 9),
            k17,
            direct,
            ("matrix symbol entries in E-tilde", "1<p<infinity"),
            "Fredholm criterion for an already admissible matrix Mellin PDO.",
            ("Does not prove the correction is a Mellin PDO.",),
        ),
        _evidence(
            "Karlovich2017HasemanSO",
            "theorem",
            "5.3",
            ("476",),
            (14, 15),
            k17,
            analogy,
            ("logarithmic star", "branchwise Mellin transport"),
            "Exact Mellin matrix representation of the Cauchy operator.",
            ("Pure-cusp geometry and localized blocks differ.",),
        ),
        _evidence(
            "Karlovich2017HasemanSO",
            "theorem",
            "5.6",
            ("479",),
            (17, 18),
            k17,
            analogy,
            ("logarithmic star", "slowly oscillating branch shift"),
            "Exact off-diagonal Mellin-PDO representation with shift.",
            ("Does not identify the paper's Wplus_12.",),
        ),
        _evidence(
            "Karlovich2017HasemanSO",
            "theorem",
            "6.2",
            ("482",),
            (20, 21, 22, 23, 24, 25),
            k17,
            analogy,
            ("global logarithmic-star auxiliaries", "matrix E-tilde symbols"),
            "Auxiliary shifts are Fredholm with index zero.",
            ("Paper auxiliaries use an independently proved weighted model.",),
        ),
        _evidence(
            "Karlovich2017HasemanSO",
            "theorem",
            "7.1",
            ("487-488",),
            (25, 26, 27),
            k17,
            analogy,
            ("Theorem 5.6 hypotheses", "coefficient in SO"),
            "Fredholm criterion for the logarithmic-star Haseman operator.",
            ("No direct conclusion for the pure-cusp Schur pivot.",),
        ),
        _evidence(
            "KKL2014TwoShifts",
            "theorem",
            "1.1",
            ("936",),
            (2, 20),
            k14,
            adaptation,
            ("unweighted Lp", "two SOS shifts"),
            "Two-shift auxiliaries are Fredholm with index zero.",
            ("Weighted pure-cusp transfer is unproved.",),
        ),
        _evidence(
            "KKL2014TwoShifts",
            "lemma",
            "2.7",
            ("940",),
            (6,),
            k14,
            adaptation,
            ("unweighted Lp", "Cauchy algebra generated by I and S"),
            "SOS shifts commute with the Cauchy algebra modulo compacts.",
            ("A weighted conjugated analogue is not established.",),
        ),
        _evidence(
            "KKL2014TwoShifts",
            "theorem",
            "3.3",
            ("942",),
            (8,),
            k14,
            adaptation,
            ("both symbols in E(R+,V(R))",),
            "Op(a)Op(b) is Op(ab) modulo compacts.",
            ("U1's multiplier is not in V when gamma1 is nontrivial.",),
        ),
        _evidence(
            "KKL2014TwoShifts",
            "lemma",
            "3.4",
            ("942",),
            (8,),
            k14,
            adaptation,
            ("a,b,c in E", "a radial-only", "c covariable-only"),
            "The separated Mellin triple product is exact.",
            ("The complete right word has not been put in this form.",),
        ),
        _evidence(
            "KKL2014TwoShifts",
            "theorem",
            "3.6",
            ("943",),
            (9,),
            k14,
            adaptation,
            ("scalar symbol in E-tilde", "boundary nonvanishing"),
            "Scalar Mellin-PDO Fredholm criterion and index formula.",
            ("Membership and nonvanishing for the correction are open.",),
        ),
        _evidence(
            "KKL2014TwoShifts",
            "lemma",
            "4.4",
            ("947", "948"),
            (13, 14),
            k14,
            adaptation,
            ("SOS shift", "specific fixed-singularity operator R_y"),
            "U_alpha R_y has an exact Mellin-PDO representation.",
            ("R_y has not been identified with Wplus_12.",),
        ),
        _evidence(
            "KKL2014TwoShifts",
            "lemma",
            "4.5",
            ("949", "950"),
            (15, 16),
            k14,
            adaptation,
            ("SOS shift", "specific fixed-singularity operator R_y"),
            "U_alpha^{+/-1}R_y has a simplified Mellin model modulo compacts.",
            ("Kernel, cutoff, weight, and branch identification remain open.",),
        ),
        _evidence(
            "KKL2014TwoShifts",
            "theorem",
            "5.2",
            ("952", "953", "954"),
            (18, 19, 20),
            k14,
            adaptation,
            ("unweighted Lp", "SOS shift", "half-line Cauchy factors"),
            "One-shift auxiliaries are Fredholm with index zero.",
            ("The paper already proves weighted auxiliary invertibility separately.",),
        ),
        _evidence(
            "Karlovich2009MellinPDO",
            "theorem",
            "4.3",
            ("88",),
            (8,),
            k09,
            adaptation,
            ("scalar operator in A_p", "invertible Fredholm symbol"),
            "Scalar Mellin-algebra Fredholm criterion and index formula.",
            ("Matrix dependency not independently verified; membership is open.",),
        ),
        _evidence(
            "KKL2014Regularization",
            "theorem",
            "5.8",
            ("206",),
            (18,),
            reg,
            adaptation,
            ("symbol in E-tilde", "all boundary and fiber values nonzero"),
            "Regularizers have form Op(b)+K in the same Mellin class.",
            ("Only R1 itself is supplied; lateral closure is not proved.",),
        ),
        _evidence(
            "Karlovich2025Cusps",
            "lemma",
            "3.1",
            ("Article 102, pp. 8-10",),
            (8, 9, 10),
            cusp,
            adaptation,
            ("pure-cusp model", "admissible power weight"),
            "Explicit transported off-diagonal remainders are compact.",
            ("Paper-specific block specialization remains to prove.",),
        ),
        _evidence(
            "Karlovich2025Cusps",
            "lemma",
            "3.2",
            ("Article 102, pp. 11-12",),
            (11, 12),
            cusp,
            adaptation,
            ("pure-cusp model", "conjugated block cases"),
            "Explicit conjugated transported remainders are compact.",
            ("Does not include R1.",),
        ),
        _evidence(
            "Karlovich2025Cusps",
            "corollary",
            "3.3",
            ("Article 102, pp. 12-13",),
            (12, 13),
            cusp,
            adaptation,
            ("cusp-star curve and weight", "bounded Cauchy operator"),
            "Transported Cauchy matrices equal explicit models modulo compacts.",
            ("Does not establish closure of the Schur correction.",),
        ),
        _evidence(
            "Karlovich2025Cusps",
            "theorem",
            "5.1",
            ("Article 102, p. 19",),
            (18, 19),
            cusp,
            adaptation,
            ("operator in the defined weighted half-line algebra",),
            "Fredholm criterion in the half-line mixed algebra.",
            ("Membership of R1 and the complete word is unproved.",),
        ),
        _evidence(
            "Karlovich2025Cusps",
            "theorem",
            "6.1",
            ("Article 102, p. 24",),
            (20, 21, 22, 23, 24),
            cusp,
            adaptation,
            ("transported generators in the source algebra",),
            "Explicit Fredholm symbols of the cusp-algebra generators.",
            ("The paper's R1 is not identified as a generator.",),
        ),
        _evidence(
            "Karlovich2025Cusps",
            "theorem",
            "6.2",
            ("Article 102, p. 25",),
            (25,),
            cusp,
            adaptation,
            ("operator belongs to cusp algebra A",),
            "The symbol map is a homomorphism and supplies a Fredholm criterion.",
            ("Algebra membership must precede use; R1 inclusion is not shown.",),
        ),
    )


def _interface(
    name: str,
    start: int,
    factors: tuple[OperatorAtom, ...],
    candidates: tuple[str, ...],
    required: tuple[str, ...],
    verified: tuple[str, ...],
    missing: tuple[str, ...],
    status: ClosureInterfaceStatus,
) -> ClosureInterfaceRecord:
    return ClosureInterfaceRecord(
        name=name,
        start_position=start,
        end_position=start + len(factors) - 1,
        factors=factors,
        candidate_results=candidates,
        required_hypotheses=required,
        verified_hypotheses=verified,
        missing_hypotheses=missing,
        status=status,
    )


def _word(
    identifier: str,
    latex_label: str,
    factors: tuple[OperatorAtom, ...],
) -> ClosureWordRecord:
    product = Product(factors)
    r1_index = factors.index(R1)
    r1_position = r1_index + 1
    b1 = factors[r1_index + 1]
    q_factors = factors[:r1_index]
    q_start = 1
    interfaces = (
        _interface(
            "Q_1 × R_1",
            q_start,
            q_factors + (R1,),
            ("KKL2014TwoShifts:2.7", "KKL2014TwoShifts:5.2"),
            ("weighted conjugated commutation or direct ordered product rule",),
            ("P^pm and U1^{-1}P^+ are exact paper factors",),
            ("no applicable weighted rule for Q_1 before R1",),
            ClosureInterfaceStatus.BLOCKED,
        ),
        _interface(
            f"R_1 × {b1.name}",
            r1_position,
            (R1, b1),
            ("KKL2014TwoShifts:3.3", "KKL2014TwoShifts:3.4"),
            ("both ordered factors represented in one admissible Mellin calculus",),
            ("R1 has an E-tilde symbol",),
            (
                "U1 multiplier is outside V when nontrivial"
                if b1 is U1
                else "paper-to-source specialization of Ghat1 multiplication",
            ),
            ClosureInterfaceStatus.SPECIALIZATION_TO_PROVE,
        ),
        _interface(
            f"{b1.name} × W_{{1,2}}^+",
            r1_position + 1,
            (b1, Wplus_12),
            (
                "Karlovich2017HasemanSO:5.6",
                "KKL2014TwoShifts:4.4",
                "Karlovich2025Cusps:3.3",
            ),
            ("kernel/symbol identification and compatible domain, weight, cutoff",),
            ("paper gives Wplus_12 exactly as a localized Wiener--Hopf block",),
            ("no Mellin-PDO or cusp-generator identification",),
            ClosureInterfaceStatus.BLOCKED,
        ),
    )
    return ClosureWordRecord(
        identifier=identifier,
        latex_label=latex_label,
        product=product,
        signature=tuple(
            build_factor_signature(position, factor)
            for position, factor in enumerate(factors, start=1)
        ),
        domain="branch 2 copy of L^p(R_+,w_delta)",
        codomain="branch 1 copy of L^p(R_+,w_delta)",
        source_branch=2,
        target_branch=1,
        current_relation="formal ordered word inside the Phase O mod-compact model",
        interfaces=interfaces,
        exact_relations=(
            "U1^{-1}P^+ and P^- come from the exact T_{1,-} expansion",
            "R1 remains one atomic regularizer factor",
        ),
        mod_compact_relations=(
            "Wplus_12 originates in the supplied off-diagonal mod-compact model",
        ),
        candidate_rules=("H1", "H2", "H3"),
        obligation_ids=("P-01", "P-02", "P-03", "P-04", "P-05", "P-06", "P-07", "P-08"),
    )


def closure_words() -> tuple[ClosureWordRecord, ...]:
    """Build the four ordered right-core words without simplification."""

    return (
        _word("K++", r"K_{++}", (U1_inverse, P_plus, R1, U1, Wplus_12)),
        _word(
            "K+Ghat",
            r"K_{+\widehat G}",
            (U1_inverse, P_plus, R1, Ghat1, Wplus_12),
        ),
        _word("K-+", r"K_{-+}", (P_minus, R1, U1, Wplus_12)),
        _word("K-Ghat", r"K_{-\widehat G}", (P_minus, R1, Ghat1, Wplus_12)),
    )


def closure_result_interface_matrix() -> tuple[ClosureResultInterfaceRow, ...]:
    """Assess the nine interfaces required by the integrated Phase P audit."""

    def row(
        interface: str,
        candidate: str,
        required: tuple[str, ...],
        verified: tuple[str, ...],
        missing: tuple[str, ...],
        status: ClosureInterfaceStatus,
    ) -> ClosureResultInterfaceRow:
        return ClosureResultInterfaceRow(
            interface,
            candidate,
            required,
            verified,
            missing,
            status,
        )

    return (
        row(
            r"P^\pm \times R_1",
            "KKL2014TwoShifts Lemma 2.7",
            ("weighted conjugated commutation in the paper's space",),
            ("P^pm are exact Mellin convolution factors",),
            ("source theorem is unweighted and concerns the Cauchy algebra",),
            ClosureInterfaceStatus.NO_RULE,
        ),
        row(
            r"U_1^{-1}P^+ \times R_1",
            "KKL2014TwoShifts Lemma 2.7 and Theorem 5.2",
            ("weighted ordered rule including the exterior dilation",),
            ("the factor order is represented exactly",),
            ("no source result covers this weighted ordered triple",),
            ClosureInterfaceStatus.BLOCKED,
        ),
        row(
            r"R_1 \times U_1",
            "exact ordered factorized pair (r1,d_gamma1)",
            ("paper Mellin quantization and normalized dilation definitions",),
            ("direct dense-core derivation and bounded extension",),
            ("membership and closure of the provisional factorized class",),
            ClosureInterfaceStatus.CERTIFIED_EXACT,
        ),
        row(
            r"R_1 \times \widehat G_1",
            "KKL2014TwoShifts Theorem 3.3",
            ("r1 and the lambda-independent Ghat1 symbol belong to E",),
            ("both symbols are in E-tilde after the paper's Phi_delta conjugation",),
            ("no exact identity is asserted; the theorem gives a compact residue",),
            ClosureInterfaceStatus.CERTIFIED_MOD_COMPACT,
        ),
        row(
            r"U_1 \times W_{1,2}^{+}",
            "KKL2014TwoShifts Lemmas 4.4--4.5",
            ("Wplus_12 is the source operator R_y after all transports",),
            ("both kernels have fixed-singularity structure",),
            ("kernels, cutoff, weight, and branch maps are different",),
            ClosureInterfaceStatus.ANALOGY_ONLY,
        ),
        row(
            r"\widehat G_1 \times W_{1,2}^{+}",
            "Karlovich2017HasemanSO Theorem 5.6",
            ("log-star block theorem specializes to the localized cusp block",),
            ("the coefficient and block order are explicit",),
            ("geometry, localization, and symbol identification",),
            ClosureInterfaceStatus.ANALOGY_ONLY,
        ),
        row(
            r"R_1B_1 \times W_{1,2}^{+}",
            "H2 via Mellin semiproduct or cusp-algebra closure",
            ("H1 and a compatible class for Wplus_12",),
            ("the suffix is contiguous in every word",),
            ("H1 and Wplus_12 class membership",),
            ClosureInterfaceStatus.BLOCKED,
        ),
        row(
            r"Q_1 \times R_1B_1W_{1,2}^{+}",
            "H3 via a weighted ordered exterior rule",
            ("H2 and a rule for each Q1",),
            ("the full right core is represented without reordering",),
            ("H2 and the weighted Q1 interface",),
            ClosureInterfaceStatus.BLOCKED,
        ),
        row(
            "complete word",
            "Mellin-PDO membership or Karlovich2025 cusp-algebra membership",
            ("all interfaces close in one verified class; residues are compact",),
            ("the AST fixes all four products and their domains",),
            ("P-01 through P-08",),
            ClosureInterfaceStatus.BLOCKED,
        ),
    )


def closure_proof_obligations() -> tuple[ClosureProofObligation, ...]:
    """Return the ordered acyclic graph P-01 through P-12."""

    def node(
        identifier: str,
        statement: str,
        depends_on: tuple[str, ...],
        sources: tuple[str, ...],
        factors: tuple[OperatorAtom, ...],
        status: ClosureObligationStatus,
        evidence_required: str,
        blocking_effect: str,
        closure_criterion: str,
    ) -> ClosureProofObligation:
        return ClosureProofObligation(
            identifier,
            statement,
            depends_on,
            sources,
            factors,
            status,
            evidence_required,
            blocking_effect,
            closure_criterion,
        )

    return (
        node(
            "P-01",
            "Identify Wplus_12 in a closure-compatible operator class.",
            (),
            ("Karlovich2025Cusps:3.3", "KKL2014TwoShifts:4.4"),
            (Wplus_12,),
            ClosureObligationStatus.BLOCKED,
            "kernel or symbol comparison including cutoff, weight, and branches",
            "blocks H2, H3, Mellin closure, and cusp-algebra membership",
            "an exact or certified mod-compact class identification",
        ),
        node(
            "P-02",
            "Represent U1 and prove the exact ordered factorized composition R1 U1.",
            (),
            (
                "paper:mellin-representation-dilation-shift",
                "phase-q:exact-right-dilation-definition-proof",
            ),
            (R1, U1),
            ClosureObligationStatus.ANALYTICALLY_PROVED,
            "direct dense-core derivation plus bounded extension",
            "removes the composition-formula blocker, not standard-class membership",
            "retain the exact ordered pair (r1,d_gamma1) without an E-tilde claim",
        ),
        node(
            "P-03",
            "Represent Ghat1 and certify the ordered R1 Ghat1 semiproduct.",
            (),
            ("paper:normalized-diagonal-symbol-admissible", "KKL2014TwoShifts:3.3"),
            (R1, Ghat1),
            ClosureObligationStatus.ANALYTICALLY_PROVED,
            "paper membership proof plus every hypothesis of Theorem 3.3",
            "discharges the coefficient case only, modulo compact operators",
            "Op(r1)M_Ghat1 is Op(r1 Ghat1) modulo compact operators",
        ),
        node(
            "P-04",
            "Place the separately controlled R1 U1 and R1 Ghat1 cases in one admissible H1 framework.",
            ("P-02", "P-03"),
            ("KKL2014TwoShifts:3.3", "KKL2014TwoShifts:3.4"),
            (R1, U1, Ghat1),
            ClosureObligationStatus.BLOCKED,
            "U1 case: provisional-class closure; Ghat1 case: already certified mod compact",
            "blocks uniform H1 and every longer candidate",
            "one named admissible framework covering both relations at their proved strengths",
        ),
        node(
            "P-05",
            "Control (R1 B1)Wplus_12 without reordering factors.",
            ("P-01", "P-04"),
            ("Karlovich2017HasemanSO:5.6", "Karlovich2025Cusps:3.3"),
            (R1, U1, Ghat1, Wplus_12),
            ClosureObligationStatus.BLOCKED,
            "compatible representation of Wplus_12 and a lateral product theorem",
            "blocks H2 and H3",
            "derive a certified model for the ordered three-factor suffix",
        ),
        node(
            "P-06",
            "Control Q1 before R1 B1 Wplus_12 without commuting it.",
            ("P-05",),
            ("KKL2014TwoShifts:2.7", "KKL2014TwoShifts:5.2"),
            (U1_inverse, P_plus, P_minus, R1),
            ClosureObligationStatus.BLOCKED,
            "weighted ordered rule for both Q1 choices",
            "blocks H3",
            "preserve the stored order and produce a controlled remainder",
        ),
        node(
            "P-07",
            "Prove compactness of every accumulated residue.",
            ("P-04", "P-05", "P-06"),
            ("Karlovich2025Cusps:3.1", "KKL2014TwoShifts:3.3"),
            (R1, Wplus_12),
            ClosureObligationStatus.BLOCKED,
            "explicit residual formulas and compactness estimates",
            "prevents any candidate from becoming certified mod compact",
            "all residuals lie in K(L^p(R_+,w_delta)) with branch typing",
        ),
        node(
            "P-08",
            "Identify the resulting symbol and its precise class.",
            ("P-04", "P-05", "P-06", "P-07"),
            ("Karlovich2009MellinPDO:4.3", "Karlovich2025Cusps:6.1"),
            (R1, U1, Ghat1, Wplus_12),
            ClosureObligationStatus.BLOCKED,
            "constructed symbol, variable dependence, and class estimates",
            "blocks algebra membership and later symbol tests",
            "one named symbol in a verified Mellin or cusp quotient class",
        ),
        node(
            "P-09",
            "Propagate the closure result to all 16 Phase O terms.",
            ("P-06", "P-07", "P-08"),
            ("phase-o:C2-T01--C2-T16",),
            (R1, Wplus_12),
            ClosureObligationStatus.BLOCKED,
            "termwise branch-compatible substitution preserving order",
            "blocks reconstruction of the complete correction",
            "all 16 terms receive a relation of identical certified strength",
        ),
        node(
            "P-10",
            "Reconstruct the complete correction C2 from the 16 controlled terms.",
            ("P-09",),
            ("phase-o:complete-correction-expansion",),
            (R1, Wplus_12),
            ClosureObligationStatus.BLOCKED,
            "finite summation with compatible compact ideals",
            "blocks any membership statement for the corrected pivot",
            "recover C2 without changing the Phase N mod-compact relation",
        ),
        node(
            "P-11",
            "Prove membership of N2first in the selected operator class or algebra.",
            ("P-10",),
            ("Karlovich2009MellinPDO:4.3", "Karlovich2025Cusps:6.2"),
            (R1, Wplus_12),
            ClosureObligationStatus.BLOCKED,
            "closure of N2 plus the reconstructed correction",
            "blocks every analytic criterion for the pivot",
            "a proved membership theorem with an identified quotient symbol",
        ),
        node(
            "P-12",
            "Only after membership, verify nonvanishing and Fredholmness conditions.",
            ("P-11",),
            ("KKL2014Regularization:5.8", "Karlovich2025Cusps:6.2"),
            (R1, Wplus_12),
            ClosureObligationStatus.BLOCKED,
            "complete boundary/fiber or determinant nonvanishing proof",
            "outside Phase P and blocks any Fredholm conclusion",
            "all hypotheses of the chosen criterion are checked",
        ),
    )


def closure_lemma_candidates() -> tuple[ClosureLemmaCandidate, ...]:
    """Build H1/H2/H3 as blocked, explicitly unproved candidates."""

    return (
        ClosureLemmaCandidate(
            MinimalLemmaChoice.H1,
            "For both B1 choices, R1 B1 = Op(c_B1)+K_B1 in one explicit admissible class.",
            (
                "R1 has its supplied E-tilde symbol",
                "each B1 is represented in the same or a proved extended calculus",
            ),
            ("KKL2014TwoShifts:3.3", "KKL2014TwoShifts:3.4"),
            (
                "R1 is a Mellin PDO",
                "R1 U1 has an exact ordered factorized representation",
                "R1 Ghat1 has a certified mod-compact semiproduct",
            ),
            ("the ordered products occur unchanged in all 16 terms",),
            (
                "U1 has multiplier gamma1^(i lambda), outside V when nontrivial",
                "the factorized U1 output has no proved closed calculus",
                "no uniform admissible framework for both B1 choices is identified",
            ),
            16,
            "A false E-class assertion for the U1 case would invalidate the lemma.",
            ClosureObligationStatus.BLOCKED,
        ),
        ClosureLemmaCandidate(
            MinimalLemmaChoice.H2,
            "For both B1 choices, R1 B1 Wplus_12 = Op(c_B1^+)+K_B1^+.",
            (
                "H1 has been discharged",
                "Wplus_12 has a compatible Mellin representation",
            ),
            (
                "Karlovich2017HasemanSO:5.6",
                "KKL2014TwoShifts:4.4",
                "Karlovich2025Cusps:3.3",
            ),
            ("the paper gives the exact localized Wiener--Hopf definition",),
            ("the suffix is contiguous in all 16 terms",),
            (
                "Wplus_12 is not identified with R_y or a Mellin PDO",
                "compact cutoff residues are not controlled in this product",
            ),
            16,
            "Depends on the currently false/unproved Wplus_12 identification.",
            ClosureObligationStatus.BLOCKED,
        ),
        ClosureLemmaCandidate(
            MinimalLemmaChoice.H3,
            "For every Q1,B1 pair, Q1 R1 B1 Wplus_12 = Op(c_Q1,B1^+)+K_Q1,B1^+.",
            (
                "H2 has been discharged",
                "both exterior Q1 factors obey a weighted ordered closure rule",
            ),
            ("KKL2014TwoShifts:2.7", "KKL2014TwoShifts:5.2"),
            ("the four Q1,B1 words are exact AST products",),
            ("the full right core is contiguous in all 16 terms",),
            (
                "H1 and H2 remain blocked",
                "the source commutation theorem is unweighted and cannot reorder here",
            ),
            16,
            "Bundles three independent missing interfaces and risks hiding them.",
            ClosureObligationStatus.BLOCKED,
        ),
    )


def minimal_closure_decision() -> MinimalClosureDecision:
    """Retain NONE until the factorized output and Wplus_12 are controlled."""

    return MinimalClosureDecision(
        decision=MinimalLemmaChoice.NONE,
        confidence=DecisionConfidence.HIGH,
        blocking_obligations=("P-01", "P-04"),
        evidence=(
            "phase-q:exact-right-dilation-definition-proof",
            "phase-q:certified-Ghat1-semiproduct",
            "paper:normalized-wh-blocks-section-six",
            "KKL2014TwoShifts:3.3",
            "KKL2014TwoShifts:4.4",
            "Karlovich2025Cusps:3.3",
        ),
        rationale=(
            "The two H1 cases now have separate certified relations, but the exact "
            "factorized U1 output has no proved closed admissible calculus. H2 and H3 "
            "also require an unproved identification of the localized Wplus_12."
        ),
        prerequisite_statement=(
            "Prove the minimum closure properties needed for the factorized pair "
            "(r1,d_gamma1), without forcing it into E-tilde; independently identify "
            "Wplus_12 before attempting H2."
        ),
    )


def build_first_pivot_closure_analysis() -> FirstPivotClosureAnalysis:
    """Build the deterministic evidence/word/obligation/candidate trace."""

    return FirstPivotClosureAnalysis(
        evidence=closure_evidence_manifest(),
        words=closure_words(),
        interface_matrix=closure_result_interface_matrix(),
        obligations=closure_proof_obligations(),
        candidates=closure_lemma_candidates(),
        decision=minimal_closure_decision(),
    )
