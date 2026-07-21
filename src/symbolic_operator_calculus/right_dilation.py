"""Exact right-dilation metadata for the Phase Q Mellin calculation.

The noncommutative operator expression remains the existing ``Product`` AST.
The types below distinguish a standard Mellin symbol from an oscillatory
right factor without claiming that their pointwise product is standard.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

import sympy as sp

from .operators import Ghat1, R1, U1, Product
from .semantics import (
    AnalyticProofObligation,
    DerivationRelationKind,
    ExactIdentity,
    ExactIdentityScope,
    ModCompactEquivalence,
    OperatorRule,
    RuleCertificationStatus,
    SemanticDerivationStep,
)

if TYPE_CHECKING:
    from .mellin import MellinSymbol


class RightDilationError(ValueError):
    """Raised when Phase Q metadata contradicts its ordered AST product."""


class StandardMellinMembership(str, Enum):
    """Explicit membership status; no status is inferred from a product."""

    CERTIFIED = "CERTIFIED"
    EXCLUDED_NONTRIVIAL_DILATION = "EXCLUDED_NONTRIVIAL_DILATION"
    NOT_DEMONSTRATED = "NOT_DEMONSTRATED"


class DilationFrequencyCase(str, Enum):
    IDENTITY = "IDENTITY"
    NONTRIVIAL = "NONTRIVIAL"


class AnalyticLemmaStatus(str, Enum):
    ANALYTICALLY_PROVED = "ANALYTICALLY_PROVED"


def _texts(value: object, field_name: str) -> tuple[str, ...]:
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
class StandardMellinSymbol:
    """Declarative reference to a symbol with separately certified membership."""

    name: str
    symbol_class: str
    membership: StandardMellinMembership
    evidence: tuple[str, ...]

@dataclass(frozen=True)
class OscillatoryMellinFactor:
    """The typed family d_gamma(lambda)=gamma**(i lambda)."""

    family_symbol: MellinSymbol
    gamma: sp.Symbol
    frequency: sp.Symbol
    case: DilationFrequencyCase
    standard_v_membership: StandardMellinMembership
    proof: tuple[str, ...]

    def __post_init__(self) -> None:
        from .mellin import MellinSymbol

        if not isinstance(self.family_symbol, MellinSymbol):
            raise TypeError("family_symbol must be a MellinSymbol.")
        if (
            self.case is DilationFrequencyCase.NONTRIVIAL
            and self.standard_v_membership
            is not StandardMellinMembership.EXCLUDED_NONTRIVIAL_DILATION
        ):
            raise RightDilationError("a nontrivial dilation must be excluded from V.")
        if (
            self.case is DilationFrequencyCase.IDENTITY
            and self.standard_v_membership is not StandardMellinMembership.CERTIFIED
        ):
            raise RightDilationError("the identity factor must have certified membership.")
        object.__setattr__(self, "proof", _texts(self.proof, "proof"))

    @property
    def expression(self) -> sp.Expr:
        if self.case is DilationFrequencyCase.IDENTITY:
            return sp.Integer(1)
        return self.family_symbol.as_expr()


@dataclass(frozen=True)
class FactorizedOscillatoryMellinSymbol:
    """Ordered pair (a,d_gamma), not a pointwise-product membership claim."""

    standard_factor: StandardMellinSymbol
    oscillatory_factor: OscillatoryMellinFactor
    ast_product: Product
    provisional_class: str
    standard_product_membership: StandardMellinMembership
    boundedness_reason: str
    pending_properties: tuple[AnalyticProofObligation, ...]

    def __post_init__(self) -> None:
        if self.ast_product != Product((R1, U1)):
            raise RightDilationError("the represented AST product must be R1 followed by U1.")
        if self.standard_product_membership is not (
            StandardMellinMembership.NOT_DEMONSTRATED
        ):
            raise RightDilationError("the factorized product is not a certified standard symbol.")
        if not self.pending_properties or not all(
            isinstance(item, AnalyticProofObligation) for item in self.pending_properties
        ):
            raise TypeError("pending_properties must contain proof obligations.")

    @property
    def ordered_symbol_factors(self) -> tuple[object, object]:
        return self.standard_factor, self.oscillatory_factor


@dataclass(frozen=True)
class RightDilationDefinitionProof:
    quantization: str
    normalized_action: str
    dense_domain: str
    substitution: str
    sign_check: str
    continuity_extension: str

@dataclass(frozen=True)
class RightDilationCompositionTrace:
    ast_product: Product
    factorized_symbol: FactorizedOscillatoryMellinSymbol
    exact_identity: ExactIdentity
    rule: OperatorRule
    semantic_steps: tuple[SemanticDerivationStep, ...]

    def __post_init__(self) -> None:
        if self.ast_product != Product((R1, U1)):
            raise RightDilationError("trace order must be R1 then U1.")
        if self.exact_identity.left != self.ast_product:
            raise RightDilationError("exact identity must use the authoritative AST product.")
        if self.exact_identity.right is not self.factorized_symbol:
            raise RightDilationError("exact identity must retain the factorized symbol object.")
        if self.rule.relation is not self.exact_identity:
            raise RightDilationError("the exact rule must reference the exact identity.")


@dataclass(frozen=True)
class CoefficientSemiproductEvidence:
    theorem: str
    printed_page: str
    pdf_page: int
    checksum: str
    paper_hypotheses: tuple[str, ...]

@dataclass(frozen=True)
class CoefficientSemiproductTrace:
    ast_product: Product
    left_symbol: StandardMellinSymbol
    right_symbol: StandardMellinSymbol
    product_symbol: StandardMellinSymbol
    relation: ModCompactEquivalence
    rule: OperatorRule

    def __post_init__(self) -> None:
        if self.ast_product != Product((R1, Ghat1)):
            raise RightDilationError("coefficient order must be R1 then Ghat1.")
        if self.relation.left != self.ast_product or self.relation.right is not self.product_symbol:
            raise RightDilationError("semiproduct relation endpoints are inconsistent.")
        if self.rule.relation is not self.relation:
            raise RightDilationError("semiproduct rule must reference its relation.")


@dataclass(frozen=True)
class DampingVariationLemma:
    name: str
    hypotheses: tuple[str, ...]
    variation_bound: sp.Expr
    conclusion: str
    proof_steps: tuple[str, ...]
    source: str
    status: AnalyticLemmaStatus

def build_oscillatory_mellin_factor(
    lam: sp.Symbol,
    gamma: sp.Symbol,
    *,
    case: DilationFrequencyCase = DilationFrequencyCase.NONTRIVIAL,
) -> OscillatoryMellinFactor:
    """Build d_gamma with an explicit identity/nontrivial/unresolved case."""

    from .mellin import build_dilation_multiplier

    symbol = build_dilation_multiplier(lam, gamma)
    membership = {
        DilationFrequencyCase.IDENTITY: StandardMellinMembership.CERTIFIED,
        DilationFrequencyCase.NONTRIVIAL: (
            StandardMellinMembership.EXCLUDED_NONTRIVIAL_DILATION
        ),
    }[case]
    proof_by_case = {
        DilationFrequencyCase.IDENTITY: "gamma=1 gives d_gamma=1 and zero variation",
        DilationFrequencyCase.NONTRIVIAL: (
            "for gamma>0, gamma!=1, |d_gamma'|=|log gamma| and total variation is infinite"
        ),
    }
    proof = (proof_by_case[case],)
    return OscillatoryMellinFactor(symbol, gamma, lam, case, membership, proof)


def _right_dilation_obligations() -> tuple[AnalyticProofObligation, ...]:
    labels = (
        ("right_gamma_sums", "Control sums carrying distinct dilation frequencies."),
        ("right_gamma_products", "Prove a product law for two factorized elements."),
        ("right_gamma_left_composition", "Analyze composition by a dilation on the left."),
        ("right_gamma_involution", "Define and control involution in the provisional class."),
        ("right_gamma_fredholm_symbol", "Construct a Fredholm symbol only after closure."),
        ("right_gamma_norm_closure", "Determine a norm and prove completeness or closure."),
        ("right_gamma_wiener_hopf", "Analyze interaction with localized Wiener--Hopf factors later."),
    )
    return tuple(AnalyticProofObligation(key, statement) for key, statement in labels)


def build_right_dilation_composition_trace(
    lam: sp.Symbol,
    gamma: sp.Symbol,
) -> RightDilationCompositionTrace:
    """Build the exact R1 U1 identity while retaining the ordered symbol pair."""

    ast_product = Product((R1, U1))
    standard = StandardMellinSymbol(
        "r_1",
        "tilde-E(R_+,V(R))",
        StandardMellinMembership.CERTIFIED,
        ("paper definition of R1", "KKL2014Regularization Theorem 5.8"),
    )
    oscillatory = build_oscillatory_mellin_factor(lam, gamma)
    factorized = FactorizedOscillatoryMellinSymbol(
        standard,
        oscillatory,
        ast_product,
        "tilde-E^(right-gamma)",
        StandardMellinMembership.NOT_DEMONSTRATED,
        "Op(r1) is bounded and the normalized right dilation is an isometry",
        _right_dilation_obligations(),
    )
    proof = RightDilationDefinitionProof(
        "paper equation: Mellin Kohn--Nirenberg kernel with (x/y)^(i lambda)",
        "Phi_delta U_gamma Phi_delta^(-1)g(x)=g(gamma x)",
        "C_c^infinity(R_+) with symmetric frequency truncation",
        "z=gamma y in the inner integral; integrations are not interchanged",
        "the log-Gaussian test gives gamma^(+i lambda)",
        "boundedness of Op(r1) and isometry of the dilation give the unique extension",
    )
    identity = ExactIdentity(
        ast_product,
        factorized,
        evidence=proof,
        scope=ExactIdentityScope.OPERATORIAL,
        hypotheses=("gamma>0", "r1 defines a bounded Mellin PDO"),
    )
    rule = OperatorRule(
        "exact_right_dilation_mellin_composition",
        DerivationRelationKind.EXACT_EQUALITY,
        identity,
        (
            "paper Mellin transform uses x^(-i lambda)",
            "paper quantization uses (x/y)^(i lambda)",
            "U_gamma is normalized by gamma^(delta+1/p)",
            "gamma>0",
        ),
        "direct derivation from the paper definitions; no external composition theorem",
        RuleCertificationStatus.CERTIFIED_EXACT,
    )
    steps = (
        SemanticDerivationStep("ordered_ast_product", DerivationRelationKind.EXACT_EQUALITY, ast_product),
        SemanticDerivationStep("factorized_symbol", DerivationRelationKind.FORMAL_SUBSTITUTION, factorized),
        SemanticDerivationStep("exact_composition", DerivationRelationKind.EXACT_EQUALITY, identity, proof),
        *(SemanticDerivationStep(item.key, item.relation_kind, item) for item in factorized.pending_properties),
    )
    return RightDilationCompositionTrace(ast_product, factorized, identity, rule, steps)


def build_coefficient_semiproduct_trace() -> CoefficientSemiproductTrace:
    """Record the verified Theorem 3.3 specialization for R1 Ghat1."""

    ast_product = Product((R1, Ghat1))
    class_name = "tilde-E(R_+,V(R))"
    r1 = StandardMellinSymbol("r_1", class_name, StandardMellinMembership.CERTIFIED, ("paper Section 6",))
    ghat = StandardMellinSymbol(
        "Ghat_1(x)", class_name, StandardMellinMembership.CERTIFIED,
        ("paper proposition: normalized diagonal symbol admissibility",),
    )
    product = StandardMellinSymbol(
        "r_1 Ghat_1", class_name, StandardMellinMembership.CERTIFIED,
        ("tilde-E is a Banach algebra",),
    )
    evidence = CoefficientSemiproductEvidence(
        "KKL2014TwoShifts Theorem 3.3",
        "942",
        8,
        "3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec",
        (
            "r1 belongs to tilde-E and hence E",
            "Ghat1(x) is lambda-independent and belongs to tilde-E",
            "after Phi_delta both operators act on L^p(R_+,dmu)",
        ),
    )
    relation = ModCompactEquivalence(
        ast_product,
        product,
        space="L^p(R_+,w_delta), equivalently L^p(R_+,dmu) after Phi_delta",
        compact_ideal="compact operators",
        evidence=evidence,
    )
    rule = OperatorRule(
        "regularizer_slowly_oscillating_coefficient_semiproduct",
        DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
        relation,
        evidence.paper_hypotheses,
        "KKL2014TwoShifts, Theorem 3.3, printed p. 942/PDF p. 8",
        RuleCertificationStatus.CERTIFIED_MOD_COMPACT,
    )
    return CoefficientSemiproductTrace(ast_product, r1, ghat, product, relation, rule)


def build_damping_variation_lemma(gamma: sp.Symbol) -> DampingVariationLemma:
    """Return the proved BV damping estimate as explicit analytic metadata."""

    if not isinstance(gamma, sp.Symbol):
        raise TypeError("gamma must be a SymPy Symbol.")
    variation_h = sp.Symbol("Var_h", nonnegative=True)
    l1_h = sp.Symbol("L1_h", nonnegative=True)
    bound = variation_h + sp.Abs(sp.log(gamma)) * l1_h
    return DampingVariationLemma(
        "oscillatory_damping_preserves_bounded_variation",
        ("gamma>0", "h is absolutely continuous", "h in V(R) intersect L1(R)"),
        bound,
        "Var(d_gamma h) <= Var(h)+|log gamma| ||h||_1 and zero endpoint limits are preserved",
        (
            "(d_gamma h)'=d_gamma(h'+i log(gamma)h) almost everywhere",
            "integrate the triangle inequality over R",
            "|d_gamma h|=|h| preserves endpoint convergence to zero",
        ),
        "direct proof from absolute continuity and total variation",
        AnalyticLemmaStatus.ANALYTICALLY_PROVED,
    )
