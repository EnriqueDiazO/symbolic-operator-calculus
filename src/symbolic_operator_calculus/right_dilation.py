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

from .operators import Ghat1, P_minus, P_plus, R1, U1, U1_inverse, Product
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
        if not isinstance(self.ast_product, Product) or not self.ast_product.factors:
            raise RightDilationError("a factorized symbol requires a nonempty AST product.")
        expected = (
            StandardMellinMembership.CERTIFIED
            if self.oscillatory_factor.case is DilationFrequencyCase.IDENTITY
            else StandardMellinMembership.NOT_DEMONSTRATED
        )
        if self.standard_product_membership is not expected:
            raise RightDilationError("standard membership contradicts the frequency case.")
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


@dataclass(frozen=True)
class RadialScalingTrace:
    original_symbol: StandardMellinSymbol
    scaled_symbol: StandardMellinSymbol
    gamma: sp.Symbol
    spatial_variable: sp.Symbol
    radial_argument: sp.Expr
    verified_contracts: tuple[str, ...]
    status: AnalyticLemmaStatus


@dataclass(frozen=True)
class LeftDilationCovarianceTrace:
    source_product: Product
    conjugated_product: Product
    scaling: RadialScalingTrace
    inverse_factor: OscillatoryMellinFactor
    factorized_symbol: FactorizedOscillatoryMellinSymbol
    left_identity: ExactIdentity
    conjugated_identity: ExactIdentity
    rules: tuple[OperatorRule, OperatorRule]
    frequency_cancellation: sp.Expr

    def __post_init__(self) -> None:
        if self.source_product != Product((U1_inverse, R1)):
            raise RightDilationError("left covariance must preserve U1_inverse R1.")
        if self.conjugated_product != Product((U1_inverse, R1, U1)):
            raise RightDilationError("conjugated covariance order is invalid.")
        if self.frequency_cancellation != 1:
            raise RightDilationError("inverse and forward frequencies must cancel.")


@dataclass(frozen=True)
class LeftCoreTrace:
    identifier: str
    ast_product: Product
    output_symbol: StandardMellinSymbol | FactorizedOscillatoryMellinSymbol
    relation: ModCompactEquivalence
    rule: OperatorRule
    output_kind: str
    compact_residual_terms: tuple[str, ...]
    weighted_conjugation: str

    def __post_init__(self) -> None:
        if self.relation.left is not self.ast_product:
            raise RightDilationError("left-core relation must retain its AST product.")
        if self.relation.right is not self.output_symbol:
            raise RightDilationError("left-core output and relation disagree.")
        if self.rule.certification_status is not RuleCertificationStatus.CERTIFIED_MOD_COMPACT:
            raise RightDilationError("every Phase R core is certified modulo compact operators.")


@dataclass(frozen=True)
class FirstPivotLeftCoreClosureTrace:
    covariance: LeftDilationCovarianceTrace
    cauchy_relations: tuple[ModCompactEquivalence, ModCompactEquivalence]
    cauchy_rules: tuple[OperatorRule, OperatorRule]
    cores: tuple[LeftCoreTrace, ...]
    compact_ideal_proof: tuple[str, ...]
    remaining_blocker: str

def build_oscillatory_mellin_factor(
    lam: sp.Symbol,
    gamma: sp.Symbol,
    *,
    case: DilationFrequencyCase = DilationFrequencyCase.NONTRIVIAL,
    inverse: bool = False,
) -> OscillatoryMellinFactor:
    """Build d_gamma with an explicit identity/nontrivial/unresolved case."""

    from .mellin import build_dilation_multiplier

    symbol = build_dilation_multiplier(lam, gamma, inverse=inverse)
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


def _standard_symbol(name: str, *evidence: str) -> StandardMellinSymbol:
    return StandardMellinSymbol(
        name, "tilde-E(R_+,V(R))", StandardMellinMembership.CERTIFIED, evidence
    )


def build_radial_scaling_trace(
    gamma: sp.Symbol,
    spatial_variable: sp.Symbol,
    *,
    symbol_name: str = "a",
) -> RadialScalingTrace:
    """Prove E-tilde invariance under x -> x/gamma for gamma>0."""

    original = _standard_symbol(symbol_name, "supplied E-tilde membership")
    scaled = _standard_symbol(
        f"{symbol_name}_gamma(x,lambda)={symbol_name}(x/gamma,lambda)",
        "direct radial-scaling invariance proof",
    )
    contracts = (
        "the C_b(V) norm is unchanged because x -> x/gamma is bijective",
        "V-valued continuity is preserved by composition",
        "cm_t(a_gamma)=cm_(t/gamma)(a) at zero and infinity",
        "covariable-translation suprema are unchanged",
        "uniform derivative-tail suprema are unchanged",
        "endpoint fibers are preserved by the radial scaling automorphism",
    )
    return RadialScalingTrace(
        original,
        scaled,
        gamma,
        spatial_variable,
        spatial_variable / gamma,
        contracts,
        AnalyticLemmaStatus.ANALYTICALLY_PROVED,
    )


def build_left_dilation_covariance_trace(
    lam: sp.Symbol,
    gamma: sp.Symbol,
    spatial_variable: sp.Symbol,
    *,
    case: DilationFrequencyCase = DilationFrequencyCase.NONTRIVIAL,
) -> LeftDilationCovarianceTrace:
    """Build both exact left-covariance identities from the paper definitions."""

    scaling = build_radial_scaling_trace(gamma, spatial_variable, symbol_name="r_1")
    source = Product((U1_inverse, R1))
    conjugated = Product((U1_inverse, R1, U1))
    inverse_factor = build_oscillatory_mellin_factor(
        lam, gamma, case=case, inverse=True
    )
    membership = (
        StandardMellinMembership.CERTIFIED
        if case is DilationFrequencyCase.IDENTITY
        else StandardMellinMembership.NOT_DEMONSTRATED
    )
    factorized = FactorizedOscillatoryMellinSymbol(
        scaling.scaled_symbol,
        inverse_factor,
        source,
        "tilde-E^(right-gamma-inverse)",
        membership,
        "Op(a_gamma) and D_(gamma^-1) are bounded",
        _right_dilation_obligations(),
    )
    proof = (
        "evaluate Op(a)f at x/gamma on C_c^infinity(R_+)",
        "((x/gamma)/y)^(i lambda)=gamma^(-i lambda)(x/y)^(i lambda)",
        "the log-Gaussian test gives gamma^(-i lambda)",
        "boundedness and the dilation isometry give the unique extension",
    )
    left_identity = ExactIdentity(
        source,
        factorized,
        evidence=proof,
        scope=ExactIdentityScope.OPERATORIAL,
        hypotheses=("gamma>0", "a belongs to E-tilde"),
    )
    conjugated_identity = ExactIdentity(
        conjugated,
        scaling.scaled_symbol,
        evidence=(*proof, "D_(gamma^-1)D_gamma=I exactly"),
        scope=ExactIdentityScope.OPERATORIAL,
        hypotheses=("gamma>0", "a belongs to E-tilde"),
    )
    preconditions = (
        "paper Mellin kernel uses (x/y)^(i lambda)",
        "D_gamma f(x)=f(gamma x)",
        "gamma>0",
    )
    rules = (
        OperatorRule(
            "exact_left_dilation_mellin_covariance",
            DerivationRelationKind.EXACT_EQUALITY,
            left_identity,
            preconditions,
            "direct derivation from the paper definitions",
            RuleCertificationStatus.CERTIFIED_EXACT,
        ),
        OperatorRule(
            "exact_dilation_conjugation_mellin_covariance",
            DerivationRelationKind.EXACT_EQUALITY,
            conjugated_identity,
            preconditions,
            "direct derivation from the paper definitions",
            RuleCertificationStatus.CERTIFIED_EXACT,
        ),
    )
    forward = build_oscillatory_mellin_factor(lam, gamma, case=case)
    cancellation = sp.powsimp(
        inverse_factor.expression * forward.expression, force=True
    )
    return LeftDilationCovarianceTrace(
        source,
        conjugated,
        scaling,
        inverse_factor,
        factorized,
        left_identity,
        conjugated_identity,
        rules,
        cancellation,
    )


def build_first_pivot_left_core_closure(
    lam: sp.Symbol,
    gamma: sp.Symbol,
    spatial_variable: sp.Symbol,
    *,
    case: DilationFrequencyCase = DilationFrequencyCase.NONTRIVIAL,
) -> FirstPivotLeftCoreClosureTrace:
    """Reduce the four ordered cores before Wplus_12, never including that factor."""

    covariance = build_left_dilation_covariance_trace(
        lam, gamma, spatial_variable, case=case
    )
    source_checksum = "3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec"
    cauchy_symbols = (
        _standard_symbol("p^+ r_1", "p^+ and r_1 belong to E-tilde"),
        _standard_symbol("p^- r_1", "p^- and r_1 belong to E-tilde"),
    )
    evidence = CoefficientSemiproductEvidence(
        "KKL2014TwoShifts Theorem 3.3",
        "942",
        8,
        source_checksum,
        (
            "p^pm are spatially constant E-tilde symbols",
            "the right symbol belongs to E-tilde",
            "both operators act on L^p(R_+,dmu) after Phi_delta",
        ),
    )
    cauchy_products = (Product((P_plus, R1)), Product((P_minus, R1)))
    cauchy_relations = tuple(
        ModCompactEquivalence(
            product,
            symbol,
            space="L^p(R_+,dmu), transported from the weighted space",
            compact_ideal="compact operators",
            evidence=evidence,
        )
        for product, symbol in zip(cauchy_products, cauchy_symbols, strict=True)
    )
    cauchy_rules = tuple(
        OperatorRule(
            f"cauchy_{sign}_mellin_semiproduct",
            DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
            relation,
            evidence.paper_hypotheses,
            f"KKL2014TwoShifts Theorem 3.3, printed 942/PDF 8, checksum {source_checksum}",
            RuleCertificationStatus.CERTIFIED_MOD_COMPACT,
        )
        for sign, relation in zip(("plus", "minus"), cauchy_relations, strict=True)
    )
    standard = {
        "L++": _standard_symbol(
            "p^+(lambda) r_1(x/gamma,lambda)",
            "radial scaling invariance of E-tilde",
        ),
        "L-Ghat": _standard_symbol(
            "p^-(lambda) r_1(x,lambda) Ghat_1(x)",
            "E-tilde Banach-algebra closure",
        ),
        "L+Ghat-base": _standard_symbol(
            "p^+(lambda) r_1(x/gamma,lambda) Ghat_1(x/gamma)",
            "radial scaling invariance of E-tilde",
        ),
        "L-+-base": cauchy_symbols[1],
    }
    forward = build_oscillatory_mellin_factor(lam, gamma, case=case)
    membership = (
        StandardMellinMembership.CERTIFIED
        if case is DilationFrequencyCase.IDENTITY
        else StandardMellinMembership.NOT_DEMONSTRATED
    )
    words = {
        "L-+": Product((P_minus, R1, U1)),
        "L++": Product((U1_inverse, P_plus, R1, U1)),
        "L-Ghat": Product((P_minus, R1, Ghat1)),
        "L+Ghat": Product((U1_inverse, P_plus, R1, Ghat1)),
    }
    outputs: dict[str, StandardMellinSymbol | FactorizedOscillatoryMellinSymbol] = {
        "L-+": FactorizedOscillatoryMellinSymbol(
            standard["L-+-base"], forward, words["L-+"],
            "tilde-E^(right-gamma)", membership,
            "bounded ordered factors", _right_dilation_obligations(),
        ),
        "L++": standard["L++"],
        "L-Ghat": standard["L-Ghat"],
        "L+Ghat": FactorizedOscillatoryMellinSymbol(
            standard["L+Ghat-base"], covariance.inverse_factor, words["L+Ghat"],
            "tilde-E^(right-gamma-inverse)", membership,
            "bounded ordered factors", _right_dilation_obligations(),
        ),
    }
    residuals = {
        "L-+": ("K_{-,r} D_gamma",),
        "L++": ("D_gamma^{-1} K_{+,r} D_gamma",),
        "L-Ghat": ("P^- K_{r,Ghat}", "K_{-,rGhat}"),
        "L+Ghat": ("D_gamma^{-1} P^+ K_{r,Ghat}", "D_gamma^{-1} K_{+,rGhat}"),
    }
    kinds = {
        "L-+": "RIGHT_FACTORIZED_MELLIN_PDO",
        "L++": "STANDARD_MELLIN_PDO",
        "L-Ghat": "STANDARD_MELLIN_PDO",
        "L+Ghat": "RIGHT_FACTORIZED_MELLIN_PDO",
    }
    cores = []
    for identifier in ("L-+", "L++", "L-Ghat", "L+Ghat"):
        relation = ModCompactEquivalence(
            words[identifier], outputs[identifier],
            space="weighted space, equivalently L^p(R_+,dmu) after Phi_delta",
            compact_ideal="compact operators",
            evidence=("exact covariance", "KKL2014TwoShifts Theorem 3.3"),
        )
        rule = OperatorRule(
            f"first_pivot_{identifier.replace('+', 'plus').replace('-', 'minus')}_closure",
            DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
            relation,
            ("gamma>0", "all standard factors belong to E-tilde"),
            f"Phase R derivation plus KKL2014TwoShifts Theorem 3.3; checksum {source_checksum}",
            RuleCertificationStatus.CERTIFIED_MOD_COMPACT,
        )
        cores.append(LeftCoreTrace(
            identifier, words[identifier], outputs[identifier], relation, rule,
            kinds[identifier], residuals[identifier],
            "Phi_delta^{-1}(unweighted relation)Phi_delta",
        ))
    return FirstPivotLeftCoreClosureTrace(
        covariance,
        cauchy_relations,
        cauchy_rules,
        tuple(cores),
        (
            "AKB is compact whenever K is compact and A,B are bounded",
            "finite sums of compact residuals remain compact",
            "Phi_delta conjugation preserves compactness",
        ),
        "right composition of each reduced core with the still-unidentified Wplus_12",
    )
