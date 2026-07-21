"""Exact structural normalization of the complete first-pivot correction.

The global correction originates in the Phase N modulo-compact pivot model.
Only the two branch normalizations are exact supplied rules. This module keeps
those logical levels separate and reuses the existing noncommutative AST.
"""

from __future__ import annotations

from dataclasses import dataclass

from .normalized_schur import (
    NormalizedFirstSchurPivotDerivation,
    build_normalized_first_schur_pivot_derivation,
)
from .operators import (
    C2,
    C2_normalized,
    G1,
    G2,
    Ghat1,
    Ghat2,
    N2,
    N2_first,
    R1,
    T1_minus,
    T2_minus,
    U1,
    U2,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    Z1_inverse,
    Z2_inverse,
    LinearCombination,
    OperatorAtom,
    Product,
    compose_ordered,
)
from .semantics import (
    DerivationRelationKind,
    ExactIdentity,
    ExactIdentityScope,
    ModCompactEquivalence,
    OperatorRule,
    RuleCertificationStatus,
    SemanticDerivationStep,
)
from .substitution import substitute_operator_subproduct


class CompleteCorrectionError(ValueError):
    """Raised when the complete-correction normalization is inconsistent."""


@dataclass(frozen=True)
class CompleteCorrectionFactorization:
    """Validated grouped form over an authoritative expanded AST expression."""

    expanded: LinearCombination
    left_difference: LinearCombination
    left_wiener_hopf: OperatorAtom
    bridge: Product
    right_difference: LinearCombination
    right_wiener_hopf: OperatorAtom
    suffix: Product

    def __post_init__(self) -> None:
        if not isinstance(self.expanded, LinearCombination):
            raise TypeError("expanded must be a LinearCombination.")
        if not all(
            isinstance(item, LinearCombination)
            for item in (self.left_difference, self.right_difference)
        ):
            raise TypeError("differences must be LinearCombination objects.")
        if not all(
            isinstance(item, OperatorAtom)
            for item in (self.left_wiener_hopf, self.right_wiener_hopf)
        ):
            raise TypeError("Wiener--Hopf factors must be OperatorAtom objects.")
        if not all(isinstance(item, Product) for item in (self.bridge, self.suffix)):
            raise TypeError("bridge and suffix must be Product objects.")
        reconstructed = compose_ordered(
            self.left_difference,
            self.left_wiener_hopf,
            self.bridge,
            self.right_difference,
            self.right_wiener_hopf,
            self.suffix,
        )
        if reconstructed != self.expanded:
            raise CompleteCorrectionError(
                "factorization does not reconstruct its expanded expression."
            )


@dataclass(frozen=True)
class CompleteCorrectionDerivation:
    """Typed trace for B1--B4 of the complete correction normalization."""

    phase_n_trace: NormalizedFirstSchurPivotDerivation
    correction_definition: ExactIdentity
    pivot_relation: ModCompactEquivalence
    left_normalization_rule: OperatorRule
    right_normalization_rule: OperatorRule
    normalized_identity: ExactIdentity
    normalized_definition: ExactIdentity
    original_factorization: object
    normalized_factorization: CompleteCorrectionFactorization
    semantic_steps: tuple[SemanticDerivationStep, ...]
    rules: tuple[OperatorRule, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.phase_n_trace, NormalizedFirstSchurPivotDerivation):
            raise TypeError("phase_n_trace must be a Phase N derivation.")
        if not isinstance(self.correction_definition, ExactIdentity):
            raise TypeError("correction_definition must be exact structural data.")
        if not isinstance(self.pivot_relation, ModCompactEquivalence):
            raise TypeError("pivot_relation must remain modulo compact operators.")
        if not all(
            isinstance(item, OperatorRule)
            for item in (self.left_normalization_rule, self.right_normalization_rule)
        ):
            raise TypeError("normalization rules must be OperatorRule objects.")
        if not isinstance(self.normalized_identity, ExactIdentity):
            raise TypeError("normalized_identity must be an ExactIdentity.")
        if not isinstance(self.normalized_definition, ExactIdentity):
            raise TypeError("normalized_definition must be an ExactIdentity.")
        if not isinstance(self.normalized_factorization, CompleteCorrectionFactorization):
            raise TypeError("normalized_factorization has an unexpected type.")
        if not self.semantic_steps or not all(
            isinstance(item, SemanticDerivationStep) for item in self.semantic_steps
        ):
            raise TypeError("semantic_steps must contain typed semantic steps.")
        if len(self.rules) != 3 or not all(
            isinstance(item, OperatorRule) for item in self.rules
        ):
            raise TypeError("rules must contain exactly three OperatorRule objects.")
        if self.rules[1:] != (
            self.left_normalization_rule,
            self.right_normalization_rule,
        ):
            raise CompleteCorrectionError("rules must retain pivot, left, right order.")
        if self.normalized_definition.right != self.normalized_factorization.expanded:
            raise CompleteCorrectionError(
                "normalized definition and factorization must share one expression."
            )


def complete_correction_expression() -> LinearCombination:
    """Return the Phase N four-term correction before exact Z-normalization."""

    trace = build_normalized_first_schur_pivot_derivation()
    return trace.mod_compact_factorization.expanded_correction


def left_normalized_difference() -> LinearCombination:
    """Return the exact normalized left difference ``U2 - Ghat2``."""

    return U2 - Ghat2


def right_normalized_difference() -> LinearCombination:
    """Return the exact normalized right difference ``U1 - Ghat1``."""

    return U1 - Ghat1


def normalized_complete_correction_expression() -> LinearCombination:
    """Apply only the four contiguous consequences of the supplied Z rules."""

    expression = complete_correction_expression()
    substitutions = (
        (Z2_inverse * Vtilde_alpha2, U2),
        (Z2_inverse * G2, Ghat2),
        (Z1_inverse * Vtilde_alpha1, U1),
        (Z1_inverse * G1, Ghat1),
    )
    for target, replacement in substitutions:
        expression = substitute_operator_subproduct(
            expression,
            target,
            replacement,
        )
    return expression


def factor_normalized_complete_correction() -> CompleteCorrectionFactorization:
    """Validate C0 against the exact four-term normalized expansion C1."""

    return CompleteCorrectionFactorization(
        expanded=normalized_complete_correction_expression(),
        left_difference=left_normalized_difference(),
        left_wiener_hopf=Wminus_21,
        bridge=Product((T1_minus, R1)),
        right_difference=right_normalized_difference(),
        right_wiener_hopf=Wplus_12,
        suffix=Product((T2_minus,)),
    )


def _left_normalization_identity() -> ExactIdentity:
    return ExactIdentity(
        left=compose_ordered(
            Z2_inverse,
            Vtilde_alpha2 - G2,
        ),
        right=left_normalized_difference(),
        evidence="exact branch-2 transport and coefficient factorizations",
        scope=ExactIdentityScope.OPERATORIAL,
        hypotheses=("Z2 is invertible",),
    )


def _right_normalization_identity() -> ExactIdentity:
    return ExactIdentity(
        left=compose_ordered(
            Z1_inverse,
            Vtilde_alpha1 - G1,
        ),
        right=right_normalized_difference(),
        evidence="exact branch-1 transport and coefficient factorizations",
        scope=ExactIdentityScope.OPERATORIAL,
        hypotheses=("Z1 is invertible",),
    )


def build_complete_correction_derivation() -> CompleteCorrectionDerivation:
    """Build B1--B4 without changing the Phase N modulo-compact strength."""

    phase_n = build_normalized_first_schur_pivot_derivation()
    original = complete_correction_expression()
    normalized = normalized_complete_correction_expression()
    correction_definition = ExactIdentity(
        left=C2,
        right=original,
        evidence="definition extracted from the Phase N correction terms",
        scope=ExactIdentityScope.STRUCTURAL,
    )
    pivot_relation = ModCompactEquivalence(
        left=N2_first,
        right=N2 + C2,
    )
    left_rule = OperatorRule(
        name="normalize_left_branch_factor",
        relation_kind=DerivationRelationKind.EXACT_EQUALITY,
        relation=_left_normalization_identity(),
        preconditions=(
            "Vtilde_alpha2 = Z2 U2",
            "G2 = Z2 Ghat2",
            "Z2 is invertible",
        ),
        source=(
            "paper sections/05_haseman_operator_with_shift.tex, "
            "normalized transport and coefficient factorizations"
        ),
        certification_status=RuleCertificationStatus.CERTIFIED_EXACT,
    )
    right_rule = OperatorRule(
        name="normalize_right_branch_factor",
        relation_kind=DerivationRelationKind.EXACT_EQUALITY,
        relation=_right_normalization_identity(),
        preconditions=(
            "Vtilde_alpha1 = Z1 U1",
            "G1 = Z1 Ghat1",
            "Z1 is invertible",
        ),
        source=(
            "paper sections/05_haseman_operator_with_shift.tex, "
            "normalized transport and coefficient factorizations"
        ),
        certification_status=RuleCertificationStatus.CERTIFIED_EXACT,
    )
    pivot_rule = OperatorRule(
        name="complete_correction_from_first_pivot",
        relation_kind=DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
        relation=pivot_relation,
        preconditions=(
            "Phase N A21 Wiener--Hopf replacement",
            "Phase N A12 Wiener--Hopf replacement",
        ),
        source="Phase N supplied off-diagonal modulo-compact relations",
        certification_status=RuleCertificationStatus.UNCERTIFIED,
    )
    normalized_identity = ExactIdentity(
        left=original,
        right=normalized,
        evidence=(left_rule, right_rule),
        scope=ExactIdentityScope.OPERATORIAL,
        hypotheses=(
            "normalize_left_branch_factor preconditions",
            "normalize_right_branch_factor preconditions",
        ),
    )
    normalized_definition = ExactIdentity(
        left=C2_normalized,
        right=normalized,
        evidence="validated normalized correction factorization",
        scope=ExactIdentityScope.STRUCTURAL,
    )
    steps = (
        SemanticDerivationStep(
            "B1_pivot_correction",
            DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
            pivot_relation,
            pivot_rule,
        ),
        SemanticDerivationStep(
            "B2_left_normalization",
            DerivationRelationKind.EXACT_EQUALITY,
            left_rule.relation,
            left_rule,
        ),
        SemanticDerivationStep(
            "B3_right_normalization",
            DerivationRelationKind.EXACT_EQUALITY,
            right_rule.relation,
            right_rule,
        ),
        SemanticDerivationStep(
            "B4_normalized_correction",
            DerivationRelationKind.EXACT_EQUALITY,
            normalized_definition,
            (left_rule, right_rule),
        ),
    )
    return CompleteCorrectionDerivation(
        phase_n_trace=phase_n,
        correction_definition=correction_definition,
        pivot_relation=pivot_relation,
        left_normalization_rule=left_rule,
        right_normalization_rule=right_rule,
        normalized_identity=normalized_identity,
        normalized_definition=normalized_definition,
        original_factorization=phase_n.mod_compact_factorization,
        normalized_factorization=factor_normalized_complete_correction(),
        semantic_steps=steps,
        rules=(pivot_rule, left_rule, right_rule),
    )
