"""Formal normalization of the first Schur pivot.

This module performs only ordered finite algebra on the existing operator AST.
The supplied inverse and normalization identities are recorded explicitly, and
the Wiener--Hopf replacements remain modulo-compact declarations.  No analytic
operator-class or Fredholm conclusion is inferred.
"""

from __future__ import annotations

from dataclasses import dataclass

from .blocks import a12_mod_compact_relation, a21_mod_compact_relation
from .operators import (
    A11_inverse,
    A12,
    A21,
    A22,
    A22_first,
    G1,
    G2,
    N2,
    N2_first,
    P_minus,
    P_plus,
    R1,
    T1_minus,
    T2_minus,
    U1_inverse,
    U2_inverse,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    Z1_inverse,
    Z2_inverse,
    LinearCombination,
    OperatorAtom,
    Product,
    Term,
    compose_ordered,
)
from .relations import ModCompactRelation
from .semantics import (
    AnalyticProofObligation,
    DerivationRelationKind,
    ExactIdentity,
    ExactIdentityScope,
    FactorClassification,
    FactorStatus,
    FormalIdentity,
    ModCompactEquivalence,
    OperatorClass,
    SemanticDerivationStep,
)
from .substitution import substitute_operator_subproduct


class NormalizedSchurDerivationError(ValueError):
    """Raised when the normalized first-pivot trace is inconsistent."""


@dataclass(frozen=True)
class NormalizedModCompactFactorization:
    """Validated factorization of the expanded normalized WH correction.

    This is metadata over the existing AST, not a nested expression tree.
    ``expanded_correction`` is the authoritative operator expression.
    """

    expanded_correction: LinearCombination
    prefix: Product
    left_difference: LinearCombination
    left_wiener_hopf: OperatorAtom
    bridge: Product
    right_difference: LinearCombination
    right_wiener_hopf: OperatorAtom
    suffix: Product

    def __post_init__(self) -> None:
        if not isinstance(self.expanded_correction, LinearCombination):
            raise TypeError("expanded_correction must be a LinearCombination.")
        if not all(
            isinstance(product, Product)
            for product in (self.prefix, self.bridge, self.suffix)
        ):
            raise TypeError("prefix, bridge, and suffix must be Products.")
        if not all(
            isinstance(difference, LinearCombination)
            for difference in (self.left_difference, self.right_difference)
        ):
            raise TypeError("factorized differences must be LinearCombinations.")
        if not all(
            isinstance(atom, OperatorAtom)
            for atom in (self.left_wiener_hopf, self.right_wiener_hopf)
        ):
            raise TypeError("Wiener--Hopf factors must be OperatorAtoms.")
        reconstructed = compose_ordered(
            self.prefix,
            self.left_difference,
            self.left_wiener_hopf,
            self.bridge,
            self.right_difference,
            self.right_wiener_hopf,
            self.suffix,
        )
        if reconstructed != self.expanded_correction:
            raise NormalizedSchurDerivationError(
                "factorization does not reconstruct the expanded correction."
            )


@dataclass(frozen=True)
class NormalizedFirstSchurPivotDerivation:
    """Deterministic formal trace for the normalized first Schur pivot."""

    reduced_definition: ExactIdentity
    pivot_definition: ExactIdentity
    normalized_expansion: ExactIdentity
    inverse_substitution: FormalIdentity
    diagonal_recognition: ExactIdentity
    exact_result: ExactIdentity
    offdiagonal_relations: tuple[ModCompactRelation, ModCompactRelation]
    mod_compact_result: ModCompactEquivalence
    mod_compact_factorization: NormalizedModCompactFactorization
    semantic_steps: tuple[SemanticDerivationStep, ...]
    factor_classifications: tuple[FactorClassification, ...]
    proof_obligations: tuple[AnalyticProofObligation, ...]

    def __post_init__(self) -> None:
        exact_relations = (
            self.reduced_definition,
            self.pivot_definition,
            self.normalized_expansion,
            self.diagonal_recognition,
            self.exact_result,
        )
        if not all(isinstance(item, ExactIdentity) for item in exact_relations):
            raise TypeError("exact trace fields must be ExactIdentity objects.")
        if not isinstance(self.inverse_substitution, FormalIdentity):
            raise TypeError("inverse_substitution must be a FormalIdentity.")
        if len(self.offdiagonal_relations) != 2 or not all(
            isinstance(item, ModCompactRelation)
            for item in self.offdiagonal_relations
        ):
            raise TypeError(
                "offdiagonal_relations must contain A21 and A12 relations."
            )
        if not isinstance(self.mod_compact_result, ModCompactEquivalence):
            raise TypeError("mod_compact_result must be a ModCompactEquivalence.")
        if not isinstance(
            self.mod_compact_factorization,
            NormalizedModCompactFactorization,
        ):
            raise TypeError(
                "mod_compact_factorization must be normalized factorization metadata."
            )
        if not self.semantic_steps or not all(
            isinstance(item, SemanticDerivationStep)
            for item in self.semantic_steps
        ):
            raise TypeError("semantic_steps must contain SemanticDerivationStep objects.")
        if not self.factor_classifications or not all(
            isinstance(item, FactorClassification)
            for item in self.factor_classifications
        ):
            raise TypeError(
                "factor_classifications must contain FactorClassification objects."
            )
        if not self.proof_obligations or not all(
            isinstance(item, AnalyticProofObligation)
            for item in self.proof_obligations
        ):
            raise TypeError(
                "proof_obligations must contain AnalyticProofObligation objects."
            )
        expanded_terms = self.mod_compact_result.right.terms
        correction = LinearCombination(expanded_terms[1:])
        if correction != self.mod_compact_factorization.expanded_correction:
            raise NormalizedSchurDerivationError(
                "modulo-compact result and factorization disagree."
            )

    @property
    def expanded_mod_compact_result(self) -> LinearCombination:
        """Return the controlled expansion of only the two WH differences."""

        return self.mod_compact_result.right


def t1_minus_expression() -> LinearCombination:
    """Return the supplied structural formula ``U1^-1 P+ + P-``."""

    return compose_ordered(U1_inverse, P_plus) + P_minus


def t2_minus_expression() -> LinearCombination:
    """Return the supplied structural formula ``U2^-1 P+ + P-``."""

    return compose_ordered(U2_inverse, P_plus) + P_minus


def a22_first_schur_exact_expression() -> LinearCombination:
    """Build ``A22 - A21 A11^(-1) A12`` in exact factor order."""

    return LinearCombination(
        (
            Term(1, A22),
            Term(-1, A21 * A11_inverse * A12),
        )
    )


def normalized_a22_first_schur_expression() -> LinearCombination:
    """Distribute ``Z2^-1`` and ``T2,-`` around the exact reduction."""

    return compose_ordered(
        Z2_inverse,
        a22_first_schur_exact_expression(),
        T2_minus,
    )


def normalized_a22_with_regularizer_expression() -> LinearCombination:
    """Formally substitute ``A11^(-1) = T1,- R1 Z1^-1``."""

    return substitute_operator_subproduct(
        normalized_a22_first_schur_expression(),
        A11_inverse,
        T1_minus * R1 * Z1_inverse,
    )


def normalized_first_schur_exact_expression() -> LinearCombination:
    """Recognize ``Z2^-1 A22 T2,-`` as ``N2`` in the formal expression."""

    return substitute_operator_subproduct(
        normalized_a22_with_regularizer_expression(),
        Z2_inverse * A22 * T2_minus,
        N2,
    )


def normalized_first_schur_mod_compact_expression() -> LinearCombination:
    """Apply exactly the two supplied off-diagonal WH replacements."""

    expression = normalized_first_schur_exact_expression()
    expression = substitute_operator_subproduct(
        expression,
        A21,
        a21_mod_compact_relation().model.expression,
    )
    return substitute_operator_subproduct(
        expression,
        A12,
        a12_mod_compact_relation().model.expression,
    )


def factor_normalized_mod_compact_correction() -> NormalizedModCompactFactorization:
    """Validate the requested grouped WH form against its four-term expansion."""

    expanded = normalized_first_schur_mod_compact_expression()
    return NormalizedModCompactFactorization(
        expanded_correction=LinearCombination(expanded.terms[1:]),
        prefix=Product((Z2_inverse,)),
        left_difference=Vtilde_alpha2 - G2,
        left_wiener_hopf=Wminus_21,
        bridge=Product((T1_minus, R1, Z1_inverse)),
        right_difference=Vtilde_alpha1 - G1,
        right_wiener_hopf=Wplus_12,
        suffix=Product((T2_minus,)),
    )


def normalized_first_schur_factor_classifications(
) -> tuple[FactorClassification, ...]:
    """Return the requested factor table as explicit caller-supplied metadata."""

    return (
        FactorClassification(
            Z2_inverse,
            OperatorClass.MULTIPLICATION_OPERATOR,
            FactorStatus.EXACT,
            "definición del transporte",
        ),
        FactorClassification(
            Vtilde_alpha2,
            OperatorClass.MULTIPLICATION_DILATION,
            FactorStatus.EXACT,
            "transporte del shift",
        ),
        FactorClassification(
            G2,
            OperatorClass.MULTIPLICATION_OPERATOR,
            FactorStatus.EXACT,
            "coeficiente branchwise",
        ),
        FactorClassification(
            Wminus_21,
            OperatorClass.LOCALIZED_WIENER_HOPF_OPERATOR,
            FactorStatus.MOD_COMPACT,
            "localización cuspídea 2025",
        ),
        FactorClassification(
            T1_minus,
            OperatorClass.AUXILIARY_SHIFT_OPERATOR,
            FactorStatus.INVERTIBLE,
            "resultado branchwise",
        ),
        FactorClassification(
            R1,
            OperatorClass.MELLIN_PDO_REGULARIZER,
            FactorStatus.REGULARIZER_MOD_COMPACT,
            "KKL 2014, Thm. 5.8",
        ),
        FactorClassification(
            Z1_inverse,
            OperatorClass.MULTIPLICATION_OPERATOR,
            FactorStatus.EXACT,
            "normalización diagonal",
        ),
        FactorClassification(
            Vtilde_alpha1,
            OperatorClass.MULTIPLICATION_DILATION,
            FactorStatus.EXACT,
            "transporte del shift",
        ),
        FactorClassification(
            G1,
            OperatorClass.MULTIPLICATION_OPERATOR,
            FactorStatus.EXACT,
            "coeficiente branchwise",
        ),
        FactorClassification(
            Wplus_12,
            OperatorClass.LOCALIZED_WIENER_HOPF_OPERATOR,
            FactorStatus.MOD_COMPACT,
            "localización cuspídea 2025",
        ),
        FactorClassification(
            T2_minus,
            OperatorClass.AUXILIARY_SHIFT_OPERATOR,
            FactorStatus.INVERTIBLE,
            "resultado branchwise",
        ),
    )


def normalized_first_schur_proof_obligations(
) -> tuple[AnalyticProofObligation, ...]:
    """Return the analytic questions deliberately left open by this phase."""

    return (
        AnalyticProofObligation(
            "justify_wh_mod_compact",
            "Justificar la sustitución Wiener--Hopf módulo compactos.",
        ),
        AnalyticProofObligation(
            "classify_complete_correction",
            "Demostrar que la corrección completa pertenece a una álgebra cerrada apropiada.",
        ),
        AnalyticProofObligation(
            "prove_product_closure",
            "Demostrar cierre bajo los productos que aparecen.",
        ),
        AnalyticProofObligation(
            "choose_analytic_framework",
            "Decidir si la salida es un Mellin PDO admisible o un elemento de la álgebra cuspídea de 2025.",
        ),
        AnalyticProofObligation(
            "apply_fredholm_criterion_later",
            "Aplicar después un criterio de Fredholmness, sin afirmar aún su conclusión.",
        ),
    )


def build_normalized_first_schur_pivot_derivation(
) -> NormalizedFirstSchurPivotDerivation:
    """Build the complete formal and modulo-compact normalized-pivot trace."""

    reduced_expression = a22_first_schur_exact_expression()
    normalized_expression = normalized_a22_first_schur_expression()
    with_regularizer = normalized_a22_with_regularizer_expression()
    exact_expression = normalized_first_schur_exact_expression()
    mod_compact_expression = normalized_first_schur_mod_compact_expression()
    pivot_word = Z2_inverse * A22_first * T2_minus
    diagonal_word = Z2_inverse * A22 * T2_minus
    inverse_word = T1_minus * R1 * Z1_inverse
    left_relation = a21_mod_compact_relation()
    right_relation = a12_mod_compact_relation()
    reduced_definition = ExactIdentity(
        left=A22_first,
        right=reduced_expression,
        evidence="definition of the first Schur pivot",
        scope=ExactIdentityScope.STRUCTURAL,
    )
    pivot_definition = ExactIdentity(
        left=N2_first,
        right=pivot_word,
        evidence="definition of normalized first Schur pivot",
        scope=ExactIdentityScope.STRUCTURAL,
    )
    normalized_expansion = ExactIdentity(
        left=pivot_word,
        right=normalized_expression,
        evidence="finite distributivity in the ordered operator AST",
        scope=ExactIdentityScope.STRUCTURAL,
    )
    inverse_substitution = FormalIdentity(
        left=normalized_expression,
        right=with_regularizer,
        justification=FormalIdentity(
            left=A11_inverse,
            right=inverse_word,
            justification="supplied formal regularizer identity",
        ),
    )
    diagonal_recognition = ExactIdentity(
        left=diagonal_word,
        right=N2,
        evidence="supplied definition of N2",
        scope=ExactIdentityScope.STRUCTURAL,
    )
    exact_result = ExactIdentity(
        left=N2_first,
        right=exact_expression,
        evidence="ordered formal substitution and diagonal recognition",
        scope=ExactIdentityScope.STRUCTURAL,
        hypotheses=(
            "A11^(-1) = T1,- R1 Z1^(-1)",
            "Z2^(-1) A22 T2,- = N2",
        ),
    )
    mod_compact_result = ModCompactEquivalence(
        left=N2_first,
        right=mod_compact_expression,
    )
    proof_obligations = normalized_first_schur_proof_obligations()
    semantic_steps = (
        SemanticDerivationStep(
            "reduced_definition",
            DerivationRelationKind.EXACT_EQUALITY,
            reduced_definition,
        ),
        SemanticDerivationStep(
            "pivot_definition",
            DerivationRelationKind.EXACT_EQUALITY,
            pivot_definition,
        ),
        SemanticDerivationStep(
            "normalized_expansion",
            DerivationRelationKind.EXACT_EQUALITY,
            normalized_expansion,
        ),
        SemanticDerivationStep(
            "inverse_substitution",
            DerivationRelationKind.FORMAL_SUBSTITUTION,
            inverse_substitution,
        ),
        SemanticDerivationStep(
            "diagonal_recognition",
            DerivationRelationKind.EXACT_EQUALITY,
            diagonal_recognition,
        ),
        SemanticDerivationStep(
            "exact_result",
            DerivationRelationKind.EXACT_EQUALITY,
            exact_result,
        ),
        SemanticDerivationStep(
            "a21_mod_compact",
            DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
            left_relation.semantic_relation,
        ),
        SemanticDerivationStep(
            "a12_mod_compact",
            DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
            right_relation.semantic_relation,
        ),
        SemanticDerivationStep(
            "normalized_mod_compact_result",
            DerivationRelationKind.MOD_COMPACT_EQUIVALENCE,
            mod_compact_result,
        ),
        *(
            SemanticDerivationStep(
                obligation.key,
                obligation.relation_kind,
                obligation,
            )
            for obligation in proof_obligations
        ),
    )
    return NormalizedFirstSchurPivotDerivation(
        reduced_definition=reduced_definition,
        pivot_definition=pivot_definition,
        normalized_expansion=normalized_expansion,
        inverse_substitution=inverse_substitution,
        diagonal_recognition=diagonal_recognition,
        exact_result=exact_result,
        offdiagonal_relations=(left_relation, right_relation),
        mod_compact_result=mod_compact_result,
        mod_compact_factorization=factor_normalized_mod_compact_correction(),
        semantic_steps=semantic_steps,
        factor_classifications=normalized_first_schur_factor_classifications(),
        proof_obligations=proof_obligations,
    )
