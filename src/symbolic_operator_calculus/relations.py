"""Minimal relation types for exact blocks and modulo-compact models."""

from __future__ import annotations

from dataclasses import dataclass, field

from .operators import LinearCombination, OperatorAtom, Product
from .semantics import (
    CertificationStatus,
    ModCompactEquivalence,
    RegularizerOperator,
)


def _validate_family(family: object) -> str:
    if not isinstance(family, str) or not family:
        raise TypeError("family must be a non-empty string.")
    return family


def _validate_positive_index(name: str, value: object) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise TypeError(f"{name} must be a positive integer.")
    if value <= 0:
        raise ValueError(f"{name} must be positive.")
    return value


@dataclass(frozen=True)
class ExactBlock:
    """Structural reference to an exact operator block, without a formula."""

    family: str
    row: int
    column: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "family", _validate_family(self.family))
        object.__setattr__(self, "row", _validate_positive_index("row", self.row))
        object.__setattr__(
            self,
            "column",
            _validate_positive_index("column", self.column),
        )


@dataclass(frozen=True)
class FormalRegularizer(RegularizerOperator):
    """Metadata associating a formal regularizer atom with an exact block.

    This association records neither an exact inverse nor an oriented
    product-identity relation.
    """

    target: ExactBlock
    operator: OperatorAtom

    def __post_init__(self) -> None:
        RegularizerOperator.__post_init__(self)
        if not isinstance(self.target, ExactBlock):
            raise TypeError("target must be an ExactBlock.")
        if not isinstance(self.operator, OperatorAtom):
            raise TypeError("operator must be an OperatorAtom.")
        if not self.operator.is_formal_regularizer:
            raise ValueError("operator must be marked as a formal regularizer.")


@dataclass(frozen=True)
class FirstSchurReduction:
    """Formal algebraic metadata for ``diagonal - left * regularizer * right``."""

    diagonal: ExactBlock
    left: ExactBlock
    regularizer: FormalRegularizer
    right: ExactBlock
    offdiagonal_product_coefficient: int = -1

    def __post_init__(self) -> None:
        if not isinstance(self.diagonal, ExactBlock):
            raise TypeError("diagonal must be an ExactBlock.")
        if not isinstance(self.left, ExactBlock):
            raise TypeError("left must be an ExactBlock.")
        if not isinstance(self.regularizer, FormalRegularizer):
            raise TypeError("regularizer must be a FormalRegularizer.")
        if not isinstance(self.right, ExactBlock):
            raise TypeError("right must be an ExactBlock.")
        if self.diagonal != ExactBlock("A", 2, 2):
            raise ValueError("diagonal must be exact block A22.")
        if self.left != ExactBlock("A", 2, 1):
            raise ValueError("left must be exact block A21.")
        if self.regularizer.target != ExactBlock("A", 1, 1):
            raise ValueError("regularizer must target exact block A11.")
        if self.right != ExactBlock("A", 1, 2):
            raise ValueError("right must be exact block A12.")
        if self.offdiagonal_product_coefficient != -1:
            raise ValueError("the first Schur reduction requires exterior sign -1.")


@dataclass(frozen=True)
class ModCompactSchurRelation:
    """Modulo-compact declaration from a Schur reduction to an AST model.

    With evidence it is marked ``EVIDENCE_SUPPLIED``, not certified by the
    program. No compact residual is inferred from the endpoint types.
    """

    exact: FirstSchurReduction
    model: LinearCombination
    space: object | None = None
    compact_ideal: object | None = None
    residual: object | None = None
    evidence: object | None = None
    semantic_relation: ModCompactEquivalence = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.exact, FirstSchurReduction):
            raise TypeError("exact must be a FirstSchurReduction.")
        if not isinstance(self.model, LinearCombination):
            raise TypeError("model must be a LinearCombination.")
        object.__setattr__(
            self,
            "semantic_relation",
            ModCompactEquivalence(
                left=self.exact,
                right=self.model,
                space=self.space,
                compact_ideal=self.compact_ideal,
                residual=self.residual,
                evidence=self.evidence,
            ),
        )

    @property
    def certification_status(self) -> CertificationStatus:
        """Expose whether external compactness evidence was supplied."""

        return self.semantic_relation.certification_status


@dataclass(frozen=True)
class WienerHopfModel:
    """Normalized Wiener-Hopf model for an exact block."""

    family: str
    row: int
    column: int
    expression: Product | LinearCombination
    normalized: bool = True

    def __post_init__(self) -> None:
        object.__setattr__(self, "family", _validate_family(self.family))
        object.__setattr__(self, "row", _validate_positive_index("row", self.row))
        object.__setattr__(
            self,
            "column",
            _validate_positive_index("column", self.column),
        )
        if not isinstance(self.expression, (Product, LinearCombination)):
            raise TypeError("expression must be a Product or LinearCombination.")
        if not isinstance(self.normalized, bool):
            raise TypeError("normalized must be a bool.")


@dataclass(frozen=True)
class ModCompactRelation:
    """Modulo-compact declaration between an exact block and a model.

    This compatibility wrapper retains the existing ``exact`` and ``model``
    endpoints while exposing an explicit :class:`ModCompactEquivalence`.
    It is uncertified by default. Evidence changes its status only to
    ``EVIDENCE_SUPPLIED`` and never proves that its residual is compact.
    """

    exact: ExactBlock
    model: WienerHopfModel
    space: object | None = None
    compact_ideal: object | None = None
    residual: object | None = None
    evidence: object | None = None
    semantic_relation: ModCompactEquivalence = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.exact, ExactBlock):
            raise TypeError("exact must be an ExactBlock.")
        if not isinstance(self.model, WienerHopfModel):
            raise TypeError("model must be a WienerHopfModel.")
        if (
            self.exact.family,
            self.exact.row,
            self.exact.column,
        ) != (
            self.model.family,
            self.model.row,
            self.model.column,
        ):
            raise ValueError("exact block and model must reference the same block.")
        object.__setattr__(
            self,
            "semantic_relation",
            ModCompactEquivalence(
                left=self.exact,
                right=self.model,
                space=self.space,
                compact_ideal=self.compact_ideal,
                residual=self.residual,
                evidence=self.evidence,
            ),
        )

    @property
    def certification_status(self) -> CertificationStatus:
        """Expose whether external compactness evidence was supplied."""

        return self.semantic_relation.certification_status
