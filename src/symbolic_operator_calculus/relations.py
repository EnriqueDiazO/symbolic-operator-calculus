"""Minimal relation types for exact blocks and modulo-compact models."""

from __future__ import annotations

from dataclasses import dataclass

from .operators import LinearCombination, Product


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
    """Relation asserting exact block is equivalent to a model modulo compacta."""

    exact: ExactBlock
    model: WienerHopfModel

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
