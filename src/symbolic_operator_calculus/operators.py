"""Small structural language for ordered noncommutative operators.

This module intentionally does not use SymPy ``Mul`` or ``Add`` to represent
operator products or ordered linear combinations. Order is stored directly in
tuples and expansion is an explicit operation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, TypeAlias

Sign: TypeAlias = Literal[-1, 1]


def _checked_sign(sign: int) -> Sign:
    if sign not in (-1, 1):
        raise ValueError("A term sign must be +1 or -1.")
    return sign  # type: ignore[return-value]


@dataclass(frozen=True)
class OperatorAtom:
    """An indivisible noncommutative operator symbol."""

    name: str
    kind: str = "operator"

    @property
    def is_commutative(self) -> bool:
        """All operator atoms in this MVP are noncommutative."""

        return False

    @property
    def is_formal_regularizer(self) -> bool:
        """Whether this atom denotes a formal regularizer."""

        return self.kind == "formal_regularizer"

    def __mul__(self, other: Factor) -> Product:
        return Product((self,)) * other

    def __rmul__(self, other: Factor) -> Product:
        return Product.from_factor(other) * self

    def __sub__(self, other: Factor) -> LinearCombination:
        return LinearCombination.from_difference(self, other)


@dataclass(frozen=True)
class Product:
    """An ordered product of atoms and unexpanded formal factors."""

    factors: tuple[Factor, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "factors", tuple(self.factors))

    @classmethod
    def from_factor(cls, factor: Factor) -> Product:
        if isinstance(factor, Product):
            return factor
        return cls((factor,))

    def __mul__(self, other: Factor) -> Product:
        if isinstance(other, Product):
            return Product(self.factors + other.factors)
        return Product(self.factors + (other,))

    def __rmul__(self, other: Factor) -> Product:
        return Product.from_factor(other) * self

    def __sub__(self, other: Factor) -> LinearCombination:
        return LinearCombination.from_difference(self, other)


@dataclass(frozen=True)
class Term:
    """A signed ordered product."""

    sign: Sign
    product: Product

    def __post_init__(self) -> None:
        object.__setattr__(self, "sign", _checked_sign(self.sign))
        object.__setattr__(self, "product", Product.from_factor(self.product))


@dataclass(frozen=True)
class LinearCombination:
    """An ordered finite linear combination of signed products."""

    terms: tuple[Term, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "terms", tuple(self.terms))

    @classmethod
    def from_factor(cls, factor: Factor) -> LinearCombination:
        if isinstance(factor, LinearCombination):
            return factor
        return cls((Term(1, Product.from_factor(factor)),))

    @classmethod
    def from_difference(cls, left: Factor, right: Factor) -> LinearCombination:
        left_terms = cls.from_factor(left).terms
        right_terms = tuple(
            Term(_checked_sign(-term.sign), term.product)
            for term in cls.from_factor(right).terms
        )
        return cls(left_terms + right_terms)

    def __mul__(self, other: Factor) -> Product:
        return Product((self,)) * other

    def __rmul__(self, other: Factor) -> Product:
        return Product.from_factor(other) * self

    def __sub__(self, other: Factor) -> LinearCombination:
        return LinearCombination.from_difference(self, other)


Factor: TypeAlias = OperatorAtom | Product | LinearCombination


def expand_ordered(expression: Factor) -> LinearCombination:
    """Expand distributively while preserving declared left-to-right order."""

    if isinstance(expression, LinearCombination):
        expanded_terms: list[Term] = []
        for term in expression.terms:
            expanded_product = expand_ordered(term.product)
            expanded_terms.extend(
                Term(_checked_sign(term.sign * expanded.sign), expanded.product)
                for expanded in expanded_product.terms
            )
        return LinearCombination(tuple(expanded_terms))

    if isinstance(expression, OperatorAtom):
        return LinearCombination.from_factor(expression)

    if isinstance(expression, Product):
        expanded_terms = (Term(1, Product(())),)
        for factor in expression.factors:
            factor_terms = expand_ordered(factor).terms
            expanded_terms = tuple(
                Term(
                    _checked_sign(left.sign * right.sign),
                    left.product * right.product,
                )
                for left in expanded_terms
                for right in factor_terms
            )
        return LinearCombination(expanded_terms)

    raise TypeError(f"Unsupported expression type: {type(expression)!r}")


Vtilde_alpha2 = OperatorAtom("Vtilde_alpha2")
G2 = OperatorAtom("G2")
Wminus_21 = OperatorAtom("Wminus_21")
R11 = OperatorAtom("R11", kind="formal_regularizer")
Vtilde_alpha1 = OperatorAtom("Vtilde_alpha1")
G1 = OperatorAtom("G1")
Wplus_12 = OperatorAtom("Wplus_12")


def main_expression() -> Product:
    """Build the unexpanded principal expression for Phase B."""

    return (
        (Vtilde_alpha2 - G2)
        * Wminus_21
        * R11
        * (Vtilde_alpha1 - G1)
        * Wplus_12
    )
