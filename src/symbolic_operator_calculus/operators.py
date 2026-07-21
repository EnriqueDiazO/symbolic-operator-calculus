"""Small structural language for ordered noncommutative operators.

This module intentionally does not use SymPy ``Mul`` or ``Add`` to represent
operator products or ordered linear combinations. Order is stored directly in
tuples and expansion is an explicit operation.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite
from typing import TypeAlias

Scalar: TypeAlias = int | float | complex | Fraction
OperatorExpression: TypeAlias = "OperatorAtom | Product | LinearCombination"


def _is_scalar(value: object) -> bool:
    return isinstance(value, (int, float, complex, Fraction)) and not isinstance(
        value,
        bool,
    )


def _is_zero_scalar(value: object) -> bool:
    return _is_scalar(value) and value == 0


def _checked_coefficient(coefficient: object) -> Scalar:
    if not _is_scalar(coefficient):
        raise TypeError(
            "Term.coefficient must be an int, float, complex, or Fraction."
        )
    if isinstance(coefficient, float) and not isfinite(coefficient):
        raise ValueError("Term.coefficient must be finite.")
    if isinstance(coefficient, complex) and not (
        isfinite(coefficient.real) and isfinite(coefficient.imag)
    ):
        raise ValueError("Term.coefficient must be finite.")
    return coefficient


def _validated_nonempty_string(
    value: object,
    *,
    class_name: str,
    field_name: str,
) -> str:
    if not isinstance(value, str):
        raise TypeError(
            f"{class_name}.{field_name} must be a nonempty str; got "
            f"{type(value).__name__}."
        )
    if not value.strip():
        raise ValueError(f"{class_name}.{field_name} must be a nonempty str.")
    return value


def _tuple_from_iterable(
    value: object,
    *,
    class_name: str,
    field_name: str,
) -> tuple[object, ...]:
    try:
        return tuple(value)
    except TypeError as exc:
        raise TypeError(
            f"{class_name}.{field_name} must be an iterable; got "
            f"{type(value).__name__}."
        ) from exc


def _unsupported_distribution() -> TypeError:
    return TypeError(
        "Products involving LinearCombination require distributivity and are "
        "not supported in this phase."
    )


@dataclass(frozen=True)
class OperatorAtom:
    """An indivisible noncommutative operator symbol."""

    name: str
    kind: str = "operator"

    def __post_init__(self) -> None:
        _validated_nonempty_string(
            self.name,
            class_name="OperatorAtom",
            field_name="name",
        )
        _validated_nonempty_string(
            self.kind,
            class_name="OperatorAtom",
            field_name="kind",
        )

    @property
    def is_commutative(self) -> bool:
        """All operator atoms in this MVP are noncommutative."""

        return False

    @property
    def is_formal_regularizer(self) -> bool:
        """Whether this atom denotes a formal regularizer."""

        return self.kind == "formal_regularizer"

    def __mul__(self, other: object) -> Product | LinearCombination:
        if _is_scalar(other):
            return LinearCombination.from_scaled_factor(other, self)
        if isinstance(other, LinearCombination):
            raise _unsupported_distribution()
        if isinstance(other, (OperatorAtom, Product)):
            return Product.from_factor(self) * other
        return NotImplemented

    def __rmul__(self, other: object) -> Product | LinearCombination:
        if _is_scalar(other):
            return LinearCombination.from_scaled_factor(other, self)
        if isinstance(other, LinearCombination):
            raise _unsupported_distribution()
        if isinstance(other, Product):
            return other * self
        return NotImplemented

    def __add__(self, other: object) -> LinearCombination:
        return LinearCombination.from_factor(self) + other

    def __radd__(self, other: object) -> LinearCombination:
        return LinearCombination.from_factor(self).__radd__(other)

    def __sub__(self, other: object) -> LinearCombination:
        return LinearCombination.from_factor(self) - other

    def __rsub__(self, other: object) -> LinearCombination:
        return LinearCombination.from_factor(self).__rsub__(other)

    def __neg__(self) -> LinearCombination:
        return LinearCombination.from_scaled_factor(-1, self)


@dataclass(frozen=True)
class Product:
    """An ordered product of operator atoms.

    ``Product(())`` is the structural composition identity (the empty word).
    It remains structurally distinct from the explicit identity atom product
    ``Product((I,))``, although both have the same action. Non-empty products
    may contain only ``OperatorAtom`` factors.
    """

    factors: tuple[OperatorAtom, ...]

    def __post_init__(self) -> None:
        factors = _tuple_from_iterable(
            self.factors,
            class_name="Product",
            field_name="factors",
        )
        for factor in factors:
            if not isinstance(factor, OperatorAtom):
                raise TypeError(
                    "Product.factors items must be OperatorAtom instances; got "
                    f"{type(factor).__name__}."
                )
        object.__setattr__(self, "factors", factors)

    @classmethod
    def from_factor(cls, factor: OperatorAtom | Product) -> Product:
        if isinstance(factor, Product):
            return factor
        if isinstance(factor, OperatorAtom):
            return cls((factor,))
        raise TypeError(
            "Product.from_factor expects an OperatorAtom or Product, got "
            f"{type(factor).__name__}."
        )

    def __mul__(self, other: object) -> Product | LinearCombination:
        if _is_scalar(other):
            return LinearCombination.from_scaled_factor(other, self)
        if isinstance(other, LinearCombination):
            raise _unsupported_distribution()
        if isinstance(other, Product):
            return Product(self.factors + other.factors)
        if isinstance(other, OperatorAtom):
            return Product(self.factors + (other,))
        return NotImplemented

    def __rmul__(self, other: object) -> Product | LinearCombination:
        if _is_scalar(other):
            return LinearCombination.from_scaled_factor(other, self)
        if isinstance(other, LinearCombination):
            raise _unsupported_distribution()
        if isinstance(other, OperatorAtom):
            return Product((other,) + self.factors)
        return NotImplemented

    def __add__(self, other: object) -> LinearCombination:
        return LinearCombination.from_factor(self) + other

    def __radd__(self, other: object) -> LinearCombination:
        return LinearCombination.from_factor(self).__radd__(other)

    def __sub__(self, other: object) -> LinearCombination:
        return LinearCombination.from_factor(self) - other

    def __rsub__(self, other: object) -> LinearCombination:
        return LinearCombination.from_factor(self).__rsub__(other)

    def __neg__(self) -> LinearCombination:
        return LinearCombination.from_scaled_factor(-1, self)


@dataclass(frozen=True)
class Term:
    """A scalar coefficient multiplying an ordered product.

    An ``OperatorAtom`` is accepted as shorthand and stored as a one-factor
    ``Product``.
    """

    coefficient: Scalar
    product: Product | OperatorAtom

    def __post_init__(self) -> None:
        object.__setattr__(self, "coefficient", _checked_coefficient(self.coefficient))
        if not isinstance(self.product, (OperatorAtom, Product)):
            raise TypeError(
                "Term.product must be an OperatorAtom or Product; got "
                f"{type(self.product).__name__}."
            )
        object.__setattr__(self, "product", Product.from_factor(self.product))


@dataclass(frozen=True)
class LinearCombination:
    """An ordered finite linear combination of scalar-weighted products.

    ``LinearCombination(())`` is the structural additive zero.
    """

    terms: tuple[Term, ...]

    def __post_init__(self) -> None:
        terms = _tuple_from_iterable(
            self.terms,
            class_name="LinearCombination",
            field_name="terms",
        )
        for term in terms:
            if not isinstance(term, Term):
                raise TypeError(
                    "LinearCombination.terms items must be Term instances; got "
                    f"{type(term).__name__}."
                )
        object.__setattr__(
            self,
            "terms",
            tuple(term for term in terms if term.coefficient != 0),
        )

    @classmethod
    def from_factor(cls, factor: object) -> LinearCombination:
        if isinstance(factor, LinearCombination):
            return factor
        if isinstance(factor, (OperatorAtom, Product)):
            return cls((Term(1, Product.from_factor(factor)),))
        if _is_zero_scalar(factor):
            return cls(())
        raise TypeError(
            "LinearCombination.from_factor expects an operator expression or "
            f"zero scalar, got {type(factor).__name__}."
        )

    @classmethod
    def from_scaled_factor(
        cls,
        coefficient: object,
        factor: OperatorAtom | Product,
    ) -> LinearCombination:
        coefficient = _checked_coefficient(coefficient)
        if coefficient == 0:
            return cls(())
        return cls((Term(coefficient, Product.from_factor(factor)),))

    @classmethod
    def from_difference(cls, left: object, right: object) -> LinearCombination:
        return cls.from_factor(left) - right

    def __add__(self, other: object) -> LinearCombination:
        if _is_zero_scalar(other):
            return self
        if isinstance(other, (OperatorAtom, Product)):
            return LinearCombination(self.terms + LinearCombination.from_factor(other).terms)
        if isinstance(other, LinearCombination):
            return LinearCombination(self.terms + other.terms)
        return NotImplemented

    def __radd__(self, other: object) -> LinearCombination:
        if _is_zero_scalar(other):
            return self
        if isinstance(other, (OperatorAtom, Product)):
            return LinearCombination.from_factor(other) + self
        return NotImplemented

    def __sub__(self, other: object) -> LinearCombination:
        if _is_zero_scalar(other):
            return self
        if isinstance(other, (OperatorAtom, Product, LinearCombination)):
            return self + (-LinearCombination.from_factor(other))
        return NotImplemented

    def __rsub__(self, other: object) -> LinearCombination:
        if _is_zero_scalar(other):
            return -self
        if isinstance(other, (OperatorAtom, Product)):
            return LinearCombination.from_factor(other) - self
        return NotImplemented

    def __neg__(self) -> LinearCombination:
        return LinearCombination(
            tuple(Term(-term.coefficient, term.product) for term in self.terms)
        )

    def __mul__(self, other: object) -> LinearCombination:
        if _is_scalar(other):
            return LinearCombination(
                tuple(Term(term.coefficient * other, term.product) for term in self.terms)
            )
        raise _unsupported_distribution()

    def __rmul__(self, other: object) -> LinearCombination:
        if _is_scalar(other):
            return self * other
        raise _unsupported_distribution()


Factor: TypeAlias = OperatorExpression


def expand_ordered(expression: Factor) -> LinearCombination:
    """Expand the current structural expression into ordered product terms."""

    if isinstance(expression, LinearCombination):
        expanded_terms: list[Term] = []
        for term in expression.terms:
            expanded_product = expand_ordered(term.product)
            expanded_terms.extend(
                Term(term.coefficient * expanded.coefficient, expanded.product)
                for expanded in expanded_product.terms
            )
        return LinearCombination(tuple(expanded_terms))

    if isinstance(expression, OperatorAtom):
        return LinearCombination.from_factor(expression)

    if isinstance(expression, Product):
        return LinearCombination((Term(1, expression),))

    raise TypeError(f"Unsupported expression type: {type(expression)!r}")


def compose_ordered(*expressions: Factor) -> LinearCombination:
    """Explicitly distribute and compose operator expressions in stored order.

    The overloaded multiplication operators intentionally reject products
    involving :class:`LinearCombination`.  This function is the opt-in path
    for the finite distributive expansion needed by formal derivations.  It
    never sorts, commutes, cancels, or combines terms.
    """

    result = LinearCombination((Term(1, Product(())),))
    for expression in expressions:
        if not isinstance(expression, (OperatorAtom, Product, LinearCombination)):
            raise TypeError(
                "compose_ordered expects operator expressions; got "
                f"{type(expression).__name__}."
            )
        right = LinearCombination.from_factor(expression)
        result = LinearCombination(
            tuple(
                Term(
                    left_term.coefficient * right_term.coefficient,
                    Product(
                        left_term.product.factors + right_term.product.factors
                    ),
                )
                for left_term in result.terms
                for right_term in right.terms
            )
        )
    return result


Vtilde_alpha2 = OperatorAtom("Vtilde_alpha2")
G2 = OperatorAtom("G2")
Wminus_21 = OperatorAtom("Wminus_21")
R11 = OperatorAtom("R11", kind="formal_regularizer")
Vtilde_alpha1 = OperatorAtom("Vtilde_alpha1")
G1 = OperatorAtom("G1")
Wplus_12 = OperatorAtom("Wplus_12")
I = OperatorAtom("I")  # noqa: E741 - mathematical identity operator
S_Rplus = OperatorAtom("S_Rplus")

# Atoms used by the normalized first-Schur-pivot derivation.  They remain
# indivisible noncommutative symbols; the ``kind`` labels are descriptive and
# do not prove membership in an operator class.
Z1_inverse = OperatorAtom("Z1_inverse", kind="multiplication_operator")
Z2_inverse = OperatorAtom("Z2_inverse", kind="multiplication_operator")
A12 = OperatorAtom("A12", kind="exact_block")
A21 = OperatorAtom("A21", kind="exact_block")
A22 = OperatorAtom("A22", kind="exact_block")
A11_inverse = OperatorAtom("A11_inverse", kind="formal_inverse")
A22_first = OperatorAtom("A22_first", kind="formal_schur_pivot")
T1_minus = OperatorAtom("T1_minus", kind="auxiliary_shift_operator")
T2_minus = OperatorAtom("T2_minus", kind="auxiliary_shift_operator")
U1_inverse = OperatorAtom("U1_inverse", kind="inverse_shift_operator")
U2_inverse = OperatorAtom("U2_inverse", kind="inverse_shift_operator")
P_plus = OperatorAtom("P_plus", kind="projection")
P_minus = OperatorAtom("P_minus", kind="projection")
R1 = OperatorAtom("R1", kind="formal_regularizer")
N2 = OperatorAtom("N2", kind="normalized_pivot")
N2_first = OperatorAtom("N2_first", kind="normalized_schur_pivot")
U1 = OperatorAtom("U1", kind="dilation_operator")
U2 = OperatorAtom("U2", kind="dilation_operator")
Ghat1 = OperatorAtom("Ghat1", kind="normalized_coefficient")
Ghat2 = OperatorAtom("Ghat2", kind="normalized_coefficient")
C2 = OperatorAtom("C2", kind="complete_first_pivot_correction")
C2_normalized = OperatorAtom(
    "C2_normalized",
    kind="normalized_complete_first_pivot_correction",
)


def main_expression() -> LinearCombination:
    """Build the ordered four-term principal expansion for the current MVP."""

    return LinearCombination(
        (
            Term(1, Vtilde_alpha2 * Wminus_21 * R11 * Vtilde_alpha1 * Wplus_12),
            Term(-1, Vtilde_alpha2 * Wminus_21 * R11 * G1 * Wplus_12),
            Term(-1, G2 * Wminus_21 * R11 * Vtilde_alpha1 * Wplus_12),
            Term(1, G2 * Wminus_21 * R11 * G1 * Wplus_12),
        )
    )
