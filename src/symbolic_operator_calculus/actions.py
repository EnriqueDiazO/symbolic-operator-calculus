"""Atomic actions of single noncommutative operator atoms.

The operator AST remains structural and independent from SymPy.  SymPy appears
here only as the scalar language for functions, symbols, and formal integrals.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import TypeAlias

import sympy as sp

from .operators import (
    G1,
    G2,
    LinearCombination,
    OperatorAtom,
    Product,
    R11,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
)

ScalarFunction: TypeAlias = Callable[..., sp.Expr]


class AtomicActionError(Exception):
    """Base class for atomic action errors."""


class MissingAtomicActionError(AtomicActionError, KeyError):
    """Raised when no action is supplied for an operator atom."""


class UnsupportedAtomicExpressionError(AtomicActionError, TypeError):
    """Raised when the atomic API receives a non-atomic operator expression."""


class UnsupportedProductExpressionError(AtomicActionError, TypeError):
    """Raised when the product API receives a non-product expression."""


class UnsupportedLinearCombinationExpressionError(AtomicActionError, TypeError):
    """Raised when the linear-combination API receives a non-combination."""


class ProductApplicationError(AtomicActionError):
    """Raised when a factor fails during product application."""

    def __init__(
        self,
        product: Product,
        factor: object,
        factor_position: int,
        factor_count: int,
    ) -> None:
        self.product = product
        self.factor = factor
        self.factor_position = factor_position
        self.factor_count = factor_count
        super().__init__(
            "Failed to apply factor "
            f"{factor_position}/{factor_count} of product {product!r}: "
            f"{factor!r}."
        )


class LinearCombinationApplicationError(AtomicActionError):
    """Raised when a term fails during linear-combination application."""

    def __init__(
        self,
        combination: LinearCombination,
        term: object,
        term_position: int,
        term_count: int,
        coefficient: object,
    ) -> None:
        self.combination = combination
        self.term = term
        self.term_position = term_position
        self.term_count = term_count
        self.coefficient = coefficient
        super().__init__(
            "Failed to apply term "
            f"{term_position}/{term_count} with coefficient {coefficient} "
            f"of linear combination {combination!r}: {term!r}."
        )


class MissingIntegrationVariableError(AtomicActionError, ValueError):
    """Raised when an integral action has no explicit integration variable."""


class UnsafeScalarSubstitutionError(AtomicActionError, ValueError):
    """Raised when this phase cannot safely substitute inside an operand."""


@dataclass(frozen=True)
class MultiplicationAction:
    """Action of a scalar multiplication operator, q(x) -> G(x) q(x)."""

    coefficient: ScalarFunction

    def apply(
        self,
        operand: sp.Expr,
        variable: sp.Symbol,
        *,
        integration_variable: sp.Symbol | None = None,
    ) -> sp.Expr:
        return self.coefficient(variable) * operand


@dataclass(frozen=True)
class TransportedShiftAction:
    """Action of rho V_alpha with alpha(x) = gamma*x."""

    rho: ScalarFunction
    gamma: sp.Expr

    def apply(
        self,
        operand: sp.Expr,
        variable: sp.Symbol,
        *,
        integration_variable: sp.Symbol | None = None,
    ) -> sp.Expr:
        shifted = _replace_free_variable(operand, variable, self.gamma * variable)
        return self.rho(variable) * shifted


@dataclass(frozen=True)
class IntegralKernelAction:
    """Action of a formal integral kernel over (0, infinity)."""

    kernel: ScalarFunction

    def apply(
        self,
        operand: sp.Expr,
        variable: sp.Symbol,
        *,
        integration_variable: sp.Symbol | None = None,
    ) -> sp.Integral:
        return _kernel_integral(self.kernel, operand, variable, integration_variable)


@dataclass(frozen=True)
class FormalRegularizerAction:
    """Kernel action of the formal regularizer R11, without inverse rules."""

    kernel: ScalarFunction

    @property
    def is_formal_regularizer(self) -> bool:
        return True

    def apply(
        self,
        operand: sp.Expr,
        variable: sp.Symbol,
        *,
        integration_variable: sp.Symbol | None = None,
    ) -> sp.Integral:
        return _kernel_integral(self.kernel, operand, variable, integration_variable)


AtomicAction: TypeAlias = (
    MultiplicationAction
    | TransportedShiftAction
    | IntegralKernelAction
    | FormalRegularizerAction
)


G1_scalar = sp.Function("G1")
G2_scalar = sp.Function("G2")
rho1 = sp.Function("rho1")
rho2 = sp.Function("rho2")
gamma1 = sp.Symbol("gamma1")
gamma2 = sp.Symbol("gamma2")
Lplus_12_kernel = sp.Function("Lplus_12")
Lminus_21_kernel = sp.Function("Lminus_21")
R11_kernel = sp.Function("R11")


def mvp_atomic_rules() -> Mapping[OperatorAtom, AtomicAction]:
    """Return the explicit atomic action mapping for the current MVP phase."""

    return MappingProxyType({
        G1: MultiplicationAction(G1_scalar),
        G2: MultiplicationAction(G2_scalar),
        Vtilde_alpha1: TransportedShiftAction(rho1, gamma1),
        Vtilde_alpha2: TransportedShiftAction(rho2, gamma2),
        Wplus_12: IntegralKernelAction(Lplus_12_kernel),
        Wminus_21: IntegralKernelAction(Lminus_21_kernel),
        R11: FormalRegularizerAction(R11_kernel),
    })


def apply_atom(
    atom: OperatorAtom,
    operand: sp.Expr,
    variable: sp.Symbol,
    rules: Mapping[OperatorAtom, AtomicAction],
    *,
    integration_variable: sp.Symbol | None = None,
) -> sp.Expr:
    """Apply exactly one operator atom using an explicit action mapping."""

    if not isinstance(atom, OperatorAtom):
        raise UnsupportedAtomicExpressionError(
            f"apply_atom expects an OperatorAtom, got {type(atom).__name__}."
        )

    try:
        action = rules[atom]
    except KeyError as exc:
        raise MissingAtomicActionError(
            f"No atomic action registered for operator atom {atom.name!r}."
        ) from exc

    return action.apply(
        operand,
        variable,
        integration_variable=integration_variable,
    )


def apply_product(
    product: Product,
    operand: sp.Expr,
    variable: sp.Symbol,
    rules: Mapping[OperatorAtom, AtomicAction],
    *,
    integration_variable: sp.Symbol | None = None,
) -> sp.Expr:
    """Apply an ordered product using the convention ``(AB)f = A(Bf)``.

    ``Product.factors`` stores factors in written left-to-right order. Applying
    a product therefore traverses those factors from right to left and delegates
    each atomic step to ``apply_atom``. This function intentionally does not
    implement ``LinearCombination``, fresh dummy variables, substitution under
    bound ``Integral`` objects, or operator simplification.

    The existing AST permits ``Product(())`` and uses it as an identity-like
    accumulator during expansion, so applying an empty product returns the
    operand unchanged.
    """

    if not isinstance(product, Product):
        raise UnsupportedProductExpressionError(
            f"apply_product expects a Product, got {type(product).__name__}."
        )

    result = operand
    factor_count = len(product.factors)
    for zero_based_position, factor in reversed(list(enumerate(product.factors))):
        try:
            result = apply_atom(
                factor,
                result,
                variable,
                rules,
                integration_variable=integration_variable,
            )
        except AtomicActionError as exc:
            raise ProductApplicationError(
                product=product,
                factor=factor,
                factor_position=zero_based_position + 1,
                factor_count=factor_count,
            ) from exc
    return result


def apply_linear_combination(
    combination: LinearCombination,
    operand: sp.Expr,
    variable: sp.Symbol,
    rules: Mapping[OperatorAtom, AtomicAction],
    *,
    integration_variable: sp.Symbol | None = None,
) -> sp.Expr:
    """Apply an ordered linear combination term by term.

    The current AST stores each term as ``Term(coefficient, Product(...))``.
    This function therefore implements ``(sum c_i A_i)f = sum c_i A_i(f)``.
    Single-atom products delegate to ``apply_atom``; longer or empty products
    delegate to ``apply_product``. There is deliberately no generic dispatcher,
    no recursive ``LinearCombination`` application, no fresh dummy variables,
    no substitution under bound ``Integral`` objects, and no simplification.
    """

    if not isinstance(combination, LinearCombination):
        raise UnsupportedLinearCombinationExpressionError(
            "apply_linear_combination expects a LinearCombination, got "
            f"{type(combination).__name__}."
        )

    scalar_terms: list[sp.Expr] = []
    term_count = len(combination.terms)
    for zero_based_position, term in enumerate(combination.terms):
        try:
            applied = _apply_linear_combination_term(
                term.product,
                operand,
                variable,
                rules,
                integration_variable=integration_variable,
            )
        except AtomicActionError as exc:
            raise LinearCombinationApplicationError(
                combination=combination,
                term=term,
                term_position=zero_based_position + 1,
                term_count=term_count,
                coefficient=term.coefficient,
            ) from exc
        scalar_terms.append(term.coefficient * applied)
    return sum(scalar_terms, sp.Integer(0))


def _apply_linear_combination_term(
    product: Product,
    operand: sp.Expr,
    variable: sp.Symbol,
    rules: Mapping[OperatorAtom, AtomicAction],
    *,
    integration_variable: sp.Symbol | None = None,
) -> sp.Expr:
    if len(product.factors) == 1 and isinstance(product.factors[0], OperatorAtom):
        return apply_atom(
            product.factors[0],
            operand,
            variable,
            rules,
            integration_variable=integration_variable,
        )
    return apply_product(
        product,
        operand,
        variable,
        rules,
        integration_variable=integration_variable,
    )


def _kernel_integral(
    kernel: ScalarFunction,
    operand: sp.Expr,
    variable: sp.Symbol,
    integration_variable: sp.Symbol | None,
) -> sp.Integral:
    if integration_variable is None:
        raise MissingIntegrationVariableError(
            "Integral kernel actions require an explicit integration variable."
        )

    substituted = _replace_free_variable(operand, variable, integration_variable)
    return sp.Integral(
        kernel(variable, integration_variable) * substituted,
        (integration_variable, 0, sp.oo),
    )


def _replace_free_variable(
    operand: sp.Expr,
    variable: sp.Symbol,
    replacement: sp.Expr,
) -> sp.Expr:
    if operand.has(sp.Integral):
        raise UnsafeScalarSubstitutionError(
            "Phase C substitutions are limited to scalar operands without "
            "bound Integral objects."
        )
    return operand.xreplace({variable: replacement})
