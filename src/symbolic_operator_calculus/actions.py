"""Atomic actions of single noncommutative operator atoms.

The operator AST remains structural and independent from SymPy.  SymPy appears
here only as the scalar language for functions, symbols, and formal integrals.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from fractions import Fraction
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
    Scalar,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
)
from .substitution import (
    SafeSubstitutionError,
    collect_bound_symbols,
    fresh_symbol,
    substitute_free_variable,
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


class UnsafeScalarSubstitutionError(AtomicActionError, SafeSubstitutionError):
    """Raised when this phase cannot safely substitute inside an operand."""


@dataclass(frozen=True)
class AppliedTerm:
    """One ordered product application with its original scalar coefficient."""

    coefficient: Scalar
    product: Product
    expression: sp.Expr

    def as_expr(self) -> sp.Expr:
        """Project this ordered applied term to a signed SymPy expression."""

        return _sympy_coefficient(self.coefficient) * self.expression


@dataclass(frozen=True)
class AppliedLinearCombination:
    """Ordered product applications before projection to a SymPy ``Add``."""

    result_terms: tuple[AppliedTerm, ...]

    @property
    def terms(self) -> tuple[AppliedTerm, ...]:
        """Alias for callers that prefer the shorter structural name."""

        return self.result_terms

    def as_expr(self) -> sp.Expr:
        """Project the ordered applied terms to their final SymPy sum."""

        return sum((term.as_expr() for term in self.result_terms), sp.Integer(0))


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
    each atomic step to ``apply_atom``. Integral factors receive deterministic
    local dummy variables when a multi-factor product needs them. This function
    intentionally does not implement ``LinearCombination`` or operator
    simplification.

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
    generated_variables: set[sp.Symbol] = set()
    first_integral_variable = integration_variable
    used_explicit_integral_variable = False
    for zero_based_position, factor in reversed(list(enumerate(product.factors))):
        factor_integration_variable = integration_variable
        if _atom_needs_integration_variable(factor, rules):
            if first_integral_variable is not None and not used_explicit_integral_variable:
                factor_integration_variable = first_integral_variable
                used_explicit_integral_variable = True
            elif factor_count > 1:
                factor_integration_variable = _fresh_integration_variable(
                    result,
                    variable,
                    generated_variables,
                )
                generated_variables.add(factor_integration_variable)
            else:
                factor_integration_variable = None
        try:
            result = apply_atom(
                factor,
                result,
                variable,
                rules,
                integration_variable=factor_integration_variable,
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
    """Apply an ordered linear combination and return its SymPy projection.

    The current AST stores each term as ``Term(coefficient, Product(...))``.
    This function therefore implements ``(sum c_i A_i)f = sum c_i A_i(f)``.
    The ordered intermediate structure is produced by
    ``apply_linear_combination_ordered`` and then projected to a scalar SymPy
    sum. There is deliberately no generic dispatcher, no recursive
    ``LinearCombination`` application, and no simplification.
    """

    return apply_linear_combination_ordered(
        combination,
        operand,
        variable,
        rules,
        integration_variable=integration_variable,
    ).as_expr()


def apply_linear_combination_ordered(
    combination: LinearCombination,
    operand: sp.Expr,
    variable: sp.Symbol,
    rules: Mapping[OperatorAtom, AtomicAction],
    *,
    integration_variable: sp.Symbol | None = None,
) -> AppliedLinearCombination:
    """Apply a linear combination while preserving ordered applied terms.

    Each result term stores the original coefficient, the original ``Product``,
    and the scalar expression obtained by applying that product to ``operand``.
    Call ``as_expr`` only when a final SymPy ``Add`` projection is needed.
    """

    if not isinstance(combination, LinearCombination):
        raise UnsupportedLinearCombinationExpressionError(
            "apply_linear_combination expects a LinearCombination, got "
            f"{type(combination).__name__}."
        )

    applied_terms: list[AppliedTerm] = []
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
        applied_terms.append(
            AppliedTerm(
                coefficient=term.coefficient,
                product=term.product,
                expression=applied,
            )
        )
    return AppliedLinearCombination(tuple(applied_terms))


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
    try:
        return substitute_free_variable(operand, variable, replacement)
    except SafeSubstitutionError as exc:
        raise UnsafeScalarSubstitutionError(str(exc)) from exc


def _atom_needs_integration_variable(
    atom: object,
    rules: Mapping[OperatorAtom, AtomicAction],
) -> bool:
    if not isinstance(atom, OperatorAtom):
        return False
    action = rules.get(atom)
    return isinstance(action, (IntegralKernelAction, FormalRegularizerAction))


def _fresh_integration_variable(
    expression: sp.Expr,
    output_variable: sp.Symbol,
    generated_variables: set[sp.Symbol],
) -> sp.Symbol:
    avoid = (
        expression.free_symbols
        | collect_bound_symbols(expression)
        | generated_variables
        | {output_variable}
    )
    return fresh_symbol(avoid, preferred_names=("y", "v", "u", "w", "z"))


def _sympy_coefficient(coefficient: Scalar) -> sp.Expr:
    if isinstance(coefficient, Fraction):
        return sp.Rational(coefficient.numerator, coefficient.denominator)
    return sp.sympify(coefficient)
