"""Kernel extraction and combined-kernel constructions for the MVP."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import sympy as sp

from .actions import (
    AppliedLinearCombination,
    AtomicAction,
    apply_atom,
    apply_linear_combination_ordered,
    mvp_atomic_rules,
)
from .blocks import a22_first_schur_correction
from .operators import (
    R11,
    LinearCombination,
    OperatorAtom,
    Product,
    Term,
)
from .substitution import collect_bound_symbols, fresh_symbol


class KernelExtractionError(ValueError):
    """Raised when a scalar action is not a supported integral action."""


class SchurFactorizationError(ValueError):
    """Raised when the canonical first Schur correction cannot be factored."""


@dataclass(frozen=True)
class KernelTerm:
    """One extracted kernel with its originating operator product."""

    coefficient: int | float | complex
    product: Product
    kernel: sp.Expr
    ordered_factors: tuple[sp.Expr, ...] | None = None

    def as_expr(self) -> sp.Expr:
        """Project this ordered kernel term to a signed SymPy expression."""

        return self.coefficient * self.kernel


@dataclass(frozen=True)
class KernelCombination:
    """Ordered extracted kernels before projection to a SymPy ``Add``."""

    terms: tuple[KernelTerm, ...]

    def as_expr(self) -> sp.Expr:
        """Project the ordered kernels to their final SymPy sum."""

        return sum((term.as_expr() for term in self.terms), sp.Integer(0))


@dataclass(frozen=True)
class FirstSchurCorrectionFactorization:
    """Specific factorization ``left * R11 * right`` of the m=2 correction."""

    correction: LinearCombination
    left: LinearCombination
    regularizer: OperatorAtom
    right: LinearCombination


R11_kernel = sp.Function("R11")


def factor_first_schur_correction(
    correction: LinearCombination | None = None,
) -> FirstSchurCorrectionFactorization:
    """Factor the canonical correction into its two sides around ``R11``."""

    correction = (
        a22_first_schur_correction() if correction is None else correction
    )
    if not isinstance(correction, LinearCombination):
        raise SchurFactorizationError("correction must be a LinearCombination.")
    if len(correction.terms) != 4:
        raise SchurFactorizationError("the first Schur correction must have four terms.")

    split_terms: list[tuple[Term, Product, Product]] = []
    for term in correction.terms:
        positions = [
            index
            for index, factor in enumerate(term.product.factors)
            if factor is R11
        ]
        if len(positions) != 1:
            raise SchurFactorizationError(
                "each correction term must contain R11 exactly once."
            )
        position = positions[0]
        prefix = Product(term.product.factors[:position])
        suffix = Product(term.product.factors[position + 1 :])
        if not prefix.factors or not suffix.factors:
            raise SchurFactorizationError(
                "R11 must have compatible nonempty factors on both sides."
            )
        split_terms.append((term, prefix, suffix))

    prefixes = _unique_products(item[1] for item in split_terms)
    suffixes = _unique_products(item[2] for item in split_terms)
    if len(prefixes) != 2 or len(suffixes) != 2:
        raise SchurFactorizationError(
            "the correction must have two unique left and right products."
        )
    if correction.terms[0].coefficient != 1:
        raise SchurFactorizationError(
            "the canonical correction must have leading coefficient one."
        )

    coefficient_by_pair = {
        (prefix, suffix): term.coefficient
        for term, prefix, suffix in split_terms
    }
    expected_pairs = tuple(
        (prefix, suffix) for prefix in prefixes for suffix in suffixes
    )
    actual_pairs = tuple((prefix, suffix) for _, prefix, suffix in split_terms)
    if actual_pairs != expected_pairs or set(coefficient_by_pair) != set(expected_pairs):
        raise SchurFactorizationError(
            "correction terms do not form the ordered Cartesian product."
        )

    left = LinearCombination(
        tuple(
            Term(coefficient_by_pair[(prefix, suffixes[0])], prefix)
            for prefix in prefixes
        )
    )
    right = LinearCombination(
        tuple(
            Term(coefficient_by_pair[(prefixes[0], suffix)], suffix)
            for suffix in suffixes
        )
    )
    factorization = FirstSchurCorrectionFactorization(
        correction=correction,
        left=left,
        regularizer=R11,
        right=right,
    )
    if _expand_first_schur_factorization(factorization) != correction:
        raise SchurFactorizationError(
            "derived factors do not reconstruct the canonical correction."
        )
    return factorization


def extract_integral_kernel(
    applied_expression: sp.Expr,
    output_variable: sp.Symbol,
    input_function: sp.FunctionClass,
    input_variable: sp.Symbol,
) -> sp.Expr:
    """Extract ``K(x, y)`` from a supported formal action ``(Tf)(x)``.

    The MVP shape is an unevaluated nested integral whose innermost integral
    has integrand ``kernel_part * f(y)`` over the input variable. This function
    removes only that final input integral, keeps the remaining integral
    variables bound, and leaves the input variable free in the extracted kernel.
    """

    if not isinstance(applied_expression, sp.Expr):
        raise KernelExtractionError("applied_expression must be a SymPy expression.")
    if not isinstance(output_variable, sp.Symbol):
        raise KernelExtractionError("output_variable must be a SymPy Symbol.")
    if not isinstance(input_variable, sp.Symbol):
        raise KernelExtractionError("input_variable must be a SymPy Symbol.")

    input_atom = input_function(input_variable)
    kernel, found = _remove_input_integral(
        applied_expression,
        input_atom,
        input_variable,
    )
    if not found:
        raise KernelExtractionError(
            "Could not find a final Integral containing the input function."
        )
    if kernel.has(input_atom):
        raise KernelExtractionError("The extracted kernel still contains f(y).")
    return _canonicalize_remaining_integrals(kernel)


def extract_applied_kernels(
    applied_combination: AppliedLinearCombination,
    output_variable: sp.Symbol,
    input_function: sp.FunctionClass,
    input_variable: sp.Symbol,
) -> KernelCombination:
    """Extract ordered kernels from an ordered applied linear combination."""

    if not isinstance(applied_combination, AppliedLinearCombination):
        raise TypeError(
            "extract_applied_kernels expects an AppliedLinearCombination, got "
            f"{type(applied_combination).__name__}."
        )

    return KernelCombination(
        tuple(
            KernelTerm(
                coefficient=term.coefficient,
                product=term.product,
                kernel=extract_integral_kernel(
                    term.expression,
                    output_variable,
                    input_function,
                    input_variable,
                ),
            )
            for term in applied_combination.result_terms
        )
    )


def m12_kernel_combination(
    left_variable: sp.Symbol,
    input_variable: sp.Symbol,
    *,
    rules: Mapping[OperatorAtom, AtomicAction] | None = None,
) -> KernelCombination:
    """Derive the ordered M12 kernel terms from the right Schur factor."""

    factorization = factor_first_schur_correction()
    return _kernel_combination_from_factor(
        factorization.right,
        left_variable,
        input_variable,
        rules,
    )


def m12_kernel(
    left_variable: sp.Symbol,
    input_variable: sp.Symbol,
    *,
    rules: Mapping[OperatorAtom, AtomicAction] | None = None,
) -> sp.Expr:
    """Project the ordered M12 kernel combination to a scalar expression."""

    return m12_kernel_combination(
        left_variable,
        input_variable,
        rules=rules,
    ).as_expr()


def m21_kernel_combination(
    output_variable: sp.Symbol,
    right_variable: sp.Symbol,
    *,
    rules: Mapping[OperatorAtom, AtomicAction] | None = None,
) -> KernelCombination:
    """Derive the ordered M21 kernel terms from the left Schur factor."""

    factorization = factor_first_schur_correction()
    return _kernel_combination_from_factor(
        factorization.left,
        output_variable,
        right_variable,
        rules,
    )


def m21_kernel(
    output_variable: sp.Symbol,
    right_variable: sp.Symbol,
    *,
    rules: Mapping[OperatorAtom, AtomicAction] | None = None,
) -> sp.Expr:
    """Project the ordered M21 kernel combination to a scalar expression."""

    return m21_kernel_combination(
        output_variable,
        right_variable,
        rules=rules,
    ).as_expr()


def c22_integrand(
    output_variable: sp.Symbol,
    outer_variable: sp.Symbol,
    middle_variable: sp.Symbol,
    input_variable: sp.Symbol,
    *,
    rules: Mapping[OperatorAtom, AtomicAction] | None = None,
) -> sp.Expr:
    """Return the scalar integrand ``M21 * R11 * M12``."""

    return (
        m21_kernel(output_variable, outer_variable, rules=rules)
        * R11_kernel(outer_variable, middle_variable)
        * m12_kernel(middle_variable, input_variable, rules=rules)
    )


def combined_kernel_c22(
    output_variable: sp.Symbol,
    input_variable: sp.Symbol,
    *,
    outer_variable: sp.Symbol | None = None,
    middle_variable: sp.Symbol | None = None,
    rules: Mapping[OperatorAtom, AtomicAction] | None = None,
) -> sp.Integral:
    """Derive the combined kernel ``C22(x, y)`` using the selected rules."""

    outer_variable, middle_variable = _resolve_c22_variables(
        output_variable,
        input_variable,
        outer_variable,
        middle_variable,
    )
    return sp.Integral(
        m21_kernel(output_variable, outer_variable, rules=rules)
        * sp.Integral(
            R11_kernel(outer_variable, middle_variable)
            * m12_kernel(middle_variable, input_variable, rules=rules),
            (middle_variable, 0, sp.oo),
        ),
        (outer_variable, 0, sp.oo),
    )


def apply_combined_kernel_c22(
    output_variable: sp.Symbol,
    input_function: sp.FunctionClass,
    input_variable: sp.Symbol,
    *,
    outer_variable: sp.Symbol | None = None,
    middle_variable: sp.Symbol | None = None,
    rules: Mapping[OperatorAtom, AtomicAction] | None = None,
) -> sp.Integral:
    """Return the selected compact action ``Integral(C22(x,y) f(y), dy)``."""

    kernel = combined_kernel_c22(
        output_variable,
        input_variable,
        outer_variable=outer_variable,
        middle_variable=middle_variable,
        rules=rules,
    )
    return sp.Integral(kernel * input_function(input_variable), (input_variable, 0, sp.oo))


def _unique_products(products: Iterable[Product]) -> tuple[Product, ...]:
    unique: list[Product] = []
    for product in products:
        if product not in unique:
            unique.append(product)
    return tuple(unique)


def _expand_first_schur_factorization(
    factorization: FirstSchurCorrectionFactorization,
) -> LinearCombination:
    return LinearCombination(
        tuple(
            Term(
                left.coefficient * right.coefficient,
                Product(
                    left.product.factors
                    + (factorization.regularizer,)
                    + right.product.factors
                ),
            )
            for left in factorization.left.terms
            for right in factorization.right.terms
        )
    )


def _kernel_expression_from_factor(
    factor: LinearCombination,
    output_variable: sp.Symbol,
    input_variable: sp.Symbol,
    rules: Mapping[OperatorAtom, AtomicAction] | None,
) -> sp.Expr:
    return _kernel_combination_from_factor(
        factor,
        output_variable,
        input_variable,
        rules,
    ).as_expr()


def _kernel_combination_from_factor(
    factor: LinearCombination,
    output_variable: sp.Symbol,
    input_variable: sp.Symbol,
    rules: Mapping[OperatorAtom, AtomicAction] | None,
) -> KernelCombination:
    input_function = sp.Function("_schur_kernel_input")
    resolved_rules = mvp_atomic_rules() if rules is None else rules
    applied = apply_linear_combination_ordered(
        factor,
        input_function(output_variable),
        output_variable,
        resolved_rules,
        integration_variable=input_variable,
    )
    extracted = extract_applied_kernels(
        applied,
        output_variable,
        input_function,
        input_variable,
    )
    ordered_terms: list[KernelTerm] = []
    for term in extracted.terms:
        leading_scalar = apply_atom(
            term.product.factors[0],
            sp.Integer(1),
            output_variable,
            resolved_rules,
        )
        remaining_kernel = term.kernel.extract_multiplicatively(leading_scalar)
        if remaining_kernel is None:
            raise KernelExtractionError(
                "could not preserve the leading scalar action in kernel order."
            )
        ordered_terms.append(
            KernelTerm(
                coefficient=term.coefficient,
                product=term.product,
                kernel=term.kernel,
                ordered_factors=(leading_scalar, remaining_kernel),
            )
        )
    return KernelCombination(tuple(ordered_terms))


def _resolve_c22_variables(
    output_variable: sp.Symbol,
    input_variable: sp.Symbol,
    outer_variable: sp.Symbol | None,
    middle_variable: sp.Symbol | None,
) -> tuple[sp.Symbol, sp.Symbol]:
    avoid = (
        output_variable.free_symbols
        | input_variable.free_symbols
        | output_variable.atoms(sp.Symbol)
        | input_variable.atoms(sp.Symbol)
        | collect_bound_symbols(output_variable)
        | collect_bound_symbols(input_variable)
    )
    outer_variable = _resolve_one_c22_variable(
        "outer_variable",
        outer_variable,
        avoid,
    )
    avoid.add(outer_variable)
    middle_variable = _resolve_one_c22_variable(
        "middle_variable",
        middle_variable,
        avoid,
    )
    return outer_variable, middle_variable


def _resolve_one_c22_variable(
    name: str,
    variable: sp.Symbol | None,
    avoid: set[sp.Symbol],
) -> sp.Symbol:
    if variable is None:
        return fresh_symbol(avoid)
    if not isinstance(variable, sp.Symbol):
        raise TypeError(f"{name} must be a SymPy Symbol.")
    if variable in avoid:
        raise ValueError(f"{name} collides with an existing symbol.")
    return variable


def _remove_input_integral(
    expression: sp.Expr,
    input_atom: sp.Expr,
    input_variable: sp.Symbol,
) -> tuple[sp.Expr, bool]:
    if isinstance(expression, sp.Integral):
        new_function, found = _remove_input_integral(
            expression.function,
            input_atom,
            input_variable,
        )
        if found:
            return sp.Integral(new_function, *expression.limits), True

        if tuple(expression.variables) == (input_variable,) and expression.function.has(
            input_atom
        ):
            return _remove_exact_input_factor(expression.function, input_atom), True
        return expression, False

    if not expression.args:
        return expression, False

    new_args: list[sp.Expr] = []
    found_count = 0
    for arg in expression.args:
        new_arg, found = _remove_input_integral(arg, input_atom, input_variable)
        new_args.append(new_arg)
        if found:
            found_count += 1

    if found_count > 1:
        raise KernelExtractionError("Found multiple candidate input integrals.")
    if found_count == 0:
        return expression, False
    return expression.func(*new_args), True


def _remove_exact_input_factor(integrand: sp.Expr, input_atom: sp.Expr) -> sp.Expr:
    factors = sp.Mul.make_args(integrand)
    matches = [index for index, factor in enumerate(factors) if factor == input_atom]
    if len(matches) != 1:
        raise KernelExtractionError(
            "The final integrand must contain exactly one explicit input factor."
        )

    remaining_factors = [
        factor for index, factor in enumerate(factors) if index != matches[0]
    ]
    return sp.Mul(*remaining_factors) if remaining_factors else sp.Integer(1)


def _canonicalize_remaining_integrals(expression: sp.Expr) -> sp.Expr:
    split = _split_single_integral_factor(expression)
    if split is None:
        return expression

    outside, outer_integral = split
    outer_function = outside * outer_integral.function
    return sp.Integral(outer_function, *outer_integral.limits)


def _split_single_integral_factor(
    expression: sp.Expr,
) -> tuple[sp.Expr, sp.Integral] | None:
    factors = sp.Mul.make_args(expression)
    integral_factors = [
        factor for factor in factors if isinstance(factor, sp.Integral)
    ]
    if len(integral_factors) != 1:
        return None

    integral = integral_factors[0]
    outside_factors = [factor for factor in factors if factor is not integral]
    outside = sp.Mul(*outside_factors) if outside_factors else sp.Integer(1)
    return outside, integral
