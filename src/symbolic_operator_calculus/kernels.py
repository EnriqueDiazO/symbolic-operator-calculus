"""Kernel extraction and combined-kernel constructions for the MVP."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .actions import AppliedLinearCombination
from .operators import Product


class KernelExtractionError(ValueError):
    """Raised when a scalar action is not a supported integral action."""


@dataclass(frozen=True)
class KernelTerm:
    """One extracted kernel with its originating operator product."""

    coefficient: int | float | complex
    product: Product
    kernel: sp.Expr

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


G1_scalar = sp.Function("G1")
G2_scalar = sp.Function("G2")
rho1 = sp.Function("rho1")
rho2 = sp.Function("rho2")
gamma1 = sp.Symbol("gamma1")
gamma2 = sp.Symbol("gamma2")
Lplus_12_kernel = sp.Function("Lplus_12")
Lminus_21_kernel = sp.Function("Lminus_21")
R11_kernel = sp.Function("R11")


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


def m12_kernel(left_variable: sp.Symbol, input_variable: sp.Symbol) -> sp.Expr:
    """Return ``M_{1,2}(left_variable, input_variable)``."""

    return rho1(left_variable) * Lplus_12_kernel(gamma1 * left_variable, input_variable) - (
        G1_scalar(left_variable) * Lplus_12_kernel(left_variable, input_variable)
    )


def m21_kernel(output_variable: sp.Symbol, right_variable: sp.Symbol) -> sp.Expr:
    """Return ``M_{2,1}(output_variable, right_variable)``."""

    return rho2(output_variable) * Lminus_21_kernel(
        gamma2 * output_variable,
        right_variable,
    ) - (
        G2_scalar(output_variable)
        * Lminus_21_kernel(output_variable, right_variable)
    )


def c22_integrand(
    output_variable: sp.Symbol,
    outer_variable: sp.Symbol,
    middle_variable: sp.Symbol,
    input_variable: sp.Symbol,
) -> sp.Expr:
    """Return the scalar integrand ``M21 * R11 * M12``."""

    return (
        m21_kernel(output_variable, outer_variable)
        * R11_kernel(outer_variable, middle_variable)
        * m12_kernel(middle_variable, input_variable)
    )


def combined_kernel_c22(
    output_variable: sp.Symbol,
    input_variable: sp.Symbol,
    *,
    outer_variable: sp.Symbol | None = None,
    middle_variable: sp.Symbol | None = None,
) -> sp.Integral:
    """Return the formal combined kernel ``C22(x, y)``."""

    outer_variable = outer_variable or sp.Symbol("u")
    middle_variable = middle_variable or sp.Symbol("v")
    return sp.Integral(
        m21_kernel(output_variable, outer_variable)
        * sp.Integral(
            R11_kernel(outer_variable, middle_variable)
            * m12_kernel(middle_variable, input_variable),
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
) -> sp.Integral:
    """Return the formal action ``Integral(C22(x, y) * f(y), dy)``."""

    kernel = combined_kernel_c22(
        output_variable,
        input_variable,
        outer_variable=outer_variable,
        middle_variable=middle_variable,
    )
    return sp.Integral(kernel * input_function(input_variable), (input_variable, 0, sp.oo))


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
