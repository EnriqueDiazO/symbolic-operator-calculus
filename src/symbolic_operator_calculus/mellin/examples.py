"""Controlled scalar Mellin examples with explicit P0-B metadata."""

from __future__ import annotations

from dataclasses import dataclass
from typing import final

import sympy as sp

from ..conditional import ConditionalIdentity
from ..domains import AssumptionContext, ComplexDomain
from ..singularities import detect_singularities
from .domains import MellinSymbolDomain
from .expressions import MellinExpression
from .symbols import MellinSymbol, MellinSymbolDependency
from .variables import (
    MellinVariable,
    MellinVariableRole,
    mellin_frequency,
    mellin_parameter,
    output_spatial_variable,
    relative_multiplicative_variable,
)


@final
@dataclass(frozen=True)
class DiagonalMellinSymbolFamily:
    """The scalar symbols z, coth(z), p_-, p_+ and two checked identities."""

    z: MellinSymbol
    s: MellinSymbol
    p_minus: MellinSymbol
    p_plus: MellinSymbol
    product_identity: ConditionalIdentity
    derivative_identity: ConditionalIdentity


def _validate_symbol(value: object, name: str) -> sp.Symbol:
    if not isinstance(value, sp.Symbol):
        raise TypeError(f"{name} must be an existing SymPy Symbol.")
    return value


def _make_domain(
    *,
    frequency: MellinVariable,
    expression: sp.Expr,
    context: AssumptionContext,
    spatial_variables: tuple[MellinVariable, ...] = (),
    relative_variable: MellinVariable | None = None,
    parameters: tuple[MellinVariable, ...] = (),
    description: str,
) -> MellinSymbolDomain:
    complex_domain = ComplexDomain.real_line(
        frequency.symbol,
        assumption_context=context,
        description="real Mellin-frequency line under explicit hypotheses",
    )
    singular_set = detect_singularities(
        expression,
        frequency.symbol,
        assumptions=context,
        domain=complex_domain,
    )
    return MellinSymbolDomain(
        frequency=frequency,
        complex_domain=complex_domain,
        spatial_variables=spatial_variables,
        relative_variable=relative_variable,
        parameters=parameters,
        assumption_context=context,
        singular_set=singular_set,
        description=description,
    )


def _make_symbol(
    expression: sp.Expr,
    *,
    frequency: MellinVariable,
    context: AssumptionContext,
    variables: tuple[MellinVariable, ...],
    parameters: tuple[MellinVariable, ...],
    dependency: MellinSymbolDependency,
    spatial_variables: tuple[MellinVariable, ...] = (),
    relative_variable: MellinVariable | None = None,
    description: str,
) -> MellinSymbol:
    domain = _make_domain(
        frequency=frequency,
        expression=expression,
        context=context,
        spatial_variables=spatial_variables,
        relative_variable=relative_variable,
        parameters=parameters,
        description=description,
    )
    formal = MellinExpression(
        expression=expression,
        domain=domain,
        variables=variables,
        parameters=parameters,
        description=description,
    )
    return MellinSymbol(formal, dependency, description)


def build_diagonal_mellin_symbols(
    lam: sp.Symbol,
    kappa: sp.Symbol,
) -> DiagonalMellinSymbolFamily:
    """Build z, coth(z), p_-/p_+ on lambda real with 0 < kappa < 1."""

    lam = _validate_symbol(lam, "lam")
    kappa = _validate_symbol(kappa, "kappa")
    frequency = mellin_frequency(lam)
    parameter = mellin_parameter(kappa)
    context = AssumptionContext((kappa > 0, kappa < 1))
    variables = (frequency,)
    parameters = (parameter,)
    z_expr = sp.pi * (lam + sp.I * kappa)
    s_expr = sp.coth(z_expr)
    p_minus_expr = (1 - s_expr) / 2
    p_plus_expr = (1 + s_expr) / 2

    def make(expression: sp.Expr, description: str) -> MellinSymbol:
        return _make_symbol(
            expression,
            frequency=frequency,
            context=context,
            variables=variables,
            parameters=parameters,
            dependency=MellinSymbolDependency.FREQUENCY_ONLY,
            description=description,
        )

    z = make(z_expr, "z(lambda) = pi(lambda + i kappa)")
    s = make(s_expr, "diagonal coth symbol")
    p_minus = make(p_minus_expr, "negative diagonal scalar projection symbol")
    p_plus = make(p_plus_expr, "positive diagonal scalar projection symbol")

    product_expression = p_minus.expression * p_plus.expression
    product = MellinSymbol(
        product_expression,
        MellinSymbolDependency.FREQUENCY_ONLY,
        "scalar product p_- p_+",
    )
    product_rhs = make(
        -1 / (4 * sp.sinh(z_expr) ** 2),
        "conditional scalar product formula",
    )
    product_identity = product.conditional_identity(
        product_rhs,
        description="conditional scalar identity for p_- p_+",
    ).symbolically_check()

    derivative_expression = s.expression.differentiate(lam)
    derivative = MellinSymbol(
        derivative_expression,
        MellinSymbolDependency.FREQUENCY_ONLY,
        "scalar derivative of coth symbol",
    )
    derivative_rhs = make(
        -sp.pi / sp.sinh(z_expr) ** 2,
        "conditional scalar derivative formula",
    )
    derivative_identity = derivative.conditional_identity(
        derivative_rhs,
        description="conditional scalar derivative identity",
    ).symbolically_check()

    return DiagonalMellinSymbolFamily(
        z=z,
        s=s,
        p_minus=p_minus,
        p_plus=p_plus,
        product_identity=product_identity,
        derivative_identity=derivative_identity,
    )


def build_exponential_sinh_symbol(
    lam: sp.Symbol,
    beta: sp.Symbol,
    kappa: sp.Symbol,
) -> MellinSymbol:
    """Build exp(beta(lambda+i kappa))/sinh(pi(lambda+i kappa))."""

    lam = _validate_symbol(lam, "lam")
    beta = _validate_symbol(beta, "beta")
    kappa = _validate_symbol(kappa, "kappa")
    frequency = mellin_frequency(lam)
    beta_parameter = mellin_parameter(beta)
    kappa_parameter = mellin_parameter(kappa)
    context = AssumptionContext(
        (beta > -sp.pi, beta < sp.pi, kappa > 0, kappa < 1)
    )
    shifted = lam + sp.I * kappa
    expression = sp.exp(beta * shifted) / sp.sinh(sp.pi * shifted)
    return _make_symbol(
        expression,
        frequency=frequency,
        context=context,
        variables=(frequency,),
        parameters=(beta_parameter, kappa_parameter),
        dependency=MellinSymbolDependency.FREQUENCY_ONLY,
        description="exponential-over-sinh scalar Mellin symbol",
    )


def build_dilation_multiplier(
    lam: sp.Symbol,
    gamma: sp.Symbol,
    *,
    inverse: bool = False,
) -> MellinSymbol:
    """Build gamma**(+/- i lambda) on the declared gamma > 0 branch."""

    lam = _validate_symbol(lam, "lam")
    gamma = _validate_symbol(gamma, "gamma")
    frequency = mellin_frequency(lam)
    parameter = mellin_parameter(gamma)
    context = AssumptionContext((gamma > 0,))
    exponent = -sp.I * lam if inverse else sp.I * lam
    expression = sp.Pow(gamma, exponent, evaluate=False)
    return _make_symbol(
        expression,
        frequency=frequency,
        context=context,
        variables=(frequency,),
        parameters=(parameter,),
        dependency=MellinSymbolDependency.FREQUENCY_ONLY,
        description=(
            "positive-base inverse-dilation multiplier"
            if inverse
            else "positive-base scalar dilation multiplier"
        ),
    )


def build_space_frequency_symbol(
    lam: sp.Symbol,
    x: sp.Symbol,
    omega: sp.Expr,
    beta: sp.Symbol,
    kappa: sp.Symbol,
) -> MellinSymbol:
    """Build exp(i omega(x) lambda) h_beta(lambda) as a scalar symbol."""

    lam = _validate_symbol(lam, "lam")
    x = _validate_symbol(x, "x")
    beta = _validate_symbol(beta, "beta")
    kappa = _validate_symbol(kappa, "kappa")
    if not isinstance(omega, sp.Expr):
        raise TypeError("omega must be an existing SymPy scalar expression.")
    if not omega.free_symbols.issubset({x}):
        raise ValueError("omega may depend only on the declared spatial variable.")
    frequency = mellin_frequency(lam)
    spatial = output_spatial_variable(x)
    beta_parameter = mellin_parameter(beta)
    kappa_parameter = mellin_parameter(kappa)
    context = AssumptionContext(
        (beta > -sp.pi, beta < sp.pi, kappa > 0, kappa < 1)
    )
    shifted = lam + sp.I * kappa
    h_beta = sp.exp(beta * shifted) / sp.sinh(sp.pi * shifted)
    expression = sp.exp(sp.I * omega * lam) * h_beta
    return _make_symbol(
        expression,
        frequency=frequency,
        context=context,
        variables=(frequency, spatial),
        parameters=(beta_parameter, kappa_parameter),
        dependency=MellinSymbolDependency.SPACE_FREQUENCY,
        spatial_variables=(spatial,),
        description="space-frequency scalar Mellin symbol",
    )


def build_relative_frequency_symbol(
    lam: sp.Symbol,
    relative: sp.Symbol,
    expression: sp.Expr,
    *,
    parameters: tuple[MellinVariable, ...] = (),
    assumption_context: AssumptionContext = AssumptionContext(),
) -> MellinSymbol:
    """Wrap an explicitly supplied q(r, lambda) without constructing a transform."""

    lam = _validate_symbol(lam, "lam")
    relative = _validate_symbol(relative, "relative")
    if not isinstance(expression, sp.Expr):
        raise TypeError("expression must be an existing SymPy scalar expression.")
    if not all(
        isinstance(item, MellinVariable)
        and item.role is MellinVariableRole.PARAMETER
        for item in parameters
    ):
        raise TypeError("parameters must contain PARAMETER MellinVariable objects.")
    frequency = mellin_frequency(lam)
    relative_variable = relative_multiplicative_variable(relative)
    allowed = {lam, relative, *(item.symbol for item in parameters)}
    if not expression.free_symbols.issubset(allowed):
        raise ValueError("relative expression contains undeclared free symbols.")
    if lam not in expression.free_symbols or relative not in expression.free_symbols:
        raise ValueError("relative-frequency example must depend on r and lambda.")
    return _make_symbol(
        expression,
        frequency=frequency,
        context=assumption_context,
        variables=(frequency, relative_variable),
        parameters=parameters,
        dependency=MellinSymbolDependency.RELATIVE_FREQUENCY,
        relative_variable=relative_variable,
        description="relative-frequency scalar Mellin symbol",
    )
