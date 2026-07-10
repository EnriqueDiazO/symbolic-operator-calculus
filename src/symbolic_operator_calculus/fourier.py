"""Fourier kernels under the MVP convention.

The project convention is explicit here:

- forward exponential: ``exp(-I*t*lambda)``
- inverse exponential: ``exp(+I*t*lambda)``
- inverse prefactor: ``1/(2*pi)``

The inverse kernels in this module are built from unevaluated SymPy integrals
first. Evaluation uses ``conds="none"`` only after checking that the decay
parameter ``d`` is declared real and positive, matching the MVP hypothesis
``d = c2 - c1 > 0``.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


class FourierEvaluationError(ValueError):
    """Raised when the Fourier calculation assumptions are insufficient."""


class FourierScalingError(ValueError):
    """Raised when a Fourier scaling lacks a positive real scale."""


@dataclass(frozen=True)
class InverseFourierIntegral:
    """Observable inverse-Fourier integral before and after evaluation."""

    symbol_expression: sp.Expr
    integrand: sp.Expr
    integration_variable: sp.Symbol
    lower_limit: sp.Expr
    upper_limit: sp.Expr
    decay_parameter: sp.Symbol
    time_variable: sp.Expr

    @property
    def limits(self) -> tuple[sp.Symbol, sp.Expr, sp.Expr]:
        """Return the single integration limit tuple."""

        return (self.integration_variable, self.lower_limit, self.upper_limit)

    @property
    def integral(self) -> sp.Integral:
        """Return the unevaluated inverse-Fourier integral."""

        return sp.Integral(self.integrand, self.limits)

    @property
    def assumptions(self) -> dict[str, bool | None]:
        """Return the assumptions relevant to convergence."""

        return {
            "decay_parameter_is_real": self.decay_parameter.is_real,
            "decay_parameter_is_positive": self.decay_parameter.is_positive,
        }

    def evaluate_default(self) -> sp.Expr:
        """Evaluate with SymPy's default condition handling."""

        return self.integral.doit()

    def evaluate(self) -> sp.Expr:
        """Evaluate after checking the explicit positive-decay hypothesis."""

        if not self.decay_parameter.is_real or not self.decay_parameter.is_positive:
            raise FourierEvaluationError(
                "Fourier kernel evaluation requires a real positive decay parameter."
            )
        return self.integral.doit(conds="none")


def frequency_variable() -> sp.Symbol:
    """Return the real Fourier frequency variable named ``lambda``."""

    return sp.Symbol("lambda", real=True)


def time_variable() -> sp.Symbol:
    """Return the real transform variable named ``t``."""

    return sp.Symbol("t", real=True)


def positive_decay_symbol() -> sp.Symbol:
    """Return the real positive decay parameter ``d``."""

    return sp.Symbol("d", real=True, positive=True)


def forward_exponential(
    time: sp.Expr | None = None,
    frequency: sp.Symbol | None = None,
) -> sp.Expr:
    """Return the forward-transform exponential ``exp(-I*t*lambda)``."""

    time = time if time is not None else time_variable()
    frequency = frequency or frequency_variable()
    return sp.exp(-sp.I * time * frequency)


def inverse_exponential(
    time: sp.Expr | None = None,
    frequency: sp.Symbol | None = None,
) -> sp.Expr:
    """Return the inverse-transform exponential ``exp(+I*t*lambda)``."""

    time = time if time is not None else time_variable()
    frequency = frequency or frequency_variable()
    return sp.exp(sp.I * time * frequency)


def inverse_prefactor() -> sp.Expr:
    """Return the inverse Fourier prefactor ``1/(2*pi)``."""

    return sp.Rational(1, 2) / sp.pi


def scaled_convolution_kernel(
    kernel: sp.Expr,
    time: sp.Symbol,
    scale: sp.Expr,
) -> sp.Expr:
    """Return ``a*h(a*t)`` for a declared positive real scale ``a``."""

    _validate_fourier_scaling(kernel, time, scale)
    return scale * kernel.xreplace({time: scale * time})


def scaled_fourier_symbol(
    symbol: sp.Expr,
    frequency: sp.Symbol,
    scale: sp.Expr,
) -> sp.Expr:
    """Return ``b(lambda/a)`` paired with ``a*h(a*t)`` by this convention."""

    _validate_fourier_scaling(symbol, frequency, scale)
    return symbol.xreplace({frequency: frequency / scale})


def _validate_fourier_scaling(
    expression: sp.Expr,
    variable: sp.Symbol,
    scale: sp.Expr,
) -> None:
    if not isinstance(expression, sp.Expr):
        raise TypeError("expression must be a SymPy expression.")
    if not isinstance(variable, sp.Symbol):
        raise TypeError("the scaling variable must be a SymPy Symbol.")
    if variable not in expression.free_symbols:
        raise FourierScalingError("the expression must depend on its scaling variable.")
    if scale.is_real is not True or scale.is_positive is not True:
        raise FourierScalingError("Fourier scaling requires a positive real scale.")


def bplus_symbol(
    decay: sp.Symbol | None = None,
    frequency: sp.Symbol | None = None,
) -> sp.Expr:
    """Return ``exp(-d*lambda)`` on the positive half-line."""

    decay = decay or positive_decay_symbol()
    frequency = frequency or frequency_variable()
    return sp.exp(-decay * frequency)


def bminus_symbol(
    decay: sp.Symbol | None = None,
    frequency: sp.Symbol | None = None,
) -> sp.Expr:
    """Return ``exp(d*lambda)`` on the negative half-line."""

    decay = decay or positive_decay_symbol()
    frequency = frequency or frequency_variable()
    return sp.exp(decay * frequency)


def positive_inverse_integral(
    time: sp.Expr | None = None,
    decay: sp.Symbol | None = None,
    frequency: sp.Symbol | None = None,
) -> InverseFourierIntegral:
    """Build the unevaluated positive inverse-Fourier integral."""

    time = time if time is not None else time_variable()
    decay = decay or positive_decay_symbol()
    frequency = frequency or frequency_variable()
    symbol_expression = bplus_symbol(decay, frequency)
    return InverseFourierIntegral(
        symbol_expression=symbol_expression,
        integrand=symbol_expression * inverse_exponential(time, frequency),
        integration_variable=frequency,
        lower_limit=sp.Integer(0),
        upper_limit=sp.oo,
        decay_parameter=decay,
        time_variable=time,
    )


def negative_inverse_integral(
    time: sp.Expr | None = None,
    decay: sp.Symbol | None = None,
    frequency: sp.Symbol | None = None,
) -> InverseFourierIntegral:
    """Build the unevaluated negative inverse-Fourier integral."""

    time = time if time is not None else time_variable()
    decay = decay or positive_decay_symbol()
    frequency = frequency or frequency_variable()
    symbol_expression = bminus_symbol(decay, frequency)
    return InverseFourierIntegral(
        symbol_expression=symbol_expression,
        integrand=symbol_expression * inverse_exponential(time, frequency),
        integration_variable=frequency,
        lower_limit=-sp.oo,
        upper_limit=sp.Integer(0),
        decay_parameter=decay,
        time_variable=time,
    )


def kplus_kernel(
    time: sp.Expr | None = None,
    decay: sp.Symbol | None = None,
    frequency: sp.Symbol | None = None,
) -> sp.Expr:
    """Evaluate ``K^+_{1,2}(t)`` from its inverse-Fourier integral."""

    calculation = positive_inverse_integral(time, decay, frequency)
    return inverse_prefactor() * calculation.evaluate()


def kminus_kernel(
    time: sp.Expr | None = None,
    decay: sp.Symbol | None = None,
    frequency: sp.Symbol | None = None,
) -> sp.Expr:
    """Evaluate ``K^-_{2,1}(t)`` from its inverse-Fourier integral."""

    calculation = negative_inverse_integral(time, decay, frequency)
    return inverse_prefactor() * calculation.evaluate()


def localized_lplus_kernel(
    output_variable: sp.Expr,
    input_variable: sp.Expr,
    *,
    decay: sp.Symbol | None = None,
) -> sp.Expr:
    """Return formal ``L^+_{1,2}(x,y)`` using ``Kplus(x-y)``."""

    chi_infinity = sp.Function("chi_infinity")
    return (
        chi_infinity(output_variable)
        * kplus_kernel(output_variable - input_variable, decay=decay)
        * chi_infinity(input_variable)
    )


def localized_lminus_kernel(
    output_variable: sp.Expr,
    input_variable: sp.Expr,
    *,
    decay: sp.Symbol | None = None,
) -> sp.Expr:
    """Return formal ``L^-_{2,1}(x,y)`` using ``Kminus(x-y)``."""

    chi_infinity = sp.Function("chi_infinity")
    return (
        chi_infinity(output_variable)
        * kminus_kernel(output_variable - input_variable, decay=decay)
        * chi_infinity(input_variable)
    )
