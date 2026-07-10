"""Exact dilation factorizations for the article's off-diagonal WH kernel."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import sympy as sp

from .fourier import scaled_convolution_kernel, scaled_fourier_symbol


class RelativeWienerHopfError(ValueError):
    """Raised when the concrete relative-dilation family is ill-defined."""


chi_plus = sp.Function("chi_plus")
chi_minus = sp.Function("chi_minus")


@dataclass(frozen=True)
class DilationConjugatedWienerHopfFactorization:
    """Two exact placements of the relative dilation for one ``W_{k,j}``."""

    k: int
    j: int
    gamma_k: sp.Expr
    gamma_j: sp.Expr
    delta_c: sp.Expr
    time_variable: sp.Symbol
    frequency_variable: sp.Symbol
    output_variable: sp.Symbol
    input_variable: sp.Symbol
    original_kernel: sp.Expr
    original_symbol: sp.Expr
    relative_scale: sp.Expr
    integral_kernel: sp.Expr
    left_convolution_kernel: sp.Expr
    left_symbol: sp.Expr
    left_reconstructed_kernel: sp.Expr
    right_convolution_kernel: sp.Expr
    right_symbol: sp.Expr
    right_reconstructed_kernel: sp.Expr

    @property
    def kernels_reconstruct_exactly(self) -> bool:
        """Whether both factorizations recover the conjugated kernel."""

        return all(
            sp.simplify(candidate - self.integral_kernel) == 0
            for candidate in (
                self.left_reconstructed_kernel,
                self.right_reconstructed_kernel,
            )
        )


def article_wiener_hopf_kernel(
    time: sp.Symbol,
    delta_c: sp.Expr,
) -> sp.Expr:
    """Return ``h_{k,j}(t) = -1/(pi*i*(t-i*delta_c))``."""

    if not isinstance(time, sp.Symbol):
        raise TypeError("time must be a SymPy Symbol.")
    return -1 / (sp.pi * sp.I * (time - sp.I * delta_c))


def article_wiener_hopf_symbol(
    k: int,
    j: int,
    delta_c: sp.Expr,
    frequency: sp.Symbol,
) -> sp.Expr:
    """Return the original, non-normalized article symbol ``b_{k,j}``."""

    _validate_indices(k, j)
    if not isinstance(frequency, sp.Symbol):
        raise TypeError("frequency must be a SymPy Symbol.")
    exponential = sp.exp(delta_c * frequency)
    if k < j:
        return 2 * chi_plus(frequency) * exponential
    return -2 * chi_minus(frequency) * exponential


def factor_dilation_conjugated_wiener_hopf(
    k: int,
    j: int,
    gamma_k: sp.Expr,
    gamma_j: sp.Expr,
    c_k: sp.Expr,
    c_j: sp.Expr,
    *,
    time: sp.Symbol | None = None,
    frequency: sp.Symbol | None = None,
    output: sp.Symbol | None = None,
    input_: sp.Symbol | None = None,
) -> DilationConjugatedWienerHopfFactorization:
    """Construct and verify both exact relative-dilation factorizations."""

    _validate_indices(k, j)
    _validate_positive_gamma("gamma_k", gamma_k)
    _validate_positive_gamma("gamma_j", gamma_j)
    time = time or sp.Symbol("t", real=True)
    frequency = frequency or sp.Symbol("lambda", real=True)
    output = output or sp.Symbol("x", real=True, positive=True)
    input_ = input_ or sp.Symbol("y", real=True, positive=True)
    if not all(isinstance(variable, sp.Symbol) for variable in (time, frequency, output, input_)):
        raise TypeError("time, frequency, output, and input_ must be SymPy Symbols.")

    delta_c = c_k - c_j
    _validate_delta_c(k, j, delta_c)
    original_kernel = article_wiener_hopf_kernel(time, delta_c)
    original_symbol = article_wiener_hopf_symbol(k, j, delta_c, frequency)
    relative_scale = gamma_k / gamma_j
    integral_kernel = gamma_j * original_kernel.xreplace(
        {time: gamma_k * output - gamma_j * input_}
    )

    left_kernel = scaled_convolution_kernel(original_kernel, time, gamma_j)
    left_symbol = _scaled_article_symbol(
        k,
        j,
        delta_c,
        frequency,
        gamma_j,
        original_symbol,
    )
    left_reconstructed = left_kernel.xreplace(
        {time: relative_scale * output - input_}
    )

    right_kernel = scaled_convolution_kernel(original_kernel, time, gamma_k)
    right_symbol = _scaled_article_symbol(
        k,
        j,
        delta_c,
        frequency,
        gamma_k,
        original_symbol,
    )
    right_reconstructed = (
        right_kernel.xreplace(
            {time: output - input_ / relative_scale}
        )
        / relative_scale
    )

    result = DilationConjugatedWienerHopfFactorization(
        k=k,
        j=j,
        gamma_k=gamma_k,
        gamma_j=gamma_j,
        delta_c=delta_c,
        time_variable=time,
        frequency_variable=frequency,
        output_variable=output,
        input_variable=input_,
        original_kernel=original_kernel,
        original_symbol=original_symbol,
        relative_scale=relative_scale,
        integral_kernel=integral_kernel,
        left_convolution_kernel=left_kernel,
        left_symbol=left_symbol,
        left_reconstructed_kernel=left_reconstructed,
        right_convolution_kernel=right_kernel,
        right_symbol=right_symbol,
        right_reconstructed_kernel=right_reconstructed,
    )
    if not result.kernels_reconstruct_exactly:
        raise RelativeWienerHopfError(
            "the relative-dilation factorizations did not reconstruct the kernel."
        )
    return result


def m2_relative_wiener_hopf_factorizations(
    gamma1: sp.Expr,
    gamma2: sp.Expr,
    decay: sp.Expr,
) -> tuple[
    DilationConjugatedWienerHopfFactorization,
    DilationConjugatedWienerHopfFactorization,
]:
    """Return the original-symbol factorizations for ``W12`` and ``W21``."""

    _validate_positive_gamma("decay", decay)
    c1 = sp.Integer(0)
    c2 = decay
    return (
        factor_dilation_conjugated_wiener_hopf(1, 2, gamma1, gamma2, c1, c2),
        factor_dilation_conjugated_wiener_hopf(2, 1, gamma2, gamma1, c2, c1),
    )


def halfline_indicator_is_scale_invariant(
    side: str,
    frequency: sp.Symbol,
    scale: sp.Expr,
) -> bool:
    """Record ``chi_±(lambda/a) = chi_±(lambda)`` for positive ``a``."""

    _validate_positive_gamma("scale", scale)
    if not isinstance(frequency, sp.Symbol) or frequency.is_real is not True:
        raise RelativeWienerHopfError("frequency must be a real SymPy Symbol.")
    if side not in {"plus", "minus"}:
        raise RelativeWienerHopfError("side must be 'plus' or 'minus'.")
    return True


def _scaled_article_symbol(
    k: int,
    j: int,
    delta_c: sp.Expr,
    frequency: sp.Symbol,
    scale: sp.Expr,
    original_symbol: sp.Expr,
) -> sp.Expr:
    raw_scaled = scaled_fourier_symbol(original_symbol, frequency, scale)
    indicator = chi_plus if k < j else chi_minus
    raw_indicator = indicator(frequency / scale)
    if not raw_scaled.has(raw_indicator):
        raise RelativeWienerHopfError("scaled symbol lost its half-line indicator.")
    halfline_indicator_is_scale_invariant(
        "plus" if k < j else "minus",
        frequency,
        scale,
    )
    coefficient = sp.Integer(2) if k < j else sp.Integer(-2)
    return coefficient * indicator(frequency) * sp.exp(delta_c * frequency / scale)


def _validate_indices(k: int, j: int) -> None:
    if not isinstance(k, int) or not isinstance(j, int) or isinstance(k, bool) or isinstance(j, bool):
        raise TypeError("k and j must be integer indices.")
    if k == j:
        raise RelativeWienerHopfError("the off-diagonal family requires k != j.")


def _validate_positive_gamma(name: str, gamma: sp.Expr) -> None:
    if not isinstance(gamma, sp.Expr):
        raise TypeError(f"{name} must be a SymPy expression.")
    if gamma.is_real is not True or gamma.is_positive is not True:
        raise RelativeWienerHopfError(f"{name} must be declared real and positive.")


def _validate_delta_c(k: int, j: int, delta_c: sp.Expr) -> None:
    expected_sign = delta_c.is_negative if k < j else delta_c.is_positive
    if delta_c.is_real is not True or expected_sign is not True:
        relation = "negative" if k < j else "positive"
        raise RelativeWienerHopfError(
            f"c_k - c_j must be declared real and {relation} for these indices."
        )


@dataclass(frozen=True)
class DilationMap:
    """A typed branch or relative dilation map."""

    kind: Literal["branch", "relative"]
    k: int
    j: int | None
    scale: sp.Expr


@dataclass(frozen=True)
class DilationOperatorModel:
    """Composition operator generated by a typed dilation map."""

    dilation: DilationMap
    inverse: bool = False


@dataclass(frozen=True)
class WienerHopfOperatorModel:
    """Wiener--Hopf operator carrying one scalar L1 symbol."""

    placement: Literal["original", "left", "right"]
    symbol: sp.Expr


RelativeOperatorFactor = DilationOperatorModel | WienerHopfOperatorModel


@dataclass(frozen=True)
class OrderedRelativeOperatorProduct:
    """Operator composition whose tuple preserves noncommutative order."""

    factors: tuple[RelativeOperatorFactor, ...]

    def __post_init__(self) -> None:
        if len(self.factors) not in {2, 3}:
            raise RelativeWienerHopfError("an operator product needs two or three factors.")
        if not all(
            isinstance(factor, (DilationOperatorModel, WienerHopfOperatorModel))
            for factor in self.factors
        ):
            raise TypeError("operator product factors must be typed operator models.")


@dataclass(frozen=True)
class RelativeWienerHopfIdentity:
    """The three exactly equal ordered operator products."""

    original: OrderedRelativeOperatorProduct
    left: OrderedRelativeOperatorProduct
    right: OrderedRelativeOperatorProduct
    exact: bool = True

    def __post_init__(self) -> None:
        if self.exact is not True:
            raise RelativeWienerHopfError("the relative Wiener--Hopf identity is exact.")


@dataclass(frozen=True)
class RelativeWienerHopfAction:
    """Common concrete integral action of all three products."""

    kernel: sp.Expr
    output_variable: sp.Symbol
    input_variable: sp.Symbol
    input_function: sp.FunctionClass
    integral: sp.Integral


@dataclass(frozen=True)
class RelativeWienerHopfDerivationTrace:
    """E2E trace built solely from one verified L1 factorization."""

    factorization: DilationConjugatedWienerHopfFactorization
    branch_k: DilationMap
    branch_j: DilationMap
    relative_dilation: DilationMap
    original_operator: OrderedRelativeOperatorProduct
    left_operator: OrderedRelativeOperatorProduct
    right_operator: OrderedRelativeOperatorProduct
    identity: RelativeWienerHopfIdentity
    action: RelativeWienerHopfAction
    symbol_correspondence_forward: sp.Expr
    symbol_correspondence_reverse: sp.Expr
    reconstruction_checks: tuple[bool, bool]


def build_relative_wiener_hopf_trace(
    k: int,
    j: int,
    gamma_k: sp.Expr,
    gamma_j: sp.Expr,
    c_k: sp.Expr,
    c_j: sp.Expr,
) -> RelativeWienerHopfDerivationTrace:
    """Build the typed exact operator trace from the verified L1 result."""

    l1 = factor_dilation_conjugated_wiener_hopf(
        k, j, gamma_k, gamma_j, c_k, c_j
    )
    branch_k = DilationMap("branch", k, None, l1.gamma_k)
    branch_j = DilationMap("branch", j, None, l1.gamma_j)
    beta = DilationMap("relative", k, j, l1.relative_scale)
    v_k = DilationOperatorModel(branch_k)
    v_j_inverse = DilationOperatorModel(branch_j, inverse=True)
    v_beta = DilationOperatorModel(beta)
    w_original = WienerHopfOperatorModel("original", l1.original_symbol)
    w_left = WienerHopfOperatorModel("left", l1.left_symbol)
    w_right = WienerHopfOperatorModel("right", l1.right_symbol)
    original = OrderedRelativeOperatorProduct((v_k, w_original, v_j_inverse))
    left = OrderedRelativeOperatorProduct((v_beta, w_left))
    right = OrderedRelativeOperatorProduct((w_right, v_beta))
    identity = RelativeWienerHopfIdentity(original, left, right)

    f = sp.Function("f")
    x, y = l1.output_variable, l1.input_variable
    action = RelativeWienerHopfAction(
        kernel=l1.integral_kernel,
        output_variable=x,
        input_variable=y,
        input_function=f,
        integral=sp.Integral(l1.integral_kernel * f(y), (y, 0, sp.oo)),
    )
    forward = l1.right_symbol
    reverse = l1.left_symbol
    checks = (
        sp.simplify(l1.left_reconstructed_kernel - l1.integral_kernel) == 0,
        sp.simplify(l1.right_reconstructed_kernel - l1.integral_kernel) == 0,
    )
    if beta.scale != l1.gamma_k / l1.gamma_j:
        raise RelativeWienerHopfError("the relative dilation does not match L1.")
    if not all(checks):
        raise RelativeWienerHopfError("an operator factorization lost its L1 kernel.")
    return RelativeWienerHopfDerivationTrace(
        factorization=l1,
        branch_k=branch_k,
        branch_j=branch_j,
        relative_dilation=beta,
        original_operator=original,
        left_operator=left,
        right_operator=right,
        identity=identity,
        action=action,
        symbol_correspondence_forward=forward,
        symbol_correspondence_reverse=reverse,
        reconstruction_checks=checks,
    )
