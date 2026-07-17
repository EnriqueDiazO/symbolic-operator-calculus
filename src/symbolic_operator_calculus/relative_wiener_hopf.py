"""Exact dilation factorizations for the article's off-diagonal WH kernel."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Literal

import sympy as sp

from .fourier import scaled_convolution_kernel, scaled_fourier_symbol
from .substitution import (
    collect_bound_symbols,
    fresh_symbol,
    substitute_free_variable,
)


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

    def __post_init__(self) -> None:
        _validate_indices(self.k, self.j)
        _validate_positive_gamma("gamma_k", self.gamma_k)
        _validate_positive_gamma("gamma_j", self.gamma_j)
        _validate_delta_c(self.k, self.j, self.delta_c)
        variables = {
            "time_variable": self.time_variable,
            "frequency_variable": self.frequency_variable,
            "output_variable": self.output_variable,
            "input_variable": self.input_variable,
        }
        for name, variable in variables.items():
            if not isinstance(variable, sp.Symbol):
                raise TypeError(f"{name} must be a SymPy Symbol.")

        _validate_positive_gamma("relative_scale", self.relative_scale)
        _require_equivalent(
            "relative_scale",
            self.relative_scale,
            self.gamma_k / self.gamma_j,
        )
        expected_original_kernel = article_wiener_hopf_kernel(
            self.time_variable,
            self.delta_c,
        )
        expected_original_symbol = article_wiener_hopf_symbol(
            self.k,
            self.j,
            self.delta_c,
            self.frequency_variable,
        )
        expected_integral_kernel = self.gamma_j * expected_original_kernel.xreplace(
            {
                self.time_variable: (
                    self.gamma_k * self.output_variable
                    - self.gamma_j * self.input_variable
                )
            }
        )
        expected_left_kernel = scaled_convolution_kernel(
            expected_original_kernel,
            self.time_variable,
            self.gamma_j,
        )
        expected_left_symbol = _scaled_article_symbol(
            self.k,
            self.j,
            self.delta_c,
            self.frequency_variable,
            self.gamma_j,
            expected_original_symbol,
        )
        expected_left_reconstructed = expected_left_kernel.xreplace(
            {
                self.time_variable: (
                    self.relative_scale * self.output_variable
                    - self.input_variable
                )
            }
        )
        expected_right_kernel = scaled_convolution_kernel(
            expected_original_kernel,
            self.time_variable,
            self.gamma_k,
        )
        expected_right_symbol = _scaled_article_symbol(
            self.k,
            self.j,
            self.delta_c,
            self.frequency_variable,
            self.gamma_k,
            expected_original_symbol,
        )
        expected_right_reconstructed = (
            expected_right_kernel.xreplace(
                {
                    self.time_variable: (
                        self.output_variable
                        - self.input_variable / self.relative_scale
                    )
                }
            )
            / self.relative_scale
        )
        expected_fields = {
            "original_kernel": expected_original_kernel,
            "original_symbol": expected_original_symbol,
            "integral_kernel": expected_integral_kernel,
            "left_convolution_kernel": expected_left_kernel,
            "left_symbol": expected_left_symbol,
            "left_reconstructed_kernel": expected_left_reconstructed,
            "right_convolution_kernel": expected_right_kernel,
            "right_symbol": expected_right_symbol,
            "right_reconstructed_kernel": expected_right_reconstructed,
        }
        for name, expected in expected_fields.items():
            _require_equivalent(name, getattr(self, name), expected)

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
    _validate_positive_index("k", k)
    _validate_positive_index("j", j)
    if k == j:
        raise RelativeWienerHopfError("the off-diagonal family requires k != j.")


def _validate_positive_index(name: str, index: object) -> None:
    if type(index) is not int:
        raise TypeError(f"{name} must be an integer index.")
    if index <= 0:
        raise RelativeWienerHopfError(f"{name} must be strictly positive.")


def _validate_positive_gamma(name: str, gamma: sp.Expr) -> None:
    if not isinstance(gamma, sp.Expr):
        raise TypeError(f"{name} must be a SymPy expression.")
    if (
        gamma.is_real is not True
        or gamma.is_positive is not True
        or gamma.is_finite is not True
    ):
        raise RelativeWienerHopfError(
            f"{name} must be declared real, positive, and finite."
        )


def _validate_choice(name: str, value: object, choices: set[str]) -> None:
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a string.")
    if value not in choices:
        allowed = ", ".join(sorted(repr(choice) for choice in choices))
        raise RelativeWienerHopfError(f"{name} must be one of {allowed}.")


def _require_equivalent(name: str, actual: object, expected: sp.Expr) -> None:
    if not isinstance(actual, sp.Expr):
        raise TypeError(f"{name} must be a SymPy expression.")
    if sp.simplify(actual - expected) != 0:
        raise RelativeWienerHopfError(
            f"{name} is inconsistent with the relative Wiener--Hopf model."
        )


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
    gamma_k: sp.Expr | None = None
    gamma_j: sp.Expr | None = None

    def __post_init__(self) -> None:
        _validate_choice("kind", self.kind, {"branch", "relative"})
        _validate_positive_index("k", self.k)
        _validate_positive_gamma("scale", self.scale)
        if self.kind == "branch":
            if self.j is not None:
                raise RelativeWienerHopfError(
                    "a branch dilation must not carry a second index."
                )
            if self.gamma_k is not None or self.gamma_j is not None:
                raise RelativeWienerHopfError(
                    "a branch dilation must not carry relative-scale metadata."
                )
            return

        _validate_positive_index("j", self.j)
        if self.k == self.j:
            raise RelativeWienerHopfError(
                "a relative dilation requires distinct indices."
            )
        _validate_positive_gamma("gamma_k", self.gamma_k)
        _validate_positive_gamma("gamma_j", self.gamma_j)
        _require_equivalent("scale", self.scale, self.gamma_k / self.gamma_j)


@dataclass(frozen=True)
class DilationOperatorModel:
    """Composition operator generated by a typed dilation map."""

    dilation: DilationMap
    inverse: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.dilation, DilationMap):
            raise TypeError("dilation must be a DilationMap.")
        if type(self.inverse) is not bool:
            raise TypeError("inverse must be a bool.")


@dataclass(frozen=True)
class WienerHopfOperatorModel:
    """Wiener--Hopf operator carrying one scalar L1 symbol."""

    placement: Literal["original", "left", "right"]
    symbol: sp.Expr
    factorization: DilationConjugatedWienerHopfFactorization

    def __post_init__(self) -> None:
        _validate_choice(
            "placement",
            self.placement,
            {"original", "left", "right"},
        )
        if not isinstance(self.symbol, sp.Expr):
            raise TypeError("symbol must be a SymPy expression.")
        if not isinstance(
            self.factorization,
            DilationConjugatedWienerHopfFactorization,
        ):
            raise TypeError(
                "factorization must be a DilationConjugatedWienerHopfFactorization."
            )
        expected_symbol = {
            "original": self.factorization.original_symbol,
            "left": self.factorization.left_symbol,
            "right": self.factorization.right_symbol,
        }[self.placement]
        _require_equivalent("symbol", self.symbol, expected_symbol)

    @property
    def k(self) -> int:
        return self.factorization.k

    @property
    def j(self) -> int:
        return self.factorization.j


RelativeOperatorFactor = DilationOperatorModel | WienerHopfOperatorModel


@dataclass(frozen=True)
class OrderedRelativeOperatorProduct:
    """Operator composition whose tuple preserves noncommutative order."""

    factors: tuple[RelativeOperatorFactor, ...]

    def __post_init__(self) -> None:
        if self.factors is None:
            raise TypeError("operator product factors must be an iterable.")
        try:
            factors = tuple(self.factors)
        except TypeError as exc:
            raise TypeError("operator product factors must be an iterable.") from exc
        object.__setattr__(self, "factors", factors)
        if len(factors) not in {2, 3}:
            raise RelativeWienerHopfError("an operator product needs two or three factors.")
        if not all(
            isinstance(factor, (DilationOperatorModel, WienerHopfOperatorModel))
            for factor in factors
        ):
            raise TypeError("operator product factors must be typed operator models.")


@dataclass(frozen=True)
class RelativeWienerHopfIdentity:
    """The three exactly equal ordered operator products."""

    original: OrderedRelativeOperatorProduct
    left: OrderedRelativeOperatorProduct
    right: OrderedRelativeOperatorProduct
    exact: bool = field(init=False, default=True)

    def __post_init__(self) -> None:
        original_model = _validate_relative_product_shape(self.original, "original")
        left_model = _validate_relative_product_shape(self.left, "left")
        right_model = _validate_relative_product_shape(self.right, "right")
        if original_model != left_model or original_model != right_model:
            raise RelativeWienerHopfError(
                "all identity products must come from the same L1 factorization."
            )


def _validate_relative_product_shape(
    product: OrderedRelativeOperatorProduct,
    expected: Literal["original", "left", "right"],
) -> DilationConjugatedWienerHopfFactorization:
    if not isinstance(product, OrderedRelativeOperatorProduct):
        raise RelativeWienerHopfError(
            f"the {expected} identity member must be an ordered product."
        )
    factors = product.factors
    if expected == "original":
        if not (
            len(factors) == 3
            and isinstance(factors[0], DilationOperatorModel)
            and isinstance(factors[1], WienerHopfOperatorModel)
            and isinstance(factors[2], DilationOperatorModel)
        ):
            raise RelativeWienerHopfError(
                "the original product must have branch, Wiener--Hopf, branch order."
            )
        first, wh_factor, last = factors
        model = wh_factor.factorization
        _validate_branch_factor(first, model.k, model.gamma_k, inverse=False)
        _validate_branch_factor(last, model.j, model.gamma_j, inverse=True)
    elif expected == "left":
        if not (
            len(factors) == 2
            and isinstance(factors[0], DilationOperatorModel)
            and isinstance(factors[1], WienerHopfOperatorModel)
        ):
            raise RelativeWienerHopfError(
                "the left product must place the relative dilation first."
            )
        relative_factor, wh_factor = factors
        model = wh_factor.factorization
        _validate_relative_factor(relative_factor, model)
    else:
        if not (
            len(factors) == 2
            and isinstance(factors[0], WienerHopfOperatorModel)
            and isinstance(factors[1], DilationOperatorModel)
        ):
            raise RelativeWienerHopfError(
                "the right product must place the relative dilation last."
            )
        wh_factor, relative_factor = factors
        model = wh_factor.factorization
        _validate_relative_factor(relative_factor, model)
    if wh_factor.placement != expected:
        raise RelativeWienerHopfError(
            f"the {expected} product has the wrong Wiener--Hopf placement."
        )
    return model


def _validate_branch_factor(
    factor: DilationOperatorModel,
    index: int,
    scale: sp.Expr,
    *,
    inverse: bool,
) -> None:
    dilation = factor.dilation
    if (
        dilation.kind != "branch"
        or dilation.k != index
        or dilation.j is not None
        or factor.inverse is not inverse
    ):
        raise RelativeWienerHopfError(
            "a branch factor has incompatible index, role, or direction."
        )
    _require_equivalent("branch scale", dilation.scale, scale)


def _validate_relative_factor(
    factor: DilationOperatorModel,
    model: DilationConjugatedWienerHopfFactorization,
) -> None:
    dilation = factor.dilation
    if (
        dilation.kind != "relative"
        or dilation.k != model.k
        or dilation.j != model.j
        or factor.inverse
    ):
        raise RelativeWienerHopfError(
            "a relative factor has incompatible indices, role, or direction."
        )
    _require_equivalent("relative gamma_k", dilation.gamma_k, model.gamma_k)
    _require_equivalent("relative gamma_j", dilation.gamma_j, model.gamma_j)
    _require_equivalent("relative scale", dilation.scale, model.relative_scale)


@dataclass(frozen=True)
class RelativeWienerHopfAction:
    """Common concrete integral action of all three products."""

    kernel: sp.Expr
    output_variable: sp.Symbol
    input_variable: sp.Symbol
    input_function: sp.FunctionClass
    integral: sp.Integral

    def __post_init__(self) -> None:
        if not isinstance(self.kernel, sp.Expr):
            raise TypeError("kernel must be a SymPy expression.")
        if not isinstance(self.output_variable, sp.Symbol):
            raise TypeError("output_variable must be a SymPy Symbol.")
        if not isinstance(self.input_variable, sp.Symbol):
            raise TypeError("input_variable must be a SymPy Symbol.")
        if not isinstance(self.input_function, sp.FunctionClass):
            raise TypeError("input_function must be a SymPy function class.")
        if not isinstance(self.integral, sp.Integral):
            raise TypeError("integral must be a SymPy Integral.")
        expected_integral = sp.Integral(
            self.kernel * self.input_function(self.input_variable),
            (self.input_variable, 0, sp.oo),
        )
        if self.integral != expected_integral:
            raise RelativeWienerHopfError(
                "integral is inconsistent with the stored relative action."
            )


@dataclass(frozen=True)
class RelativeProductAction:
    """Direct and normalized actions computed from one canonical product."""

    product: OrderedRelativeOperatorProduct
    operand: sp.Expr
    direct_expression: sp.Integral
    normalized_expression: sp.Integral
    output_variable: sp.Symbol
    input_variable: sp.Symbol
    kernel: sp.Expr
    direct_integration_variable: sp.Symbol

    def __post_init__(self) -> None:
        _classify_relative_product(self.product)
        if not isinstance(self.operand, sp.Expr):
            raise TypeError("operand must be a SymPy expression.")
        if not isinstance(self.direct_expression, sp.Integral):
            raise TypeError("direct_expression must be a SymPy Integral.")
        if not isinstance(self.normalized_expression, sp.Integral):
            raise TypeError("normalized_expression must be a SymPy Integral.")
        if not isinstance(self.output_variable, sp.Symbol):
            raise TypeError("output_variable must be a SymPy Symbol.")
        if not isinstance(self.input_variable, sp.Symbol):
            raise TypeError("input_variable must be a SymPy Symbol.")
        if not isinstance(self.kernel, sp.Expr):
            raise TypeError("kernel must be a SymPy expression.")
        if not isinstance(self.direct_integration_variable, sp.Symbol):
            raise TypeError(
                "direct_integration_variable must be a SymPy Symbol."
            )
        expected = sp.Integral(
            self.kernel * self.operand,
            (self.input_variable, 0, sp.oo),
        )
        if self.normalized_expression != expected:
            raise RelativeWienerHopfError(
                "normalized_expression is inconsistent with its kernel and operand."
            )
        if self.direct_integration_variable not in collect_bound_symbols(
            self.direct_expression
        ):
            raise RelativeWienerHopfError(
                "the direct integration variable must be bound by the direct action."
            )


@dataclass(frozen=True)
class RelativeProductActionVerification:
    """Derived evidence that all three canonical product actions agree."""

    identity: RelativeWienerHopfIdentity
    original: RelativeProductAction
    left: RelativeProductAction
    right: RelativeProductAction
    actions_verified: bool = field(init=False, default=True)

    def __post_init__(self) -> None:
        if not isinstance(self.identity, RelativeWienerHopfIdentity):
            raise TypeError("identity must be a RelativeWienerHopfIdentity.")
        actions = (self.original, self.left, self.right)
        if not all(isinstance(action, RelativeProductAction) for action in actions):
            raise TypeError("all verification members must be RelativeProductAction.")
        if tuple(action.product for action in actions) != (
            self.identity.original,
            self.identity.left,
            self.identity.right,
        ):
            raise RelativeWienerHopfError(
                "verification actions must correspond to the stored identity."
            )
        output_variables = {action.output_variable for action in actions}
        input_variables = {action.input_variable for action in actions}
        if len(output_variables) != 1 or len(input_variables) != 1:
            raise RelativeWienerHopfError(
                "verification actions must use common output and input variables."
            )
        if any(
            sp.simplify(self.original.kernel - action.kernel) != 0
            for action in (self.left, self.right)
        ):
            raise RelativeWienerHopfError(
                "the canonical relative product kernels are not equal."
            )
        model = _validate_relative_product_shape(self.identity.original, "original")
        expected = model.integral_kernel.xreplace(
            {
                model.output_variable: self.original.output_variable,
                model.input_variable: self.original.input_variable,
            }
        )
        if any(sp.simplify(action.kernel - expected) != 0 for action in actions):
            raise RelativeWienerHopfError(
                "a canonical product kernel does not match the conjugated L1 kernel."
            )


@dataclass(frozen=True)
class RelativeSymbolCorrespondenceVerification:
    """Independent Fourier-scaling evidence for both relative placements."""

    factorization: DilationConjugatedWienerHopfFactorization
    pair: tuple[int, int]
    frequency_variable: sp.Symbol
    left_scale: sp.Expr
    right_scale: sp.Expr
    computed_left_symbol: sp.Expr
    stored_left_symbol: sp.Expr
    computed_right_symbol: sp.Expr
    stored_right_symbol: sp.Expr
    computed_left_kernel: sp.Expr
    stored_left_kernel: sp.Expr
    computed_right_kernel: sp.Expr
    stored_right_kernel: sp.Expr
    left_symbol_verified: bool = field(init=False, default=True)
    right_symbol_verified: bool = field(init=False, default=True)
    left_kernel_verified: bool = field(init=False, default=True)
    right_kernel_verified: bool = field(init=False, default=True)
    correspondences_verified: bool = field(init=False, default=True)

    def __post_init__(self) -> None:
        model = _validated_symbol_factorization(self.factorization)
        if (
            type(self.pair) is not tuple
            or len(self.pair) != 2
            or any(type(index) is not int for index in self.pair)
        ):
            raise TypeError("pair must be a two-item tuple of integer indices.")
        if self.pair != (model.k, model.j):
            raise RelativeWienerHopfError(
                "pair is inconsistent with the stored L1 factorization."
            )
        _validate_verification_frequency(self.frequency_variable)
        _require_algebraically_equivalent(
            "left_scale",
            self.left_scale,
            model.gamma_j,
        )
        _require_algebraically_equivalent(
            "right_scale",
            self.right_scale,
            model.gamma_k,
        )

        expected_stored_left_symbol = substitute_free_variable(
            model.left_symbol,
            model.frequency_variable,
            self.frequency_variable,
        )
        expected_stored_right_symbol = substitute_free_variable(
            model.right_symbol,
            model.frequency_variable,
            self.frequency_variable,
        )
        expected_computed_left_symbol = _reconstruct_scaled_relative_symbol(
            model,
            self.frequency_variable,
            self.left_scale,
        )
        expected_computed_right_symbol = _reconstruct_scaled_relative_symbol(
            model,
            self.frequency_variable,
            self.right_scale,
        )
        expected_computed_left_kernel = scaled_convolution_kernel(
            model.original_kernel,
            model.time_variable,
            self.left_scale,
        )
        expected_computed_right_kernel = scaled_convolution_kernel(
            model.original_kernel,
            model.time_variable,
            self.right_scale,
        )
        comparisons = (
            (
                "stored_left_symbol",
                self.stored_left_symbol,
                expected_stored_left_symbol,
            ),
            (
                "stored_right_symbol",
                self.stored_right_symbol,
                expected_stored_right_symbol,
            ),
            (
                "computed_left_symbol",
                self.computed_left_symbol,
                expected_computed_left_symbol,
            ),
            (
                "computed_right_symbol",
                self.computed_right_symbol,
                expected_computed_right_symbol,
            ),
            (
                "stored_left_kernel",
                self.stored_left_kernel,
                model.left_convolution_kernel,
            ),
            (
                "stored_right_kernel",
                self.stored_right_kernel,
                model.right_convolution_kernel,
            ),
            (
                "computed_left_kernel",
                self.computed_left_kernel,
                expected_computed_left_kernel,
            ),
            (
                "computed_right_kernel",
                self.computed_right_kernel,
                expected_computed_right_kernel,
            ),
            (
                "left symbol correspondence",
                self.computed_left_symbol,
                self.stored_left_symbol,
            ),
            (
                "right symbol correspondence",
                self.computed_right_symbol,
                self.stored_right_symbol,
            ),
            (
                "left kernel correspondence",
                self.computed_left_kernel,
                self.stored_left_kernel,
            ),
            (
                "right kernel correspondence",
                self.computed_right_kernel,
                self.stored_right_kernel,
            ),
        )
        for name, actual, expected in comparisons:
            _require_algebraically_equivalent(name, actual, expected)


def verify_relative_symbol_correspondences(
    factorization_or_identity: (
        DilationConjugatedWienerHopfFactorization
        | RelativeWienerHopfIdentity
    ),
    *,
    frequency_variable: sp.Symbol | None = None,
) -> RelativeSymbolCorrespondenceVerification:
    """Reconstruct and verify both L1 symbol/kernel scaling correspondences."""

    if isinstance(factorization_or_identity, RelativeWienerHopfIdentity):
        model = _validate_relative_product_shape(
            factorization_or_identity.original,
            "original",
        )
    elif isinstance(
        factorization_or_identity,
        DilationConjugatedWienerHopfFactorization,
    ):
        model = factorization_or_identity
    else:
        raise TypeError(
            "factorization_or_identity must be a relative factorization or identity."
        )
    model = _validated_symbol_factorization(model)
    frequency = (
        model.frequency_variable
        if frequency_variable is None
        else frequency_variable
    )
    _validate_verification_frequency(frequency)

    computed_left_symbol = _reconstruct_scaled_relative_symbol(
        model,
        frequency,
        model.gamma_j,
    )
    computed_right_symbol = _reconstruct_scaled_relative_symbol(
        model,
        frequency,
        model.gamma_k,
    )
    computed_left_kernel = scaled_convolution_kernel(
        model.original_kernel,
        model.time_variable,
        model.gamma_j,
    )
    computed_right_kernel = scaled_convolution_kernel(
        model.original_kernel,
        model.time_variable,
        model.gamma_k,
    )
    return RelativeSymbolCorrespondenceVerification(
        factorization=model,
        pair=(model.k, model.j),
        frequency_variable=frequency,
        left_scale=model.gamma_j,
        right_scale=model.gamma_k,
        computed_left_symbol=computed_left_symbol,
        stored_left_symbol=substitute_free_variable(
            model.left_symbol,
            model.frequency_variable,
            frequency,
        ),
        computed_right_symbol=computed_right_symbol,
        stored_right_symbol=substitute_free_variable(
            model.right_symbol,
            model.frequency_variable,
            frequency,
        ),
        computed_left_kernel=computed_left_kernel,
        stored_left_kernel=model.left_convolution_kernel,
        computed_right_kernel=computed_right_kernel,
        stored_right_kernel=model.right_convolution_kernel,
    )


def _validated_symbol_factorization(
    factorization: object,
) -> DilationConjugatedWienerHopfFactorization:
    if not isinstance(
        factorization,
        DilationConjugatedWienerHopfFactorization,
    ):
        raise TypeError(
            "factorization must be a DilationConjugatedWienerHopfFactorization."
        )
    required_fields = (
        "k",
        "j",
        "gamma_k",
        "gamma_j",
        "time_variable",
        "frequency_variable",
        "original_kernel",
        "original_symbol",
        "left_convolution_kernel",
        "left_symbol",
        "right_convolution_kernel",
        "right_symbol",
    )
    for name in required_fields:
        if not hasattr(factorization, name):
            raise RelativeWienerHopfError(
                f"factorization is missing required field {name!r}."
            )
    return factorization


def _validate_verification_frequency(frequency: object) -> None:
    if not isinstance(frequency, sp.Symbol):
        raise TypeError("frequency_variable must be a SymPy Symbol.")
    if frequency.is_real is not True:
        raise RelativeWienerHopfError(
            "frequency_variable must be a real SymPy Symbol."
        )


def _reconstruct_scaled_relative_symbol(
    model: DilationConjugatedWienerHopfFactorization,
    frequency: sp.Symbol,
    scale: sp.Expr,
) -> sp.Expr:
    original_symbol = substitute_free_variable(
        model.original_symbol,
        model.frequency_variable,
        frequency,
    )
    scaled_symbol = scaled_fourier_symbol(
        original_symbol,
        frequency,
        scale,
    )
    indicator = chi_plus if model.k < model.j else chi_minus
    raw_indicator = indicator(frequency / scale)
    if not scaled_symbol.has(raw_indicator):
        raise RelativeWienerHopfError(
            "the reconstructed symbol lost its half-line indicator."
        )
    halfline_indicator_is_scale_invariant(
        "plus" if model.k < model.j else "minus",
        frequency,
        scale,
    )
    return scaled_symbol.xreplace({raw_indicator: indicator(frequency)})


def _require_algebraically_equivalent(
    name: str,
    actual: object,
    expected: object,
) -> None:
    if not isinstance(actual, sp.Expr) or not isinstance(expected, sp.Expr):
        raise TypeError(f"{name} must compare two SymPy expressions.")
    if actual == expected:
        return
    try:
        equivalent = sp.simplify(actual - expected) == 0
    except (TypeError, ValueError):
        equivalent = False
    if not equivalent:
        raise RelativeWienerHopfError(
            f"{name} is inconsistent with the independently scaled L1 model."
        )


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

    def __post_init__(self) -> None:
        if not isinstance(
            self.factorization,
            DilationConjugatedWienerHopfFactorization,
        ):
            raise TypeError(
                "factorization must be a DilationConjugatedWienerHopfFactorization."
            )
        if not all(
            isinstance(item, DilationMap)
            for item in (self.branch_k, self.branch_j, self.relative_dilation)
        ):
            raise TypeError("trace dilations must be DilationMap instances.")
        model = self.factorization
        _validate_branch_map(self.branch_k, model.k, model.gamma_k)
        _validate_branch_map(self.branch_j, model.j, model.gamma_j)
        _validate_relative_map(self.relative_dilation, model)
        if not isinstance(self.identity, RelativeWienerHopfIdentity):
            raise TypeError("identity must be a RelativeWienerHopfIdentity.")
        if (
            self.original_operator != self.identity.original
            or self.left_operator != self.identity.left
            or self.right_operator != self.identity.right
        ):
            raise RelativeWienerHopfError(
                "trace products must be the products stored by its identity."
            )
        identity_model = _validate_relative_product_shape(
            self.original_operator,
            "original",
        )
        if identity_model != model:
            raise RelativeWienerHopfError(
                "trace products must come from its stored L1 factorization."
            )
        if not isinstance(self.action, RelativeWienerHopfAction):
            raise TypeError("action must be a RelativeWienerHopfAction.")
        _require_equivalent("action kernel", self.action.kernel, model.integral_kernel)
        if (
            self.action.output_variable != model.output_variable
            or self.action.input_variable != model.input_variable
        ):
            raise RelativeWienerHopfError(
                "trace action variables must match its L1 factorization."
            )
        _require_equivalent(
            "symbol_correspondence_forward",
            self.symbol_correspondence_forward,
            model.right_symbol,
        )
        _require_equivalent(
            "symbol_correspondence_reverse",
            self.symbol_correspondence_reverse,
            model.left_symbol,
        )
        if (
            type(self.reconstruction_checks) is not tuple
            or len(self.reconstruction_checks) != 2
            or any(type(check) is not bool for check in self.reconstruction_checks)
        ):
            raise TypeError("reconstruction_checks must be a pair of bool values.")
        if self.reconstruction_checks != (True, True):
            raise RelativeWienerHopfError(
                "both stored L1 kernel reconstructions must hold."
            )


def _validate_branch_map(
    dilation: DilationMap,
    index: int,
    scale: sp.Expr,
) -> None:
    if dilation.kind != "branch" or dilation.k != index or dilation.j is not None:
        raise RelativeWienerHopfError(
            "a trace branch dilation has incompatible index or role."
        )
    _require_equivalent("branch scale", dilation.scale, scale)


def _validate_relative_map(
    dilation: DilationMap,
    model: DilationConjugatedWienerHopfFactorization,
) -> None:
    if (
        dilation.kind != "relative"
        or dilation.k != model.k
        or dilation.j != model.j
    ):
        raise RelativeWienerHopfError(
            "the trace relative dilation has incompatible indices or role."
        )
    _require_equivalent("relative gamma_k", dilation.gamma_k, model.gamma_k)
    _require_equivalent("relative gamma_j", dilation.gamma_j, model.gamma_j)
    _require_equivalent("relative scale", dilation.scale, model.relative_scale)


def apply_relative_factor(
    factor: RelativeOperatorFactor,
    operand: sp.Expr,
    *,
    output_variable: sp.Symbol,
    integration_variables: Iterable[sp.Symbol] = (),
) -> sp.Expr:
    """Apply one typed relative-domain factor to a SymPy operand."""

    if not isinstance(operand, sp.Expr):
        raise TypeError("operand must be a SymPy expression.")
    if not isinstance(output_variable, sp.Symbol):
        raise TypeError("output_variable must be a SymPy Symbol.")
    variables = _validated_integration_variables(integration_variables)
    if isinstance(factor, DilationOperatorModel):
        if variables:
            raise RelativeWienerHopfError(
                "a dilation action does not use integration variables."
            )
        scale = factor.dilation.scale
        replacement = (
            output_variable / scale
            if factor.inverse
            else scale * output_variable
        )
        return substitute_free_variable(operand, output_variable, replacement)
    if not isinstance(factor, WienerHopfOperatorModel):
        raise TypeError("factor must be a typed relative operator factor.")
    if len(variables) > 1:
        raise RelativeWienerHopfError(
            "one Wiener--Hopf factor accepts at most one integration variable."
        )
    kernel = _wiener_hopf_kernel(factor)
    integration_variable = _select_integration_variable(
        operand,
        kernel,
        output_variable,
        variables[0] if variables else None,
    )
    kernel_at_difference = substitute_free_variable(
        kernel,
        factor.factorization.time_variable,
        output_variable - integration_variable,
    )
    integrated_operand = substitute_free_variable(
        operand,
        output_variable,
        integration_variable,
    )
    return sp.Integral(
        kernel_at_difference * integrated_operand,
        (integration_variable, 0, sp.oo),
    )


def apply_ordered_relative_product(
    product: OrderedRelativeOperatorProduct,
    operand: sp.Expr,
    *,
    output_variable: sp.Symbol,
    input_variable: sp.Symbol | None = None,
    integration_variables: Iterable[sp.Symbol] = (),
) -> RelativeProductAction:
    """Apply one canonical product right-to-left and normalize its action."""

    if not isinstance(product, OrderedRelativeOperatorProduct):
        raise TypeError("product must be an OrderedRelativeOperatorProduct.")
    if not isinstance(operand, sp.Expr):
        raise TypeError("operand must be a SymPy expression.")
    if not isinstance(output_variable, sp.Symbol):
        raise TypeError("output_variable must be a SymPy Symbol.")
    placement, model = _classify_relative_product(product)
    source_variable = model.input_variable if input_variable is None else input_variable
    if not isinstance(source_variable, sp.Symbol):
        raise TypeError("input_variable must be a SymPy Symbol.")
    variables = _validated_integration_variables(integration_variables)
    integral_factor_count = sum(
        isinstance(factor, WienerHopfOperatorModel) for factor in product.factors
    )
    if variables and len(variables) != integral_factor_count:
        raise RelativeWienerHopfError(
            "integration_variables must supply one distinct symbol per integral factor."
        )

    evaluation_variable = _evaluation_variable(
        operand,
        source_variable,
        output_variable,
    )
    expression = substitute_free_variable(
        operand,
        source_variable,
        evaluation_variable,
    )
    integral_index = integral_factor_count - 1
    for factor in reversed(product.factors):
        if isinstance(factor, WienerHopfOperatorModel):
            selected = (variables[integral_index],) if variables else ()
            integral_index -= 1
            expression = apply_relative_factor(
                factor,
                expression,
                output_variable=evaluation_variable,
                integration_variables=selected,
            )
        else:
            expression = apply_relative_factor(
                factor,
                expression,
                output_variable=evaluation_variable,
            )
    direct_expression = substitute_free_variable(
        expression,
        evaluation_variable,
        output_variable,
    )
    if not isinstance(direct_expression, sp.Integral):
        raise RelativeWienerHopfError(
            "a canonical relative product must produce an integral action."
        )
    direct_integration_variable = direct_expression.variables[0]
    normalized_operand, normalized_input = _normalized_operand(
        operand,
        source_variable,
        output_variable,
    )
    kernel = _normalized_relative_kernel(
        product,
        placement,
        output_variable,
        normalized_input,
    )
    normalized_expression = sp.Integral(
        kernel * normalized_operand,
        (normalized_input, 0, sp.oo),
    )
    return RelativeProductAction(
        product=product,
        operand=normalized_operand,
        direct_expression=direct_expression,
        normalized_expression=normalized_expression,
        output_variable=output_variable,
        input_variable=normalized_input,
        kernel=kernel,
        direct_integration_variable=direct_integration_variable,
    )


def verify_relative_product_actions(
    identity: RelativeWienerHopfIdentity,
    operand: sp.Expr | None = None,
    *,
    output_variable: sp.Symbol | None = None,
    input_variable: sp.Symbol | None = None,
) -> RelativeProductActionVerification:
    """Derive and compare the three canonical actions of one identity."""

    if not isinstance(identity, RelativeWienerHopfIdentity):
        raise TypeError("identity must be a RelativeWienerHopfIdentity.")
    model = _validate_relative_product_shape(identity.original, "original")
    output = model.output_variable if output_variable is None else output_variable
    input_ = model.input_variable if input_variable is None else input_variable
    if not isinstance(output, sp.Symbol):
        raise TypeError("output_variable must be a SymPy Symbol.")
    if not isinstance(input_, sp.Symbol):
        raise TypeError("input_variable must be a SymPy Symbol.")
    applied_operand = sp.Function("f")(input_) if operand is None else operand
    if not isinstance(applied_operand, sp.Expr):
        raise TypeError("operand must be a SymPy expression.")
    actions = tuple(
        apply_ordered_relative_product(
            product,
            applied_operand,
            output_variable=output,
            input_variable=input_,
        )
        for product in (identity.original, identity.left, identity.right)
    )
    return RelativeProductActionVerification(
        identity=identity,
        original=actions[0],
        left=actions[1],
        right=actions[2],
    )


def _validated_integration_variables(
    integration_variables: Iterable[sp.Symbol],
) -> tuple[sp.Symbol, ...]:
    if integration_variables is None:
        raise TypeError("integration_variables must be an iterable of Symbols.")
    try:
        variables = tuple(integration_variables)
    except TypeError as exc:
        raise TypeError(
            "integration_variables must be an iterable of Symbols."
        ) from exc
    if not all(isinstance(variable, sp.Symbol) for variable in variables):
        raise TypeError("integration_variables must contain only SymPy Symbols.")
    names = tuple(variable.name for variable in variables)
    if len(set(names)) != len(names):
        raise RelativeWienerHopfError(
            "integration_variables must not repeat symbol names."
        )
    return variables


def _select_integration_variable(
    operand: sp.Expr,
    kernel: sp.Expr,
    output_variable: sp.Symbol,
    proposed: sp.Symbol | None,
) -> sp.Symbol:
    avoid = (
        operand.atoms(sp.Symbol)
        | kernel.atoms(sp.Symbol)
        | collect_bound_symbols(operand)
        | {output_variable}
    )
    if proposed is not None and proposed.name not in {
        symbol.name for symbol in avoid
    }:
        return proposed
    if proposed is not None:
        avoid.add(proposed)
    return fresh_symbol(avoid, preferred_names=("y", "u", "v", "w", "z"))


def _evaluation_variable(
    operand: sp.Expr,
    input_variable: sp.Symbol,
    output_variable: sp.Symbol,
) -> sp.Symbol:
    operand_symbols = operand.atoms(sp.Symbol) | collect_bound_symbols(operand)
    if output_variable != input_variable and output_variable not in operand_symbols:
        return output_variable
    avoid = operand_symbols | {input_variable, output_variable}
    return fresh_symbol(avoid, preferred_names=("q", "r", "s", "w", "v"))


def _normalized_operand(
    operand: sp.Expr,
    input_variable: sp.Symbol,
    output_variable: sp.Symbol,
) -> tuple[sp.Expr, sp.Symbol]:
    bound = collect_bound_symbols(operand)
    if input_variable != output_variable and input_variable not in bound:
        return operand, input_variable
    avoid = operand.atoms(sp.Symbol) | bound | {input_variable, output_variable}
    normalized_input = fresh_symbol(
        avoid,
        preferred_names=("y", "u", "v", "w", "z"),
    )
    return (
        substitute_free_variable(operand, input_variable, normalized_input),
        normalized_input,
    )


def _classify_relative_product(
    product: OrderedRelativeOperatorProduct,
) -> tuple[
    Literal["original", "left", "right"],
    DilationConjugatedWienerHopfFactorization,
]:
    if not isinstance(product, OrderedRelativeOperatorProduct):
        raise TypeError("product must be an OrderedRelativeOperatorProduct.")
    if len(product.factors) == 3:
        placement: Literal["original", "left", "right"] = "original"
    elif isinstance(product.factors[0], DilationOperatorModel):
        placement = "left"
    else:
        placement = "right"
    return placement, _validate_relative_product_shape(product, placement)


def _wiener_hopf_kernel(factor: WienerHopfOperatorModel) -> sp.Expr:
    model = factor.factorization
    kernel = {
        "original": model.original_kernel,
        "left": model.left_convolution_kernel,
        "right": model.right_convolution_kernel,
    }.get(factor.placement)
    if not isinstance(kernel, sp.Expr):
        raise RelativeWienerHopfError(
            "the Wiener--Hopf factor has no supported convolution kernel."
        )
    return kernel


def _normalized_relative_kernel(
    product: OrderedRelativeOperatorProduct,
    placement: Literal["original", "left", "right"],
    output_variable: sp.Symbol,
    input_variable: sp.Symbol,
) -> sp.Expr:
    if placement == "original":
        left_dilation, wh_factor, right_dilation = product.factors
        kernel = _wiener_hopf_kernel(wh_factor)
        gamma_k = left_dilation.dilation.scale
        gamma_j = right_dilation.dilation.scale
        argument = gamma_k * output_variable - gamma_j * input_variable
        jacobian = gamma_j
    elif placement == "left":
        relative_dilation, wh_factor = product.factors
        kernel = _wiener_hopf_kernel(wh_factor)
        argument = (
            relative_dilation.dilation.scale * output_variable - input_variable
        )
        jacobian = sp.Integer(1)
    else:
        wh_factor, relative_dilation = product.factors
        kernel = _wiener_hopf_kernel(wh_factor)
        relative_scale = relative_dilation.dilation.scale
        argument = output_variable - input_variable / relative_scale
        jacobian = 1 / relative_scale
    return sp.simplify(
        jacobian
        * substitute_free_variable(
            kernel,
            wh_factor.factorization.time_variable,
            argument,
        )
    )


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
    beta = DilationMap(
        "relative",
        k,
        j,
        l1.relative_scale,
        gamma_k=l1.gamma_k,
        gamma_j=l1.gamma_j,
    )
    v_k = DilationOperatorModel(branch_k)
    v_j_inverse = DilationOperatorModel(branch_j, inverse=True)
    v_beta = DilationOperatorModel(beta)
    w_original = WienerHopfOperatorModel("original", l1.original_symbol, l1)
    w_left = WienerHopfOperatorModel("left", l1.left_symbol, l1)
    w_right = WienerHopfOperatorModel("right", l1.right_symbol, l1)
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
