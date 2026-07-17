"""Exact dilation factorizations for the article's off-diagonal WH kernel."""

from __future__ import annotations

from dataclasses import dataclass, field
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
