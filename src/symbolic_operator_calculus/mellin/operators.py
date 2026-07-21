"""Typed operators and exact dilation covariance for Mellin calculations.

The noncommutative :class:`~symbolic_operator_calculus.operators.Product`
remains the authoritative operator AST.  The models in this module associate
analytic metadata with its atoms; they never make SymPy expressions act as
operator products and never infer Fredholm-algebra membership.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import sympy as sp

from ..domains import AssumptionContext, ConsistencyStatus
from ..operators import OperatorAtom, Product
from ..semantics import (
    DerivationRelationKind,
    ExactIdentity,
    ExactIdentityScope,
    FormalIdentity,
    ModCompactEquivalence,
    OperatorRule,
    RegularizerOperator,
    RuleCertificationStatus,
)
from ..singularities import detect_singularities
from .domains import MellinDomainRoleError, MellinSymbolDomain
from .expressions import MellinExpression
from .symbols import MellinSymbol, MellinSymbolDependency
from .variables import MellinVariableRole, mellin_parameter


FOURIER_CONVENTION = (
    "Ff(lambda)=int_R f(t) exp(-i t lambda) dt; "
    "F^-1g(t)=(2 pi)^-1 int_R g(lambda) exp(i t lambda) dlambda"
)
MELLIN_QUANTIZATION = (
    "Op_M(a)f(x)=(2 pi)^-1 int_R int_R+ "
    "a(x,lambda)(x/y)^(i lambda)f(y)dy/y dlambda"
)


class MellinOperatorError(ValueError):
    """Raised when typed operator metadata is insufficient or inconsistent."""


class DilationConvention(str, Enum):
    """Supported pullback convention for a positive constant dilation."""

    ARGUMENT_MULTIPLICATION = "f(x) -> gamma**nu * f(gamma*x)"


class WeightNormConvention(str, Enum):
    """The article's multiplicative-weight norm convention."""

    MULTIPLICATIVE_POWER_WEIGHT = "integral |x**delta f(x)|**p dx"


class RegularizerRepresentationStatus(str, Enum):
    """Strength of an explicitly supplied Mellin representation."""

    ASSUMED = "CONDITIONAL_ON_REGULARIZER_REPRESENTATION"
    CERTIFIED_EXACT = "CERTIFIED_EXACT"
    CERTIFIED_MOD_COMPACT = "CERTIFIED_MOD_COMPACT"


class AlgebraMembershipStatus(str, Enum):
    """Proof status of a declared operator-algebra membership claim."""

    PROVED = "PROVED"
    UNPROVED = "UNPROVED_MIXED_ALGEBRA_MEMBERSHIP"
    DISPROVED = "DISPROVED"


def _nonempty(value: object, field_name: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{field_name} must be a non-empty str.")
    if not value.strip():
        raise ValueError(f"{field_name} must be a non-empty str.")
    return value


def _items(value: object, field_name: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{field_name} must be an iterable, not text.")
    try:
        return tuple(value)  # type: ignore[arg-type]
    except TypeError as exc:
        raise TypeError(f"{field_name} must be an iterable.") from exc


def _sympy_scalar(value: object, field_name: str) -> sp.Expr:
    if isinstance(value, bool):
        raise TypeError(f"{field_name} must be a scalar SymPy expression.")
    try:
        expression = sp.sympify(value)
    except (sp.SympifyError, TypeError) as exc:
        raise TypeError(f"{field_name} must be a scalar SymPy expression.") from exc
    if not isinstance(expression, sp.Expr):
        raise TypeError(f"{field_name} must be a scalar SymPy expression.")
    return expression


def _context_has(context: AssumptionContext, proposition: sp.Basic) -> bool:
    return context.contains_exact(proposition)


def _is_explicitly_positive(
    value: sp.Expr,
    context: AssumptionContext,
) -> bool:
    if value.is_positive is True:
        return True
    try:
        return _context_has(context, sp.StrictGreaterThan(value, 0))
    except TypeError:
        return False


def _require_positive(
    value: sp.Expr,
    context: AssumptionContext,
    field_name: str,
) -> None:
    if value.is_nonpositive is True:
        raise MellinOperatorError(f"{field_name} must be strictly positive.")
    if not _is_explicitly_positive(value, context):
        raise MellinOperatorError(
            f"{field_name}>0 must be provable from the SymPy declaration or "
            "stored exactly in the AssumptionContext."
        )


def _equivalent(left: sp.Expr, right: sp.Expr) -> bool:
    return sp.simplify(left - right) == 0


def _scaled_atom(atom: OperatorAtom, suffix: str) -> OperatorAtom:
    return OperatorAtom(f"{atom.name}_{suffix}", kind=atom.kind)


@dataclass(frozen=True)
class WeightedLpSpace:
    """One scalar weighted half-line realization used by the article.

    ``weight_exponent=delta`` means
    ``||f||**p = integral |x**delta f(x)|**p dx``.  It does not mean that
    ``x**delta`` is a measure density.
    """

    label: str
    p: sp.Expr
    weight_exponent: sp.Expr
    assumptions: AssumptionContext
    source: str
    norm_convention: WeightNormConvention = (
        WeightNormConvention.MULTIPLICATIVE_POWER_WEIGHT
    )
    measure: str = "dx on (0,infinity)"

    def __post_init__(self) -> None:
        object.__setattr__(self, "label", _nonempty(self.label, "label"))
        p = _sympy_scalar(self.p, "p")
        delta = _sympy_scalar(self.weight_exponent, "weight_exponent")
        if not isinstance(self.assumptions, AssumptionContext):
            raise TypeError("assumptions must be an AssumptionContext.")
        if self.assumptions.consistency_status is ConsistencyStatus.INCONSISTENT:
            raise MellinOperatorError("a weighted space cannot use inconsistent assumptions.")
        if p.is_real is False or p.is_finite is False:
            raise MellinOperatorError("p must be real and finite.")
        if p.is_number:
            if sp.StrictGreaterThan(p, 1) is not sp.true:
                raise MellinOperatorError("p must satisfy 1<p<infinity.")
        elif not (
            (p - 1).is_positive is True
            or _context_has(self.assumptions, sp.StrictGreaterThan(p, 1))
        ):
            raise MellinOperatorError("the explicit assumption p>1 is required.")
        if delta.is_real is False or delta.is_finite is False:
            raise MellinOperatorError("weight_exponent must be real and finite.")
        if not isinstance(self.norm_convention, WeightNormConvention):
            raise TypeError("norm_convention must be a WeightNormConvention.")
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(self, "measure", _nonempty(self.measure, "measure"))
        object.__setattr__(self, "p", p)
        object.__setattr__(self, "weight_exponent", delta)

    @property
    def normalization_exponent(self) -> sp.Expr:
        """Return the exponent derived from the actual norm."""

        return self.weight_exponent + 1 / self.p

    @property
    def kappa(self) -> sp.Expr:
        """Alias the article's normalized-dilation exponent ``delta+1/p``."""

        return self.normalization_exponent


@dataclass(frozen=True)
class WeightedDilationOperator:
    """A positive pullback dilation with every normalization made explicit."""

    atom: OperatorAtom
    gamma: sp.Expr
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    assumptions: AssumptionContext
    source: str
    normalization_power: sp.Expr | None = None
    convention: DilationConvention = DilationConvention.ARGUMENT_MULTIPLICATION
    evidence: tuple[object, ...] = ()
    bounded: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.domain, WeightedLpSpace) or not isinstance(
            self.codomain, WeightedLpSpace
        ):
            raise TypeError("domain and codomain must be WeightedLpSpace objects.")
        if self.domain != self.codomain:
            raise MellinOperatorError("a scalar dilation must preserve one weighted space.")
        if not isinstance(self.assumptions, AssumptionContext):
            raise TypeError("assumptions must be an AssumptionContext.")
        combined = self.domain.assumptions.combine(self.assumptions)
        if combined.consistency_status is ConsistencyStatus.INCONSISTENT:
            raise MellinOperatorError("dilation assumptions are inconsistent.")
        gamma = _sympy_scalar(self.gamma, "gamma")
        _require_positive(gamma, combined, "gamma")
        power = (
            self.domain.normalization_exponent
            if self.normalization_power is None
            else _sympy_scalar(self.normalization_power, "normalization_power")
        )
        if not isinstance(self.convention, DilationConvention):
            raise TypeError("convention must be a DilationConvention.")
        if not isinstance(self.bounded, bool):
            raise TypeError("bounded must be a bool.")
        object.__setattr__(self, "gamma", gamma)
        object.__setattr__(self, "normalization_power", power)
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))

    @property
    def normalization_factor(self) -> sp.Expr:
        return sp.Pow(self.gamma, self.normalization_power, evaluate=False)

    @property
    def norm_factor(self) -> sp.Expr:
        """Return ``||V_gamma||`` from the direct norm calculation."""

        return sp.Pow(
            self.gamma,
            self.normalization_power - self.domain.normalization_exponent,
            evaluate=False,
        )

    @property
    def is_isometric_normalization(self) -> bool:
        return _equivalent(
            self.normalization_power,
            self.domain.normalization_exponent,
        )

    def inverse(self, atom: OperatorAtom | None = None) -> WeightedDilationOperator:
        """Return the exact inverse using the same normalization power."""

        inverse_gamma = sp.simplify(1 / self.gamma)
        inverse_assumptions = self.assumptions
        combined = self.domain.assumptions.combine(inverse_assumptions)
        if not _is_explicitly_positive(inverse_gamma, combined):
            inverse_assumptions = inverse_assumptions.with_assumption(
                sp.StrictGreaterThan(inverse_gamma, 0)
            )
        return WeightedDilationOperator(
            atom=atom or OperatorAtom(f"{self.atom.name}_inverse", kind=self.atom.kind),
            gamma=inverse_gamma,
            domain=self.codomain,
            codomain=self.domain,
            assumptions=inverse_assumptions,
            source=f"inverse of {self.source}",
            normalization_power=self.normalization_power,
            convention=self.convention,
            evidence=self.evidence,
            bounded=self.bounded,
        )


@dataclass(frozen=True)
class MultiplicationOperator:
    atom: OperatorAtom
    symbol: sp.Expr
    radial_variable: sp.Symbol
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    source: str
    evidence: tuple[object, ...] = ()
    bounded: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.symbol, sp.Expr):
            raise TypeError("symbol must be a SymPy expression.")
        if not isinstance(self.radial_variable, sp.Symbol):
            raise TypeError("radial_variable must be a SymPy Symbol.")
        _validate_operator_spaces(self.domain, self.codomain)
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))
        if not isinstance(self.bounded, bool):
            raise TypeError("bounded must be a bool.")


@dataclass(frozen=True)
class WienerHopfOperator:
    atom: OperatorAtom
    symbol: sp.Expr
    frequency_variable: sp.Symbol
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    symbol_class: str
    source: str
    frequency_scaling_stable: bool
    stability_evidence: tuple[object, ...]
    localized: bool = False
    localization_controlled: bool = False
    fourier_convention: str = FOURIER_CONVENTION
    evidence: tuple[object, ...] = ()
    bounded: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.symbol, sp.Expr):
            raise TypeError("symbol must be a SymPy expression.")
        if not isinstance(self.frequency_variable, sp.Symbol):
            raise TypeError("frequency_variable must be a SymPy Symbol.")
        if (
            self.symbol.free_symbols
            and self.frequency_variable not in self.symbol.free_symbols
        ):
            raise MellinOperatorError("a Wiener--Hopf symbol must use its frequency.")
        _validate_operator_spaces(self.domain, self.codomain)
        object.__setattr__(self, "symbol_class", _nonempty(self.symbol_class, "symbol_class"))
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(
            self,
            "fourier_convention",
            _nonempty(self.fourier_convention, "fourier_convention"),
        )
        stability = _items(self.stability_evidence, "stability_evidence")
        if self.frequency_scaling_stable and not stability:
            raise MellinOperatorError(
                "frequency-scaling stability requires explicit evidence."
            )
        object.__setattr__(self, "stability_evidence", stability)
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))
        if not all(
            isinstance(value, bool)
            for value in (
                self.frequency_scaling_stable,
                self.localized,
                self.localization_controlled,
                self.bounded,
            )
        ):
            raise TypeError("Wiener--Hopf status fields must be bool values.")


@dataclass(frozen=True)
class MellinConvolutionOperator:
    atom: OperatorAtom
    symbol: MellinSymbol
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    symbol_class: str
    source: str
    evidence: tuple[object, ...] = ()
    bounded: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.symbol, MellinSymbol):
            raise TypeError("symbol must be a MellinSymbol.")
        if self.symbol.dependency not in {
            MellinSymbolDependency.FREQUENCY_ONLY,
            MellinSymbolDependency.PARAMETRIC_CONSTANT,
            MellinSymbolDependency.CONSTANT,
        }:
            raise MellinOperatorError(
                "a Mellin convolution symbol cannot depend on the radial variable."
            )
        _validate_operator_spaces(self.domain, self.codomain)
        object.__setattr__(self, "symbol_class", _nonempty(self.symbol_class, "symbol_class"))
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))
        if not isinstance(self.bounded, bool):
            raise TypeError("bounded must be a bool.")


@dataclass(frozen=True)
class MellinPseudodifferentialOperator:
    atom: OperatorAtom
    symbol: MellinSymbol
    domain: WeightedLpSpace
    codomain: WeightedLpSpace
    symbol_class: str
    source: str
    radial_scaling_stable: bool
    stability_evidence: tuple[object, ...]
    quantization: str = MELLIN_QUANTIZATION
    evidence: tuple[object, ...] = ()
    bounded: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.atom, OperatorAtom):
            raise TypeError("atom must be an OperatorAtom.")
        if not isinstance(self.symbol, MellinSymbol):
            raise TypeError("symbol must be a MellinSymbol.")
        if self.symbol.dependency is not MellinSymbolDependency.SPACE_FREQUENCY:
            raise MellinOperatorError(
                "a radial Mellin PDO must have SPACE_FREQUENCY dependency."
            )
        if len(self.symbol.spatial_variables) != 1:
            raise MellinOperatorError("a scalar radial Mellin PDO needs one spatial variable.")
        _validate_operator_spaces(self.domain, self.codomain)
        object.__setattr__(self, "symbol_class", _nonempty(self.symbol_class, "symbol_class"))
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(self, "quantization", _nonempty(self.quantization, "quantization"))
        stability = _items(self.stability_evidence, "stability_evidence")
        if self.radial_scaling_stable and not stability:
            raise MellinOperatorError(
                "radial-scaling stability requires explicit evidence."
            )
        for item in stability:
            try:
                hash(item)
            except TypeError as exc:
                raise TypeError(
                    "Mellin-PDO stability_evidence items must be hashable."
                ) from exc
        object.__setattr__(self, "stability_evidence", stability)
        object.__setattr__(self, "evidence", _items(self.evidence, "evidence"))
        if not isinstance(self.radial_scaling_stable, bool) or not isinstance(
            self.bounded, bool
        ):
            raise TypeError("Mellin PDO status fields must be bool values.")


def _validate_operator_spaces(
    domain: WeightedLpSpace,
    codomain: WeightedLpSpace,
) -> None:
    if not isinstance(domain, WeightedLpSpace) or not isinstance(
        codomain, WeightedLpSpace
    ):
        raise TypeError("domain and codomain must be WeightedLpSpace objects.")


@dataclass(frozen=True)
class RadialScalingEvidence:
    """Internal evidence recording the orientation of a radial scale."""

    gamma: sp.Expr
    radial_variable: sp.Symbol
    original_expression: sp.Expr
    scaled_expression: sp.Expr
    symbol_class: str
    stability_evidence: tuple[object, ...]
    orientation: str = "V_gamma Op(a) V_gamma_inverse gives a(gamma*x,lambda)"


def scale_radial_symbol(
    symbol: MellinSymbol,
    gamma: sp.Expr,
    assumptions: AssumptionContext,
    *,
    symbol_class: str,
    stability_evidence: tuple[object, ...],
) -> MellinSymbol:
    """Return ``a(gamma*x,lambda)`` without violating variable roles.

    ``MellinExpression.substitute`` deliberately rejects ``x -> gamma*x``
    because that generic operation mixes a spatial role with a parameter.
    This function reconstructs the typed domain explicitly: ``x`` remains the
    sole spatial variable and every new free symbol in ``gamma`` is declared
    as a parameter.
    """

    if not isinstance(symbol, MellinSymbol):
        raise TypeError("symbol must be a MellinSymbol.")
    if not isinstance(assumptions, AssumptionContext):
        raise TypeError("assumptions must be an AssumptionContext.")
    symbol_class = _nonempty(symbol_class, "symbol_class")
    stability = _items(stability_evidence, "stability_evidence")
    if not stability:
        raise MellinOperatorError(
            "radial symbol scaling requires explicit symbol-class stability evidence."
        )
    for item in stability:
        try:
            hash(item)
        except TypeError as exc:
            raise TypeError("stability_evidence items must be hashable.") from exc
    gamma = _sympy_scalar(gamma, "gamma")
    context = symbol.domain.assumption_context.combine(assumptions)
    if context.consistency_status is ConsistencyStatus.INCONSISTENT or (
        sp.simplify(context.as_boolean_expression()) is sp.false
    ):
        raise MellinOperatorError(
            "radial scaling assumptions are explicitly inconsistent."
        )
    _require_positive(gamma, context, "gamma")
    spatial = symbol.spatial_variables
    if symbol.dependency is not MellinSymbolDependency.SPACE_FREQUENCY or len(spatial) != 1:
        raise MellinOperatorError(
            "radial scaling requires one SPACE_FREQUENCY Mellin symbol."
        )
    radial = spatial[0]
    if radial.role is not MellinVariableRole.OUTPUT_SPATIAL:
        raise MellinOperatorError("radial scaling requires an output spatial variable.")
    if radial.symbol in gamma.free_symbols or symbol.frequency.symbol in gamma.free_symbols:
        raise MellinOperatorError("gamma cannot depend on a spatial or frequency variable.")

    declared = {item.symbol: item for item in symbol.domain.declared_variables}
    parameters = list(symbol.domain.parameters)
    for parameter_symbol in sorted(gamma.free_symbols, key=sp.srepr):
        existing = declared.get(parameter_symbol)
        if existing is not None:
            if existing.role is not MellinVariableRole.PARAMETER:
                raise MellinDomainRoleError(
                    "each free symbol of gamma must have PARAMETER role."
                )
            continue
        parameters.append(
            mellin_parameter(
                parameter_symbol,
                description="positive constant radial-scaling parameter",
            )
        )

    scaled_expression = symbol.as_expr().xreplace(
        {radial.symbol: gamma * radial.symbol}
    )
    scaling_evidence = RadialScalingEvidence(
        gamma,
        radial.symbol,
        symbol.as_expr(),
        scaled_expression,
        symbol_class,
        stability,
    )
    domain = MellinSymbolDomain(
        frequency=symbol.domain.frequency,
        complex_domain=symbol.domain.complex_domain,
        spatial_variables=symbol.domain.spatial_variables,
        relative_variable=symbol.domain.relative_variable,
        parameters=tuple(parameters),
        assumption_context=context,
        singular_set=symbol.domain.singular_set,
        description=symbol.domain.description,
        evidence=(*symbol.domain.evidence, scaling_evidence),
    )
    detected = detect_singularities(
        scaled_expression,
        domain.frequency.symbol,
        assumptions=domain.assumption_context,
        domain=domain.complex_domain,
    )
    domain = domain.with_singularities(detected)
    expression = MellinExpression(
        expression=scaled_expression,
        domain=domain,
        variables=symbol.expression.variables,
        parameters=tuple(parameters),
        dummy_variables=symbol.expression.dummy_variables,
        verification_status=symbol.verification_status,
        evidence=(*symbol.evidence, scaling_evidence),
        description=(
            f"Radial scaling of {symbol.description}"
            if symbol.description is not None
            else "Radially scaled Mellin symbol"
        ),
    )
    return MellinSymbol(
        expression,
        symbol.dependency,
        expression.description,
    )


@dataclass(frozen=True)
class DilationCovarianceEvidence:
    """Definition-level calculation supporting one covariance identity."""

    operator_family: str
    convention: str
    substitution: str
    normalization_check: str
    orientation_check: str
    operator_source: str
    symbol_class: str | None
    stability_evidence: tuple[object, ...]


@dataclass(frozen=True)
class DilationCovarianceTrace:
    """One exact covariance identity tied to authoritative AST products."""

    source_product: Product
    result_product: Product
    dilation: WeightedDilationOperator
    inverse_dilation: WeightedDilationOperator
    original_operator: object
    transformed_operator: object
    assumptions: AssumptionContext
    exact_identity: ExactIdentity
    rule: OperatorRule

    def __post_init__(self) -> None:
        if not isinstance(self.assumptions, AssumptionContext):
            raise TypeError("assumptions must be an AssumptionContext.")
        if not isinstance(self.dilation, WeightedDilationOperator) or not isinstance(
            self.inverse_dilation, WeightedDilationOperator
        ):
            raise TypeError(
                "dilation and inverse_dilation must be WeightedDilationOperator objects."
            )
        original_atom = getattr(self.original_operator, "atom", None)
        transformed_atom = getattr(self.transformed_operator, "atom", None)
        if not isinstance(original_atom, OperatorAtom) or not isinstance(
            transformed_atom, OperatorAtom
        ):
            raise TypeError(
                "stored covariance operators must expose authoritative atoms."
            )
        expected_source = Product(
            (self.dilation.atom, original_atom, self.inverse_dilation.atom)
        )
        expected_result = Product((transformed_atom,))
        if self.source_product != expected_source:
            raise MellinOperatorError(
                "covariance source AST does not match its stored operators."
            )
        if self.result_product != expected_result:
            raise MellinOperatorError(
                "covariance result AST does not match its transformed operator."
            )
        if self.exact_identity.left != self.source_product:
            raise MellinOperatorError("covariance identity lost its source AST.")
        if self.exact_identity.right != self.result_product:
            raise MellinOperatorError("covariance identity lost its result AST.")
        if self.rule.relation is not self.exact_identity:
            raise MellinOperatorError("covariance rule must reference its identity.")
        if (
            self.rule.relation_kind is not DerivationRelationKind.EXACT_EQUALITY
            or self.rule.certification_status
            is not RuleCertificationStatus.CERTIFIED_EXACT
        ):
            raise MellinOperatorError(
                "a dilation covariance trace requires a certified exact rule."
            )


def _require_inverse_pair(
    dilation: WeightedDilationOperator,
    inverse: WeightedDilationOperator,
) -> tuple[WeightedLpSpace, AssumptionContext]:
    if not isinstance(dilation, WeightedDilationOperator) or not isinstance(
        inverse, WeightedDilationOperator
    ):
        raise TypeError("dilation and inverse must be WeightedDilationOperator objects.")
    if not _equivalent(dilation.gamma * inverse.gamma, sp.Integer(1)):
        raise MellinOperatorError("the two dilation scales are not inverse.")
    if dilation.domain != inverse.domain or dilation.codomain != inverse.codomain:
        raise MellinOperatorError("inverse dilations must act on the same weighted space.")
    if dilation.convention is not inverse.convention:
        raise MellinOperatorError("inverse dilations use different action conventions.")
    if not _equivalent(dilation.normalization_power, inverse.normalization_power):
        raise MellinOperatorError("inverse dilations use different normalizations.")
    if not dilation.bounded or not inverse.bounded:
        raise MellinOperatorError("covariance requires bounded dilation factors.")
    context = (
        dilation.domain.assumptions.combine(dilation.assumptions).combine(
            inverse.assumptions
        )
    )
    if context.consistency_status is ConsistencyStatus.INCONSISTENT or (
        sp.simplify(context.as_boolean_expression()) is sp.false
    ):
        raise MellinOperatorError(
            "inverse dilation assumptions are explicitly incompatible."
        )
    return dilation.domain, context


def _require_bounded_operator(operator: object, family: str) -> None:
    if getattr(operator, "bounded", None) is not True:
        raise MellinOperatorError(
            f"exact {family} dilation covariance requires a bounded operator."
        )


def _covariance_trace(
    *,
    family: str,
    dilation: WeightedDilationOperator,
    operator_atom: OperatorAtom,
    inverse: WeightedDilationOperator,
    transformed_atom: OperatorAtom,
    original_operator: object,
    transformed_operator: object,
    assumptions: AssumptionContext,
    substitution: str,
    orientation: str,
) -> DilationCovarianceTrace:
    source = Product((dilation.atom, operator_atom, inverse.atom))
    result = Product((transformed_atom,))
    evidence = DilationCovarianceEvidence(
        family,
        dilation.convention.value,
        substitution,
        "gamma**nu and gamma**(-nu) cancel exactly",
        orientation,
        getattr(original_operator, "source", "direct definition-level calculation"),
        getattr(original_operator, "symbol_class", None),
        tuple(getattr(original_operator, "stability_evidence", ())),
    )
    identity = ExactIdentity(
        source,
        result,
        evidence=evidence,
        scope=ExactIdentityScope.OPERATORIAL,
        hypotheses=(
            *assumptions.assumptions,
            "all displayed operators are bounded",
            "the stated domain/codomain and normalization metadata agree",
        ),
    )
    rule = OperatorRule(
        name=f"exact_{family}_dilation_covariance",
        relation_kind=DerivationRelationKind.EXACT_EQUALITY,
        relation=identity,
        preconditions=(
            "the two dilations are inverse",
            "domain, p, weight, normalization, and convention agree",
            "gamma>0",
        ),
        source="direct calculation from repository Fourier/Mellin definitions",
        certification_status=RuleCertificationStatus.CERTIFIED_EXACT,
    )
    return DilationCovarianceTrace(
        source,
        result,
        dilation,
        inverse,
        original_operator,
        transformed_operator,
        assumptions,
        identity,
        rule,
    )


def conjugate_multiplication_by_dilation(
    dilation: WeightedDilationOperator,
    operator: MultiplicationOperator,
    inverse: WeightedDilationOperator,
) -> DilationCovarianceTrace:
    """Certify ``V_gamma M_a V_gamma^-1=M_{a(gamma*x)}``."""

    space, context = _require_inverse_pair(dilation, inverse)
    if not isinstance(operator, MultiplicationOperator):
        raise TypeError("operator must be a MultiplicationOperator.")
    _require_bounded_operator(operator, "multiplication")
    _require_same_operator_space(space, operator.domain, operator.codomain)
    scaled = operator.symbol.xreplace(
        {operator.radial_variable: dilation.gamma * operator.radial_variable}
    )
    transformed = MultiplicationOperator(
        _scaled_atom(operator.atom, "radially_scaled"),
        scaled,
        operator.radial_variable,
        operator.domain,
        operator.codomain,
        f"{operator.source}; exact dilation covariance of multiplication",
        (*operator.evidence, "direct pointwise calculation"),
        operator.bounded,
    )
    return _covariance_trace(
        family="multiplication",
        dilation=dilation,
        operator_atom=operator.atom,
        inverse=inverse,
        transformed_atom=transformed.atom,
        original_operator=operator,
        transformed_operator=transformed,
        assumptions=context,
        substitution="evaluate a at gamma*x",
        orientation="a_gamma(x)=a(gamma*x)",
    )


def conjugate_wiener_hopf_by_dilation(
    dilation: WeightedDilationOperator,
    operator: WienerHopfOperator,
    inverse: WeightedDilationOperator,
) -> DilationCovarianceTrace:
    """Certify ``V_gamma W(b) V_gamma^-1=W(b(lambda/gamma))``."""

    space, context = _require_inverse_pair(dilation, inverse)
    if not isinstance(operator, WienerHopfOperator):
        raise TypeError("operator must be a WienerHopfOperator.")
    _require_bounded_operator(operator, "Wiener--Hopf")
    _require_same_operator_space(space, operator.domain, operator.codomain)
    if operator.fourier_convention != FOURIER_CONVENTION:
        raise MellinOperatorError(
            "Wiener--Hopf covariance requires the repository Fourier convention."
        )
    if operator.localized:
        raise MellinOperatorError(
            "localized Wiener--Hopf covariance is not certified: the current "
            "type cannot represent the required cutoff rescaling chi(x) -> "
            "chi(gamma*x)."
        )
    if not operator.frequency_scaling_stable:
        raise MellinOperatorError("the declared Wiener--Hopf class is not scale-stable.")
    lam = operator.frequency_variable
    scaled_symbol = operator.symbol.xreplace({lam: lam / dilation.gamma})
    transformed = WienerHopfOperator(
        atom=_scaled_atom(operator.atom, "frequency_scaled"),
        symbol=scaled_symbol,
        frequency_variable=lam,
        domain=operator.domain,
        codomain=operator.codomain,
        symbol_class=operator.symbol_class,
        source=(
            f"{operator.source}; Fourier scaling "
            "h_gamma(t)=gamma*h(gamma*t)"
        ),
        frequency_scaling_stable=True,
        stability_evidence=operator.stability_evidence,
        localized=operator.localized,
        localization_controlled=operator.localization_controlled,
        fourier_convention=operator.fourier_convention,
        evidence=(*operator.evidence, "positive half-line is scale invariant"),
        bounded=operator.bounded,
    )
    return _covariance_trace(
        family="wiener_hopf",
        dilation=dilation,
        operator_atom=operator.atom,
        inverse=inverse,
        transformed_atom=transformed.atom,
        original_operator=operator,
        transformed_operator=transformed,
        assumptions=context,
        substitution="y=gamma*z; h_gamma(t)=gamma*h(gamma*t)",
        orientation="b_gamma(lambda)=b(lambda/gamma)",
    )


def conjugate_mellin_convolution_by_dilation(
    dilation: WeightedDilationOperator,
    operator: MellinConvolutionOperator,
    inverse: WeightedDilationOperator,
) -> DilationCovarianceTrace:
    """Certify exact invariance of a radial-independent Mellin convolution."""

    space, context = _require_inverse_pair(dilation, inverse)
    if not isinstance(operator, MellinConvolutionOperator):
        raise TypeError("operator must be a MellinConvolutionOperator.")
    _require_bounded_operator(operator, "Mellin-convolution")
    _require_same_operator_space(space, operator.domain, operator.codomain)
    return _covariance_trace(
        family="mellin_convolution",
        dilation=dilation,
        operator_atom=operator.atom,
        inverse=inverse,
        transformed_atom=operator.atom,
        original_operator=operator,
        transformed_operator=operator,
        assumptions=context,
        substitution="y=gamma*z with dy/y=dz/z",
        orientation="c(lambda) is independent of the radial variable",
    )


def conjugate_mellin_pdo_by_dilation(
    dilation: WeightedDilationOperator,
    operator: MellinPseudodifferentialOperator,
    inverse: WeightedDilationOperator,
    assumptions: AssumptionContext,
) -> DilationCovarianceTrace:
    """Certify ``V_gamma Op_M(a) V_gamma^-1=Op_M(a(gamma*x,lambda))``."""

    space, pair_context = _require_inverse_pair(dilation, inverse)
    if not isinstance(operator, MellinPseudodifferentialOperator):
        raise TypeError("operator must be a MellinPseudodifferentialOperator.")
    _require_bounded_operator(operator, "Mellin-PDO")
    if not isinstance(assumptions, AssumptionContext):
        raise TypeError("assumptions must be an AssumptionContext.")
    context = pair_context.combine(assumptions).combine(operator.symbol.assumptions)
    if context.consistency_status is ConsistencyStatus.INCONSISTENT or (
        sp.simplify(context.as_boolean_expression()) is sp.false
    ):
        raise MellinOperatorError(
            "Mellin-PDO covariance assumptions are explicitly incompatible."
        )
    _require_same_operator_space(space, operator.domain, operator.codomain)
    if operator.quantization != MELLIN_QUANTIZATION:
        raise MellinOperatorError(
            "Mellin-PDO covariance requires the repository Mellin quantization."
        )
    if not operator.radial_scaling_stable:
        raise MellinOperatorError("the Mellin symbol class is not radially scale-stable.")
    scaled_symbol = scale_radial_symbol(
        operator.symbol,
        dilation.gamma,
        context,
        symbol_class=operator.symbol_class,
        stability_evidence=operator.stability_evidence,
    )
    transformed = MellinPseudodifferentialOperator(
        atom=_scaled_atom(operator.atom, "radially_scaled"),
        symbol=scaled_symbol,
        domain=operator.domain,
        codomain=operator.codomain,
        symbol_class=operator.symbol_class,
        source=f"{operator.source}; exact Mellin-PDO dilation covariance",
        radial_scaling_stable=True,
        stability_evidence=operator.stability_evidence,
        quantization=operator.quantization,
        evidence=(*operator.evidence, "dy/y is invariant under y=gamma*z"),
        bounded=operator.bounded,
    )
    return _covariance_trace(
        family="mellin_pdo",
        dilation=dilation,
        operator_atom=operator.atom,
        inverse=inverse,
        transformed_atom=transformed.atom,
        original_operator=operator,
        transformed_operator=transformed,
        assumptions=context,
        substitution="evaluate at gamma*x and set y=gamma*z; dy/y=dz/z",
        orientation="a_gamma(x,lambda)=a(gamma*x,lambda)",
    )


def _require_same_operator_space(
    expected: WeightedLpSpace,
    domain: WeightedLpSpace,
    codomain: WeightedLpSpace,
) -> None:
    if domain != expected or codomain != expected:
        raise MellinOperatorError("operator domain/codomain do not match the dilation space.")


@dataclass(frozen=True)
class RegularizerMellinRepresentation:
    """A separate, never implicit representation of an abstract regularizer."""

    regularizer: RegularizerOperator
    mellin_operator: MellinPseudodifferentialOperator
    status: RegularizerRepresentationStatus
    hypotheses: tuple[object, ...]
    source: str
    evidence: tuple[object, ...] = ()
    space: WeightedLpSpace | None = None
    compact_ideal: object | None = None
    semantic_relation: ExactIdentity | FormalIdentity | ModCompactEquivalence = field(
        init=False
    )

    def __post_init__(self) -> None:
        if not isinstance(self.regularizer, RegularizerOperator):
            raise TypeError("regularizer must be a RegularizerOperator.")
        if not isinstance(self.regularizer.operator, OperatorAtom):
            raise TypeError(
                "the abstract regularizer must retain an authoritative OperatorAtom."
            )
        if not isinstance(self.mellin_operator, MellinPseudodifferentialOperator):
            raise TypeError("mellin_operator must be a MellinPseudodifferentialOperator.")
        if not isinstance(self.status, RegularizerRepresentationStatus):
            raise TypeError("status must be a RegularizerRepresentationStatus.")
        hypotheses = _items(self.hypotheses, "hypotheses")
        evidence = _items(self.evidence, "evidence")
        if self.status is not RegularizerRepresentationStatus.ASSUMED and not evidence:
            raise MellinOperatorError("a certified representation requires evidence.")
        if self.status is RegularizerRepresentationStatus.ASSUMED and not hypotheses:
            raise MellinOperatorError("an assumed representation requires a hypothesis.")
        source = _nonempty(self.source, "source")
        if self.space is not None and self.space != self.mellin_operator.domain:
            raise MellinOperatorError("representation space and Mellin operator disagree.")
        relation: ExactIdentity | FormalIdentity | ModCompactEquivalence
        if self.status is RegularizerRepresentationStatus.CERTIFIED_EXACT:
            relation = ExactIdentity(
                self.regularizer.operator,
                self.mellin_operator.atom,
                evidence=evidence,
                scope=ExactIdentityScope.OPERATORIAL,
                hypotheses=hypotheses,
            )
        elif self.status is RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT:
            relation = ModCompactEquivalence(
                self.regularizer.operator,
                self.mellin_operator.atom,
                space=self.space or self.mellin_operator.domain,
                compact_ideal=self.compact_ideal or "compact operators",
                evidence=evidence,
            )
        else:
            relation = FormalIdentity(
                self.regularizer.operator,
                self.mellin_operator.atom,
                justification=(source, hypotheses, evidence),
            )
        object.__setattr__(self, "hypotheses", hypotheses)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "semantic_relation", relation)

    @property
    def relation_kind(self) -> DerivationRelationKind:
        if self.status is RegularizerRepresentationStatus.CERTIFIED_MOD_COMPACT:
            return DerivationRelationKind.MOD_COMPACT_EQUIVALENCE
        if self.status is RegularizerRepresentationStatus.ASSUMED:
            return DerivationRelationKind.FORMAL_SUBSTITUTION
        return DerivationRelationKind.EXACT_EQUALITY


@dataclass(frozen=True)
class AlgebraMembershipClaim:
    """A membership claim whose status is never inferred from factor types."""

    expression: object
    algebra: str
    status: AlgebraMembershipStatus
    assumptions: tuple[object, ...]
    source: str
    evidence: tuple[object, ...] = ()

    def __post_init__(self) -> None:
        if self.expression is None:
            raise ValueError("expression must be supplied.")
        object.__setattr__(self, "algebra", _nonempty(self.algebra, "algebra"))
        if not isinstance(self.status, AlgebraMembershipStatus):
            raise TypeError("status must be an AlgebraMembershipStatus.")
        assumptions = _items(self.assumptions, "assumptions")
        evidence = _items(self.evidence, "evidence")
        if self.status is AlgebraMembershipStatus.PROVED and not evidence:
            raise MellinOperatorError("PROVED algebra membership requires evidence.")
        object.__setattr__(self, "assumptions", assumptions)
        object.__setattr__(self, "source", _nonempty(self.source, "source"))
        object.__setattr__(self, "evidence", evidence)


__all__ = [
    "AlgebraMembershipClaim",
    "AlgebraMembershipStatus",
    "DilationConvention",
    "DilationCovarianceEvidence",
    "DilationCovarianceTrace",
    "FOURIER_CONVENTION",
    "MELLIN_QUANTIZATION",
    "MellinConvolutionOperator",
    "MellinOperatorError",
    "MellinPseudodifferentialOperator",
    "MultiplicationOperator",
    "RadialScalingEvidence",
    "RegularizerMellinRepresentation",
    "RegularizerRepresentationStatus",
    "WeightNormConvention",
    "WeightedDilationOperator",
    "WeightedLpSpace",
    "WienerHopfOperator",
    "conjugate_mellin_convolution_by_dilation",
    "conjugate_mellin_pdo_by_dilation",
    "conjugate_multiplication_by_dilation",
    "conjugate_wiener_hopf_by_dilation",
    "scale_radial_symbol",
]
