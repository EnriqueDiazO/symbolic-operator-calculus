"""LaTeX presentation of the structured first Schur derivation trace."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import sympy as sp
from sympy.printing.latex import LatexPrinter

from .actions import AppliedLinearCombination, PrincipalValue
from .conditional import ConditionalIdentity
from .complete_correction import CompleteCorrectionDerivation
from .derivations import FirstSchurDerivationTrace
from .domains import AssumptionContext, ComplexDomain, DomainRegionKind
from .kernels import KernelCombination
from .normalized_schur import (
    NormalizedFirstSchurPivotDerivation,
    NormalizedModCompactFactorization,
)
from .operators import (
    G1,
    G2,
    I,
    LinearCombination,
    OperatorAtom,
    Product,
    Scalar,
    Term,
)
from .relations import ExactBlock, FirstSchurReduction
from .relative_wiener_hopf import RelativeWienerHopfDerivationTrace
from .semantics import (
    ApproximateEquality,
    ExactIdentity,
    FormalIdentity,
    KernelAnnotatedExpression,
    ModCompactEquivalence,
)
from .singularities import SingularSet


OPERATOR_LATEX: dict[str, str] = {
    "A11_inverse": r"A_{1,1}^{(-1)}",
    "A12": r"A_{1,2}",
    "A21": r"A_{2,1}",
    "A22": r"A_{2,2}",
    "A22_first": r"A_{2,2}^{(1)}",
    "C2": r"\mathcal C_{2}",
    "C2_normalized": r"\mathcal C_{2}^{\mathrm{norm}}",
    "G1": r"G_{1}\,I",
    "G2": r"G_{2}\,I",
    "Ghat1": r"\widehat G_{1}",
    "Ghat2": r"\widehat G_{2}",
    "I": r"I",
    "N2": r"N_{2}",
    "N2_first": r"N_{2}^{(1)}",
    "P_minus": r"P^{-}",
    "P_plus": r"P^{+}",
    "R1": r"R_{1}",
    "R11": r"R_{1,1}",
    "S_Rplus": r"S_{\mathbb R_+}",
    "T1_minus": r"T_{1,-}",
    "T2_minus": r"T_{2,-}",
    "U1": r"U_{1}",
    "U2": r"U_{2}",
    "U1_inverse": r"U_{1}^{-1}",
    "U2_inverse": r"U_{2}^{-1}",
    "Vtilde_alpha1": r"\widetilde V_{\alpha_1}",
    "Vtilde_alpha2": r"\widetilde V_{\alpha_2}",
    "Wminus_21": r"W^-_{2,1}",
    "Wplus_12": r"W^+_{1,2}",
    "Z1_inverse": r"Z_{1}^{-1}",
    "Z2_inverse": r"Z_{2}^{-1}",
}

SCALAR_FUNCTION_LATEX: dict[str, str] = {
    "G1": r"G_{1}",
    "G2": r"G_{2}",
    "Lminus_21": r"L^-_{2,1}",
    "Lplus_12": r"L^+_{1,2}",
    "R11": r"R_{1,1}",
    "chi_infinity": r"\chi_\infty",
    "rho1": r"\rho_1",
    "rho2": r"\rho_2",
    "chi_plus": r"\chi_+",
    "chi_minus": r"\chi_-",
}

SCALAR_SYMBOL_LATEX: dict[str, str] = {
    "gamma1": r"\gamma_1",
    "gamma2": r"\gamma_2",
}


class LatexRenderingError(ValueError):
    """Raised when an object has no supported K7B LaTeX representation."""


def _semantic_endpoint_latex(value: object) -> str:
    if isinstance(value, sp.Basic):
        return _StructuredScalarLatexPrinter().doprint(value)
    escaped = r"\_".join(str(value).split("_"))
    return rf"\text{{{escaped}}}"


def render_assumption_context_latex(context: AssumptionContext) -> str:
    """Render every stored condition without claiming it has been proved."""

    if not isinstance(context, AssumptionContext):
        raise TypeError("context must be an AssumptionContext.")
    if not context.assumptions:
        return r"\text{no explicit assumptions}"
    conditions = r" \\ ".join(sp.latex(item) for item in context.assumptions)
    return rf"\substack{{{conditions}}}"


def render_complex_domain_latex(domain: ComplexDomain) -> str:
    """Render the declared region, assumptions, and exclusions."""

    if not isinstance(domain, ComplexDomain):
        raise TypeError("domain must be a ComplexDomain.")
    variable = sp.latex(domain.variable)
    regions: list[str] = []
    for region in domain.regions:
        if region.kind is DomainRegionKind.COMPLEX_PLANE:
            regions.append(rf"{variable}\in\mathbb{{C}}")
        elif region.kind is DomainRegionKind.REAL_LINE:
            regions.append(rf"{variable}\in\mathbb{{R}}")
        elif region.kind is DomainRegionKind.UPPER_HALF_PLANE:
            regions.append(rf"\operatorname{{Im}}{variable}>0")
        elif region.kind is DomainRegionKind.LOWER_HALF_PLANE:
            regions.append(rf"\operatorname{{Im}}{variable}<0")
        elif region.kind is DomainRegionKind.VERTICAL_STRIP:
            regions.append(
                rf"{sp.latex(region.lower)}<\operatorname{{Re}}{variable}"
                rf"<{sp.latex(region.upper)}"
            )
        elif region.kind is DomainRegionKind.HORIZONTAL_STRIP:
            regions.append(
                rf"{sp.latex(region.lower)}<\operatorname{{Im}}{variable}"
                rf"<{sp.latex(region.upper)}"
            )
        elif region.condition is not None:
            regions.append(sp.latex(region.condition))
    if domain.exclusions:
        for exclusion in domain.exclusions:
            if isinstance(exclusion, sp.Set):
                regions.append(rf"{variable}\notin{sp.latex(exclusion)}")
            else:
                regions.append(rf"{variable}\ne{sp.latex(exclusion)}")
    regions.append(render_assumption_context_latex(domain.assumption_context))
    return r"\;\land\;".join(regions)


def render_singular_set_latex(singular_set: SingularSet) -> str:
    """Render stored singularities, conditions, and origin labels."""

    if not isinstance(singular_set, SingularSet):
        raise TypeError("singular_set must be a SingularSet.")
    if not singular_set.singularities:
        return r"\operatorname{Sing}_{\mathrm{declared}}=\varnothing"
    entries: list[str] = []
    for singularity in singular_set.singularities:
        conditions = render_assumption_context_latex(singularity.conditions)
        order = "" if singularity.order is None else rf",\;m={singularity.order}"
        location = (
            rf"{sp.latex(singularity.variable)}\in{sp.latex(singularity.location)}"
            if isinstance(singularity.location, sp.Set)
            else rf"{sp.latex(singularity.variable)}={sp.latex(singularity.location)}"
        )
        entries.append(
            rf"\text{{{' '.join(singularity.kind.value.split('_'))}}}:\;"
            rf"{location}"
            rf"{order}\;[\text{{{' '.join(singularity.origin.value.split('_'))}}};\;"
            rf"{conditions}]"
        )
    return r"\substack{" + r" \\ ".join(entries) + "}"


def render_conditional_identity_latex(identity: ConditionalIdentity) -> str:
    """Render a conditional scalar equality with all validity metadata visible."""

    if not isinstance(identity, ConditionalIdentity):
        raise TypeError("identity must be a ConditionalIdentity.")
    lhs = _semantic_endpoint_latex(identity.lhs)
    rhs = _semantic_endpoint_latex(identity.rhs)
    conditions = render_assumption_context_latex(identity.assumption_context)
    domain = render_complex_domain_latex(identity.domain)
    singularities = render_singular_set_latex(identity.singular_set)
    status = " ".join(identity.verification_status.value.split("_"))
    return (
        rf"\left.{lhs}={rhs}\;\right|_{{{conditions}}}"
        rf"\quad\substack{{\text{{conditional scalar identity}} \\ "
        rf"\text{{scope: {' '.join(identity.scope.value.split('_'))}}} \\ "
        rf"\text{{status: {status}}} \\ {domain} \\ {singularities}}}"
    )


def render_semantic_relation_latex(relation: object) -> str:
    """Visually distinguish the semantic relation types introduced in P0."""

    if isinstance(relation, ConditionalIdentity):
        return render_conditional_identity_latex(relation)
    if isinstance(relation, ExactIdentity):
        return (
            rf"{_semantic_endpoint_latex(relation.left)}"
            rf"\overset{{\text{{exact: {relation.scope.value}}}}}{{=}}"
            rf"{_semantic_endpoint_latex(relation.right)}"
        )
    if isinstance(relation, FormalIdentity):
        return (
            rf"{_semantic_endpoint_latex(relation.left)}"
            rf"\overset{{\text{{formal}}}}{{=}}"
            rf"{_semantic_endpoint_latex(relation.right)}"
        )
    if isinstance(relation, ModCompactEquivalence):
        return (
            rf"{_semantic_endpoint_latex(relation.left)}"
            rf"\equiv_{{\mathrm{{mod}}\,\mathcal K}}"
            rf"{_semantic_endpoint_latex(relation.right)}"
        )
    if isinstance(relation, ApproximateEquality):
        return (
            rf"{_semantic_endpoint_latex(relation.left)}"
            rf"\overset{{\text{{approximate}}}}{{\approx}}"
            rf"{_semantic_endpoint_latex(relation.right)}"
        )
    raise LatexRenderingError("unsupported semantic relation type.")


class _StructuredScalarLatexPrinter(LatexPrinter):
    """Localized scalar printer with the project's mathematical notation."""

    def _print_Function(self, expr: sp.Function, exp: str | None = None) -> str:
        name = expr.func.__name__
        notation = SCALAR_FUNCTION_LATEX.get(name)
        if notation is None:
            return super()._print_Function(expr, exp)
        if exp is not None:
            notation = rf"{notation}^{{{exp}}}"
        arguments = ",".join(self._print(argument) for argument in expr.args)
        return rf"{notation}\!\left({arguments}\right)"

    def _print_Symbol(self, expr: sp.Symbol) -> str:
        notation = SCALAR_SYMBOL_LATEX.get(expr.name)
        return notation if notation is not None else super()._print_Symbol(expr)

    def _print_PrincipalValue(self, expr: PrincipalValue) -> str:
        return rf"\operatorname{{p.v.}}\,{self._print(expr.args[0])}"


@dataclass(frozen=True)
class RenderedDerivationStep:
    """One stable presentation step with separate title and LaTeX formula."""

    key: str
    title: str
    latex: str

    def __post_init__(self) -> None:
        if not self.key or not self.title or not self.latex:
            raise LatexRenderingError("rendered step fields must be non-empty.")


@dataclass(frozen=True)
class RenderedFirstSchurDerivation:
    """Ordered LaTeX presentation of the first reduced m=2 derivation."""

    steps: tuple[RenderedDerivationStep, ...]

    def __post_init__(self) -> None:
        if not self.steps or not all(
            isinstance(step, RenderedDerivationStep) for step in self.steps
        ):
            raise LatexRenderingError(
                "steps must be a non-empty tuple of rendered derivation steps."
            )
        keys = tuple(step.key for step in self.steps)
        if len(set(keys)) != len(keys):
            raise LatexRenderingError("rendered step keys must be unique.")

    def as_latex(self) -> str:
        """Concatenate the already-rendered formulas without document markup."""

        return "\n\n".join(step.latex for step in self.steps)


@dataclass(frozen=True)
class RenderedNormalizedFirstSchurPivot:
    """Ordered LaTeX presentation of the normalized first-pivot trace."""

    steps: tuple[RenderedDerivationStep, ...]
    proof_obligations: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.steps or not all(
            isinstance(step, RenderedDerivationStep) for step in self.steps
        ):
            raise LatexRenderingError(
                "steps must contain rendered normalized-pivot steps."
            )
        keys = tuple(step.key for step in self.steps)
        if len(set(keys)) != len(keys):
            raise LatexRenderingError("rendered step keys must be unique.")
        if not self.proof_obligations or not all(
            isinstance(item, str) and item.strip()
            for item in self.proof_obligations
        ):
            raise LatexRenderingError(
                "proof_obligations must contain non-empty strings."
            )

    def as_latex(self) -> str:
        """Concatenate formulas without adding document-level markup."""

        return "\n\n".join(step.latex for step in self.steps)


@dataclass(frozen=True)
class RenderedCompleteCorrectionDerivation:
    """Ordered presentation of exact normalization inside a mod-compact model."""

    steps: tuple[RenderedDerivationStep, ...]

    def __post_init__(self) -> None:
        if not self.steps or not all(
            isinstance(item, RenderedDerivationStep) for item in self.steps
        ):
            raise LatexRenderingError(
                "steps must contain RenderedDerivationStep objects."
            )
        keys = tuple(step.key for step in self.steps)
        if len(set(keys)) != len(keys):
            raise LatexRenderingError("complete-correction step keys must be unique.")

    def as_latex(self) -> str:
        return "\n\n".join(step.latex for step in self.steps)


@dataclass(frozen=True)
class RenderedRelativeWienerHopfTrace:
    """Ordered LaTeX presentation of one exact relative-WH trace."""

    steps: tuple[RenderedDerivationStep, ...]

    def as_latex(self) -> str:
        return "\n\n".join(step.latex for step in self.steps)


def render_relative_wiener_hopf_trace_latex(
    trace: RelativeWienerHopfDerivationTrace,
) -> RenderedRelativeWienerHopfTrace:
    """Render stored L1/L2 trace data without recomputing mathematics."""

    if not isinstance(trace, RelativeWienerHopfDerivationTrace):
        raise TypeError("trace must be a RelativeWienerHopfDerivationTrace.")
    l1 = trace.factorization
    k, j = l1.k, l1.j
    symbol = render_scalar_latex(l1.original_symbol)
    left_symbol = render_scalar_latex(l1.left_symbol)
    right_symbol = render_scalar_latex(l1.right_symbol)
    kernel = render_scalar_latex(l1.integral_kernel)
    scale = render_scalar_latex(trace.relative_dilation.scale)
    integral = render_scalar_latex(trace.action.integral)
    steps = (
        RenderedDerivationStep(
            "original_operator", "Original operator",
            rf"\mathcal W_{{{k},{j}}}=V_{{\alpha_{k}}}W\!\left({symbol}\right)V_{{\alpha_{j}}}^{{-1}}",
        ),
        RenderedDerivationStep(
            "conjugated_kernel", "Conjugated kernel",
            rf"\mathcal K_{{{k},{j}}}(x,y)={kernel}",
        ),
        RenderedDerivationStep(
            "relative_dilation", "Relative dilation",
            rf"\beta_{{{k},{j}}}=\alpha_{j}^{{-1}}\circ\alpha_{k},\qquad "
            rf"\beta_{{{k},{j}}}(x)={scale}x",
        ),
        RenderedDerivationStep(
            "left_symbol", "Left symbol",
            rf"b_{{{k},{j}}}^{{\mathrm L}}(\lambda)={left_symbol}",
        ),
        RenderedDerivationStep(
            "left_factorization", "Left factorization",
            rf"\mathcal W_{{{k},{j}}}=V_{{\beta_{{{k},{j}}}}}W\!\left(b_{{{k},{j}}}^{{\mathrm L}}\right)",
        ),
        RenderedDerivationStep(
            "right_symbol", "Right symbol",
            rf"b_{{{k},{j}}}^{{\mathrm R}}(\lambda)={right_symbol}",
        ),
        RenderedDerivationStep(
            "right_factorization", "Right factorization",
            rf"\mathcal W_{{{k},{j}}}=W\!\left(b_{{{k},{j}}}^{{\mathrm R}}\right)V_{{\beta_{{{k},{j}}}}}",
        ),
        RenderedDerivationStep(
            "exact_identity", "Exact identity within validated model",
            rf"\mathcal W_{{{k},{j}}}=V_{{\beta_{{{k},{j}}}}}W\!\left(b_{{{k},{j}}}^{{\mathrm L}}\right)"
            rf"=W\!\left(b_{{{k},{j}}}^{{\mathrm R}}\right)V_{{\beta_{{{k},{j}}}}}",
        ),
        RenderedDerivationStep(
            "symbol_correspondence", "Symbol correspondence",
            rf"b_{{{k},{j}}}^{{\mathrm R}}(\lambda)=b_{{{k},{j}}}^{{\mathrm L}}"
            rf"\!\left(\frac{{\gamma_{j}}}{{\gamma_{k}}}\lambda\right),\qquad "
            rf"b_{{{k},{j}}}^{{\mathrm L}}(\lambda)=b_{{{k},{j}}}^{{\mathrm R}}"
            rf"\!\left(\frac{{\gamma_{k}}}{{\gamma_{j}}}\lambda\right)",
        ),
        RenderedDerivationStep(
            "common_action", "Common action",
            rf"\left(\mathcal W_{{{k},{j}}}f\right)(x)=\int_0^\infty "
            rf"\mathcal K_{{{k},{j}}}(x,y)f(y)\,dy={integral}",
        ),
    )
    return RenderedRelativeWienerHopfTrace(steps)


def render_operator_atom_latex(atom: OperatorAtom) -> str:
    """Render one known noncommutative operator atom."""

    if not isinstance(atom, OperatorAtom):
        raise TypeError("atom must be an OperatorAtom.")
    try:
        return OPERATOR_LATEX[atom.name]
    except KeyError as exc:
        raise LatexRenderingError(
            f"no LaTeX notation is registered for operator {atom.name!r}."
        ) from exc


def render_product_latex(product: Product) -> str:
    """Render a product in its stored noncommutative factor order."""

    if not isinstance(product, Product):
        raise TypeError("product must be a Product.")
    if not product.factors:
        return render_operator_atom_latex(I)
    return r"\,".join(
        _render_operator_product_factor_latex(atom) for atom in product.factors
    )


def render_term_latex(term: Term, *, include_sign: bool = True) -> str:
    """Render one AST term without changing its product structure."""

    if not isinstance(term, Term):
        raise TypeError("term must be a Term.")
    sign, magnitude = _coefficient_sign_and_magnitude(term.coefficient)
    if magnitude == 0:
        return "0"
    body = render_product_latex(term.product)
    coefficient = _coefficient_magnitude_latex(magnitude)
    rendered = body if coefficient == "1" else rf"{coefficient}\,{body}"
    if include_sign and sign == -1:
        return rf"- {rendered}"
    return rendered


def render_linear_combination_latex(combination: LinearCombination) -> str:
    """Render terms and signs in their stored AST order."""

    if not isinstance(combination, LinearCombination):
        raise TypeError("combination must be a LinearCombination.")
    pieces: list[str] = []
    for index, term in enumerate(combination.terms):
        sign, _ = _coefficient_sign_and_magnitude(term.coefficient)
        body = render_term_latex(term, include_sign=False)
        if index == 0:
            pieces.append(rf"- {body}" if sign == -1 else body)
        else:
            pieces.append(rf" - {body}" if sign == -1 else rf" + {body}")
    return "".join(pieces) if pieces else "0"


def render_normalized_first_schur_pivot_exact_latex(
    trace: NormalizedFirstSchurPivotDerivation,
) -> str:
    """Render the final exact formal identity in noncommutative order."""

    if not isinstance(trace, NormalizedFirstSchurPivotDerivation):
        raise TypeError("trace must be a NormalizedFirstSchurPivotDerivation.")
    return (
        f"{_render_normalized_endpoint_latex(trace.exact_result.left)} = "
        f"{_render_normalized_endpoint_latex(trace.exact_result.right)}"
    )


def render_original_complete_correction_latex(
    trace: CompleteCorrectionDerivation,
) -> str:
    """Render the Phase N factorized definition of ``C2``."""

    if not isinstance(trace, CompleteCorrectionDerivation):
        raise TypeError("trace must be a CompleteCorrectionDerivation.")
    factorization = trace.original_factorization
    prefix = _render_normalized_endpoint_latex(factorization.prefix)
    left = _render_normalized_endpoint_latex(factorization.left_difference)
    left_wh = _render_normalized_atom_latex(factorization.left_wiener_hopf)
    bridge = _render_normalized_endpoint_latex(factorization.bridge)
    right = _render_normalized_endpoint_latex(factorization.right_difference)
    right_wh = _render_normalized_atom_latex(factorization.right_wiener_hopf)
    suffix = _render_normalized_endpoint_latex(factorization.suffix)
    return (
        r"\mathcal C_{2} = "
        rf"{prefix}\,\left({left}\right)\,{left_wh}\,{bridge}\,"
        rf"\left({right}\right)\,{right_wh}\,{suffix}"
    )


def render_normalized_complete_correction_latex(
    trace: CompleteCorrectionDerivation,
) -> str:
    """Render C0 from the validated normalized factorization metadata."""

    if not isinstance(trace, CompleteCorrectionDerivation):
        raise TypeError("trace must be a CompleteCorrectionDerivation.")
    factorization = trace.normalized_factorization
    left = _render_normalized_endpoint_latex(factorization.left_difference)
    left_wh = _render_normalized_atom_latex(factorization.left_wiener_hopf)
    bridge = _render_normalized_endpoint_latex(factorization.bridge)
    right = _render_normalized_endpoint_latex(factorization.right_difference)
    right_wh = _render_normalized_atom_latex(factorization.right_wiener_hopf)
    suffix = _render_normalized_endpoint_latex(factorization.suffix)
    return (
        r"\mathcal C_{2}^{\mathrm{norm}} = "
        rf"\left({left}\right)\,{left_wh}\,{bridge}\,"
        rf"\left({right}\right)\,{right_wh}\,{suffix}"
    )


def render_complete_correction_derivation_latex(
    trace: CompleteCorrectionDerivation,
) -> RenderedCompleteCorrectionDerivation:
    """Render B1--B4 with distinct exact and modulo-compact symbols."""

    if not isinstance(trace, CompleteCorrectionDerivation):
        raise TypeError("trace must be a CompleteCorrectionDerivation.")
    left_rule = trace.left_normalization_rule.relation
    right_rule = trace.right_normalization_rule.relation
    steps = (
        RenderedDerivationStep(
            "B1_pivot_correction",
            "Correction inherited from the Phase N modulo-compact model",
            r"N_{2}^{(1)} \simeq N_{2} + \mathcal C_{2}",
        ),
        RenderedDerivationStep(
            "B1_correction_definition",
            "Factorized complete correction",
            render_original_complete_correction_latex(trace),
        ),
        RenderedDerivationStep(
            "B2_left_normalization",
            "Certified exact left branch normalization",
            _render_prefixed_difference_identity_latex(left_rule),
        ),
        RenderedDerivationStep(
            "B3_right_normalization",
            "Certified exact right branch normalization",
            _render_prefixed_difference_identity_latex(right_rule),
        ),
        RenderedDerivationStep(
            "B4_normalized_correction",
            "Exact normalized correction inside the supplied model",
            render_normalized_complete_correction_latex(trace),
        ),
    )
    return RenderedCompleteCorrectionDerivation(steps)


def render_normalized_first_schur_pivot_mod_compact_latex(
    trace: NormalizedFirstSchurPivotDerivation,
    *,
    expand_differences: bool = False,
) -> str:
    """Render the factorized or controlled four-term modulo-compact result."""

    if not isinstance(trace, NormalizedFirstSchurPivotDerivation):
        raise TypeError("trace must be a NormalizedFirstSchurPivotDerivation.")
    if not isinstance(expand_differences, bool):
        raise TypeError("expand_differences must be a bool.")
    right = (
        _render_normalized_endpoint_latex(trace.expanded_mod_compact_result)
        if expand_differences
        else _render_normalized_mod_compact_factorization_latex(
            trace.mod_compact_factorization
        )
    )
    return (
        f"{_render_normalized_endpoint_latex(trace.mod_compact_result.left)} "
        rf"\simeq {right}"
    )


def render_normalized_first_schur_pivot_latex(
    trace: NormalizedFirstSchurPivotDerivation,
) -> RenderedNormalizedFirstSchurPivot:
    """Render every algebraic step of the normalized-pivot trace."""

    if not isinstance(trace, NormalizedFirstSchurPivotDerivation):
        raise TypeError("trace must be a NormalizedFirstSchurPivotDerivation.")
    factorization = trace.mod_compact_factorization
    left_relation, right_relation = trace.offdiagonal_relations
    inverse_rule = trace.inverse_substitution.justification
    if not isinstance(inverse_rule, FormalIdentity):
        raise LatexRenderingError(
            "inverse substitution must retain its supplied formal rule."
        )
    left_model = (
        r"- \left("
        + _render_normalized_endpoint_latex(factorization.left_difference)
        + r"\right)\,"
        + _render_normalized_atom_latex(factorization.left_wiener_hopf)
    )
    right_model = (
        r"\left("
        + _render_normalized_endpoint_latex(factorization.right_difference)
        + r"\right)\,"
        + _render_normalized_atom_latex(factorization.right_wiener_hopf)
    )
    steps = (
        _normalized_relation_step(
            "reduced_definition",
            "Exact first Schur pivot definition",
            trace.reduced_definition.left,
            "=",
            trace.reduced_definition.right,
        ),
        _normalized_relation_step(
            "pivot_definition",
            "Exact normalized pivot definition",
            trace.pivot_definition.left,
            "=",
            trace.pivot_definition.right,
        ),
        _normalized_relation_step(
            "normalized_expansion",
            "Ordered normalization and distribution",
            trace.normalized_expansion.left,
            "=",
            trace.normalized_expansion.right,
        ),
        _normalized_relation_step(
            "inverse_rule",
            "Supplied formal regularizer identity",
            inverse_rule.left,
            r"\overset{\mathrm{formal}}{=}",
            inverse_rule.right,
        ),
        _normalized_relation_step(
            "inverse_substitution",
            "Formal inverse substitution",
            trace.inverse_substitution.left,
            r"\overset{\mathrm{formal}}{=}",
            trace.inverse_substitution.right,
        ),
        _normalized_relation_step(
            "diagonal_recognition",
            "Exact recognition of the diagonal term",
            trace.diagonal_recognition.left,
            "=",
            trace.diagonal_recognition.right,
        ),
        RenderedDerivationStep(
            "exact_result",
            "Final exact formal identity",
            render_normalized_first_schur_pivot_exact_latex(trace),
        ),
        RenderedDerivationStep(
            "offdiagonal_models",
            "Supplied off-diagonal equivalences modulo compact operators",
            r"\begin{aligned}"
            + _render_exact_block_latex(left_relation.exact)
            + rf" &\simeq {left_model} \\ "
            + _render_exact_block_latex(right_relation.exact)
            + rf" &\simeq {right_model}"
            + r"\end{aligned}",
        ),
        RenderedDerivationStep(
            "mod_compact_result",
            "Factorized equivalence modulo compact operators",
            render_normalized_first_schur_pivot_mod_compact_latex(trace),
        ),
        RenderedDerivationStep(
            "controlled_expansion",
            "Controlled expansion of the two differences",
            render_normalized_first_schur_pivot_mod_compact_latex(
                trace,
                expand_differences=True,
            ),
        ),
    )
    return RenderedNormalizedFirstSchurPivot(
        steps,
        tuple(obligation.statement for obligation in trace.proof_obligations),
    )


def render_normalized_factor_classification_markdown(
    trace: NormalizedFirstSchurPivotDerivation,
) -> str:
    """Render the trace's declarative factor classification as Markdown."""

    if not isinstance(trace, NormalizedFirstSchurPivotDerivation):
        raise TypeError("trace must be a NormalizedFirstSchurPivotDerivation.")
    lines = [
        "| Factor | Clase operatorial conocida | Estatus | Fuente matemática |",
        "|---|---|---|---|",
    ]
    lines.extend(
        "| $"
        + _render_normalized_atom_latex(row.factor)
        + "$ | "
        + row.operator_class.value
        + " | "
        + row.status.value
        + " | "
        + row.mathematical_source
        + " |"
        for row in trace.factor_classifications
    )
    return "\n".join(lines)


def render_normalized_proof_obligations_markdown(
    trace: NormalizedFirstSchurPivotDerivation,
) -> str:
    """Render every still-pending analytic proof obligation."""

    if not isinstance(trace, NormalizedFirstSchurPivotDerivation):
        raise TypeError("trace must be a NormalizedFirstSchurPivotDerivation.")
    return "\n".join(
        f"{index}. {obligation.statement}"
        for index, obligation in enumerate(trace.proof_obligations, start=1)
    )


def _render_prefixed_difference_identity_latex(identity: ExactIdentity) -> str:
    left = identity.left
    if not isinstance(left, LinearCombination) or len(left.terms) != 2:
        raise LatexRenderingError(
            "branch normalization left side must have exactly two terms."
        )
    first, second = left.terms
    if (first.coefficient, second.coefficient) != (1, -1):
        raise LatexRenderingError("branch normalization must have signs (+,-).")
    if not first.product.factors or not second.product.factors:
        raise LatexRenderingError("branch normalization terms must be nonempty.")
    prefix = first.product.factors[0]
    if second.product.factors[0] is not prefix:
        raise LatexRenderingError("branch normalization terms must share a prefix.")
    first_tail = Product(first.product.factors[1:])
    second_tail = Product(second.product.factors[1:])
    inner = LinearCombination((Term(1, first_tail), Term(-1, second_tail)))
    return (
        _render_normalized_atom_latex(prefix)
        + r"\,\left("
        + _render_normalized_endpoint_latex(inner)
        + r"\right) = "
        + _render_normalized_endpoint_latex(identity.right)
    )


def _normalized_relation_step(
    key: str,
    title: str,
    left: object,
    relation_symbol: str,
    right: object,
) -> RenderedDerivationStep:
    return RenderedDerivationStep(
        key,
        title,
        f"{_render_normalized_endpoint_latex(left)} {relation_symbol} "
        f"{_render_normalized_endpoint_latex(right)}",
    )


def _render_normalized_endpoint_latex(value: object) -> str:
    if isinstance(value, OperatorAtom):
        return _render_normalized_atom_latex(value)
    if isinstance(value, Product):
        return r"\,".join(
            _render_normalized_atom_latex(atom) for atom in value.factors
        ) or "I"
    if isinstance(value, LinearCombination):
        pieces: list[str] = []
        for index, term in enumerate(value.terms):
            sign, magnitude = _coefficient_sign_and_magnitude(term.coefficient)
            body = r"\,".join(
                _render_normalized_atom_latex(atom)
                for atom in term.product.factors
            ) or "I"
            coefficient = _coefficient_magnitude_latex(magnitude)
            if coefficient != "1":
                body = rf"{coefficient}\,{body}"
            if index == 0:
                pieces.append(("- " if sign == -1 else "") + body)
            else:
                pieces.append((" - " if sign == -1 else " + ") + body)
        return "".join(pieces) if pieces else "0"
    raise LatexRenderingError(
        f"unsupported normalized-pivot endpoint {type(value).__name__}."
    )


def _render_normalized_atom_latex(atom: OperatorAtom) -> str:
    if atom is G1:
        return r"G_{1}"
    if atom is G2:
        return r"G_{2}"
    return render_operator_atom_latex(atom)


def _render_normalized_mod_compact_factorization_latex(
    factorization: NormalizedModCompactFactorization,
) -> str:
    prefix = _render_normalized_endpoint_latex(factorization.prefix)
    left = _render_normalized_endpoint_latex(factorization.left_difference)
    left_wh = _render_normalized_atom_latex(factorization.left_wiener_hopf)
    bridge = _render_normalized_endpoint_latex(factorization.bridge)
    right = _render_normalized_endpoint_latex(factorization.right_difference)
    right_wh = _render_normalized_atom_latex(factorization.right_wiener_hopf)
    suffix = _render_normalized_endpoint_latex(factorization.suffix)
    return (
        r"N_{2} + "
        rf"{prefix}\,\left({left}\right)\,{left_wh}\,{bridge}\,"
        rf"\left({right}\right)\,{right_wh}\,{suffix}"
    )


def render_scalar_latex(
    expression: sp.Expr | KernelAnnotatedExpression,
) -> str:
    """Render a scalar expression and expose retained kernel semantics."""

    if isinstance(expression, KernelAnnotatedExpression):
        scalar = _StructuredScalarLatexPrinter().doprint(expression.expression)
        status_text = ", ".join(
            " ".join(status.value.split("_"))
            for status in expression.semantic_statuses
        )
        evidence_count = sum(
            evidence is not None for evidence in expression.evidences
        )
        return (
            rf"\underbrace{{{scalar}}}_{{\substack{{"
            rf"\text{{kernel status: {status_text}}} \\ "
            rf"\text{{caller hypotheses: {len(expression.hypotheses)}; "
            rf"evidence objects: {evidence_count}}}}}}}"
        )
    if not isinstance(expression, sp.Expr):
        raise TypeError("expression must be a SymPy expression.")
    return _StructuredScalarLatexPrinter().doprint(expression)


def render_kernel_combination_latex(combination: KernelCombination) -> str:
    """Render ordered kernel terms without projecting them to a SymPy Add."""

    if not isinstance(combination, KernelCombination):
        raise TypeError("combination must be a KernelCombination.")
    printer = _StructuredScalarLatexPrinter({"order": "none"})
    pieces: list[str] = []
    for index, term in enumerate(combination.terms):
        sign, magnitude = _coefficient_sign_and_magnitude(term.coefficient)
        scalar = (
            printer.doprint(term.kernel)
            if term.ordered_factors is None
            else " ".join(
                printer.doprint(factor) for factor in term.ordered_factors
            )
        )
        coefficient = _coefficient_magnitude_latex(magnitude)
        body = scalar if coefficient == "1" else rf"{coefficient}\,{scalar}"
        if index == 0:
            pieces.append(rf"- {body}" if sign == -1 else body)
        else:
            pieces.append(rf" - {body}" if sign == -1 else rf" + {body}")
    return "".join(pieces) if pieces else "0"


def render_combined_kernel_action_latex(
    action: sp.Expr | KernelAnnotatedExpression,
    output_variable: sp.Symbol,
    input_function: sp.FunctionClass,
    input_variable: sp.Symbol,
) -> str:
    """Render the selected C22 action and its existing expanded expression."""

    if not isinstance(action, (sp.Expr, KernelAnnotatedExpression)):
        raise TypeError("action must be a SymPy or annotated kernel expression.")
    if not isinstance(output_variable, sp.Symbol) or not isinstance(
        input_variable,
        sp.Symbol,
    ):
        raise TypeError("output_variable and input_variable must be SymPy Symbols.")
    output = _scalar_symbol(output_variable)
    input_ = _scalar_symbol(input_variable)
    operand = render_scalar_latex(input_function(input_variable))
    expanded = render_scalar_latex(action)
    return (
        r"\begin{aligned}"
        rf"\left(C_{{2,2}}^{{(1)}}f\right)\!\left({output}\right)"
        rf" &= \int_0^\infty C_{{2,2}}^{{(1)}}\!\left({output},{input_}\right)"
        rf" {operand}\,d{input_} \\ "
        rf"&= {expanded}"
        r"\end{aligned}"
    )


def render_first_schur_derivation_latex(
    trace: FirstSchurDerivationTrace,
) -> RenderedFirstSchurDerivation:
    """Render the nine ordered presentation steps from one K7A trace."""

    if not isinstance(trace, FirstSchurDerivationTrace):
        raise TypeError("trace must be a FirstSchurDerivationTrace.")
    steps = (
        _render_exact_definition_step(trace),
        _render_offdiagonal_models_step(trace),
        _render_correction_sign_step(trace),
        _render_reduced_model_step(trace),
        _render_correction_expansion_step(trace),
        _render_left_kernel_step(trace),
        _render_right_kernel_step(trace),
        _render_combined_kernel_step(trace),
        _render_compact_action_step(trace),
    )
    return RenderedFirstSchurDerivation(steps)


def _render_exact_definition_step(
    trace: FirstSchurDerivationTrace,
) -> RenderedDerivationStep:
    reduction = trace.reduced_relation.exact
    latex = (
        rf"{_render_reduced_block_latex(reduction)} = "
        rf"{_render_exact_block_latex(reduction.diagonal)}"
        rf" - {_render_exact_block_latex(reduction.left)}\,"
        rf"{render_operator_atom_latex(reduction.regularizer.operator)}\,"
        rf"{_render_exact_block_latex(reduction.right)}."
    )
    return RenderedDerivationStep(
        "exact_definition",
        "Formal algebraic definition",
        latex,
    )


def _render_offdiagonal_models_step(
    trace: FirstSchurDerivationTrace,
) -> RenderedDerivationStep:
    sides = (trace.factorization.left, trace.factorization.right)
    relation_lines: list[str] = []
    for relation, side in zip(trace.offdiagonal_relations, sides, strict=True):
        relation_sign = _factorized_relation_sign(
            relation.model.expression,
            side,
        )
        prefix = "- " if relation_sign < 0 else ""
        relation_lines.append(
            rf"{_render_exact_block_latex(relation.exact)} &\simeq "
            rf"{prefix}{_render_factored_side_latex(side)}"
        )
    relations = r" \\ ".join(relation_lines)
    latex = rf"\begin{{aligned}}{relations}\end{{aligned}}"
    return RenderedDerivationStep(
        "offdiagonal_models",
        "Off-diagonal models",
        latex,
    )


def _render_correction_sign_step(
    trace: FirstSchurDerivationTrace,
) -> RenderedDerivationStep:
    left = _render_factored_side_latex(trace.factorization.left)
    right = _render_factored_side_latex(trace.factorization.right)
    signs = trace.sign_trace
    sign_line = (
        rf"{_unit_sign_latex(signs.schur_sign)}"
        rf"\left[{_unit_sign_latex(signs.leading_a21_sign)}B_{{2,1}}\right]"
        rf"\,R_{{1,1}}\,B_{{1,2}}"
        rf" = {_unit_sign_latex(signs.correction_sign)}"
        rf"B_{{2,1}}\,R_{{1,1}}\,B_{{1,2}}"
    )
    latex = (
        r"\begin{aligned}"
        rf"B_{{2,1}} &:= {left} \\ "
        rf"B_{{1,2}} &:= {right} \\ "
        rf"{sign_line}"
        r"\end{aligned}"
    )
    return RenderedDerivationStep(
        "correction_sign",
        "Correction sign",
        latex,
    )


def _render_reduced_model_step(
    trace: FirstSchurDerivationTrace,
) -> RenderedDerivationStep:
    reduction = trace.reduced_relation.exact
    sign = _binary_sign_latex(trace.sign_trace.correction_sign)
    latex = (
        rf"{_render_reduced_block_latex(reduction)} \simeq "
        rf"{_render_exact_block_latex(reduction.diagonal)} {sign} "
        rf"B_{{2,1}}\,R_{{1,1}}\,B_{{1,2}}."
    )
    return RenderedDerivationStep("reduced_model", "Reduced model", latex)


def _render_correction_expansion_step(
    trace: FirstSchurDerivationTrace,
) -> RenderedDerivationStep:
    lines = _render_operator_term_lines(trace.correction)
    latex = (
        r"\begin{aligned}"
        r"B_{2,1}\,R_{1,1}\,B_{1,2} &={}& "
        rf"{lines[0]} \\ "
        + r" \\ ".join(rf"&&&{line}" for line in lines[1:])
        + r".\end{aligned}"
    )
    return RenderedDerivationStep(
        "correction_expansion",
        "Correction expansion",
        latex,
    )


def _render_left_kernel_step(
    trace: FirstSchurDerivationTrace,
) -> RenderedDerivationStep:
    latex = (
        rf"M_{{2,1}}\!\left({_scalar_symbol(trace.output_variable)},"
        rf"{_scalar_symbol(trace.outer_variable)}\right)"
        rf" = {render_kernel_combination_latex(trace.m21_combination)}."
    )
    return RenderedDerivationStep("left_kernel", "Left kernel", latex)


def _render_right_kernel_step(
    trace: FirstSchurDerivationTrace,
) -> RenderedDerivationStep:
    latex = (
        rf"M_{{1,2}}\!\left({_scalar_symbol(trace.middle_variable)},"
        rf"{_scalar_symbol(trace.input_variable)}\right)"
        rf" = {render_kernel_combination_latex(trace.m12_combination)}."
    )
    return RenderedDerivationStep("right_kernel", "Right kernel", latex)


def _render_combined_kernel_step(
    trace: FirstSchurDerivationTrace,
) -> RenderedDerivationStep:
    output = _scalar_symbol(trace.output_variable)
    input_ = _scalar_symbol(trace.input_variable)
    outer = _scalar_symbol(trace.outer_variable)
    middle = _scalar_symbol(trace.middle_variable)
    left_lines = _render_kernel_term_lines(trace.m21_combination)
    right_lines = _render_kernel_term_lines(trace.m12_combination)
    latex = (
        r"\begin{aligned}"
        rf"C_{{2,2}}^{{(1)}}\!\left({output},{input_}\right)"
        rf" &={{}}& \int_0^\infty M_{{2,1}}\!\left({output},{outer}\right) \\ "
        rf"&&&\quad\times\Biggl[\int_0^\infty R_{{1,1}}\!\left({outer},{middle}\right)"
        rf" M_{{1,2}}\!\left({middle},{input_}\right)\,d{middle}\Biggr]"
        rf"\,d{outer} \\[0.6em] "
        rf"C_{{2,2}}^{{(1)}}\!\left({output},{input_}\right)"
        rf" &={{}}& \int_0^\infty\Bigl[{left_lines[0]} \\ "
        rf"&&&\qquad {left_lines[1]}\Bigr] \\ "
        rf"&&&\quad\times\Biggl\{{\int_0^\infty"
        rf" R_{{1,1}}\!\left({outer},{middle}\right)\Bigl[{right_lines[0]} \\ "
        rf"&&&\qquad\qquad {right_lines[1]}\Bigr]"
        rf"\,d{middle}\Biggr\}}\,d{outer}."
        r"\end{aligned}"
    )
    return RenderedDerivationStep(
        "combined_kernel",
        f"Combined kernel ({_kernel_semantic_summary(trace.combined_kernel)})",
        latex,
    )


def _render_compact_action_step(
    trace: FirstSchurDerivationTrace,
) -> RenderedDerivationStep:
    reduction = trace.compact_action.relation.exact
    output = _scalar_symbol(trace.output_variable)
    input_ = _scalar_symbol(trace.input_variable)
    diagonal_lines = _render_applied_term_lines(trace.compact_action.diagonal)
    diagonal_block = _render_exact_block_latex(reduction.diagonal)
    model_block = _render_reduced_model_block_latex(reduction)
    latex = (
        r"\begin{gathered}\begin{aligned}"
        rf"{_render_reduced_block_latex(reduction)} &\simeq "
        rf"{model_block}, \\ "
        rf"\left({model_block}"
        rf" f\right)\!\left({output}\right)"
        rf" &={{}}& \left({diagonal_block}f\right)\!\left({output}\right) \\ "
        rf"&&&\quad + \int_0^\infty C_{{2,2}}^{{(1)}}"
        rf"\!\left({output},{input_}\right)"
        rf" f\!\left({input_}\right)\,d{input_}."
        r"\end{aligned}\\[0.8em]\begin{aligned}"
        rf"\left({diagonal_block}f\right)\!\left({output}\right)"
        rf" &={{}}& {diagonal_lines[0]} \\ "
        + r" \\ ".join(rf"&&&{line}" for line in diagonal_lines[1:])
        + r".\end{aligned}\end{gathered}"
    )
    return RenderedDerivationStep(
        "compact_model_action",
        "Compact model action declaration "
        f"({_kernel_semantic_summary(trace.compact_action.correction)})",
        latex,
    )


def _kernel_semantic_summary(expression: KernelAnnotatedExpression) -> str:
    statuses = ", ".join(status.value for status in expression.semantic_statuses)
    evidence_count = sum(evidence is not None for evidence in expression.evidences)
    return (
        f"kernel status: {statuses}; caller hypotheses: {len(expression.hypotheses)}; "
        f"evidence objects: {evidence_count}"
    )


def _render_applied_combination_latex(
    combination: AppliedLinearCombination,
) -> str:
    pieces: list[str] = []
    for index, term in enumerate(combination.terms):
        sign, magnitude = _coefficient_sign_and_magnitude(term.coefficient)
        scalar = render_scalar_latex(term.expression)
        coefficient = _coefficient_magnitude_latex(magnitude)
        body = (
            scalar
            if coefficient == "1"
            else rf"{coefficient}\,\left({scalar}\right)"
        )
        if index == 0:
            pieces.append(rf"- {body}" if sign == -1 else body)
        else:
            pieces.append(rf" - {body}" if sign == -1 else rf" + {body}")
    return "".join(pieces)


def _render_operator_term_lines(combination: LinearCombination) -> tuple[str, ...]:
    lines: list[str] = []
    for index, term in enumerate(combination.terms):
        sign, _ = _coefficient_sign_and_magnitude(term.coefficient)
        body = render_term_latex(term, include_sign=False)
        lines.append(_semantic_line_sign(index, sign) + body)
    return tuple(lines)


def _render_kernel_term_lines(combination: KernelCombination) -> tuple[str, ...]:
    printer = _StructuredScalarLatexPrinter({"order": "none"})
    lines: list[str] = []
    for index, term in enumerate(combination.terms):
        sign, magnitude = _coefficient_sign_and_magnitude(term.coefficient)
        factors = term.ordered_factors or (term.kernel,)
        scalar = r"\,".join(printer.doprint(factor) for factor in factors)
        coefficient = _coefficient_magnitude_latex(magnitude)
        body = scalar if coefficient == "1" else rf"{coefficient}\,{scalar}"
        lines.append(_semantic_line_sign(index, sign) + body)
    return tuple(lines)


def _render_applied_term_lines(
    combination: AppliedLinearCombination,
) -> tuple[str, ...]:
    lines: list[str] = []
    for index, term in enumerate(combination.terms):
        sign, magnitude = _coefficient_sign_and_magnitude(term.coefficient)
        scalar = render_scalar_latex(term.expression)
        coefficient = _coefficient_magnitude_latex(magnitude)
        body = scalar if coefficient == "1" else rf"{coefficient}\,\left({scalar}\right)"
        lines.append(_semantic_line_sign(index, sign) + body)
    return tuple(lines)


def _semantic_line_sign(index: int, sign: int | None) -> str:
    if index == 0:
        return "- " if sign == -1 else ""
    return r"{}- " if sign == -1 else r"{}+ "


def _render_exact_block_latex(block: ExactBlock) -> str:
    return rf"A_{{{block.row},{block.column}}}"


def _render_operator_product_factor_latex(atom: OperatorAtom) -> str:
    rendered = render_operator_atom_latex(atom)
    if atom is G1 or atom is G2:
        return rf"\left({rendered}\right)"
    return rendered


def _render_factored_side_latex(side: LinearCombination) -> str:
    if len(side.terms) != 2:
        raise LatexRenderingError("a Schur side must contain exactly two terms.")
    suffix = side.terms[0].product.factors[-1]
    if any(term.product.factors[-1] is not suffix for term in side.terms):
        raise LatexRenderingError("Schur-side terms must share one right factor.")
    inner_terms = LinearCombination(
        tuple(
            Term(
                term.coefficient,
                Product(term.product.factors[:-1]),
            )
            for term in side.terms
        )
    )
    inner = _render_ungrouped_linear_combination_latex(inner_terms)
    return rf"\left({inner}\right)\,{render_operator_atom_latex(suffix)}"


def _render_ungrouped_linear_combination_latex(
    combination: LinearCombination,
) -> str:
    pieces: list[str] = []
    for index, term in enumerate(combination.terms):
        sign, magnitude = _coefficient_sign_and_magnitude(term.coefficient)
        body = r"\,".join(
            render_operator_atom_latex(atom) for atom in term.product.factors
        )
        coefficient = _coefficient_magnitude_latex(magnitude)
        if coefficient != "1":
            body = rf"{coefficient}\,{body}"
        if index == 0:
            pieces.append(rf"- {body}" if sign == -1 else body)
        else:
            pieces.append(rf" - {body}" if sign == -1 else rf" + {body}")
    return "".join(pieces)


def _factorized_relation_sign(
    model: LinearCombination,
    side: LinearCombination,
) -> int:
    if model == side:
        return 1
    if model == -side:
        return -1
    raise LatexRenderingError("off-diagonal model does not match its Schur side.")


def _render_reduced_block_latex(reduction: FirstSchurReduction) -> str:
    return rf"{_render_exact_block_latex(reduction.diagonal)}^{{(1)}}"


def _render_reduced_model_block_latex(reduction: FirstSchurReduction) -> str:
    return rf"{_render_exact_block_latex(reduction.diagonal)}^{{(1),\mathrm{{model}}}}"


def _scalar_symbol(symbol: sp.Symbol) -> str:
    return _StructuredScalarLatexPrinter().doprint(symbol)


def _coefficient_sign_and_magnitude(
    coefficient: Scalar,
) -> tuple[int | None, Scalar]:
    if isinstance(coefficient, complex):
        if coefficient.imag != 0:
            if coefficient.real == 0 and coefficient.imag < 0:
                return -1, -coefficient
            return None, coefficient
        coefficient = _real_complex_component(coefficient.real)
    if coefficient < 0:
        return -1, -coefficient
    return 1, coefficient


def _coefficient_magnitude_latex(coefficient: Scalar) -> str:
    if isinstance(coefficient, Fraction):
        return _StructuredScalarLatexPrinter().doprint(
            sp.Rational(coefficient.numerator, coefficient.denominator)
        )
    if isinstance(coefficient, complex):
        rendered = _StructuredScalarLatexPrinter().doprint(
            _complex_scalar_expression(coefficient)
        )
        if coefficient.real != 0 and coefficient.imag != 0:
            return rf"\left({rendered}\right)"
        return rendered
    return _StructuredScalarLatexPrinter().doprint(sp.sympify(coefficient))


def _real_complex_component(value: float) -> int | float:
    return int(value) if value.is_integer() else value


def _complex_scalar_expression(value: complex) -> sp.Expr:
    real = _float_scalar_expression(value.real)
    imaginary = _float_scalar_expression(value.imag)
    return real + imaginary * sp.I


def _float_scalar_expression(value: float) -> sp.Expr:
    return sp.Integer(int(value)) if value.is_integer() else sp.Float(value)


def _unit_sign_latex(sign: Scalar) -> str:
    if sign == 1:
        return "+"
    if sign == -1:
        return "-"
    return _coefficient_magnitude_latex(sign)


def _binary_sign_latex(sign: Scalar) -> str:
    if sign == 1:
        return "+"
    if sign == -1:
        return "-"
    raise LatexRenderingError("the correction sign must be +1 or -1.")
