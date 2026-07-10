"""Structured trace of the first reduced Schur derivation for m=2."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import sympy as sp

from .actions import AtomicAction, mvp_atomic_rules
from .blocks import (
    a12_mod_compact_relation,
    a21_mod_compact_relation,
    a22_first_schur_correction,
    a22_first_schur_model,
    a22_first_schur_mod_compact_relation,
    a22_first_schur_reduction,
)
from .kernels import (
    FirstSchurCorrectionFactorization,
    KernelCombination,
    combined_kernel_c22,
    factor_first_schur_correction,
    m12_kernel,
    m12_kernel_combination,
    m21_kernel,
    m21_kernel_combination,
)
from .operators import LinearCombination, OperatorAtom, Scalar, main_expression
from .relations import ModCompactRelation, ModCompactSchurRelation
from .schur import (
    FirstSchurCompactModelAction,
    apply_a22_first_schur_model_compact,
)


class DerivationTraceError(ValueError):
    """Raised when the structured first-Schur trace is inconsistent."""


@dataclass(frozen=True)
class FirstSchurCorrectionSignTrace:
    """Signs that explain the positive first Schur correction."""

    schur_sign: Scalar
    leading_a21_sign: Scalar
    correction_sign: Scalar

    def __post_init__(self) -> None:
        if self.correction_sign != self.schur_sign * self.leading_a21_sign:
            raise DerivationTraceError(
                "correction_sign must be the product of its two source signs."
            )


@dataclass(frozen=True)
class FirstSchurDerivationTrace:
    """Complete structured derivation trace for the first reduced m=2 model."""

    offdiagonal_relations: tuple[ModCompactRelation, ModCompactRelation]
    sign_trace: FirstSchurCorrectionSignTrace
    reduced_relation: ModCompactSchurRelation
    correction: LinearCombination
    factorization: FirstSchurCorrectionFactorization
    m21_combination: KernelCombination
    m12_combination: KernelCombination
    m21: sp.Expr
    m12: sp.Expr
    combined_kernel: sp.Integral
    compact_action: FirstSchurCompactModelAction
    output_variable: sp.Symbol
    input_variable: sp.Symbol
    outer_variable: sp.Symbol
    middle_variable: sp.Symbol

    def __post_init__(self) -> None:
        if len(self.offdiagonal_relations) != 2 or not all(
            isinstance(relation, ModCompactRelation)
            for relation in self.offdiagonal_relations
        ):
            raise DerivationTraceError(
                "offdiagonal_relations must contain the A21 and A12 relations."
            )
        if not isinstance(self.reduced_relation, ModCompactSchurRelation):
            raise DerivationTraceError(
                "reduced_relation must be a ModCompactSchurRelation."
            )
        if self.factorization.correction != self.correction:
            raise DerivationTraceError(
                "factorization must reference the trace correction."
            )
        if self.m21 != self.m21_combination.as_expr():
            raise DerivationTraceError("M21 must project its ordered combination.")
        if self.m12 != self.m12_combination.as_expr():
            raise DerivationTraceError("M12 must project its ordered combination.")
        if self.compact_action.relation != self.reduced_relation:
            raise DerivationTraceError(
                "compact action and trace must share the reduced relation."
            )
        variables = (
            self.output_variable,
            self.input_variable,
            self.outer_variable,
            self.middle_variable,
        )
        if not all(isinstance(variable, sp.Symbol) for variable in variables):
            raise DerivationTraceError("all trace variables must be SymPy Symbols.")
        if len(set(variables)) != 4:
            raise DerivationTraceError("all four trace variables must be distinct.")


def build_first_schur_derivation_trace(
    operand: sp.Expr,
    output_variable: sp.Symbol,
    *,
    input_variable: sp.Symbol,
    outer_variable: sp.Symbol,
    middle_variable: sp.Symbol,
    rules: Mapping[OperatorAtom, AtomicAction] | None = None,
) -> FirstSchurDerivationTrace:
    """Build and verify the structured m=2 derivation using existing APIs."""

    _validate_explicit_variables(
        output_variable,
        input_variable,
        outer_variable,
        middle_variable,
    )
    resolved_rules = mvp_atomic_rules() if rules is None else rules
    left_relation = a21_mod_compact_relation()
    right_relation = a12_mod_compact_relation()
    exact_reduction = a22_first_schur_reduction()
    reduced_relation = a22_first_schur_mod_compact_relation()
    correction = a22_first_schur_correction()
    sign_trace = FirstSchurCorrectionSignTrace(
        schur_sign=exact_reduction.offdiagonal_product_coefficient,
        leading_a21_sign=left_relation.model.expression.terms[0].coefficient,
        correction_sign=(
            exact_reduction.offdiagonal_product_coefficient
            * left_relation.model.expression.terms[0].coefficient
        ),
    )
    factorization = factor_first_schur_correction(correction)
    m21_combination = m21_kernel_combination(
        output_variable,
        outer_variable,
        rules=resolved_rules,
    )
    m12_combination = m12_kernel_combination(
        middle_variable,
        input_variable,
        rules=resolved_rules,
    )
    m21 = m21_combination.as_expr()
    m12 = m12_combination.as_expr()
    combined_kernel = combined_kernel_c22(
        output_variable,
        input_variable,
        outer_variable=outer_variable,
        middle_variable=middle_variable,
        rules=resolved_rules,
    )
    compact_action = apply_a22_first_schur_model_compact(
        operand,
        output_variable,
        input_variable=input_variable,
        outer_variable=outer_variable,
        middle_variable=middle_variable,
        rules=resolved_rules,
    )
    trace = FirstSchurDerivationTrace(
        offdiagonal_relations=(left_relation, right_relation),
        sign_trace=sign_trace,
        reduced_relation=reduced_relation,
        correction=correction,
        factorization=factorization,
        m21_combination=m21_combination,
        m12_combination=m12_combination,
        m21=m21,
        m12=m12,
        combined_kernel=combined_kernel,
        compact_action=compact_action,
        output_variable=output_variable,
        input_variable=input_variable,
        outer_variable=outer_variable,
        middle_variable=middle_variable,
    )
    _validate_derived_trace(
        trace,
        operand,
        resolved_rules,
    )
    return trace


def _validate_explicit_variables(
    output_variable: sp.Symbol,
    input_variable: sp.Symbol,
    outer_variable: sp.Symbol,
    middle_variable: sp.Symbol,
) -> None:
    variables = (
        output_variable,
        input_variable,
        outer_variable,
        middle_variable,
    )
    if not all(isinstance(variable, sp.Symbol) for variable in variables):
        raise DerivationTraceError("all trace variables must be SymPy Symbols.")
    if len(set(variables)) != 4:
        raise DerivationTraceError("all four trace variables must be distinct.")


def _validate_derived_trace(
    trace: FirstSchurDerivationTrace,
    operand: sp.Expr,
    rules: Mapping[OperatorAtom, AtomicAction],
) -> None:
    if trace.reduced_relation.exact != a22_first_schur_reduction():
        raise DerivationTraceError("reduced relation has an unexpected exact side.")
    if trace.reduced_relation.model != a22_first_schur_model():
        raise DerivationTraceError("reduced relation has an unexpected model side.")
    if trace.correction != a22_first_schur_correction():
        raise DerivationTraceError("trace correction differs from the K4B API.")
    if trace.correction != trace.sign_trace.correction_sign * main_expression():
        raise DerivationTraceError("derived sign does not control the correction.")
    if trace.m21 != m21_kernel(
        trace.output_variable,
        trace.outer_variable,
        rules=rules,
    ):
        raise DerivationTraceError("M21 is inconsistent with its variables.")
    if trace.m21_combination != m21_kernel_combination(
        trace.output_variable,
        trace.outer_variable,
        rules=rules,
    ):
        raise DerivationTraceError("ordered M21 combination is inconsistent.")
    if trace.m12 != m12_kernel(
        trace.middle_variable,
        trace.input_variable,
        rules=rules,
    ):
        raise DerivationTraceError("M12 is inconsistent with its variables.")
    if trace.m12_combination != m12_kernel_combination(
        trace.middle_variable,
        trace.input_variable,
        rules=rules,
    ):
        raise DerivationTraceError("ordered M12 combination is inconsistent.")
    if trace.combined_kernel != combined_kernel_c22(
        trace.output_variable,
        trace.input_variable,
        outer_variable=trace.outer_variable,
        middle_variable=trace.middle_variable,
        rules=rules,
    ):
        raise DerivationTraceError("combined kernel is inconsistent with the trace.")
    expected_action = apply_a22_first_schur_model_compact(
        operand,
        trace.output_variable,
        input_variable=trace.input_variable,
        outer_variable=trace.outer_variable,
        middle_variable=trace.middle_variable,
        rules=rules,
    )
    if trace.compact_action != expected_action:
        raise DerivationTraceError("compact action is inconsistent with the trace.")
