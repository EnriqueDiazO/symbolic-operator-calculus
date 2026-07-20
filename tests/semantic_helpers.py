"""Explicit caller-supplied semantic objects used by integration tests."""

import sympy as sp

from symbolic_operator_calculus import (
    IntegrationDomain,
    KernelRepresentation,
    KernelRepresentationStatus,
    a11_formal_regularizer,
)


def explicit_r11_kernel_representation() -> KernelRepresentation:
    """Return test evidence for a deliberately formal R11 kernel assumption."""

    output, input_ = sp.symbols("_r11_output _r11_input")
    return KernelRepresentation(
        operator=a11_formal_regularizer(),
        kernel_expression=sp.Function("R11")(output, input_),
        output_variable=output,
        input_variable=input_,
        integration_domain=IntegrationDomain(0, sp.oo, "positive half-line"),
        semantic_status=KernelRepresentationStatus.FORMAL,
        hypotheses=("R11 has the displayed formal kernel representation",),
        evidence=None,
    )


def regularizer_rules():
    """Return default atomic rules with the explicit test representation."""

    from symbolic_operator_calculus import mvp_atomic_rules

    return mvp_atomic_rules(
        regularizer_kernel=explicit_r11_kernel_representation()
    )
