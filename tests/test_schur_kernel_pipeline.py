from dataclasses import replace

import pytest
import sympy as sp

import symbolic_operator_calculus.kernels as kernel_module
from symbolic_operator_calculus import (
    G1,
    G2,
    R11,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    FirstSchurCorrectionFactorization,
    LinearCombination,
    Product,
    SchurFactorizationError,
    Term,
    a22_first_schur_correction,
    apply_combined_kernel_c22,
    apply_linear_combination_ordered,
    c22_integrand,
    combined_kernel_c22,
    explicit_wiener_hopf_rules,
    extract_applied_kernels,
    factor_first_schur_correction,
    kminus_kernel,
    kplus_kernel,
    m12_kernel,
    m21_kernel,
    mvp_atomic_rules,
    positive_decay_symbol,
)


LEFT_FACTORS = (
    (Vtilde_alpha2, Wminus_21),
    (G2, Wminus_21),
)
RIGHT_FACTORS = (
    (Vtilde_alpha1, Wplus_12),
    (G1, Wplus_12),
)


def coefficients(expression):
    return tuple(term.coefficient for term in expression.terms)


def factors(expression):
    return tuple(term.product.factors for term in expression.terms)


def reconstruct(factorization):
    return LinearCombination(
        tuple(
            Term(
                left.coefficient * right.coefficient,
                Product(
                    left.product.factors
                    + (factorization.regularizer,)
                    + right.product.factors
                ),
            )
            for left in factorization.left.terms
            for right in factorization.right.terms
        )
    )


def nested_integrand(kernel):
    outer_factors = sp.Mul.make_args(kernel.function)
    inner = next(factor for factor in outer_factors if isinstance(factor, sp.Integral))
    outer = sp.Mul(*(factor for factor in outer_factors if factor is not inner))
    return outer * inner.function


def test_pipeline_factorization_starts_from_k4b_api(monkeypatch):
    canonical = a22_first_schur_correction()
    calls = []

    def tracked_correction():
        calls.append(True)
        return canonical

    monkeypatch.setattr(kernel_module, "a22_first_schur_correction", tracked_correction)

    result = factor_first_schur_correction()

    assert calls == [True]
    assert result.correction == canonical


def test_factorization_has_exact_left_regularizer_and_right_structure():
    factorization = factor_first_schur_correction()

    assert isinstance(factorization, FirstSchurCorrectionFactorization)
    assert factorization.regularizer is R11
    assert coefficients(factorization.left) == (1, -1)
    assert coefficients(factorization.right) == (1, -1)
    assert factors(factorization.left) == LEFT_FACTORS
    assert factors(factorization.right) == RIGHT_FACTORS
    assert all(
        term.product.factors.count(R11) == 1
        for term in factorization.correction.terms
    )


def test_factorization_reconstructs_the_canonical_correction_exactly():
    factorization = factor_first_schur_correction()

    assert reconstruct(factorization) == factorization.correction
    assert reconstruct(factorization) == a22_first_schur_correction()


@pytest.mark.parametrize(
    "invalid",
    (
        LinearCombination(a22_first_schur_correction().terms[:3]),
        LinearCombination(
            (
                Term(1, Vtilde_alpha2 * Wminus_21 * Vtilde_alpha1 * Wplus_12),
                *a22_first_schur_correction().terms[1:],
            )
        ),
        LinearCombination(
            (
                Term(
                    1,
                    Vtilde_alpha2
                    * Wminus_21
                    * R11
                    * R11
                    * Vtilde_alpha1
                    * Wplus_12,
                ),
                *a22_first_schur_correction().terms[1:],
            )
        ),
        LinearCombination(
            (
                a22_first_schur_correction().terms[1],
                a22_first_schur_correction().terms[0],
                *a22_first_schur_correction().terms[2:],
            )
        ),
        LinearCombination(
            (
                *a22_first_schur_correction().terms[:3],
                Term(-1, a22_first_schur_correction().terms[3].product),
            )
        ),
    ),
)
def test_factorization_rejects_incompatible_canonical_structure(invalid):
    with pytest.raises(SchurFactorizationError):
        factor_first_schur_correction(invalid)


def test_formal_m21_and_m12_are_derived_with_ast_signs():
    x, y, u, v = sp.symbols("x y u v")
    gamma1, gamma2 = sp.symbols("gamma1 gamma2")
    rho1, rho2 = sp.Function("rho1"), sp.Function("rho2")
    g1, g2 = sp.Function("G1"), sp.Function("G2")
    lplus, lminus = sp.Function("Lplus_12"), sp.Function("Lminus_21")

    assert m21_kernel(x, u) == (
        rho2(x) * lminus(gamma2 * x, u) - g2(x) * lminus(x, u)
    )
    assert m12_kernel(v, y) == (
        rho1(v) * lplus(gamma1 * v, y) - g1(v) * lplus(v, y)
    )


def test_m_kernel_signs_follow_factorization_not_a_manual_formula(monkeypatch):
    x, y, u, v = sp.symbols("x y u v")
    original = factor_first_schur_correction()
    reversed_left = LinearCombination(
        tuple(Term(-term.coefficient, term.product) for term in original.left.terms)
    )
    reversed_right = LinearCombination(
        tuple(Term(-term.coefficient, term.product) for term in original.right.terms)
    )

    monkeypatch.setattr(
        kernel_module,
        "factor_first_schur_correction",
        lambda: replace(original, left=reversed_left),
    )
    assert m21_kernel(x, u) == -kernel_module._kernel_expression_from_factor(
        original.left, x, u, None
    )

    monkeypatch.setattr(
        kernel_module,
        "factor_first_schur_correction",
        lambda: replace(original, right=reversed_right),
    )
    assert m12_kernel(v, y) == -kernel_module._kernel_expression_from_factor(
        original.right, v, y, None
    )


def test_formal_c22_uses_derived_m_kernels_in_required_integral_order():
    x, y, u, v = sp.symbols("x y u v")
    result = combined_kernel_c22(
        x,
        y,
        outer_variable=u,
        middle_variable=v,
    )
    inner = next(iter(result.atoms(sp.Integral) - {result}))

    assert result.limits == ((u, 0, sp.oo),)
    assert inner.limits == ((v, 0, sp.oo),)
    assert result.function == m21_kernel(x, u) * inner
    assert inner.function == sp.Function("R11")(u, v) * m12_kernel(v, y)
    assert not result.has(sp.Rational(1, 2), sp.Integer(2), sp.Integer(-2))


def test_explicit_pipeline_uses_normalized_fourier_kernels_without_formal_names():
    x, y, u, v = sp.symbols("x y u v")
    gamma1, gamma2 = sp.symbols("gamma1 gamma2")
    d = positive_decay_symbol()
    chi = sp.Function("chi_infinity")
    rules = explicit_wiener_hopf_rules(decay=d)

    m21 = m21_kernel(x, u, rules=rules)
    m12 = m12_kernel(v, y, rules=rules)
    c22 = combined_kernel_c22(
        x,
        y,
        outer_variable=u,
        middle_variable=v,
        rules=rules,
    )

    assert m21.has(
        chi(gamma2 * x),
        kminus_kernel(gamma2 * x - u, d),
        chi(x),
        kminus_kernel(x - u, d),
    )
    assert m12.has(
        chi(gamma1 * v),
        kplus_kernel(gamma1 * v - y, d),
        chi(v),
        kplus_kernel(v - y, d),
    )
    assert not c22.has(sp.Function("Lplus_12"), sp.Function("Lminus_21"))


def test_automatic_c22_variables_are_fresh_and_distinct_under_name_collisions():
    u, v = sp.symbols("u v")

    output_collision = combined_kernel_c22(u, sp.Symbol("y"))
    input_collision = combined_kernel_c22(sp.Symbol("x"), v)

    assert set(output_collision.variables).isdisjoint({u})
    assert output_collision.variables[0] not in output_collision.free_symbols
    assert input_collision.variables[0] != v
    assert output_collision.variables[0] != next(
        integral.variables[0]
        for integral in output_collision.atoms(sp.Integral)
        if integral is not output_collision
    )


def test_explicit_valid_variables_are_respected_and_collisions_rejected():
    x, y, u, v = sp.symbols("x y u v")

    result = combined_kernel_c22(
        x,
        y,
        outer_variable=u,
        middle_variable=v,
    )

    assert result.variables == [u]
    assert any(integral.variables == [v] for integral in result.atoms(sp.Integral))
    with pytest.raises(ValueError):
        combined_kernel_c22(x, y, outer_variable=x)
    with pytest.raises(ValueError):
        combined_kernel_c22(x, y, outer_variable=u, middle_variable=u)


def test_end_to_end_factorized_c22_expands_to_four_extracted_kernels():
    x, y, u, v = sp.symbols("x y u v")
    f = sp.Function("f")
    correction = a22_first_schur_correction()
    applied = apply_linear_combination_ordered(
        correction,
        f(x),
        x,
        mvp_atomic_rules(),
    )
    extracted = extract_applied_kernels(applied, x, f, y)
    expanded = sp.Add.make_args(sp.expand(c22_integrand(x, u, v, y)))
    expected = tuple(
        term.coefficient * nested_integrand(term.kernel)
        for term in extracted.terms
    )

    assert coefficients(correction) == (1, -1, -1, 1)
    assert len(expanded) == len(expected) == 4
    assert all(
        any(sp.simplify(candidate - term) == 0 for candidate in expanded)
        for term in expected
    )


def test_compact_action_uses_same_api_for_formal_and_explicit_rules():
    x, y, u, v = sp.symbols("x y u v")
    f = sp.Function("f")
    explicit_rules = explicit_wiener_hopf_rules(decay=positive_decay_symbol())

    formal = apply_combined_kernel_c22(
        x,
        f,
        y,
        outer_variable=u,
        middle_variable=v,
    )
    explicit = apply_combined_kernel_c22(
        x,
        f,
        y,
        outer_variable=u,
        middle_variable=v,
        rules=explicit_rules,
    )

    assert formal.limits == explicit.limits == ((y, 0, sp.oo),)
    assert formal.has(sp.Function("Lplus_12"), sp.Function("Lminus_21"))
    assert not explicit.has(sp.Function("Lplus_12"), sp.Function("Lminus_21"))
    assert explicit.has(sp.Function("chi_infinity"))
