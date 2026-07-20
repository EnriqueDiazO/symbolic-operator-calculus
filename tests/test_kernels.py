import sympy as sp

from symbolic_operator_calculus import (
    G1,
    G2,
    R11,
    Vtilde_alpha1,
    Vtilde_alpha2,
    Wminus_21,
    Wplus_12,
    KernelAnnotatedExpression,
    KernelCombination,
    KernelTerm,
    apply_combined_kernel_c22 as _apply_combined_kernel_c22,
    apply_linear_combination_ordered,
    c22_integrand as _c22_integrand,
    collect_bound_symbols,
    combined_kernel_c22 as _combined_kernel_c22,
    expand_ordered,
    extract_applied_kernels,
    extract_integral_kernel,
    m12_kernel,
    m21_kernel,
    main_expression,
)
from semantic_helpers import explicit_r11_kernel_representation, regularizer_rules


def c22_integrand(*args, **kwargs):
    kwargs["regularizer_kernel"] = explicit_r11_kernel_representation()
    return _c22_integrand(*args, **kwargs).expression


def combined_kernel_c22(*args, **kwargs):
    kwargs["regularizer_kernel"] = explicit_r11_kernel_representation()
    return _combined_kernel_c22(*args, **kwargs).expression


def apply_combined_kernel_c22(*args, **kwargs):
    kwargs["regularizer_kernel"] = explicit_r11_kernel_representation()
    return _apply_combined_kernel_c22(*args, **kwargs).expression


REQUIRED_FACTORS = (
    (Vtilde_alpha2, Wminus_21, R11, Vtilde_alpha1, Wplus_12),
    (Vtilde_alpha2, Wminus_21, R11, G1, Wplus_12),
    (G2, Wminus_21, R11, Vtilde_alpha1, Wplus_12),
    (G2, Wminus_21, R11, G1, Wplus_12),
)


def scalar_symbols():
    x, y, u, v = sp.symbols("x y u v")
    return x, y, u, v, sp.Function("f")


def applied_terms():
    x, _, _, _, f = scalar_symbols()
    return apply_linear_combination_ordered(
        expand_ordered(main_expression()),
        f(x),
        x,
        regularizer_rules(),
    )


def kernel_terms():
    x, y, _, _, f = scalar_symbols()
    return extract_applied_kernels(applied_terms(), x, f, y)


def nested_integrand(kernel):
    assert isinstance(kernel, sp.Integral)
    outer_factors = sp.Mul.make_args(kernel.function)
    inner_factors = [
        factor for factor in outer_factors if isinstance(factor, sp.Integral)
    ]
    assert len(inner_factors) == 1
    inner = inner_factors[0]
    assert isinstance(inner, sp.Integral)
    outer_scalar = sp.Mul(
        *(factor for factor in outer_factors if factor is not inner),
    )
    return outer_scalar * inner.function


def nested_integral_parts(kernel):
    assert isinstance(kernel, sp.Integral)
    outer_factors = sp.Mul.make_args(kernel.function)
    inner_factors = [
        factor for factor in outer_factors if isinstance(factor, sp.Integral)
    ]
    assert len(inner_factors) == 1
    inner = inner_factors[0]
    outer_scalar = sp.Mul(
        *(factor for factor in outer_factors if factor is not inner),
    )
    return outer_scalar, inner


def expected_kernel_integrands():
    x, y, u, v, _ = scalar_symbols()
    gamma1, gamma2 = sp.symbols("gamma1 gamma2")
    rho1 = sp.Function("rho1")
    rho2 = sp.Function("rho2")
    g1 = sp.Function("G1")
    g2 = sp.Function("G2")
    lplus = sp.Function("Lplus_12")
    lminus = sp.Function("Lminus_21")
    r11 = sp.Function("R11")

    return (
        rho2(x) * lminus(gamma2 * x, u) * r11(u, v) * rho1(v) * lplus(gamma1 * v, y),
        rho2(x) * lminus(gamma2 * x, u) * r11(u, v) * g1(v) * lplus(v, y),
        g2(x) * lminus(x, u) * r11(u, v) * rho1(v) * lplus(gamma1 * v, y),
        g2(x) * lminus(x, u) * r11(u, v) * g1(v) * lplus(v, y),
    )


def signed_expected_kernel_integrands():
    return tuple(
        coefficient * integrand
        for coefficient, integrand in zip((1, -1, -1, 1), expected_kernel_integrands())
    )


def add_terms(expression):
    return sp.Add.make_args(sp.expand(expression))


def test_i1_extract_integral_kernel_removes_only_input_integral():
    x, y, _, _, f = scalar_symbols()
    applied = applied_terms().result_terms[0].expression

    kernel = extract_integral_kernel(applied, x, f, y)

    assert isinstance(kernel, KernelAnnotatedExpression)
    assert len(applied.atoms(sp.Integral)) == 3
    assert len(kernel.atoms(sp.Integral)) == 2
    assert y not in collect_bound_symbols(kernel.expression)
    assert kernel.kernel_representations == applied.kernel_representations


def test_i2_input_function_disappears_from_extracted_kernel():
    _, y, _, _, f = scalar_symbols()
    kernels = kernel_terms()

    assert all(not term.kernel.has(f(y)) for term in kernels.terms)


def test_i3_input_variable_remains_free_in_extracted_kernel():
    _, y, _, _, _ = scalar_symbols()
    kernels = kernel_terms()

    assert all(y in term.kernel.free_symbols for term in kernels.terms)


def test_i4_outer_and_middle_integrals_remain_bound():
    _, _, u, v, _ = scalar_symbols()
    kernels = kernel_terms()

    assert tuple(collect_bound_symbols(term.kernel) for term in kernels.terms) == (
        {u, v},
        {u, v},
        {u, v},
        {u, v},
    )


def test_i5_k1_contains_rho2_and_rho1():
    x, _, _, v, _ = scalar_symbols()
    k1 = kernel_terms().terms[0].kernel

    assert k1.has(sp.Function("rho2")(x))
    assert k1.has(sp.Function("rho1")(v))


def test_i6_k2_contains_rho2_and_g1():
    x, _, _, v, _ = scalar_symbols()
    k2 = kernel_terms().terms[1].kernel

    assert k2.has(sp.Function("rho2")(x))
    assert k2.has(sp.Function("G1")(v))


def test_i7_k3_contains_g2_and_rho1():
    x, _, _, v, _ = scalar_symbols()
    k3 = kernel_terms().terms[2].kernel

    assert k3.has(sp.Function("G2")(x))
    assert k3.has(sp.Function("rho1")(v))


def test_i8_k4_contains_g2_and_g1():
    x, _, _, v, _ = scalar_symbols()
    k4 = kernel_terms().terms[3].kernel

    assert k4.has(sp.Function("G2")(x))
    assert k4.has(sp.Function("G1")(v))


def test_i9_kernel_coefficients_keep_required_signs():
    kernels = kernel_terms()

    assert tuple(term.coefficient for term in kernels.terms) == (1, -1, -1, 1)


def test_i10_each_kernel_term_keeps_originating_product():
    kernels = kernel_terms()

    assert isinstance(kernels, KernelCombination)
    assert all(isinstance(term, KernelTerm) for term in kernels.terms)
    assert tuple(term.product.factors for term in kernels.terms) == REQUIRED_FACTORS


def test_i11_m12_has_exactly_two_scalar_terms():
    _, y, _, v, _ = scalar_symbols()

    assert len(sp.Add.make_args(m12_kernel(v, y))) == 2


def test_i12_m21_has_exactly_two_scalar_terms():
    x, _, u, _, _ = scalar_symbols()

    assert len(sp.Add.make_args(m21_kernel(x, u))) == 2


def test_i13_expanded_m21_r11_m12_has_four_terms():
    x, y, u, v, _ = scalar_symbols()

    expanded_terms = add_terms(c22_integrand(x, u, v, y))

    assert len(expanded_terms) == 4


def test_i14_expanded_m_products_match_signed_kernel_integrands():
    x, y, u, v, _ = scalar_symbols()
    expanded_terms = add_terms(c22_integrand(x, u, v, y))

    expected_terms = signed_expected_kernel_integrands()

    assert all(
        any(sp.simplify(expanded - expected) == 0 for expanded in expanded_terms)
        for expected in expected_terms
    )


def test_i15_c22_contains_two_nested_integrals():
    x, y, _, _, _ = scalar_symbols()

    c22 = combined_kernel_c22(x, y)

    assert len(c22.atoms(sp.Integral)) == 2


def test_i16_c22_integrand_contains_m21_r11_and_m12_factors():
    x, y, u, v, _ = scalar_symbols()
    c22 = combined_kernel_c22(x, y, outer_variable=u, middle_variable=v)
    outer_scalar, inner = nested_integral_parts(c22)
    inner_factors = set(sp.Mul.make_args(inner.function))

    assert outer_scalar == m21_kernel(x, u)
    assert sp.Function("R11")(u, v) in inner_factors
    assert m12_kernel(v, y) in inner_factors


def test_i17_final_action_adds_one_outer_input_integral():
    x, y, _, _, f = scalar_symbols()

    action = apply_combined_kernel_c22(x, f, y)

    assert isinstance(action, sp.Integral)
    assert tuple(action.variables) == (y,)


def test_i18_final_action_contains_three_total_integrals():
    x, y, _, _, f = scalar_symbols()

    action = apply_combined_kernel_c22(x, f, y)

    assert len(action.atoms(sp.Integral)) == 3


def test_i19_structural_kernel_terms_do_not_depend_on_sympy_mul_order():
    kernels = kernel_terms()

    assert not isinstance(kernels, sp.Add)
    assert tuple(term.product.factors for term in kernels.terms) == REQUIRED_FACTORS
    assert tuple(nested_integrand(term.kernel) for term in kernels.terms) == (
        expected_kernel_integrands()
    )
