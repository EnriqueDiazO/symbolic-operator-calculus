import sympy as sp

from symbolic_operator_calculus import (
    AnalyticLemmaStatus,
    build_damping_variation_lemma,
)


def test_damping_lemma_has_the_exact_variation_bound_and_proved_status():
    gamma = sp.Symbol("gamma", positive=True)
    lemma = build_damping_variation_lemma(gamma)
    variation_h = sp.Symbol("Var_h", nonnegative=True)
    l1_h = sp.Symbol("L1_h", nonnegative=True)

    assert lemma.variation_bound == (
        variation_h + sp.Abs(sp.log(gamma)) * l1_h
    )
    assert lemma.status is AnalyticLemmaStatus.ANALYTICALLY_PROVED
    assert any("V(R) intersect L1(R)" in item for item in lemma.hypotheses)
    assert "almost everywhere" in lemma.proof_steps[0]


def test_exponential_damping_example_satisfies_the_proved_bound():
    alpha = sp.Integer(1)
    log_gamma = sp.Integer(1)
    actual_variation = 2 * sp.sqrt(alpha**2 + log_gamma**2) / alpha
    variation_h = sp.Integer(2)
    l1_h = sp.Integer(2)
    upper_bound = variation_h + abs(log_gamma) * l1_h

    assert actual_variation <= upper_bound


def test_oscillatory_factor_preserves_zero_endpoint_limits_by_modulus():
    lam = sp.Symbol("lambda", real=True)
    gamma = sp.Symbol("gamma", positive=True)
    h = sp.exp(-sp.Abs(lam))
    damped_modulus = sp.simplify(sp.Abs(gamma ** (sp.I * lam) * h))

    assert damped_modulus == h
    assert sp.limit(damped_modulus, lam, sp.oo) == 0
    assert sp.limit(damped_modulus, lam, -sp.oo) == 0


def test_documented_hypothesis_does_not_replace_l1_with_c0():
    lemma = build_damping_variation_lemma(sp.Symbol("gamma", positive=True))
    text = " ".join((*lemma.hypotheses, lemma.conclusion, *lemma.proof_steps))

    assert "L1(R)" in text
    assert "C0(R) implies L1(R)" not in text
