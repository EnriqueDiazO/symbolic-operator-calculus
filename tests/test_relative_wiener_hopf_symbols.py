from dataclasses import FrozenInstanceError, replace

import pytest
import sympy as sp

import symbolic_operator_calculus as soc
import symbolic_operator_calculus.relative_wiener_hopf as rwh


@pytest.fixture(scope="module")
def parameters():
    return sp.symbols(
        "gamma1 gamma2 d",
        real=True,
        positive=True,
        finite=True,
    )


@pytest.fixture(scope="module")
def models(parameters):
    return rwh.m2_relative_wiener_hopf_factorizations(*parameters)


@pytest.fixture(scope="module")
def traces(parameters):
    gamma1, gamma2, decay = parameters
    return (
        rwh.build_relative_wiener_hopf_trace(
            1, 2, gamma1, gamma2, sp.Integer(0), decay
        ),
        rwh.build_relative_wiener_hopf_trace(
            2, 1, gamma2, gamma1, decay, sp.Integer(0)
        ),
    )


def _canonical_scaled_symbol(model, frequency, scale):
    original = rwh.substitute_free_variable(
        model.original_symbol,
        model.frequency_variable,
        frequency,
    )
    scaled = rwh.scaled_fourier_symbol(original, frequency, scale)
    indicator = rwh.chi_plus if model.k < model.j else rwh.chi_minus
    return scaled.xreplace(
        {indicator(frequency / scale): indicator(frequency)}
    )


def test_symbol_verification_api_is_publicly_exported():
    assert soc.RelativeSymbolCorrespondenceVerification is (
        rwh.RelativeSymbolCorrespondenceVerification
    )
    assert soc.verify_relative_symbol_correspondences is (
        rwh.verify_relative_symbol_correspondences
    )


@pytest.mark.parametrize("model_index", [0, 1])
def test_both_symbol_placements_are_reconstructed_from_original(models, model_index):
    model = models[model_index]

    verification = rwh.verify_relative_symbol_correspondences(model)

    assert verification.computed_left_symbol == _canonical_scaled_symbol(
        model, model.frequency_variable, model.gamma_j
    )
    assert verification.computed_right_symbol == _canonical_scaled_symbol(
        model, model.frequency_variable, model.gamma_k
    )
    assert sp.simplify(
        verification.computed_left_symbol - verification.stored_left_symbol
    ) == 0
    assert sp.simplify(
        verification.computed_right_symbol - verification.stored_right_symbol
    ) == 0
    assert verification.left_symbol_verified is True
    assert verification.right_symbol_verified is True


@pytest.mark.parametrize("model_index", [0, 1])
def test_both_kernel_placements_are_reconstructed_from_original(models, model_index):
    model = models[model_index]

    verification = rwh.verify_relative_symbol_correspondences(model)

    assert verification.computed_left_kernel == rwh.scaled_convolution_kernel(
        model.original_kernel, model.time_variable, model.gamma_j
    )
    assert verification.computed_right_kernel == rwh.scaled_convolution_kernel(
        model.original_kernel, model.time_variable, model.gamma_k
    )
    assert sp.simplify(
        verification.computed_left_kernel - verification.stored_left_kernel
    ) == 0
    assert sp.simplify(
        verification.computed_right_kernel - verification.stored_right_kernel
    ) == 0
    assert verification.left_kernel_verified is True
    assert verification.right_kernel_verified is True


@pytest.mark.parametrize("model_index", [0, 1])
def test_symbol_and_kernel_routes_share_each_placement_scale(models, model_index):
    model = models[model_index]

    verification = rwh.verify_relative_symbol_correspondences(model)

    assert verification.left_scale == model.gamma_j
    assert verification.right_scale == model.gamma_k
    assert verification.pair == (model.k, model.j)
    assert verification.correspondences_verified is True


@pytest.mark.parametrize("name", ["xi", "eta"])
@pytest.mark.parametrize("model_index", [0, 1])
def test_alternative_real_frequency_is_used_consistently(
    models,
    model_index,
    name,
):
    model = models[model_index]
    frequency = sp.Symbol(name, real=True)

    verification = rwh.verify_relative_symbol_correspondences(
        model,
        frequency_variable=frequency,
    )

    assert verification.frequency_variable == frequency
    assert model.frequency_variable not in verification.computed_left_symbol.free_symbols
    assert model.frequency_variable not in verification.computed_right_symbol.free_symbols
    assert verification.computed_left_symbol == _canonical_scaled_symbol(
        model, frequency, model.gamma_j
    )
    assert verification.computed_right_symbol == _canonical_scaled_symbol(
        model, frequency, model.gamma_k
    )


def test_original_factors_signs_and_halfline_supports_are_preserved(models):
    positive = rwh.verify_relative_symbol_correspondences(models[0])
    negative = rwh.verify_relative_symbol_correspondences(models[1])
    positive_frequency = positive.frequency_variable
    negative_frequency = negative.frequency_variable

    for symbol in (positive.computed_left_symbol, positive.computed_right_symbol):
        assert symbol.has(rwh.chi_plus(positive_frequency))
        assert not symbol.has(rwh.chi_minus(positive_frequency))
        assert symbol.as_coeff_Mul()[0] == 2
    for symbol in (negative.computed_left_symbol, negative.computed_right_symbol):
        assert symbol.has(rwh.chi_minus(negative_frequency))
        assert not symbol.has(rwh.chi_plus(negative_frequency))
        assert symbol.as_coeff_Mul()[0] == -2


def test_identity_input_and_action_evidence_remain_separate(traces):
    for trace in traces:
        symbol_evidence = rwh.verify_relative_symbol_correspondences(trace.identity)
        action_evidence = rwh.verify_relative_product_actions(trace.identity)

        assert len(trace.identity.exact_relations) == 2
        assert action_evidence.actions_verified is True
        assert symbol_evidence.correspondences_verified is True


def test_scaling_functions_are_invoked_independently(monkeypatch, models):
    symbol_calls = []
    kernel_calls = []
    real_symbol_scaling = rwh.scaled_fourier_symbol
    real_kernel_scaling = rwh.scaled_convolution_kernel

    def tracked_symbol_scaling(symbol, frequency, scale):
        symbol_calls.append((symbol, frequency, scale))
        return real_symbol_scaling(symbol, frequency, scale)

    def tracked_kernel_scaling(kernel, time, scale):
        kernel_calls.append((kernel, time, scale))
        return real_kernel_scaling(kernel, time, scale)

    monkeypatch.setattr(rwh, "scaled_fourier_symbol", tracked_symbol_scaling)
    monkeypatch.setattr(rwh, "scaled_convolution_kernel", tracked_kernel_scaling)

    verification = rwh.verify_relative_symbol_correspondences(models[0])

    assert {call[2] for call in symbol_calls} == {models[0].gamma_j, models[0].gamma_k}
    assert {call[2] for call in kernel_calls} == {models[0].gamma_j, models[0].gamma_k}
    assert verification.computed_left_symbol == _canonical_scaled_symbol(
        models[0], models[0].frequency_variable, models[0].gamma_j
    )


def test_invalid_domain_object_and_frequency_are_rejected(models):
    with pytest.raises(TypeError, match="factorization or identity"):
        rwh.verify_relative_symbol_correspondences(object())
    with pytest.raises(TypeError, match="frequency_variable"):
        rwh.verify_relative_symbol_correspondences(
            models[0],
            frequency_variable=sp.Integer(1),
        )
    with pytest.raises(ValueError, match="real SymPy Symbol"):
        rwh.verify_relative_symbol_correspondences(
            models[0],
            frequency_variable=sp.Symbol("z"),
        )


def test_foreign_or_swapped_stored_correspondences_are_rejected(models):
    positive = rwh.verify_relative_symbol_correspondences(models[0])
    negative = rwh.verify_relative_symbol_correspondences(models[1])

    with pytest.raises(ValueError, match="stored_left_symbol"):
        replace(positive, stored_left_symbol=negative.stored_left_symbol)
    with pytest.raises(ValueError, match="stored_left_symbol"):
        replace(positive, stored_left_symbol=positive.stored_right_symbol)
    with pytest.raises(ValueError, match="stored_left_kernel"):
        replace(positive, stored_left_kernel=positive.stored_right_kernel)
    with pytest.raises(ValueError, match="left_scale"):
        replace(positive, left_scale=positive.right_scale)


def test_controlled_incomplete_factorization_is_rejected():
    incomplete = object.__new__(rwh.DilationConjugatedWienerHopfFactorization)

    with pytest.raises(ValueError, match="missing required field"):
        rwh.verify_relative_symbol_correspondences(incomplete)


def test_verification_is_frozen_hashable_and_booleans_are_not_configurable(models):
    verification = rwh.verify_relative_symbol_correspondences(models[0])

    assert hash(verification) == hash(verification)
    with pytest.raises(FrozenInstanceError):
        verification.correspondences_verified = False
    with pytest.raises(TypeError):
        rwh.RelativeSymbolCorrespondenceVerification(
            factorization=models[0],
            pair=(1, 2),
            frequency_variable=models[0].frequency_variable,
            left_scale=models[0].gamma_j,
            right_scale=models[0].gamma_k,
            computed_left_symbol=models[0].left_symbol,
            stored_left_symbol=models[0].left_symbol,
            computed_right_symbol=models[0].right_symbol,
            stored_right_symbol=models[0].right_symbol,
            computed_left_kernel=models[0].left_convolution_kernel,
            stored_left_kernel=models[0].left_convolution_kernel,
            computed_right_kernel=models[0].right_convolution_kernel,
            stored_right_kernel=models[0].right_convolution_kernel,
            correspondences_verified=False,
        )
