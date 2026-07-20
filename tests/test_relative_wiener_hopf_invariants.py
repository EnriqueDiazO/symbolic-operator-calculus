from dataclasses import replace
import inspect

import pytest
import sympy as sp

from symbolic_operator_calculus import (
    DilationMap,
    DilationOperatorModel,
    OrderedRelativeOperatorProduct,
    RelativeWienerHopfIdentity,
    WienerHopfOperatorModel,
    build_relative_wiener_hopf_trace,
    factor_dilation_conjugated_wiener_hopf,
    m2_relative_wiener_hopf_factorizations,
)


def _trace(
    k: int = 1,
    j: int = 2,
    gamma_k: sp.Expr | None = None,
    gamma_j: sp.Expr | None = None,
):
    gamma_k = gamma_k or sp.Symbol(f"gamma_{k}", positive=True)
    gamma_j = gamma_j or sp.Symbol(f"gamma_{j}", positive=True)
    decay = sp.Symbol(f"d_{k}_{j}", positive=True)
    c_k, c_j = (sp.Integer(0), decay) if k < j else (decay, sp.Integer(0))
    return build_relative_wiener_hopf_trace(
        k,
        j,
        gamma_k,
        gamma_j,
        c_k,
        c_j,
    )


@pytest.mark.parametrize("k,j", [(1, 2), (2, 1), (1, 3), (3, 2)])
def test_valid_pairs_are_stored_and_all_domain_objects_are_hashable(k, j):
    trace = _trace(k, j)

    assert (trace.factorization.k, trace.factorization.j) == (k, j)
    assert (trace.relative_dilation.k, trace.relative_dilation.j) == (k, j)
    assert trace == _trace(k, j)
    objects = (
        trace.factorization,
        trace.branch_k,
        trace.branch_j,
        trace.relative_dilation,
        *trace.original_operator.factors,
        trace.original_operator,
        trace.left_operator,
        trace.right_operator,
        trace.identity,
        trace.action,
        trace,
    )
    assert len(set(objects)) <= len(objects)
    assert {item: position for position, item in enumerate(objects)}[trace] == (
        len(objects) - 1
    )


@pytest.mark.parametrize("bad_index", [None, True, False, 1.0, "1"])
def test_branch_map_rejects_non_integer_indices(bad_index):
    with pytest.raises(TypeError):
        DilationMap("branch", bad_index, None, sp.Integer(1))


@pytest.mark.parametrize("bad_index", [0, -1, -5])
def test_branch_map_rejects_non_positive_indices(bad_index):
    with pytest.raises(ValueError):
        DilationMap("branch", bad_index, None, sp.Integer(1))


@pytest.mark.parametrize("bad_index", [None, True, False, 1.0, "1"])
@pytest.mark.parametrize("field", ["k", "j"])
def test_relative_map_rejects_non_integer_indices(field, bad_index):
    relative = _trace().relative_dilation
    with pytest.raises(TypeError):
        replace(relative, **{field: bad_index})


@pytest.mark.parametrize("bad_index", [0, -1, -5])
@pytest.mark.parametrize("field", ["k", "j"])
def test_relative_map_rejects_non_positive_indices(field, bad_index):
    relative = _trace().relative_dilation
    with pytest.raises(ValueError):
        replace(relative, **{field: bad_index})


@pytest.mark.parametrize("pair", [(1, 1), (2, 2)])
def test_relative_map_and_l1_factory_reject_diagonal_pairs(pair):
    relative = _trace().relative_dilation
    with pytest.raises(ValueError):
        replace(relative, k=pair[0], j=pair[1])
    with pytest.raises(ValueError):
        factor_dilation_conjugated_wiener_hopf(
            pair[0],
            pair[1],
            sp.Integer(2),
            sp.Integer(3),
            sp.Integer(0),
            sp.Integer(1),
        )


@pytest.mark.parametrize("bad_index", [None, True, False, 1.0, "1"])
@pytest.mark.parametrize("field", ["k", "j"])
def test_l1_factory_rejects_non_integer_indices(field, bad_index):
    values = {"k": 1, "j": 2, field: bad_index}
    with pytest.raises(TypeError):
        factor_dilation_conjugated_wiener_hopf(
            values["k"],
            values["j"],
            sp.Integer(2),
            sp.Integer(3),
            sp.Integer(0),
            sp.Integer(1),
        )


@pytest.mark.parametrize("bad_index", [0, -1, -5])
@pytest.mark.parametrize("field", ["k", "j"])
def test_l1_factory_rejects_non_positive_indices(field, bad_index):
    values = {"k": 1, "j": 2, field: bad_index}
    with pytest.raises(ValueError):
        factor_dilation_conjugated_wiener_hopf(
            values["k"],
            values["j"],
            sp.Integer(2),
            sp.Integer(3),
            sp.Integer(0),
            sp.Integer(1),
        )


@pytest.mark.parametrize(
    "scale",
    [
        sp.Integer(1),
        sp.Integer(2),
        sp.Rational(3, 2),
        sp.Symbol("gamma", positive=True),
    ],
)
def test_valid_positive_scales_are_accepted(scale):
    branch = DilationMap("branch", 1, None, scale)
    assert branch.scale == scale
    assert hash(branch)


@pytest.mark.parametrize("scale", [0, -1])
def test_scale_requires_a_sympy_expression(scale):
    with pytest.raises(TypeError):
        DilationMap("branch", 1, None, scale)


@pytest.mark.parametrize(
    "scale",
    [
        sp.Integer(0),
        sp.Integer(-1),
        sp.I,
        sp.oo,
        -sp.oo,
        sp.nan,
        sp.Symbol("gamma"),
        sp.Symbol("gamma", real=True),
        sp.Symbol("gamma", negative=True),
        sp.Symbol("gamma", zero=True),
    ],
)
def test_branch_map_rejects_non_positive_real_finite_scales(scale):
    with pytest.raises(ValueError):
        DilationMap("branch", 1, None, scale)


@pytest.mark.parametrize(
    "scale",
    [
        sp.Integer(0),
        sp.Integer(-1),
        sp.I,
        sp.oo,
        -sp.oo,
        sp.nan,
        sp.Symbol("gamma"),
        sp.Symbol("gamma", real=True),
        sp.Symbol("gamma", negative=True),
    ],
)
def test_l1_factory_rejects_non_positive_real_finite_scales(scale):
    with pytest.raises(ValueError):
        factor_dilation_conjugated_wiener_hopf(
            1,
            2,
            sp.Integer(2),
            scale,
            sp.Integer(0),
            sp.Integer(1),
        )


@pytest.mark.parametrize("kind", ["wrong", ""])
def test_dilation_map_rejects_unknown_kinds(kind):
    with pytest.raises(ValueError):
        DilationMap(kind, 1, None, sp.Integer(1))


@pytest.mark.parametrize("kind", [None, 1])
def test_dilation_map_rejects_non_string_kinds(kind):
    with pytest.raises(TypeError):
        DilationMap(kind, 1, None, sp.Integer(1))


def test_branch_and_relative_map_roles_are_coherent():
    with pytest.raises(ValueError):
        DilationMap("branch", 1, 2, sp.Integer(1))

    relative = _trace().relative_dilation
    assert relative.kind == "relative"
    assert relative.j == 2


def test_relative_scale_is_the_exact_ratio_for_numeric_and_symbolic_scales():
    numeric = _trace(
        gamma_k=sp.Integer(2),
        gamma_j=sp.Integer(3),
    ).relative_dilation
    symbolic = _trace().relative_dilation

    assert numeric.scale == sp.Rational(2, 3)
    assert sp.simplify(
        symbolic.scale
        - symbolic.gamma_k / symbolic.gamma_j
    ) == 0
    with pytest.raises(ValueError):
        replace(numeric, scale=sp.Rational(3, 2))


def test_dilation_operator_requires_a_map_and_a_boolean_direction():
    trace = _trace()
    with pytest.raises(TypeError):
        DilationOperatorModel(None)
    with pytest.raises(TypeError):
        replace(trace.original_operator.factors[0], inverse="wrong")


@pytest.mark.parametrize("placement", ["wrong", ""])
def test_wiener_hopf_operator_rejects_unknown_placements(placement):
    operator = _trace().original_operator.factors[1]
    with pytest.raises(ValueError):
        replace(operator, placement=placement)


@pytest.mark.parametrize("placement", [None, 1])
def test_wiener_hopf_operator_rejects_non_string_placements(placement):
    operator = _trace().original_operator.factors[1]
    with pytest.raises(TypeError):
        replace(operator, placement=placement)


def test_wiener_hopf_operator_rejects_a_symbol_from_another_placement():
    trace = _trace()
    original = trace.original_operator.factors[1]
    with pytest.raises(ValueError):
        replace(original, symbol=trace.factorization.left_symbol)
    with pytest.raises(TypeError):
        replace(original, symbol=None)


def test_product_normalizes_iterables_preserves_order_and_is_hashable():
    factors = _trace().left_operator.factors
    product = OrderedRelativeOperatorProduct(factor for factor in factors)

    assert isinstance(product.factors, tuple)
    assert product.factors == factors
    assert hash(product)


@pytest.mark.parametrize("factors", [None, (None, None), (object(), object())])
def test_product_rejects_none_and_foreign_factors(factors):
    with pytest.raises(TypeError):
        OrderedRelativeOperatorProduct(factors)


@pytest.mark.parametrize("factors", [(), (DilationOperatorModel,)])
def test_product_rejects_noncanonical_lengths(factors):
    with pytest.raises((TypeError, ValueError)):
        OrderedRelativeOperatorProduct(factors)


def test_factory_products_have_the_three_canonical_ordered_shapes():
    trace = _trace()
    original = trace.original_operator.factors
    left = trace.left_operator.factors
    right = trace.right_operator.factors

    assert [type(factor) for factor in original] == [
        DilationOperatorModel,
        WienerHopfOperatorModel,
        DilationOperatorModel,
    ]
    assert [type(factor) for factor in left] == [
        DilationOperatorModel,
        WienerHopfOperatorModel,
    ]
    assert [type(factor) for factor in right] == [
        WienerHopfOperatorModel,
        DilationOperatorModel,
    ]
    assert [original[1].placement, left[1].placement, right[0].placement] == [
        "original",
        "left",
        "right",
    ]
    assert original[0].dilation.k == left[0].dilation.k == right[1].dilation.k
    assert original[2].dilation.k == left[0].dilation.j == right[1].dilation.j
    assert not original[0].inverse
    assert original[2].inverse
    assert not left[0].inverse
    assert not right[1].inverse
    assert all(hash(product) for product in (trace.original_operator, trace.left_operator, trace.right_operator))


def test_exact_is_derived_and_cannot_be_supplied_by_callers():
    trace = _trace()
    signature = inspect.signature(RelativeWienerHopfIdentity)

    assert "exact" not in signature.parameters
    assert len(trace.identity.exact_relations) == 2
    with pytest.raises(TypeError):
        RelativeWienerHopfIdentity(
            trace.original_operator,
            trace.left_operator,
            trace.right_operator,
            exact=True,
        )
    with pytest.raises(TypeError):
        RelativeWienerHopfIdentity(
            trace.original_operator,
            trace.left_operator,
            trace.right_operator,
            exact=False,
        )


def test_identity_rejects_swapped_duplicated_and_mixed_products():
    trace12 = _trace(1, 2)
    trace21 = _trace(2, 1)
    invalid = (
        (trace12.original_operator, trace12.right_operator, trace12.left_operator),
        (trace12.original_operator, trace12.left_operator, trace12.left_operator),
        (trace12.original_operator, trace12.left_operator, trace21.right_operator),
    )

    for products in invalid:
        with pytest.raises(ValueError):
            RelativeWienerHopfIdentity(*products)


def test_identity_rejects_wrong_factor_order_and_placement():
    trace = _trace()
    reversed_left = OrderedRelativeOperatorProduct(
        tuple(reversed(trace.left_operator.factors))
    )
    reversed_right = OrderedRelativeOperatorProduct(
        tuple(reversed(trace.right_operator.factors))
    )

    with pytest.raises(ValueError):
        RelativeWienerHopfIdentity(
            trace.original_operator,
            reversed_left,
            trace.right_operator,
        )
    with pytest.raises(ValueError):
        RelativeWienerHopfIdentity(
            trace.original_operator,
            trace.left_operator,
            reversed_right,
        )


def test_identity_rejects_changed_scale_ratio_and_symbol():
    trace = _trace(gamma_k=sp.Integer(2), gamma_j=sp.Integer(3))
    with pytest.raises(ValueError):
        replace(trace.relative_dilation, scale=sp.Rational(3, 2))
    with pytest.raises(ValueError):
        replace(
            trace.left_operator.factors[1],
            symbol=trace.factorization.right_symbol,
        )


def test_identity_rejects_a_factor_from_another_pair():
    trace12 = _trace(1, 2)
    trace13 = _trace(1, 3)
    mixed_original = OrderedRelativeOperatorProduct(
        (
            trace12.original_operator.factors[0],
            trace12.original_operator.factors[1],
            trace13.original_operator.factors[2],
        )
    )

    with pytest.raises(ValueError):
        RelativeWienerHopfIdentity(
            mixed_original,
            trace12.left_operator,
            trace12.right_operator,
        )


def test_l1_metadata_rejects_invalid_indices_scale_ratio_and_symbol():
    factorization = _trace().factorization
    with pytest.raises(ValueError):
        replace(factorization, k=0)
    with pytest.raises(ValueError):
        replace(factorization, gamma_k=sp.Integer(0))
    with pytest.raises(ValueError):
        replace(factorization, relative_scale=sp.Integer(1))
    with pytest.raises(ValueError):
        replace(factorization, left_symbol=factorization.right_symbol)


def test_trace_rejects_foreign_maps_products_actions_and_correspondences():
    trace12 = _trace(1, 2)
    trace21 = _trace(2, 1)

    with pytest.raises(ValueError):
        replace(trace12, branch_k=trace21.branch_k)
    with pytest.raises(ValueError):
        replace(trace12, right_operator=trace21.right_operator)
    with pytest.raises(ValueError):
        replace(trace12, action=trace21.action)
    with pytest.raises(ValueError):
        replace(
            trace12,
            symbol_correspondence_forward=trace21.symbol_correspondence_forward,
        )


def test_m2_factory_preserves_both_directions_and_structural_ratios():
    gamma1, gamma2, decay = sp.symbols(
        "gamma1 gamma2 decay",
        positive=True,
    )
    forward, reverse = m2_relative_wiener_hopf_factorizations(
        gamma1,
        gamma2,
        decay,
    )

    assert (forward.k, forward.j) == (1, 2)
    assert (reverse.k, reverse.j) == (2, 1)
    assert sp.simplify(forward.relative_scale - gamma1 / gamma2) == 0
    assert sp.simplify(reverse.relative_scale - gamma2 / gamma1) == 0
