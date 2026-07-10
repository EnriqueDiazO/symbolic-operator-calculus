import pytest

from symbolic_operator_calculus import (
    ExactBlock,
    LinearCombination,
    ModCompactRelation,
    OperatorAtom,
    Term,
    WienerHopfModel,
    a12_mod_compact_relation,
    a12_wh_model,
)


def test_exact_block_preserves_structured_indices():
    block = ExactBlock("A", 1, 2)

    assert block.family == "A"
    assert block.row == 1
    assert block.column == 2


def test_exact_block_rejects_zero_or_negative_indices():
    with pytest.raises(ValueError):
        ExactBlock("A", 0, 2)
    with pytest.raises(ValueError):
        ExactBlock("A", 1, -2)


def test_wiener_hopf_model_preserves_row_and_column():
    model = a12_wh_model()

    assert model.family == "A"
    assert model.row == 1
    assert model.column == 2
    assert model.normalized is True


def test_mod_compact_relation_is_distinct_type():
    relation = a12_mod_compact_relation()

    assert isinstance(relation, ModCompactRelation)
    assert not isinstance(relation, ExactBlock)
    assert not isinstance(relation, WienerHopfModel)


def test_mod_compact_relation_keeps_exact_and_model_distinct():
    relation = a12_mod_compact_relation()

    assert relation.exact is not relation.model
    assert relation.exact != relation.model
    assert isinstance(relation.exact, ExactBlock)
    assert isinstance(relation.model, WienerHopfModel)


def test_mod_compact_relation_has_no_automatic_rewrite():
    relation = a12_mod_compact_relation()

    assert not hasattr(relation, "rewrite")
    assert not hasattr(relation, "as_expr")
    assert relation.exact == ExactBlock("A", 1, 2)


def test_mod_compact_relation_rejects_mismatched_block_indices():
    A = OperatorAtom("A")
    expression = LinearCombination((Term(1, A),))
    model = WienerHopfModel("A", 2, 1, expression)

    with pytest.raises(ValueError):
        ModCompactRelation(ExactBlock("A", 1, 2), model)
