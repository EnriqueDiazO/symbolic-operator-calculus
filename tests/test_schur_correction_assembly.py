from __future__ import annotations

from dataclasses import replace

import pytest
import sympy as sp

from p1d_fixtures import build_p1d_case
from symbolic_operator_calculus import (
    AlgebraFactorEvidence,
    AlgebraMembershipNotCertified,
    Branch,
    BranchedDiagonalOperator,
    CutoffAuditRecord,
    CutoffReplacementNotCertified,
    CutoffSide,
    CutoffTreatmentNotCertified,
    CutoffTreatmentStatus,
    ExactIdentity,
    ExactIdentityScope,
    FactorCoverageClassification,
    FactorSymbol,
    FactoredRegularizerStatus,
    FirstSchurPivotCertificate,
    ModCompactEquivalence,
    NoncommutativeProductRuleCertificate,
    NoncommutativeProductRuleMissing,
    OffDiagonalBlockModel,
    OperatorAtom,
    PivotFredholmStatus,
    Product,
    RelativeDilationNotResolved,
    RelativeDilationStatus,
    SchurAlgebraMembershipStatus,
    SchurAssemblyStatus,
    SchurCorrectionSymbolStatus,
    SchurSignConvention,
    article_factor_membership_evidence,
    assemble_first_schur_correction,
    audit_common_algebra_membership,
    audit_relative_dilation_pair,
    audit_schur_cutoffs,
    certify_first_schur_pivot,
    conjugate_cutoff_by_dilation,
    derive_schur_correction_symbol,
    require_certified_schur_symbol,
    require_cutoff_replacement,
)


def _assembly(case=None):
    case = case or build_p1d_case()
    assembly = assemble_first_schur_correction(
        left_factor=case.left,
        diagonal_regularizer=case.central,
        right_factor=case.right,
        source_blocks=(case.left_source, case.right_source),
        assumptions=case.assumptions,
    )
    return case, assembly


def _actual_algebra(case, assembly):
    cutoffs = audit_schur_cutoffs(assembly)
    factors = article_factor_membership_evidence(assembly)
    certificate = audit_common_algebra_membership(
        assembly=assembly,
        factor_evidence=factors,
        assumptions=case.assumptions,
        compact_ideal=case.ideal,
        cutoff_audits=cutoffs,
        evidence=("read-only literature audit",),
    )
    return cutoffs, factors, certificate


def _certified_algebra(case, assembly, *, extension=False):
    algebra = "Synthetic sourced mixed algebra"
    factors = tuple(
        AlgebraFactorEvidence(
            factor=factor,
            factor_label=f"factor-{index}",
            candidate_algebra=algebra,
            classification=FactorCoverageClassification.EXACT_THEOREM,
            source="synthetic complete theorem fixture",
            hypotheses=case.assumptions.assumptions,
            evidence=(f"membership theorem {index}",),
            compact_ideal=case.ideal,
        )
        for index, factor in enumerate(assembly.expanded_expression.factors)
    )
    certificate = audit_common_algebra_membership(
        assembly=assembly,
        factor_evidence=factors,
        assumptions=case.assumptions,
        common_algebra=algebra,
        compact_ideal=case.ideal,
        closure_theorem="ordered-product closure theorem",
        inverse_closedness_theorem="inverse-closedness and Fredholm theorem",
        symbol_homomorphism="ordered multiplicative symbol homomorphism",
        explicit_extension=extension,
        cutoff_audits=audit_schur_cutoffs(assembly),
        evidence=("complete synthetic algebra package",),
    )
    return algebra, factors, certificate


def _no_pair_audit(case):
    return audit_relative_dilation_pair(
        case.u2,
        case.u1,
        source="typed branch-dilation comparison",
        evidence=("gamma1 and gamma2 are independent positive parameters",),
    )


def test_exterior_factors_retain_orientation_weight_differences_and_wh_metadata():
    case = build_p1d_case()

    assert case.left.branch_from == Branch(1)
    assert case.left.branch_to == Branch(2)
    assert case.right.branch_from == Branch(2)
    assert case.right.branch_to == Branch(1)
    assert case.left.domain == case.left.codomain == case.space
    assert case.left.weight == case.space.weight_exponent
    assert [term.coefficient for term in case.left.difference.expanded_expression.terms] == [
        1,
        -1,
    ]
    assert len(case.left.localizers) == 2
    assert len(case.right.localizers) == 2
    assert case.left.wiener_hopf_symbols
    assert case.right.wiener_hopf_symbols
    assert case.left.relative_dilations == (case.u2,)
    assert case.right.relative_dilations == (case.u1,)
    assert len(case.left.multipliers) == 2


def test_exact_assembly_matches_grouped_source_and_seven_factor_words():
    case, assembly = _assembly()
    candidate = case.central.candidate_regularizer

    assert assembly.assembly_status is (
        SchurAssemblyStatus.CERTIFIED_EXACT_SCHUR_ASSEMBLY
    )
    assert case.central.status is (
        FactoredRegularizerStatus.CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER
    )
    assert assembly.grouped_expression == Product(
        (case.left.atom, candidate.atom, case.right.atom)
    )
    assert assembly.source_expression == Product(
        (
            case.left_source.atom,
            candidate.mellin_regularizer.atom,
            case.right_source.atom,
        )
    )
    assert assembly.expanded_expression == Product(
        (
            case.left.difference.atom,
            case.left.localized_wiener_hopf.atom,
            candidate.left_interface.atom,
            candidate.mellin_regularizer.atom,
            candidate.right_interface_inverse.atom,
            case.right.difference.atom,
            case.right.localized_wiener_hopf.atom,
        )
    )
    assert isinstance(assembly.assembly_relation, ExactIdentity)
    assert isinstance(assembly.source_relation, ExactIdentity)
    assert isinstance(assembly.correction_definition, ExactIdentity)
    assert assembly.correction_definition.left == assembly.correction_atom
    assert assembly.correction_definition.right == assembly.grouped_expression
    assert assembly.word_ledger.grouped_regularizer_word == assembly.grouped_expression
    assert assembly.word_ledger.expanded_word == assembly.expanded_expression
    assert assembly.word_ledger.raw_schur_word is None
    assert assembly.domain.branch == assembly.codomain.branch == Branch(2)


def test_grouped_word_keeps_the_p1c_regularizer_opaque_and_does_not_distribute():
    case, assembly = _assembly()

    assert assembly.grouped_expression.factors[1] == (
        case.central.candidate_regularizer.atom
    )
    assert len(assembly.grouped_expression.factors) == 3
    assert len(assembly.expanded_expression.factors) == 7
    assert case.left.difference.expanded_expression not in (
        assembly.expanded_expression.factors
    )
    assert case.right.difference.expanded_expression not in (
        assembly.expanded_expression.factors
    )


def test_missing_source_factor_fails_closed_without_a_relation():
    case = build_p1d_case()
    assembly = assemble_first_schur_correction(
        left_factor=case.left,
        diagonal_regularizer=case.central,
        right_factor=case.right,
        source_blocks=(case.left_source,),
        assumptions=case.assumptions,
    )

    assert assembly.assembly_status is SchurAssemblyStatus.BLOCKED_BY_MISSING_FACTOR
    assert assembly.assembly_relation is None
    assert assembly.source_relation is None
    assert assembly.correction_definition is None


def test_reversed_source_blocks_are_not_recovered_from_their_names():
    case = build_p1d_case()
    assembly = assemble_first_schur_correction(
        left_factor=case.left,
        diagonal_regularizer=case.central,
        right_factor=case.right,
        source_blocks=(case.right_source, case.left_source),
        assumptions=case.assumptions,
    )

    assert assembly.assembly_status is (
        SchurAssemblyStatus.BLOCKED_BY_SOURCE_ORDER_MISMATCH
    )
    assert assembly.assembly_relation is None


def test_wrong_branch_orientation_blocks_before_source_matching():
    case = build_p1d_case()
    wrong_right = replace(
        case.right,
        branch_from=Branch(1),
        branch_to=Branch(2),
    )
    assembly = assemble_first_schur_correction(
        left_factor=case.left,
        diagonal_regularizer=case.central,
        right_factor=wrong_right,
        source_blocks=(case.left_source, case.right_source),
        assumptions=case.assumptions,
    )

    assert assembly.assembly_status is (
        SchurAssemblyStatus.BLOCKED_BY_DOMAIN_OR_CODOMAIN
    )
    assert assembly.assembly_relation is None


def test_assembly_status_vocabulary_is_exactly_the_seven_allowed_values():
    assert {item.value for item in SchurAssemblyStatus} == {
        "CERTIFIED_EXACT_SCHUR_ASSEMBLY",
        "CERTIFIED_MOD_COMPACTS_SCHUR_ASSEMBLY",
        "CONDITIONAL_SCHUR_ASSEMBLY",
        "BLOCKED_BY_SOURCE_ORDER_MISMATCH",
        "BLOCKED_BY_DOMAIN_OR_CODOMAIN",
        "BLOCKED_BY_SIGN_CONVENTION",
        "BLOCKED_BY_MISSING_FACTOR",
    }


def test_cutoff_audit_enumerates_four_retained_localizers_in_order():
    _, assembly = _assembly()
    records = audit_schur_cutoffs(assembly)

    assert len(records) == 4
    assert [record.side for record in records] == [
        CutoffSide.LEFT,
        CutoffSide.RIGHT,
        CutoffSide.LEFT,
        CutoffSide.RIGHT,
    ]
    assert [record.branch for record in records] == [
        Branch(2),
        Branch(1),
        Branch(1),
        Branch(2),
    ]
    assert all(
        record.status is CutoffTreatmentStatus.NO_REPLACEMENT_NEEDED
        for record in records
    )
    assert all(record.scale == 1 for record in records)


def test_exact_cutoff_covariance_does_not_certify_replacement():
    case, assembly = _assembly()
    record = audit_schur_cutoffs(assembly)[2]
    trace = conjugate_cutoff_by_dilation(
        case.u1,
        record.cutoff,
        case.u1_inverse,
    )
    covariance_only = CutoffAuditRecord(
        cutoff=record.cutoff,
        cutoff_function=record.cutoff_function,
        argument=record.argument,
        scale=record.scale,
        branch=record.branch,
        side=record.side,
        adjacent_operator=record.adjacent_operator,
        support_claim=record.support_claim,
        source=record.source,
        status=CutoffTreatmentStatus.EXACT_COVARIANCE_ONLY,
        replacement_cutoff=trace.scaled_cutoff,
        covariance_trace=trace,
        evidence=(trace.exact_identity,),
    )

    assert sp.simplify(trace.scaled_cutoff.radial_scale - case.gamma1) == 0
    with pytest.raises(CutoffReplacementNotCertified, match="not a replacement"):
        require_cutoff_replacement(covariance_only)


def test_a_separate_mod_compact_cutoff_theorem_is_required_and_returned():
    case, assembly = _assembly()
    record = audit_schur_cutoffs(assembly)[0]
    replacement = replace(
        record.cutoff,
        atom=OperatorAtom("M_scaled_chi"),
        radial_scale=case.gamma2,
    )
    relation = ModCompactEquivalence(
        record.cutoff.atom,
        replacement.atom,
        space=case.space,
        compact_ideal=case.ideal,
        evidence=("specific replacement theorem",),
    )
    certified = replace(
        record,
        status=CutoffTreatmentStatus.MOD_COMPACTS_REPLACEMENT_CERTIFIED,
        replacement_cutoff=replacement,
        semantic_relation=relation,
        evidence=("specific replacement theorem",),
    )

    assert require_cutoff_replacement(certified) is relation


def test_relative_dilation_audit_cancels_only_an_adjacent_inverse_pair():
    case = build_p1d_case()
    audit = audit_relative_dilation_pair(
        case.u1_inverse,
        case.u1,
        source="exact inverse dilation theorem",
        evidence=("typed inverse pair",),
    )

    assert audit.are_inverse
    assert audit.adjacent
    assert audit.status is RelativeDilationStatus.EXACT_ADJACENT_CANCELLATION
    assert isinstance(audit.semantic_relation, ExactIdentity)


def test_actual_inverse_pair_is_blocked_by_the_ordered_intervening_factors():
    case = build_p1d_case()
    candidate = case.central.candidate_regularizer
    between = (
        OperatorAtom("P_plus"),
        candidate.mellin_regularizer.atom,
        candidate.right_interface_inverse.atom,
        case.right.difference.coefficient.atom,
    )
    audit = audit_relative_dilation_pair(
        case.u1_inverse,
        case.u1,
        intervening_factors=between,
        source="article seven-factor word audit",
        evidence=("no commutation or distribution theorem",),
    )

    assert audit.are_inverse
    assert not audit.adjacent
    assert audit.intervening_factors == between
    assert audit.status is RelativeDilationStatus.BLOCKED_BY_INTERVENING_FACTORS
    assert audit.semantic_relation is None


def test_noninverse_pair_and_exact_covariance_have_distinct_statuses():
    case = build_p1d_case()
    no_pair = _no_pair_audit(case)
    covariance = ExactIdentity(
        Product((case.u1.atom, OperatorAtom("M"), case.u1_inverse.atom)),
        Product((OperatorAtom("M_scaled"),)),
        evidence=("exact covariance theorem",),
        scope=ExactIdentityScope.OPERATORIAL,
    )
    reduced = audit_relative_dilation_pair(
        case.u1,
        case.u1_inverse,
        intervening_factors=(OperatorAtom("M"),),
        covariance_relation=covariance,
        source="exact multiplication covariance",
        evidence=(covariance,),
    )

    assert no_pair.status is RelativeDilationStatus.NO_RELATIVE_DILATION_PAIR
    assert reduced.status is RelativeDilationStatus.EXACT_COVARIANCE_REDUCTION
    assert reduced.semantic_relation is covariance


def test_algebra_audit_records_the_real_common_source_gap_factor_by_factor():
    case, assembly = _assembly()
    _, factors, algebra = _actual_algebra(case, assembly)

    assert len(factors) == 7
    assert factors[1].classification is FactorCoverageClassification.EXACT_THEOREM
    assert factors[3].candidate_algebra == "KKL-2014 radial Mellin PDO algebra"
    assert factors[5].classification is FactorCoverageClassification.SOURCE_GAP
    assert algebra.status is (
        SchurAlgebraMembershipStatus.BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP
    )
    assert algebra.missing_generator
    assert algebra.missing_closure_property
    assert algebra.missing_cutoff_theorem
    assert algebra.missing_dilation_membership
    assert algebra.missing_inverse_closedness
    assert algebra.missing_symbol_homomorphism
    assert not algebra.is_certified


def test_common_algebra_cannot_be_certified_by_mixing_component_algebras():
    case, assembly = _assembly()
    factors = article_factor_membership_evidence(assembly)
    algebra = audit_common_algebra_membership(
        assembly=assembly,
        factor_evidence=factors,
        assumptions=case.assumptions,
        common_algebra="Karlovich-2025 half-line algebra",
        compact_ideal=case.ideal,
        closure_theorem="component closure only",
        inverse_closedness_theorem="component inverse theorem only",
        symbol_homomorphism="component symbol only",
        evidence=("incomplete mixed claim",),
    )

    assert algebra.status is (
        SchurAlgebraMembershipStatus.BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP
    )


@pytest.mark.parametrize("extension", [False, True])
def test_synthetic_complete_source_package_certifies_common_or_explicit_extension(
    extension,
):
    case, assembly = _assembly()
    _, factors, algebra = _certified_algebra(case, assembly, extension=extension)

    expected = (
        SchurAlgebraMembershipStatus.CERTIFIED_MEMBERSHIP_IN_EXPLICIT_EXTENSION
        if extension
        else SchurAlgebraMembershipStatus.CERTIFIED_COMMON_ALGEBRA_MEMBERSHIP
    )
    assert algebra.status is expected
    assert algebra.is_certified
    assert tuple(item.factor for item in factors) == assembly.expanded_expression.factors


def test_failed_assembly_prevents_algebra_attempt():
    case = build_p1d_case()
    blocked = assemble_first_schur_correction(
        left_factor=case.left,
        diagonal_regularizer=case.central,
        right_factor=case.right,
        source_blocks=(case.left_source,),
        assumptions=case.assumptions,
    )
    algebra = audit_common_algebra_membership(
        assembly=blocked,
        factor_evidence=(),
        assumptions=case.assumptions,
    )

    assert algebra.status is (
        SchurAlgebraMembershipStatus.NOT_ATTEMPTED_BECAUSE_ASSEMBLY_FAILED
    )


def test_algebra_status_vocabulary_is_exactly_the_nine_allowed_values():
    assert {item.value for item in SchurAlgebraMembershipStatus} == {
        "CERTIFIED_COMMON_ALGEBRA_MEMBERSHIP",
        "CERTIFIED_MEMBERSHIP_IN_EXPLICIT_EXTENSION",
        "CONDITIONAL_ALGEBRA_MEMBERSHIP",
        "BLOCKED_BY_CUTOFF_REPLACEMENT",
        "BLOCKED_BY_RELATIVE_DILATION",
        "BLOCKED_BY_OFF_DIAGONAL_WH_MEMBERSHIP",
        "BLOCKED_BY_CENTRAL_REGULARIZER_MEMBERSHIP",
        "BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP",
        "NOT_ATTEMPTED_BECAUSE_ASSEMBLY_FAILED",
    }


def test_actual_symbol_gate_stops_at_algebra_membership_without_a_formula():
    case, assembly = _assembly()
    cutoffs, _, algebra = _actual_algebra(case, assembly)
    result = derive_schur_correction_symbol(
        assembly=assembly,
        algebra_certificate=algebra,
        cutoff_audits=cutoffs,
    )

    assert result.symbol_status is (
        SchurCorrectionSymbolStatus.BLOCKED_BY_ALGEBRA_MEMBERSHIP
    )
    assert result.symbol is None
    with pytest.raises(AlgebraMembershipNotCertified):
        require_certified_schur_symbol(result)


def test_certified_algebra_still_requires_all_four_cutoff_records():
    case, assembly = _assembly()
    _, _, algebra = _certified_algebra(case, assembly)
    result = derive_schur_correction_symbol(
        assembly=assembly,
        algebra_certificate=algebra,
        cutoff_audits=(),
        dilation_audits=(_no_pair_audit(case),),
    )

    assert result.symbol_status is SchurCorrectionSymbolStatus.BLOCKED_BY_CUTOFFS
    with pytest.raises(CutoffTreatmentNotCertified):
        require_certified_schur_symbol(result)


def test_symbol_gate_rejects_an_unresolved_relative_dilation():
    case, assembly = _assembly()
    algebra_name, _, algebra = _certified_algebra(case, assembly)
    unresolved = audit_relative_dilation_pair(
        case.u1_inverse,
        case.u1,
        intervening_factors=(OperatorAtom("R1"),),
        source="ordered word",
        evidence=("intervening regularizer",),
    )
    rule = NoncommutativeProductRuleCertificate(
        algebra=algebra_name,
        ordered=True,
        multiplicative=True,
        source="synthetic rule",
        evidence=("complete product theorem",),
    )
    result = derive_schur_correction_symbol(
        assembly=assembly,
        algebra_certificate=algebra,
        product_rule=rule,
        cutoff_audits=audit_schur_cutoffs(assembly),
        dilation_audits=(unresolved,),
    )

    assert result.symbol_status is (
        SchurCorrectionSymbolStatus.BLOCKED_BY_NONCOMMUTATIVE_PRODUCT
    )
    with pytest.raises(RelativeDilationNotResolved):
        require_certified_schur_symbol(result)


def test_symbol_gate_requires_an_ordered_noncommutative_product_rule():
    case, assembly = _assembly()
    _, _, algebra = _certified_algebra(case, assembly)
    result = derive_schur_correction_symbol(
        assembly=assembly,
        algebra_certificate=algebra,
        cutoff_audits=audit_schur_cutoffs(assembly),
        dilation_audits=(_no_pair_audit(case),),
    )

    assert result.symbol_status is (
        SchurCorrectionSymbolStatus.BLOCKED_BY_NONCOMMUTATIVE_PRODUCT
    )
    with pytest.raises(NoncommutativeProductRuleMissing):
        require_certified_schur_symbol(result)


def test_fully_certified_synthetic_symbol_preserves_factor_order():
    case, assembly = _assembly()
    algebra_name, _, algebra = _certified_algebra(case, assembly)
    factor_symbols = tuple(
        FactorSymbol(
            factor=factor,
            symbol=("sigma", index, factor.name),
            algebra=algebra_name,
            source="synthetic factor-symbol theorem",
            evidence=(f"symbol theorem {index}",),
        )
        for index, factor in enumerate(assembly.expanded_expression.factors)
    )
    rule = NoncommutativeProductRuleCertificate(
        algebra=algebra_name,
        ordered=True,
        multiplicative=True,
        source="synthetic ordered product theorem",
        evidence=("multiplicative homomorphism",),
    )
    result = derive_schur_correction_symbol(
        assembly=assembly,
        algebra_certificate=algebra,
        factor_symbols=factor_symbols,
        product_rule=rule,
        cutoff_audits=audit_schur_cutoffs(assembly),
        dilation_audits=(_no_pair_audit(case),),
    )
    symbol = require_certified_schur_symbol(result)

    assert result.symbol_status is (
        SchurCorrectionSymbolStatus.CERTIFIED_SCHUR_CORRECTION_SYMBOL
    )
    assert symbol.ordered_factor_symbols == factor_symbols
    assert symbol.expression == tuple(item.symbol for item in factor_symbols)
    assert not isinstance(symbol.expression, sp.Expr)


def _pivot_data(case, assembly):
    diagonal = BranchedDiagonalOperator(
        atom=OperatorAtom("A22"),
        branch=Branch(2),
        domain=case.space,
        codomain=case.space,
        source="article section 06:27-42",
        evidence=("A22 exact definition",),
    )
    left_block = OffDiagonalBlockModel(
        block_atom=OperatorAtom("A21"),
        branch_from=Branch(1),
        branch_to=Branch(2),
        domain=case.space,
        codomain=case.space,
        normalized_exterior=case.left,
        model_sign=-1,
        compact_ideal=case.ideal,
        source="article section 06:86-95",
        evidence=("A21 is minus L21 modulo compact",),
    )
    right_block = OffDiagonalBlockModel(
        block_atom=OperatorAtom("A12"),
        branch_from=Branch(2),
        branch_to=Branch(1),
        domain=case.space,
        codomain=case.space,
        normalized_exterior=case.right,
        model_sign=1,
        compact_ideal=case.ideal,
        source="article section 06:78-85",
        evidence=("A12 is plus R12 modulo compact",),
    )
    sign = SchurSignConvention(
        raw_schur_sign=-1,
        left_model_sign=-1,
        right_model_sign=1,
        correction_sign_in_model_pivot=1,
        source="article sections 06:281-354",
        evidence=("external Schur minus and normalized A21 minus",),
    )
    return diagonal, left_block, right_block, sign


def test_first_pivot_keeps_raw_minus_and_positive_model_correction_separate():
    case, assembly = _assembly()
    diagonal, left_block, right_block, sign = _pivot_data(case, assembly)
    pivot = certify_first_schur_pivot(
        diagonal_block=diagonal,
        correction=assembly,
        left_block=left_block,
        right_block=right_block,
        sign=sign,
    )

    assert isinstance(pivot, FirstSchurPivotCertificate)
    assert pivot.raw_schur_word == Product(
        (
            left_block.block_atom,
            case.central.candidate_regularizer.atom,
            right_block.block_atom,
        )
    )
    assert [term.coefficient for term in pivot.source_expression.terms] == [1, -1]
    assert [term.coefficient for term in pivot.pivot_expression.terms] == [1, 1]
    assert isinstance(pivot.pivot_relation, ModCompactEquivalence)
    assert pivot.pivot_relation.compact_ideal == case.ideal
    assert pivot.fredholm_status is PivotFredholmStatus.NOT_CERTIFIED
    assert pivot.symbol_status is SchurCorrectionSymbolStatus.NOT_ATTEMPTED
    assert "\\simeq" in pivot.latex


def test_wrong_sign_convention_blocks_the_model_pivot_and_fredholm_promotion():
    case, assembly = _assembly()
    diagonal, left_block, right_block, _ = _pivot_data(case, assembly)
    wrong = SchurSignConvention(
        raw_schur_sign=-1,
        left_model_sign=1,
        right_model_sign=1,
        correction_sign_in_model_pivot=-1,
        source="intentionally wrong normalized sign",
        evidence=("negative test",),
    )
    pivot = certify_first_schur_pivot(
        diagonal_block=diagonal,
        correction=assembly,
        left_block=left_block,
        right_block=right_block,
        sign=wrong,
    )

    assert pivot.assembly_status is SchurAssemblyStatus.BLOCKED_BY_SIGN_CONVENTION
    assert pivot.pivot_expression is None
    assert pivot.pivot_relation is None
    assert pivot.fredholm_status is PivotFredholmStatus.NOT_CERTIFIED


def test_symbol_status_vocabulary_is_exactly_the_six_allowed_values():
    assert {item.value for item in SchurCorrectionSymbolStatus} == {
        "CERTIFIED_SCHUR_CORRECTION_SYMBOL",
        "CONDITIONAL_SCHUR_CORRECTION_SYMBOL",
        "BLOCKED_BY_ALGEBRA_MEMBERSHIP",
        "BLOCKED_BY_CUTOFFS",
        "BLOCKED_BY_NONCOMMUTATIVE_PRODUCT",
        "NOT_ATTEMPTED",
    }


def test_all_eight_required_open_obligations_remain_visible():
    _, assembly = _assembly()

    assert {item.key for item in assembly.proof_obligations} == {
        "full_regularizer_algebra_membership",
        "cutoff_replacement_mod_compacts",
        "off_diagonal_wh_membership",
        "relative_dilation_membership",
        "wh_mellin_wh_sandwich_membership",
        "fredholm_algebra_membership",
        "schur_correction_symbol",
        "first_schur_pivot_fredholmness",
    }
