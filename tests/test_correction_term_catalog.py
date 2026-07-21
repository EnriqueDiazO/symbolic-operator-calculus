from symbolic_operator_calculus import (
    P_plus,
    R1,
    U1_inverse,
    DerivationRelationKind,
    InterfaceStatus,
    MotifStatus,
    TermOperatorClass,
    build_complete_correction_expansion_trace,
    render_correction_interface_matrix_markdown,
    render_correction_motif_frequency_markdown,
    render_correction_term_catalog_latex,
    render_correction_term_catalog_markdown,
)


def test_catalog_assigns_stable_ids_coefficients_and_ast_terms():
    trace = build_complete_correction_expansion_trace()

    assert tuple(record.identifier for record in trace.terms) == tuple(
        f"C2-T{index:02d}" for index in range(1, 17)
    )
    assert tuple(record.coefficient for record in trace.terms) == (
        1,
        -1,
        -1,
        1,
    ) * 4
    assert tuple(record.term for record in trace.terms) == trace.c3.terms
    assert all(
        record.origin_relation is DerivationRelationKind.EXACT_EQUALITY
        for record in trace.terms
    )


def test_signatures_preserve_concrete_factor_position_and_declared_semantics():
    trace = build_complete_correction_expansion_trace()
    first = trace.terms[0]

    assert tuple(item.factor for item in first.signature) == (
        first.term.product.factors
    )
    assert tuple(item.position for item in first.signature) == tuple(
        range(1, len(first.signature) + 1)
    )
    assert first.signature[0].operator_class is TermOperatorClass.DILATION_OPERATOR
    assert first.signature[1].operator_class is (
        TermOperatorClass.LOCALIZED_WIENER_HOPF_MINUS
    )
    assert next(item for item in first.signature if item.factor is R1).operator_class is (
        TermOperatorClass.MELLIN_PDO_REGULARIZER
    )
    assert next(item for item in first.signature if item.factor is R1).relation_kind is (
        DerivationRelationKind.FORMAL_SUBSTITUTION
    )


def test_every_motif_is_an_exact_slice_and_no_commutation_motif_is_invented():
    trace = build_complete_correction_expansion_trace()

    for record in trace.terms:
        for motif in record.motifs:
            assert record.term.product.factors[
                motif.start_position - 1 : motif.end_position
            ] == motif.factors
    assert all(
        motif.status is not MotifStatus.FORBIDDEN_WITHOUT_COMMUTATION
        for record in trace.terms
        for motif in record.motifs
    )
    assert all(
        any(motif.name == "mellin_regularizer_core" for motif in record.motifs)
        for record in trace.terms
    )


def test_motif_frequency_is_deterministic_and_matches_all_terms():
    trace = build_complete_correction_expansion_trace()
    frequencies = {
        (item.name, item.status): item.count for item in trace.motif_frequency
    }

    assert frequencies == {
        ("left_auxiliary_piece", MotifStatus.EXACT_RECOGNIZED): 16,
        ("mellin_regularizer_core", MotifStatus.UNRESOLVED): 16,
        ("mixed_cauchy_dilation", MotifStatus.EXACT_RECOGNIZED): 16,
        ("normalized_left_shift", MotifStatus.SUPPLIED_MOD_COMPACT): 16,
        ("normalized_right_shift", MotifStatus.SUPPLIED_MOD_COMPACT): 16,
        ("right_auxiliary_piece", MotifStatus.EXACT_RECOGNIZED): 16,
    }


def test_open_interfaces_exclude_only_the_known_exact_auxiliary_subword():
    trace = build_complete_correction_expansion_trace()
    first = trace.terms[0]

    assert all(
        (item.left, item.right) != (U1_inverse, P_plus)
        for item in first.open_interfaces
    )
    assert any(item.right is R1 for item in first.open_interfaces)
    assert all(record.proof_obligations for record in trace.terms)
    assert all(
        obligation.relation_kind
        is DerivationRelationKind.ANALYTIC_PROOF_OBLIGATION
        for record in trace.terms
        for obligation in record.proof_obligations
    )


def test_interface_matrix_has_required_rows_counts_and_conservative_statuses():
    trace = build_complete_correction_expansion_trace()
    rows = {row.interface: row for row in trace.interface_matrix}

    assert rows["multiplicación × Wiener--Hopf"].term_count == 16
    assert rows["dilatación × Wiener--Hopf"].term_count == 16
    assert rows["Wiener--Hopf × dilatación"].term_count == 16
    assert rows["Wiener--Hopf × Cauchy factor"].term_count == 16
    assert rows["Cauchy factor × R_1"].term_count == 16
    assert rows["R_1 × multiplicación"].term_count == 8
    assert rows["R_1 × dilatación"].term_count == 8
    assert rows["Wiener--Hopf × auxiliar"].term_count == 8
    assert rows["Wiener--Hopf × R_1"].status is InterfaceStatus.NOT_APPLICABLE
    assert rows["R_1 × Wiener--Hopf"].status is InterfaceStatus.NOT_APPLICABLE
    assert all(
        row.status
        not in (
            InterfaceStatus.CERTIFIED_EXACT,
            InterfaceStatus.CERTIFIED_MOD_COMPACT,
        )
        for row in trace.interface_matrix
        if row.term_count
    )


def test_catalog_motif_and_interface_renderers_are_stable_and_complete():
    first = build_complete_correction_expansion_trace()
    second = build_complete_correction_expansion_trace()
    markdown = render_correction_term_catalog_markdown(first)
    latex = render_correction_term_catalog_latex(first)

    assert markdown == render_correction_term_catalog_markdown(second)
    assert latex == render_correction_term_catalog_latex(second)
    assert markdown.count("\n") == 17
    assert latex.count("C2-T") == 16
    assert r"\begin{longtable}" in latex
    assert render_correction_motif_frequency_markdown(first).count("\n") == 7
    assert render_correction_interface_matrix_markdown(first).count("\n") == 12
