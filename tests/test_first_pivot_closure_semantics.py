from symbolic_operator_calculus import (
    ClosureInterfaceStatus,
    ClosureObligationStatus,
    DecisionConfidence,
    MinimalLemmaChoice,
    build_first_pivot_closure_analysis,
    render_closure_evidence_manifest_yaml,
    render_closure_graph_markdown,
    render_closure_interface_matrix_markdown,
    render_closure_obligations_yaml,
    render_minimal_closure_decision_yaml,
)


def test_manifest_has_complete_verified_provenance_without_private_paths():
    analysis = build_first_pivot_closure_analysis()
    rendered = render_closure_evidence_manifest_yaml(analysis.evidence)

    assert len(analysis.evidence) == 21
    assert all(record.status == "verified" for record in analysis.evidence)
    assert all(record.printed_pages for record in analysis.evidence)
    assert all(record.pdf_pages for record in analysis.evidence)
    assert all(len(record.source_checksum) == 64 for record in analysis.evidence)
    assert all(record.hypotheses_used for record in analysis.evidence)
    assert all(record.limitations for record in analysis.evidence)
    assert rendered.count("  - bibkey:") == 21
    assert "/home/" not in rendered


def test_interface_matrix_contains_required_rows_and_no_false_certification():
    matrix = build_first_pivot_closure_analysis().interface_matrix
    rendered = render_closure_interface_matrix_markdown(matrix)

    assert len(matrix) == 9
    assert tuple(row.status for row in matrix) == (
        ClosureInterfaceStatus.NO_RULE,
        ClosureInterfaceStatus.BLOCKED,
        ClosureInterfaceStatus.BLOCKED,
        ClosureInterfaceStatus.SPECIALIZATION_TO_PROVE,
        ClosureInterfaceStatus.ANALOGY_ONLY,
        ClosureInterfaceStatus.ANALOGY_ONLY,
        ClosureInterfaceStatus.BLOCKED,
        ClosureInterfaceStatus.BLOCKED,
        ClosureInterfaceStatus.BLOCKED,
    )
    assert "CERTIFIED_EXACT" not in rendered
    assert "CERTIFIED_MOD_COMPACT" not in rendered


def test_obligation_graph_is_ordered_and_has_no_analytic_proof_claim():
    obligations = build_first_pivot_closure_analysis().obligations
    known: set[str] = set()
    for obligation in obligations:
        assert set(obligation.depends_on) <= known
        known.add(obligation.identifier)

    assert tuple(item.identifier for item in obligations) == tuple(
        f"P-{index:02d}" for index in range(1, 13)
    )
    assert all(
        item.status is not ClosureObligationStatus.ANALYTICALLY_PROVED
        for item in obligations
    )
    assert "ANALYTICALLY_PROVED" not in render_closure_obligations_yaml(obligations)
    assert "P_11 --> P_12" in render_closure_graph_markdown(obligations)


def test_candidates_are_open_and_decision_is_exactly_none_high_confidence():
    analysis = build_first_pivot_closure_analysis()
    decision = analysis.decision

    assert tuple(candidate.identifier for candidate in analysis.candidates) == (
        MinimalLemmaChoice.H1,
        MinimalLemmaChoice.H2,
        MinimalLemmaChoice.H3,
    )
    assert all(candidate.affected_term_count == 16 for candidate in analysis.candidates)
    assert all(
        candidate.status is ClosureObligationStatus.BLOCKED
        for candidate in analysis.candidates
    )
    assert decision.decision is MinimalLemmaChoice.NONE
    assert decision.confidence is DecisionConfidence.HIGH
    assert decision.blocking_obligations == ("P-01", "P-02", "P-04")
    rendered = render_minimal_closure_decision_yaml(decision)
    assert rendered.startswith("decision: NONE\nconfidence: high\n")
    assert "Fredholm" not in decision.prerequisite_statement
