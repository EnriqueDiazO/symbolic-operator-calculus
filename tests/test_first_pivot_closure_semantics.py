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


def test_interface_matrix_contains_two_phase_q_certifications_only():
    matrix = build_first_pivot_closure_analysis().interface_matrix
    rendered = render_closure_interface_matrix_markdown(matrix)

    assert len(matrix) == 9
    assert tuple(row.status for row in matrix) == (
        ClosureInterfaceStatus.CERTIFIED_MOD_COMPACT,
        ClosureInterfaceStatus.CERTIFIED_MOD_COMPACT,
        ClosureInterfaceStatus.CERTIFIED_EXACT,
        ClosureInterfaceStatus.CERTIFIED_MOD_COMPACT,
        ClosureInterfaceStatus.ANALOGY_ONLY,
        ClosureInterfaceStatus.ANALOGY_ONLY,
        ClosureInterfaceStatus.BLOCKED,
        ClosureInterfaceStatus.BLOCKED,
        ClosureInterfaceStatus.BLOCKED,
    )
    assert rendered.count("CERTIFIED_EXACT") == 1
    assert rendered.count("CERTIFIED_MOD_COMPACT") == 3


def test_obligation_graph_marks_only_the_two_phase_q_nodes_proved():
    obligations = build_first_pivot_closure_analysis().obligations
    known: set[str] = set()
    for obligation in obligations:
        assert set(obligation.depends_on) <= known
        known.add(obligation.identifier)

    assert tuple(item.identifier for item in obligations) == tuple(
        f"P-{index:02d}" for index in range(1, 13)
    )
    assert {
        item.identifier
        for item in obligations
        if item.status is ClosureObligationStatus.ANALYTICALLY_PROVED
    } == {"P-02", "P-03", "P-04", "P-06"}
    assert render_closure_obligations_yaml(obligations).count(
        "ANALYTICALLY_PROVED"
    ) == 4
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
    assert decision.blocking_obligations == ("P-01", "P-05")
    rendered = render_minimal_closure_decision_yaml(decision)
    assert rendered.startswith("decision: NONE\nconfidence: high\n")
    assert "Fredholm" not in decision.prerequisite_statement
