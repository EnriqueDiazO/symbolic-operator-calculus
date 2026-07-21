# First-pivot result × interface matrix

This matrix evaluates only contiguous products in the stored noncommutative
words. “Verified” means that the stated fact was found either in the paper
definition or in the checked source result; it does not mean that every
hypothesis needed for application has been proved. No row changes a
modulo-compact relation into an exact equality.

| Interface | Candidate result | Required hypotheses | Verified | Missing | Status |
|---|---|---|---|---|---|
| $P^\pm \times R_1$ | KKL2014TwoShifts Lemma 2.7 | weighted conjugated commutation in the paper's space | P^pm are exact Mellin convolution factors | source theorem is unweighted and concerns the Cauchy algebra | `NO_RULE` |
| $U_1^{-1}P^+ \times R_1$ | KKL2014TwoShifts Lemma 2.7 and Theorem 5.2 | weighted ordered rule including the exterior dilation | the factor order is represented exactly | no source result covers this weighted ordered triple | `BLOCKED` |
| $R_1 \times U_1$ | KKL2014TwoShifts Theorem 3.3 or a proved extended calculus | both symbols belong to one composition calculus | R1 has the supplied E-tilde symbol | gamma1^(i lambda) is outside V when gamma1 is nontrivial | `BLOCKED` |
| $R_1 \times \widehat G_1$ | KKL2014TwoShifts Theorem 3.3 | r1 and the lambda-independent Ghat1 symbol belong to E | r1 is in E-tilde; Ghat1 is an SO multiplication coefficient | paper-to-source symbol-class specialization | `SPECIALIZATION_TO_PROVE` |
| $U_1 \times W_{1,2}^{+}$ | KKL2014TwoShifts Lemmas 4.4--4.5 | Wplus_12 is the source operator R_y after all transports | both kernels have fixed-singularity structure | kernels, cutoff, weight, and branch maps are different | `ANALOGY_ONLY` |
| $\widehat G_1 \times W_{1,2}^{+}$ | Karlovich2017HasemanSO Theorem 5.6 | log-star block theorem specializes to the localized cusp block | the coefficient and block order are explicit | geometry, localization, and symbol identification | `ANALOGY_ONLY` |
| $R_1B_1 \times W_{1,2}^{+}$ | H2 via Mellin semiproduct or cusp-algebra closure | H1 and a compatible class for Wplus_12 | the suffix is contiguous in every word | H1 and Wplus_12 class membership | `BLOCKED` |
| $Q_1 \times R_1B_1W_{1,2}^{+}$ | H3 via a weighted ordered exterior rule | H2 and a rule for each Q1 | the full right core is represented without reordering | H2 and the weighted Q1 interface | `BLOCKED` |
| $complete word$ | Mellin-PDO membership or Karlovich2025 cusp-algebra membership | all interfaces close in one verified class; residues are compact | the AST fixes all four products and their domains | P-01 through P-08 | `BLOCKED` |

## Source-result assessment

- **KKL2014TwoShifts, Theorem 3.3 (printed p. 942; PDF p. 8).** Its
  semiproduct conclusion requires both symbols to lie in
  \(\mathcal E(\mathbb R_+,V(\mathbb R))\). This supports a
  specialization to prove for \(R_1\widehat G_1\), after recording the
  lambda-independent multiplication symbol. It does not apply to
  \(R_1U_1\), because the paper explicitly notes that
  \(\gamma_1^{i\lambda}\notin V(\mathbb R)\) for a nontrivial dilation.
- **KKL2014TwoShifts, Lemma 3.4 (printed p. 942; PDF p. 8).** The exact
  separated triple product requires all three symbols to lie in the source
  class, with radial-only and covariable-only outer symbols. The present word
  has not been put in that form; the lemma is not applicable yet.
- **KKL2014TwoShifts, Lemmas 4.4–4.5 (printed pp. 947–950; PDF pp. 13–16).**
  These concern the particular fixed-singularity operator \(R_y\).
  \(W_{1,2}^{+}\) has not been identified with \(R_y\), so the comparison
  is an unproved adaptation only.
- **KKL2014Regularization, Theorem 5.8 (printed p. 206; PDF p. 18).** It
  supplies a Mellin-PDO regularizer after its nonvanishing hypotheses. It
  supports the declared status of \(R_1\), not closure of products on either
  side of \(R_1\).
- **Karlovich2017HasemanSO, Theorem 5.6 (printed p. 479; PDF pp. 17–18).**
  The exact shifted off-diagonal representation belongs to a logarithmic-star
  model. It is an analogy for the localized pure-cusp block.
- **Karlovich2017HasemanSO, Theorem 6.2 (printed p. 482; PDF pp. 20–25).**
  Its global logarithmic-star auxiliary operators are not the weighted
  branchwise auxiliaries of the paper; no closure row is certified by it.
- **Karlovich2025Cusps, Lemmas 3.1–3.2 and Corollary 3.3 (article pp. 8–13).**
  They support pure-cusp localization and compact block remainders after a
  paper-specific specialization. They do not place \(R_1\) or the complete
  word in the 2025 algebra.
- **Karlovich2025Cusps, Theorems 5.1, 6.1, and 6.2 (article pp. 19–25).**
  Algebraic closure, symbols, and Fredholm criteria are available only after
  membership in the source algebra has been proved. That antecedent is open.
- **Karlovich2009MellinPDO, Theorem 4.3 (printed p. 88; PDF p. 8).** Only
  its scalar scope is independently verified here. It is deliberately not
  applied before Mellin-algebra membership.

## Matrix conclusion

The only near-term source specialization is the
\(R_1\widehat G_1\) row. It does not settle the paired
\(R_1U_1\) case needed by H1. Every row containing
\(W_{1,2}^{+}\) remains either analogy-only or blocked.
