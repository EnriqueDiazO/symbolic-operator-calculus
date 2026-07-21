# Decision on the minimal first-pivot closure lemma

## Candidate status

These are research candidates, not established theorems. In every formula,
\(K\) denotes a compact residue only if its compactness is separately proved.

### H1

For each \(B_1\in\{U_1,\widehat G_1\}\), seek a symbol \(c_{B_1}\) in one
explicitly named admissible class such that

\[
 R_1B_1=\operatorname{Op}(c_{B_1})+K_{B_1}.
\]

- **Hypotheses:** the supplied Mellin representation of \(R_1\); membership of
  both choices of \(B_1\) in the same verified calculus; an ordered
  semiproduct theorem; compactness of the residue.
- **Sources:** KKL2014TwoShifts, Theorem 3.3 and Lemma 3.4.
- **Established:** Phase Q proves the exact ordered factorized identity for
  \(R_1U_1\), and certifies the \(R_1\widehat G_1\) semiproduct modulo
  compact operators under the paper hypothesis.
- **Formal:** the contiguous subword occurs without reordering in all 16
  correction terms.
- **Open:** the nontrivial dilation multiplier remains outside \(V\), the
  factorized output has no proved closed calculus, and no one admissible
  framework for both choices has been identified.
- **Would cover:** the \(R_1B_1\) subword in all 16 terms, but not the following
  \(W_{1,2}^{+}\) or the exterior \(Q_1\).
- **Risk:** an unjustified \(V\)-membership assertion for \(U_1\).

Status: `BLOCKED`.

### H2

For each \(B_1\), seek \(c_{B_1}^{+}\) such that

\[
 R_1B_1W_{1,2}^{+}
 =\operatorname{Op}(c_{B_1}^{+})+K_{B_1}^{+}.
\]

- **Hypotheses:** H1; a compatible representation of
  \(W_{1,2}^{+}\); lateral closure; compact cutoff residues.
- **Sources:** Karlovich2017HasemanSO, Theorem 5.6;
  KKL2014TwoShifts, Lemmas 4.4–4.5; Karlovich2025Cusps,
  Corollary 3.3.
- **Established:** the paper's exact localized Wiener--Hopf definition.
- **Formal:** this contiguous suffix occurs in all 16 terms.
- **Open:** H1 and P-01/P-05.
- **Would cover:** the suffix in all 16 terms, but not the exterior \(Q_1\).
- **Risk:** silently identifying \(W_{1,2}^{+}\) with \(R_y\), a Mellin PDO,
  or a 2025 generator.

Status: `BLOCKED`.

### H3

For every pair
\((Q_1,B_1)\in\{U_1^{-1}P^+,P^-\}\times\{U_1,\widehat G_1\}\), seek
\(c_{Q_1,B_1}^{+}\) such that

\[
 Q_1R_1B_1W_{1,2}^{+}
 =\operatorname{Op}(c_{Q_1,B_1}^{+})+K_{Q_1,B_1}^{+}.
\]

- **Hypotheses:** H2; a weighted ordered rule for both \(Q_1\) choices;
  compactness and identified symbol class.
- **Sources:** KKL2014TwoShifts, Lemma 2.7 and Theorem 5.2, only as unproved
  weighted adaptations.
- **Established:** the four full right cores are exact AST products.
- **Formal:** the four patterns account for the recurring core of all 16 terms.
- **Open:** P-01 through P-08.
- **Would cover:** the complete right core in all 16 terms.
- **Risk:** bundling three unrelated missing interfaces and hiding the first
  failed hypothesis.

Status: `BLOCKED`.

## Explicit criteria

| Criterion | H1 | H2 | H3 | Consequence |
|---|---|---|---|---|
| Appears in all 16 terms | yes | yes | yes | no discriminator |
| Preserves factor order | yes | yes | yes | all are formally safe |
| Has currently verifiable hypotheses | two separate short relations, but no uniform class | no | no | no candidate is defensible now |
| Reduces the most obligations | local reduction | larger suffix | whole core | breadth cannot replace missing hypotheses |
| Avoids a false identification | yes if class is enlarged and proved | no, currently depends on P-01 | no, inherits P-01 | excludes H2/H3 |
| Does not anticipate Fredholmness | yes | yes | yes | all are scoped below P-12 |

## Machine-readable decision

```yaml
decision: NONE
confidence: high
blocking_obligations:
  - "P-01"
  - "P-04"
evidence:
  - "phase-q:exact-right-dilation-definition-proof"
  - "phase-q:certified-Ghat1-semiproduct"
  - "paper:normalized-wh-blocks-section-six"
  - "KKL2014TwoShifts:3.3"
  - "KKL2014TwoShifts:4.4"
  - "Karlovich2025Cusps:3.3"
rationale: "The two H1 cases now have separate certified relations, but the exact factorized U1 output has no proved closed admissible calculus. H2 and H3 also require an unproved identification of the localized Wplus_12."
prerequisite_statement: "Prove the minimum closure properties needed for the factorized pair (r1,d_gamma1), without forcing it into E-tilde; independently identify Wplus_12 before attempting H2."
```

## Selected statement

The selected option remains **NONE**. Consequently, no uniform H1/H2/H3
equality-plus-compact statement is promoted to a lemma after Phase Q.

The next statement to formulate and prove is a prerequisite, not a disguised
version of H1:

> Prove just the closure properties actually needed for the ordered factorized
> pair \((r_1,d_{\gamma_1})\), without forcing its pointwise product into
> \(\widetilde{\mathcal E}\). Independently identify \(W_{1,2}^{+}\) before
> attempting H2.

The exact \(R_1U_1\) composition itself is no longer open, and the coefficient
semiproduct is no longer merely plausible. What remains open is the common
closure framework required by H1 and every longer word.

## Mathematical status and impact

Confidence remains **high** because the two resolved cases are now recorded at
their distinct logical strengths, while the uniform closure antecedent is
still absent. The four words and their occurrence in the 16 terms are formally
certified. None of those 16 complete terms thereby acquires Mellin-PDO
membership, cusp-algebra membership, an identified final symbol, or a compact
remainder. Reconstruction of
\(\mathcal C_2\), membership of \(N_2^{(1)}\), and Fredholm analysis remain at
P-10, P-11, and P-12 respectively.

## Preserved history

Phase P selected `NONE` because \(R_1U_1\) had no proved ordered composition
rule and the coefficient case was only a specialization to prove. Phase Q
discharges precisely those two short calculations. It does not retroactively
claim that Phase P had proved them, and it does not modify the Phase P finding
that \(W_{1,2}^{+}\) is unidentified.
