# P1-C adversarial review

## Initial review

| Challenge | Status | Finding |
|---|---|---|
| Was \(\mathcal R_1\) confused with the full regularizer? | `PASS` | The candidate retains \(\mathcal T_{1,-}\) and \(Z_1^{-1}\) as ordered exterior factors |
| Was either interface lost? | `PASS` | Both are present in every candidate AST |
| Was the order reversed? | `PASS` | The sole candidate order is \(\mathcal T_{1,-}\mathcal R_1Z_1^{-1}\) |
| Was a different source link substituted? | `PASS` | The link is derived from the article's \(\mathcal A=Z\widehat{\mathcal A}\) and \(\mathcal N=\widehat{\mathcal A}\mathcal T\) identities |
| Was bilateral core regularity merely assumed from one side? | `PASS` | Article section 06 writes both oriented relations; each full product consumes the matching one |
| Were compact ideals mixed? | `PASS_WITH_ASSUMPTION` | Certification requires one typed \(\mathcal K(X)\); runtime mismatch checks remain to be implemented |
| Was compactness multiplied by an unbounded operator? | `PASS_WITH_ASSUMPTION` | The article proves bounded invertibility; typed runtime checks remain to be implemented |
| Was exact invertibility confused with Calkin invertibility? | `PASS` | \(Z_1\) and \(\mathcal T_{1,-}\) use exact inverse certificates; \(\mathcal R_1\) remains modulo compact |
| Was \(Z_1\) identified with KKL 2017's auxiliary \(V\)? | `PASS` | The bibliography audit explicitly classifies this as `NOT_A_MATCH` |
| Was the project formula attributed to a precedent source? | `PASS` | KKL 2016/2017 are precedent only; the exact factorization is article-internal |
| Were all factors absorbed into one Mellin PDO? | `PASS` | The full word remains a `Product`; only \(\mathcal R_1\)'s transported core has a Mellin-PDO representation |
| Was common-algebra membership asserted? | `PASS` | It remains `UNPROVED` |
| Was the weighted space lost? | `PASS_WITH_ASSUMPTION` | Every planned interface carries the same `WeightedLpSpace`; runtime validation remains to be implemented |
| Were domain or codomain altered? | `PASS_WITH_ASSUMPTION` | All article factors preserve \(X\); fail-closed checks remain to be implemented |
| Was the left proof reused as the right proof? | `PASS` | The second derivation independently uses \(\mathcal T_{1,-}^{-1}\) and \(\mathcal R_1\mathcal H_1\simeq I\) |
| Is the proposed status stronger than the evidence? | `PASS` | It matches the article's general proposition and its explicit \(k=1\) specialization |

## Refutation attempts

1. Dropping \(Z_1^{-1}\) fails because the exact link begins with
   \(Z_1^{-1}\mathcal A_{1,1}\), not \(\mathcal A_{1,1}\) alone.
2. Dropping \(\mathcal T_{1,-}\) fails because the Mellin core is defined as
   \(\widehat{\mathcal A}_{1,1}\mathcal T_{1,-}\).
3. Swapping the factors fails structurally in the noncommutative AST; no
   compact-commutator theorem for these concrete interfaces was found.
4. Using only \(\mathcal H_1\mathcal R_1\simeq I\) proves only the product with
   \(\mathcal A_{1,1}\) on the left. The other product needs the distinct
   \(\mathcal R_1\mathcal H_1\simeq I\) evidence.
5. Calling \(\mathcal R_1\) an exact inverse fails: both source relations are
   explicitly modulo \(\mathcal K(X)\).
6. Calling \(B_1\) a single Mellin PDO fails for lack of a theorem absorbing
   the exact interfaces into the Mellin calculus.

## Gate

There is no mathematical `FAIL` in the source-backed derivation. Items marked
`PASS_WITH_ASSUMPTION` are implementation gates: P1-C must remain fail-closed
until it verifies boundedness, spaces, exact inverse evidence, and one common
bilateral compact ideal.

## Post-implementation review

The runtime gates are now implemented and directly exercised:

| Attempted overclaim or malformed input | Final status |
|---|---|
| Wrong link order or duplicated role | `PASS`: rejected by `DomainOrOrderError` |
| Incompatible domain/codomain or weight | `PASS`: rejected or returned as `BLOCKED_BY_DOMAIN_OR_ORDER` before derivation |
| Different compact-ideal label or space | `PASS`: raises `CompactIdealMismatch` before residues are combined |
| Non-bilateral compact ideal | `PASS`: rejected before propagation |
| Unbounded exterior factor | `PASS`: rejected before applying the ideal property |
| Multiplier with a certified positive-half-line zero | `PASS`: exact inverse interface rejected |
| Auxiliary inverse supplied as a bare shift | `PASS`: type rejected |
| Missing inverse/link/regularizer evidence | `PASS`: fail-closed error |
| Only \(\mathcal H_1\mathcal R_1\simeq I\) | `PASS`: only the left product is certified |
| Only \(\mathcal R_1\mathcal H_1\simeq I\) | `PASS`: only the right product is certified |
| Conditional link or core regularizer | `PASS`: result remains `CONDITIONAL_TWO_SIDED_FACTORED_REGULARIZER` |
| Unverified link or core regularizer | `PASS`: explicit blocking status, no reconstructed relation |
| Schur insertion of a unilateral result | `PASS`: rejected |
| Expansion of the certified center during Schur insertion | `PASS`: no expanded expression is produced |
| Automatic common algebra, single PDO, or Schur symbol | `PASS`: fixed respectively at `UNPROVED`, `UNPROVED`, and `NOT_CERTIFIED` |

The left trace uses the `H1 R1` relation and the bounded pair
`Z1, Z1_inverse`. The right trace is a different tuple of steps and uses the
`R1 H1` relation and the bounded pair `T1_minus, T1_minus_inverse`. There is no
code path that copies the left semantic conclusion into the right field.

No final adversarial item is classified `FAIL`. The theorem assumptions stay
visible and therefore the result is certified under, rather than independently
of, the article's Fredholm and weighted-space hypotheses.
