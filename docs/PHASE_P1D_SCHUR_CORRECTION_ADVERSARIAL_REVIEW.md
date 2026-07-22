# P1-D adversarial review

## Source and word integrity

| Challenge | Verdict | Required behavior |
|---|---|---|
| Confuse \(\mathcal R_1\) with the full P1-C regularizer | `PASS` | grouped center is exactly \(\mathcal T_{1,-}\mathcal R_1Z_1^{-1}\) |
| Reuse the P1-B literal five-factor interface | `PASS` | reject duplication of \(\mathcal T_{1,-}\) and \(Z_1^{-1}\) |
| Swap an exterior factor across the center | `PASS` | typed AST order mismatch blocks assembly |
| Match factors only by printed names | `PASS` | compare typed objects, roles, branches, spaces, and exact ASTs |
| Distribute either \(\widetilde V-\mathcal G\) difference | `PASS` | retain the two linear combinations as subtrees |
| Drop any of the four cutoffs | `PASS` | actual expanded word retains all four |
| Treat source exactness as pivot exactness | `PASS` | correction assembly exact; original/model pivot relation only modulo compact |

## Sign attacks

1. The raw Schur definition carries one external minus.
2. The normalized model for \(\mathcal A_{2,1}\) carries one internal minus.
3. The normalized model for \(\mathcal A_{1,2}\) carries a plus.
4. Therefore the model correction enters the pivot with a plus.

Changing either normalized sign, attaching the Schur sign to the correction
atom itself, or claiming an exact relation between original and model pivot
must fail. The correction remains intrinsically the positive word (G); sign
bookkeeping is owned by the pivot certificate.

## Cutoff attacks

| Attempt | Expected result |
|---|---|
| Keep each localized block intact | `NO_REPLACEMENT_NEEDED` |
| Use exact dilation covariance to replace \(\chi(\gamma x)\) by \(\chi(x)\) | `CutoffReplacementNotCertified` |
| Claim compact equivalence from support geometry alone | `CutoffReplacementNotCertified` |
| Transport the cutoff by an explicit exact covariance identity | `EXACT_COVARIANCE_ONLY`; this still is not a replacement |
| Provide a sourced common-ideal compact certificate | `REPLACEMENT_MOD_COMPACTS` |
| Supply only an assumption | `REPLACEMENT_ASSUMED`, which cannot support a certified symbol |

Exact covariance transports a cutoff; it never deletes it.

## Dilation attacks

| Attempt | Expected result |
|---|---|
| Cancel a typed adjacent reciprocal pair | `EXACT_CANCELLATION` |
| Apply a fully evidenced exact covariance trace | `EXACT_COVARIANCE_REDUCTION` |
| Cancel through \(\mathcal R_1\), \(Z_1^{-1}\), a projection, multiplier, or sum | `BLOCKED_BY_INTERVENING_FACTORS` |
| Infer a pair merely because two branch scales occur | `NO_RELATIVE_DILATION_PAIR` |
| Expand sums to expose a desirable pair | rejected; no distributive search |

The actual word contains no certified cancellation opportunity.

## Algebra attacks

The following implications are invalid and must remain unimplemented:

- every factor is bounded, therefore the word belongs to the desired algebra;
- factors occur in several well-studied algebras, therefore their product
  belongs to one of them;
- exact dilation covariance implies membership of the dilation in the base
  algebra;
- invertibility of \(\mathcal T_{1,-}\) implies inverse-closedness of a mixed
  algebra containing it;
- the 2025 half-line algebra's treatment of \(\mathcal N_1\) automatically
  contains the chosen general radial regularizer \(\mathcal R_1\);
- naming a crossed-product extension proves its closure or Fredholm theorem;
- DLS 1995 permits deletion of both concrete localization cutoffs;
- a Fredholm theorem for a different exact algebra applies to this word.

The article instance must expose the missing generator, closure, cutoff,
dilation, inverse-closedness, and symbol-homomorphism fields, with status
`BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP`.

## Symbol attacks

`derive_schur_correction_symbol` must raise or return a typed blocker before
forming any formula if a prerequisite is absent. Dedicated error paths cover:

- `AlgebraMembershipNotCertified`;
- `CutoffTreatmentNotCertified`;
- `RelativeDilationNotResolved`;
- `NoncommutativeProductRuleMissing`;
- `SchurSymbolNotAvailable` for a blocked result requested as certified.

A successful synthetic path requires factor symbols in the exact expanded
order and a certified noncommutative product rule. It returns a structured
ordered-symbol object, not an uncontrolled commutative SymPy product.

## Fredholm overclaim attack

The article proves that the full two-branch operator is Fredholm if and only
if the Schur pivot is Fredholm under the stated hypotheses. It does not prove
that the pivot is Fredholm. The P1-D pivot certificate must therefore use
`NOT_CERTIFIED` even though the correction assembly is exact.

## Implementation gate

The post-implementation review must test malformed domains, reversed branch
maps, source-block disagreement, regularizer statuses weaker than P1-C's
bilateral certificate, compact-ideal mismatch, cutoff replacement requests,
intervening dilation factors, incomplete algebra evidence, and missing symbol
rules. No blocking result may carry a certified semantic relation or a
candidate symbol.

## Post-implementation review

| Refutation attempt | Final classification | Runtime evidence |
|---|---|---|
| Reverse the two article source blocks | `PASS` | `BLOCKED_BY_SOURCE_ORDER_MISMATCH`, with no semantic relation |
| Omit one source factor | `PASS` | `BLOCKED_BY_MISSING_FACTOR` |
| Reverse a branch orientation | `PASS` | `BLOCKED_BY_DOMAIN_OR_CODOMAIN` before source matching |
| Lose either P1-C interface | `PASS` | the seven-factor AST is compared exactly with both nested source definitions |
| Expand the two differences distributively | `PASS` | both differences remain authoritative atoms with separate exact linear-combination definitions |
| Delete or merge cutoffs | `PASS` | four audit records, all `NO_REPLACEMENT_NEEDED` |
| Promote exact cutoff covariance to replacement | `PASS` | `require_cutoff_replacement` raises `CutoffReplacementNotCertified` |
| Cancel the nonadjacent \(U_1^{-1}/U_1\) pair | `PASS` | `BLOCKED_BY_INTERVENING_FACTORS`, no relation |
| Mix Karlovich-2025 and KKL-2014 component coverage | `PASS` | common algebra remains `BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP` |
| Request a symbol after that algebra result | `PASS` | `BLOCKED_BY_ALGEBRA_MEMBERSHIP`, `symbol is None` |
| Omit one of the four cutoff audits in a synthetic certified algebra | `PASS` | `BLOCKED_BY_CUTOFFS` |
| Omit the ordered product rule | `PASS` | `BLOCKED_BY_NONCOMMUTATIVE_PRODUCT` |
| Supply a complete synthetic sourced algebra package | `PASS` | the positive path returns an ordered structured symbol, demonstrating that the gate is reachable |
| Change the normalized \(A_{2,1}\) sign | `PASS` | pivot status `BLOCKED_BY_SIGN_CONVENTION` |
| Infer pivot Fredholmness from assembly | `PASS` | fixed `NOT_CERTIFIED`, even in the successful assembly path |

The focused suite contains 28 tests in addition to all pre-existing P1-B/P1-C
regressions. No adversarial item is classified `FAIL`. The remaining source
gaps are represented as proof obligations rather than assumptions.
