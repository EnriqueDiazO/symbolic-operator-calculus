# P1-B adversarial review

## Review outcome

The implementation review has no `FAIL`.  It has `OPEN` mathematical inputs,
all of which are exposed as obligations or blocking statuses.  Therefore the
software phase is reproducible, while the article-level factorization is not
claimed complete.

| # | Audit target | Status | Adversarial finding |
|---:|---|---|---|
| 1 | Order of the word | `PASS` | `original_word` is the supplied five-atom `Product`; expanded and factorized candidates are separate fields. |
| 2 | Dilation orientation | `PASS` | Both scales are typed, positive, normalized, and checked as inverse. Same-orientation pairs block. |
| 3 | Cutoff rescaling | `PASS` | Exact covariance produces `chi(gamma x)` and a new radial scale. It never equates this cutoff with `chi(x)`. |
| 4 | Coefficient rescaling | `PASS` | Coefficients remain typed multiplication operators; no coefficient is crossed in the full-word recognizer. The exact multiplication covariance remains separately available. |
| 5 | Wiener--Hopf rescaling | `PASS` | The prior exact frequency-rescaling rule is reused only in its supported scope. Localized WH objects are distinct typed objects. |
| 6 | Position of `P^pm` | `PASS` | `T_{1,-}=U_1^-1P^+ +P^-` remains a `LinearCombination`. A projection between an exterior and dilation blocks the pattern. |
| 7 | `Phi^-1 Op_M(r) Phi` interface | `PASS` | `TransportedMellinCore` records all three atoms. A bare or assumed representation is blocked. |
| 8 | Compact ideal propagation | `PASS` | Certified modulo-compact inputs must share one explicit ideal; an assumed cutoff remains formal. |
| 9 | Matrix compatibility | `PASS` | Existing `ExactBlock(A,2,1)` and `ExactBlock(A,1,2)` metadata are reused by the cancellation bridge. |
| 10 | Sign convention | `PASS` | Raw product, signed contribution, and pivot are separate types; embedded raw minus signs are rejected. |
| 11 | Domain/codomain | `PASS` | Every accepted contiguous factor preserves the same full `WeightedLpSpace`, hence the same `p`, weight, assumptions, and norm convention. |
| 12 | Overcertification | `PASS` | Blocked results cannot carry a factorized word or semantic relation; status is not inferred from names. |

## Source-interface attack

The most dangerous input uses the article names while following the prompt's
five-argument grouping.  Typed metadata says that the left `B^-` contains the
left auxiliary and the right `B^+` contains the right auxiliary.  Expansion
then shows

```text
... T1_minus T1_minus R1 Z1_inverse Z1_inverse ...
```

and returns `BLOCKED_BY_REGULARIZER_INTERFACE`.  A test with coincidentally
similar atom names but without the typed roles does **not** trigger this
classification; it is instead rejected by the actual operator pattern.  This
guards against string-based certification.

## Cutoff attacks

- passing an ordinary multiplication operator where a cutoff is required:
  rejected by type;
- passing a boolean cutoff: rejected by type;
- omitting claims: `BLOCKED_BY_CUTOFF_CONTROL`;
- supplying one unilateral claim for a two-sided localization: still blocked;
- supplying an `UNPROVED` claim: still blocked;
- supplying an `ASSUMED` bilateral claim: candidate only as `FormalIdentity`;
- supplying a `CERTIFIED` bilateral claim: at most modulo compact operators in
  its named ideal.

No compactness is inferred from overlapping or separated-looking support
descriptions.

## Order and interface attacks

- a `P^+` or `P^-` between exterior and dilation: blocked;
- an extra factor on either core: blocked;
- a mismatched weighted space: blocked;
- two forward or two inverse dilation orientations: blocked;
- a bare `RegularizerMellinRepresentation` without explicit `Phi` transport:
  blocked;
- an assumed transported regularizer: blocked;
- localized exterior supplied to the old cancellation API: never delegated;
- no actual adjacent pair: no cancellation call.

## Open mathematical rows

| Question | Status | Consequence |
|---|---|---|
| Correct nonduplicating primitive interface for the prompt | `OPEN` | The literal requested word remains blocked. |
| Actual left-core summand factorization | `OPEN` | No article-level `L_{2,1}V_gamma Phi^-1` is asserted. |
| Actual right-core inverse orientation | `OPEN` | No article-level `Phi V_gamma^-1 R_{1,2}` is asserted. |
| Scaled-cutoff replacement in the mixed word | `OPEN` | Verified two-cutoff localization is retained. |
| Common WH--Mellin--dilation algebra | `OPEN` | No closure, inverse-closedness, or Fredholm symbol is inferred. |
| Full correction/pivot symbol | `OPEN` | No final determinant or index formula is asserted. |

## Final adversarial decision

The software correctly demonstrates all six statuses on controlled typed
inputs and correctly stops on the article-labelled five-factor input.  The
single article classification is
`BLOCKED_BY_REGULARIZER_INTERFACE`.  The pivot is not declared solved.
