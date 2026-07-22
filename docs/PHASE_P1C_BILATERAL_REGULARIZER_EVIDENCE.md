# P1-C evidence matrix

## Classification vocabulary

The audit uses only `ARTICLE_INTERNAL_EXACT`,
`ARTICLE_INTERNAL_MOD_COMPACTS`, `EXACT_THEOREM`,
`DIRECT_ALGEBRAIC_CONSEQUENCE`, `PLAUSIBLE_ADAPTATION`, and `SOURCE_GAP`.

| Claim | Classification | Evidence and scope |
|---|---|---|
| \(\widetilde V_1=Z_1U_1\) | `ARTICLE_INTERNAL_EXACT` | Article section 05, lines 339-408 |
| \(\mathcal A_{1,1}=Z_1\widehat{\mathcal A}_{1,1}\) without commuting \(Z_1\) and \(P^\pm\) | `ARTICLE_INTERNAL_EXACT` | Article section 05, lines 429-507 |
| \(\mathcal T_{1,-}=U_1^{-1}P^++P^-\) | `ARTICLE_INTERNAL_EXACT` | Article section 05, lines 515-529 |
| Exact invertibility of \(\mathcal T_{1,-}\) | `EXACT_THEOREM` | Article theorem and proof, section 05, lines 532-782 |
| \(\mathcal H_1:=\mathcal N_1=\widehat{\mathcal A}_{1,1}\mathcal T_{1,-}\) | `ARTICLE_INTERNAL_EXACT` | Article sections 05:792-900 and 06:109-132; `H1` is explicitly an API alias |
| \(Z_1^{-1}\mathcal A_{1,1}\mathcal T_{1,-}=\mathcal H_1\) | `DIRECT_ALGEBRAIC_CONSEQUENCE` | Exact substitution and cancellation from the preceding article identities |
| \(\mathcal R_1=\Phi^{-1}\operatorname{Op}(r_{1,1})\Phi\) | `ARTICLE_INTERNAL_EXACT` | Article section 06, lines 161-168 |
| \(\mathcal H_1\mathcal R_1\simeq I\) | `ARTICLE_INTERNAL_MOD_COMPACTS` | Article section 06, lines 134-171; regularizer input attributed there to KKL 2014, Theorem 5.8 |
| \(\mathcal R_1\mathcal H_1\simeq I\) | `ARTICLE_INTERNAL_MOD_COMPACTS` | Same independent, oppositely oriented relation |
| \(\mathcal A_{1,1}B_1\simeq I\) | `DIRECT_ALGEBRAIC_CONSEQUENCE` | Article section 05:1473-1493 and specialization 06:209-227 |
| \(B_1\mathcal A_{1,1}\simeq I\) | `DIRECT_ALGEBRAIC_CONSEQUENCE` | Article section 05:1495-1529 and specialization 06:229-252 |
| Full word belongs to one common operator algebra | `SOURCE_GAP` | Neither article nor audited bibliography proves membership for the complete \(\mathcal T_{1,-}\mathcal R_1Z_1^{-1}\) word |
| Full word is one Mellin PDO | `SOURCE_GAP` | False as an automatic inference; no source proves absorption of both interfaces |
| A Schur-correction symbol follows now | `SOURCE_GAP` | Requires the still-open mixed-algebra and localization work |

## Directed bibliography audit

The read-only source is
`haseman-literature-index/reports/P1C_SOURCE_INGESTION_AND_REGULARIZER_AUDIT.md`,
cross-checked against `catalog/theorem_index.yaml` and `catalog/uses.yaml`.

| Source result | Verification | Permitted use in P1-C |
|---|---|---|
| KKL 2014, Regularization Theorem 5.8 | Statement and proof marked verified in the theorem index | Direct source for the two-sided Mellin-PDO regularizer invoked by article section 06, once the article's symbol hypotheses are assumed |
| KKL 2016, Theorem 7.1 | Verified | Factorized-regularizer precedent only; not the project link identity |
| KKL 2017, Theorem 3.3 | Verified, including proof | The factors `F` and `V` are only a structural precedent. `T1_minus` versus `F` is a plausible match; `Z1` versus `V` is not a match |
| KKL 2017, Corollary 3.4 | Verified | General transport of regularizers in that source, not the formula proved in this article |
| Karlovich 2022, Theorems 5.2--5.4 | Verified | Calkin/common-algebra conclusions only after exact algebra membership; no permission to absorb \(T_1\) or \(Z_1\) |
| Karlovich 2021 QC Mellin-PDO results | Verified in the index | Conditional symbol-class context; not a proof about the full factorized word |
| Kravchenko--Litvinchuk 1994 | Scanned/book-level record | Historical precedent, not direct certified applicability here |
| Karlovich 2006 algebra results | Catalogued | No audited theorem establishes the complete P1-C word's membership |

The article's own factorization and the two full-product calculations are
internal exact/modulo-compact arguments. They are not attributed to KKL 2017.

## Hypotheses retained

The bilateral conclusion is conditional on the Fredholm hypotheses stated in
`sections/06_shur_reduction.tex:101-107`, including the symbol assumptions
needed to invoke the KKL 2014 regularization theorem. It also retains the
weighted-space assumptions ensuring bounded exact invertibility of \(Z_1\)
and \(\mathcal T_{1,-}\). The software records these hypotheses; it does not
prove them from raw coefficient data.

## Claims deliberately left open

- `full_regularizer_algebra_membership`
- `single_mellin_pdo_representation`
- `cutoff_replacement_mod_compacts`
- `wh_mellin_wh_sandwich_membership`
- `fredholm_algebra_membership`
- `schur_correction_symbol`
