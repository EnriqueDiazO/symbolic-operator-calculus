# Audit of the regularizer–coefficient semiproduct

## Conclusion

Under the paper's explicit hypothesis \(G_1\in SO(\mathbb R_+)\), the ordered
coefficient case is

\[
 \boxed{
 \operatorname{Op}(r_1)M_{\widehat G_1}
 \simeq
 \operatorname{Op}(r_1\widehat G_1)
 }
 \qquad\text{on }L^p(\mathbb R_+,d\mu).
\]

After conjugation this gives the corresponding statement for
\(R_1\widehat{\mathcal G}_1\) on the weighted space. Its status is
`CERTIFIED_MOD_COMPACT`, not exact equality. This audit is independent of the
exact, factorized treatment of \(R_1U_1\).

## Verified evidence

| Result | Printed/PDF pages | Checksum | Use |
|---|---|---|---|
| `KKL2014TwoShifts`, Theorem 3.3 | 942 / 8 | `3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec` | If \(a,b\in\mathcal E(\mathbb R_+,V(\mathbb R))\), then \(\operatorname{Op}(a)\operatorname{Op}(b)\simeq\operatorname{Op}(ab)\). |
| `KKL2014TwoShifts`, Lemma 3.4 | 942 / 8 | `3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec` | Exact separated triple product; audited but not used for this two-factor formula. |
| `KKL2014Regularization`, Theorem 5.8 | 206 / 18 | `7ed69d096f9e31983cb9b4701e8b486fd15736d82163144f27f9bae244ce41f6` | Provenance for the regularizer symbol in \(\widetilde{\mathcal E}\) under its stated regularization hypotheses. |

These records were consumed through `litctl`; no Fredholm criterion from the
sources is used in Phase Q.

## Hypothesis-by-hypothesis check

1. The paper constructs
   \(R_1=\Phi_\delta^{-1}\operatorname{Op}(r_1)\Phi_\delta\), with
   \(r_1\in\widetilde{\mathcal E}(\mathbb R_+,V(\mathbb R))\), hence
   \(r_1\in\mathcal E\).
2. The paper's normalized coefficient satisfies
   \(\widehat G_1=\zeta_1^{-1}G_1\in SO(\mathbb R_+)\) when
   \(G_1\in SO(\mathbb R_+)\).
3. The paper explicitly proves that
   \((x,\lambda)\mapsto\widehat G_1(x)\) belongs to
   \(\widetilde{\mathcal E}\): it is independent of \(\lambda\), so its
   covariable derivative and translation difference vanish.
4. Multiplication by \(\widehat G_1\) commutes with the multiplication
   transport \(\Phi_\delta\). Thus both factors act on the unweighted model
   \(L^p(\mathbb R_+,d\mu)\) in the form required by Theorem 3.3.
5. The theorem is applied in the stored order \(r_1\) followed by
   \(\widehat G_1\); no factors are commuted.
6. Since \(\widetilde{\mathcal E}\) is a Banach algebra, the pointwise product
   \(r_1\widehat G_1\) retains the declared standard symbol membership.

Every required membership, space, weight transport, dependence, and order
hypothesis has therefore been checked. Lemma 3.4 is unnecessary and is not
used to upgrade the result to an equality.

## Registered rule

| Field | Value |
|---|---|
| name | `regularizer_slowly_oscillating_coefficient_semiproduct` |
| logical type | `MOD_COMPACT_EQUIVALENCE` |
| preconditions | both symbols in the verified Mellin class; paper's \(\Phi_\delta\) transport; \(G_1\in SO\) |
| source | `KKL2014TwoShifts`, Theorem 3.3 |
| pages/checksum | printed 942; PDF 8; checksum above |
| status | `CERTIFIED_MOD_COMPACT` |

The compact residue is supplied by the cited semiproduct theorem. No statement
about \(W_{1,2}^{+}\), membership of the longer word, or Fredholmness follows.
