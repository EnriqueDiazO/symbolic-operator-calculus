# Audit of the left Cauchy-factor Mellin semiproduct

## Conclusion

If \(a\in\widetilde{\mathcal E}(\mathbb R_+,V(\mathbb R))\), then after
the paper's weighted conjugation,

\[
 \boxed{
 P^\pm\operatorname{Op}(a)
 \simeq
 \operatorname{Op}(p^\pm a)
 }
 \qquad\bmod\mathcal K.
\]

The status is `CERTIFIED_MOD_COMPACT`. The order is always \(p^\pm a\); no
factor is commuted.

## Exact Cauchy representation

The paper defines

\[
 P^\pm=\frac12(I\pm S_{\mathbb R_+})
\]

and proves exactly

\[
 \Phi_\delta P^\pm\Phi_\delta^{-1}
 =\operatorname{Op}(p^\pm),
\]

where

\[
 p^\pm(\lambda)=\frac12\left(1\pm
 \coth(\pi(\lambda+i\kappa_\delta))\right).
\]

These are Mellin convolution symbols independent of \(x\). They belong to
\(V(\mathbb R)\); spatial constancy gives bounded radial continuity and zero
radial oscillation. Translation continuity in the \(V\)-norm follows from
absolute continuity, uniform continuity, and continuity of translations of
their \(L^1\) derivatives. Their explicitly exponentially decaying
derivatives give the uniform tail condition. Hence
\(p^\pm\in\widetilde{\mathcal E}\subset\mathcal E\).

This does not make \(P^\pm\) complementary projections. The paper proves

\[
 P^-P^+=P^+P^-\ne0.
\]

No zero-product simplification is used.

## Source result and provenance

`KKL2014TwoShifts`, Theorem 3.3 states that, on
\(L^p(\mathbb R_+,d\mu)\),

\[
 \operatorname{Op}(b)\operatorname{Op}(a)
 \simeq\operatorname{Op}(ba)
\]

when both symbols belong to
\(\mathcal E(\mathbb R_+,V(\mathbb R))\).

| Field | Value |
|---|---|
| printed page | 942 |
| PDF page | 8 |
| checksum | `3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec` |
| application kind | direct specialization after paper membership checks |
| logical strength | modulo compact operators |

All hypotheses are satisfied for \(b=p^\pm\) and the standard symbols used
in Phase R. Since `E-tilde` is a Banach algebra, \(p^\pm a\) remains in that
strengthened class.

Lemma 3.4 from the same page was also audited. Its exact triple-product result
is unnecessary here and is not applied to any oscillatory dilation factor.
`Karlovich2017HasemanSO`, Theorem 3.5 was inspected as required but is a
Fredholm criterion; it is not used.

## Registered rules

The public trace records `cauchy_plus_mellin_semiproduct` and
`cauchy_minus_mellin_semiproduct`, each with:

- logical type `MOD_COMPACT_EQUIVALENCE`;
- the verified symbol, space, and order preconditions;
- Theorem 3.3, printed/PDF pages and checksum;
- status `CERTIFIED_MOD_COMPACT`.

No Fredholmness, right Wiener--Hopf composition, or standard membership for a
symbol containing \(d_\gamma\) follows from these rules.
