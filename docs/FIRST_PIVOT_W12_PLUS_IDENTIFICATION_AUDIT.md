# Audit of the identification of \(W_{1,2}^{+}\)

## Result

```yaml
conclusion: NOT_IDENTIFIED
relation_to_R_y: NONE_PROVED
mellin_membership: OPEN
cusp_algebra_membership: OPEN
```

This audit compares definitions, kernels, transports, cutoffs, spaces, and
source hypotheses. A resemblance between fixed-singularity kernels is not used
as an operator identity or a modulo-compact equivalence.

## 1. Exact paper definition

After the paper's branchwise normalization, the positive off-diagonal model is

\[
 W_{1,2}^{+}
 =\chi_\infty
   W\!\left(\chi_+(\cdot)e^{(c_1-c_2)(\cdot)}\right)
   \chi_\infty I.
\]

This is an exact definition of the model operator. The replacement of the
transported off-diagonal block by this model is a separate relation modulo
compact operators.

## 2. Kernel

Before localization, the paper's Fourier inversion gives a kernel proportional
to

\[
 \frac{1}{\pi i}\frac{1}{\tau-x+i(c_1-c_2)}.
\]

The actual model includes multiplication by \(\chi_\infty\) at both sides.
Consequently, neither the unlocalized kernel alone nor its singularity location
determines the localized operator.

## 3. Wiener--Hopf symbol

The supplied Wiener--Hopf multiplier is

\[
 \chi_+(\xi)e^{(c_1-c_2)\xi}.
\]

This is the symbol in the paper's additive/Fourier Wiener--Hopf
representation. It has not been converted into an admissible symbol of the
Mellin class needed by KKL2014TwoShifts, Theorem 3.3.

## 4. Localization

The cutoff is the paper's \(\chi_\infty\), placed on both sides of the
Wiener--Hopf operator. Any proposed Mellin transformation must account for
both multipliers and for commutator or remainder terms; dropping them is not a
formal simplification allowed by the AST.

## 5. Domain and codomain

The operator maps the branch-2 copy of
\(L^p(\mathbb R_+,w_\delta)\) to the branch-1 copy of the same weighted
half-line space, where
\(1<p<\infty\) and \(-1/p<\delta<1-1/p\). The two branch labels are material:
source results on a single unweighted half-line do not automatically typecheck
this block.

## 6. Relation with the transported block

The paper supplies a modulo-compact replacement of its transported
off-diagonal Cauchy block by \(W_{1,2}^{+}\). That statement is retained as
`CERTIFIED_MOD_COMPACT` only at the already supplied block-model boundary. It
does not imply that products such as \(R_1B_1W_{1,2}^{+}\) have a Mellin
symbol.

## 7. Comparison with \(R_y\) in KKL 2014

KKL2014TwoShifts, Lemmas 4.4--4.5, use

\[
 (R_yf)(t)=\frac{1}{\pi i}\int_0^\infty
 \left(\frac{t}{\tau}\right)^{1/y-1/p}
 \frac{f(\tau)}{\tau+t}\,d\tau,
\]

whose Mellin multiplier is
\(1/\sinh(\pi(\lambda+i/y))\). Compared with \(W_{1,2}^{+}\), this has a
different denominator, multiplicative power, representation, cutoff status,
weight bookkeeping, and branch typing. No exact or modulo-compact
identification was found. Lemmas 4.4--4.5 therefore remain
`UNPROVED_ADAPTATION`.

## 8. Comparison with Karlovich 2017 blocks

Karlovich2017HasemanSO, Theorem 5.6 (printed p. 479; PDF pp. 17--18), gives an
exact Mellin-PDO representation of shifted off-diagonal blocks on a
logarithmic star. It suggests a method, but the pure-cusp geometry, weighted
transport, and localized block in the present paper have not been shown to be
that operator. The comparison is `ANALOGY`.

## 9. Comparison with Karlovich 2025 generators

Karlovich2025Cusps, Lemmas 3.1--3.2 and Corollary 3.3, give localized pure-cusp
matrix models modulo compacts. They are the closest geometric evidence for the
paper's block. A paper-specific specialization is still required to identify
the exact cutoff, branch transport, weight, and algebra generator. The source
does not place the additional regularizer \(R_1\) in that algebra.

## 10. Mellin-PDO feasibility

A Mellin representation is plausible only as a research route, not as a
current result. Nothing verified permits declaring the localized
Wiener--Hopf operator itself to be a Mellin PDO or multiplying a hypothetical
symbol by \(r_1\).

## 11. Additional transformation required

A viable identification proof would have to:

1. conjugate the weighted branch operator through \(\Phi_\delta\);
2. express the additive localized kernel in logarithmic/Mellin coordinates;
3. retain both cutoff multipliers;
4. identify an explicit symbol class, including all variables and limits;
5. prove the required composition theorem or compact commutator estimates;
6. verify source and target branch types.

## 12. Conclusion

`NOT_IDENTIFIED`. The current evidence supplies neither
\(W_{1,2}^{+}=R_y\), nor \(W_{1,2}^{+}\simeq R_y\), nor membership of
\(W_{1,2}^{+}\) in a Mellin or 2025 cusp algebra. This blocks H2 and H3, but it
is not the only obstruction: the uniform H1 candidate is already blocked by
the \(R_1U_1\) interface.
