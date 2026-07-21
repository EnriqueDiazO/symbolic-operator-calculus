# First-pivot closure before the localized Wiener–Hopf block

## Scope and hypotheses

This report stops immediately before \(W_{1,2}^{+}\). Assume the paper's
weighted model, \(\gamma_1>0\),

\[
 r_1,\widehat G_1,p^\pm
 \in\widetilde{\mathcal E}(\mathbb R_+,V(\mathbb R)),
\]

with \(\widehat G_1\) viewed as a lambda-independent symbol, and use the
already certified coefficient semiproduct. The exact relations are:

\[
 \Phi_\delta U_1^{\pm1}\Phi_\delta^{-1}=D_{\gamma_1^{\pm1}},
 \qquad
 \Phi_\delta P^\pm\Phi_\delta^{-1}=\operatorname{Op}(p^\pm),
 \qquad
 \Phi_\delta R_1\Phi_\delta^{-1}=\operatorname{Op}(r_1).
\]

No complementarity of \(P^\pm\) is used.

## Four unweighted reductions

Write

\[
 a_{+r,\gamma_1}(x,\lambda)
 :=p^+(\lambda)r_1(x/\gamma_1,\lambda)
\]

and

\[
 a_{+rG,\gamma_1}(x,\lambda)
 :=p^+(\lambda)r_1(x/\gamma_1,\lambda)
   \widehat G_1(x/\gamma_1).
\]

Then:

\[
\begin{aligned}
 L_{-+}
 &=P^-\operatorname{Op}(r_1)D_{\gamma_1}\\
 &\simeq
 \operatorname{Op}_{\mathrm{right}\text{-}\gamma_1}
 (p^-r_1,d_{\gamma_1}),\\[1mm]
 L_{++}
 &=D_{\gamma_1}^{-1}P^+\operatorname{Op}(r_1)D_{\gamma_1}\\
 &\simeq \operatorname{Op}(a_{+r,\gamma_1}),\\[1mm]
 L_{-\widehat G}
 &=P^-\operatorname{Op}(r_1)M_{\widehat G_1}\\
 &\simeq \operatorname{Op}(p^-r_1\widehat G_1),\\[1mm]
 L_{+\widehat G}
 &=D_{\gamma_1}^{-1}P^+\operatorname{Op}(r_1)M_{\widehat G_1}\\
 &\simeq
 \operatorname{Op}_{\mathrm{right}\text{-}\gamma_1^{-1}}
 (a_{+rG,\gamma_1},d_{\gamma_1^{-1}}).
\end{aligned}
\]

Every displayed \(\simeq\) is a certified equivalence modulo compact
operators. The covariance steps inside the derivations are exact.

## Weighted translation

Conjugation by \(\Phi_\delta^{-1}\) yields, in the original factor order,

\[
\begin{aligned}
 P^-R_1U_1
 &\simeq\Phi_\delta^{-1}
 \operatorname{Op}_{\mathrm{right}\text{-}\gamma_1}
 (p^-r_1,d_{\gamma_1})\Phi_\delta,\\
 U_1^{-1}P^+R_1U_1
 &\simeq\Phi_\delta^{-1}
 \operatorname{Op}(a_{+r,\gamma_1})\Phi_\delta,\\
 P^-R_1\widehat G_1
 &\simeq\Phi_\delta^{-1}
 \operatorname{Op}(p^-r_1\widehat G_1)\Phi_\delta,\\
 U_1^{-1}P^+R_1\widehat G_1
 &\simeq\Phi_\delta^{-1}
 \operatorname{Op}_{\mathrm{right}\text{-}\gamma_1^{-1}}
 (a_{+rG,\gamma_1},d_{\gamma_1^{-1}})\Phi_\delta.
\end{aligned}
\]

The AST stores these four left-hand sides unchanged and keeps one atomic
\(R_1\) in each.

## Logical types, symbol classes, and residues

| Core | Final relation | Output class | Compact residual |
|---|---|---|---|
| \(L_{-+}\) | `MOD_COMPACT_EQUIVALENCE` | right-factorized \((p^-r_1,d_{\gamma_1})\) | \(K_{-,r}D_{\gamma_1}\) |
| \(L_{++}\) | `MOD_COMPACT_EQUIVALENCE` | standard `E-tilde` | \(D_{\gamma_1}^{-1}K_{+,r}D_{\gamma_1}\) |
| \(L_{-\widehat G}\) | `MOD_COMPACT_EQUIVALENCE` | standard `E-tilde` | \(P^-K_{r,G}+K_{-,rG}\) |
| \(L_{+\widehat G}\) | `MOD_COMPACT_EQUIVALENCE` | right-factorized \((a_{+rG,\gamma_1},d_{\gamma_1^{-1}})\) | \(D_{\gamma_1}^{-1}(P^+K_{r,G}+K_{+,rG})\) |

Each residual is compact because bounded left/right multiplication preserves
the compact ideal, finite sums of compact operators remain compact, and
conjugation by the isometric isomorphism preserves compactness.

For the two factorized outputs, standard membership of the first component is
certified, but standard membership of the complete oscillatory product is not
asserted. For \(\gamma_1=1\), both oscillatory factors reduce to one and radial
scaling becomes the identity.

## Evidence and remaining blocker

The Cauchy semiproduct uses `KKL2014TwoShifts`, Theorem 3.3, printed p. 942,
PDF p. 8, checksum
`3f1b91576171681ff3b927685664c169134948dda515dec3737a72fc7a7ae1ec`.
Left covariance and radial-scaling invariance are direct proofs from the
paper's definitions.

P-04 and P-06 are discharged at the strengths above. P-01 and P-05 remain
blocked. The only structural factor common to all four unfinished words is
now their rightmost \(W_{1,2}^{+}\). This report neither identifies that
operator nor derives a right-composition rule, algebra membership, or any
Fredholm conclusion.
