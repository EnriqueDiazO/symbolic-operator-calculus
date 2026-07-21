# conditional draft — not ready for insertion into the article

The following is a conditional notation draft only.  It is not a proved
article lemma and must not be inserted while the listed interface, cutoff,
and algebra obligations remain open.

```latex
\paragraph{Audit of the full first-Schur word (conditional draft).}
The notation used in the present manuscript gives
\begin{align*}
 \mathcal B_{2,1}^{-}
 &=
 (\widetilde V_{\alpha_2}-\mathcal G_2)
 W_{2,1}^{-}\mathcal T_{1,-},\\
 \mathcal B_{1,2}^{+}
 &=
 Z_1^{-1}(\widetilde V_{\alpha_1}-\mathcal G_1)W_{1,2}^{+}.
\end{align*}
Consequently, the literal word
\[
 \mathcal B_{2,1}^{-}\mathcal T_{1,-}\mathcal R_1
 Z_1^{-1}\mathcal B_{1,2}^{+}
\]
contains two copies of both $\mathcal T_{1,-}$ and $Z_1^{-1}$ after
expansion.  It is not the manuscript identity
\[
 C_{2,2}^{(1)}
 =\mathcal B_{2,1}^{-}\mathcal R_1\mathcal B_{1,2}^{+}.
\]
No idempotence or cancellation of the repeated factors is assumed.

For a normalized dilation
$(V_\gamma f)(x)=\gamma^\nu f(\gamma x)$ one has exactly
\[
 V_\gamma M_\chi V_{\gamma^{-1}}
 =M_{\chi_\gamma},
 \qquad
 \chi_\gamma(x)=\chi(\gamma x).
\]
In particular, this identity does not imply
$M_{\chi_\gamma}\simeq M_\chi$.

Assume, as additional hypotheses not proved here, that there are
order-preserving representations in one weighted operator space
\begin{align*}
 \mathcal B_{2,1}^{\mathrm{prim}}\mathcal T_{1,-}
 &\mathrel{\diamond}
 L_{2,1}V_\gamma\Phi_\delta^{-1},\\
 Z_1^{-1}\mathcal B_{1,2}^{\mathrm{prim}}
 &\mathrel{\diamond}
 \Phi_\delta V_{\gamma^{-1}}R_{1,2},\\
 \mathcal R_1
 &=\Phi_\delta^{-1}\operatorname{Op}_M(r_{1,1})\Phi_\delta,
\end{align*}
where $\diamond$ is either equality or equality modulo one explicitly named
compact ideal.  Assume also every required replacement of a scaled cutoff is
proved on its stated side.  Then the middle dilation pair is genuinely
adjacent and the exact covariance identity gives the candidate reduction
\[
 L_{2,1}V_\gamma
 \Phi_\delta^{-1}\operatorname{Op}_M(r_{1,1})\Phi_\delta
 V_{\gamma^{-1}}R_{1,2}
 \mathrel{\diamond}
 L_{2,1}\Phi_\delta^{-1}
 \operatorname{Op}_M(r_{1,1}^{\gamma})
 \Phi_\delta R_{1,2},
\]
with
\[
 r_{1,1}^{\gamma}(x,\lambda)
 =r_{1,1}(\gamma x,\lambda).
\]

The preceding conditional calculation proves neither the required exterior
factorizations nor the cutoff replacements, mixed-algebra membership,
Fredholmness, or the symbol of the full Schur correction.
```

## Open labels before any article use

```text
left_core_factorization
right_core_factorization
cutoff_replacement_mod_compacts
regularizer_mellin_representation
wh_mellin_wh_sandwich_membership
fredholm_algebra_membership
schur_correction_symbol
```

Current article-level status:
`BLOCKED_BY_REGULARIZER_INTERFACE`.
