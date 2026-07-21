# P1-B left-core factorization

## Classification

`BLOCKED`

The requested left core is

\[
  \mathcal L_{2,1}=\mathcal B^-_{2,1}\mathcal T_{1,-}.
\]

Under the article's definition of `B^-_{2,1}` this expands exactly to

\[
 (\widetilde V_{\alpha_2}-\mathcal G_2)
 M_{\chi_\infty}W(b^-_{2,1})M_{\chi_\infty}
 (U_1^{-1}P^+ +P^-)(U_1^{-1}P^+ +P^-).
\]

The first auxiliary is already inside `B^-`; the second is introduced by the
requested interface.  The AST must retain both.  This is not the article's
left factor in `sections/06_shur_reduction.tex:382-390`.

## Exact information that may be used

The following transformations are exact when applied to a typed contiguous
subword with compatible domains:

- `V_gamma M_a V_gamma^-1=M_{a(gamma dot)}`;
- dilation/Wiener--Hopf covariance with frequency rescaling;
- `V_gamma M_chi V_gamma^-1=M_{chi(gamma dot)}`;
- Mellin-PDO radial covariance;
- cancellation of an adjacent typed inverse-dilation pair.

None turns a sum into a product.  In particular,
`T_{1,-}=U_1^-1 P^+ +P^-` is an ordered linear combination and the article
notes that the half-line `P^+` and `P^-` are not complementary projections
(`sections/05_haseman_operator_with_shift.tex:321-334`).  It is therefore
invalid to replace either auxiliary by `U_1^-1`, to pass `U_1^-1` across
`P^+`, or to declare the square idempotent.

The leftmost transported shift also expands as
`Vtilde_2=M_{rho_2}V_{gamma_2}`.  Its dilation scale is `gamma_2`, while the
auxiliary uses `gamma_1^-1`; combining them would produce a relative scale
only inside particular distributed summands and only after order-sensitive
covariance.  The two cutoffs of `W^-_21` then become scaled cutoffs and cannot
be replaced by the originals automatically.

## Why the target form is unavailable

The desired form

\[
  \mathcal L_{2,1}=L_{2,1}V_\gamma\Phi_\delta^{-1}
\]

does not occur as a contiguous factorization of the actual AST:

1. no `Phi_delta^-1` occurs in this core;
2. every candidate dilation is separated from the core boundary by `P^+`, or
   occurs only in one summand;
3. there is an extra full `T_{1,-}` in the literal requested word;
4. moving through the localized Wiener--Hopf block rescales both its WH
   symbol and its two cutoffs;
5. no verified cutoff-replacement theorem identifies the scaled cuts with
   the original ones in this mixed word.

## Minimal obligations

- `left_core_factorization`: state a corrected, nonduplicating interface and
  prove an exact or modulo-compact decomposition of every distributed summand.
- `cutoff_replacement_mod_compacts`: prove the required one- or two-sided
  compactness for each scaled cutoff difference in this precise word.
- `projection_dilation_interface`: prove the necessary order-preserving rule
  for `U_1^-1 P^+`; Mellin-convolution invariance alone does not authorize an
  unrecorded rewrite.
- `left_core_domain_compatibility`: retain the branch source/codomain and the
  weighted/Mellin transport boundaries.

No `L_{2,1}` is certified by the audited sources.
