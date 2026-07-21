# P1-B right-core factorization

## Classification

`BLOCKED`

The requested right core is

\[
  \mathcal R_{1,2}=Z_1^{-1}\mathcal B^+_{1,2}.
\]

But the article defines

\[
  \mathcal B^+_{1,2}
  =Z_1^{-1}(\widetilde V_{\alpha_1}-\mathcal G_1)W^+_{1,2}
\]

at `sections/06_shur_reduction.tex:391-399`.  The requested core therefore
has the exact AST

\[
  Z_1^{-1}Z_1^{-1}
  (\widetilde V_{\alpha_1}-\mathcal G_1)
  M_{\chi_\infty}W(b^+_{1,2})M_{\chi_\infty}.
\]

The repeated `Z_1^-1` is not the article's right factor and cannot be removed.

## Orientation and order

The transported shift is

\[
  \widetilde V_{\alpha_1}=M_{\rho_1}V_{\gamma_1}
  =Z_1U_1,
\]

with a forward, orientation-preserving scale `gamma_1>0`; see
`sections/05_haseman_operator_with_shift.tex:86-160,339-408`.  The target
right factor instead asks for `Phi_delta V_{gamma^-1}R_{1,2}`.  No inverse
dilation occurs in this core.  Expanding the difference produces one forward
dilation summand and one pure multiplication summand.  Cancelling one of the
two `Z_1^-1` factors against the `Z_1` inside only the first summand would
still require a distributive derivation and leaves the `-G_1` summand; it
cannot justify deleting the other leading `Z_1^-1`.

The two cutoff factors in `W^+_12` also remain.  Moving a dilation across this
localized block requires simultaneous rescaling of its multiplier cutoffs and
its Wiener--Hopf symbol.  The exact rescaling calculation in
`sections/05_haseman_operator_with_shift.tex:2172-2278` preserves the scaled
cuts; it does not replace them modulo compact operators.

## Why the target form is unavailable

There is no contiguous `Phi_delta V_{gamma^-1}` prefix in the AST.  Creating
one would require all of the following unproved operations: remove a repeated
normalization multiplier, extract an inverse dilation from a forward
transported shift, cross coefficient multipliers and a noncommutative sum,
and control the two cutoff replacements.  The left and right obstructions are
therefore not symmetric.

## Minimal obligations

- `right_core_factorization`: state a corrected, nonduplicating interface and
  derive its order-preserving summand representation.
- `right_inverse_dilation_orientation`: identify a genuine inverse relative
  dilation in the actual word; it cannot be introduced by notation.
- `cutoff_replacement_mod_compacts`: prove the exact side-specific compact
  difference needed after rescaling.
- `right_core_domain_compatibility`: verify the weighted and Mellin spaces at
  every `Phi_delta` boundary.

No `R_{1,2}` of the requested target form is certified.
