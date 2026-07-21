# P1-B Schur sign convention

## Three distinct objects

The implementation must distinguish:

1. `RawSchurCorrection`: the ordered product
   `A21 A11^(-1) A12`, with coefficient `+1` intrinsic to the raw product;
2. `SignedSchurContribution`: a separate sign, normally `-1`, applied exactly
   once to the raw correction in the algebraic pivot;
3. `SchurPivot`: `A22` plus that signed contribution.

Thus the structural algebraic definition is

\[
  A_{2,2}^{(1)}=A_{2,2}-
  \underbrace{A_{2,1}A_{1,1}^{(-1)}A_{1,2}}_{\text{raw correction}}.
\]

This is the definition in
`sections/06_shur_reduction.tex:281-305`.

## Article model convention

The normalized off-diagonal model has

\[
 A_{1,2}\simeq +(\widetilde V_{\alpha_1}-G_1)W^+_{1,2},
 \qquad
 A_{2,1}\simeq -(\widetilde V_{\alpha_2}-G_2)W^-_{2,1};
\]

see `sections/06_shur_reduction.tex:44-95`.  Substitution into the pivot's
outer minus gives a positive model term (`:310-348`).  The article names this
positive expression `C_{2,2}^{(1)}`:

\[
 C_{2,2}^{(1)}=
 (\widetilde V_{\alpha_2}-G_2)W^-_{2,1}A_{1,1}^{(-1)}
 (\widetilde V_{\alpha_1}-G_1)W^+_{1,2}.
\]

Accordingly,

\[
 A_{2,2}^{(1)}\simeq A_{2,2}+C_{2,2}^{(1)},
 \qquad
 A_{2,1}A_{1,1}^{(-1)}A_{1,2}\simeq-C_{2,2}^{(1)}.
\]

The exact equality later asserted for `C_{2,2}^{(1)}` concerns the chosen
model correction, while its relation to the original pivot remains modulo
compact operators (`sections/06_shur_reduction.tex:359-380`).

## Fail-closed rule

A raw correction may be constructed only with intrinsic coefficient `+1`.
A signed contribution accepts only the external signs `+1` or `-1` and must
reject a raw input already marked negative.  A pivot consumes exactly one
signed contribution.  This prevents the normalized minus in `A21` and the
Schur-complement minus from being silently applied twice.

No code or document in P1-B changes the article's naming convention.
