# Phase L1: relative dilation and the Wiener–Hopf symbol

## 1. Original problem

For `k != j`, start from the exact half-line operator

\[
\mathcal W_{k,j}=V_{\alpha_k}W(b_{k,j})V_{\alpha_j}^{-1}
\]

and the convolution kernel

\[
h_{k,j}(t)=-\frac{1}{\pi i}\frac{1}{t-i(c_k-c_j)}.
\]

The corresponding integral kernel is required to be

\[
\mathcal K_{k,j}(x,y)
=-\frac{\gamma_j}{\pi i}
\frac{1}{\gamma_kx-\gamma_jy-i(c_k-c_j)}.
\]

The question is whether the two outer dilations can be replaced by one
relative dilation and a Wiener–Hopf operator with a new symbol.

## 2. Definitions used

There is no weight in the operator considered here:

\[
(V_\alpha f)(x)=f(\alpha(x)),\qquad
\alpha_k(x)=\gamma_kx,\qquad \gamma_k>0.
\]

Thus this is `V_{alpha_k}`, not
`\widetilde V_{\alpha_k}=\rho_kV_{\alpha_k}`. On the half-line,

\[
(W(b)f)(x)=\int_0^\infty h(x-y)f(y)\,dy,
\qquad h=\mathcal F^{-1}b.
\]

Positive dilations preserve `(0,infinity)` exactly.

## 3. Fourier convention

The project uses

\[
(\mathcal Ff)(\lambda)=\int_{\mathbb R}f(t)e^{-it\lambda}\,dt,
\qquad
(\mathcal F^{-1}g)(t)=\frac1{2\pi}
\int_{\mathbb R}g(\lambda)e^{it\lambda}\,d\lambda.
\]

No SymPy transform convention is substituted for these definitions.

## 4. Kernel of the conjugated operator

Apply the three operators from right to left:

\[
\begin{aligned}
(V_{\alpha_k}W(b_{k,j})V_{\alpha_j}^{-1}f)(x)
&=\int_0^\infty
h_{k,j}(\gamma_kx-z)f(z/\gamma_j)\,dz.
\end{aligned}
\]

With `z=gamma_j y`, positivity preserves the limits and
`dz=gamma_j dy`. Hence

\[
(\mathcal W_{k,j}f)(x)=\int_0^\infty
\gamma_jh_{k,j}(\gamma_kx-\gamma_jy)f(y)\,dy.
\]

Substitution of `h_{k,j}` gives exactly the stated `mathcal K_{k,j}`.

## 5. Origin of the factor gamma_j

The prefactor `gamma_j` is the Jacobian of `z=gamma_j y`. It is not a
Fourier normalization and cannot be removed. There is exactly one such
factor.

## 6. Relative shift

For composition operators,

\[
V_\alpha V_\eta=V_{\eta\circ\alpha}.
\]

Consequently

\[
V_{\alpha_k}V_{\alpha_j}^{-1}
=V_{\alpha_j^{-1}\circ\alpha_k}.
\]

The relative map is therefore

\[
\beta_{k,j}=\alpha_j^{-1}\circ\alpha_k,
\qquad
\beta_{k,j}(x)=r_{k,j}x,
\qquad
r_{k,j}=\frac{\gamma_k}{\gamma_j}.
\]

The inverse belongs to `alpha_j`, not `alpha_k`. Scalar dilations commute as
point maps, but the displayed composition is fixed by the operator product.

## 7. Shift on the left

Factor

\[
\gamma_kx-\gamma_jy-i(c_k-c_j)
=\gamma_j\left(r_{k,j}x-y-i\frac{c_k-c_j}{\gamma_j}\right).
\]

Define

\[
h_{k,j}^{\mathrm L}(t)=\gamma_jh_{k,j}(\gamma_jt).
\]

Then

\[
\mathcal K_{k,j}(x,y)=h_{k,j}^{\mathrm L}(r_{k,j}x-y),
\]

so

\[
\mathcal W_{k,j}
=V_{\beta_{k,j}}W(b_{k,j}^{\mathrm L}).
\]

## 8. Left convolution kernel

Explicitly,

\[
h_{k,j}^{\mathrm L}(t)
=-\frac{\gamma_j}{\pi i}
\frac1{\gamma_jt-i(c_k-c_j)}.
\]

There is no `gamma_j^2`.

## 9. Left symbol

The Fourier scaling law below gives

\[
b_{k,j}^{\mathrm L}(\lambda)
=b_{k,j}\!\left(\frac{\lambda}{\gamma_j}\right).
\]

## 10. Shift on the right

Define

\[
h_{k,j}^{\mathrm R}(t)=\gamma_kh_{k,j}(\gamma_kt).
\]

For `r=r_{k,j}`, changing variables `y=rz` in
`W(b^{\mathrm R})V_\beta f` gives

\[
\begin{aligned}
(W(b^{\mathrm R})V_\beta f)(x)
&=\int_0^\infty h^{\mathrm R}(x-z)f(rz)\,dz\\
&=\int_0^\infty \frac1r
h^{\mathrm R}\!\left(x-\frac yr\right)f(y)\,dy.
\end{aligned}
\]

Since `(1/r)gamma_k=gamma_j` and `gamma_k/r=gamma_j`, this kernel is
`gamma_j h(gamma_k x-gamma_j y)`. Therefore

\[
\mathcal W_{k,j}=W(b_{k,j}^{\mathrm R})V_{\beta_{k,j}}.
\]

## 11. Right convolution kernel and symbol

The kernel and symbol are

\[
h_{k,j}^{\mathrm R}(t)=\gamma_kh_{k,j}(\gamma_kt),
\qquad
b_{k,j}^{\mathrm R}(\lambda)
=b_{k,j}\!\left(\frac{\lambda}{\gamma_k}\right).
\]

## 12. Fourier scaling law

Let `a>0` and `h_a(t)=a h(at)`. With `u=at`, `dt=du/a`:

\[
\begin{aligned}
(\mathcal Fh_a)(\lambda)
&=\int_{\mathbb R}a h(at)e^{-it\lambda}\,dt\\
&=\int_{\mathbb R}h(u)e^{-iu\lambda/a}\,du\\
&=(\mathcal Fh)(\lambda/a).
\end{aligned}
\]

This proves both symbol formulas with the project's forward exponential and
introduces no `2 pi` factor.

## 13. Independent operator verification

The conjugation identity following from the same scaling law is

\[
V_aW(b)V_a^{-1}=W\!\left(b(\cdot/a)\right).
\]

Using `V_{\alpha_k}=V_{\beta_{k,j}}V_{\alpha_j}` gives the left
factorization. Inserting
`V_{\alpha_k}^{-1}V_{\alpha_k}` around `W(b)` gives the right one. These
operator manipulations reproduce the relative scale and symbols found from
the kernels.

## 14. Equality of the two forms

Both reconstructed kernels simplify identically to

\[
\gamma_jh_{k,j}(\gamma_kx-\gamma_jy).
\]

Therefore

\[
V_{\beta_{k,j}}W(b_{k,j}^{\mathrm L})
=W(b_{k,j}^{\mathrm R})V_{\beta_{k,j}}
=\mathcal W_{k,j}.
\]

## 15. General cases

Because `gamma_k,gamma_j>0`, half-line indicators satisfy
`chi_±(lambda/gamma)=chi_±(lambda)`.

If `k<j`, then `Delta c=c_k-c_j<0` and

\[
\begin{aligned}
b_{k,j}^{\mathrm L}(\lambda)
&=2\chi_+(\lambda)
e^{(c_k-c_j)\lambda/\gamma_j},\\
b_{k,j}^{\mathrm R}(\lambda)
&=2\chi_+(\lambda)
e^{(c_k-c_j)\lambda/\gamma_k}.
\end{aligned}
\]

If `k>j`, then `Delta c=c_k-c_j>0` and

\[
\begin{aligned}
b_{k,j}^{\mathrm L}(\lambda)
&=-2\chi_-(\lambda)
e^{(c_k-c_j)\lambda/\gamma_j},\\
b_{k,j}^{\mathrm R}(\lambda)
&=-2\chi_-(\lambda)
e^{(c_k-c_j)\lambda/\gamma_k}.
\end{aligned}
\]

## 16. The four m=2 symbols

Let `d=c_2-c_1>0`.

\[
\begin{array}{c|c|c}
&\text{shift left}&\text{shift right}\\ \hline
\mathcal W_{1,2}
&2\chi_+(\lambda)e^{-d\lambda/\gamma_2}
&2\chi_+(\lambda)e^{-d\lambda/\gamma_1}\\[0.4em]
\mathcal W_{2,1}
&-2\chi_-(\lambda)e^{d\lambda/\gamma_1}
&-2\chi_-(\lambda)e^{d\lambda/\gamma_2}
\end{array}
\]

These are the original article symbols. No block factor `1/2` or normalized
K1/K2 symbol has been applied.

## 17. Exact versus modulo compact operators

Under the stated starting identity
`mathcal W_{k,j}=V_{alpha_k}W(b_{k,j})V_{alpha_j}^{-1}`, both
factorizations are exact equalities. No `chi_infinity`, block approximation,
or compact remainder occurs in L1. If a later article step localizes or
replaces this exact operator, `simeq` begins at that later step, not here.

## 18. Symbolic and numerical verification

The implementation records the original, left, and right kernels and symbols
in immutable metadata. SymPy simplification verifies both reconstructed
kernels against the original conjugated kernel. Tests also verify the known
inverse transform pairs `2 kplus=h12` and `-2 kminus=h21`, the positive-scale
law, indicator stability, all four m=2 symbols, and absence of spurious
normalizations.

A numerical check uses
`gamma_k=3/2`, `gamma_j=5/3`, `c_k=1`, `c_j=4`, `x=7/5`, and `y=4/5`.
The original, left-reconstructed, and right-reconstructed complex kernel
values agree within `1e-13`.

## 19. Limitations

The implementation is specific to positive scalar dilations and the concrete
off-diagonal family. It is not a general kernel-factorization CAS, does not
add a generic Wiener–Hopf operator to the noncommutative AST, and does not
modify normalized blocks, Schur reduction, rendering, or export.

## 20. Recommended article continuation

Immediately after the displayed formula for `mathcal K_{k,j}`, insert the
following verified statement:

```latex
Set
\[
\beta_{k,j}=\alpha_j^{-1}\circ\alpha_k,
\qquad
\beta_{k,j}(x)=\frac{\gamma_k}{\gamma_j}x.
\]
Then the operator admits the two exact factorizations
\[
\mathcal W_{k,j}
=V_{\beta_{k,j}}W\!\left(b_{k,j}^{\mathrm L}\right)
=W\!\left(b_{k,j}^{\mathrm R}\right)V_{\beta_{k,j}},
\]
where
\[
b_{k,j}^{\mathrm L}(\lambda)
=b_{k,j}\!\left(\frac{\lambda}{\gamma_j}\right),
\qquad
b_{k,j}^{\mathrm R}(\lambda)
=b_{k,j}\!\left(\frac{\lambda}{\gamma_k}\right).
\]
Indeed, with
\(h_{k,j}=\mathcal F^{-1}b_{k,j}\), the two convolution kernels are
\(h_{k,j}^{\mathrm L}(t)=\gamma_jh_{k,j}(\gamma_jt)\) and
\(h_{k,j}^{\mathrm R}(t)=\gamma_kh_{k,j}(\gamma_kt)\); both forms have
the integral kernel
\(\gamma_jh_{k,j}(\gamma_kx-\gamma_jy)\).
```
