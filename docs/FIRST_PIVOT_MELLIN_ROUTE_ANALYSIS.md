# Mellin route for the first-pivot right core

## Outcome

```yaml
route: MELLIN
status: BLOCKED
first_uniform_blocker: R_1_times_U_1
secondary_blocker: W_1_2_plus_not_identified
final_symbol_generated: false
```

The analysis follows the stored order
\(Q_1R_1B_1W_{1,2}^{+}\). It never moves \(R_1\), commutes adjacent factors,
or invokes a symbol theorem before checking membership.

## F1. The product \(R_1B_1\)

### Common data

The paper defines

\[
 R_1=\Phi_\delta^{-1}\operatorname{Op}(r_1)\Phi_\delta,
 \qquad
 r_1\in\widetilde{\mathcal E}(\mathbb R_+,V(\mathbb R)).
\]

KKL2014TwoShifts, Theorem 3.3 (printed p. 942; PDF p. 8), gives
semiproduct closure modulo compacts when **both** symbols belong to
\(\mathcal E(\mathbb R_+,V(\mathbb R))\). Lemma 3.4 on the same page gives
an exact separated triple product only when all three source-class hypotheses
and the outer-variable restrictions hold.

### Case \(B_1=U_1\)

Under the paper's Mellin transport, the normalized constant dilation has
multiplier

\[
 d_{\gamma_1}(\lambda)=\gamma_1^{i\lambda}.
\]

For \(\gamma_1\ne1\), the paper explicitly records
\(d_{\gamma_1}\notin V(\mathbb R)\). Thus the second symbol hypothesis of
Theorem 3.3 fails in the verified calculus. Neither Theorem 3.3 nor Lemma 3.4
proves

\[
 R_1U_1\simeq\operatorname{Op}(r_1d_{\gamma_1}).
\]

That displayed formula is therefore a prohibited conclusion at this stage,
not a derived result. The missing step is an ordered composition theorem in an
explicitly enlarged calculus containing the oscillatory dilation multiplier.

Status: `BLOCKED`.

### Case \(B_1=\widehat G_1\)

Assuming the paper's coefficient hypothesis, \(\widehat G_1\) is a
multiplication operator with an \(SO(\mathbb R_+)\) coefficient. A
lambda-independent symbol is compatible in shape with the source Mellin
calculus. To apply Theorem 3.3 one must still record the conjugated operator,
variable dependence, and membership in the exact
\(\mathcal E(\mathbb R_+,V(\mathbb R))\) class used by the theorem.

Only after that specialization is proved would the theorem yield the ordered
semiproduct relation modulo compacts. No equality is claimed here.

Status: `SPECIALIZATION_TO_PROVE` (P-03 and the \(\widehat G_1\) part of
P-04).

### F1 conclusion

H1 is quantified over **both** values of \(B_1\). The plausible
\(\widehat G_1\) specialization does not repair the blocked \(U_1\) case.
Therefore H1 is not yet a defendible uniform lemma in a named verified class.

## F2. The product \(R_1B_1W_{1,2}^{+}\)

The F1 result is already incomplete. Independently,
\(W_{1,2}^{+}\) is defined as a localized Wiener--Hopf operator and has not
been proved to belong to the same Mellin class. Theorem 3.3 cannot be applied
to an unclassified right factor. Lemma 3.4 cannot be applied because no
three-symbol representation satisfying its source-class and variable
hypotheses exists.

KKL2014TwoShifts, Lemmas 4.4--4.5, do not solve this interface: their \(R_y\)
has not been identified with \(W_{1,2}^{+}\). Karlovich2017HasemanSO,
Theorem 5.6, is an analogy on a different geometry. Hence no
\(c_{B_1}^{+}\) is constructed.

Status: `BLOCKED` (P-01 and P-05).

## F3. The product \(Q_1R_1B_1W_{1,2}^{+}\)

The two exact exterior choices are

\[
 Q_1\in\{U_1^{-1}P^+,P^-\}.
\]

The paper gives exact Mellin symbols for \(P^\pm\), but that fact does not
authorize moving them through \(R_1\) or discarding the exterior dilation.
KKL2014TwoShifts, Lemma 2.7, is an unweighted commutation statement for SOS
shifts and the Cauchy algebra; its weighted transported specialization has not
been proved. Theorem 5.2 concerns different unweighted auxiliaries.

Since F2 is blocked and no weighted ordered exterior rule is available, no
\(c_{Q_1,B_1}^{+}\) is constructed.

Status: `BLOCKED` (P-06).

## Deferred criteria

KKL2014TwoShifts, Theorem 3.6; Karlovich2009MellinPDO, Theorem 4.3 in its
verified scalar scope; and KKL2014Regularization, Theorem 5.8 may be relevant
only after the appropriate Mellin membership, symbol, and nonvanishing
hypotheses have been proved. None is applied to the first pivot in Phase P.

## Mellin-route conclusion

No final symbol is generated. The first minimum research task on this route is
to define and prove ordered right-dilation closure for
\(\operatorname{Op}(r_1)U_1\) in a calculus that genuinely contains
\(\gamma_1^{i\lambda}\). The localization and class identification of
\(W_{1,2}^{+}\) remains the next independent obstruction.
