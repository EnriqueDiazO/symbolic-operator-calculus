# Cusp-algebra route for the first-pivot right core

## Outcome

```yaml
route: KARLOVICH_2025_CUSP_ALGEBRA
conclusion: NOT_SUPPORTED
membership_of_complete_words: OPEN
symbol_map_applied: false
fredholm_criterion_applied: false
```

`NOT_SUPPORTED` means that the verified source does not yet prove membership
of the four present words. It does not mean that a future specialization is
impossible.

## Factor-by-factor audit

| Factor or interface | Current fact | Missing fact for the 2025 route | Status |
|---|---|---|---|
| \(P^\pm\) | Exact Mellin convolution/Cauchy factors after \(\Phi_\delta\) | Match to the precise transported generator and weight conventions of the 2025 algebra | `SPECIALIZATION_TO_PROVE` |
| \(U_1\) | Exact normalized branchwise dilation; invertible | Membership as a generator or multiplier of the 2025 algebra | `NO_RULE` |
| \(\widehat G_1\) | Exact multiplication by an \(SO\) coefficient under the paper hypothesis | Match to the 2025 coefficient algebra on the transported branch | `SPECIALIZATION_TO_PROVE` |
| \(R_1\) | Mellin-PDO regularizer with \(r_1\in\widetilde{\mathcal E}\) | A theorem placing this variable-coefficient regularizer in the 2025 cusp algebra | `BLOCKED` |
| \(W_{1,2}^{+}\) | Exact localized Wiener--Hopf model; supplied block relation modulo compacts | Identification with a named 2025 off-diagonal generator, including cutoff, weight, and branches | `BLOCKED` |
| Branch 2 \(\to\) branch 1 | Typed in the Phase P records | Matrix placement and transported domain/codomain specialization in the source algebra | `SPECIALIZATION_TO_PROVE` |
| Complete ordered product | Exact AST word | Membership of every factor and closure under this ordered product | `BLOCKED` |

## Localization evidence

Karlovich2025Cusps, Lemmas 3.1--3.2 and Corollary 3.3 (article pp. 8--13),
prove compact localization remainders and explicit matrices for the source
pure-cusp model. These results are geometrically closer than the logarithmic
star results of 2017. Their use for the present block still requires an
explicit specialization of the cutoff, branch transport, kernel, and power
weight. They do not include the additional operator \(R_1\).

Accordingly, the already supplied off-diagonal block replacement in the paper
remains a relation modulo compacts, but it does not become a proof that
\(W_{1,2}^{+}\) or \(R_1B_1W_{1,2}^{+}\) belongs to the source algebra.

## Algebra, compact ideal, and symbol map

Karlovich2025Cusps, Theorem 5.1, provides a criterion inside its defined
weighted half-line mixed algebra. Theorems 6.1--6.2 identify source-generator
symbols and supply a homomorphism and Fredholm criterion for operators already
in the cusp algebra. These are conditional resources, not membership theorems
for arbitrary Mellin PDOs.

For the four right-core words, the following antecedents are still missing:

1. \(R_1\) is a member or a controlled multiplier of the 2025 algebra;
2. \(W_{1,2}^{+}\) is the required localized generator with correct branch
   typing;
3. the products with \(U_1\) and \(\widehat G_1\) remain inside the algebra;
4. every remainder created by the specialization lies in the compact ideal;
5. the complete ordered product has a source-defined symbol.

Until these are established, the homomorphism cannot be used to multiply
symbols, and no nonvanishing or Fredholm test is available.

## Conclusion

The cusp route is `NOT_SUPPORTED` at the current evidence boundary. Its main
potential advantage is that it treats cusp localization natively; its decisive
disadvantage is that no verified result places \(R_1\) in the algebra. This
route should be revisited only after either proving that inclusion or replacing
\(R_1\) by a source-admissible representative through a separately justified
modulo-compact statement.
