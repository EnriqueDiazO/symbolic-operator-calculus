# P1-D algebra audit: the mixed first-correction word

## Verdict

The exact correction assembly and membership in a Fredholm algebra are
separate questions. The literature audit classifies the actual word as

```text
BLOCKED_BY_COMMON_ALGEBRA_SOURCE_GAP
```

and its requested Fredholm symbol as

```text
BLOCKED_BY_ALGEBRA_MEMBERSHIP
```

This is a positive fail-closed result: exact operator composition is
certified, while the missing common-algebra theorem is represented explicitly
instead of being inferred from factorwise precedents.

## Factor coverage matrix

| Factor or feature | Best audited coverage | Classification for the complete word |
|---|---|---|
| \(\mathcal W^-_{2,1}\), \(\mathcal W^+_{1,2}\), with both cutoffs | Exact project identification; Karlovich 2025 half-line Wiener--Hopf algebra | `EXACT_FACTOR_MEMBERSHIP` for each concrete block |
| Diagonal Mellin core \(\mathcal N_1\) | Karlovich 2025 half-line algebra | `EXACT_FACTOR_MEMBERSHIP` for the core, not its chosen regularizer |
| \(\mathcal R_1=\Phi^{-1}\operatorname{Op}_M(r_{1,1})\Phi\) | KKL 2014 radial Mellin-PDO calculus | `COVERED_ONLY_IN_DIFFERENT_ALGEBRA` |
| \(\mathcal T_{1,-}=U_1^{-1}P^++P^-\) | article proves bounded exact invertibility | `BOUNDED_OPERATOR_ONLY`; no membership theorem for the target algebra |
| \(\widetilde V_{\alpha_k}=M_{\rho_k}V_{\gamma_k}\) | exact project formula and covariance | `DILATION_GENERATOR_NOT_COVERED` |
| Coefficient multipliers and projections | covered in several component calculi | `EXACT_FACTOR_MEMBERSHIP` only within their respective calculi |
| Ordered products mixing all rows above | no audited common source | `MIXED_PRODUCT_CLOSURE_NOT_COVERED` |
| Four retained cutoff multipliers | exact in the concrete localized blocks | no theorem allowing their deletion or transporting every factor into one algebra |

Factorwise membership statements cannot be unioned into membership in a
single algebra. The common certificate must name one algebra, one
representation on the same weighted space, a product-closure theorem, a
common compact ideal, and the relevant inverse-closed/Fredholm-symbol result.
No audited source provides that tuple.

## Audited sources and their boundaries

The following reports in the read-only repository
`/home/enriquedo/PersonalProjects/haseman-literature-index/reports` record the
source audit:

- `MIXED_WIENER_HOPF_MELLIN_DILATION_AUDIT.md`: no source simultaneously
  covers the localized Wiener--Hopf factors, a general radial Mellin PDO,
  nontrivial dilations, the ordered mixed products, compact ideal, inverse
  closedness, and a Fredholm symbol.
- `FIRST_PIVOT_WIENER_HOPF_BLOCK_IDENTIFICATION.md`: identifies the concrete
  doubly localized Wiener--Hopf blocks without deleting their cutoffs.
- `FIRST_PIVOT_REGULARIZER_IN_CUSP_ALGEBRA.md`: the chosen \(\mathcal R_1\)
  is not reduced to a Mellin convolution covered by the 2025 half-line
  algebra.
- `FIRST_PIVOT_PRODUCT_MEMBERSHIP_2025.md`: component coverage does not prove
  closure for the complete correction product.
- `DLS1995_CUTOFF_AND_LOCALIZATION_AUDIT.md`: the source does not authorize
  removal of the current double cutoffs and does not supply the needed
  half-line Wiener--Hopf/Mellin/dilation algebra.
- `MIXED_LOCAL_CUSP_ALGEBRA_ROUTE_DECISION.md`: the current case is the route
  requiring a new mixed algebra, not an application of an existing common
  one.
- `MINIMAL_MIXED_CUSP_LEMMA_REQUIRED.md`: even the minimal bridge requires a
  dilation-stable crossed-product algebra, product closure, a faithful orbit
  symbol, and a Fredholm theorem; these are unproved obligations.

Duduchava's commutative quotient framework does not silently solve the
dilation problem: a nontrivial dilation has a noncompact commutator with
appropriate multipliers, so it cannot simply be inserted into that
commutative quotient. The BKM inverse-closed and Fredholm results apply to
their precisely defined algebra, not to an unnamed extension containing the
present dilation and radial Mellin-PDO generators.

## Missing theorem fields

The runtime `CommonAlgebraCertificate` must preserve at least these blockers:

| Field | Actual P1-D value |
|---|---|
| `missing_generator` | general transported radial Mellin PDO and nontrivial dilation-containing interfaces in one algebra |
| `missing_closure_property` | closure under the exact seven-factor ordered product and its linear-combination subfactors |
| `missing_cutoff_theorem` | one theorem coordinating the retained double cutoffs with every other generator; no deletion is attempted |
| `missing_dilation_membership` | membership of \(V_{\gamma_k}\) and \(U_1^{-1}P^++P^-\) in the proposed common algebra |
| `missing_inverse_closedness` | inverse-closedness for the actual mixed extension, not merely for a component algebra |
| `missing_symbol_homomorphism` | a multiplicative symbol on the actual noncommutative mixed algebra plus a Fredholm criterion |

An algebra created as a Python label is not a mathematical extension theorem.
The implementation may certify an explicitly supplied synthetic algebra only
when every factor, closure property, ideal, and symbol rule is separately
evidenced. The article instance supplies no such certificate.

## Symbol gate

A correction symbol cannot be obtained by blindly multiplying factor symbols.
The derivation requires all of the following:

1. exact correction assembly;
2. certified membership of every factor in the same algebra or a sourced
   extension;
3. a common compact ideal and compatible representation;
4. certified treatment of every cutoff;
5. resolved relative-dilation effects;
6. a multiplicative rule for the ordered noncommutative product;
7. a factor symbol for every expanded factor.

The actual word fails item 2 first. Its symbol result must therefore contain
no candidate formula and use `BLOCKED_BY_ALGEBRA_MEMBERSHIP`. This does not
preclude later work proving the missing mixed-algebra lemma; it prevents P1-D
from presenting that future theorem as already available.

