# P0-B: explicit assumptions and declared complex domains

P0-B adds semantic containers that retain mathematical conditions before any
future Mellin work.  These containers record declarations; they are not a
theorem prover and do not establish analytic or operatorial properties.

## Assumption contexts

`AssumptionContext` is an immutable, hashable, presentation-ordered collection
of explicit SymPy Boolean propositions.  Exact structural duplicates are
removed, and literal `True` is discarded because it adds no hypothesis. No
other hypothesis is discarded because it appears redundant.

Its consistency status has deliberately narrow meaning:

- `CONSISTENT`: the empty context, literal truth, or a system of elementary
  numeric bounds for symbols that the internal rules can certify;
- `INCONSISTENT`: a literal or elementary direct contradiction was found;
- `UNDETERMINED`: the limited rules cannot decide consistency.

In particular, `UNDETERMINED` does not mean consistent.  The implementation
never calls `bool()` on a symbolic relation and does not use failure by SymPy
to find a contradiction as evidence of consistency.

## Complex domains

`ComplexDomain` records a principal SymPy symbol, a conservative intersection
of elementary regions, an `AssumptionContext`, explicit exclusions, an
optional description, and separately retained external evidence.  Supported
regions are the complex plane, real line, upper and lower half-planes,
vertical and horizontal open strips, and regions defined by an explicit
Boolean condition.

Point membership and structural domain containment return `YES`, `NO`, or
`UNDETERMINED`.  Only results established by the limited internal rules are
reported as `YES` or `NO`.  Intersections add restrictions; they never widen a
domain.  Obvious empty intersections, such as the real line with an open
upper half-plane, are rejected explicitly.

The domain's evidence status reuses P0-A's `UNCERTIFIED` and
`EVIDENCE_SUPPLIED` distinction.  Supplying a domain or evidence does not
prove convergence, boundedness, compactness, or existence of an operator.

## Singular sets

`SingularSet` contains immutable `Singularity` records and conditioned
`SingularityAvoidance` records. Each singularity retains its location or set,
kind, order when known, relevant variable, conditions, origin, evidence, and
description. Origins distinguish user declarations, limited internal rules,
and external evidence. Evidence is recorded but never promoted to an
internally verified proof.

The internal detector is intentionally not universal. It has structural rules
only for:

- poles of `coth(pi*z)` at `z = I*n`, `n` integer;
- poles of negative powers of `sinh(pi*z)`;
- the specialization `z = lambda + I*kappa`, where a real-line pole requires
  `lambda = 0` and integer `kappa`;
- principal-branch sensitivity of logarithms and noninteger powers.

For an explicit numeric interval containing no integer, including
`0 < kappa < 1` and `-1 < kappa < 0`, the absence of an integer `kappa` is
stored as a limited conditioned avoidance rule. The underlying conditional
pole record remains; it is not erased. At integer values such as `kappa = -1`,
`kappa = 0`, or `kappa = 1`, `lambda = 0` is reported as singular. No numerical
sampling is used to declare a pole or branch absent. If the shifted principal
variable is not explicitly real, the detector retains the full shifted
integer lattice instead of applying the real-line specialization. A supplied
`ComplexDomain.real_line` is also explicit real-variable information, even
when the underlying SymPy symbol has no `real=True` property.

Negative powers also receive a conservative general exclusion at zeros of
their base when no more specific supported pole rule applies. This prevents a
symbolic cancellation such as `x/x` from erasing the required condition
`x != 0`; it does not classify arbitrary zeros or their orders.

For `gamma**(I*lambda)`, a positive-real-base treatment is available only when
the explicit context contains `gamma > 0`. A SymPy symbol property by itself
does not add that proposition to the context. The principal-branch records are
retained together with the conditioned positive-base avoidance.

## Conditional scalar identities

`ConditionalIdentity` is distinct from P0-A's `ExactIdentity`. It stores its
left and right scalar expressions, `AssumptionContext`, `ComplexDomain`,
`SingularSet`, scalar scope, verification status, evidence, and description.
It cannot be constructed with operatorial scope and `require_exact_identity`
rejects it.

Verification statuses mean:

- `DECLARED`: conditions and the equality are caller declarations;
- `SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS`: the concrete
  `sympy_simplify_difference` routine returned a residual whose `is_zero` is
  explicitly `True`, while every condition remained attached;
- `EVIDENCE_SUPPLIED`: an external evidence object was provided.

`SYMBOLICALLY_CHECKED_UNDER_ASSUMPTIONS` is not an operatorial theorem,
analytic continuation result, or proof of boundedness. Supplying evidence does
not automatically change a `DECLARED` identity's status and never turns
external evidence into an internal check.

`require_applicable_identity` accepts a conditional identity only when the
available context contains every required hypothesis, consistency is decided,
the available domain is provably no broader than the declared domain, and the
stored singularities are provably avoided. A missing condition, a broader or
unrelated domain, a singular point, or an undetermined check causes an
explicit exception. Callers can additionally request a concrete verification
status; external evidence does not satisfy a request for an internal symbolic
check.

Substitution updates expressions, assumptions, domain, and singularities;
values such as `kappa = 0` are rejected because they violate the stored open
interval. Scalar differentiation is available only after the original scalar
identity has passed the internal check and is applicable on its declared
domain; it preserves the domain and assumptions and unions newly detected
singularities with the old set. Combining identities unions contexts,
intersects domains, unions singular sets, retains evidence as separate
records, and never strengthens verification status.

Calling `.as_expr()` is a deliberate one-way projection to a SymPy equality.
The result has lost the semantic metadata and no API reconstructs a stronger
status from it.

## Reference family

`build_coth_conditional_identities(lambda, kappa)` constructs the scalar
functions

\[
s(\lambda)=\coth\!\bigl(\pi(\lambda+i\kappa)\bigr),\qquad
p_\pm(\lambda)=\frac{1\pm s(\lambda)}{2},
\]

and conditioned, scalar-only identities for

\[
p_-p_+=-\frac{1}{4\sinh^2\!\bigl(\pi(\lambda+i\kappa)\bigr)},
\]

\[
s'(\lambda)=-\frac{\pi}{\sinh^2\!\bigl(\pi(\lambda+i\kappa)\bigr)},
\qquad
p_\pm'(\lambda)=\mp\frac{\pi}{2\sinh^2\!\bigl(\pi(\lambda+i\kappa)\bigr)}.
\]

Every identity retains `lambda in R`, `0 < kappa < 1`, its limited pole
records, scalar scope, and the exact algebraic-check record.

## Current boundary

P0-B prepares semantic metadata for later Mellin work; it does not implement
Mellin transforms or convolution. It also does not implement weighted `L^p`
spaces, `V(R)`, `C_0(R)`, normalized shifts, total variation, general operator
matrices, determinants, Fredholmness, indices, compactness proofs, analytic
regularizers, boundary fibers, or a general inequality/theorem solver. SymPy
alone does not provide a complete theory of branches or poles.
