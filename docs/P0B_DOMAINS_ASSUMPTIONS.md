# P0-B: explicit assumptions and declared complex domains

P0-B adds semantic containers that retain mathematical conditions before any
future Mellin work.  These containers record declarations; they are not a
theorem prover and do not establish analytic or operatorial properties.

## Assumption contexts

`AssumptionContext` is an immutable, hashable, presentation-ordered collection
of explicit SymPy Boolean propositions.  Exact structural duplicates are
removed, but no hypothesis is discarded because it appears redundant.

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

## Current boundary

This first P0-B layer does not implement Mellin transforms, branch analysis,
pole detection, conditional identities, weighted spaces, Fredholm theory, or
automatic inequality solving.  Singular sets and conditional scalar
identities are introduced only in the next guarded part of P0-B.
