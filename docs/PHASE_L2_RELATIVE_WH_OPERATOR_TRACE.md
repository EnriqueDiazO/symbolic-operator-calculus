# Phase L2: relative Wiener--Hopf operator trace

## Objective and L1 provenance

This phase adds a typed, immutable end-to-end representation of the exact
relative-dilation Wiener--Hopf identity. The public factory
factor_dilation_conjugated_wiener_hopf is its single L1 construction boundary.
The resulting DilationConjugatedWienerHopfFactorization supplies all symbols,
kernels, variables, relative scale and reconstruction results. Fourier scaling
remains exclusively in scaled_convolution_kernel and scaled_fourier_symbol.

## Parametrized operators

DilationMap, DilationOperatorModel, WienerHopfOperatorModel and
OrderedRelativeOperatorProduct form a small domain model. Products store typed
factors in tuples and do not use SymPy Mul. Dynamic OperatorAtom names were not
used because they would encode indices, maps and scalar symbols as strings.

The original product has three ordered factors. Each factorization has two.
RelativeWienerHopfIdentity marks the relation as exact and has no
modulo-compact semantics.

## Traceability, rendering, and action

build_relative_wiener_hopf_trace retains the complete L1 result and places its
same symbol and kernel objects into the operator and action models. It validates
both kernel reconstructions and both symbol correspondences. The common action
is one half-line integral using the stored conjugated kernel.

render_relative_wiener_hopf_trace_latex accepts only the completed trace and
returns ten ordered steps. It calculates no symbols or kernels.

The notebook relative_wiener_hopf.ipynb builds both directions of the m=2 case
through public APIs and displays renderer output. It has no stored outputs.

## Tests and limitations

Tests cover immutability, typed order, L1 provenance, exactness, kernel
reconstruction, common action, symbol correspondence, both m=2 directions and
rendering. This is not a generic operator AST, arbitrary Wiener--Hopf evaluator,
matrix model, compact-equivalence engine or LaTeX parser. It covers one
off-diagonal pair and retains the assumptions enforced by L1.
