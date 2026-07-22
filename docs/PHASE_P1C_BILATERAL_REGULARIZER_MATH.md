# P1-C mathematics: the bilateral factored diagonal regularizer

## Data and order

Let

\[
 X=L^p(\mathbb R_+,w_\delta),\qquad
 B_1=\mathcal T_{1,-}\mathcal R_1Z_1^{-1}.
\]

All displayed products are ordered noncommutative products. The exact link is

\[
 Z_1^{-1}\mathcal A_{1,1}\mathcal T_{1,-}=\mathcal H_1,
 \tag{L}
\]

where \(\mathcal H_1\) is the P1-C alias for the article's \(\mathcal N_1\).
Assume the source-backed exact inverse identities for \(Z_1\) and
\(\mathcal T_{1,-}\), and the two independent core relations

\[
 \mathcal H_1\mathcal R_1-I\in\mathcal K(X),\qquad
 \mathcal R_1\mathcal H_1-I\in\mathcal K(X).
 \tag{R}
\]

## Product with \(\mathcal A_{1,1}\) on the left

Exact left multiplication of (L) by \(Z_1\), followed by exact cancellation,
gives

\[
 \mathcal A_{1,1}\mathcal T_{1,-}=Z_1\mathcal H_1.
\]

Therefore, by association only,

\[
\begin{aligned}
 \mathcal A_{1,1}B_1
 &=\mathcal A_{1,1}\mathcal T_{1,-}\mathcal R_1Z_1^{-1}\\
 &=Z_1\mathcal H_1\mathcal R_1Z_1^{-1},\\
 \mathcal A_{1,1}B_1-I
 &=Z_1(\mathcal H_1\mathcal R_1-I)Z_1^{-1}
 \in\mathcal K(X).
\end{aligned}
\]

This proof uses the \(\mathcal H_1\mathcal R_1\) relation and is recorded as
the `left_product` trace. It is not reused as the opposite product proof.

## Product with \(\mathcal A_{1,1}\) on the right

Exact right multiplication of (L) by \(\mathcal T_{1,-}^{-1}\), followed by
exact cancellation, gives

\[
 Z_1^{-1}\mathcal A_{1,1}
 =\mathcal H_1\mathcal T_{1,-}^{-1}.
\]

Hence

\[
\begin{aligned}
 B_1\mathcal A_{1,1}
 &=\mathcal T_{1,-}\mathcal R_1Z_1^{-1}\mathcal A_{1,1}\\
 &=\mathcal T_{1,-}\mathcal R_1\mathcal H_1
   \mathcal T_{1,-}^{-1},\\
 B_1\mathcal A_{1,1}-I
 &=\mathcal T_{1,-}(\mathcal R_1\mathcal H_1-I)
   \mathcal T_{1,-}^{-1}
 \in\mathcal K(X).
\end{aligned}
\]

This separate proof needs the exact inverse of \(\mathcal T_{1,-}\) and the
oppositely oriented \(\mathcal R_1\mathcal H_1\) relation.

## Compact-ideal transport

For a bilateral operator ideal \(\mathcal K(X)\), bounded \(L,R\in\mathcal
B(X)\), and \(X_0-Y_0\in\mathcal K(X)\),

\[
 LX_0R-LY_0R=L(X_0-Y_0)R\in\mathcal K(X).
\]

The first trace applies this with \((L,R)=(Z_1,Z_1^{-1})\); the second uses
\((L,R)=(\mathcal T_{1,-},\mathcal T_{1,-}^{-1})\). The implementation must
reject a non-bilateral ideal, an unbounded exterior factor, a space mismatch,
or differently labelled compact ideals.

## Certification boundary

Under the stated assumptions, both residues belong to the same
\(\mathcal K(X)\), so the main status is
`CERTIFIED_TWO_SIDED_FACTORED_REGULARIZER`.

The following do not follow from this calculation:

- membership of \(B_1\) in a common Mellin/Wiener--Hopf algebra;
- representation of \(B_1\) as one Mellin PDO;
- any symbol for the complete Schur correction;
- any cutoff replacement or new dilation cancellation in the Schur word.

The safe Schur integration consequently retains \(B_1\) as one structured
central object inside \(L_{2,1}B_1R_{1,2}\).

## Implemented certificate

`src/symbolic_operator_calculus/diagonal_regularizer.py` implements this
derivation without adding a second AST or semantic-relation hierarchy:

- `DiagonalLinkIdentity` accepts only the typed orientation
  `Z1_inverse A11 T1_minus = H1` and checks all five roles;
- `InvertibleAuxiliaryOperator` preserves the `LinearCombination` defining
  `T1_minus` and records both exact inverse products;
- `InvertibleMultiplicationOperator` records both exact multiplier products,
  requires nonvanishing evidence, and rejects zeros certified on
  \((0,\infty)\);
- `MellinCoreRegularizerCertificate` records the two core products separately
  and never promotes a unilateral relation;
- `propagate_compact_ideal` implements the bilateral-ideal step with explicit
  boundedness, space, and ideal checks;
- `FactoredTwoSidedRegularizerCertificate` retains both original products,
  two three-step traces, their semantic relations, evidence, assumptions,
  blockers, open obligations, and LaTeX;
- `SchurRegularizerInsertion` inserts the certified center as one opaque atom
  and fixes the common-algebra and single-PDO statuses at `UNPROVED` and the
  Schur symbol at `NOT_CERTIFIED`.

The public certification call is:

```python
result = certify_factored_two_sided_diagonal_regularizer(
    diagonal_operator=A11,
    left_interface=T1_minus,
    mellin_operator=H1,
    mellin_regularizer=R1,
    right_interface=Z1,
    right_interface_inverse=Z1_inverse,
    link_identity=link,
    assumptions=context,
)
```

Its strongest result is available only when both independently oriented core
relations survive every fail-closed check with one common bilateral ideal.
