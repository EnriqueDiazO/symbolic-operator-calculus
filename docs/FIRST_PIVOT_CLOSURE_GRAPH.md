# First-pivot closure graph

The nodes below are generated from the public Phase P API. An edge means that
the source obligation must be discharged before the target can close. Status
`SOURCE_VERIFIED` records supporting source evidence, and
`FORMALLY_REDUCED` records a symbolic reduction; neither means an analytic
proof. No node is `ANALYTICALLY_PROVED`.

```mermaid
flowchart TD
    P_01["P-01: BLOCKED"]
    P_02["P-02: FORMALLY_REDUCED"]
    P_03["P-03: SOURCE_VERIFIED"]
    P_04["P-04: BLOCKED"]
    P_02 --> P_04
    P_03 --> P_04
    P_05["P-05: BLOCKED"]
    P_01 --> P_05
    P_04 --> P_05
    P_06["P-06: BLOCKED"]
    P_05 --> P_06
    P_07["P-07: BLOCKED"]
    P_04 --> P_07
    P_05 --> P_07
    P_06 --> P_07
    P_08["P-08: BLOCKED"]
    P_04 --> P_08
    P_05 --> P_08
    P_06 --> P_08
    P_07 --> P_08
    P_09["P-09: BLOCKED"]
    P_06 --> P_09
    P_07 --> P_09
    P_08 --> P_09
    P_10["P-10: BLOCKED"]
    P_09 --> P_10
    P_11["P-11: BLOCKED"]
    P_10 --> P_11
    P_12["P-12: BLOCKED"]
    P_11 --> P_12
```

| ID | Statement | Status | Closure criterion |
|---|---|---|---|
| `P-01` | Identify Wplus_12 in a closure-compatible operator class. | `BLOCKED` | an exact or certified mod-compact class identification |
| `P-02` | Represent U1 in a calculus compatible with the R1 semiproduct. | `FORMALLY_REDUCED` | a proved ordered composition rule for Op(r1) followed by U1 |
| `P-03` | Represent Ghat1 as an admissible Mellin multiplication symbol. | `SOURCE_VERIFIED` | prove the lambda-independent Ghat1 symbol belongs to E(R+,V(R)) |
| `P-04` | Apply a valid ordered semiproduct rule to R1 B1 for both B1 choices. | `BLOCKED` | produce an identified symbol class and compact remainder for both B1 |
| `P-05` | Control (R1 B1)Wplus_12 without reordering factors. | `BLOCKED` | derive a certified model for the ordered three-factor suffix |
| `P-06` | Control Q1 before R1 B1 Wplus_12 without commuting it. | `BLOCKED` | preserve the stored order and produce a controlled remainder |
| `P-07` | Prove compactness of every accumulated residue. | `BLOCKED` | all residuals lie in K(L^p(R_+,w_delta)) with branch typing |
| `P-08` | Identify the resulting symbol and its precise class. | `BLOCKED` | one named symbol in a verified Mellin or cusp quotient class |
| `P-09` | Propagate the closure result to all 16 Phase O terms. | `BLOCKED` | all 16 terms receive a relation of identical certified strength |
| `P-10` | Reconstruct the complete correction C2 from the 16 controlled terms. | `BLOCKED` | recover C2 without changing the Phase N mod-compact relation |
| `P-11` | Prove membership of N2first in the selected operator class or algebra. | `BLOCKED` | a proved membership theorem with an identified quotient symbol |
| `P-12` | Only after membership, verify nonvanishing and Fredholmness conditions. | `BLOCKED` | all hypotheses of the chosen criterion are checked |

## Critical paths

- H1 depends on P-02, P-03, P-04, and P-07.
- H2 additionally depends on P-01 and P-05.
- H3 additionally depends on P-06.
- Propagation to the correction begins only at P-09.
- Membership and any later Fredholm test are strictly ordered as P-11 then
  P-12.
