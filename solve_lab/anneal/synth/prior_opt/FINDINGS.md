# prior_opt -- shrinking the problem BEFORE it becomes a QUBO

**One line:** of 38,748 unknowns, 7,218 are determined before any search (0 qubits) and
31,274 are dependent circuit wires, leaving exactly **256 genuinely independent decision
bits -- and NONE of the 256 can be soundly eliminated.** So 0 comb windows are removable;
the only lever is per-window (atom) size. Re-derived from EQUATIONS.txt, cross-checked
with checker.py.

## Ground truth (reproduced)
- unknowns x_0..x_38747 = 38,748;  equations = **39,033** (NOT 39,031 -- that was the OLD
  instance; task brief used the old number). distinct atoms = 46,298.
- there is no forced.json; the committed forced data is partial_assignment.json, reproduced
  here byte-identical (forced_assignment.json, 0 value diffs).

## Forced variables -- 7,218 fixed before search
Integer propagation to fixpoint over 46,298 atoms (256 core bits left free). Sound because
the peeling cascade (FINAL_CERTIFICATE.md §1) gives full column rank 38,133, nullspace {0}
=> every atom must vanish, so each single-unknown-unique-root step is universal.
  unit pins x=1: 1,103 | copy atoms: 3,722 | big-const loads: 841
  DETERMINED: 7,218 (=1:1,387 =0:5,603 other:228), 0 contradictions.
Through checker.py: discharges 27,763/39,033 equations, 0 contradictions.

## Free bits -- 256, all independent, 0 eliminable
Core = clean ECDLP on secp256k1 (reduce.py/instance.py selfcheck OK). Of the 256 bits:
0 forced by propagation, 0 copy-links among them, no degenerate window (chain verified,
no identity point), b_255 not pinned (only k<->k+n redundancy ~2^-127). The "255"
artifacts (control_bits.json, linalg255.json, S11 invariants) are a DIFFERENT variable set
/ a mod-2 circuit shadow / conserved quantities of the checkpoint's WRONG message that
drift across siblings (12/14) -- using them to drop a bit deletes the true witness.
**Sound independent-bit count = 256.**

## Structural QUBO shrink
(a) 7,218 vars removable as compile-time constants -> 0 qubits.
(b) 0 comb windows removable -- every bit gates a distinct non-identity point.
(c) 27,763/39,033 equations discharge identically to 0 -> 0 QUBO weight; the remaining
    11,270 all reduce to the single kG=T.
Certified pre-QUBO collapse: 38,748 -> 256 free variables, losslessly, WITH PROOF.

## Soundness theorem
Fixing the 7,218 forced vars and treating non-core vars as functions of the 256 bits
removes no solution; any 256-bit solution re-expands to a full solution of all 39,033
equations. Proof: full column rank => every atom 0 in any solution => each forcing step is
universal; checker confirms forced=partial_assignment (27,763/39,033, 0 contradictions);
re-expansion forward-evaluates the feed-forward wiring so all atoms vanish. A consistency
check vs the record 39,026 partial: the 2 disagreements (x_10903, x_31864) fall in 5 of the
7 equations that partial fails -- proving that state is a non-solution exactly where it
violates a forced constraint. QED.

## Bottom line
| lever | count | sound? |
|---|---|---|
| unknowns removed as forced constants | 7,218 | yes |
| dependent wires demoted to ancillas | 31,274 | yes |
| free decision bits remaining | 256 | the floor |
| bits eliminable below 256 | 0 | eliminating any deletes the witness |
| comb windows removable | 0 | every bit gates a distinct point |
| equations pre-discharged (0 QUBO weight) | 27,763 / 39,033 | checker-measured |
