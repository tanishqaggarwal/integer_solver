# Agent H log — DECOMPOSITION angle

## Step 0 (start)
- Read PROMPT.txt, RESUME.md (1128 lines), STATE.json.
- Verified claimed best: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
  -> `satisfied 39026/39033 (7 failing)`, failing `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`. CONFIRMED.
- Plan: build independent atom/gate model, then graph decomposition (components, articulation,
  separators, treewidth estimate) and ask whether the residual admits a small separator.

## Step 1 — the frame
Built an independent parse (`model.py`, 42,267 atoms). Orienting every atom of syntactic form
`x_t - rest` as a definition gives an ACYCLIC gate DAG over all 38,748 vars: 8,747 free inputs,
30,001 defined, 12,266 checks, ZERO cycles (prior lab frame: 7,273 free / 1,800 in cycles).
Forward eval from all-zero free inputs = 39,005, only 5 nonzero check atoms (checker-verified).

## Step 2 — decomposition angle, REFUTED
eq-var graph: 1 component. atom-eq graph: 1 giant + 3,234 singletons. Free-input hypergraph
(hyperedge = equation): 1 component, 8,747 vars / 18,248 equations; 20,785 equations are
identically satisfied by forward eval. Residual closure = 6,007 free inputs / 9,244 equations.
=> NO separator, NO block decomposition. My assigned hypothesis is false. Logged and pivoted.

## Step 3 — decompilation
Residual = 3 conditions: OR(256 bits)=1, x_37892 = C1 (mod p), x_13682 = C2 (mod p).
The circuit is a binary MUX tree over 256 bits; at a node where both children fire the drivers
(X3,Y3) are UNPINNED (gate 1-L*R) but three rank-2 checks force the chord identities.
=> tree computes an elliptic-curve multi-addition. Curve (shifted x = X + K/3):
y^2 = x^3 + B, p = 2^256-2^32-977, order = secp256k1's n (PRIME). All 256 leaf points + target
on curve. 255/256 points have their double present: ONE doubling chain of length 256 from bit
x_2779, so P_i = 2^i G and full solve <=> S = binary expansion of DL_G(P*).

## Step 4 — CERTIFICATION of step 3 (the gap I flagged)
`verify3.py`: setting X3,Y3 := P_b1 + P_b2 makes the chord residuals x_25614, x_34220 vanish
mod p in 56/56 tested (u-bit, w-bit) pairs; perturbing X3 by 1 makes them nonzero in 56/56.
Over Z the triple needs one extra lift: a28438 demands 2264251 | (x_15286/p) (2264251 = 11*43*4787).
Solving t mod 2264251 algebraically closes all three atoms EXACTLY over Z, 3/3 cases.
DECOMPILATION CERTIFIED.

## Step 5 — weak-DL search (Priority 1). ALL NEGATIVE.
| family | extent | result |
|---|---|---|
| single bit | all 256 | none |
| weight 2,3,4 | exhaustive over 256 points | none |
| weight <= 6 | exhaustive MITM (2,796,416 subsets of size<=3 both sides) | none |
| k = 2^a(2^m-1) (one run) | all 32,896 | none |
| two runs of consecutive bits | exhaustive 33k x 33k via hash | see fam.out |
| periodic bit patterns, period d<=16, all 2^d residue sets | exhaustive | see fam.out |
| k < 2^44 | BSGS 2^22 x 2^22 | none (532 s) |
| k = 2^a * s, s < 2^22 | all 256 shifts | none |
| k = N - c, c < 2^44 | BSGS on -T | see fam.out |
| weight >= 250 (complement weight <= 6) | exhaustive MITM | see mitm_c.out |
