# Fundamental hardness: ~2^256, no square-root speedup

## Answer to "is it 2^108? is there a sqrt speedup?"
**No to both.** Measured, not estimated:

### Effective search dimension n = 256 (not 41, not 108)
- 2366 boolean bits = 1522 gate outputs + 844 free inputs.
- Of the 844 free: **256 SELECTORS** (the codeword) + 588 inert don't-cares
  (toggling them changes zero equations).
- Full 256-toggle scan: **none is inert** — each breaks 6-131 equations. All live,
  all independent; **k = 0 bits removable** by linear/GF(2) elimination.
- agentA sets exactly 2 selectors {x_2081, x_24601} — a "2-of-256" near-solution.
- `x_15298 = OR(178 tree-1)·OR(78 tree-2)`, trees disjoint, union = 256. The 41
  "core-relevant"/load bits and all 4 residue bits live in tree-2; tree-1 has none.

### Naive exponent = 2^(n-k) = 2^256
The gap obstruction is controlled by only 2 selectors (a=x_2081, b=x_4287) → 4 channels,
3 distinct (M1,M2) values, none matching the response. The OR-trees produce ZERO extra
message choices at the gap. So the obstruction proves the neighborhood is dead but leaves
the full 256-bit codeword to invert.

### Why every square-root route fails
| Route | Verdict |
|---|---|
| Subset-sum MITM on the 41 additive load bits | Executable at 2^20.5 but **vacuous** — contribution matrix rank-1 mod p, bits are mod-p-INERT for the obstruction (verified: flipping any leaves G1,G2 unchanged mod p) |
| Two-tree product MITM `OR(178)·OR(78)` | Trees disjoint BUT tree-1 is pure activation (feeds no message/load wire); cross-tree triples are degree ≥3 entangled → **no clean additive interface to split on** |
| Reciprocal-lock baby-step/giant-step (r1·r2 ≡ 1 mod p, verified exactly) | The lock is a **tautology of the rank-1 coupling**: r2 is *defined* as 1/r1 by the 1-dim null space. The join condition is the LINEAR Q = G2 − r1·G1 = 0, not A(x)·B(y)≡1 over independent halves. No collision instance, nothing to search |
| ISD / syndrome decoding | Wrong tool — dense ℤ/p rank-1 certificate, no low-weight GF(2) target |
| Degree-2 MQ / bilinear MITM (2^128) | **Refuted: the verifier is degree > 2.** All 84 core-bit triples violate degree-2 |

### Why 2^108 specifically has no basis
2^108 would be a MITM over a ~216-item interface (216/2). No such interface exists:
k=0 (no linear reduction) and there is no degree-2 structure to interpolate. A 2^108/2^128
MITM *is* real arithmetic for the integer-load codeword — but that quantity is
**mod-p-inert**, so it cannot close the obstruction. Right cost, wrong problem.

### Structure: over-determined high-degree Boolean inversion
- Message LOAD is additive/subset-sum (the only linear structure) — deliberately
  **firewalled** from the verifier.
- VERIFIER is high-degree (>MQ), p-adically rigid. 21,814/39,033 equations are
  codeword-dependent (~85 constraints per selector bit) ⇒ **planted witness is
  essentially UNIQUE** (needle-in-haystack).
- Contrast: the PREVIOUS instance (small prime 2^31-1) fell to exactly a 2^(n/2)
  subset-sum-collision MITM because its matching interface was two clean subset-sums.
  Here that interface is blocked.

## The recursive-cascade confirmation (measured this session)
The obstruction is not a single wall but a **designed recursive cascade** — healing one
layer surfaces an isomorphic copy at the next, with the SAME residual value:
1. Layer 1 (gap): G1 = x_7068−x_2099−7376877·x_642, G2 = x_4432−x_19964−x_28730 → 11 fails.
2. Aligning the response wires → 16 leaf-ripple fails, reducing to
   L1 = `x_2964−x_26756−x_579` (x_579=p·x_19569), L2 = `9367949·(x_24548−x_25442)−x_7927`
   (x_7927=p·x_11052) — **identical residuals mod p to G1,G2** (617050…, 333101…).
3. Solving L1,L2 (x_2964:=x_26756, x_24548:=x_25442) → **12 fails**
   [4339,7740,8915,9549,10093,10882,13192,18937,23176,25225,28881,32535], reducing to
   `28669: x_23238−x_24292−10937191·x_33421` and `28671: x_36462−x_12985−x_34942`
   — the SAME algebraic form again.
Each layer's absorbers are side-effect-free p-multipliers; each heal is exact; the
obstruction simply relocates. This is the structural reason local repair cannot terminate.

## Bottom line
Fundamental hardness ≈ **2^256** needle-in-haystack inversion of an over-determined,
high-degree, p-adically-rigid Boolean circuit with an essentially unique witness.
**Not 2^108. No square-root speedup exists** — the one linear handle (additive load) is
firewalled from the verifier, and the verifier's degree >2 entanglement blocks every
collision interface.
