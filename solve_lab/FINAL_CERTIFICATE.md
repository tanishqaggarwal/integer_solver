# Final structural certificate

## 1. The peeling certificate (replaces the erroneous "square matrix" claim)
The equation x atom matrix is **39033 x 38133 — NOT square** (the earlier squareness
claim was wrong). But the conclusion it supported is true and now proven directly:

Seeded by the 3 single-atom equations **18843, 19066, 20807**, the cascade
"an equation with exactly one not-yet-zero atom forces that atom to zero"
propagates to **all 38133 atoms**. Replayed and verified step by step; pivot
coefficients all lie in [1, 75].

Consequences:
- Full COLUMN rank 38133 over Q and over GF(q) for every prime q > 75.
- **Null space = {0}** -> no atom can be nonzero and cancel against others.
- Elimination only divides by integers <= 75, and the checker works over Z, so
  `M*a = 0, a in Z^38133  =>  a = 0`.
Therefore every atom must vanish exactly, including `p - x_26064` and `a2 - x_24453`.

## 2. Both pins are rigid
- `p - x_26064` (atom 36568) occurs in 13 equations: 8429, 11166, 11915, 12594,
  23869, 25313, 26785, 31400, 32300, 36106, 36767, 37257, 37666. They contain ONLY
  boolean atoms (x-x^2, x, 1-x) and copy atoms (x_a-x_b) — no large absorber.
- `a2 - x_24453` occurs in exactly ONE equation, 27494 — the only atom in the whole
  system with equation-multiplicity 1. Its 10 co-atoms are all boolean.
- Empirically: setting the p-class to 1 (or 0) breaks exactly the 13 pins;
  x_24453 := 0 breaks 5179 equations (a2 feeds 383 gadgets).

## 3. All 220 copy links are enforced by EQUATIONS, not by the harness
Rebuilt the class using only the 3558 pure `x_a - x_b` atoms extracted from
EQUATIONS.txt: the class of x_26064 has exactly 220 members with an explicit
219-edge spanning tree; every edge-atom occurs in >= 9 equations. All five slack
multipliers reduce to p by short explicit chains, e.g.
`x_28599 = x_9325 = x_1692 = x_26064 = p`. Nothing is merely assumed by forward-eval.

## 4. A real harness defect (premise vindicated, but the freedom is useless)
**164 of the harness's 32006 gates impose equalities the equations do NOT imply** —
each traces to an atom `x_out - x_in - x_slack` where the harness silently pins
`x_slack = 0`. So `H.forward()` really was searching a strict submanifold, and the
original off-manifold premise in OFFMANIFOLD.md was correct (its later "correction",
which cited the bogus square matrix, was wrong to withdraw it).
However the residual freedom is exactly **p-granular** (3707 absorber atoms
`x_s = p*x_f`), while at best_agentA exactly 2 of 38133 atoms are nonzero and BOTH
are nonzero mod p. p-granular slack shifts sums only by multiples of p, so it cannot
close either gap. Route closed.

## 5. The previous instance's loophole, exactly identified
The OLD instance
has the same 220-member copy class but it is **unpinned**, and the old witness assigns
**1 to all 220 members** (220/220). Every absorber becomes `1*slack`, so every
`target - accumulated - slack = 0` check is vacuously satisfiable. That is how the
previous instance was solved WITHOUT any cryptanalysis — its witness has 255 zeros and
a single 1 across the 256 selectors, and contains no occurrence of p or a2.

## 6. The patch is exactly the two pins
With constants abstracted, OLD vs NEW atom-shape multisets differ by exactly:
NEW-only `p - x_26064` and `a2 - x_24453` (plus one re-randomised blinding factor).
BIG constants (>2^100): **514 shared, 0 old-only, exactly 2 new-only**. Atoms 38133 vs
38131; equations 39033 vs 39031. The 512 leaf constants and the target T are identical,
so the hidden witness is the same in both instances.

## Bottom line
The hardness of this instance rests entirely on those two pins, and both are provably
unbreakable over Z. Every structural escape (atom cancellation, copy-class splitting,
under-constrained wires, integer-vs-mod-p slack, degeneracy, backward solve from the
target, weak division encodings) is closed with evidence. What remains is the genuine
256-bit accumulator-subset problem, of which the search sweeps covered ~6e-61.

Verified states: best_agentA_39022.json = 39022/39033 (record);
commonmode_39021.json = 39021 (both gaps closed exactly over Z, residue = 2 atoms);
nondegen_gapclosed_39013.json = 39013.
