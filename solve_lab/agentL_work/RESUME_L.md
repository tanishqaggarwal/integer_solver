# RESUME_L — agent L.  Angle: the same-OR-group double-leaf case (F's one open case).

## 0. Baseline re-verified BY ME
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
-> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
Nothing of mine beats it yet.

## 1. THE OPEN CASE IS SETTLED — NO SUM.  (exp1.py, exact re-propagation)
Test node **x23153**, children = free leaves **x19326** and **x28825** (both pin their own wire pairs).
Measured mux output wires (x35346, x10824) by exact forward propagation:

| leaf bits | sel_a | sel_b | sel_ab | slot wire holds |
|---|---|---|---|---|
| a=1,b=0 | 1 | 0 | 0 | **A** = (C10424, C27436) mod p exactly |
| a=0,b=1 | 0 | 1 | 0 | **B** = (C20930, C30632) mod p exactly |
| a=1,b=1 (SAME OR-GROUP) | 0 | 0 | 1 | **chordK(A,B)** exactly — NOT A+B, NOT A, NOT B |

All 10 local residual atoms vanish mod p in all three cases.  The 3 stage checks have a unique
solution for the chord-output pair and the 3rd check is consistent (over-determined, satisfied).

**Why there is no sum** (the mechanism, measured not assumed):
* each free leaf L has 3 residual atoms: `L*(L-1)` (boolean) and two pins `L*(w - C) - p*h`.
* the SAME wires w are the mux value wires; their off-guard is **`(1 - L)*w = p*h`** — guarded by
  the LEAF BIT, not by the selector.  So when two leaves in one group are ON, both pins DO fire and
  both value wires DO hold their constants.
* but the mux is `out = sel_a*va + sel_b*vb + sel_ab*vab` with
  `sel_a = a(1-b)`, `sel_b = b(1-a)`, `sel_ab = ab` — so both firing pins are multiplied by **0**.
  Exactly one term survives, always.  The sum is a sum of three terms of which two are annihilated.
=> **the fold law is unchanged; the tree model stands; F's caveat is refuted, not confirmed.**

## 2. CORRECTION TO THE CONFIGURATION COUNT: 2^178 - 1, not 2^256 - 1
Full OR-tree rebuilt from scratch (ortree2.py/census.py): **two perfect binary OR-trees of depth 7**
rooted at x8599 and x21839, 127 internal nodes each = **254 internal nodes, 256 leaves**.
Of the 256 leaves: **178 are FREE booleans, 78 are DEFINED AS THE CONSTANT 0** (dead leaves, their
value wires are also literal 0).  178+78 is exactly the split F read as the root's "178|78 leaf
support" and E read as "178 | 41+21+16 channels" — it is the LIVE/DEAD leaf split.
All 254 nodes have exactly the (2,2,2) three-way, two-coordinate mux (253 exact, 1 with a 3rd
ab-term).  => **reachable configurations = 2^178 - 1 non-empty live-leaf subsets.**

## 3. Other measured facts (my own knob set = ALL 8,747 free vars)
* **Every handle var appears in exactly ONE atom.**  3,681 handle-only free vars (their variation
  changes every atom by a multiple of p), 5,066 value vars, 0 unclassified.  handles2.py.
* |Z| = 7,202 wires that are 0 mod p for every assignment.
* The **deliverable has exactly ONE leaf ON: x24601.**  It is the trivial single-leaf pass-through.
  Its 4 nonzero atoms are the root both-live guards `(1-sel_root)*x8731`, `(1-sel_root)*x9118` and
  the two root slot links `x4432-x19964`, `x7068-x2099`.  i.e. the deliverable pays 7 equations for
  *not* firing the root stage.
* Coordinate orientation at x23153: **mux coord2 is the x-coordinate, coord1 is y.**

## 4. Files (all in agentL_work/)
`trace.py ortree.py ortree2.py census.py wire.py link.py crux.py onset.py fail7.py handles2.py
exp1.py` + pickles `ortree2.pkl nodes.pkl outwires.pkl handles.pkl`.

## 5. NEXT
Build the complete fold evaluator over all 254 nodes + 178 leaf constants (F only had 47/72 stages
wired; the OR-tree route wires 254/254).  Then extract the target and invert.
