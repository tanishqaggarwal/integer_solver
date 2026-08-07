# RESUME_L — agent L.  Angle: the same-OR-group double-leaf case (F's one open case).
Integers, congruences and polynomials only.  No curve/group framing anywhere.

--------------------------------------------------------------------------------------------
## 0. SCORES
* Shared baseline **39,026 / 39,033**, re-verified BY ME:
  `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
  -> `satisfied 39026/39033 (7 failing)`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
* My own best: **39,018 / 39,033** = `agentL_work/assign_L1.json`, checker-verified
  (`satisfied 39018/39033 (15 failing)`).  **Below baseline — do not adopt it as the deliverable.**
  It is not a search result; it is the *canonical* assignment my constructor produces for one
  chosen leaf, and its value is that it leaves exactly the two target congruences nonzero and
  nothing else (see §3).  Its cost of 15 is simply what those two atoms happen to price at.
* **No infeasibility is claimed anywhere in this file.**

--------------------------------------------------------------------------------------------
## 1. THE ASSIGNED QUESTION IS SETTLED — THERE IS NO SUM.  (`exp1.py`, exact re-propagation)
F's open case: "two leaves ON in the SAME OR-group fire two pins on different wires and the
slot may then see a SUM rather than a single value."  **Refuted, by measurement.**

Test node **x23153**, children = free leaves **x19326** and **x28825**, each pinning its own
wire pair.  Mux output wires (x35346, x10824) measured by exact forward propagation:

| leaf bits | sel_a | sel_b | sel_ab | slot wire actually holds |
|---|---|---|---|---|
| a=1,b=0 | 1 | 0 | 0 | **A** = (C10424, C27436) mod p, exactly |
| a=0,b=1 | 0 | 1 | 0 | **B** = (C20930, C30632) mod p, exactly |
| a=1,b=1 **(same OR-group)** | 0 | 0 | 1 | **chordK(A,B)** exactly — not A+B, not A, not B |

All 10 local residual atoms vanish mod p in all three cases; the 3 stage checks have a unique
solution for the chord-output pair and the over-determined 3rd check is satisfied.

**Mechanism (measured, not assumed).**  Both pins DO fire — that part of F's worry is correct.
Each live leaf L carries `L*(L-1)` plus two pins `L*(w - C) - M*p*h`, and those same wires w are
the mux value wires.  Their off-guard is **`(1-L)*w = p*h`, guarded by the LEAF BIT, not by the
selector**.  So with two leaves of a group ON, both value wires hold their constants.  But the
mux is `out = sel_a*va + sel_b*vb + sel_ab*vab` with `sel_a=a(1-b)`, `sel_b=b(1-a)`,
`sel_ab=ab`, so both firing pins are multiplied by **0**.  The slot is a sum of three terms of
which two are always annihilated.  **The fold law is unchanged; the tree model stands.**

**Generalised beyond the one node.**  The single-node measurement above is the clean isolated
case, but the general statement is verified by §3: full assignments built for ON-sets of size
1, 2, 5, 73, 200 and 256 leaves — which necessarily contain many same-OR-group multi-leaf
configurations at every depth — each collapse to exactly the same two residual atoms and
nothing else.  If any same-group case produced a sum, those constructions would break.

**Not the same phenomenon as E's saturation.**  E's "a channel contributes at most once however
many of its bits are on" is the *pass-through* branch (exactly one child live -> value passes
through unchanged).  The double-leaf case is the *other* branch: both children live -> the
value is chordK(A,B), which differs from both A and B.  Two distinct branches of the same mux.

--------------------------------------------------------------------------------------------
## 2. F'S 2^256 - 1 COUNT IS CORRECT.  (I briefly had 2^178; that was wrong and is retracted.)
`global.py` enumerates **every** OR node in the circuit (383 of them) and finds the unique one
that is nobody's child: **root = x9274**.  Under it: **384 leaves — 256 free booleans and 128
literal `:= 0`** (verified: every dead leaf's definition is exactly the constant 0).  My first
pass built the OR tree upward from x8599/x21839 only, which is a *sub*-forest (254 nodes,
178 live / 78 dead leaves) — hence the wrong count.  **Reachable configurations = 2^256 - 1.**

Note where the lab's recurring "178 | 78" actually comes from: it is the **live-leaf split at
the root**, `x1711` (178 live) vs `x17215` (78 live).  Same number F read as root leaf support
and E read as channels 178 | 41+21+16.

--------------------------------------------------------------------------------------------
## 3. THE COMPLETE REDUCTION (this is the main result)
A fully calibrated model of all 383 nodes and all 256 live leaves now exists and is validated:

* `buildall.py` -> `full_model.pkl` : 383 nodes, **all 383 with a clean two-coordinate 3-way
  mux**; 256 live leaves, **all 256 with both pin constants extracted numerically, 0 conflicts**;
  128 dead leaves all literal 0.
* `calib2.py` -> `calib2.pkl` : per-node parent/child **coordinate alignment** (from the slot
  links; 4 recovered numerically) and per-node **chord orientation**, obtained by solving each
  node's own 3 stage checks and matching against chordK with the universal
  `K = 97553848499418123410591666447050222001188385549510401465815187079080512838891`.
  **383/383 calibrated, 0 failures** (188 orient=1, 67 orient=0, 128 dead).
* `mkassign2.py` : ON-set -> full integer assignment (fold mod p, then handles lifted to Z).

**Measured, for every ON-set tried (|S| = 1, 2, 5, 73, 200, 256 — random and extreme):
the whole 39,033-atom system collapses to exactly the SAME 4 atoms, and to exactly 2 once the
two target wires are pinned:**

    (x18956 - x37892) - p*h                      ->  x37892 == T2  (mod p)
    (x24468 - x13682) - 12354891*p*h             ->  x13682 == T1  (mod p)

x13682 and x37892 are the root mux output wires.  So:

> **MOD P**, the 39,033-equation system reduces to: some non-empty subset S of the 256 live
> leaves folds, through the calibrated tree, to the pair
>
>     TARGET (root coord order) =
>       ( 44859544763832475231923253825569092119321525945631045653619508440821028887,
>         36200939269128454586076546451607958467047992891178506183612554289882454126226 )

Everything else closes mod p for any S.  **Over Z there are 927 further conditions** (see the
c > 1 bullet below); they are sparse but they are NOT discharged for free, so the criterion
above is necessary, and sufficient only once those 927 are also solved.  Verified over Z end to
end for |S| = 1 only: `assign_L1.json`, exact checker, 15 failing equations, all attributable to
the two TARGET atoms.

Supporting facts I measured with **my own knob set (all 8,747 free vars)**:
* `handles2.py` / **`hcheck.py` — CORRECTED, agent P was right.**  My "appears in exactly one
  atom" was measured on the **free cofactor u**, not on the P-multiple h.  Measured: my 3,681
  "handle" vars are all FREE, appear in **0 residual atoms directly** and in **exactly 1
  definition** — they are the cofactors.  Cofactor-appears-once is necessary for freedom but
  **not sufficient**: the atom must also be solvable over Z, and `c*p | R` is strictly stronger
  than `R == 0 mod p` when c > 1.  Splitting my 3,681 atoms by the measured multiplier c:
  **c == 1 for 2,747** (condition collapses to the mod-p congruence) and **c > 1 for 927**
  (real extra integer condition), 7 with zero slope.  **My 927 matches P's 927 exactly**, from
  an independent decomposition.  **So the lift is NOT free, and S3 must be read mod p.**
  |Z| (wires == 0 mod p for every assignment) = 7,202.
* **How binding are the 927?  (`repaircount.py`, my own measurement, which P does not have.)**
  They are real but *sparsely* violated, and they grow with the ON-set:
      |S| = 1  : 2 of 927 violated, 2 distinct atoms, greedy repair discharges ALL -> 0 left
      |S| = 2  : 4 distinct atoms violated, greedy repair leaves 1 undischarged
      |S| = 17 : 36 distinct atoms violated, greedy repair leaves 8 undischarged
  Every violated atom is inside the c > 1 set (verified).  My repair is a greedy round-robin
  that shifts one value wire by a multiple of p per atom; it cycles for larger |S| (hits the
  60-round cap).  **What is needed is a simultaneous CRT solve over the ~766 shift parameters,
  not a round-robin** — that is exactly P's rank question, and these numbers are the data for it.
* `slopes.py`: every handle's slope is divisible by p (2,747 are exactly +-p, the rest +-M*p).
* **The deliverable's ON-set: TWO selector bits, ONE propagating leaf** (`onset_deliv.py`,
  `rootcheck.py`, `delivsite.py`).  My earlier "one leaf ON" was read off a STALE partial model
  and I retract it as stated.  Measured from the deliverable's own JSON — all 256 leaves are
  free vars, so the file value IS the bit, no inference:
  **exactly x2081 and x24601 are set to 1**, everything else 0.  M, K, P and R are right.
  LCA(2081,24601) = the ROOT x9274; 24601 sits under the root's a-child (178-side), 2081 under
  the b-child (78-side); root sel_ab = 1.  Root cause of my error: my stale sub-forest covered
  only the 178-side, so it could not see 2081 — same root cause as my retracted 2^178.
  **What the fleet was disagreeing about**: I was reading what *reaches* the root, which is one
  leaf; M/K/P/R were reading the selector configuration, which is two.  Q's census settles it
  the same way.  Both are true of different objects, and the distinction is the mechanism:

**MEASURED MECHANISM OF THE DELIVERABLE (`rootcheck.py`, hard numbers):**
    root va input = root vb input, EXACTLY, in both coordinates
        = (37841415183514949237467304684128824427406379377151921996714091976892367869714,
           82007976112976807461901870199198737303514020147647909878034348606308756230357)
        = my model's value for leaf 24601 transported to the root frame
    root vab wires x30213, x22162 = TARGET, EXACTLY, the pair I derived independently in S3.
  So the 2081 branch is overwritten to carry 24601's value; the root then sees two equal inputs
  and its own checks stop constraining its output; the output is set straight to TARGET and
  everything above closes.  **The deliverable's independently-set root wires holding exactly my
  independently-derived TARGET is the strongest single cross-check of S3 in this file.**
  Cut site (`delivsite.py`): child node **x27994** (the node holding leaf 2081), parent
  **x4971.va**.  Broken atoms: 2 guards on x27994's vab wires (x8731, x9118) + 2 slot links on
  x4971's va wires (x4432, x7068) = 4 atoms = 7 equations.

--------------------------------------------------------------------------------------------
## 4. WHY 39,026 IS NOT EASY TO BEAT (priced, not asserted)
`price.py` prices every atom by how many equations contain it, and every **slot-link pair**
(breaking a pair lets the fold hit the target above it while only those 2 atoms fail):
* cheapest slot-link pair anywhere in the tree: **cost 11** (node x4971 va, wires 4432/7068;
  and node x23242 va) — worse than the deliverable's 7.
* the root pair my constructor breaks: cost 15 (hence 39,018).
* cost histogram over 378 pairs: 11:2, 12:18, 13:32, 14:59, 15:72, 16:79, 17:69, 18:24, 19:18, 20:3, 21:2.
The deliverable reaches 7 with **four** nonzero atoms, i.e. it exploits **cancellation** between
atoms inside shared equations — a cheaper cut than any 2-atom cut in this family.
**Knob set for that statement: the 378 parent-slot-wire pairs of the 383-node tree, at any
selector configuration** (the pricing is configuration-independent — it is incidence only).
No claim is made about other break sets.

--------------------------------------------------------------------------------------------
## 5. SEARCH STATUS
`subsearch.py` (LCA + cumulative coordinate swaps, O(1) per fold; validated 300/300 against the
full circuit-driven fold) and `s3.py` (validated 200/200):
* **|S| = 1 : 256 folds, no hit.**
* **|S| = 2 : 32,640 folds, no hit.**
* **|S| = 3 : 2,763,520 folds, NO HIT** (`s3b.py` / `s3b.log`, 143 s).
* **|S| = 4 : 174,792,640 folds — RUNNING** (`s4.py` / `s4.log`, measured rate 400k / 56 s => ~6.8 h).  If it is
  still running or was killed, relaunch with `setsid nohup python3 s4.py > s4.log 2>&1 &`.
* **No degeneracy anywhere** (`degen.py`): 0 of 32,640 pairs and 0 of ~370 random sets of size
  3..256 produce a zero chord denominator, so the fold is well defined on everything tested —
  which is what the 2^256-1 count needs.
`s3b.py` uses **Montgomery batch inversion** (`batch_inv`): 200k inversions in 1.4 s vs 849 us
each naively — a 120x speedup, and the whole cost of a fold is its inversions.  Rate measured:
400,000 triples / 22 s.  At that rate **|S| = 4 (174,792,640 sets) is ~2.7 h** — the next run to
launch, and it needs only a generic "merge the pair with the deepest LCA" loop over 4 items.
Root MITM split is 178 | 78 live leaves, so a root-level meet-in-the-middle is 2^78 — out of
reach.  12 of 383 nodes have live-leaf support > 24; 371 are <= 24.

--------------------------------------------------------------------------------------------
## 6. NEXT, IN ORDER
1. Finish |S|=3; then |S|=4 (174M) with batched modular inversion (Montgomery trick) — the
   per-fold cost is currently ~1 inversion per chord and that is the whole budget.
2. Invert the target *down* the 78-side: `invchord` is in `fastfold.py` and is exact.  At each
   node on the 78-side spine the sibling's value is a function only of that sibling's subtree,
   so the target inverts to a required value at a node of support <= 24, where forward
   enumeration is 2^24.  The blocker is that the 178-side value must be guessed first — unless
   the 178 side is entirely OFF, which is one cheap case worth running to completion.
3. Cancellation-aware cut search: find break sets of 3-5 atoms whose equations cancel down to
   < 7.  My pricing only covers 2-atom cuts.  **Tried and failed the naive version**
   (`cut.py`): inverting the target down to a cheap node and putting the required value on that
   node's parent slot wires leaves 6 nonzero atoms and 39-47 failing equations, far worse than
   the incidence bound of 11-13 — the extra atoms come from divisibility repair spilling.  The
   deliverable reaches 7 with FOUR nonzero atoms, i.e. by arranging *value* cancellation inside
   shared equations, not by minimising incidence.  Any attempt to beat 7 has to search
   cancellation, not support.

--------------------------------------------------------------------------------------------
## 6b. CANCELLATION SEARCH — instrument BUILT and CALIBRATED, family MIS-SPECIFIED (open)
`cancel.py` / `cansearch.py`.

**What works.**  An **exact in-memory scorer** now exists: `CK.load_equations()` once (~91-153 s)
then ~1.1 s per candidate, so hundreds of placements can be priced exactly instead of by
incidence.  Calibrated on two known points: **deliverable -> 7** and **my `assign_L1` -> 15**,
both matching `checker.py` exactly.  This is the instrument the lab was missing, and it is the
thing to reuse.  **Do not price with `E.score`** — it reports 13 for the deliverable where the
truth is 7, so every incidence/F-model number in S4 and `cut.py` is inflated and only ordinally
useful.

**What does not work yet.**  My generalisation of the deliverable's cut — ON = {G, L}; overwrite
the L branch at a chosen site with G's value so LCA(G,L) sees two equal inputs; then drive that
node's vab wires to TARGET — **does not reproduce the deliverable even at the deliverable's own
site**:
    G=24601 L=2081 csite=27994 set_vab=True  -> 9 nonzero atoms, EXACT failing 49
    G=24601 L=2081 csite=27994 set_vab=False -> 7 nonzero atoms, EXACT failing 47
against the deliverable's 4 atoms / 7 equations.  So the family as I specified it is wrong: it
breaks 5 extra atoms.  **The next step is diagnostic, not a sweep** — diff my 9 broken atoms
against the deliverable's known 4
    ((x7075*x8731)+x31864), ((5113045*(x7075*x9118))-x29854),
    ((x4432-x19964)-x28730),  ((x7068-x2099)-(7376877*x642))
and find which wires the deliverable leaves alone that I am overwriting.  Prime suspects: the
top slot links x24468/x18956 -> x13682/x37892, which I pin to T1/T2 while also rewriting the
root vab wires, and the un-converged divisibility repair (see the c>1 bullet).  Only once the
constructor reproduces 7 at the known site is a sweep over sites meaningful.

**Caveat inherited from the c > 1 finding:** any candidate whose greedy repair does not converge
carries extra nonzero atoms that have nothing to do with the cut, which will masquerade as a bad
placement.  Fix the repair (simultaneous CRT) before trusting a sweep.

## 7. FILES (all in `agentL_work/`)
Code: `trace.py ortree.py ortree2.py census.py wire.py link.py crux.py onset.py fail7.py
handles.py handles2.py exp1.py model.py model2.py calib.py fold.py fold2.py global.py
buildall.py calib2.py fold3.py lift.py slopes.py mkassign.py mkassign2.py dbg.py target.py
fastfold.py subsearch.py val2.py price.py s3.py`
Data: `full_model.pkl calib2.pkl pins.pkl handles.pkl slopes.pkl price.pkl target.pkl
ortree2.pkl nodes.pkl outwires.pkl ors.pkl assign_L1.json s3.log slopes.log`
Rebuild from cold: `python3 handles2.py; python3 global.py; python3 buildall.py;
python3 calib2.py; python3 slopes.py; python3 target.py` (~5 min total).
