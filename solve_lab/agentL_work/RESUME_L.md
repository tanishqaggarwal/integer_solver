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
  that shifts one value wire by a multiple of p per atom; it cycles for larger |S|.
  **`repairfix.py`: reordering the shifts bottom-up (deepest wire first) does NOT help** —
  identical results to round-robin at every size, so the residue is genuinely a simultaneous
  system, not an ordering artefact.  Undischarged count scales with |S|:
      |S| =  1 -> 0      |S| =  2 -> 1      |S| = 17 -> 9      |S| = 40 -> 21
  **A simultaneous CRT solve over the ~766 shift parameters is required** — that is exactly P's
  rank question, and these numbers are the data for it.
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
## 6b. CANCELLATION SEARCH — instrument BUILT, family FIXED, cost is a VALUE choice (DIAGNOSED)
`cancel.py` / `cansearch.py`.

**What works.**  An **exact in-memory scorer** now exists: `CK.load_equations()` once (~91-153 s)
then ~1.1 s per candidate, so hundreds of placements can be priced exactly instead of by
incidence.  Calibrated on two known points: **deliverable -> 7** and **my `assign_L1` -> 15**,
both matching `checker.py` exactly.  This is the instrument the lab was missing, and it is the
thing to reuse.  **Do not price with `E.score`** — it reports 13 for the deliverable where the
truth is 7, so every incidence/F-model number in S4 and `cut.py` is inflated and only ordinally
useful.

**FAMILY MISSPECIFICATION — FOUND AND FIXED (`diffcut.py`, `cansearch2.py`).**
The bug: I injected the forged value at ONE site and never re-propagated it up the branch, so
every slot wire above the cut still carried leaf 2081's value.  `diffcut.py` showed it directly
— the deliverable carries leaf 24601's value on the WHOLE 2081 branch (x28505.va, x16102.vb,
x23131.vb, x13976.va, x17215.vb, x9274.vb) while mine carried 2081's.  Hence my 5 extra atoms:
2 slot links above the cut, and the 3 ROOT stage checks (x15298), which broke because the root's
two inputs were still unequal.  `build2()` in `cansearch2.py` propagates from the cut site up to
m and **now breaks EXACTLY the deliverable's four atoms** — identical support.

**THE RESULT THAT MATTERS.  Same support, different cost:**
    my build2, same 4 atoms as the deliverable  -> EXACT failing 13
    the deliverable, same 4 atoms               -> EXACT failing  7
    my build2 with vab left at 0, only 2 atoms  -> EXACT failing 11
So (a) **cancellation is a VALUE property, not a support property** — proven, because the
support is now byte-identical and the cost still differs by 6; and (b) **more broken atoms can
cost FEWER equations** (4 atoms -> 11 vs 2 atoms -> 13... in fact 4->13 and 2->11 for mine, and
4->7 for the deliverable), so minimising atom count is the wrong objective.

**WHERE THE 6 EQUATIONS LIVE (`valdiff.py`).**  With support identical, only **12 variables**
differ mod p between my assignment and the deliverable, and they are all cofactors/handles:
    x105, x1329, x3387, x5081, x5676, x9413, x10903, x11436, x14393, x14768, x17325, x22820
The deliverable sets them to specific nonzero integers; **my constructor leaves them at 0**,
because `relift` skips exactly the atoms that are nonzero mod p, so their handles never get set.
**That is the whole cancellation degree of freedom.**  The site fixes WHICH atoms break; the
handle values fix HOW MANY equations they cost.  The search is therefore
**site x handle-values**, not site alone — which is precisely what M's primitive prices.

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

--------------------------------------------------------------------------------------------
## 6c. THE USEFUL SEARCH SPACE IS 15 ATOMS (`refilter2.py`, `atomfilter.py`, `emit.py`)
M's filter — *a site can help only if its corrupted atoms appear in the equations that fail at
the uncorrupted baseline* — applied exactly.

**First: my incidence map was wrong and I rebuilt it.**  The map I priced with in S4/S6b came
from F's `E.eqres` and sees only 12-13 of the 25 target equations.  Rebuilt exactly instead:
every residual atom has exactly ONE free cofactor u (3,681, verified in `hcheck.py`) and u
occurs nowhere else, so **equation e contains atom a  <=>  u_a in vars(e)** — read straight off
`checker.load_equations()`'s varsets, no model of the equation algebra needed.  `eqmap_exact.pkl`.

**BASELINE DISCREPANCY — FLAGGED, NOT RESOLVED.**  M's baseline fails **25** equations; mine
fails **13**, and **my 13 are a strict subset of M's 25** (`B \ M25 = {}`).  The 12 M sees and
I do not: 5324, 9041, 11226, 15558, 21000, 22534, 22997, 28929, 29330, 32026, 35512, 38051.
My baseline = the deliverable with its 16 tuned handle/cofactor vars zeroed.  I filtered on the
**UNION (25)**, not the intersection — the intersection is my 13 and would over-discard anything
that could only fix one of M's extra 12.  Whoever reconciles these should use the union until
the two baselines are the same object.

**RESULT: 15 incident atoms -- CORRECTED BY T TO 18.**  My census was scoped to slot-link
guards; T's 33 extra genuine p-handles (guards = stage checks and leaf pins) include **3 more
incident** ones (u=x10422, u=x15120, u=x35531, all hitting 12231/12350/14584/29125).  My
*criterion* is sound -- T verified `eqs(u)==eqs(atom_u)` at 3,681/3,681, zero violations -- what
was too narrow is the FAMILY it ranged over.  Use the full p-handle family (3,707 / 3,714).
All 3,666 others provably cannot change any target equation, whatever value they are given.
Grouped by the node they guard (h = the P-multiple to corrupt, u = its free cofactor):

    node x27994 sel_ab (its 3 STAGE CHECKS + 2 guards)     5 atoms, total rt 41
        rt 10  h=x23754  u=x6947    ((6122989*(x21279*x2239))-x23754)     <- NOT corrupted by
        rt  9  h=x35619  u=x33168   ((x21279*x31731)+x35619)                 the deliverable
        rt  8  h=x31864  u=x10903   ((x7075*x8731)+x31864)                <- deliverable
        rt  8  h=x9629   u=x950     ((x21279*x9106)-(13523997*x9629))     <- NOT corrupted
        rt  6  h=x29854  u=x1329    ((5113045*(x7075*x9118))-x29854)      <- deliverable
    slot links (x4971.va, x4971.vb, x36871.vb)             5 atoms, total rt 26
        rt 10  h=x28730  u=x9413    ((x4432-x19964)-x28730)               <- deliverable
        rt  9  h=x642    u=x17325   ((x7068-x2099)-(7376877*x642))        <- deliverable
        rt  5  h=x37413  u=x11099   ((x15324-x37254)-(8481759*x37413))    <- x4971.vb, the OTHER side
        rt  1  h=x34113  u=x23110   ((x22043-x4264)-x34113)
        rt  1  h=x28355  u=x32349   ((179131*(x27500-x12143))-x28355)
    node x35155 sel_ab                                     3 atoms, total rt 10
        rt  4  h=x1844   u=x21574   ((x1956*x17065)-x1844)
        rt  3  h=x29305  u=x1613    ((x1956*x23318)-(8235511*x29305))
        rt  3  h=x2892   u=x6090    ((x1956*x14199)+x2892)
    node x14803 sel_ab                                     2 atoms, total rt  2
        rt  1  h=x23822  u=x22526   ((x16495*x6247)-x23822)
        rt  1  h=x7945   u=x34868   ((7720481*(x16495*x1504))-x7945)

**So the ENTIRE useful corruption space is the subsets of these 15 handles = 32,768 candidates**
— enumerable outright at M's stated throughput, not a heuristic shortlist.  The deliverable is
the 4-subset {642, 28730, 29854, 31864} = 39,026.  **The strongest lead: the deliverable does
not touch x23754 / x35619 / x9629 — the three stage checks of the very node it cuts at — and
they carry the highest incidence in the whole system (rt 10/9/8).**  Next after those, x37413
(rt 5), the sibling slot x4971.vb.  `emit_for_M.json` carries the 15 atoms and all 2,048
supersets of the deliverable's set in descending total rt.
**Caveat on the ordering only:** total-rt is a sum of incidences and S6b proved cost is a value
property, so treat the order as an enumeration convenience, not a prediction.

--------------------------------------------------------------------------------------------
## 6d. REALIZABILITY OF THE 15 (`roles32.py`, `realiz.py`) — answers "which subsets are cuts?"
The 15 incident atoms sit on exactly FOUR nodes, in three mechanistic classes:

    x27994 vab guards      x31864 x29854      driven only when sel_ab(x27994) == 0
    x27994 stage checks    x23754 x35619 x9629 driven only when sel_ab(x27994) == 1
    x4971.va slot links    x28730 x642        driven always
    x4971.vb slot link     x37413             driven always
    x36871.vb slot links   x34113 x28355      driven always
    x35155 stage checks    x1844 x29305 x2892  NEVER driven  (see below)
    x14803 stage checks    x23822 x7945        NEVER driven  (see below)

**HARD RESULT 1 — five of the 15 are permanently vacuous.**  A node's stage checks are gated by
sel_ab, and sel_ab == 1 requires BOTH child subtrees to contain a live leaf.  Live-leaf counts:
    x27994: 1 | 1  -> sel_ab CAN be 1        x4971 : 2 | 1  -> CAN be 1
    x36871: 2 | 1  -> CAN be 1
    x35155: 1 | 0  -> **sel_ab is 0 in ALL 2^256 configurations**
    x14803: 1 | 0  -> **sel_ab is 0 in ALL 2^256 configurations**
So x1844, x29305, x2892, x23822, x7945 can never carry a circuit-derived value in any reachable
configuration whatsoever.

**HARD RESULT 2 — at x27994 the guards and the stage checks are mutually exclusive.**
Guards are `(1-sel_ab)*vab`; checks are `sel_ab*(...)`.  At sel_ab=0 the guards bite and the
checks are vacuous; at sel_ab=1 the reverse.  The deliverable runs at sel_ab(x27994)=0
(measured, x21279=0), which is **why it does not corrupt x23754/x35619/x9629 — they are vacuous
there, not overlooked.**

**BUT — DO NOT READ THIS AS "M IS PRICING NOISE".**  Freeing a handle always works formally: at
sel_ab=0 the check `((x21279*x9106)-(13523997*x9629))` collapses to `-13523997*x9629`, so
demoting x9629 supplies a **pure unconstrained additive term** in whichever equations contain it.
That is arguably the *best* cancellation knob available, precisely because it is uncoupled from
circuit values.  So every one of the 32,768 subsets is runnable and worth pricing.  What changes
is the INTERPRETATION: a hit on the vacuous atoms is a **cancellation knob**, not a structural
cut, and it will not generalise to other sites or other ON-sets the way a cut would.
The realizability filter therefore does NOT reduce M's 35,960 — I could not cut it, and saying so
is the answer.  What it gives instead is a partition: 10 of the 15 can be circuit-driven, 5 can
only ever be knobs, and at x27994 no assignment drives guards and checks at once.

--------------------------------------------------------------------------------------------
## 6e. THE ALIAS SLACK IS p x (FREE VAR) — Q'S GATE, ANSWERED  (`slack.py`, `slack2.py`)
Q asked: is anything forcing x_4116 and its five sibling shared factors to zero?

**The question dissolves: x_4116 is not a variable that could be nonzero.  IT IS THE CONSTANT p.**
Evaluated from their definitions alone, with no free variable anywhere in the chain:

    x4116 = x16153 = x1962 = x12682 = x19049 = x15616
          = 115792089237316195423570985008687907853269984665640564039457584007908834671663 = p

They are 6 of the **220 constant-p wires** in the instance (every constant multiple of p present
has multiplier exactly 1).  **That is why they carry no unary pin — constants do not need
pinning**, and it is why they are "shared": it is the same constant reused, 66 times for x4116.

**GENERAL, NOT ANECDOTAL (`slack2.py`).**  Of the 12,232 wires defined as a product of two wires,
the 3,681 that appear as slack in a residual atom (one free cofactor) are
**3,681 / 3,681 of the form (constant multiple of p) x (free variable).  ZERO exceptions.**
The other 7,697 products are the selector products sel*value, a different population.

**CONSEQUENCE**
    slack == 0  (mod p)  in EVERY assignment, unconditionally
    => parent_input == mux_out  (mod p)  exactly
    => **the coordinate hand-off follows the measured tree, unconditionally, mod p**
    over Z the slack is a free multiple of p — the lift freedom.

This corroborates my `slopes.py` result (all 3,681 handle slopes divisible by p, 0 exceptions)
from a completely independent direction: that was numerical, this is structural.

**WHAT THIS DOES AND DOES NOT CLOSE.**  It closes the hand-off **mod p** — Q's worry that the
parent input is "the mux output plus an unpinned amount" is answered: the amount is p x (free
var), which is 0 mod p.  It does **NOT** close the integer statement.  Over Z that slack is
genuinely free, and its residue is exactly the **927 c>1 divisibility conditions** of the
c > 1 bullet, which remain open (P's rank question).  So: **the mod-p reduction now closes on
measurement; the Z statement still needs the 927.**  Anyone reporting the existence result as
unconditional must say "mod p" or must first discharge the 927.

**T's corrections to S6b, accepted:** the gap is **5, not 6**; the far side is **12, not 13**
(the 13 was my own `build2`'s score, inflated by the un-converged repair I flagged); and only
**4 of the 12 cofactors move anything** — x1329 (+3), x9413 (+4), x10903 (+3), x17325 (+4) —
so **the cofactor freedom is 4-dimensional, not 12**.  The core claim (cancellation is a value
property, support byte-identical) stands; my list was longer than the effect.

--------------------------------------------------------------------------------------------
## 6f. WHY THE ROUND-ROBIN CANNOT WORK: THE SHIFT SYSTEM IS NONLINEAR (`crt.py`)
Inherited P's method (factor, prime-by-prime, CRT) and both guards.  **P's three files are NOT
runnable here** — `plift5.py`, `prank.py`, `pcompose2.py` all load `model4.pkl`, `slp.pkl`,
`blocks.pkl`, `leaves.pkl` from `agentP_work/` and `import pfold`, none of which were copied.
I did not reach into that directory.  I took the method onto my own model instead.

**Measured shape of the stuck conditions** (run the greedy repair to its fixpoint, then probe
each surviving condition at shift t = 0, 1, 2 and check whether d(2) == 2*d(1)):

    |S| = 2   1 stuck:  c = 6672769 (prime),  6 influencing wires,  **NONLINEAR**
    |S| = 17  8 stuck:  c = 15194385 = 3^4*5*37517   NONLINEAR
                        c = 3849267  = 3*19*67531    NONLINEAR
                        c = 6672769  (prime)         NONLINEAR
                        c = 10696593 = 3*3565531     NONLINEAR
                        c = 10353929 = 127*81527     NONLINEAR
                        c = 2264251  = 11*43*4787    NONLINEAR
                        c = 10937191 = 449*24359     **LINEAR**, wires (23238,+1) (2964,-1)
                        c = 13040669 = 19*199*3449   **LINEAR**, wires (10261,+1) (27156,-1)

**This is the actual reason the round-robin fails, and it is not an ordering artefact** (S3
already ruled ordering out).  Two things are going on at once:
 1. The two LINEAR conditions have d/p = +-1 on a *shared* wire, so greedy fixes one and the
    next fix re-breaks it — pure simultaneity.
 2. The other six are genuinely **nonlinear** in the shifts: a shift enters the chord law
    through a product, so after dividing by p a term p*t_w*t_v survives mod c.  A linear
    solve over the shift parameters cannot express them.

**This independently corroborates P's own expansion** from a different model: P's `n1'` carries
`P*(E*b^2 + 2*a*A*b - d^2) + P^2*a*b^2` — quadratic and cubic in the shift parameters.  P had
this right; my measurement confirms the nonlinearity is real and not specific to P's block.
**So "a simultaneous CRT solve over the ~766 shift parameters" — my own phrasing, twice — is
not by itself enough: the system is polynomial, not linear.**

**CONCRETE NEXT STEP (cheap, and the right one).**  Do not brute-force: c is up to ~1.5e7 and
E.run is ~0.07 s.  Instead **fit and solve exactly**: for a chosen wire the atom is a polynomial
in t of degree <= 3 (P's expansion bounds it), so evaluate at t = 0,1,2,3, interpolate the
coefficients exactly, then root-find mod each prime factor of c and CRT — c factors into small
primes in 7 of the 8 cases above (3, 5, 11, 19, 43, 127, 199, 449, 3449, 4787, ...), where root
finding is trivial.  Then **verify by direct recomputation**, per P's second guard.  The two
linear conditions should be solved jointly with the rest, not greedily.

--------------------------------------------------------------------------------------------
## 6g. FIT-AND-SOLVE WORKS PER CONDITION; |S|=2 CLOSES OVER Z  (`solve927.py`)
Ran the recipe.  Per condition: evaluate R(t) at t=0..4, exact Newton forward-difference fit,
root-find mod each prime power of c, CRT, then **verify by direct recomputation** (P's guard).

**RESULT 1 — |S| = 2 CLOSES COMPLETELY OVER Z.**
    round 0: 1 stuck -> SOLVED c=6672769 (prime) deg=2 wire x24908 t=2990790, verified
    FINAL: **0 undischarged, 2 nonzero atoms** — and those 2 are the target congruences.
This is the first ON-set beyond |S|=1 for which every one of the 927 conditions is discharged.

**RESULT 2 — P's degree bound is CONFIRMED on a second model.**  I fitted to degree 4 and
recorded the top nonzero degree every time.  Observed degrees: **1, 2 and 3.  A degree-4 term
never appeared, in any condition, in either ON-set.**  P's expansion bounds it at 3; measured
independently here, it holds.

**RESULT 3 — cost is driven by the largest prime factor of c, not by c.**  c = 6672769 is prime
and takes ~59 s (6.6M evaluations of the *fitted polynomial* — cheap, unlike E.run);
c = 15194385 = 3^4*5*37517 takes ~1 s.  Brute-forcing the fit is fine; brute-forcing E.run is not.

**REMAINING GAP, PRECISELY LOCATED — the solver is still per-condition, so it oscillates on
SHARED wires.**  |S| = 17 goes 8 stuck -> 3 -> 3 -> 3.  Every individual condition is solved and
verified exactly, but the same wires recur — **x23238 and x10261 each carry two different
conditions**, so clearing one re-breaks the other.  This is the same simultaneity that defeated
the greedy round-robin, now one level up: I replaced "linear solve" with "exact polynomial solve"
but not "one at a time" with "jointly".

**THE FIX, AND IT IS SMALL.**  For a wire w carrying conditions {a_1..a_k}, do not solve them
one at a time: fit each a_j as a polynomial in t_w, root-find each mod its own c_j, and take the
**INTERSECTION of the root sets** via CRT across the distinct c_j — a t_w that clears all of them
at once, or a proof that none exists, in which case fall back to a different wire for one of
them.  Only 2-3 wires are contended, so the intersection is over 2 conditions at a time.
Everything needed is already in `solve927.py`; it is `solve_one` that needs to become
`solve_group`.

--------------------------------------------------------------------------------------------
## 6h. JOINT SOLVE: |S|=17 GOES 8 -> 2, AND THE RESIDUE IS BIVARIATE  (`solve927g.py`)
`solve_one` -> `solve_group`: for a contended wire, intersect root sets prime-power by
prime-power, CRT, then verify by direct recomputation.

**It works, and same-wire contention is gone.**  Round 0 cleared **2 conditions jointly on one
wire** (x23238, t=79784602390776, verified) -- exactly the case that oscillated before.
Undischarged went **8 -> 2** and stayed at 2.

**THE RESIDUE IS A DIFFERENT OBSTRUCTION, NOT OSCILLATION FROM GREED.**  The surviving pair sits
on **two different wires, x9776 and x10261**.  Each is individually solvable and verified every
round -- but clearing one on x9776 breaks the other on x10261 and vice versa, so the loop cycles
    round 2: x9776 t=1890710 , x10261 t=1550230
    round 3: x9776 t=6051501 , x10261 t=1345905
    round 4: x9776 t=4302428 , x10261 t=11694764   ... stable 2-cycle
Per-wire grouping cannot reach this: the coupling is **across wires**, so it is a genuine
**bivariate** system, not a contended single wire.

**SUPERSEDED BY S6i — the bivariate reading was WRONG.  See below.**

**(original, retained for the record) THE FIX, STILL BOUNDED.**  Solve the pair simultaneously in (t1,t2): each atom is degree <= 3
in each variable (bound confirmed in S6g), so for each prime power q^e of the two moduli, **loop
t1 over q^e and root-find t2 from the resulting univariate polynomial** -- one loop, not a double
loop, so ~10^7 cheap polynomial evaluations, comparable to the 59 s already spent on a single
prime c.  Then CRT across prime powers and verify by recomputation.  Generally: `solve_group`
must range over wire SETS, taken from the connected components of the "shares a condition"
graph -- here exactly one component of size 2.

--------------------------------------------------------------------------------------------
## 6i. THE COUPLING IS WITH ALREADY-SATISFIED CONDITIONS, NOT AMONG THE STUCK ONES (`bivar.py`)
**I was wrong in S6h and the component computation says so.**

    COMPONENT SIZES at |S|=17 of the shares-a-condition graph:  [1, 1]

The two residual conditions are in **separate components — they share no influencing wire at
all**.  Measured: c1 = 1707229 = 43*39703 and c2 = 5930437 (prime), **coprime**; and on the
wire pair (x5460,x5616) atom1 has degree (2,1) while atom2 has degree **(0,0)** — those wires do
not touch atom2.  So there is no bivariate coupling between them.  (Degree-<=3 per variable does
hold separately: (2,1) observed, nothing above 3 anywhere.)

**So why does it cycle?**  Because `solve_group` groups only the currently-STUCK atoms on a wire
and ignores the c>1 atoms that wire influences which are currently SATISFIED.  Clearing a stuck
condition on wire w silently breaks an already-discharged one elsewhere, which then reappears as
stuck next round.  **The 2-cycle is between the stuck set and the satisfied set, not within the
stuck set.**  That is also why the residual pair is not stable across runs: at S6h it was
(x9776,x10261) with different moduli, here it is a different pair — the residue is
path-dependent, which a genuine structural obstruction would not be.

**THE ACTUAL FIX.**  For wire w, `solve_group` must range over **every c>1 atom that w
influences**, not just the violated ones, and require:
    violated atoms  -> t in their root set   (clear them)
    satisfied atoms -> t in their root set   (PRESERVE them; note t=0 is always in it)
i.e. intersect root sets over the whole influenced set, which is exactly the machinery already
written — only the atom list passed in is wrong.  This is a one-line change to how `ats` is
built in `solve927g.py`, plus keeping the existing direct-recomputation guard.

## 6ii. METRIC + ARTIFACT DISCIPLINE (T's audit, accepted in full)
* **Report NONZERO ATOMS (of 9,032), not the "stuck" count.**  My stuck list was
  `[a for a in relift(vv) if r[...] % p == 0]`, which **silently drops bad-list entries whose
  residual is not 0 mod p**.  At |S|=2 the two dropped entries are exactly the target
  congruences (benign); at |S|=17 the two metrics genuinely differ.  The complete number is the
  nonzero-atom count.  Every earlier "0 undischarged" in this file should be read with that
  caveat; the companion nonzero-atom figures quoted alongside them are the sound ones.
* **ALWAYS dump the assignment and run `checker.py`.**  `solve927.py` dumped nothing, so the
  |S|=2 closure was model-internal until T reproduced it, dumped it and checked it:
  `39018/39033`, 2 nonzero atoms = the two target congruences, and their 15-equation footprint
  is EXACTLY the checker's failing set.  The closure is real, but I should have produced that
  artifact myself.  `closeS.py` now dumps `close_<tag>.json` for every run.
* T also confirmed the degree bound a third time (deg <= 6, 8, 10 re-fits, same top degrees) and
  established it **bounds cost, not correctness** — a wrong bound can only lose a solution, never
  admit a false one, because the recomputation guard rejects a bad root.

## 6j. PROCESS RULE (learned twice, flagged, now recorded)
* **Never `pkill -f <pattern>` where the pattern can match this shell.**  `pkill -f solve927g`
  matched the wrapping shell and killed it (exit 144), twice.  Record the PID at launch
  (`echo $!` / `os.getpid()`) and kill that, or match on a token unique to the job.
* **Never split a source file on a literal that the file itself contains.**  Splitting
  `solve927g.py` on `"if __name__"` cut inside its own `.split("if __name__")` call.  Use an
  explicit `#MAINSTART` marker comment — three occurrences of this bug this session.

**PERFORMANCE BLOCKER IN `closeS.py` (found the hard way, fix is easy).**  The S6i fix needs,
for each wire w, the set of c>1 atoms w influences.  I computed it by *probing* — `influences()`
calls `E.run` twice per (atom,wire), so 927 atoms x ~0.14 s = ~130 s **per wire**, and the run
did not finish.  **Build that map once, structurally**, from `atomvalvars` / `vars_of` (which is
already in memory) instead of probing: wire -> [atoms mentioning it].  Probing is then only
needed to confirm a nonzero derivative on the handful of surviving candidates.  With that change
the |S| = 3, 5, 8, 17 sweep is minutes, not hours.

--------------------------------------------------------------------------------------------
## 6k. STRUCTURAL MAP FIXED THE FIRST BLOCKER; A SECOND ONE IS NOW THE BINDING COST (`closeS2.py`)
The influence map is now built once, structurally, from `vars_of`/`atomvalvars` — no probing:

    structural influence map: 1901 wires, mean 1.4 c>1 atoms/wire, **max 3**

So the groups are tiny and the S6i fix is cheap, exactly as predicted.  **The 130 s/wire probing
cost is gone.**

**But the sweep still does not finish, and the reason is a DIFFERENT bottleneck: `rootset_pp`.**
It computes the FULL root set by enumerating `t` over `q^e`.  For a large prime modulus
(c = 5930437, c = 6672769) that is ~6M `peval` calls, ~30-60 s — and the S6i fix now calls it
**once per atom per wire, including every `keep` atom**, so the cost multiplied rather than fell.

**THE FIX, AND IT IS THE RIGHT SHAPE ANYWAY.**  Never enumerate a `keep` atom's root set — for
those we only need to TEST a candidate, not enumerate.  So:
    1. enumerate the root set of the VIOLATED atom only (one large-prime enumeration, unavoidable);
    2. for each candidate t in it, TEST every `keep` atom by evaluating its fitted polynomial —
       O(1) per test instead of O(q^e);
    3. first t passing all tests wins; then the existing direct-recomputation guard.
That turns the per-wire cost from (#atoms x q^e) into (q^e + #candidates x #atoms).  Also cache
`fit(vv,i,w)` per (atom,wire) — it is recomputed many times across the outer rounds.

**STATUS OF THE QUESTION THIS WAS MEANT TO ANSWER:** still open.  I have **no |S| = 3/5/8 data**,
so "does the integer lift close for small |S| only, or generally?" is unanswered.  |S| = 2
remains the only ON-set beyond a single leaf verified closed over Z (by T, `39018/39033`,
2 nonzero atoms = the target congruences).

--------------------------------------------------------------------------------------------
## 6l. MEASURED: COST IS FINE (186 s), CORRECTNESS IS NOT (control FAILS)  (`closeS3.py`)
Implemented the fix as specified — enumerate only the VIOLATED atom's roots; preserve every
`keep` atom by forcing `t == 0 mod c_keep`; `fit()` cached per (atom,wire,generation);
direct-recomputation guard retained.  Ran `|S| = 2` alone, timed, as the control.

    |S| = 2   NONZERO ATOMS = 8 of 9032   WALL CLOCK = 186.2 s   -> close_S2.json

**COST IS AFFORDABLE.**  186 s per configuration.  |S| = 3/5/8 would be ~10 min total.  My
earlier fear that this approach was too expensive was WRONG.

**BUT THE CONTROL FAILS.**  It must reproduce T's 2 nonzero atoms (the target congruences) and
instead gives **8**.  The extra six:
    ((x24908-x17601)+x5201)              slot link
    ((6788513*(x16742-x19083))-x9254)    root slot link
    ((x12186-x23927)-x25758)             root slot link
    ((537773*(x15298*x37758))-x35605)    ROOT stage check
    ((x15298*x11150)+x4007)              ROOT stage check
    ((x18956-x37892)-x32237)             target congruence
This is the signature of a value being shifted that must NOT move: `x15298` is the root sel_ab,
so the root's stage checks are breaking, and the root slot links are breaking with them.
**Diagnosis: `t == 0 mod c_keep` is the wrong preservation constraint.**  It preserves the keep
atom's DIVISIBILITY but the shift still moves the wire by `p*t`, and for a wire feeding the root
mux that changes a value the stage checks pin.  The keep set must include the atoms a wire
affects *structurally*, not only the `c > 1` ones — `W2A` was built over `CGT2` (the 927) alone,
so every `c == 1` atom the wire touches was invisible to the guard.
**FIX: build `W2A` over ALL atoms with a handle (3,681), not just the 927.**  The recomputation
guard then rejects any `t` that breaks a `c == 1` atom, which is exactly what leaked here.

**STATUS: the line is NOT closed, and NOT closed on cost.**  It is blocked on a correctness
regression with a named cause and a one-line fix.  No |S| = 3/5/8 data exists.

**RETRACTION, and it nearly went into the record as a result.**  I first reported this run as
"exceeded 13 minutes and did not finish" and wrote a measured-stop conclusion on that basis.
It was false: the run had finished in 186 s.  I had checked liveness with `pgrep -f closeS3`,
which matched **my own shell's command line** — the same bug as S6j, for the THIRD time, and the
first time it produced a wrong empirical claim rather than just a dead shell.
**Never test a process's liveness with a pattern that appears in the testing command.**

--------------------------------------------------------------------------------------------
## 6m. ROUND 5: GUARD REWRITTEN CORRECTLY; CONTROL PRODUCED NO RESULT  (`closeS4.py`)
**The scoping fix was the wrong fix and I replaced it with the right one.**  Widening `W2A` to
the 3,681 handle-carrying atoms would still have missed the ~5,351 atoms with **no** handle —
those cannot absorb anything and must stay exactly zero, and `((x24908-x17601)+x5201)` in the
failing list is one of them.  So `closeS4.py` drops scoping entirely and uses a **global guard**:
accept a shift only if the total nonzero-atom count strictly decreases, verified by direct
recomputation.  That subsumes every scoping question — `c > 1`, `c == 1`, and handle-less atoms
alike — and it optimises the metric we actually report.

**NO CONTROL RESULT.**  The run exited after ~110 s having printed only its two header lines: no
result line, no traceback, and `close_S2.json` was left at its earlier timestamp, so `close()`
never returned.  Checked by **recorded PID** (`c4.pid`), not by pattern — that part worked.
Cause not established; most likely process lifetime across the session rather than a code fault,
but **I have no evidence either way and am not going to guess.**  Re-run is
`setsid nohup python3 closeS4.py S2 2 > c4_S2.log 2>&1 & echo $! > c4.pid` and it needs to be
given several minutes before the log is read.

**STATE OF THE QUESTION AFTER FIVE ATTEMPTS: OPEN, WITH NO DATA.**
`|S| = 2` remains the only ON-set beyond a single leaf verified closed over Z, and that
verification is T's.  There is no |S| = 3/5/8 measurement.  What IS established: cost is
**186 s/configuration** (measured), the algorithm is correct in shape, and the remaining defect
was a guard-scope error whose replacement is written but unverified.

## 6n. PROCESS RULE, GENERALISED (three failures, one root cause)
**Never identify a process by matching a command line — record the PID at launch and check that.**
`pkill -f closeS3` killed my own shell twice; `pgrep -f "closeS3.py S2"` matched my own shell a
third time and produced a **false empirical claim** ("exceeded 13 minutes") that I filed as a
measured result and had to retract.  A rule naming `pkill` does not generalise: it recurred as
`pgrep`.  The rule is about **command-line matching as an identification method**, not about any
particular tool.  Launch with `... & echo $! > job.pid`; test with `kill -0 $(cat job.pid)`.

## 7. FILES (all in `agentL_work/`)
Code: `trace.py ortree.py ortree2.py census.py wire.py link.py crux.py onset.py fail7.py
handles.py handles2.py exp1.py model.py model2.py calib.py fold.py fold2.py global.py
buildall.py calib2.py fold3.py lift.py slopes.py mkassign.py mkassign2.py dbg.py target.py
fastfold.py subsearch.py val2.py price.py s3.py`
Data: `full_model.pkl calib2.pkl pins.pkl handles.pkl slopes.pkl price.pkl target.pkl
ortree2.pkl nodes.pkl outwires.pkl ors.pkl assign_L1.json s3.log slopes.log`
Rebuild from cold: `python3 handles2.py; python3 global.py; python3 buildall.py;
python3 calib2.py; python3 slopes.py; python3 target.py` (~5 min total).
