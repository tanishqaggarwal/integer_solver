# RESUME_P — Agent P (independent auditor of the tree claim)

**Everything below was derived from `EQUATIONS.txt` alone.** No file from `agentF_work/`
(or any other agent dir) was read, imported, or executed. FLEET.md was read only for the
statement of the claim under audit.

## 0. Baseline re-verified myself
`python3 solve_lab/agentE_work/verifyE.py solve_lab/best/new_instance_partial_39026.json`
→ **39026/39033**, failing `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`. Confirmed.
My own parser independently reproduces **exactly those 7** indices.

## 1. My parse (files in this dir, run in this order)
| file | what it does |
|---|---|
| `pparse.py`/`pparse3.py` | recursive-descent parser, polynomial algebra, `peel()` for `scalar*L^k` |
| `pparse4.py` → `model4.pkl` | **atom decomposition**: 39,033 eqs, 39,277 distinct atoms, all deg ≤2, ≤3 vars |
| `porder.py` → `order.pkl` | intra-equation atom order ⇒ the 39,277 atoms form a **single Hamiltonian path** |
| `pslp.py` → `slp.pkl` | orientation + forward evaluator (reproduces 38,744/38,748 vars of the deliverable) |
| `pzero.py` → `wire3.pkl` | P-alias class, vars ≡0 mod P (handles), mod-P equality classes |
| `pextract.py` → `blocks.pkl` | template match of the 383 law-blocks |
| `pgraph.py`/`ptopo.py` → `graph.pkl`,`topo.pkl` | stage wiring + topology |
| `pfold.py` | **independent fold evaluator** (validated: 376/382 stages match the deliverable) |
| `paudit.py` | the audit checks |

Facts: every equation is `scalar · L^k = 0`, k∈{1,2,4}; L is a linear combination of 3..24
atoms. The atoms tile one straight-line program of 39,277 gates; equations are contiguous
windows of it. 38,748 vars.

## 2. What the instance is (my derivation)
* One variable is pinned to a 256-bit constant **P** = `115792089237316195423570985008687907853269984665640564039457584007908834671663`
  (atom at SLP position 19242; 220 alias copies). "Handle" variables enter only as `P·u`,
  so **almost every constraint is a congruence mod P**. 7,323 vars are provably ≡0 mod P.
* A second variable is pinned to a 256-bit constant **Q** = `97553848499418123410591666447050222001188385549510401465815187079080512838891`
  (var `x_24453`, SLP position 822). It occurs in **383 atoms, all of the identical shape
  `E = D + Q`**, spaced *exactly 43 apart* in the SLP (positions 2578 … 19004).
* **256 selector variables**, each boolean, each attached to **exactly 2 coordinate
  variables**, each of those pinned to **1 of 512 distinct 287–296-bit constants**:
  `s·(x − K) ≡ 0` and `(1−s)·x ≡ 0` (mod P) ⇒ `x ≡ s·K`. Extracted by pure regex from the
  raw file: 512 triples, `B→#K = {1:512}`, `A→#B = {2:256}`. **256 leaves confirmed.**
* Two final congruences pin the root pair to a **target** `T = (T_x, T_y)`.

## 3. THE LAW — confirmed, and stated exactly
Each of the 383 blocks computes, from six inputs `i1..i6`:
```
A = i1 - i2          B = i4 - i3          E = i1 + i2 + i5 + Q
N1 = E·A² - B²                    (degree 3)
N2 = A·(i3 + i6) - B·(i2 - i5)    (degree 2)
```
and imposes three congruences `c_k1·N1 + c_k2·N2 ≡ 0 (mod P)`, k=1..3, with a
**different random 3×2 integer matrix in every block** (382 distinct) but rank 2, hence
equivalent to `N1 ≡ N2 ≡ 0`.
With `X=(i2,i3)`, `Y=(i1,i4)`, `Z=(i5,i6)` and `λ=(Yy−Xy)/(Yx−Xx)`:
```
Zx = λ² − Xx − Yx − Q ,   Zy = λ(Xx − Zx) − Xy
```
* **382/383 blocks matched the template byte-for-byte** (`pextract.py`), with identical
  signs (`Q_off=+1, sN1=+1, sN2=−1`) — **zero exceptions**. The 383rd (q=19004) is the root;
  it obeys the same law with the two operands swapped (the law is commutative), verified
  term by term by hand from the SLP dump.
* **Invertible: 300/300** random triples — given `X` and `Z`, `Y` is recovered in closed form
  (`λ=(Zy+Xy)/(Xx−Zx)`, `Yx=λ²−Xx−Zx−Q`, `Yy=Xy+λ(Yx−Xx)`).
* **Commutative 200/200.** **Associative 300/300 *on the leaf set*** (0/400 on random pairs —
  associativity holds only on the locus below).
* All 256 leaf pairs `(x,y)` (in one consistent coordinate order) satisfy the *same* relation
  `y² − (x + Q/3)³ ≡ b (mod P)`, `b = 64019533680030876408443198762210829058751700634554282185987325820393598524794`.
  256/256, one common value. The law is **closed** on that locus (300/300) and the **target
  also lies on it**. So the fold is a sum in an abelian group.

## 4. The mux — confirmed
Every block ends with two 3-way multiplexers `out = (1−a)b·X + a(1−b)·Y + ab·Z` over the
liveness bits, and computes `ab−a−b` (= −(a OR b)) as the output liveness.
**"a gate fires only when both inputs are live; with one live input the value passes
straight through" — CONFIRMED, 381/381.**
Liveness is **fully determined by the 256 selectors**: of the 762 liveness slots,
607 alias the child's OR output, 127 are pinned to 0, 28 alias a dead subtree's OR.
**No extra boolean freedom: the configuration space is exactly 2^256.**

## 5. WHERE I DISAGREE WITH F — counts and depth
| quantity | F / FLEET.md | **measured by me** |
|---|---|---|
| stages | **96** | **383 law-blocks** = 255 merges + 101 pass-throughs + 27 dead |
| depth | **6** | **9** (178-side depth 8, 78-side depth 7, + root) |
| leaves | 256 | **256 ✔** |
| "72 fully-determined + 24 leaf-adjacent" | — | 165 stage⊕stage, 89 leaf⊕leaf, 77 leaf⊕dead, 23 stage⊕dead, 27 dead⊕dead |
| root split 178 \| 78 | 178 \| 78 | **178 \| 78 ✔** (stage 381 support 178, stage 380 support 78) |
| law + universal constant | one, zero exceptions | **✔ confirmed 383/383** |
| invertible | yes | **✔ 300/300** |
| 2^256 − 1 configurations | yes | **✔** (2^256 selector settings; empty set leaves the root not live) |

The **255 merge nodes over 256 leaves** is exactly a binary tree, so F's picture is right in
kind; the numbers 96 and 6 are not. F apparently never resolved the 56 "undecoded slot
pairs" — **I resolved all 27 of mine: they are the outputs of the 27 dead stages, provably
≡0 mod P (their leaf support is empty). Nothing in them breaks the tree model.**

## 6. THE FINDING THAT MATTERS — the law is VACUOUS on the diagonal
If a merge stage receives **two equal live input pairs** (`X == Y`), then `A = B = 0`, so
`N1 = E·0 − 0 = 0` and `N2 = 0·H − 0·J = 0` **identically**, and the stage's output
`Z = (i5,i6)` is **completely unconstrained mod P**. Everything above it is then free, and
the root can be driven to the target by inverting the law up the chain (closed form, §3).

**The 39,026 deliverable is a live instance of exactly this.** I measured it:
* its fold matches my evaluator at **376 of 382 stages**;
* the 6 mismatches are the chain `277 → 330 → 357 → 370 → 377 → 380` (supports 3,6,10,21,37,78);
* it forces `stage380 output == stage381 output` (I checked: both equal
  `(82007976…230357, 37841415…869714)` mod P), which makes the root free, and then sets the
  root to `T` — the root congruences are satisfied.
* The cost of forging that equality is exactly its 4 corrupted variables
  (`x_642, x_28730, x_29854, x_31864`) and its **7 failing equations**.

**Consequence for the campaign:** the satisfying set is *not* only "subsets whose group sum
is T". It is
> `{ S : Σ_{i∈S} L_i = T }  ∪  { S : some merge stage sees two equal live inputs }`,
and the second family is a *collision* condition available at each of the 255 merge stages.
No optimality/"coding optimum" argument in this lab accounted for it.
Also **excluded**: any S for which an intermediate sum has no affine representative
(two inputs with equal x but unequal y ⇒ `N1 = −B² ≠ 0`, infeasible) — that is a hard
constraint on the search, not a freedom.

## 7. Honest limits of what I proved
* Everything above is **mod P**. The full problem is over ℤ; the handles absorb the
  P-multiples, but each congruence also carries a *small* multiplier `c` (e.g. 12630605) on
  the handle side, so the real condition is `c·P | R`, not just `P | R`. I have **not** built
  the integer lift, so "solve the subset-sum ⇒ full solution" is a **conjecture**, not proved.
  (Evidence for it: the deliverable is such a lift for its own configuration.)
* I did not search for a collision; I only established that collisions are sufficient.
* Knob set for every "determined/free" statement above: **the 256 selector variables**,
  with liveness derived; configuration: **arbitrary** unless stated (the 376/382 measurement
  is at the deliverable's configuration, ON-set = leaves {21, 167} in my index order).

## 8. Best verified score
**39,026 / 39,033** — the existing deliverable, re-verified. I did not beat it.

## 9. Highest-EV next experiments
1. **Build the integer lift** (`config → full 38,748-var assignment`). Choose every
   coordinate variable's residue in [0,P) and tune the lift `r + kP` so each `c·P | R`
   divisibility holds (CRT over the small `c`'s). This converts the problem into *pure*
   subset-sum and gives a one-shot emitter for whoever finds the subset. It is also the
   only way to test whether a non-target configuration can beat 39,026.
2. **Collision search at the root** — MITM over the 78-side (`2^78` sums) against the
   178-side. Infeasible head-on, but the 78-side is the smallest such handle in the
   instance and worth profiling before dismissing.
3. **Defect accounting**: measure, for each atom, how many equations it appears in and
   whether pairs of nonzero atoms can cancel inside an equation. The deliverable achieves
   7 failing equations from 7 nonzero atoms; establish whether a 2-nonzero-atom forgery
   (e.g. corrupt only the root pair) yields fewer than 7 failing equations. **This is the
   only route I see to beating 39,026 without solving the subset-sum.**

## 10. New measurement: the defect-window screen (why 7 is hard to beat)
The 39,277 atoms lie on a single path, and **36,079 of 39,033 equations use an exactly
consecutive run of atoms** on it. Each atom sits in 13.1 equations on average (min 1, max 22).
The deliverable's 7 nonzero atoms are the *consecutive* SLP positions **36291–36297**, a
low-coverage pocket (per-atom coverage 6–7); they touch only **12** equations, of which the
7 defect values cancel **5** — leaving exactly the 7 failures.

Sliding-window screen `min_p ( |equations touching [p,p+w-1]| − w )` over all p ≥ 1700:
the minimum is **4**, for every w from 1 to 12, attained at **SLP position 5497**.
So an *optimistic* floor for a window-localised defect is 4 failing equations
(→ 39,029), versus the deliverable's 7. The screen is optimistic: it assumes all w atom
values in the window are independently settable, which the deliverable's own window shows
they are not (it got 5 cancellations from 7 atoms, not 7). Still, **the region around SLP
position 5497 has never been probed by this lab and is a strictly better carrier than the
deliverable's 36291–36297 pocket on this measure.** That is experiment #3 in §9, made concrete.

Caveats, stated per the standing rules: this screen was computed over **all 39,277 atoms**
(knob set = every atom, unrestricted), and it is **configuration-independent** — it is a
property of the equation/atom incidence alone and says nothing about which defect values are
algebraically reachable at any selector configuration.

---
# CHECK-IN 8 — the 5497 carrier is dead. Prediction withdrawn.

## 11. PRIMARY RESULT: I withdraw the "deficit 4 at SLP 5497 → 39,029" prediction.
It was never a score and it is not even a candidate. **Measured, not argued:**

* SLP position 5497 is the atom `x29741 − x17440 + x27926` — the `A = i1 − i2` wire **inside
  law-block index 68** (`+Q` gate at 5502, leaf-support 2, both inputs leaves). Footprint 5.
* Perturbing it changes that block's `N1` and `N2`, which forces the block's **three
  law-congruence atoms at SLP 34872 / 34874 / 34876** (footprints 13 / 12 / 11) to become
  nonzero. The realizable defect is therefore **4 atoms touching 22 equations**, not 1 atom
  touching 5.
* The only way to avoid that is a configuration where block 68's gate is off — and then the
  perturbation never reaches the root, so it buys nothing toward the target while still
  costing its 5 equations.

**So the footprint screen is dead**, and agent C's check-in-3 conclusion stands: no function
of incidence alone prices a defect. My screen's error term is exactly the *coupling* — an
atom's true cost includes the congruence atoms its perturbation forces downstream, which is
residue content, not incidence. `min_p(|eqs touching [p,p+w−1]| − w) = 4` remains true as an
incidence statement and is **useless** as a cost.

## 12. Coordinator item 3 (run first, as instructed): the 5497 window is NOT a decoy
K's decoy explanation does **not** apply here, so this is a third failure mode, not K's.
* In my parse, **1,158 atoms have footprint < 7**, of which **1,152 are idempotency atoms**
  (`x − x²`; 1,145 of them footprint 1) and 6 are not. So K's *kind* of claim is corroborated
  at much larger scale — low footprint overwhelmingly means idempotency decoy. (K's count of
  12 vs my 1,158 is a decomposition difference; atom indices/counts are not comparable across
  directories, as you noted.)
* But **every atom in the 5490–5511 window is genuine law-block arithmetic** — products,
  squares, three-term linear combinations. `isidempotency = False` for all of them. The
  window sits in the 6-atom non-decoy remainder.
* Net: the screen found a *real* low-footprint atom and still failed, because of coupling.
  Both explanations are needed; neither subsumes the other.

## 13. Coordinator item 2 — corrected placement cost table (short version)
I did **not** redo the enumeration K already did. What I did compute, because it was already
built: for all 382 blocks, the 3 law-congruence atoms and their handle-defs, and the union of
equations touched by breaking 2 of 3.
* Cheapest **live** merges: block 279 (`S,S`, supp 3) touches **10**; block 2 (`L,L`, supp 2),
  block 151 (`L,L`), block 193, block 311, block 330 all touch **11**; the deliverable's own
  site touches **12** with 7 atoms. Worst is 39.
* But a live merge corrupted to hit the target directly gives **three** broken congruences,
  not two: with `Z` forced to `T`, `N1` and `N2` are determined and no rank-2 combination
  `c_k1·N1 + c_k2·N2` vanishes. That is why the deliverable instead pays for the *vacuous*
  route. **Answer to your question: 7 is the price of the placement, not of the degeneracy
  as such** — but I did not find a placement priced below it, and per your instruction I
  stopped rather than enumerate further.

## 14. On K's unreachability argument — corroboration, and one specific hole
K's mechanism and mine agree, reached from different decompositions; I am glad to have it
corroborated. On the impossibility half I can neither confirm nor refute it, because it is
stated in a vocabulary my parse does not produce. Two things I can say precisely:

1. **It rests on a premise about how the instance was built** (that the 256 leaves correspond
   one-to-one to distinct exponents in an order-N arithmetic). Nothing in my parse yields
   that; my parse yields only 256 constant pairs and a degree-3 identity. That premise is the
   withdrawn lens, so the argument inherits whatever caution attaches to it.
2. **A concrete hole worth someone's time.** The step "`|x − y| < 2ⁿ < N`, so the difference
   cannot vanish" and "at the root `x − y = ±N`" both need the relevant modulus to *exceed*
   the largest signed subset difference, so that only one wrap (`k = ±1`) has to be excluded.
   If the modulus that actually governs coordinate-pair equality is smaller than 2²⁵⁶, the
   enumeration is over `x − y = kN` for several `k`, and a carry walk covering only `k = ±1`
   is incomplete. **Whoever owns that argument should state which modulus bounds the walk and
   check that it exceeds 2²⁵⁶.** That is the one place I would push.

## 15. Adjudication with S: in my decomposition it is subset selection, not one-hot
Decisive and cheap to re-run: **zero atoms in the entire file touch two or more distinct
selector variables.** Each of the 256 selectors appears in only **5–6 atoms**, all local to
its own coordinate load and its own liveness fan-out. There is no cardinality atom, no
one-hot tie, no cross-selector coupling anywhere. So the residual I am describing is
**free independent subset selection over 256 booleans, 2²⁵⁶ configurations.**
If S's decomposition saturates to one-hot, S and I are describing different objects, and this
one-line test (`any atom containing two selector vars` → none) is the thing to run against
S's parse to find out which.

## 16. Vocabulary and scope corrections to §3 above
Restated neutrally, with no reading of what the instance "is": each block imposes the two
polynomial identities `E·A² − B² ≡ 0` and `A·(i3+i6) − B·(i2−i5) ≡ 0 (mod P)` with
`A = i1−i2`, `B = i4−i3`, `E = i1+i2+i5+Q`. The properties I verified are properties of
**these identities**: solving for `(i5,i6)` given the other four is a closed-form rational
map (300/300); the map is symmetric in its two operand pairs (200/200); iterating it is
order-independent **on the 256 constant pairs** (300/300), all of which satisfy the common
identity `y² − (x + Q/3)³ ≡ b (mod P)`, as does the target. **I draw no conclusion about
solvability from any of this, and my own results point the other way** — the degeneracy is a
*second* family of satisfying assignments, not a barrier.

## 17. Standing caveats, unchanged
* Everything is **mod P**. The integer condition is `c·P | R`; the lift is **still unbuilt**,
  so "solve the residual ⇒ full solution" remains a **conjecture**, evidenced only by the
  deliverable being such a lift for its own configuration.
* Knob set for every determined/free claim: **the 256 selector variables**, liveness derived.
  §11–13 measurements are configuration-independent incidence/coupling facts; the 376/382
  fold agreement is at the deliverable's configuration (ON-set = leaves {21, 167} in my index).

## 18. Best verified score — unchanged, and I did not beat it
**39,026 / 39,033.** I produced no new assignment this check-in. The carrier I predicted does
not exist, and I am reporting that as plainly as I would have reported a win.

---
# CHECK-IN 11 — the lift is NOT free. L's property is false in my parse.

## 19. Item 1, tested directly. **0 of 3,707 handle variables appear in exactly one atom.**
Handle = a variable defined by an atom of the form `h − (P-alias)·u`. My parse finds **3,707**
of them. Occurrence count, exhaustive, no exceptions:

| variable | count | atoms it appears in |
|---|---|---|
| handle `h` (the P-multiple) | **3,707** | **exactly 2 — every single one** |
| cofactor `u` (the free multiplier) | 3,707 | **exactly 1 — every single one** |

So L's load-bearing property **is false of the handles** and **true of the cofactors**. Given
L reports 3,681 against my 3,707 — a decomposition difference, not a disagreement about the
file — I think the likeliest reading is that **L counted `u` and concluded about `h`.**

## 20. Why "appears in exactly one atom" does not imply free
The two atoms a handle sits in are
```
  h − u·P = 0        (definition)      ⇒  P | h
  R − c·h = 0        (congruence)      ⇒  h = R/c
```
`u` does appear exactly once, and that is necessary for freedom but **not sufficient**: you
must also be able to *solve that atom over ℤ* for the variable, and `h − u·P = 0` is solvable
for `u` only when `P | h`. Composing the pair gives exactly the condition I have been
flagging, **`c·P | R`**, which is strictly stronger than the mod-P statement `R ≡ 0 (mod P)`
whenever `c > 1`.

## 21. How much is actually free — measured
Splitting all 3,707 handles by the multiplier `c` on `h` in the congruence atom
(the def-side coefficient is `±1` for all 3,707, so it contributes nothing):

* **2,780 handles have `c = 1`.** For these the integer condition collapses to `P | R`, which
  *is* the mod-P congruence. **These are genuinely free — L is right about 75% of them.**
* **927 handles have `c > 1`.** For these `c·P | R` is a real extra integer condition that
  mod-P reasoning does not deliver.

**So the reduction is about three-quarters closed and one-quarter open. It is not
unconditional, and per your "if and only if" I am not stating it in unconditional form.**

## 22. Empirical corroboration that the def atom is a real constraint, not an absorber
All four variables the 39,026 deliverable corrupts — `x642, x28730, x29854, x31864` — **are
handle variables**, and they are **exactly the four (of 3,707) for which `P` does not divide
the value**; the other 3,703 all satisfy `P | h`. Its entire 7-equation deficit is the two
atoms of a handle failing to be simultaneously satisfiable. If handles absorbed freely, that
deficit would not exist. Note the scope: this shows the pair is a *joint* constraint; it does
**not** show an honest configuration must fail it — the deliverable broke those handles on
purpose, to plant the degeneracy.

## 23. What would actually close it, and the number that worries me
Not an argument — a construction. Build the integer lift for **one** configuration and verify.
The available freedom is the integer lift of each free coordinate variable (`r + kP`), which
moves `R/P` and so can tune `R/P mod c`. The count to check first:
**~766 lift parameters (2 per law-block) against 927 conditions with `c > 1`.** Fewer knobs
than conditions is not fatal — the conditions are moduli of ~7 bits and one parameter can
serve several by CRT if the coefficients are invertible — but **nobody has counted the rank**,
and that rank is the whole question. The one encouraging data point: the deliverable satisfies
all 927 of them at its own configuration, so a lift exists for at least one configuration.

## 24. Status of my standing caveat: **unchanged, now quantified**
Everything I have reported remains mod P. "Solve the residual ⇒ full solution" is **still a
conjecture**, and I now know precisely which part is unproved: the 927 `c > 1` divisibilities.
Knob set unchanged (the 256 selectors, liveness derived). Best verified score unchanged at
**39,026 / 39,033**; I produced no new assignment.

---
# CHECK-IN 13 — I did not count the rank. Reporting the failure plainly.

## 25. Item 1: **NOT ANSWERED.** The rank is still uncounted.
I did not deliver the measurement I named. Stating it without softening: the quarter of the
reduction that was open before this check-in is **still open**, and I have not moved it.

What I built instead was the machine that would let me count it (`plift2.py`): seed the 256
selectors, the 512 leaf coordinates and the 764 block law outputs, then walk the SLP with a
worklist, solving each atom for its single remaining unknown **over ℤ** and recording every
division that does not come out exactly. Counting the rank needs the derivative system
`∂(R/P)/∂t_v mod c` over the 927 `c > 1` conditions, and that needs a *complete* lift to
differentiate around. I do not have one.

## 26. Item 2: partial lift, and the obstruction is **in my tooling, not in the instance**
This distinction is the whole point of reporting it, so I am putting it first: what stopped
me is a limitation of my propagation, **not** a property of the file. Nothing below is
evidence about whether the instance's lift is free.

| configuration | undetermined vars | integer-division obstructions | nonzero atoms | equations failing | score |
|---|---|---|---|---|---|
| all selectors OFF | 18,417 | **0** | 194 | 2,645 | **36,388 / 39,033** |
| one leaf ON | 18,417 | 18 | 234 | 2,815 | **36,218 / 39,033** |

**Both are far below the 39,026 baseline and neither is a competitive assignment.** They are
incomplete constructions in which 18,417 variables were never determined and defaulted to 0;
the failures they show are that default, not the instance resisting.

The stall is diagnosed but not fixed: the first nonzero atom is SLP 10834, the `+Q` gate of
block 192 (`x38494 = x11478 + Q`), and it is nonzero only because `x11478` was still
undetermined when the worklist drained. The propagation stalls wherever an atom holds **two**
unknowns at once — a handle and its cofactor, or a mod-P copy target and its handle — and my
seeding does not cover the copy targets. That is a fixable gap in `plift2.py`, roughly:
also seed every mod-P copy target to equal its source exactly over ℤ, which forces those
handles to 0 and unblocks the cascade. I ran out of budget before doing it.

## 27. The one real datum, with its scope stated
For the all-off configuration the constructor completed with **0 integer-division
obstructions** — every divisibility it was asked to perform came out exactly. **This is not
evidence that the lift is free.** That configuration drives essentially every residual `R` to
0, so it exercises the 927 `c > 1` conditions only trivially. It shows the constructor is
sound on a degenerate input; nothing more.

## 28. Status of the reduction: **still conditional. I am not stating it unconditionally.**
Item 1 did not come back clean — it did not come back at all — so the "if and only if" is
unmet and I am holding the line. Unchanged from check-in 11: 2,780 of 3,707 handles are
genuinely free at `c = 1`; **927 carry the strictly stronger integer condition `c·P | R`
whose satisfiability is unproved**, and the rank that would decide it is uncounted. Knob set
unchanged (256 selectors, liveness derived); everything else remains mod P.

## 29. Best verified score: **39,026 / 39,033**, unchanged, not mine.
I produced no assignment that beats or approaches it this check-in.

## 30. Next step, concretely
Fix the seeding gap in §26 (seed mod-P copy targets to their sources over ℤ), confirm a
complete lift at the all-off and one-leaf configurations, then differentiate: for each of the
927 `c > 1` conditions compute `∂(R/P)/∂t_v mod c` against the ~766 lift parameters and take
the rank modulo each prime power dividing the `c`'s. That is the measurement, and it is one
working lift away.

---
# CHECK-IN 16 — constructor complete. First rank measurement, on an admittedly wrong parameter set.

## 31. The seeding fix landed, plus a second fix I had not diagnosed
`plift5.py` is a working, complete lift constructor. Two changes over `plift2.py`:
1. **The seeding fix you quoted back to me**: when an atom holds two unknowns and one is a
   handle `h` (a variable defined by `h − (P-alias)·u`), set `h = 0`. That is the general form
   of "seed every copy target to equal its source exactly over ℤ" — it forces `u = 0` and makes
   the copy exact. Nonzero atoms dropped from 194 → 3.
2. **A second stall I had not predicted: ordering.** With a plain worklist a variable could be
   solved from a *downstream* constraint before its own definition atom fired (e.g. `x23927`
   was being back-solved out of the target congruence). Replacing the deque with a heap keyed
   by SLP position fixed it. **A third fix on top:** the `h = 0` rule must *not* fire at `h`'s
   own definition atom, or it pre-empts the very divisibility test I am trying to measure.

## 32. Verified end-to-end behaviour of the constructor
| configuration | undetermined | division obstructions | nonzero atoms | failing eqs |
|---|---|---|---|---|
| all off | 9,040 | 2 | 3 | 27 |
| one leaf ON | 9,040 | **2** | **2 — exactly SLP 39273 & 39275, the two target congruences** | 17 |
| two leaves, live merge at block 2 | 9,040 | 4 | 4 — one block-2 law-congruence pair + the two targets | 27 |

The one-leaf row is the clean one: **the only atoms that fail are the two target congruences,
nothing else.** That is the constructor behaving exactly as a correct lift should.
These are diagnostics, **not score attempts and not competitive assignments**; the remaining
9,040 undetermined variables are still defaulted to 0.

Cross-check against L's table: at `|S| = 2` L reports 4 distinct `c > 1` atoms violated; I get
one violated `c > 1` pair at the single live block plus the two targets. Same order, same
places, from unshared decompositions.

## 33. **I was wrong about the parameter count, and it changes the answer**
I told you ~766 lift parameters (2 per law-block). **That is an undercount.** The 512 leaf
coordinates are pinned only *modulo* P: the leaf atom is `s·(x − K) = c·h` with `h = u·P`, so
`x = K + c·u·P` is legal for any `u`. I had been fixing `x = K` exactly and treating leaves as
rigid. The correct parameter set is **~1,278** (766 block outputs + 512 leaf coordinates),
and at a leaf⊕leaf merge the local count is **6**, not 2.

## 34. The rank measurement I did get — scoped to what it actually covers
Per live block the system is closed-form. Writing `n1 = N1/P`, `n2 = N2/P`, shifting
`i5 → i5 + P·t5`, `i6 → i6 + P·t6`:
```
 d(n1) = t5·A²        d(n2) = −(t5·B + t6·A)        condition: c_k | c_k1·n1 + c_k2·n2
```
At block 2 the three moduli are `c = (1, 1, 7038713)`, and `7038713 = 11·23·43·647`, each
prime dividing exactly one modulus — so one equation in two unknowns per prime.
**Result: q = 11, 43, 647 solvable; q = 23 NOT solvable** (the shift directions are degenerate
mod 23). So with the 2-parameter-per-block model **the rank is not full.**

**But that measurement is on the wrong parameter set** (§33): it omits the four leaf-coordinate
lifts that this block's own inputs carry. With them, `A` and `B` themselves move and the mod-23
degeneracy very likely lifts. **So I am not reporting "the rank is deficient" — I am reporting
that the rank is deficient in a 2-parameter model I now know to be too small, and that the
6-parameter computation is the one that decides it.**

## 35. Handover — what remains between this and the rank, precisely
Everything below is one step, on top of `plift5.py`:
1. Extend the shift model to the full parameter set: for each live block, parameters
   `t5, t6` (step `P`) **and** the lifts of whichever of `i1..i4` are leaf coordinates
   (step `c·P`, with `c` that leaf atom's own multiplier) or are parent-block outputs.
2. Re-derive `d(n1), d(n2)` including the `A`- and `B`-motion terms — `A = i1 − i2` and
   `B = i4 − i3` are no longer constants, so `n1 = (E·A² − B²)/P` picks up terms in
   `t_{i1}, t_{i2}, t_{i3}, t_{i4}`; all still polynomial and exactly computable.
3. Solve prime-by-prime (`fac(c_k)`, small systems mod `q^e`, CRT across primes) — the
   loop in this check-in already does this correctly and runs in under a second once the
   moduli are factored. **Do not brute-force over `lcm(c_k)`; I burned a run doing that.**
4. Blocks couple: block `j`'s output lift shifts its parent's inputs, so the full system is
   block-coupled, not block-diagonal. For a configuration with one live merge it *is*
   decoupled, which is why `|S| = 2` is the right place to start.

## 36. Status of the reduction: **still conditional. Unchanged.**
2,780 of 3,707 handles free at `c = 1`; **927 carrying `c·P | R`, satisfiability unproved.**
L's independent 927 (from `c = 1` for 2,747 plus 7 zero-slope) matching mine raises my
confidence in the *count*, not in the *satisfiability*. The rank is still uncounted on the
correct parameter set, and I am not stating the reduction unconditionally.
Knob set: 256 selectors, liveness derived; everything else mod P.

## 37. Best verified score: **39,026 / 39,033 — unchanged, not mine.**

---
# CHECK-IN 20 — mod 23 lifts. The deficiency was an artifact of the small model.

## 38. Result
At **block 2, |S| = 2, with the full 6-parameter model**, all three `c·P | R` conditions are
**simultaneously satisfiable**. The mod-23 failure I reported last round was an artifact of the
2-parameter model, exactly as I suspected. It is **not** an obstruction in the instance.

Constructed, not argued:
```
 legal mu-steps (i1..i6)      : [4373213, 7633471, 1, 1, 1, 1]
 three conditions (c1,c2,c_k) : [(5788325,9395331,1), (9705029,4851321,1), (10233687,4279357,7038713)]
 7038713 = 11 * 23 * 43 * 647 ; root pairs (t1,t2) mod q: 20 / 22 / 41 / 645
 CRT solution                 : t1 = 0 , t2 = 383619
 shifted point still mod-P valid                : True
 condition 3 (c = 7038713) now divides exactly  : True
 conditions 1 and 2 have modulus 1              : vacuous
```

## 39. Two things that make me trust this one
1. **The verdict is confirmed by direct recomputation, not by my expansion.** I rebuild the
   shifted integers `i_k + P·mu_k`, recompute `N1, N2` from scratch, and check `P | N1`,
   `P | N2` and `c_k·P | R`. The expansion is only used to *search*.
2. **That check caught a real bug in my own expansion.** I had written the `n2` shift as
   `A*h2` where the algebra gives `B*h2`. The first run reported "condition 3 divides: False"
   from the direct recomputation while my expansion said the residual was 0 — the disagreement
   is what exposed it. After the fix the expansion is **verified exact against direct
   recomputation on 5 random shifts**, and the CRT solution verifies directly.
   Had I trusted the expansion alone I would have reported the opposite result last round.

**Robustness note:** my step derivation for `i3, i4` returned 1, which I could not confirm and
which would be *optimistic* if their true step is `c > 1`. It does not matter here — the CRT
solution has `t1 = 0` and uses only parameter `i2`, whose step (7633471) is properly derived
from its own leaf atom. The verdict rests on one correctly-derived parameter.

## 40. Scope — one block, |S| = 2. I am not generalising it.
This settles **block 2 at |S| = 2**, the configuration where I had established the system is
genuinely decoupled. It is one block out of 255 merges and one of the 927 conditions.
It says nothing about:
* whether the other 926 conditions lift;
* whether they lift **simultaneously** — at |S| > 2 several blocks are live at once and the
  parameters couple through parent inputs, so the per-block CRT above does not compose;
* larger `|S|`, which I did not touch this round.

A useful incidental: at block 2, **two of the three conditions have modulus `c = 1` and are
vacuous**; only one carries `c > 1`. That ratio is consistent with the global 927 / 3707 ≈ 25%.

## 41. Reduction status — **unchanged and still conditional.**
2,780 of 3,707 handles free at `c = 1`; **927 carrying `c·P | R`, satisfiability unproved.**
One condition settled favourably at one block is not 927 conditions settled, and I am not
softening the status on the strength of it. L's matching 927 raises confidence in the count,
not the satisfiability. Knob set: 256 selectors, liveness derived; everything else mod P.

## 42. Handover — the next step is coupling, not more single blocks
`prank.py` runs the whole computation end to end and is left in a state where the next agent
runs it rather than rebuilds it. To extend:
1. **Do not brute-force `lcm(c_k)`** (I burned a run on that) and **do not trust an expansion
   without the direct-recomputation check** (it cost me a wrong sign this round). Both guards
   are already in `prank.py`.
2. Generalise `step[]`: my leaf-multiplier lookup failed for `i3, i4`. Fix it before any
   configuration where those parameters are load-bearing.
3. The real open question is **composition**: at `|S| > 2` multiple blocks are live, block `j`'s
   output lift shifts its parent's inputs, and the system stops being block-diagonal. The right
   next measurement is `|S| = 3` or `|S| = 4` with two live merges in a parent/child relation —
   that is the smallest case where coupling actually bites.

## 43. Best verified score: **39,026 / 39,033 — unchanged, not mine.** No score attempted.

---
# CHECK-IN 25 — i3/i4 resolved; composition NOT settled. Handover.

## 44. The i3/i4 step weakness I flagged: **resolved, and it was never a lookup failure**
```
i3/i4 step audit at block 2:
   i1 = x30632   step=4373213    atom found=True
   i2 = x27436   step=7633471    atom found=True
   i3 = x10424   step=1          atom found=True
   i4 = x20930   step=1          atom found=True
```
`stepof()` **did** find the atom for all four; the y-coordinate leaf atoms genuinely carry
multiplier **1**. So step 1 is the correct value, not a default masking a miss, and
`stepof()` is correct as written — it now returns a `found` flag so the two cases can never be
confused again. My check-in-20 worry that this was "optimistic if wrong" is **discharged**:
it was right.

## 45. Composition: **NOT SETTLED. I did not get an answer.** Saying so plainly.
`pcompose.py` is written and finds the parent/child candidates, but the run did not complete.
The blocker is mine and it is mundane: **`conds()` re-scans all 39,277 atoms to look up each
condition's modulus on every single call, and it is called inside a `q²` inner loop.** Hoist
the three `(c_k1, c_k2, c_k)` triples out — they are per-block constants — and the run becomes
fast. Nothing about the instance blocked this; my inner loop did.

## 46. An unverified hypothesis I want on the record as a hypothesis
While wiring `pcompose.py` I noticed that a parent's input is joined to its child's output by a
**mod-P copy congruence** `x_a − x_b = c·h`, which itself carries a free lift of step `c·P`.
**If** that holds, the parent's conditions could be satisfied using the copy-edge lifts without
touching the child's, and composition would largely decouple — which would be the favourable
answer. **I did not verify this and it is not a result.** It is precisely what step 5 of
`pcompose.py` was built to test, and it is untested. Do not cite it as anything else.

## 47. HANDOVER — exactly what to run first
State of the tooling, all in `solve_lab/agentP_work/`:
* `plift5.py` — working integer lift constructor (seeds copy targets via the `h = 0` rule,
  processes in SLP order, does not pre-empt the divisibility test at a handle's own def atom).
  At one leaf ON, the **only** nonzero atoms are the two target congruences.
* `prank.py` — the 6-parameter rank computation, runs end to end, verified verdict at block 2.
* `pcompose.py` — composition test, written, **needs the §45 hoist before it will run**.

Two guards, both already in the files, both learned the hard way:
1. **Never brute-force over `lcm(c_k)`** — factor and go prime-by-prime, then CRT.
2. **Never trust a symbolic expansion without direct recomputation** — rebuild the shifted
   integers and recheck `P | N1`, `P | N2`, `c_k·P | R` from scratch. This caught a real sign
   bug in my own algebra (`A*h2` where it should be `B*h2`) and would have inverted my
   check-in-20 result had I skipped it.

**The one configuration to run first:** `cands[0]` in `pcompose.py` — the parent/child pair
with the smallest leaf support in which the child is a leaf⊕leaf merge and the parent's other
slot is live. That is the smallest configuration with two live merges in a parent/child
relation, and therefore the first one that can distinguish "each condition lifts individually"
from "they lift simultaneously". Simultaneity is what the reduction needs; individual lifting
is all anyone has shown, including me.

## 48. Reduction status — unchanged. Fifth time, and still not softened.
**2,780 of 3,707 handles free at `c = 1`; 927 carrying `c·P | R`; satisfiability OPEN.**
One block settled favourably at `|S| = 2` (check-in 20) is one of 927, in the configuration
where the system is decoupled — i.e. in exactly the case that cannot test simultaneity.
Knob set: 256 selectors, liveness derived; everything else mod P.

## 49. Best verified score: **39,026 / 39,033 — unchanged, not mine.** No score attempted.

---
# CHECK-IN 28 (final) — composition holds at cands[0], jointly verified.

## 50. Result
`pcompose2.py` (hoist applied: all per-variable lift steps and per-block condition triples
computed once — 9,164 and 382 of them — instead of re-scanning 39,277 atoms inside the loop).

```
cands[0] : parent block 193  <-  child block 2 ,  |S| = 3 ,  2 live blocks
child  conditions : [(0,1), (0,1), (6982445, 7038713)]      -> 1 non-vacuous, c = 11*23*43*647
parent conditions : [(0,1), (1449394, 1599077), (0,1)]      -> 1 non-vacuous, c = 59*27103
child  lift steps : [4373213, 7633471, 1, 1, 1, 1]
parent lift steps : [1, 4373107, 1, 8854455, 1, 1]   (slots (1,2) carry the child's output)

child  conditions solvable with the child's own parameters   : True
parent conditions solvable with the copy-edge lifts ALONE    : True
parent conditions solvable with all six parent parameters    : True
```
**And then jointly verified rather than inferred:**
```
child  CRT shift vector : [6297909, 1145837, 0, 0, 0, 0]   mod 7038713
parent CRT shift vector : [0, 0, 1229442, 0, 0, 0]         mod 1599077
both shifts applied SIMULTANEOUSLY, recomputed directly from the shifted integers:
  CHILD  : P|N1 True   P|N2 True   c*P|R True
  PARENT : P|N1 True   P|N2 True   c*P|R True
  parent shift touches only slot 2 (a copy edge) -> child variables untouched : True
  JOINT VERIFICATION : True
```
So at this pair, **individual lifting and simultaneous lifting both hold** — the first data
point in this campaign that distinguishes the two, and it is favourable.

## 51. Why I ran the joint check, and what it changes
My first verdict line was a logical AND of two separate searches — an *inference* from
parameter disjointness, which is exactly the shape of thing my second guard exists to catch.
So I CRT'd each block's per-prime root vectors coordinate-wise, applied **both** shifts at
once, and recomputed `N1, N2` and the divisibilities from the shifted integers. It passed.
Had it not, the disjointness reasoning would have been wrong and I would have reported the
opposite. The verdict above is verified, not inferred.

## 52. My check-in-25 hypothesis: **confirmed at one pair. Not established.**
The mechanism is the one I guessed and explicitly refused to claim: the parent's input is
joined to the child's output by a mod-P copy congruence carrying its own free lift, so the
parent's condition is discharged by a copy-edge lift that never touches the child. That is now
**observed**, at one pair, with the child verified untouched. It is not proved in general and
should be cited as "confirmed at cands[0]", nothing broader.

## 53. Scope — I am not generalising this either
* **One** parent/child pair, `|S| = 3`, two live merges.
* Only **2 of the 6** conditions present in this configuration are non-vacuous; four have
  `c = 1`. So this settles **2 of the 927**.
* Deeper chains (grandparent/parent/child), sibling live merges feeding one parent, and every
  `|S| > 3` remain untested. Nothing here says the copy-edge lift is always available or always
  independent.

## 54. Reduction status — unchanged. Sixth time.
**2,780 of 3,707 handles free at `c = 1`; 927 carrying `c·P | R`; satisfiability OPEN.**
Two favourable conditions at one pair is not 927. Knob set: 256 selectors, liveness derived;
everything else mod P.

## 55. Final state of the durable output
* `plift5.py` — working integer lift constructor; at one leaf ON the only nonzero atoms are the
  two target congruences.
* `prank.py` — 6-parameter rank computation, verified at block 2, `|S| = 2`.
* `pcompose2.py` — composition test, **hoisted and runnable**; reproduces §50 end to end.
* Two guards, both learned by being bitten: **never brute-force over `lcm(c_k)`** (factor,
  prime-by-prime, CRT); **never trust a symbolic expansion or a disjointness argument without
  direct recomputation** — it caught a sign bug at check-in 20 and it is what turned §50 from
  inferred into verified.
* Next configuration to run: `cands[1..3]` (all support 3) to see whether the copy-edge
  independence repeats, then the first three-deep chain.

## 56. Best verified score: **39,026 / 39,033 — unchanged, not mine.** No score attempted.
