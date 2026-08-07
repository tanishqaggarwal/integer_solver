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
