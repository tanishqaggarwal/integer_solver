# RESUME_S — agent S (lattice methods), self-contained handoff

## 0. Verification rules
- Baseline `../best/new_instance_partial_39026.json` = **39,026/39,033**, re-verified by me with
  plain `solve_lab/checker.py` (failing `[12231,12270,12350,14584,18673,22044,29125]`).
- States with >4,300-digit integers must be verified with `../agentE_work/verifyE.py`, NOT
  `checker.py` (which raises ValueError on *parse*). Always say which was used.
- **Every "cannot move X" claim below states its knob set AND its selector configuration.**
  I found and corrected exactly the failure the brief warned about: agent E's channel
  measurement filtered to boolean knobs, which hid the single-row knobs. See §2.

## 1. Setup (rebuild in ~5 s)
Work from `solve_lab/agentS_work/`. Symlinks to E's `orient.pkl`, `users.pkl`,
`triple8_seed.json` are already here; `common.py` loads agent E's forward engine read-only.
    python3 -c "import sys;sys.path.insert(0,'.');import common"    # ~3 s
| script | what it does |
|---|---|
| `meas.py` | unfiltered delta table, all cone knobs x all 5 cluster rows |
| `table.py` | exact delta supports over ALL atoms -> `table_cfg0.pkl` |
| `handles.py`/`handles2.py` | global pure-handle census (structural prefilter + exact step) |
| `hcheck.py` | per-atom pure-handle search inside that atom's cone |
| `local1.py`,`local2.py`,`local3.py`,`fix1.py` | the deliverable's own local system |
| `cascade.py` | exact cascade repair on a full assignment, scored at equation level |
| `sat.py` | saturation test measuring atom values directly |
| `image.py`,`bfs.py`,`bfs2.py` | reachable image of the selector->residual map |
| `lattice.py`,`lat2.py`,`lat3.py`,`lat5.py` | the exact integer lattice solve |
| `basin2.py`,`basin3.py`,`sacrifice.py` | minimum-equation-cost sacrifice-set search |
| `final.py` | configuration-dependence re-check of the a20215 step |
| `reach.py` | is the BFS closure local or global? (ANSWER: local) |
| `reach2.py` | can affine knobs move a20215 mod p? (ANSWER: yes, by +-1) |
| `reach3.py` | joint solve at random configs — INCONCLUSIVE, see sec 8.3 |
| `kernel.py` | sec 8.3: displace along the affine kernel, re-measure the obstruction |
| `kernel2.py` | sec 8.3 complement: obstruction at each BFS image point (moves the mod-p class) |
| `trade.py` | trade-knob walk — VACUOUS, checks span membership up front (sec 6g) |
| `relax.py` | relaxed selectors — genuine structural move, other rows infeasible (sec 6h) |
| `structural.py` | row-deficiency analysis over all 72 logged configurations (sec 6i) |
| `structural2.py` | finds off-shape full-row-rank configurations — the valid test cases (sec 6i) |
| `dirsearch.py`,`combine.py` | deficiency-directed sweep; post-solve-class independence (sec 6j) |
| `degen2.py` | degeneracy discriminator: zero AND unresponsive |
| `selcouple.py`,`selcouple2.py` | selector-coupling census + classification |

## 2. MAIN RESULT — the endgame condition, derived exactly (not a filtered barrier)
**1,815 free variables are pure handles**: each moves exactly one atom, affinely, by exactly
**±p** (p = 2^256 - 2^32 - 977) — `handles2.py`. An atom is therefore satisfiable iff its
residual is ≡0 modulo its handle step. Measured steps at cfg0: a20215 p, a20212 p (also
12354891p), a7389 p, a10187 **7942211*p**, a747 p, a32148 p, a28035 p, a28037 p,
a20652 6672769*p, a28033 8640431*p. **a28647, a30787, a26958, a40306, a726 have NO handle.**
The trade knobs are exactly 1-for-1 and were the ones E's boolean filter excluded:
x_14853 (a20212 −1, a28647 +1), x_6083 (a7389 +1, a28647 −1), x_31339 (a10187 +1, a20215 −1),
x_18956 (a747 +8863713, a20215 +1), plus the pure handles x_22820, x_11436, x_26489, x_37012.

`lat3.py` at cfg0 (= `triple8_seed.json`, selectors x_1530 = x_1603 = 1; residual a20215, a28647):
build `D n = -R` over the **complete affine knob set** — 54 knobs, every affine knob in the
5-row cluster cone *including every single-row knob*, plus the pure handle of every atom in
play, plus the cone free vars of the no-handle atoms. 49 atoms in play, only a726 and a28647
unhandled. Eliminate the two target rows, solve the other 47 exactly over Z (FEASIBLE, kernel
dim 7), then compute the image lattice on (a20215, a28647):

> **the reachable lattice is exactly p·Z², and the required offset is NOT in it.**

So the endgame condition at cfg0 is exactly `a20215 ≡ 0 (mod p)` **and** `a28647 ≡ 0 (mod p)`,
and both are nonzero (…263210356252 and …695745279353).

## 3. ⚠️ RETRACTED — I refuted this myself. Read this section before using anything below it.
**Retraction 1 — the closure is LOCAL, not global** (`reach.py`). I described the §3 BFS as
"terminated by exhaustion". It is exhaustive only over the closure reachable *from cfg0 under
single flips*. I tested it against configurations it could not reach — random subsets of every
weight, each with all 256 selectors set explicitly rather than a few flips off cfg0.
**148 of 300 landed OUTSIDE the 48-tuple image, producing 14 new tuples.** Even |S| = 1 landed
outside (cfg0 has x_1530 = x_1603 = 1, so "only this selector on" is 3 flips away, not 1).
So "the image closed at 48 tuples" does **not** bind the instance.

**Retraction 2 — "a20215 mod p is never 0" is true but is NOT an obstruction** (`reach2.py`).
p is prime, so any knob moving a20215 by a step ≢ 0 (mod p) makes every residue reachable.
At cfg0 **and** at random |S| = 17/64/128/200 configurations outside the closure, exactly 3 affine
knobs move a20215 and **2 of them move it by ±1**: x_18956 (step +1) and x_31339 / x_30213
(step −1). So a20215 ≡ 0 (mod p) is reachable *on its own*, trivially. My §3 measured the image of
the **selector** map with affine knobs held fixed, and I over-read that as an obstruction.

**What actually survives is §2, and only in its joint form.** x_18956 moves a20215 by 1 but pays
8863713 into a747, whose only handle steps by p; keeping a747 satisfied forces
8863713·n ≡ 0 (mod p), hence n ≡ 0 (mod p) since gcd(8863713, p) = 1 — so a20215 moves by
multiples of p *only once the other rows must stay satisfied*. That is exactly the p·Z² image §2
computed. §2 and §3 were never independent; §3 was a weaker restatement that dropped the "other
rows stay satisfied" clause, which is the entire content.

**Consequence for the tension with Q's existence result: there is no contradiction, and the fault
was mine.** My result is a cfg0-*local* joint statement, the reachable space is demonstrably
larger than I mapped, and nothing I measured forbids a satisfying assignment elsewhere.
Do not cite §3 against any existence claim.

## 3 (original, now scoped down to what it actually shows). Half the condition is reachable
`bfs.py` — BFS over configurations, move set = **all 256 cluster booleans + the switch x_30163**,
non-boolean affine knobs at cfg0 values (justified: §2 shows they act only in p·Z² on the pair).
It terminated with no new tuples under single flips, giving 48 distinct mod-p 5-tuples.
**This is a local closure only — see the retraction above.**

| row | distinct values mod p | reaches 0 ? |
|---|---|---|
| a7389 | 4 | **yes** |
| a10187 | 4 | **yes** |
| a20212 | 2 | **yes** |
| **a20215** | **2** | **NO** |
| a28647 | 4 | **yes** |

`bfs2.py` (same, plus the second switch x_11559) reached 94 tuples before its timeout and found
configurations with **four of the five rows simultaneously 0 mod p** — but a20215 mod p was
still only ever 44859544763832475231923253825569092119321525945631045653619508440821028887
or 22981624690591324143788809642515852940280603493270692712106986169263210356252.
`final.py` re-measured at two such configurations: a20215's handle step is still exactly p there,
and the other rows go infeasible.
**Correct reading of all of the above:** a20215 mod p = 0 is the one thing the *selector* map never
reached with affine knobs held fixed. That is a fact about the selector map's image and nothing
more. `reach2.py` shows affine knobs move a20215 by ±1 at every configuration tested, so a20215
mod p = 0 is reachable outright; and `reach.py` shows the 48-tuple image is a local closure.
Scope of the BFS, stated plainly: move set = 256 cluster booleans + the switches, affine knobs
pinned at triple8_seed values, single flips from cfg0. **Not the whole space.**

### 3b. `lat5.py` corroboration (independent of the BFS mod-p reading)
**Status: 22 of 48 configurations completed, 0 feasible; a20215 was in the bad set in 22 of 22.** The remainder was sharded across three
`lat5p.py` workers but the box is heavily contended (7+ other agents' python jobs); the workers
were getting roughly 6% CPU each, so the sweep was left running rather than completed. Anyone
resuming should check `runs_lat5b.log` (configs 0-21, serial) and `runs_lat5p_{0,1,2}.log`
(configs 18-47, sharded) before re-running anything.

`lat5.py` re-measures knobs, handles AND targets at each configuration and runs the full exact
integer solve. Across every configuration it processed, **a20215 was in the bad set in 100% of
them** — and it was never the *reported* blocking row, because it does have a handle (step p) and
so is individually satisfiable; the infeasibility always surfaces at some other row first
(commonest: rows 21548, 11789, 26960, 10187, or "core infeasible"). This is a second, independent
route to the same conclusion as §3: a20215 nonzero is the invariant, not an artifact of reading
the image mod p.

## 4. SATURATION — reproduced, and STRONGER than previously recorded
`sat.py` measures the atom values directly (not bad-atom-dict deltas). At cfg0 the 219 moving
booleans fall into 2 classes (178 / 41) on (a20215, a28647) mod p, and **every tested subset —
within-class AND CROSS-class — reproduces the delta of a single member.** Cross-class saturation
is new. Both non-boolean "knobs" x_30163 and x_11559 are likewise **switches**: every nonzero
value, from 1 to 10^20, gives an identical delta. Consequence: the selector→residual map is a
priority / one-hot selection, not a sum. **It is not a subset-sum, so LLL on the bit vector is
the wrong tool** — an independent second reason the earlier LLL attempt failed.

## 5. The deliverable is a CANCELLATION CODE — its 7 is not 7 broken atoms
`local1.py`/`fix1.py`: the 39,026 assignment has 8 nonzero atoms whose footprints union to
**12** equations; 5 of the 12 are satisfied by cancellation. Its residual atoms are all
`x_A − p·x_private` handle atoms, with x_17499 = x_22665 = x_28599 = x_28961 = p exactly;
a36663 is `x_31864`, a36662 is `x_7075 * x_8731`.
Setting x_31864 = 0 and x_10903 = 0 zeroes a36661 and a36663 with **provably** no side effects
(occ = 2 and occ = 1) — and makes the score **worse: 11 fails / 39,022**. Depth-4 exact cascade
repair (`cascade.py`) found nothing better than 39,026.

**Objective note — CORRECTED, read the correction before using this.** 769 atoms have
equation-footprint exactly 1 and cover 769 *distinct* equations, so a state whose only nonzero
atoms are ≤6 footprint-1 atoms fails ≤6 equations and scores ≥39,027. That arithmetic is sound,
but I checked what those equations actually contain and it is **incidence structure only, not a
realisability claim** (exactly the caution in FLEET check-in 3):
- **0 of the 769 equations contain only one atom.** 768 contain exactly 2 atoms, 1 contains 11.
- So no footprint-1 atom owns its equation outright — each shares it with a partner atom.
- Cuts both ways, and the second way is the interesting one: because these are 2-atom equations,
  **cancellation is available there too**, so a nonzero footprint-1 atom can cost *zero*
  equations if its partner cancels it. That makes the target more attainable, not less — but
  nothing here shows any such state is reachable, and I did not reach one.

At cfg0 the minimum footprint of any atom in play is **10**, and the
trade graph does not reach a footprint-1 atom: the only knob linking the cluster to one
(x_11559 → a40306, nf 1) is a switch, not a continuous trade (`route.py`).

## 6. Sacrifice-set search (minimum-equation-cost residual placement)
`basin2.py`/`basin3.py` solve the exact integer atom system with a small set of atoms *dropped*
(allowed to carry residual), preferring tiny footprints. In the deliverable's own basin
(E-forward of its free vars: bad = a23618, a34120, a36660, a36661, a36662; 25 fails) the best
over all 1-, 2-, 3- and 4-atom drops is **39,019** (bad = a3939, a3940 — footprints 13 and 14
overlapping in 13 equations). Not competitive with the code's 7.

## 6b. DEGENERACY (P's family) — tested directly in my own parse, `degen2.py`
P's degenerate family = a merge sees two equal live inputs, A = B = 0, both congruences vanish
*identically*. The observable signature is NOT "the row is zero" (that is just the congruence
being satisfied) but **"the row is zero AND has stopped responding to every knob in its cone"** —
an identically-vanishing constraint cannot be moved. Measured over all 48 converged BFS
configurations x 333 cone knobs:

| row | zero & frozen (DEGENERATE) | zero & responsive | nonzero & responsive | nonzero & frozen |
|---|---|---|---|---|
| a7389 | **0** | 8 | 40 | 0 |
| a10187 | **0** | 8 | 40 | 0 |
| a20212 | **0** | 24 | 24 | 0 |
| a20215 | **0** | 0 | 48 | 0 |
| a28647 | **0** | 12 | 36 | 0 |

**Zero degenerate cases.** Every row is responsive at every configuration (responsiveness never
falls below 2, though it collapses from 222 to 3–4 — that collapse is the saturation of §4, not
degeneracy). a20215 is additionally the only row never *exactly* zero (0 of 48; the others hit
exact zero 8/8/24/12 times).

## 6c. SELECTOR COUPLING — P's one-line test, run on my parse (`selcouple.py`, `selcouple2.py`)
**Translation first.** My selector set = the 256 boolean free vars of the cluster cone. My
atom-local occurrence count per selector is **min 5, max 9, mean 5.89** (distribution
5:93, 6:112, 7:40, 8:8, 9:3). P reports "each of the 256 selectors appears in only 5–6 atoms."
That independent statistic matching is what licenses comparing the two — the objects are the same.

**Result, literal:** P reports *zero* atoms touching ≥2 distinct selectors. My parse finds **48**.
**Result, after classification (`selcouple2.py`):** of those 48 — **1** is a booleanity
certificate (a31143, a weighted sum of `2*x*(1-x)` and `x*x - x` terms over 6 selectors, which
constrains each selector separately and restricts no subset), and **47** are *bundled*: each
selector sits in its own additive term (`... + 19*(x_1058 - x_24365) + -31*(x_11540 - x_16586) + ...`).
**No atom anywhere multiplies two distinct selectors** (60 selector self-products, all booleanity).

So there is **no genuine cross-selector coupling in my parse either** — no cardinality atom, no
one-hot tie. The count differs from P's only because my decomposition is coarser and packs
independent per-selector load constraints into one additive atom; this is the same
decomposition-difference P noted for its 1,158 vs K's 12. **P's conclusion stands on my data.**
I am NOT reporting a conflict here.

## 6d. Reconciling "one-hot" (mine) with "free independent subset selection" (P's)
They describe **opposite ends of the same map**, and both are right.
- P's statement is about the **domain**: the 256 selectors are freely and independently
  choosable, 2^256 configurations, no atom couples them. §6c confirms this in my parse.
- My statement is about the **image**: the map from that free domain to the residual
  (a7389, a10187, a20212, a20215, a28647) is massively non-injective. `bfs.py` closed the image
  at **48 points**, and `sat.py` shows every subset — within-class and cross-class — reproduces
  the delta of a single member.

The two are compatible, and the compatibility is the load-bearing point: **a free independent
subset domain does not make the residual a subset-SUM.** If it were a sum, 2^256 inputs would
give ~2^256 residues and LLL on the bit vector would be the right tool. It gives 48. So my §4
refutation of "the residual is a subset-sum" is a claim about the *map*, not about the domain,
and it does not contradict P. Knob set and configuration for my version, as always: 256 cluster
booleans + both switches, base = triple8_seed, non-boolean affine knobs justified by §2.

## 6e. K's unreachability claim — consistency with my data
Consistent, and my data is independent positive evidence for it, with one honest limit.
§6b measured exactly K's predicted observable and found **0 degenerate cases in 48x333 probes**:
no row ever froze, so no stage was ever made degenerate by any configuration I reached. a20215's
never being exactly zero in 48 of 48 is the sharpest instance. But my evidence is *empirical
non-observation over a reached space*, not a proof over all 2^256 — my BFS converged on the
**image** (48 points) rather than enumerating the domain, so a degenerate configuration that is
isolated in the domain but never adjacent to cfg0 under single flips would be invisible to me.
So I corroborate K, I do not close it — and I specifically do NOT rely on K's carry-walk argument,
which P reports covers only one wrap (k = +-1) and is complete only if the governing modulus
exceeds the largest signed subset difference. My corroboration is independent of that step.

## 6f. §8.3 — is the joint p·Z² obstruction configuration-independent? (`kernel.py`, `kernel2.py`)
Knob set: the 54 affine knobs `lat2.system` measures at cfg0 (every single-row knob, plus each
atom's pure handle). Base configuration: cfg0 = triple8_seed, x_1530 = x_1603 = 1.

**Method.** `reach3.py`'s random selector scrambles were useless (66–467 bad atoms, blocking rows
nowhere near the cluster). Instead displace along the **affine kernel** — directions that hold
every non-target row at its satisfied value by construction — so every test configuration is
still a near-solution. Kernel dim at cfg0 is 7 over 54 knobs.

**Result A — the affine model is EXACT, not merely locally linear.** Every kernel displacement
tried (particular solution alone; 6 random small-coefficient combinations; 4 with coefficients up
to ±10⁶; unit vectors) landed on **bad atoms = exactly {a20215, a28647}, 28 fails, 39,005** —
identical to cfg0. Not one displacement broke the other rows. The measured structure was also
invariant: **54 knobs, 47 other rows, kernel dim 7** at every single one.

**Result B — the obstruction survived every kernel displacement**: membership of the residual in
the reachable lattice was **NO** every time, with the same p·Z² lattice and the same residual
class (a20215 ≡ 22981624690591…356252, a28647 ≡ 44159679639019…279353 mod p).

**⚠️ But state the limitation, because it is severe and it is structural.** Motion along the
kernel changes a20215 only by *multiples of p* — that IS the p·Z² result. So the residual's
**mod-p class is invariant along the kernel by construction**, and the membership answer *cannot*
change unless the measured knob set changes. Result B is therefore **not** evidence that the
obstruction is configuration-independent; it is a test of structural stability, and it passed.
Anyone citing Result B as "the obstruction is real" is making the §3 mistake again.

Exact counts (`runs_kernel.log`): **18 of 18 displacements blocked, 0 dissolved, 0 broke the
other rows.** NOTE: `kernel.py`'s own auto-generated closing line said this was "evidence it is a
statement about the instance". **That line was wrong** — it is the §3 error again. I have fixed
the script and appended a correction to the log so it cannot mislead a later reader.

**`kernel2.py` supplies the part that does move the class**: each BFS image point sits at a
different mod-p 5-tuple, so solve the other rows exactly there and re-test membership. This is
`lat3.analyse`, the correct formulation — `lat5.py` demanded *all* rows including the targets,
which is strictly stronger and hence a weaker probe of this question.
**Result: 14 image points, 5 distinct (a20215, a28647) mod-p classes covered —
blocked = 2, solved = 0, other-rows-infeasible = 12.**

**Honest verdict: the question is still open, and the experiment was starved.** Only **2 of 14**
image points were valid test cases; the other 12 could not be brought to a near-solution at all,
so they say nothing about the obstruction. Two blocked data points is not enough to call the
p·Z² obstruction configuration-independent, and Result B (18/18) does not help because kernel
motion cannot change the mod-p class by construction. **Do not upgrade any of this to a claim
about the instance.** What is established is only: (i) the affine model is exact under kernel
motion, and (ii) the obstruction is not dissolved at the 2 image points where it could be tested.

**The experiment worth running next**, and the reason this one starved: valid test cases need a
configuration that *both* moves the mod-p class *and* admits a near-solution, and BFS image points
mostly fail the second. Generate them the other way round — start from the cfg0 near-solution and
apply the **selector moves that change the class** (the trade knobs x_14853, x_6083, x_31339,
x_18956 are 1-for-1 and preserve solvability far better than selector flips), re-solving the other
rows after each move. That searches near-solutions by class rather than sampling classes and
hoping they are near-solutions.

## 6g. The trade-knob walk CANNOT answer the question — structural, checked up front (`trade.py`)
The proposed fix for §6f's starvation was: start from the cfg0 near-solution and move by the
1-for-1 trade knobs (x_14853, x_6083, x_31339, x_18956), re-solving the other rows each time.
**It is vacuous, and `trade.py` checks this before doing the walk rather than after:**

    54-knob affine set size = 54 ; trade knobs present in it: [14853, 6083, 31339, 18956] ; absent: []

**All four trade knobs are already inside the span `lat3.analyse` optimises over.** The analyse
step already explores every integer combination of them, so displacing along one and re-solving
cannot change the membership answer. Confirmed empirically, not just argued:
- every displaced point re-measures to the **identical** system — 54 knobs, 47 other rows,
  kernel dim 7 — and
- the **post-solve** residual is the identical class every time
  (a20215 ≡ 22981624690591…356252 mod p, one distinct value across the whole walk).

**⚠️ Counting trap, and I fell into it for one iteration.** The walk *does* report "VALID" cases
under the obvious criterion (class moved, other rows solvable) — because the class is measured at
the displaced point *before* the re-solve. The re-solve then washes the displacement out. **A case
is only an independent test if the POST-SOLVE residual class differs.** By that criterion the
whole trade walk is **one** test repeated, not N. Anyone reporting the raw VALID count from
`runs_trade.log` will overstate the evidence by roughly its length. Same shape as §3 and §6f
Result B: a quantity that looks like a measurement but is fixed by construction.

## 6h. Relaxed selectors (agent R's lead) — tested in my parse, `relax.py`
This *is* structurally outside the span: selectors are non-affine, so taking one off {0,1} changes
the measured system itself — the only thing that can move the membership answer. Confirmed: at
x_12714 ∈ {2, −1} the system re-measures to **53 knobs, 52 other rows, kernel dim 5** (vs
54/47/7) and the target class moves. So the move is genuine, unlike §6g.
**But the other rows go infeasible**, so it is not yet a valid test case either.
Counts when I stopped (runs were still going, see `runs_relax.log`): **5 attempts, 5
other-rows-infeasible, 0 valid cases.** Starvation rate so far 100%.

**Verdict on §8.3 overall: this line is closed for the trade knobs and starving for relaxed
selectors.** Valid-case counts across all three generators I tried:
| generator | attempts | valid (independent) cases | blocked | solved |
|---|---|---|---|---|
| BFS image points (`kernel2.py`) | 14 | 2 | 2 | 0 |
| kernel displacement (`kernel.py`) | 18 | **0** (class fixed by construction) | — | — |
| trade knobs (`trade.py`) | 14 | **1** (same test repeated; 1 distinct post-solve class) | 1 | 0 |
| relaxed selectors (`relax.py`) | 5 | **0** (all other-rows-infeasible) | 0 | 0 |
Three independent data points total, all blocked. **That is not enough to claim
configuration-independence and I am not claiming it.** The honest state: every cheap generator of
test configurations either cannot move the post-solve class (trade knobs, kernel) or destroys
solvability (relaxed selectors, random selector flips, BFS image points 12/14). The question is
not obviously answerable by sampling at all, and that itself is the finding.

**One divergence from R's lead, reported as a divergence and not a refutation** (different parses,
R's atom indices are not comparable to mine and I imported nothing): R reports that relaxing a
selector does not force its mux atoms nonzero, only the boolean-ness atoms. In my parse relaxing
x_12714 broke **6** atoms — `[10569, 20212, 20649, 20652, 32148, 32628]`, 74 fails — which
includes cluster/mux atoms a20212, a20649, a20652, a32148, not only a booleanity atom. Either the
parses decompose differently or the claim needs narrowing to particular selectors. Worth R
re-checking on the specific selectors it has in mind before the floor of 39,027 is relied on.

## 6i. ANSWER to "is there a direction that moves the post-solve class AND preserves solvability?"
**YES. Such directions exist, the infeasibility is NOT intrinsic to leaving the span, and it is
characterisable structurally** (`structural.py`, `structural2.py` — these read only my own run
logs, no new sampling, because my own starvation table says sampling cannot answer this).

**The mechanism is row deficiency.** Every `lat3.analyse` line records the shape of the
other-rows system: `knobs=K other-rows=M kernel-dim=d`, so rank = K − d and
**deficiency = M − (K − d)**. Recovered across **72 logged configurations**:

| deficiency | feasible | infeasible |
|---|---|---|
| 0 | **47** | 4 |
| > 0 | **0** | **21** |

- **deficiency > 0 ⟹ infeasible, 21 of 21, no exceptions.** This is why leaving the span usually
  kills solvability: breaking atoms adds *rows* faster than it adds *knobs*, the system becomes
  over-determined, and it dies. That is the structural answer to why relaxed selectors starved —
  x_12714 goes to (53, 52, 5), deficiency 4; x_16348 to deficiency 5.
- **Deficiency 0 is necessary but NOT sufficient** — 4 of 51 zero-deficiency systems were still
  infeasible. Full row rank over Q does not give solvability over Z; the residue conditions
  (`rhs % modulus != 0`) still bite.

**⚠️ CORRECTION to the table above — I inflated it, the same way the trade walk inflated its VALID
count.** The "47 feasible at deficiency 0" counts *log lines*, and **46 of them are cfg0's shape
(54/47/7) repeated** across the kernel and trade runs — the same configuration re-measured, not 47
configurations. Restricted to genuinely **distinct** configurations (the 22 image points swept):

| deficiency | feasible | infeasible |
|---|---|---|
| 0 | **2** (img0, img4) | **6** (img1, img6, img9, img11, img18, img19) |
| > 0 | **0** | **14** |

So the honest reading: **deficiency > 0 ⟹ infeasible survives (14/14 here, 21/21 overall)** and is
the real mechanism. But **deficiency 0 ⟹ feasible is only 2 of 8 (25%)**, not the ~92% the inflated
table implied. Deficiency 0 is a necessary condition and a weak predictor, not a generator.
That is the third time a count of repeated identical tests has masqueraded as independent
evidence in my own work (§3 image closure, §6g VALID count, here). **Check for repeats before
reporting any rate.**

**The existence proof is `img4`** (`runs_kernel2.log`), verified directly from the log rather than
from my scraper:

    [img4(|on|=1)] knobs=62 other-rows=54 kernel-dim=8 : FEASIBLE      -> rank 54 = rows, deficiency 0
    MEMBERSHIP of -residual in reachable lattice: NO
    row a20215 residual mod p = 84623865150894944922022514250283073331537809300837942433997946904552394898251
    row a28647 residual mod p = 47440525290535544708674620248717482764057423781520898897748888929433376949781

It **left cfg0's shape** (62/54/8 vs 54/47/7), **kept full row rank**, was **FEASIBLE**, and its
**post-solve class differs from cfg0's** on both coordinates. That is exactly the valid,
independent test case §6f was starved of — and it was **blocked**.

**This refines, and partly corrects, my §6g/§6h conclusion.** "The question may not be answerable
by sampling at all" is right about *blind* sampling and wrong as a general statement: it is
answerable with a **deficiency-directed generator**. The recipe follows from the table: look for
selector settings that add knobs at least as fast as they add rows (img4 gained +8 knobs against
+7 rows). Of 26 logged off-shape configurations, 5 had deficiency 0 and 1 of those was feasible —
a ~4% hit rate blind, but a targeted search optimising (K − d) − M should do far better than that.

**Status of the endgame condition: still open, but now with a working generator.** Independent
data points remain **2** (img0 at cfg0's class, img4 at a different class), both blocked. Two is
not configuration-independence and I am still not claiming it. But the residual side is no longer
blocked on "we cannot make test cases" — it is blocked on running a deficiency-directed search,
which is a concrete, bounded next step rather than an open-ended one.

## 6j. THE DEFICIENCY-DIRECTED SEARCH — run, and it STARVED (`dirsearch.py`, `combine.py`)
Pool: all 48 BFS image points (img4, the one existence proof, came from here). Sharded 3 ways.
Independence criterion applied throughout is **post-solve** class distinctness; the pre-solve
class is the trap and appears nowhere in `dirsearch.py` or `combine.py`.

    image points analysed          : 24 of 48
    other rows infeasible          : 22   (not test cases at all)
    other rows solvable            :  2   -> blocked 2, SOLVED 0
    INDEPENDENT test cases (distinct POST-SOLVE class): 2
       a20215=22981624690591... a28647=44159679639019...  img0   blocked
       a20215=84623865150894... a28647=47440525290535...  img4   blocked
    starvation rate: 92%

**The directed search starved, and I can now say why with a measured mechanism rather than a
suspicion.** Feasibility requires deficiency 0 (necessary, 21/21 no exceptions); deficiency 0
requires that breaking selectors adds knobs at least as fast as rows; and that in turn only
happens at **very low weight** — both feasible configurations found in the entire campaign are
|on| = 0 and |on| = 1. Every |on| ≥ 2 image point analysed so far is infeasible. The reachable
low-weight pool is tiny (the BFS enumerates ~7 configurations at |on| ≤ 1) and is already
exhausted, so there is no supply of independent test cases to be had from this pool at all.

**Verdict: 2 independent test cases, both blocked. That is NOT configuration-independence and I
am not claiming it.** It is the same 2 I had before the directed search; the search added 10 more
analysed configurations and 0 new independent cases. The line closes with a measured reason —
a 92% starvation rate and a structural account of it — rather than with a stretched claim.

**If anyone wants to reopen it**, the only untried supply of low-deficiency configurations is
*outside* the BFS-reachable pool: configurations at low weight that cfg0 cannot reach by flips
(§3's retraction showed the reachable closure is local, so such configurations exist). Generating
them needs a constructive method, not sampling — sampling is what starved here.

## 7. Scores
- Best verified: **39,026** — the existing deliverable, re-verified by me with `checker.py`.
- **I did not beat it.** E's 39,005 reproduced exactly. Best of my own constructions: 39,019.

## 8. Next experiments, priority order
1. **`lat5.py` — see §3b for the result.** Serial `lat5.py` is slow (~45 s/config); use
   `lat5p.py` instead, which takes `WK`/`NW` env vars and shards the configuration list across
   workers (`for w in 0 1 2; do WK=$w NW=3 python3 lat5p.py & done`). Note `lat5p.py` as written
   skips `i<18` — remove that guard to sweep from the start.
2. **Do NOT try to close the BFS.** `reach.py` settled it: the closure is local, half of random
   configurations land outside, and no amount of extra BFS budget makes it a global statement.
   The BFS was the wrong instrument — it explores the *domain* by adjacency when the thing that
   matters is the *joint lattice condition*, which `lat2.system` + `sparse.solve_sparse` tests
   directly at any configuration in ~40 s. Sample configurations, don't walk them.
3. **The single live question is whether the joint p·Z² obstruction is configuration-independent.**
   **`reach3.py` was my attempt and it FAILED to test it — do not cite it as corroboration.**
   It ran the full exact solve at 4 random configurations outside the BFS closure; all 4 came
   back infeasible, but that result is worthless for this question:
   - those configurations have **66 / 168 / 316 / 467 bad atoms** (cfg0 has 2), so they are
     wrecked states, nowhere near a solution;
   - the blocking rows were **4956, 1050, 364, 364** — nothing to do with the 5-row cluster.
   All it showed is that random selector settings break the instance everywhere, which was never
   in doubt. A real test needs configurations that are *outside cfg0's BFS closure but still
   near-solutions* (a handful of bad atoms). Constructing those is the open problem: start from a
   known-good state and move by the affine kernel directions (`lat3.py` gives kernel dim 7 at
   cfg0) rather than by scrambling selectors, since kernel moves preserve the other rows by
   construction. That is the experiment I would run next and did not get to.
4. The cancellation code is where the score actually lives. 7 is the code's optimum for the
   deliverable's atom family; the untried question is whether a **different** residual family
   (there are many, one per basin) admits a weight-6 code. `basin2.py` generalises to any basin —
   feed it other seeds.
5. ~~a726~~ **ANSWERED.** a726 (`x_24195 * x_19097`) is unhandled but **never bad**: across every
   configuration `lat5.py` has processed it is nonzero in 0 of them while appearing in the
   no-handle list. So it is not a forced-nonzero obstruction — it is a *satisfied but rigid* row:
   it constrains the solve (must stay 0) and could not be repaired if broken. Treat it as a
   side-condition on any repair, not as a target. a28647 is the only genuinely unhandled *bad* atom.
