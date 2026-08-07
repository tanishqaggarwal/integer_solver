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

## 3. Half the condition is reachable; the other half is not (scoped)
`bfs.py` — BFS over configurations, move set = **all 256 cluster booleans + the switch x_30163**,
non-boolean affine knobs at cfg0 values (justified: §2 shows they act only in p·Z² on the pair).
It **terminated by exhaustion** (gen6 produced 0 new), giving exactly **48 distinct mod-p
5-tuples**:

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
`final.py` re-measured at two such configurations: a20215's handle step is still exactly p
there, and the other rows go infeasible. **a20215 mod p = 0 is the one thing nothing reached.**
Scope, stated plainly: move set = 256 cluster booleans + both switches + 75 affine cone knobs;
base = triple8_seed. This is the whole cone of the 5 rows, so no free variable outside it can
affect them — but bfs2 did not converge, so this is *measured*, not proved.

### 3b. `lat5.py` corroboration (independent of the BFS mod-p reading)
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

## 7. Scores
- Best verified: **39,026** — the existing deliverable, re-verified by me with `checker.py`.
- **I did not beat it.** E's 39,005 reproduced exactly. Best of my own constructions: 39,019.

## 8. Next experiments, priority order
1. **`lat5.py` was still running at 26/48 configurations, 0 feasible.** Finish it
   (`python3 lat5.py`, ~40 s/config). It re-measures knobs, handles AND targets per
   configuration and does the full exact solve — the cleanest remaining falsification test.
2. **Let `bfs2.py` converge** (it had not; raise the time limit and the `frontier=nf[:40]` cap).
   If it closes with a20215 mod p still taking 2 values, §3 becomes a closed statement over the
   entire cone rather than a measurement.
3. **Attack a20215 directly rather than the cluster.** a20215 is `x_24530 - x_5647 * x_24908`
   with 274 cone free vars. Everything above says the whole endgame is the single condition
   `a20215 ≡ 0 (mod p)`; that is now a much smaller target than "solve the instance".
4. The cancellation code is where the score actually lives. 7 is the code's optimum for the
   deliverable's atom family; the untried question is whether a **different** residual family
   (there are many, one per basin) admits a weight-6 code. `basin2.py` generalises to any basin —
   feed it other seeds.
5. ~~a726~~ **ANSWERED.** a726 (`x_24195 * x_19097`) is unhandled but **never bad**: across every
   configuration `lat5.py` has processed it is nonzero in 0 of them while appearing in the
   no-handle list. So it is not a forced-nonzero obstruction — it is a *satisfied but rigid* row:
   it constrains the solve (must stay 0) and could not be repaired if broken. Treat it as a
   side-condition on any repair, not as a target. a28647 is the only genuinely unhandled *bad* atom.
