# Agent N — LOG.  Re-ranking by integer reachability; full coverage of the placement layers.

Baseline re-verified myself: `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
→ `satisfied 39026/39033 (7 failing) [12231,12270,12350,14584,18673,22044,29125]`.

---

## Step 1 — the stage-B integer test was incomplete, and it mattered

`agentH_work/stageB.py:solve_int` decides `M_I t = -b_I` over Z by fraction-free elimination,
then sets every non-pivot coordinate to 0 and requires each pivot division to be exact. Systems
solvable with nonzero free coordinates are reported unsolvable → **false negatives** in the very
quantity H's closing note identifies as binding.

`zsolve.py` replaces it with a complete test: `-b_I` lies in the integer column lattice of `M_I`
iff it lies in the integer row lattice of `M_I^T`; take the row HNF (python-flint) and reduce.
`max_zero_rows` then finds the exact maximum number of simultaneously integrally-zeroable rows by
max-clique-style DFS over the downward-closed family of solvable subsets — not capped at
`min(n, rank, 8)` as H's enumeration was.

Validation: `python3 zsolve.py` — 300 random systems against brute force in both directions, plus
200 random `max_zero_rows` instances against exhaustive subset scan. Calibration: the witness
placement reproduces `|R|=12, OPT=5, score 39026`.

**Size of the correction.** H's stage B reported the rank-8 detach sets zeroing exactly ONE row.
The complete test gives **OPT = 5** for them (score 39,025, not 39,020). Over the 1,147 handles,
**706 of 1,147** scores move upward once stage B is actually run.

---

## Step 2 — full coverage of the detach layers (H: singletons+pairs+4.0% of triples)

`optN.py` (exact, ~0.19 s/set) + `sweepF.py` (fast path, see step 3).

| layer | sets | coverage | OPT | best score |
|-------|------|----------|-----|-----------|
| singletons | 65 | **100%** | 5 for all 65 | 39,026 (`{28730}` only) |
| pairs | 2,080 | **100%** | 5 for all 2,080 | 39,026 (the 64 containing 28730) |
| triples | 43,680 | **100%** (H: 1,744 = 4.0%) | **5 for all 43,680** | 39,026 (the 2,016 = C(64,2) containing 28730) |
| quadruples | 677,040 | **100%** | 5 for all 677,040 | 39,026 (the 41,664 = C(64,3) containing 28730) |

Every run exhaustive (`exh=True`), no node-cap truncation. `|R| = 12` exactly when `28730 ∈ D`,
`|R| = 13` otherwise; `outside = 0` throughout. So `failing = |R| - 5`, and the entire detach axis
is 39,025 / 39,026 with nothing above.

## Step 3 — why the sweep collapses, and how full coverage became affordable

Measured: across placements the region `R` takes only **two** values —
`[2554,6816,8124,8680,9123,9421,12231,12270,12350,14584,18673,22044,29125]` and the same set minus
`8680` (the witness's) — and the zero-collateral knob list and knob-response matrix `M` depend only
on `R`. Only the row-target vector `b` varies, and over all 45,825 placements of size ≤ 3 it takes
**15 distinct values**. `sweepF.py` therefore prices by `(R, b)` and caches, with a random audit
that re-does the full build on cache hits and checks the knob list and `M` really match.
**Audits: 212 on triples, 35 on pairs, 20 on singletons — 0 mismatches.** Fast path cross-checked
set-by-set against the exact path on all 65 singletons: 0 mismatches.

The reason for the collapse is measured exactly in step 9: `make(D)` gives detached pool members
their **witness** values, and for 61 of the 65 pool variables the witness value already equals the
gate value, so detaching them changes nothing.  Only `{642, 28730, 29854, 31864}` differ, and
`x_28730` alone accounts for `|R| = 12`.

## Step 4 — the other two carrier classes, re-ranked exactly (H ran stage B only on top scorers)

- **20 cascade pins**, all of them (`runs/pins.jsonl`): best 39,018, unchanged. **17 of 20 have
  OPT = 0** — not one row of their region is integrally zeroable. `a26729`/`a26731` reach OPT = 13
  of |R| = 28.
- **1,147 handles**, all of them (`runs/handles_0.jsonl`): 1,146 at 39,017 and one at 39,011.
  `|R| - OPT = 16` **exactly, for all 1,146** — regardless of |R| ∈ {16,27,28,29,30,31}. All exhaustive.

Cost per carrier class (`|R| - OPT + outside`): witness placement **7**; other detach placements 8;
cascade pins 15–30; handles 16 (one 22).

## Step 5 — what actually blocks the sixth row

`obstruct.py`. For the witness region (|R| = 12, the 7 zero-collateral knobs of `Frame(POOL)`,
selector configuration = the witness's own, loaded from `best/new_instance_partial_39026.json`):

1. **All 12 rows are simultaneously solvable over Q.** The barrier is not rank, not consistency —
   it is entirely integrality. A mod-q consistency pass rejects **no** subset the search visits.
2. **All C(12,6) = 924 six-row subsets are integrally blocked, and in 924/924 the obstruction
   denominator is divisible by**
   `p = 115792089237316195423570985008687907853269984665640564039457584007908834671663`.

So a sixth row requires a knob acting on the region with granularity finer than `p`. Rational rank
cannot see that — which is exactly why H's `rank > deficit` criterion was over-generous.

## Step 6 — H's zero-collateral filter was NOT the binding restriction

`widen.py`. Every knob in H's pipeline had to disturb nothing outside the region individually.
That is only the coordinate-axis part of the collateral kernel; an integer *combination* of knobs
that individually break outside equations can have zero net collateral.

Wide candidate set = all 49 free inputs that move any atom of any region equation (they disturb
139 outside equations). `ker_Z(collateral map)` has **dimension 14** and rank **8** on the region —
strictly larger than H's rank-7 coordinate lattice. **OPT is still 5.** Same for `{28730}`,
`{17499}`, `{642,28730,31864}`.

## Step 7 — the region IS fully zeroable; the price is 69 equations

With the wide knob set and **no** collateral constraint, all 12 region rows are integrally zeroable
(OPT = 12 = |R|; 13 = |R| for the |R|=13 regions). `tradeoff.py` constructs it: `t0` has 18 nonzero
knobs with entries up to **4,917 bits**. Predicted collateral 69 outside rows; **actual score after
exact re-evaluation 38,964 — the linear model was exact at 4,917-bit knob values**. (For `{17499}`
the model predicted 73 and the truth was 74, so the model is good but not always exact — any hit
must be re-evaluated, never trusted from the model.)

A ±1 descent over the 37-dimensional `ker_Z(M_region)` did not reduce the 69.

## Step 8 — how much does collateral buy?  Exhaustively, for a budget of 0, 1 and 2: nothing

`drop.py`.  Let `W` be the outside equations we allow to break; admissible moves are then
`ker_Z(M_{outside\W})`, and `g(W)` = max region rows integrally zeroable inside that lattice:

        failing = (|R| - g(W)) + |W| ,   beating 7 needs   g(W) >= (|R|-6) + |W|.

| budget | subsets swept | coverage | max g | needed | verdict |
|--------|---------------|----------|-------|--------|---------|
| \|W\|=0 | 1 | exhaustive | 5 | 6 | no |
| \|W\|=1 | 139 | **exhaustive** | **5** | 7 | no |
| \|W\|=2 | 9,591 | **exhaustive** | **5** | 8 | no |

Each row of collateral bought would have to buy strictly more than one region row; in fact it buys
**zero**.  `g` does not move off 5 anywhere in the swept budget.  Over `|W| <= 2` the local model
cannot leave fewer than 7 failing equations at the witness placement.

Separately, `mus.py` extracted a minimal integrally-unsolvable subset of the 151-row local model:
**one MUS of size 31** (all 12 region rows plus 19 outside rows), after which the remaining 120
rows are solvable.  So the obstruction is a single large coupled block, not many small independent
ones — which is why the disjoint-MUS lower bound is only 1 and why the budget sweep above is the
informative measurement.

**Scope of this claim** (standing rule): knob set = all 49 free inputs of `Frame(POOL)` that move
any atom of any region equation — strictly larger than the 7 zero-collateral knobs H used;
selector configuration = the witness's own, from `best/new_instance_partial_39026.json` via
`optN.BASEFV`.  It is a statement about this frame and this configuration, at collateral budget
<= 2.  It is NOT a claim that nothing can move these rows in any frame.

---

## Confirmed / refuted

- **CONFIRMED** (H): integer reachability, not rational rank, is the binding quantity — and now
  with a reason: the region is rationally unobstructed and the obstruction is a single prime `p`.
- **CONFIRMED** (H): 39,026 survives. Best over 45,825 detach placements at 100% coverage for
  k ≤ 3 is exactly 39,026.
- **REFUTED**: H's stage-B integer test. It under-counts; 706 of 1,147 handle scores change.
- **REFUTED**: that the zero-collateral knob filter was what limited the search. Widening to the
  full net-zero-collateral integer kernel (rank 8 > 7) leaves OPT at 5.
- **REFUTED (as a reading of H's tables)**: "the rank>deficit sets zero exactly one row." They
  zero five.

## Step 9 — the detach axis is CLOSED, not sampled

`make(D)` gives detached pool members their **witness** values.  Comparing the witness value of
each pool variable against its gate value in the all-re-attached state: only **4 of the 65** differ,
and they are exactly `{642, 28730, 29854, 31864}` — the witness set.  Detaching any of the other
61 leaves every variable unchanged.  Therefore `make(D)` depends only on
`D ∩ {642,28730,29854,31864}`, the whole `2^65` detach lattice has exactly **16 states**, and that
matches the 16 distinct `(R,b)` signatures the 722,865-placement sweep measured.

All 16 priced exactly (`runs/detach_closure.json`): OPT = 5 for all 16, `outside = 0` for all 16,
`|R| = 12` iff `28730 ∈ D`.  Best over the entire lattice: **39,026**, attained by the 8 states
containing `28730`.  Verified: `make(POOL)` reproduces the witness assignment bit-for-bit
(score 39,026), `make([])` scores 39,020.

**The detach axis is finished.**  Sweeping k = 5, 6, ... is provably redundant.

## Quadruple layer (for the record, run before the closure was found)
677,040 sets, **100%**, 0.016 s/set, 1,337 audits with 0 mismatches: OPT = 5 for all 677,040,
`|R| = 12` for exactly 41,664 = C(64,3) (those containing 28730) → 39,026, the rest 39,025.
All exhaustive.

## Coverage table, final
| layer | sets priced | coverage | OPT | best |
|-------|-------------|----------|-----|------|
| detach singletons | 65 | 100% | 5 | 39,026 |
| detach pairs | 2,080 | 100% | 5 | 39,026 |
| detach triples | 43,680 | **100%** (H: 4.0%) | 5 | 39,026 |
| detach quadruples | 677,040 | **100%** | 5 | 39,026 |
| whole detach lattice | 2^65 | **100% exactly (16 states)** | 5 | **39,026** |
| cascade pins | 20 | 100% (H: top scorers only) | 0–13 | 39,018 |
| handles | 1,147 | 100% (H: top scorers only) | 0–17 | 39,017 |
| collateral budget \|W\| ≤ 2 | 9,731 | exhaustive | g = 5 | 39,026 |

---

# Follow-up (coordinator check-in 64): O's Lemma cross-check, and re-orientation

## Step 10 — O's Lemma is a SECOND, independent obstruction, not the one I measured

`tcheck.py`.  Confirmed in my frame: `eq_terms[8680] = (m=1, sq=True, [(1, 37887)])` — eq8680 is a
pure square of the single check atom 37887, whose source is
`(x_4432 - x_19964 - x_28730 + 6*(...) + ...)`, so syntactically `dT/dx_4432 = +1` and
`dT/dx_28730 = -1` as O states.  **`optN.inner` returns the INNER form, never its square**, so my
linear model already carried `T` rather than `T²` — no correction needed.

Then the cross-check, and the answer is **no, T = 0 is not what the 924/924 measures**:

- `T = 0` **already holds at the witness** (`T = 0` at `make(POOL)` and at `make([28730])`).
- eq8680 is **exactly the one equation detaching `x_28730` buys**: `fixed by detaching 28730: [8680]`,
  broken: none.  So O's Lemma *is* the 39,025 → 39,026 step, precisely.
- The witness region is the 12 rows **excluding 8680**.  `T = 0` is therefore not among the
  constraints the 924 six-row subsets are asked to satisfy.  The p-obstruction is **a second,
  independent obstruction** on the remaining 12 rows.

And a new result that ties the two together.  In the 13-row region (where `T ≠ 0`):
row 8680 is **not individually integrally zeroable** with the zero-collateral knobs, and the
**max rows zeroable subject to 8680 being zeroed is 0**.  So the knobs cannot reach `T = 0` at all —
the *only* way to satisfy O's mandatory constraint is to detach `x_28730`, which is exactly what the
witness does.  O's Lemma says `T = 0` is compulsory; my measurement says the frame can obtain it in
exactly one way.

## Step 11 — the frame-depth axis is saturated: same 49 knobs at every depth

`deepen.py`.  `pool.py` stops two levels above the region atoms.  Detaching further up adds free
inputs without changing the witness assignment, so it is a strict enlargement of the knob set at a
fixed state.

| depth | pool | frame free inputs | \|R\| | wide knobs | narrow knobs | OPT wide | OPT narrow | failing |
|---|---|---|---|---|---|---|---|---|
| 2 | 65 | 8,812 | 12 | 49 | 7 | 12 | **5** | 7 |
| 3 | 81 | 8,828 | 12 | 49 | 7 | 12 | **5** | 7 |
| 4 | 95 | 8,842 | 12 | 49 | 7 | 12 | **5** | 7 |
| 5 | 111 | 8,858 | 12 | 49 | 7 | 12 | **5** | 7 |
| 6 | 114 | 8,861 | 12 | 49 | 7 | 12 | **5** | 7 |
| 8, 12 | 116 (saturated) | 8,863 | 12 | 49 | 7 | 12 | **5** | 7 |

The pool saturates at 116 variables and the frame at 8,863 free inputs, and **the knob set on the
region is 49 wide / 7 narrow at every depth**.  Deeper detachment adds free inputs, none of which
touch the region.  Score 39,026 throughout.

## Step 12 — re-orientation, executed: every legal move is realizable in the current frame, and worse

`orient.py` + `reorient2.py`.  Census: 13,332 of 42,267 atoms (31.5%) admit more than one legal
`x_t - rest` reading; 10,956 of the 30,001 definition atoms do.

For the region the picture is exact.  Every legal unit target of every region atom is **already a
free input of `Frame(POOL)` with measured response exactly ±1**:

| atom | value at witness | legal unit targets (all already free in the frame) |
|---|---|---|
| 22229 | nonzero | `x_7068` (+1), `x_2099` (−1) |
| 22230 | nonzero | `x_28730` (+1) |
| 22231 | 0 | `x_4432` (+1), `x_19964` (−1), `x_28730` (−1) |
| 35758 | nonzero | `x_29854` (+1) |
| 35759 | nonzero | `x_29854` (−1) |
| 35760 | nonzero | `x_31864` (+1) |
| 35761 | nonzero | `x_31864` (+1) |
| 35762 | nonzero | `x_642` (+1) |
| **37887 (= T)** | 0 | **NONE** |

**Why that settles it.** Orienting a check atom `x_v − rest` into a definition forces the atom to 0
for every free-input choice; but where `x_v` is already free, "force the atom to 0" and "choose the
value of `x_v` that zeroes it" describe the **same set of assignments**.  More generally, re-orienting
an atom from target `x_u` to target `x_w` makes `x_u` free and turns `x_w`'s old definer into a check —
which is exactly what detaching `x_u` does in `frameB`.  **Re-orientation is detachment**, so the
detach closure and the depth saturation already cover it.

Executed anyway, as knob settings, and measured with the real scorer:

- zero atom 22229 via `x_7068` → **39,008**; via `x_2099` → **39,007**
- zero 22230 via `x_28730` → **39,021**;  35758 via `x_29854` → **39,023**;  35759 → **39,023**
- zero 35760 via `x_31864` → **39,022**;  35761 → **39,022**;  35762 via `x_642` → **39,021**
- **all 127 combinations of the 7 available moves: best is the empty combination, 39,026.**

Note the last row of the table: **atom 37887 = T has no legal unit target at all**, so `T` can never
be oriented into a definition.  `T = 0` must be obtained by value, and (step 10) the only value move
that obtains it is detaching `x_28730`.

**Residual, stated honestly.** This closes re-orientation *for the region*.  A global re-orientation
of the 10,956 re-orientable definition atoms elsewhere in the circuit would change which equations
are auto-satisfied outside the region; I did not rebuild `fwd2.pkl` wholesale, because the argument
above shows each such swap is a detachment, and the score is decided in the region.

## Confirmed / refuted (follow-up)
- **CONFIRMED** (O): eq8680 = T², T linear with `dT/dx_4432 = +1`, `dT/dx_28730 = −1`; and my model
  already carried T, not T².
- **REFUTED** (as an explanation of my result): that `T = 0` is what the 924/924 obstruction
  measures.  `T = 0` already holds at the witness and 8680 is not in the witness region — the two
  are **independent obstructions**, and O's is exactly the 39,025 → 39,026 step.
- **REFUTED**: that re-orientation is a new axis.  It is detachment, every legal region move is
  realizable in the existing frame, and all 127 combinations are ≤ 39,026.
- **REFUTED**: that frame depth is a lever.  The region's knob set is 49/7 at every depth to
  saturation.

---

# Follow-up 2 (check-in 74): T's S⁴ flag — a real correction — and the wholesale re-orientation

## Step 13 — CORRECTION and RETRACTION.  T's note was right and it was load-bearing

`atom_src[37887]` parses as `BinOp(Mult)` with the two operands **textually identical**, so
`av[37887] = S²` and `eq_terms[8680] = (1, sq=True, [(1,37887)])` gives **eq8680 = S⁴**.
`optN.inner` strips one level and returns `S²` — a quadratic. T's diagnosis is exactly right.

**Blast radius, measured not assumed.** `sqaudit.py` scans every row of every model for top-level
square atoms AND tests affineness numerically (probe each knob at steps 1, 2, 3; require
`value(t) − value(0) = t·(value(1) − value(0))`). Result: **exactly two non-affine rows anywhere** —
eq 8680 (atom 37887) and eq 13985 (atom 39967). Every other row of every model is exactly affine
against all 49/53 knobs. After stripping both rows to their linear cores, **0 non-affine pairs remain**.

**What survives.** The witness region (|R| = 12) contains *neither* row — 8680 and 13985 are both
outside it and both have core zero there. So the 924/924 p-obstruction, OPT = 5 at |R| = 12, and
the exhaustive detach closure at |R| = 12 are **unaffected**.

**What changes — and I am retracting a claim I sent you.** In the |R| = 13 regions row 8680 *is* in
the region, and it was mis-modelled. Corrected:

| | mis-modelled (S²) | corrected (S) |
|---|---|---|
| OPT at \|R\|=13 | 5 | **6** |
| 8680 individually zeroable | **no** | **yes** |
| max rows zeroable subject to 8680 = 0 | **0** | **6** |
| score | 39,025 | **39,026** |

**RETRACTED**: my step-10 statement that "row 8680 is not individually integrally zeroable and max
rows zeroable subject to 8680 being zeroed is 0, so detaching `x_28730` is the only way to reach
`T = 0`." That was computed on the quadratic row. **The knobs CAN reach `S = 0`.** In every |R| = 13
state the optimal 6 rows are `[2554, 6816, 8124, 8680, 9123, 9421]` — 8680 among them.

Constructed and **independently verified with `solve_lab/checker.py`**:
`agentN_work/N_r13_39026.json`, built from `D = []` (no detachment at all), largest variable 909
digits — `satisfied 39026/39033`, failing `[12231,12270,12350,14584,18673,22044,29125]`.
So there are (at least) two routes to `S = 0`: detach `x_28730`, or set the knobs. Both land on
39,026 with the identical failing set.

**Corrected closure table**: all **16** detach states score **39,026**, not the 39,025 / 39,026 split
I reported. The `|R| = 12` vs `|R| = 13` distinction dissolves once the row is modelled correctly.
The conclusion (39,026, not beaten) is unchanged; the intermediate numbers were wrong.

O's Lemma still holds and is still the reason 8680 must be zeroed — what changes is *how many ways*
the frame has to do it. Please relay the retraction to O and T.

## Step 14 — T's identity, adopted

T's formulation is better than mine and I am adopting it: a pool variable `v` is defined by an atom
`(v − RHS)`, so **witness(v) ≠ gate(v) ⟺ that defining atom is nonzero at the deliverable.** The
deliverable has exactly 7 nonzero atoms, and exactly 4 pool variables have a nonzero defining atom —
`{642, 28730, 29854, 31864}`. That makes the 4-of-65 reduction checkable without my frame. T also
closed the gap my argument left (0 of the 61 reach a witness variable directly, transitively in the
pool, or anywhere in the full 30,001-definition DAG), so the 2⁶⁵ lattice has **exactly 16 states by
proof**, not by enumeration.

## Step 15 — the wholesale re-orientation: the 7 survive

`fwd5.py` rebuilds the frame from scratch under a chosen target rule (fix a target per atom, then
run `fwd2`'s propagation), then forward-evaluates from the witness values on that frame's free inputs.
Definition semantics need no re-parsing: an atom is `s·x_v + F` with `s = ±1`, so `x_v ← x_v − s·atom`.

| rule | defs | checks | free | score | failing | nonzero checks (whole instance) | region atoms left nonzero |
|---|---|---|---|---|---|---|---|
| **fwd2 baseline** | 30,001 | 12,266 | 8,747 | **39,020** | 13 | — | — |
| first | 30,970 | 11,297 | 7,778 | 38,996 | 37 | 6 | 35759, 35760 |
| last | 23,170 | 19,097 | 15,578 | 39,006 | 27 | 4 | 22230, 35759, 35760 |
| lowvar | 25,863 | 16,404 | 12,885 | 39,005 | 28 | 5 | 35759, 35761, 22231, 37887 |
| highvar | 25,878 | 16,389 | 12,870 | 38,999 | 34 | 7 | 22230, 35758, 35760 |
| random/1 | 25,384 | 16,883 | 13,364 | 38,984 | 49 | 12 | 22230, 35758, 35761 |
| random/2 | 25,103 | 17,164 | 13,645 | 39,005 | 28 | 5 | 35758, 35760, 22231, 37887 |
| random/3 | 25,136 | 17,131 | 13,612 | 38,955 | 78 | 16 | 35759, 35760, 37887 |
| random/4 | 25,265 | 17,002 | 13,483 | 39,006 | 27 | 4 | 22230, 35759, 35761 |
| random/5 | 25,351 | 16,916 | 13,397 | 39,005 | 28 | 5 | 35758, 35760, 22231, 37887 |
| prefer (aimed at the region) | 30,965 | 11,302 | 7,783 | **39,020** | **13** | 5 | 22229, 35759, 35760, 22231, 37887 |

**No orientation beats the baseline.** The best alternative (`prefer`, built specifically to turn
region atoms into definitions) ties at 39,020 with **the identical 13 failing equations**.

**The 7 survive.** In every orientation the failing equations reduce to nonzero atoms drawn from the
same nine `{22229, 22230, 22231, 35758, 35759, 35760, 35761, 35762, 37887}` — 3 to 5 of them are left
nonzero, never none. Which ones varies; the carrier set does not.

Two structural facts fall out:
- **Atom 37887 is a CHECK in all 10 orientations.** It has no legal unit target, so no orientation can
  ever force it to zero. `S = 0` is always a value condition, never a structural one — the global
  version of what I found locally.
- **Fewer nonzero atoms is not better.** `last` and `random/4` leave only **4** nonzero check atoms in
  the whole instance — fewer than the baseline — yet score 39,006 vs 39,020, because those atoms sit
  in more equations. Atom count is the wrong objective; equation incidence is the right one.

**Scope.** This is the forward score from witness free-input values in each frame (the all-attached
analogue, which is 39,020 in the baseline). It is not each frame's optimum after knob optimisation.
The reason I did not compute those: re-orientation is detachment (step 12), the region's knob set is
49 at every frame depth to saturation (step 11), and the detach lattice is closed at 16 states by
T-verified proof (step 14) — so the region optimum is 39,026 regardless of orientation.

## Cross-link
`x_28730` is simultaneously one of my 4 witness variables, one of the h-wires in L's cancellation
set, and the variable entering O's `S` with `dS/dx_28730 = −1`. Three threads, one wire.

## Confirmed / refuted (follow-up 2)
- **CONFIRMED** (T): eq8680 is S⁴, not S²; `inner` strips one level too few. Real defect, 2 rows.
- **RETRACTED** (mine): "the knobs cannot reach `S = 0`; detaching `x_28730` is the only way."
  Corrected: 8680 is zeroable, OPT at \|R\|=13 is 6 not 5, and all 16 detach states score 39,026.
- **CONFIRMED** (T): the 4-of-65 reduction, now by proof via the nonzero-defining-atom identity.
- **REFUTED**: that a different global orientation changes the carriers. Over 10 orientations the
  failing equations always reduce to the same nine region atoms, and 37887 is never a definition.

## Step 16 — the corrected model re-run end to end, and it reproduces

`fixrows.py` rebuilds the witness model with BOTH square rows stripped to their linear cores and
re-runs everything that depended on them:

```
|R|=12 outside=139 knobs=49 ; non-affine (row,knob) pairs AFTER the fix: 0
outside rows nonzero at the witness after stripping: 0
region OPT with ALL 49 knobs, no collateral limit: 12 of 12 (exh=True)
|W|=0: kernel dim 14, g=5, need 6 -> no
|W|=1:   139 subsets exhaustive, max g=5, need 7 -> no
|W|=2: 9,591 subsets exhaustive, max g=5, need 8 -> no
=> corrected model, |W| <= 2: best failing = 7 (score 39026)
```

Identical to the pre-correction numbers, and now on a model with **zero** non-affine rows. The
witness-placement results were never at risk (neither square row is in that region); this makes that
explicit rather than assumed.
