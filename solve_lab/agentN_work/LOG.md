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
