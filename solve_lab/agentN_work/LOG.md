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
| quadruples | 677,040 | in progress | — | 39,026 so far |

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

The reason for the collapse: `make(D)` gives the detached pool members their **witness** values.
For 64 of the 65 pool variables the witness value already equals the gate value, so detaching them
changes nothing. `x_28730` is the exception, and it alone accounts for `|R| = 12`.

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
