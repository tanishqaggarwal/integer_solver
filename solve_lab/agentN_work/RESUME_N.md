# Agent N — RESUME.  Re-ranking by INTEGER REACHABILITY; triple-layer coverage.

## Best verified score
**39,026 / 39,033** — `solve_lab/best/new_instance_partial_39026.json`, re-verified by me with
`solve_lab/checker.py` (30 s load, `satisfied 39026/39033`, failing
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`).  Not beaten so far.

## What I corrected in agent H's pipeline
`agentH_work/stageB.py:solve_int` is **incomplete**.  After fraction-free elimination it sets
every non-pivot coordinate to 0 and requires each pivot division to be exact.  Systems that are
integrally solvable with nonzero free coordinates are reported unsolvable.  That is a FALSE
NEGATIVE in exactly the quantity H's closing note calls binding.

`agentN_work/zsolve.py` replaces it: `M_I t = -b_I` is solvable over Z iff `-b_I` lies in the
integer column lattice of `M_I`, tested by row-HNF (python-flint) + reduction.  Self-tested
against brute force on 500 random systems (both directions).  `max_zero_rows` then finds the
exact maximum number of simultaneously integrally-zeroable rows by max-clique-style DFS over
the downward-closed family of solvable row subsets (not capped at rank, unlike H's).

**Effect of the correction**: H's stage B reported the rank-8 detach sets zeroing exactly ONE
row.  The complete test gives **OPT = 5** for them.  H's scores for those sets were too low by
up to 5 equations.

## Established (measured)
1. **The region is rationally unobstructed; the entire barrier is integrality.**  For the
   witness region (|R| = 12, 7 zero-collateral knobs) all 12 rows are simultaneously solvable
   over Q, yet the integer optimum is 5.  A mod-q consistency pass rejects **no** row subset
   the search visits.
2. **The obstruction is a single prime.**  All C(12,6) = **924** six-row subsets of the witness
   region are integrally blocked, and in **924/924** the obstruction denominator is divisible by
   p = 115792089237316195423570985008687907853269984665640564039457584007908834671663.
   Knob set: the 7 zero-collateral free inputs of `Frame(POOL)` on the witness region;
   selector configuration: the witness's own, loaded from `best/new_instance_partial_39026.json`
   via `optN.BASEFV`.  So a 6th row needs a knob acting on the region with granularity finer
   than p — rational rank cannot see this, which is why H's criterion failed.
3. **Cascade pins, all 20, re-ranked exactly** (`runs/pins.jsonl`): best 39,018, unchanged.
   17 of 20 have OPT = **0** — not one row of their region is integrally zeroable.  The two
   best (a26729, a26731) reach OPT = 13 of |R| = 28.
4. **All 65 detach singletons**: OPT = 5 for every one; |R| = 13 for 64 of them and 12 for
   `{28730}` alone.  Score 39,025 for 64, 39,026 for `{28730}`.  The differentiator is |R|,
   not OPT.

## Re-entry
    cd solve_lab/agentN_work
    python3 optN.py            # calibration: witness -> |R|=12 OPT=5 score 39026
    python3 zsolve.py          # self-test of the integer solvability oracle
    python3 obstruct.py        # the p-obstruction measurement above
    python3 sweepN.py runs/tri_1.jsonl 3 1 3      # resumable; skips sets already in the file

## Files
- `zsolve.py`   complete integer solvability + exact max-integrally-zeroable-rows
- `optN.py`     corrected stage B over detach sets (~0.19 s/set uncontended)
- `carrierN.py` corrected stage B over the 20 cascade pins and the 1,147 handles
- `sweepN.py`   resumable sweep driver (JSONL, one record per set)
- `obstruct.py` localises the integrality obstruction to the prime p
- `runs/`       singles.jsonl, pairs.jsonl, tri_*.jsonl, pins.jsonl, handles_*.jsonl
