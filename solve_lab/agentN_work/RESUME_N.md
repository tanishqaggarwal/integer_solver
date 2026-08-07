# Agent N — RESUME.  Re-ranking by INTEGER REACHABILITY; the detach axis closed exactly.

## Best verified score
**39,026 / 39,033** — `solve_lab/best/new_instance_partial_39026.json`, re-verified by me twice with
`solve_lab/checker.py`: `satisfied 39026/39033 (7 failing)`, failing
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`.  **I did not beat it.**
Nothing I produced needs `verifyE.py`: my states stay inside `checker.py`'s 4,300-digit limit.

## The correction (H's stage B was wrong, and it mattered)
`agentH_work/stageB.py:solve_int` sets every non-pivot coordinate to 0 and demands exact pivot
divisions, so integrally solvable systems get reported unsolvable — **false negatives** in the one
quantity H's closing note calls binding.  `zsolve.py` replaces it with a complete test (`-b_I` in
the integer column lattice of `M_I`, decided by row-HNF over Z via python-flint) plus an exact
max-clique-style search for the largest integrally zeroable row subset, not capped at
`min(n, rank, 8)`.  Self-tested against brute force on 500 random systems in both directions.
- H reported the rank-8 detach sets zeroing **1** row; the truth is **5**.
- **706 of 1,147** handle scores move once stage B is actually run.

## Coverage delivered (H: singletons, pairs, 4.0% of triples)
| layer | sets | coverage |
|-------|------|----------|
| singletons | 65 | 100% |
| pairs | 2,080 | 100% |
| triples | 43,680 | **100%** |
| quadruples | 677,040 | **100%** |
| **whole detach lattice** | **2^65** | **100%, exactly — see below** |
| cascade pins | 20 | 100% |
| handles | 1,147 | 100% |

**OPT = 5 for every one of the 722,865 detach placements**, every run exhaustive, `outside = 0`
throughout.  `|R| = 12` iff `28730 ∈ D`, else 13; so `failing = |R| - 5` and the axis is
39,025 / 39,026 with nothing above.

## Why the axis closes exactly (this retires it, not just samples it)
`make(D)` gives detached pool members their **witness** values.  Only **4 of the 65** pool
variables have a witness value different from their gate value, and they are exactly
`{642, 28730, 29854, 31864}` — the witness set.  Detaching any of the other 61 is a literal no-op.
So `make(D)` depends only on `D ∩ {642,28730,29854,31864}`: the whole 2^65 detach lattice has
**16 distinct states**, matching the 16 distinct `(R,b)` signatures the sweeps measured.
All 16 are priced exactly in `runs/detach_closure.json`: OPT = 5 for all 16, best **39,026**.
**The detach axis is finished; do not sweep k = 5, 6, ... — it is provably the same 16 states.**

## What blocks the sixth row (`obstruct.py`)
For the witness region (|R| = 12; knob set = the 7 zero-collateral free inputs of `Frame(POOL)`;
selector configuration = the witness's own, from `best/new_instance_partial_39026.json`):
1. all 12 rows are simultaneously solvable **over Q** — the barrier is purely integrality;
2. all **924/924** six-row subsets are integrally blocked, and in **924/924** the obstruction
   denominator is divisible by
   `p = 115792089237316195423570985008687907853269984665640564039457584007908834671663`.

## Two of H's restrictions lifted, both negative
- **Zero-collateral filter lifted** (`widen.py`): H only ever used knobs that individually disturb
  nothing outside the region.  Using all 49 region-touching free inputs and taking the integer
  kernel of the collateral map gives a **rank-8** admissible lattice (vs H's rank 7).  OPT still 5.
- **Collateral allowed** (`drop.py`): with `W` outside equations permitted to break,
  `failing = (|R| - g(W)) + |W|`.  Exhaustive at `|W| = 0` (1 subset), `|W| = 1` (139) and
  `|W| = 2` (9,591): **max g = 5 at every budget**.  Collateral buys zero region rows.
- With the wide knobs and no collateral limit the region *is* fully zeroable (`tradeoff.py`,
  `t0` up to 4,917 bits) but costs 69 equations → 38,964.  The affine model was exact there.

## Re-entry
    cd solve_lab/agentN_work
    python3 zsolve.py                       # self-test of the integer oracle
    python3 optN.py                         # calibration: witness -> |R|=12 OPT=5 score 39026
    python3 obstruct.py                     # the p-obstruction measurement
    python3 drop.py 2                       # collateral-budget sweep
    python3 sweepF.py runs/x.jsonl 3        # resumable layer sweep (skips sets already present)

## Single next experiment
Every axis measured here is inside one frame and one selector configuration.  The detach axis is
now closed exactly, so the next lever is the **frame** itself: `fwd2`'s orientation makes 61 of the
65 pool variables gate-consistent with the witness, which is why they are inert.  Re-orient — pick
a different acyclic orientation of the atom set so that a *different* group of variables carries
witness-vs-gate disagreement — and re-run `optN.price` on the resulting 16-state closure.
Concretely: rebuild `fwd2.pkl` choosing, for each atom with several `x_t - rest` readings, a
different target `x_t`, then recompute the pool and its closure.  That is the only way to change
`b` beyond the 16 values reachable now, and `b` is the only input to OPT that ever varied.
