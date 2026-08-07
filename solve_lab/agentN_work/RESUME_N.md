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

## Follow-up: O's Lemma cross-check, frame depth, re-orientation (all closed)
- **O's Lemma confirmed, and it is a SECOND obstruction, not mine.**  `eq_terms[8680] =
  (1, True, [(1, 37887)])`, and `optN.inner` returns the inner form, never its square, so my model
  already carried `T`.  `T = 0` **already holds at the witness**, and eq8680 is **exactly the one
  equation detaching `x_28730` buys** — O's Lemma *is* the 39,025 → 39,026 step.  The witness region
  excludes 8680, so the 924/924 p-obstruction is independent of it.  New: in the 13-row region the
  max rows zeroable **subject to 8680 being zeroed is 0**, so the knobs cannot reach `T = 0` at all;
  detaching `x_28730` is the only way.  (`tcheck.py`)
- **Frame depth is not a lever.**  Deepening the pool from 65 to its saturation at 116 variables
  (frame free inputs 8,812 → 8,863) leaves the region's knob set at **49 wide / 7 narrow and OPT at
  5 at every depth**.  (`deepen.py`, `runs/deepen.json`)
- **Re-orientation is detachment, and every legal region move is worse.**  31.5% of atoms admit
  another `x_t - rest` reading.  For the region, **every legal unit target is already a free input
  of `Frame(POOL)` with measured response ±1**, so each re-orientation is a knob setting, not a new
  frame.  Executed: best single move 39,023, worst 39,007; **all 127 combinations ≤ 39,026, the best
  being the empty one**.  Atom 37887 (= T) has **no** legal unit target, so `T` can never be
  structurally forced to zero.  (`orient.py`, `reorient2.py`)

## Single next experiment
Everything I can vary inside this frame is now closed: the detach lattice exactly (16 states),
frame depth to saturation, collateral to budget 2, and re-orientation of the region.  The region's
7 nonzero atoms are `{22229, 22230, 35758, 35759, 35760, 35761, 35762}` and its knob set is 49,
permanently.  The only untested thing left is a **global** re-orientation: the 10,956 definition
atoms elsewhere in the circuit that admit another reading decide which equations are auto-satisfied
OUTSIDE the region.  Rebuild `fwd2.pkl` wholesale under a different target rule, then check whether
the 7 failing equations still reduce to the same 7 nonzero atoms.  If they do, 39,026 is the frame's
ceiling under every orientation and the search must move to a different `checks` decomposition
entirely.
