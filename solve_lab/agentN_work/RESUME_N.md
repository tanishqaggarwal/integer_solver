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

## CORRECTION (agent T's audit, check-in 74) — read this before reusing my numbers
`atom_src[37887]` is `S*S` with textually identical operands, so **eq8680 = S⁴**, and `optN.inner`
strips one level and returns `S²` — a quadratic.  `sqaudit.py` measures the blast radius rather than
assuming it: a numeric affineness test (probe every knob at steps 1, 2, 3) finds **exactly two
non-affine rows anywhere in any model** — eq 8680 and eq 13985 — and **zero** after both are stripped
to their linear cores.
- **Unaffected**: the witness region (|R|=12) contains neither row, so the 924/924 p-obstruction,
  OPT = 5, the exhaustive |W| ≤ 2 collateral sweep and the 39,026 optimum all stand (re-derived on
  the corrected model: kernel dim 14, g = 5 at |W| = 0 and |W| = 1 exhaustive, 0 non-affine pairs).
- **Corrected**: at |R| = 13 row 8680 *is* in the region and was mis-modelled.  OPT is **6**, not 5;
  8680 **is** integrally zeroable; **all 16 detach states score 39,026**, not a 39,025/39,026 split.
- **RETRACTED**: my claim that "the knobs cannot reach S = 0, so detaching `x_28730` is the only
  way".  They can.  `agentN_work/N_r13_39026.json` is built from `D = []` (no detachment), verified
  independently by `solve_lab/checker.py` at **39026/39033** with the identical failing set.
O's Lemma is untouched — what changed is how many ways the frame has to satisfy it.

## The reduction, in agent T's checkable form (adopted)
A pool variable `v` is defined by an atom `(v − RHS)`, so **witness(v) ≠ gate(v) ⟺ that defining
atom is nonzero at the deliverable.**  The deliverable has 7 nonzero atoms; exactly 4 pool variables
have a nonzero defining atom — `{642, 28730, 29854, 31864}`.  T further verified that 0 of the other
61 reach a witness variable directly, transitively in the pool, or anywhere in the full 30,001-node
definition DAG — so the 2⁶⁵ lattice has **exactly 16 states by proof**, not by enumeration.

## Follow-up: O's Lemma cross-check, frame depth, re-orientation (all closed)
- **O's Lemma confirmed, and it is a SECOND obstruction, not mine.**  `eq_terms[8680] =
  (1, True, [(1, 37887)])`, and `optN.inner` returns the inner form, never its square, so my model
  already carried `T`.  `T = 0` **already holds at the witness**, and eq8680 is **exactly the one
  equation detaching `x_28730` buys** — O's Lemma *is* the 39,025 → 39,026 step.  The witness region
  excludes 8680, so the 924/924 p-obstruction is independent of it.  (`tcheck.py`)  [The extra
  claim I made here — that the knobs cannot reach `S = 0` — was wrong; see the CORRECTION above.]
- **Frame depth is not a lever.**  Deepening the pool from 65 to its saturation at 116 variables
  (frame free inputs 8,812 → 8,863) leaves the region's knob set at **49 wide / 7 narrow and OPT at
  5 at every depth**.  (`deepen.py`, `runs/deepen.json`)
- **Re-orientation is detachment, and every legal region move is worse.**  31.5% of atoms admit
  another `x_t - rest` reading.  For the region, **every legal unit target is already a free input
  of `Frame(POOL)` with measured response ±1**, so each re-orientation is a knob setting, not a new
  frame.  Executed: best single move 39,023, worst 39,007; **all 127 combinations ≤ 39,026, the best
  being the empty one**.  Atom 37887 (= T) has **no** legal unit target, so `T` can never be
  structurally forced to zero.  (`orient.py`, `reorient2.py`)

## Wholesale re-orientation (`fwd5.py`, `runs/fwd5.json`) — the 7 survive
Rebuilt the frame from scratch under 10 global target rules (first / last / lowvar / highvar / 5
randoms / one aimed at the region).  **No orientation beats the baseline forward score of 39,020**;
the best alternative (`prefer`) ties it with **the identical 13 failing equations**.  In every
orientation the failing equations reduce to nonzero atoms drawn from the same nine
`{22229, 22230, 22231, 35758, 35759, 35760, 35761, 35762, 37887}` — 3 to 5 left nonzero, never none.
**Atom 37887 is a check in all 10**: it has no legal unit target, so no orientation can ever force
`S = 0` structurally.  And fewer nonzero atoms is not better — `last` and `random/4` leave only 4 in
the whole instance yet score 39,006, because those atoms sit in more equations.

## Single next experiment
Everything I can vary inside this frame is now closed: the detach lattice exactly (16 states),
frame depth to saturation, collateral to budget 2, and re-orientation of the region.  The region's
7 nonzero atoms are `{22229, 22230, 35758, 35759, 35760, 35761, 35762}` and its knob set is 49,
permanently.  Global re-orientation is now done too and the carriers survive it.  Everything reachable by
choosing a frame, a detach set, knob values, or collateral up to budget 2 is closed at 39,026.
The one thing none of it touches: every model here is **linear in the knobs**, and the two rows that
were not (8680, 13985) turned out to be the only place a genuinely nonlinear response lives.  The
next experiment is to stop linearising: solve the region's 12 rows as a **polynomial** system in the
49 knobs — the atoms carry products like `x_26874*x_6947`, so the true response surface is quadratic
in places and the affine model can only ever find affine solutions.  Concretely: take the witness
region, keep the exact atom expressions instead of finite differences, and hand the resulting
integer polynomial system to a Groebner/`msolve` computation over the 7 zero-collateral knobs.  That
is the only unexplored direction left on this thread.
