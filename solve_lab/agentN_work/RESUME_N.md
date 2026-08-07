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

---

# POST-RESTART UPDATE (step 17): the polynomial direction is CLOSED, and the barrier now has a certificate

## Environment after the restart
`fwd2.pkl` and all site-packages were gone. Rebuilt `fwd2.pkl` from `fwd2.py` (30,001 defs / 12,266
checks / 8,747 free — identical). Installed: **sympy 1.14.0, python-flint 0.9.0, Singular 4.3.2**
(`apt install singular`). **NOT available: msolve, Sage, Macaulay2, Magma, CoCoA, PARI/gp.**
Deliverable re-verified: **39,026/39,033**, failing `[12231,12270,12350,14584,18673,22044,29125]`.

## SIZE FIRST (as ordered), then the solve

| knob set | unknowns | rows | max total degree | max terms/row | max coef bits |
|---|---|---|---|---|---|
| 7 narrow | 7 | 12 | **1** | 8 | 2,435 |
| 49 wide (widen.py) | 49 | 12+139 | **2** | 48 | 6,083 |
| **68 complete** | **68** | **12+231** | **4** | **665** | **6,083** |

The 7 narrow knobs touch **zero** downstream DAG variables, so that model is affine **by proof**.
Everything computed on it (the 924/924 p-obstruction, OPT = 5) is exact.

## The result: the polynomial system is a linear system in disguise
`poly68b.py` iterates to a fixed point, using a collateral row as a linear constraint only when its
**exact** restriction is linear: rank 68 -> 37 -> 15. **Exactly one nonlinear generator survives in
the whole system — eq 8680 — and on the rank-15 lattice it is a single monomial `s^2` with
coefficient 1.** Singular: `dim 14, radical_dim 14, radical_gbsize 1, radical_gens_are_linear = 1,
one component of dim 14 and generator degree 1`. **The ideal's radical is generated by linear forms**,
and over Z a perfect square is the same condition as its base. All 12 region rows are degree 1 there.

`solve68.py`: **OPT = 5 of 12, exhaustive, 0 of 924 six-subsets integrally solvable, score 39,026.**
Carrying the products changes nothing. **There is no assignment to dump — none beats the deliverable.**

## A mod-p CERTIFICATE replaces the search
With `p = 115792089237316195423570985008687907853269984665640564039457584007908834671663`:
on the exact rank-14 variety **791 of 792 five-subsets, 924 of 924 six-subsets, 792 of 792
seven-subsets** are **inconsistent mod p**; the single consistent 5-subset is exactly the OPT set
`[2554,6816,8124,9123,9421]`. Cross-checked against zsolve: **the certificate settles 924/924 —
p is the whole story.** `p = 2,3,5,1000003,2^61-1` each certify **0 of 924**. `b_i mod p = 0` in
exactly the 5 zeroable rows and `!= 0` in exactly the 7 failing ones.

## Corrections, with blast radius measured
1. **widen.py's knob set was incomplete.** Its step-1 filter **missed 18 knobs**
   `[2239,4068,4339,8173,9106,9325,11368,13502,14466,18822,24490,24559,26064,26874,28969,30095,31731,34660]`.
   Exact polynomial support: 67 of 68 candidates genuinely move the region.
2. **widen.py's kernel was the right rank for the wrong lattice.** Its secant matrix differs from
   the true linear part in exactly **28** of 7,399 (row,knob) entries (the 28 squares). Both
   lattices have rank 14 but **their union has rank 15**, so neither contains the other.
3. **Blast radius: ZERO on every conclusion.** OPT = 5 on the exact model, on the complete knob set,
   at |W| = 0 and |W| = 1. 39,026 stands.

## A barrier, stated with its knob set and configuration
On the **zero-collateral rank-14 lattice** (Frame(POOL), D = `[642,28730,29854,31864]`, selectors
from `best/new_instance_partial_39026.json`, knob set = all 68 free inputs supporting the region),
eq **29125**'s entire knob response is `= 0 mod p` while `b != 0 mod p` — unzeroable there.
**This is NOT ambient**: over the full 68-knob space, `pambient.py` finds **0 of 12 rows unzeroable
mod p**. The barrier belongs to the lattice, not to the row.

## Collateral budgets (need `g >= W + 6` to beat 39,026)
- **|W| = 0**: g = 5, exhaustive, mod-p certified. REFUTED.
- **|W| = 1**: all **231** live collateral rows dropped in turn on the complete knob set; 217 exact,
  14 by relaxation (relaxing can only raise the optimum). **max g = 5.** REFUTED.
- **|W| = 2**: 26,565 pairs, **priced at ~63 CPU-hours / ~11 wall-hours** under current contention.
  **1,554 done, max g = 5.** Resumable: `python3 budget68b.py w2 <shard> 4`, appends to
  `runs/budget68_w2_<shard>.jsonl`, skips completed pairs.

**Why the rest is predictable** (`pgrow.py`): paying collateral raises `rk_Q(M)` 7 -> 9, but
`rk_p(M)` and `rk_p([M|b])` rise **in lockstep** (3->4 and 4->5), so the inconsistency gap stays
**exactly 1** in all 15 lattice-enlarging |W| = 1 cases. And `rk_Q(M) = 7` on a rank-14 lattice —
**7 of 14 directions already do not move the region**. The region is not dimension-starved, it is
p-starved. (Measured at |W| in {0,1}; an indicator for |W| >= 2, not a proof.)

## The carriers: the polynomial framing gives NO purchase
Of the nine `{22229,22230,22231,35758..35762,37887}`, only **37887** is structurally forced to be a
check (a square, so no bare variable, so no orientation can make it a definition, degree 4, 1
equation). The other eight admit **1-3 legal unit targets** each and sit at degree 1-2 — and **3 of
33 region atoms have no unit target**, so "no legal target" does not pick out the carrier set.
**REFUTED as an explanation.** Their only shared distinction here is high equation incidence (6-10
equations each), which points back at my own "price in equations, not atoms" rule.

## Re-entry (post-restart)
    cd solve_lab/agentN_work
    python3 fwd2.py                  # ONLY if fwd2.pkl is missing again
    python3 polyexact.py narrow      # degree-1 proof for the 7 narrow knobs
    python3 polyfull.py              # exact polys + verification against the frame
    python3 poly68b.py               # saturation loop + Singular dim/degree
    python3 solve68.py               # the exact integer solve: OPT 5, 0/924
    python3 pcert.py                 # the mod-p certificate
    python3 pgrow.py                 # why collateral budget does not touch it
    python3 budget68b.py w2 <s> 4    # resumable |W|=2 shard

## The rank-gap experiment: RUN, and the gap is an INVARIANT (step 18)

Post-restart integrity first: `*.pkl` is globally gitignored, so **both `fwd2.pkl` and `model.pkl`
were wiped**. My chain **hard-failed** on the missing pickle before any measurement, so nothing I
reported came from a partial cache. Rebuilt from `EQUATIONS.txt` and verified faithful, not merely
runnable: 42,267 atoms / 39,033 equations; 30,001 defs / 12,266 checks / 8,747 free; `frameB.py`
reloads a known state to **score 39,026** with the exact 7 nonzero atoms, the identical failing set
and **0 vars differing**; `optN.py` calibration reproduces `|R|=12 OPT=5 score=39026`.

**All 16 detach states — the whole 2^65 lattice by proof — priced exactly (`pgap.py`):**

| class | \|R\| | knobs | lattice | rk_Q(M) | gap_Q | rk_p(M) | **gap_p** | OPT | score |
|---|---|---|---|---|---|---|---|---|---|
| 8 states **with** `28730` | 12 | 68 | 14 | 7 | **0** | 3 | **1** | 5 | **39,026** |
| 8 states **without** `28730` | 13 | 76 | 15 | 8 | **0** | 4 | **1** | 6 | **39,026** |

**`gap_Q = 0` and `gap_p = 1` in all 16.** With `pgrow.py` (gap 1 across all 15 lattice-enlarging
`|W|=1` drops), the gap is invariant across the entire detach axis and across collateral budget 1.
Reproduces the corrected step-16 table exactly, which is the cross-check that the model is right.

**Rule 9 caught a real error in my own new code**: the first `pgap.py` run reported the
`28730 ∉ D` states at 39,025 with `gap_Q = 1` — the *pre-correction* number — because `price()`
kept only the constant and linear parts of each region row and so **truncated eq 8680's square
instead of rooting it**. Fixed by replacing a single-square-atom row with its `square_base`. Not
reported before it was checked.

**The placement axis is exactly those 16.** `best/new_instance_partial_39026.json`,
`N_r13_39026.json` and `H_frameB_39026.json` all load at frame score 39,026, `|R| = 12`, the
identical 7 nonzero atoms — **three independent artifacts, one configuration**. Foreign assignments
(`best_partial_3901*.json`) land at score 36,761 with `|R| = 2,273`: they are not in this frame's
coordinates and pricing them would measure my re-derivation, not their configuration.

## Single next experiment
The gap is invariant on everything reachable from this frame, so **the target is now a configuration
this frame cannot reach**. The gap is a property of `p` and the region response; `p` enters through
the frame's constants, and the one axis never varied is the **selector setting** — every
configuration I priced inherits the witness's selectors. So: find a selector setting whose region
response is not rank-deficient mod `p`. It is still one rank computation per configuration, and it
is the only remaining input to `rk_p(M)` that I have held fixed throughout.

## Dropped deliberately
`|W| = 2` (26,565 pairs, ~63 CPU-hours under contention) — its `|W| = 1` sibling is refuted over all
231 rows and `pgrow.py` explains why the budget cannot touch the gap. **Left resumable at
2,004/26,565, max g = 5**: `python3 budget68b.py w2 <shard> 4`.

---

# POST-RESTART-2 (step 19): environment rebuilt and re-verified; THE SELECTOR AXIS

## Environment after the second restart — checked before anything was measured
`*.pkl` wiped again **and site-packages wiped again** (no sympy, no python-flint). Reinstalled
`sympy 1.14.0` + `python-flint`; rebuilt `fwd2.pkl` and `model.pkl` from `EQUATIONS.txt`.
Verified FAITHFUL, not merely runnable, to the same five-point standard:

| check | value | matches published |
|---|---|---|
| `model.get()` | 42,267 atoms / 39,033 equations | yes |
| `fwd2.py` | 30,001 defs / 12,266 checks / 8,747 free | yes |
| `frameB.py` known state | 39,026, nz `[22229,22230,35758..35762]`, failing `[12231,12270,12350,14584,18673,22044,29125]`, **0 vars differing** | yes |
| `optN.py` calibration | `|R|=12 |S|=8 knobs=7 rank=7 z0=5 OPT=5 outside=0 failing=7 score=39026 exhaustive=True lin=True` | yes |
| `solve_lab/checker.py best/new_instance_partial_39026.json` | **39026/39033**, identical failing set | yes |

My chain again **hard-failed** on the first missing artifact (`kerquad.py` → `runs/polyfull.pkl`)
rather than degrading — `ikc.py` lifts `int_kernel_columns` VERBATIM out of `kerquad.py` so the
rank measurement does not need that pickle.

## New tooling
- `seltree.py` → `runs/seltree.json`: the selector hierarchy recovered **from my own frame**, with
  no dependence on another agent's chain. Each defined variable's Frame support restricted to the
  256 selectors gives 784 distinct blocks; they are 99.5% laminar (1,618 crossing pairs of 306,936).
- `psel.py`: structural configuration generator + size probe (`runs/psel_size.json`).
- `pselrank.py`: the rank measurement, exact linear part by Newton interpolation.

## The selector axis — construction (two independent regimes, both self-derived)
The 256 selectors are pure free inputs of `Frame(POOL)`. The deliverable has **exactly two live**,
`{2081, 24601}`. `runs/seltree.json` shows the hierarchy is **253 selectors meeting at one root
plus THREE outliers `{2081, 4287, 13195}` that never join it** — so the deliverable's live set is
structurally a MIXED pair: one outlier plus one leaf at depth 10 inside the tree.

- **regime "pinned"** (`psel.state_for`): selectors set, every other free input left at the
  witness. The leaf pins of a newly live leaf then break by construction.
- **regime "consistent"** (`pselc.state_for`): every leaf pin is `sel*(w - C) - m*z` with **512 of
  512 wires `w` FREE** and **512 of 512 `z = a*b` carrying a FREE factor**, so setting `w := C` on
  the live leaves and zeroing one factor of every `z` satisfies **all 512 pins simultaneously for
  ANY live set** — verified `badpins = 0` in all 51 configurations. *The pins never obstruct a
  selector setting.* (U's pin result, reached independently from my own parse.)

Knob set in both regimes = every free input syntactically supporting an atom of the region, so the
**knobs re-tune everything else**; only the selectors are pinned. `pselrank.py` takes the EXACT
linear part (Newton interpolation on t = 0..6 wherever the second difference is nonzero — a plain
secant is wrong in exactly the entries where the response is quadratic), and reports both the
AMBIENT ranks (all knobs free) and the LATTICE ranks (integer kernel of the linear collateral
response). The lattice is a RELAXATION of `pgap.py`'s exact saturation lattice, so
`gap_p >= 1` there implies `gap_p >= 1` exactly — the safe direction for a hunt.
**Calibration: the deliverable's own setting gives `LAT dim=15 rk_Q=8 rk_p=4 gap_Q=0 gap_p=1`,
i.e. exactly `pgap.py`'s `gap_p = 1` with every rank shifted +1, which is what a relaxation must do.**

## In flight at time of writing (resumable, PIDs in `runs/pselrank.pid`)
    python3 pselrank.py pselrank_s<S>.jsonl <S> 4               # pinned regime, 4 shards
    python3 pselrank.py pselrankC_s<S>.jsonl consistent <S> 2   # consistent regime, 2 shards
    python3 pselsum.py [pinned|consistent]                      # the gap_p distribution
Both skip any tag already recorded IN THEIR OWN REGIME, so a restart resumes for free.
51 structural configurations per regime (`psel.configs()`), ~200 s each under fleet contention.
