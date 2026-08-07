# RESUME_O — agent O.  Representative sweep from the 106 proven bits; the channel tension.

## 0. Verification rule used here
All my states have values ≤ ~2,450 bits (≈740 digits), so **plain `solve_lab/checker.py`
parses them** — no need for `verifyE.py`.  Every score below was produced by
`checker.py <file>` or by `engine.eqfails`, and the two agree wherever both were run.

## 1. Scores
- Deliverable re-verified by me: `checker.py solve_lab/best/new_instance_partial_39026.json`
  → **satisfied 39026/39033**, failing `[12231,12270,12350,14584,18673,22044,29125]`.
- **My best: 39,026** — I did not beat it.  Two independent reconstructions of that score
  from my own region model, both checker-verified:
  `agentO_work/grow23618_39026.json` and `agentO_work/region_opt_39026.json`
  (the second uses *different* values for all seven region variables — 3013b vs 2400b for
  x_642 etc. — so it is a genuinely distinct point at the same score).

## 2. What I ran (all files in agentO_work/)
| script | what |
|---|---|
| `sweep.py`,`sweep_run.py` | complete singles sweep: all **106** proven pin-solvable bits as channel representatives at E's cfg0, through E's `channels.simsolve` → `runs/singles.jsonl` |
| `pairsweep.py` | (a-bit, b-bit) pairs through E's `full11.solve_pair` → `runs/pairs.jsonl` (68 of 2,800 done before I stopped it — the answer was already constant) |
| `simO.py`,`lazy.py` | simultaneous solve with the **unfiltered** knob set (E's `affine_cols` rejects a knob if it is non-affine on *any* atom; I drop non-affine (knob,row) pairs instead) + lazy row activation |
| `diag3/5/6/7/8.py` | per-row knob-reachability, deliverable-vs-E-frame diffs, DAG-divergence roots |
| `coset.py`,`regionopt.py` | the witness's residual as 12 equations over 8 atoms; exhaustive subset optimisation |
| `regiongrow.py`,`growlist.py`,`grow23618.py`,`qbound.py`,`snf.py`,`fullreg.py`,`consts.py` | region growth, rational bound, Smith/echelon obstruction, full-system integral solvability |

## 3. Answers to my two assigned questions

### 3a. Representative sweep — RUN AND COMPLETE, negative
All 106 proven pin-solvable bits, one at a time, at cfg0, through E's §16 simultaneous solve
(`maxr=3, maxv=2000`, knob set = E's affine-filtered closure):
**no representative beats the empty set.**  Best per group: channel-0 (a-tree) 38,989
(x_22492); channel-1 38,995 (x_13710); channel-2 38,977 (x_490); the 16 "inert" bits
38,928–38,953.  x_1530/x_1603 return 39,005 because they are already on.
So "E only tried 2 representatives per channel" is **not** the explanation for monotonicity.

### 3b. The monotonicity TENSION — resolved, with a mechanism
Three separate facts, each measured:

1. **cfg0 is not "no channels live"; it is the (0,1) branch.**  At `triple8_seed`,
   `x_7715 = 0`, `x_34554 = 1`, so `x_15298 = x_7715*x_34554 = 0`.  Atoms
   `20649 = x_15298*x_11150`, `20652 = x_15298*x_25739 - 6672769*x_29804`,
   `32148 = 537773*(x_15298*x_37758) - x_35605` are therefore **vacuous** at cfg0.
   Turning on *any* a-tree bit sets `x_7715 = 1`, hence `x_15298 = 1`, and those three rows
   go live at once.  That single gate is the whole monotone cost — it is not per-bit
   accumulation.  This is why every a-tree single leaves `{20649,20652,32148}` bad.
2. **The three activated rows are knob-starved, and that IS a real measurement.**
   At the (1,1) pair state (a=24601, b=2081) with knob set = **all 3,545 non-boolean free
   variables in the 6-round cone closure**, the only knobs with nonzero delta on those rows
   are `x_22162`, `x_30213` (coefficients coprime to p) and `x_5146`, `x_2936` (coefficients
   ≡ 0 mod p).  Three conditions, two mod-p degrees of freedom — and `x_22162`, `x_30213`
   are already consumed by the MUX closure.  Knob set and configuration stated.
3. **The branch mismatch.** The 106 bits' pin systems were proved solvable by `bitfeas2`
   at base `{18956:C, bit:1}` with `CORE={20212,20215,24403,747}` excluded — i.e. in the
   **single-tree branches** (1,0) for the a-scan and (0,1) for the b-scan.  Plugging them
   into cfg0 / (1,1) is a different branch, which is why their pins are never repaired there.

### 3c. …and the deeper reason E's line cannot reach the deliverable
**E's forward map cannot represent the 39,026 witness.**  Feed the witness's own free-variable
values through `engine.forward` and you get a *different* state: 23 variables differ, and the
divergence has exactly four roots — `x_31864`, `x_29854`, `x_642`, `x_28730` — whose defining
atoms (36663, 36659, 36664, 23617) the witness deliberately **violates**.  E's forward always
satisfies a defining atom, so it cannot express the witness; that state scores 25 failing
equations under E's orientation, not 7.
Consequently the repair the witness uses is **outside the dependency cone of the residual**:
cone(20649,20652,32148) has 277 free variables and **none** of the witness's six carriers
(`x_1329, x_8731, x_9118, x_9413, x_10903, x_17325`) is among them; each carrier's +1 probe
has *zero* delta on all three rows.  No cone-generated closure — at any `maxr`/`maxv` — can
reach it.  So the enumeration was attacking a residual whose repair its knob generator cannot
express.  Raising `maxr`/`maxv` is not the fix (tested: `maxr=6, maxv=8000`, 3,511–3,552 knobs,
lazy row activation over 30 iterations — the residual moves but never drops below E's 28).

## 4. The new object: the witness's residual is a 13-equation, 8-unknown integer system
- The witness's 8 bad atoms `{23616,23617,36659,36660,36661,36662,36663,36664}` are touched by
  exactly **12 equations**; adding atom `a23618` makes it 13 (the extra one is **eq8680**).
- **Seven variables are private to the 8-atom region** (occur in no atom outside it):
  `x_642, x_1329, x_9413, x_10903, x_17325, x_29854, x_31864`; `a23618` frees an eighth,
  `x_28730`.  Private ⇒ moving them cannot break any other equation, so
  **failing equations = |E(R)| − maxsat(R)** exactly (verified against `eqfails`).
- Atoms 36662 and 36663 have **identical coefficient columns in all 12 equations**, so the
  residual matrix has rank 7, not 8.
- **Exhaustive** over all 4,095 subsets of the 12 (`regionopt.py`): exactly **one** is
  integrally satisfiable — the witness's own `{2554,6816,8124,9123,9421}`.  With `a23618`
  absorbed, maxsat rises 5 → 6 but |E| rises 12 → 13: a **1-for-1 trade**, independently
  reproducing agent H's eq8680 result in a different formalism.
- **Over Q the system is fully solvable**: rank 8 = #unknowns, the 5 dependent rows are
  exactly consistent (0 = 0), so the rational solution is unique and satisfies all 13.
- **Over Z exactly four divisibilities block it** (`snf.py`): three by **p** (the 256-bit
  constant) and one by a 279-bit pivot (gcd 2 with its rhs → a 278-bit modulus).
- **Why growth cannot help** (`fullreg.py`, `consts.py`): `x_17499 = x_22665 = x_28961 =
  x_28599 = p` **exactly**, and every adjacent atom that frees a new variable has the form
  `x_t − p·x_new`, so that knob's whole column is divisible by p.  All **39 single** and all
  **741 double** growths fail on the *same* row, eq29125, by the same factor of p.

  *(Knob set: the ≤10 variables private to R0 ∪ {a23618} ∪ {≤2 adjacent atoms}.  Selector
  configuration: the witness's, `x_2081` and `x_24601` on, all other cluster booleans off,
  all non-private variables at the witness's values.  Outside that knob set I claim nothing.)*

## 4b. THE RATE — the configuration scan is dead twice over, but the conditions INVERT
Coordinator asked for the expected hit rate before spending cores.  Here it is, measured.

**Correction to §4: there are FIVE blocking coordinates, not four** (`rate.py`).  The unique
rational solution has denominators `2458959, p, p, p, 2458959·p` on
`x_642, x_1329, x_9413, x_10903, x_17325` — i.e. **4 conditions mod p and 2 mod 2458959**
(= 3 × 819653; note 7376877 = 3 × 2458959 is the literal in atom 23616).

### Kill 1 — the rate (`hitrate.py`)
Admissible boundary changes form a coset `δ0 + Λ0`.  Measured period in each direction with an
exact solvability oracle: **p** for `const(a23618)`, `const(a36660)`, `const(a36662)`; for
`const(a23616)` the period exceeds every modulus tested (up to 2458959·p).  So
`[Z⁴ : Λ0] ≥ 2^768` and the **scan hit rate is ≈ 2⁻⁷⁶⁷**.
Expected hits in 2,800 configurations ≈ **10⁻²²⁷**; in all 13,884 ≈ 10⁻²²⁶.  Not a search.

### Kill 2 — zero variance, which is worse than a bad rate (`bscan.py`)
The four quantities a configuration would have to move are
`K1 = x_7068−x_2099`, `L = x_4432−x_19964`, `K2 = 5113045·x_9118`, `J = x_7075·x_8731`.
Across **35 configurations** (empty, 12 a-bits, 12 b-bits, 10 (a,b) pairs) in E's frame all four
are **identically 0** — 1 distinct value out of 35, for every one of them.  The scan would have
measured a single point 2,800 times.  **They are assignment knobs, not configuration knobs**:
the witness has `J = 2428 bits` and `K2 ≠ 0` because it *assigns* the free variables
x_8731 / x_9118, not because of its selector pair.

### The positive result — the conditions are invertible, and I have the target
`invert.py` / `tunable.py` / `target.py`: put the boundary shift δ into the unknown vector and
solve `A z + B δ = b0` over Z.
- δ = 0: unsolvable (the witness).  **0 of 9** single supports, **0 of 36** pairs, **0 of 84**
  triples work; **12 of 126** quadruples do — including `{a23616, a23618, a36660, a36662}`,
  i.e. **exactly the four constants that are not p-multipliers**.
- Applying that δ0 makes **all 13 region equations hold** (verified end-to-end, `target.log`).
- Shifts required: 2440 / 2419 / 2428 / 2429 bits.  The a36660 shift must be divisible by the
  carrier factor 5113045; `gcd(5113045, p) = 1`, and the period there is p, so a CRT lift in the
  same class mod p always exists (computed: `x_9118 += ` a 2406-bit value).  Exact values in
  `target.json`.
- **Two of the four shifts are free**: `a36662` is carried by x_8731 and `a36660` by x_9118 —
  precisely agent H's zero-collateral knobs.  The remaining two, `K1` and `L`, are carried by
  ordinary derived variables whose collateral **must be evaluated in frame B**; my atom-level
  model holds non-private variables fixed and cannot express re-derivation, which is exactly
  why it sees 8 knobs where H sees 9.
- So the honest statement is: *if* the K1 and L shifts can be realised without disturbing
  anything outside E(R), the result is all 39,033 equations.  That "if" is the open question —
  it is not established, and I am not claiming it.

Why a move-based search would never find this: δ0 is an exact 2,429-bit lattice target.  H's
70,008 single moves and 576 zero-collateral pairs over the same nine knobs could not have
stumbled on it.  That is the difference between sampling the conditions and inverting them.

## 4c. Reconciling with agent M on eq29125 — two different claims, no contradiction
M withdrew a divisibility obstruction on eq29125: single-row solvability is `gcd(coef) | s0`,
and eq29125's gcd is 1, so all seven failing rows pass.  **That is correct and it does not
conflict with my barrier**, because we are testing different things:
- **M**: is eq29125 satisfiable *on its own*, over M's full affine knob set?  Yes — gcd 1.
- **Me**: eq29125 is the row at which the *simultaneous* elimination of all 13 region equations
  fails, over the ≤10 variables **private to my region**.  The factor p is not in the row as
  written; it appears only *after* the other twelve are eliminated, because `x_17499 = p`
  exactly, leaving coefficient p against an rhs that is not ≡ 0 mod p.
A row can be individually satisfiable and jointly infeasible; gcd-1 is a statement about the row
before elimination, my p is a statement about the pivot after it.  **My barrier is the one with
the tighter knob set, so its scope must always be carried with it** — it says nothing about
M's knob set, and nothing outside the region-private variables.
My §4b result independently *supports* M's direction: because the obstruction lives only in the
joint elimination and not in the row, changing the rhs clears it — which is exactly what δ0 does.

## 5. Single highest-value next experiment (REVISED after §4b)
**DO NOT run the 2,800-configuration scan** — I proposed it, then measured it dead twice (§4b:
rate 2⁻⁷⁶⁷, and the four boundary quantities are configuration-invariant, 1 distinct value
across 35 configurations).  It would burn 1.5 hours to sample one point repeatedly.

Instead, the target is now explicit.  In **agent H's frame B** (`frameB.Frame([642, 28730,
29854, 31864])`, which reproduces the witness bit-for-bit — *not* E's `forward`, which
provably cannot represent it):
1. Apply the two free shifts: `x_8731 += δ0(a36662)` and `x_9118 += (CRT-lifted δ0(a36660))
   / 5113045` — both are H's zero-collateral knobs, so these cost nothing.  Values in
   `target.json`.
2. The whole question then reduces to the remaining two shifts,
   `K1 = x_7068 − x_2099 += δ0(a23616)` and `L = x_4432 − x_19964 += δ0(a23618)`.
   Measure their collateral **in frame B** (each shift is free to move by multiples of p in the
   `a23618` direction, so there is a p-indexed family of representatives to choose from —
   pick one whose collateral is zero, if any exists).
3. If a zero-collateral representative exists for both, all 39,033 equations follow.  If none
   does, the cost of the cheapest representative is the new bound, and it is a *measured* bound
   rather than a search outcome.

Second choice, if that closes negative: look for regions elsewhere in the instance whose
private variables have columns **coprime to p**.  Every p-multiplier atom in this region
(`a23617, a36659, a36661, a36664`, all with multiplier exactly p) is what forces the mod-p
conditions; a region built on unit-coefficient atoms would not have them.

## 6. Do not redo
- Singles over all 106 representatives at cfg0 (done, negative, `runs/singles.jsonl`).
- Deeper closure on E's residual (`maxr` 3→6, `maxv` 2000→8000, knob set unfiltered, lazy row
  activation): never below 28 failing.  The knobs are not the problem — the *cone* is.
- Single- and double-atom growth of the witness region: all fail on eq29125 mod p.
- (a,b) pairs through `full11.solve_pair`: the first 68 all give **39,013** with residual
  exactly `{20649,20652,32148}`, independent of b.  Low information at that granularity.
