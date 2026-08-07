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

## 5. Single highest-value next experiment
The four blocking residues are `(≈2,420-bit numerator) mod p`.  The numerators depend on the
selector configuration; **p does not**.  So: for each of the 2,800 (a-bit, b-bit) pairs drawn
from the 106 proven bits, construct the witness analogue — agent H's `frameB.Frame([642,
28730,29854,31864])` reproduces the witness bit-for-bit and is the right tool, *not* E's
`forward` — and test whether the 13-equation region system is integrally solvable there.
That test is ~2 s per configuration once the state is built.  A configuration in which the
four numerators are ≡ 0 mod p gives **all 39,033 equations**, because the private variables
cannot break anything outside E(R).

Second choice: enlarge the private set beyond the cone — look for regions elsewhere in the
instance whose private variables have columns **coprime to p** (all four of this region's
p-multiplier atoms are the obstruction; a region built on unit-coefficient atoms would not be).

## 6. Do not redo
- Singles over all 106 representatives at cfg0 (done, negative, `runs/singles.jsonl`).
- Deeper closure on E's residual (`maxr` 3→6, `maxv` 2000→8000, knob set unfiltered, lazy row
  activation): never below 28 failing.  The knobs are not the problem — the *cone* is.
- Single- and double-atom growth of the witness region: all fail on eq29125 mod p.
- (a,b) pairs through `full11.solve_pair`: the first 68 all give **39,013** with residual
  exactly `{20649,20652,32148}`, independent of b.  Low information at that granularity.
