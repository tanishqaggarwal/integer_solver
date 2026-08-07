# RESUME_W — agent W.  IN PROGRESS.

**Best verified score: 39,026 / 39,033** — `solve_lab/best/new_instance_partial_39026.json`,
re-verified from cold with `solve_lab/checker.py`, failing
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`.  **I did not beat it.**

Inherited agent O's closed thread.  Two tasks: finish O's frame-B budget (TASK 1) and settle
T's frame-B flag on `|K| = 34` (TASK 2).

## Environment I built (agentW_work only; agentH_work never written to)
`model.py`, `fwd2.py`, `frameB.py` are byte-copies of agent H's; run from `agentW_work` they
build `model.pkl`/`fwd2.pkl` **here**, so H's directory stays untouched (H had no `.pkl` at all —
importing `frameB` from `agentH_work` would have created them there).  My `model.pkl` and
`fwd2.pkl` are **md5-identical to agent N's**, so this is H's model, not a re-derivation.
`PYTHONDONTWRITEBYTECODE=1` throughout.  No git commands.

## TASK 2 — SETTLED.  O's `|K| = 34` is CORRECT in frame B.
`python3 w_K.py` (and `w_K_default.py` for the contrast).

Frame B = `frameB.Frame([642, 28730, 29854, 31864])` reproduces the witness bit-for-bit:
score **39,026**, same 7 failures, **0 of 38,748 variables differing**, 7 nonzero check atoms
`[22229, 22230, 35758, 35759, 35760, 35761, 35762]`.

| orientation | free inputs | witness score | nonzero atoms | \|U\| | \|C\| | \|U∩C\| | **\|K\|** |
|---|---|---|---|---|---|---|---|
| **frame B** `[642,28730,29854,31864]` | 8,751 | 39,026 | 7 | **15** | **26** | **7** | **34** |
| default (no detach) | 8,747 | **39,020** | 5 (incl. a37887) | 30 | 26 | 26 | 30 |

O's three numbers — 15 free inputs reaching a region atom, 26 free carriers of `S`, union 34 —
**all reproduce exactly in frame B**.  T's rebuild (12 / 11 / 23, overlap 0) is a fact about
**F's parse in the default orientation**, which is a different point: the witness's free values
pushed through the default DAG score **39,020**, not 39,026, and its nonzero atom set is a
different set of 5 (including `a37887` itself, i.e. eq8680 fails there).  So T's numbers and O's
are measurements of different objects, exactly as O's scope line said.
**Ledger row → CONDITIONAL, scope: agent H's model (`model.py`/`fwd2.py`), frame B's
orientation.  Not a defect.**  Note `|C| = 26` is the same in both orientations — the
frame-dependence is entirely in `U`.

## TASK 1 — in progress.  See below.  Setup reproduces O's exactly:
`python3 w_setup.py` → knobs **34**, rows **175**, satisfied **168**, failing **7**,
nonlinear rows dropped **16**, reachable checks **64**, reachable equations **190**,
S-row support **17**, `S0 = 0`.  Identical to `runs/fb_j3.log`'s header.

### Correction to O's own budget statement
O's `T_COMPENSATION.md` says the 14 completed triples are "every triple containing eq12231".
There are **C(6,2) = 15** such triples; O's list has **14**.  `[12231, 22044, 29125]` was **not
reached** either.  So the unfinished set is 21 triples, of which 20 avoid eq12231 and one
contains it.  Does not change any O conclusion.

### The solver oracle is sound (`w_oracle.py`)
`solve_sparse` returns `None` for five distinct reasons, two of which are give-ups
(`core too large`, `coefficient blowup`) and would NOT be negatives.  Over 527 solves spanning
all 127 bought-sets and a random b=2 sample, **every** `None` was `core infeasible` — i.e. FLINT
HNF proving no integer solution.  O's negatives rest on a sound oracle.

## TASK 1 — RESULTS

### (a) The 21 unreached triples at j=3, b<=2 — same method as O (`w_j3.py`, `w_j3.log`)
Running / complete; see `w_j3.log`.  Exact integer `solve_sparse`, ~14,198 solves per triple,
b=0 then all 168 b=1 then all C(168,2)=14,028 b=2.  **No improvement in any triple reached.**

### (b) The linear model's pricing is EXACT, verified OUTSIDE the model
`w_dump.py` dumps the "buy eq12231, break eq2554" solution as a full 38,748-variable
assignment; `solve_lab/checker.py` returns **39,026/39,033 failing
`[2554, 12270, 12350, 14584, 18673, 22044, 29125]`** — exactly what the model predicted
(`w_trade_12231_break2554.json`, 27 variables differ from the witness).  So frame-B model
predictions of score and failing set are trustworthy at this radius.

### (c) AUDIT of O's collateral accounting — CONFIRMED, with one number corrected
`w_audit.py` / `w_audit.json`.  Each of the 7 failing rows was bought against each of the 6
essential satisfied rows and the result priced **exactly through `frameB.State`**, not through
the linear model.  **32 of the 42 combinations are feasible and every one lands on exactly
39,026 with exactly the predicted failing set** — one row in, one row out, no collateral
anywhere.  O's collateral accounting is sound.

**But O's "every purchase costs exactly eq8680" is too strong.**  The trade is not seven-way,
it is **32-way**, and eq8680 is only one of six possible prices:

| bought | prices that work (1-for-1, all exactly 39,026) |
|---|---|
| eq12231 | 2554, 6816, 8124, 9123, 9421, **8680** |
| eq12270 | 2554, 6816, 8124, 9123, **8680** |
| eq12350 | 2554, 6816, 8124, 9123, 9421, **8680** |
| eq14584 | 2554, 6816, 8124, 9123, 9421, **8680** |
| eq18673 | 2554, 6816, 8124, 9123, 9421, **8680** |
| eq22044 | 2554, **8680** |
| eq29125 | **8680** only |

Only **eq29125** has eq8680 as its unique price.  This also explains T's unaudited flag that
"7 knobs move a failing row with `dS = 0`": there are trades that never touch `S` at all.
*The conclusion — every trade is 1-for-1, no gain — is untouched.  The mechanism claim was
narrower than the truth.*  **O's Lemma (`S = 0` forced) is completely untouched by this**; it
is about which assignments exist, not about which row you pay with.

### (d) THE STRUCTURAL RESULT — the frame-B obstruction is ARITHMETIC, not linear
`w_struct.py`, `w_essential.py`.  Exact rational arithmetic over the 175x34 system:

- **rank([A | b]) over all 175 rows = 28, and the rhs column is NOT a pivot** — so the FULL
  system, *including all seven failing rows*, is **consistent over ℚ**.  Over ℚ one can buy all
  7 for free.  **Every part of the frame-B negative is integrality, none of it is rank.**
  (This is the local, frame-B counterpart of O's §4 "over ℚ unique solution, over ℤ five
  coordinates blocked", now measured over all 175 rows rather than the 13 region rows.)
- rank(A_SAT) = 26, so the deltas keeping all 168 satisfied rows satisfied form a **rank-8**
  lattice; the 7 failing rows add only **2** further dimensions.
- The 168 satisfied rows are **homogeneous** (rhs = 0 for every one).  Therefore the admissible
  integer deltas for a kept-set KEEP are exactly `ker_Z(A_KEEP) = Z^34 ∩ ker_Q(A_KEEP)`, so
  **breaking rows that do not drop rank_Q changes nothing at any budget.**
- Exactly **6 SAT rows are essential** (single deletion drops rank 26 -> 25):
  **`{2554, 6816, 8124, 9123, 9421, S}`** — the region's own equations plus the `S` row.
  The other 162 are individually redundant.

### (e) EXHAUSTIVE over the essential-row break family, all j = 1..7 (`w_exhaust.py`)
All 2^6 = 64 subsets of the essential rows x all 127 bought-sets, exact integer oracle:

> **minbreak(P) = |P| exactly, for every P with 1 <= |P| <= 6.  GAIN = 0 everywhere.
> All seven together are NOT buyable at any b <= 6.**

Perfect k-for-k at every k.  30 seconds, vs 33 minutes for O's 14 triples.

### (f) HONEST CORRECTION to my own step (e)
I first claimed redundant-row breaks were worthless outright.  **False.**  `w_pack.py` packs
only **t = 2** disjoint rank-20 subsets out of the 162 redundant rows, and the coordinate census
shows functionals supported on as few as **2** rows — so small cocircuits *do* live inside the
redundant rows and deleting 2..6 of them CAN enlarge the lattice.  (Only rank-dropping
deletions matter; that part stands.)  Hence the cocircuit enumeration in (g).

### (g) Cocircuit enumeration — how far the closure actually reaches (`w_cocirc.py`, `w_close.py`)
The code `c -> (row_e . c)_e` has rank 26 and length 168.  Fix an information set `I` (26 rows,
identity block).  Any break-set `T = N(x)` with `|T| <= 6` is a **union of minimal cocircuits**
(proof in `w_close.py`'s docstring: if `T` properly contained the union `T'` of the minimal
cocircuits inside it, the two kernels would coincide and `x` would satisfy `T \ T'`).

- `s = |T ∩ I| = 1` and `s = 2` are enumerated **EXACTLY** (no window, no degeneracy: `s=1` is
  the basis codeword itself; `s=2` is the most-frequent-ratio over all 142 outside columns).
  A cocircuit of size <= 2 has `|T ∩ I| <= 2`, so **all cocircuits of size <= 2 are complete.**
- `s = 3..6` uses an adaptive column window and is a **bounded search, not exhaustion**:
  3,070,206 subsets were skipped as degenerate.

Result: 5,419 candidate supports, **70 minimal** — sizes `{1:6, 2:1, 3:2, 4:5, 5:14, 6:42}`.
The size-1s are the 6 essential rows.  The **only** minimal cocircuit of size 2 is
**`{22563, 8687}`**, which contains no essential row — this is the concrete witness that (f)'s
correction was necessary.  Union closure to size <= 6: **520 candidate break-sets**.

## SCOPE, stated as budget not exhaustion
| claim | status |
|---|---|
| j=1, b=0 | **exhaustive** (O, reconfirmed) |
| j=2, b<=1 | **exhaustive** — every break-set of size <= 1 that drops rank is an essential row, and all 6 were tested |
| **j=3, b<=2, all 35 triples** | **exhaustive** — all cocircuits of size <= 2 are enumerated exactly, so every rank-dropping break-set of size <= 2 was tested.  This supersedes O's 14 triples *and* my brute-force 21 |
| j=3, b<=2 by O's own brute force | O: 14 triples; **W: the remaining 21** (`w_j3.log`), independent agreement |
| j=1..7, breaks restricted to the 6 essential rows | **exhaustive**: minbreak(P) = \|P\| exactly, gain 0 (`w_exhaust.py`) |
| j=4..7, b<=6, general breaks | **budget, not exhaustion**: 520 union-of-minimal-cocircuit break-sets tested (`w_close.py`); the size-3..6 cocircuit search skipped 3.07M degenerate subsets |
| everything above | scope: **34 of 8,751 free inputs, frame B, agent H's model** |

## Do not redo
- The ℚ relaxation of the frame-B system.  All 175 rows are simultaneously ℚ-satisfiable, so
  every ℚ feasibility test is vacuously YES and no ℚ/LP relaxation can ever prune anything here.
- Brute-forcing break-sets that do not drop `rank_Q`.  162 of the 168 satisfied rows are
  individually redundant; deleting them changes `ker_Z` not at all.  O's C(168,2)=14,028 inner
  loop per triple is doing ~14,000 solves where **16** distinct lattices exist.
- Reading a greedy net-zero as a negative (O's warning, still correct).

### (h) RESULT of the closing test (`w_close.py`, `w_close.json`)
520 candidate break-sets (all unions of minimal cocircuits, size <= 6) x all 127 bought-sets,
exact integer oracle, **6,806 solves in 49 s**:

> **BEST GAIN = 0.**  For every bought-set P and every rank-dropping break-set B of size < |P|
> tested, `(SAT \ B) + P` is integer-infeasible.  Nothing beats 39,026.

Combined with (e): `minbreak(P) = |P|` exactly at every size 1..6, and the full seven is not
buyable at any b <= 6.  **The frame-B region is a perfect k-for-k trade at every k.**

## Re-entry
```
cd solve_lab/agentW_work
python3 w_K.py          # TASK 2: |K| = 34 in frame B                          (~10 s)
python3 w_setup.py      # O's setup reproduced: 34/175/168/7/16                (~5 s)
python3 w_audit.py      # the 32-way 1-for-1 trade, priced exactly             (~2 min)
python3 w_essential.py  # the 6 essential rows                                 (~3 min)
python3 w_exhaust.py    # exhaustive over essential breaks, all j              (~30 s)
python3 w_cocirc.py     # cocircuit enumeration                                (~4 min)
python3 w_close.py      # the closing test                                     (~1 min)
python3 w_j3.py         # the 21 triples O never reached, O's own method       (~30 min)
```
Artifacts: `w_trade_12231_break2554.json` (checker-verified 39,026, a *new* point, 27 vars from
the witness), `w_K.json`, `w_struct.json`, `w_essential.json`, `w_exhaust.json`,
`w_cocirc_raw.json`, `w_close.json`, `w_audit.json`.

## What is left, honestly
1. The `s = 3..6` cocircuit search skipped 3.07M degenerate window subsets.  Closing that —
   e.g. multiple information sets and randomised windows, or a proper minimum-weight-codeword
   routine on a rank-26 length-168 code — would upgrade j=4..7 from budget to exhaustive.
2. Everything here is **inside K**.  Per O's §9 and confirmed by (d), any improvement must come
   from the other **8,717** free inputs.  O's Lemma constrains none of them.
