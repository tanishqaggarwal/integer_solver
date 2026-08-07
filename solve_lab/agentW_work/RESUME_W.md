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
