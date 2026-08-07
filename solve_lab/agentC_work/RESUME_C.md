# Agent C — RESUME.  Best verified of mine: 39,013 (agentC_work/BEST_39013.json, checker-confirmed).

Setup: `cd solve_lab/s9 && python3 atomize.py && poly.py && gates.py && fwd.py` (57 s) rebuilds the
atom/gate caches.  Installed (absent before): z3-solver, python-sat, cvc5, python-flint, ortools,
sympy, numpy.  Everything here is stated over the raw integer/polynomial content of EQUATIONS.txt.

## Model (agentC_work/supp2.py, supp3.py, lib2.py, ort.py)
* All 900 nontrivial SCCs of the definition graph are size-2 duplicate equalities, so the system is
  triangular: **8,173 free variables**, 30,575 defined ones, every defining atom linear in its own
  variable with coefficient 1.  Forward substitution is exact and costs 24 ms.
* Free variables = 0  ->  0 defining atoms violated, **only 6 of 10,792 remaining atoms nonzero**,
  score 39,005.  Those 6 collapse to three conditions (see LOG.md Step 2), all three of which are met
  EXACTLY OVER Z by setting four free variables; closing the rest reaches **39,013** in 2.3 s.

## Engines
`close.py` -> `close2.py` (re-solvable repairs) -> `close3.py` (realise a value down the definition
DAG) -> `close4.py` (frame DETACHMENT: drop a defining atom so its variable becomes a control input).
`close4` protects the 220 variables equal to the literal p from being altered during realisation.

## Placement result, and the honest state of it
* `exact10513.py`: the cheapest cluster by equation-count (built on x_10513) has a 12 x 8 incidence
  matrix of **rank 8, kernel zero**.  Five of its eight atoms are SHADOWS — fixed linear combinations
  of the other three, each alone in its own squared equation — so they inflate the atom count without
  adding freedom.  A predicted 39,027 there is **REFUTED**; true cost >= 11, measured 38,988.
* **eq20538 = 30*a8427** is a single-atom equation forcing that atom to zero, structurally identical
  to **eq29125 = a22230** in the deliverable's cluster.  Every cheap cluster found so far carries one.
* `minweight.py` replaces my earlier "settable atom" classifier with a classifier-free quantity: the
  minimum weight of `M v` over nonzero rational `v`.  **CALIBRATION: it returns 5 on the deliverable's
  own cluster, whose true cost is 7.**  It is therefore a sound LOWER bound but loose by exactly the
  congruence term: every handle in the file is `p*(free)`, so each atom value is confined to a fixed
  residue class mod p, and a rational kernel direction is generally not realisable.
  **=> No purely structural classifier can reproduce 7.**  The missing 2 is not a property of the
  incidence matrix; it depends on which literal constants land in the residue vector.
  Consequence: defect placement is **not closed** — but the bound is also very loose (4 against a
  true >= 11 on x_10513), so ranking clusters structurally has low predictive value.

## Highest-value next experiment (for whoever takes this on)
Make the pricing residue-aware instead of structural: for a candidate cluster compute the residue
vector `v0 = (atom values mod p)` from an actual construction, then a subset `T` of the cluster's
equations can vanish only if `M_T v0 = 0 (mod p)`; maximise `|T|` subject to that plus the integer
solve over the p-lattice.  That is the quantity that equals 7 on the deliverable, and it is the only
version worth scanning.  It needs one construction per cluster, so pick candidates by
`|E| - |S|` first (agentC_work/truecost.json, clusters.json) and price only the top few dozen.

## Commands
python3 checker.py agentC_work/BEST_39013.json    -> 39013/39033
python3 agentC_work/exact10513.py                 -> the 12 x 8 matrix and its rank
python3 agentC_work/minweight.py                  -> calibration + classifier-free scan
