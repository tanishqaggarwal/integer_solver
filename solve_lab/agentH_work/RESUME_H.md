# Agent H — RESUME.  Integer/polynomial analysis of EQUATIONS.txt.

## Best verified score
**39,026 / 39,033** — `solve_lab/best/new_instance_partial_39026.json`, re-verified by me with
`solve_lab/checker.py` (failing [12231,12270,12350,14584,18673,22044,29125]).  I did not beat it.
I *did* reconstruct it exactly inside my own frame (see frame B below), which makes it attackable
with 1,478 free inputs the lab's earlier frame did not expose — and that attack came back negative.
My own from-scratch best is 39,018 (`scan1.py`, any single selector).

## Durable assets (all kept)
- `model.py` — independent parse, 42,267 atoms / 39,033 equations.
- `fwd2.py` — **the frame**: orient every atom of syntactic form `x_t - rest` as `x_t := rest`.
  ACYCLIC over all 38,748 variables; 8,747 free inputs, 30,001 defined, 12,266 checks, NO cycles.
  Forward eval from all-zero free inputs = **39,005** with only 5 nonzero check atoms.
- `frameB.py` — parametrised frame; `Frame([642,28730,29854,31864])` reproduces the 39,026 witness
  bit-for-bit (0 variables differing) with 8,751 free inputs.
- `fast.py` / `frameB.State` — incremental exact evaluator, 1.4 ms per move.
- `close2.py` — bottom-up cascade closer: closes the entire pin tree in ~32 exact assignments with
  NO search (smallest-atom-support first, freeze each assigned variable).  `close3.py` adds an
  indirect phase.  This constructs rather than searches, and strictly dominates beam search (39,012).
- `DECOMPOSITION.json` — the graph census.

## Established, in raw terms
1. **Decomposition hypothesis REFUTED** (my own angle): eq-var graph is one component; the
   free-input hypergraph is one component (8,747 vars / 18,248 equations); 20,785 equations are
   identically satisfied by forward eval; the residual closure is 6,007 free inputs / 9,244
   equations.  No separator, no block decomposition exists.
2. **Residual = three conditions**: `x_9274 = 1`, `x_37892 = C1 (mod p)`, `x_13682 = C2 (mod p)`,
   p = 115792089237316195423570985008687907853269984665640564039457584007908834671663.
3. **Combination law of the 512 load-pin constants, MEASURED**: for two active selectors the
   delivered pair satisfies, over Z/p,
   `(X1+X2+X3+K)(X2-X1)^2 = (Y2-Y1)^2` and `(Y3+Y1)(X2-X1) = (Y2-Y1)(X1-X3)`,
   K = 97553848499418123410591666447050222001188385549510401465815187079080512838891.
   56/56 confirmations with a 1-off control failing 56/56.  Over Z one extra obligation appears,
   `2264251 | (x_15286/p)` with 2264251 = 11*43*4787; solving it for one lift parameter closes all
   three checks EXACTLY over Z (3/3).  **The constants combine nonlinearly — no linear relation
   among them can express the residual.**
4. **Integer-relation search over the constant table: negative.** 512 constants, all distinct mod p,
   gcd 1; no single, no pairwise sum, no pairwise difference equals C1 or C2 (or C1±K, C2±K) mod p.
   The additive formulation would sit at density 1.000, where lattice relation-finding does not apply.
5. **39,026 is locally rigid, confirmed from an independent orientation**: region = 12 equations,
   5 satisfied; exactly one atom (22231) lies wholly inside it and no zero-collateral knob moves it;
   exactly nine zero-collateral knobs exist {642,1329,8731,9118,9413,10903,17325,29854,31864}.
   70,008 single moves and 576 zero-collateral pairs: no improvement.

## Re-entry
cd solve_lab/agentH_work
python3 model.py; python3 fwd2.py; python3 support.py     # rebuild caches (~1 min)
python3 frameB.py     # reproduce the 39,026 witness exactly in the detached frame
python3 close2.py <u-selector> <w-selector>   # construct + close a branch with no search

## 7. Rank-raising sweep — RUN, terminates at collateral budget 1 (LOG.md Step 10)
Crossover: one extra knob raises rank by <=1, so it pays only if it drags in ZERO new equations.
Census over all 8,751 free inputs: the +0 class is EXACTLY the 9 knobs already in the lattice;
the next cheapest is x_28730 at +1 equation.  Only 5 free inputs move a22231 at all, and all
16,806 integer combinations of them (coeffs -3..3, plus p-scaled pairs) give zero zero-collateral
directions and best failing 7.  a37887 lives in exactly ONE equation, eq 8680; of the 17 free
inputs that move it, only x_28730 moves nothing else outside the region, and it also moves a22231.
**No direction moves a22231 without moving a37887 -> a22231 buys 1 row and costs eq 8680, exactly.**

## 8. PER-PLACEMENT SWEEP OF THE CASCADE - the whole cascade is now priced (LOG.md Step 11)
20 pin atoms extracted by `chain.py`; each made the defect carrier by `sweep.py` and priced with
the five-stage pipeline.  Full table in LOG.md.  Headline:

    WITNESS  |R|=12 |S|= 8 deficit= 4 knobs=9 rank=7 zeroable=5 failing= 7  score=39026  rank > deficit
    best cascade pin (a26731) |R|=28 |S|=14 deficit=14 knobs=9 rank=6 zeroable=13 failing=15  score=39018
    11 pins tie at 39,018; the rest fall to 39,003-39,011.

**The witness placement is the ONLY placement in the cascade whose knob-image rank exceeds its
balance deficit.  By the criterion I derived, 7 is the floor for the whole cascade, not just for
one region.**  Cascade pins have only 2-3 zero-collateral knobs because a pin's free input is
consumed closing the pin itself; the witness sits at the one place where NINE free inputs act with
zero collateral.

## Single next experiment
The sweep covers carriers that are cascade PINS.  The one carrier class left unpriced is a defect
placed on a HANDLE rather than a pin - the p-quantised quotient variables (x_7497, x_11436,
x_22820, x_14393 and their ~1,240 siblings).  A handle carries the defect in multiples of p rather
than in a residue, so its region and its knob image are different objects from anything in the
table above.  Run `sweep.py` with the carrier set replaced by the handle atoms.
