# Agent D log — part 3

(line withdrawn per user instruction -- see FLEET.md)
(line withdrawn per user instruction -- see FLEET.md)
residual first). All numbers below are mine, computed in agentD_work only.

* A's closed form carries an offset K = x_24453 = C_A. That offset is the `a2` coefficient of a
(line withdrawn per user instruction -- see FLEET.md)
  x3 = lam^2 - a2 - x1 - x2, which is exactly my measured
(line withdrawn per user instruction -- see FLEET.md)
  because they assumed a2 = 0.
* Fitting a4, a6 from TWO table entries with a2 = C_A:
      a4 = 114170008767671698752186727197936107864370654164657728518655355473804451402762
      a6 = 77755683306591771556999954628254672912734268662742093169295805431582354953490
  **256/256 table point pairs fit, and the pinned target (x3,y3) = (C_C mod p, C_B mod p) fits.**
* Depressing x -> x + a2/3 gives the short form y^2 = x^3 + A x + B with **A = 0 exactly** and
      B = 64019533680030876408443198762210829058751700634554282185987325820393598524794
(line withdrawn per user instruction -- see FLEET.md)
(line withdrawn per user instruction -- see FLEET.md)
* Group order: of the six candidate orders from 4p = L^2 + 27 M^2, the annihilating one is
  n = 115792089237316195423570985008687907852837564279074904382605163141518161494337
(line withdrawn per user instruction -- see FLEET.md)
(line withdrawn per user instruction -- see FLEET.md)
  255 of 256 points P have 2P also in the table, and following the chain from the unique root
  visits all 256.  So table = { 2^i * P_0 }.
(line withdrawn per user instruction -- see FLEET.md)
(line withdrawn per user instruction -- see FLEET.md)
(line withdrawn per user instruction -- see FLEET.md)

### Did my pipeline ever move (x3,y3)?  YES — and it does not open the instance.
x_22162 (x3) and x_30213 (y3) are FREE INPUTS (ancestor cone size 1). `ecsolve2.py` moved them by
arbitrary residues mod p
(dx3 = 28179459175168346424568155963891666812586454418931334303510798603201179609960,
 dy3 = 84170865372055296292320112314180874192322915282344166253532298893879157282939)
and reached a checker-verified 39,017 with A = B = 0 exactly mod p.  But (x3,y3) is the CLAIMED
output, asserted equal to the public target by a688 (y3 == C_B) and a1618 (x3 == C_C); moving it
(line withdrawn per user instruction -- see FLEET.md)
So the only real knobs are the 256 selector bits — quantitatively, `scanAB.py` perturbed ALL
7,273 free inputs with full advice re-solve and found **zero** cost-free movers of (A,B), and only
**10 distinct (dA:dB) directions** exist, the cheapest costing 10 equations.

## 39,026 is EXACTLY optimal for its placement (independent re-derivation)
`gens26b.py` scans all 38,748 variables with the witness's five deliberately-broken GATE atoms
BLOCKED (without blocking, a plain ripple silently repairs them — this is the trap that makes
naive search collapse from 39,026 to ~39,008).  Cost-free generator lattice = 9 generators:
  x_642 (a22229 in steps of -7376877, a35762 +1), x_1329, x_8731 (a35761 in steps of 1),
  x_9118 (steps of 5113045), x_9413, x_10903, x_17325, x_29854, x_31864.
Only a22231 is completely frozen.  `lat26b.py` + `intsolve.py` (column-HNF integer solver) then
enumerates ALL 2^12 subsets of the 12 equations and solves each exactly over Z:
sizes 12..6 are ALL integrally infeasible; size 5 is solvable ({2554,6816,8124,9123,9421}, k=0).
**MAX satisfiable = 5, failing = 7, score = 39,026.**  Prior claim CONFIRMED, with a strictly
larger generator set than the previous derivation used.

## Alternative placements I priced exactly
| residual placement | |E| | free (out=0) atoms | measured/derived failing | score |
|---|---|---|---|---|
| the deliverable's 7 atoms | 12 | 8 (9 generators) | 7 (proved) | 39,026 |
| {a688, a1618, a40608} = the (x3,y3) pins | 16 | 4 | 16 (measured) | 39,017 |
| {a31670, a31672} = bit x_24601's pins | 12 | 4 | ~10 (n - c = 2 satisfiable) | ~39,023 |
| {a3576, a3578} = bit x_2081's pins | 15 | - | - | worse |
| D_adv's EC-identity atoms | 20 | 9 | ~13-17 | <= 39,020 |
Minimum union over all independent knob PAIRS that span (dA,dB): 12 equations (x_22152+x_33462).
