# Agent I — RESUME

Everything here is stated as integer / polynomial linear algebra over `EQUATIONS.txt`.

## Score status
- Deliverable re-verified by me with `solve_lab/checker.py`: **39,026 / 39,033**,
  failing `[12231, 12270, 12350, 14584, 18673, 22044, 29125]`. CONFIRMED.
- Nothing of mine beats it yet. No new best assignment written.

## My independent model (do not redo)
- `parse.py` -> `atoms.pkl` : 39,033 equations, **40,885 distinct atoms**. Each equation
  is `outer(CORE)` with `outer` one of `c*`, `P*P`, `P^3`, `c1*P + c2*P`, and
  `CORE = a0 + c1*a1 + c2*a2 + ...`, a left-nested integer combination of atoms.
- `poly.py` -> `polys.pkl` : every atom is a polynomial of **degree <= 2 with <= 3 terms**.
- `model.py <assign.json>` : exact integer scorer, **0.1 s** (vs 16 s for checker.py),
  reproduces `checker.py` line for line.
- `prop.py` : exact propagation over Z from the empty assignment. Forces **5,624
  variables with zero conflicts**: 4,109 zeros, 1,286 ones, 220 variables equal to
  `p = 2^256 - 2^32 - 977`, 9 equal to
  `K = 97553848499418123410591666447050222001188385549510401465815187079080512838891`.
  Those are the only large constants pinned anywhere in the instance.

## The best technical idea I produced: EFFECTIVE-SUPPORT REDUCTION
Pure linear algebra, no framing needed:

> **A variable whose coefficient vanishes under the current partial assignment is not
> an unknown.** Reduce each atom against the partial assignment FIRST, drop the
> monomials whose coefficient became 0, and only then count unknowns.

`3,707` atoms have the shape `X - Y*H` where `Y` is one of the 220 variables pinned to
`p`. Work mod `p` and the coefficient of `H` dies, so the atom is the **unit clause
`X == 0 (mod p)`**, not a 2-unknown constraint. Every handle `H` occurs in exactly one
atom (verified), so those 3,707 atoms are exactly 3,707 independent congruences with a
free quotient each — no information is lost by reducing mod p.

Measured effect (`fp.py` + `boolscore.py`): propagation goes from 13,490 determined
variables to **28,701 of 38,748, from 1,156 boolean decisions, in ~1 second**, and the
whole 39,033-equation system collapses to **exactly 3 violated atoms**:
```
a17810: X2287  -  8272701 * X35389
a17813: X21889 -  8646263 * X35389
a17816: X25156 - 10159099 * X35389
```
Three different multipliers, and each `X_k` is an independently derived multiple of
`X6671`, so the 3 rows have rank 2 in `(X35389, X6671)` and force **X35389 = X6671 = 0
mod p**. By exact symbolic back-substitution of the instance's own atoms (sympy, an
identity check, not numeric agreement) these two are the polynomials
```
X35389 = (x2-x1)^2 * (x3 + x1 + x2 + K) - (y2-y1)^2
X6671  = (y3+y1)*(x2-x1) - (y2-y1)*(x1-x3)
x1=X12186  y1=X16742  x2=X14853  y2=X24908  x3=X22162  y3=X30213  K=X24453
```
These seven variables are pinned through 256 boolean selector bits (each selector, when
on, pins two variables to explicit 296-bit constants through atoms
`b*(x - C) - m*h`). All 512 such handles `h` are forced to 0 mod p (verified), so an
active selector fixes its value exactly mod p. Single-bit flips: 123 measured, **all
inert** on `X35389` and `X6671`.

## Cut-cost measurements (mine, reproducible)
Releasing atoms (allowing them nonzero) and re-solving mod p, complete conflict-free
mod-p states exist with only **two** atoms released — free `a40368` (`X24453 - K`,
which occurs in exactly **1** equation) and solve `X35389 = 0` for `X24453`, then free
one of the seven coordinate variables and solve `X6671 = 0` for it:
```
{a40368, a29331}  13 equations   {a40368, a10066}  14
{a40368, a10067}  15             {a40368, a26748}  15
{a40368, a10065}  15             {a40368, a29334}  15   {a40368, a29333} 15
```
`mincut.py` enumerated every atom in the derivation chain of the seven variables: **no
second release costs fewer than 12 equations**, so the cheapest mod-p cut costs 13
equations => 39,020. `cutlocal.py` (421 atoms on the defect path) gives the same floor.

**The deliverable beats all of these not by having a smaller defect but by CANCELLATION:
its 7 nonzero atoms span 12 equations of which 5 cancel.**

## THE LIVE EXPERIMENT — min |{e : sum_a c_{e,a} v_a != 0}|
`cancel.py`, `search.py`, `exact.py`, `cascade.py` (all done). Parameterisation is realisable **by
construction**: a *knob* is a variable occurring only in atoms of the chosen support, so
perturbing it leaves every atom outside the support exactly zero and every equation
outside `E(S)` exactly satisfied. For a support `S`:
`core_e(d) = base_e + sum_k Mat[e][k]*d_k`, exact integers, linear in the knobs
(linearity re-checked numerically for every knob).

**Result so far.** Support = the 9 atoms nonzero at the deliverable
`{a23432, a23433, a36225, a36226, a36227, a36228, a36229, a37537, a37538}` plus the free
compensator `a23434 = X4432 - X19964 - X28730`:
```
knobs (occur only inside the support): X642 X1329 X9413 X10903 X17325 X28730 X29854 X31864
equations touched: 13  [2554,6816,8124,8680,9123,9421,12231,12270,12350,14584,18673,22044,29125]
MAX simultaneously satisfiable (exact, Smith normal form over ALL 2^13 subsets): 6
=> failing 7  => SCORE 39,026
```
So the compensator is worth exactly nothing, and I can now say precisely why:
`a23434` adds one new equation, **eq8680**, whose row is `[0,0,0,0,0,-1,0,0]` with base 0
— it pins the new knob `d28730` to zero. Without that pin, `eq29125` would be
satisfiable: its row is `[0,0,-p,0,0,-30,0,0]` and `gcd(p,30) = 1`, so
`p*d9413 + 30*d28730 = base` is solvable for ANY base. **eq8680 is the single row that
costs the deliverable its seventh equation.** Its other 17 atoms all occur in 11-16
equations each, which is why a naive compensator for it does not pay.

## Support search — the full result (`search.py` -> `search.log`, `search_result.json`)
Only **26 atoms** in the whole instance are adjacent to `E(S)` (share an equation with
the deliverable's residual). New-equation cost of the cheapest ones: 1, 2, 2, 3, 3, 4, ...
Evaluated **all 26 singles, all 28 pairs among the 8 cheapest, all 20 triples among the
6 cheapest**:

> **min failing = 7 for every single one of them.**

|E| grows by exactly what the new knobs buy, every time — over 74 distinct supports with
|E| from 12 to 21 and knob counts from 7 to 10. `exact.py` upgrades this from greedy
(an upper bound on minfail) to an exhaustive proof for the tightest supports by
enumerating every sacrificed set of size <= 6.

Independently re-derived with my own fixed HNF solver (a genuine bug found and fixed in
`intsolve.col_hnf`: the column operation had its roles swapped, so it looped forever):
for the deliverable's own support the exact minimum is **7**, and the optimal sacrificed
set is exactly `[12231, 12270, 12350, 14584, 18673, 22044, 29125]` — the deliverable's
own failing lines, recovered from scratch by integer linear algebra.

## THE eq8680 HUNT — the campaign's remaining question, run to ground
`eq8680.py` (census + exact branch-and-bound), `hunt.py` (prioritised), `cascade_rand.py`.

**Why the candidate list is complete.** To change eq8680's core, some atom *in eq8680*
must change value. eq8680 has exactly 18 atoms. To change an atom's value, a knob must
move it; a knob is a variable all of whose atoms are in the support. So the complete
candidate set is `{v2a[x] : x a variable of an eq8680 atom}` — **43 groups**, of which
**30 have nonzero net effect on eq8680's core** (the other 13 move it by exactly 0 and
provably cannot compensate). That census is in `eq8680.log`; it is an enumeration, not
a sample.

Supporting census: **only 28 atoms in the entire instance share an equation with E(S)**,
and their imported-equation counts run 1, 2, 2, 3, 3, 4, 4, 5, ... The cheapest
compensator of all — `X19964`, whose group is just `{a1631, a23434}`, ONE new atom, net
effect exactly −1 on eq8680, the perfect counterweight to `d28730` — was missed by my
earlier adjacency search because `a1631`'s own equations do not touch E(S) at all.

**Result: ALL 22 single knob-groups give minfail > 6 — zero hits, zero timeouts.**
(66 PAIRS of groups running in the background -> `hunt.log`; every one completed so far
is also > 6. Re-check with `grep -E 'BEATS|minfail = ' hunt.log` — an empty result means
no pair beats 39,026 either.)
```
X19964 add=(1631,)              |E|=27 (+14)  knobs=9   minfail > 6
X4432  add=(2427,22331,33706)   |E|=40 (+27)  knobs=9   minfail > 6
X6947  add=(23435,)             |E|=14 (+ 1)  knobs=9   minfail > 6
X33168 add=(23437,)             |E|=16 (+ 3)  knobs=9   minfail > 6
X10422 add=(11772,)             |E|=22 (+ 9)  knobs=9   minfail > 6
X11099 add=(11774,)             |E|=24 (+11)  knobs=9   minfail > 6
X22526 add=(11776,)             |E|=25 (+12)  knobs=9   minfail > 6
X34868 add=(11778,)             |E|=27 (+14)  knobs=9   minfail > 6
X950   add=(20290,)             |E|=18 (+ 5)  knobs=9   minfail > 6
X15120 add=(20292,)             |E|=20 (+ 7)  knobs=9   minfail > 6
X35531 add=(20294,)             |E|=21 (+ 8)  knobs=9   minfail > 6
X18253 add=(20292,20293)        |E|=21 (+ 8)  knobs=10  minfail > 6
X37720 add=(20294,20295)        |E|=22 (+ 9)  knobs=10  minfail > 6
X23822 add=(11776,11777)        |E|=26 (+13)  knobs=10  minfail > 6
X7945  add=(11778,11779)        |E|=28 (+15)  knobs=10  minfail > 6
X9629  add=(20290,20291)        |E|=19 (+ 6)  knobs=10  minfail > 6
X37413 add=(11774,11775)        |E|=25 (+12)  knobs=10  minfail > 6
X37254 add=(11775,11917)        |E|=42 (+29)  knobs=9   minfail > 6
X35619 add=(23437,23438,35830)  |E|=17 (+ 4)  knobs=10  minfail > 6
X23642 / X30108 / X15324 / X4432 ...  all minfail > 6
```
The reason is uniform and visible in the numbers: every knob that moves eq8680 moves an
atom that lives in 11-16 equations, and the imported equations all sit at base 0 and are
moved off zero by that same knob. **Each candidate buys one row and pays for 5-27.**

**Coordinator item 1 — the cascade CHOICE is not a free parameter.** `cascade_rand.py`
randomises which dependent variable each absorbed atom consumes: **8/8 seeds give the
identical outcome** — 1,817 atoms absorbed, 1,817 dependents, and exactly the same 7
knobs `[642, 1329, 9413, 10903, 17325, 29854, 31864]`. The closure is a confluent fixed
point; the choice was never the limit.

**Second performance bug found and fixed**: `minfail_bnb` sorted already-satisfied rows
FIRST, so the DFS burned its whole sacrifice budget on trivially keepable rows before
reaching the hard ones. Reversing the order made it ~20x faster (300 s timeouts became
2-36 s exact answers). This is why the first run appeared to hang, not a hard instance.

## Cascade (two-level) closure — the move class a one-level analysis cannot see
`cascade.py`. A variable with ONE atom outside the support is still usable if that atom
can be held at zero by re-solving it for another variable whose own atoms are inside.
Such atoms cost **no equations**. Closing that rule to a fixed point from the
deliverable's support absorbs **1,817 atoms** as repairable — and the knob count stays
at **exactly 7**. The closure consumes one variable as a dependent for every one it
frees. `minfail >= 7` re-proved exhaustively with the cascade knobs. So the ripple /
two-level move class buys nothing here.

## Probe of the "off-branch" structure — negative, and worth recording
`1,853` atoms occur in exactly one equation. They are **926 pairs plus one**: every one
of 926 equations carries exactly two of them, always a product/difference atom `P` and a
bare-variable atom `X`. That looked like free slack (`X` set to cancel `P` at zero
cost). It is not: **every one of those 926 variables occurs in exactly 3 atoms** — the
bare atom, a definition `X - Y*Z`, and the pairing atom `P + X` — so the trade is local
and self-cancelling. Number of equations with a genuinely free absorber (bare variable,
unit coefficient, variable in no other atom): **0**.

## Reproduce
```
cd /home/user/integer_solver/solve_lab/agentI_work
python3 parse.py && python3 poly.py && python3 dag.py   # caches, ~1 min
python3 model.py ../best/new_instance_partial_39026.json
python3 prop.py           # Z propagation from empty
python3 boolscore.py wit  # 1 s mod-p re-solve, prints the 3 violated atoms
python3 cancel.py         # the exact cancellation optimum for the deliverable's support
python3 beam.py           # support search (needs the faster solvability test, see below)
```

## THE GAP — the reason this experiment matters
"All atoms zero" is **sufficient** for all 39,033 equations but **not necessary**: every
equation is an integer combination of 3-24 atoms, **1,853 atoms occur in exactly one
equation**, and the deliverable itself carries 9 nonzero atoms. Every optimality
argument this campaign has produced — mine included — is computed inside the
all-atoms-zero branch. The cancellation formulation above is the probe of what is
outside it.

## Single highest-value next step
`beam.py` is correct but too slow: `intsolve.solve_int` (column HNF) suffers coefficient
blow-up on the 300-bit entries. Replace the solvability test with (i) a fraction-free
Bareiss rank test over Q as a fast filter, then (ii) the Smith-normal-form integrality
test only on the survivors, and test directly for "can we fail <= 6" by enumerating the
**sacrificed** set of size <= 6 rather than searching from the top. Then run the support
search over all atoms adjacent to `E(S)`, and in particular look for a compensator for
**eq8680** whose own equations are already inside `E(S)`.
