# Agent I — LOG  (FINAL)

## 0. The gap I flagged is closed (by agent F)
I headlined that "all atoms zero" is sufficient but not necessary, leaving every
optimality result in this lab conditional on one branch. Agent F computed
`rank(M) = 39,033`, `dim ker(M) = 0` for the equation-atom incidence matrix (peeling
certificate + Wiedemann, pivots non-divisible by any odd prime, so over Z). Any
assignment satisfying all 39,033 equations therefore makes ALL atoms exactly zero.
The model is an equivalence, not a restriction; everything below is unconditional.

Reconciliation: my decomposition has **40,885** atoms, so my own incidence matrix has
`dim ker >= 1,852` by arithmetic. Those directions are exactly the **926 self-cancelling
single-equation pairs** I measured in 7c; F's coarser M quotients them out. The results
agree — only the granularity differs. Match decompositions before comparing kernels.


All entries are integer / polynomial linear algebra over `EQUATIONS.txt`.

## 1. Baseline
Re-verified `best/new_instance_partial_39026.json` with `solve_lab/checker.py`:
39,026/39,033, failing `[12231,12270,12350,14584,18673,22044,29125]`. Confirmed.

## 2. Independent parse
`parse.py` -> 39,033 equations, **40,885 distinct atoms**; every atom is degree <= 2 with
<= 3 terms. Equation = `outer(CORE)`, CORE a left-nested integer combination of atoms.
`model.py` scores an assignment exactly in 0.1 s and reproduces `checker.py` exactly.

## 3. Z propagation from empty (`prop.py`)
Forces 5,624 variables, **zero conflicts**: 4,109 zeros, 1,286 ones, 220 variables equal
to `p = 2^256-2^32-977`, 9 equal to a second 256-bit constant `K`. Those are the only
large pinned constants in the instance.

## 4. Effective-support reduction — the main technical result
A variable whose coefficient vanishes under the current partial assignment is not an
unknown. 3,707 atoms have the shape `X - Y*H` with `Y` pinned to `p`; every handle `H`
occurs in exactly one atom, so mod `p` each is the unit clause `X == 0 (mod p)` and no
information is lost. Reducing atoms against the partial assignment BEFORE counting
unknowns took mod-p propagation from 13,490 to **28,701 of 38,748 variables from 1,156
boolean decisions in ~1 s**, and reduced the whole system to **3 violated atoms**
`X2287 - 8272701*X35389`, `X21889 - 8646263*X35389`, `X25156 - 10159099*X35389`
(rank 2 in `(X35389, X6671)` => both forced to 0 mod p).

Symbolic back-substitution of the instance's own atoms (exact identity check, sympy):
```
X35389 = (x2-x1)^2*(x3+x1+x2+K) - (y2-y1)^2
X6671  = (y3+y1)*(x2-x1) - (y2-y1)*(x1-x3)
x1=X12186 y1=X16742 x2=X14853 y2=X24908 x3=X22162 y3=X30213 K=X24453
```

## 5. Failed / abandoned lines (with reasons)
- **No-good loop over the choice of free inputs** (`fploop.py`, `loop2.py`): correct but
  60x too slow before the effective-support fix; superseded by `boolscore.py`.
- **Structural triangulation** (`orient*.py`): stalls at 24,376/38,748 because gates
  `X - Y*Z` with two unknowns cannot be oriented forward; not useful.
- **Single-flip scan of the 1,156 boolean decisions** (`flipscan.py`): 123 measured,
  **all inert** on `X35389` and `X6671`. Abandoned.
- **Blind cut screen over all 1,865 atoms in <= 6 equations** (`cutscan.py`): ~9 h under
  CPU contention; replaced by the targeted `cutlocal.py` (421 atoms on the defect path).

## 6. Cut-cost measurements
Complete conflict-free mod-p states exist with only **two** atoms released: free
`a40368 = X24453 - K` (which occurs in exactly **1** equation) and solve `X35389 = 0`
for `X24453`, then free one coordinate variable and solve `X6671 = 0` for it. Cheapest
total cost measured: **13 equations** (`{a40368, a29331}`) => 39,020.
`mincut.py`: no second release costs fewer than 12 equations. `cutlocal.py`: same floor.
=> the deliverable wins by **cancellation**, not by a smaller defect.

## 7. The cancellation experiment (`cancel.py`) — exact
Support = the 9 atoms nonzero at the deliverable + the free compensator
`a23434 = X4432 - X19964 - X28730`. Knobs = variables occurring only inside the support
(`X642 X1329 X9413 X10903 X17325 X28730 X29854 X31864`), so every point of the
parameterisation is realisable by construction. 13 equations touched; exact
Smith-normal-form enumeration over all subsets gives **max 6 simultaneously satisfiable
=> 7 failing => 39,026**, i.e. the compensator is worth exactly zero.

Why, precisely: `a23434` adds one equation **eq8680** whose row is `[0,...,-1,...,0]`
with base 0, pinning the new knob `d28730 = 0`. Without that pin `eq29125` is
satisfiable — its row is `[0,0,-p,0,0,-30,0,0]` and `gcd(p,30)=1`. **eq8680 is the single
row that costs the deliverable its seventh equation.**

## 7b. Support search — the full negative result
`search.py`. Only **26 atoms** in the instance are adjacent to the deliverable's 12
equations; cheapest new-equation costs 1, 2, 2, 3, 3, 4, ... Evaluated all 26 singles,
all 28 pairs among the 8 cheapest, all 20 triples among the 6 cheapest — 74 supports
with |E| from 12 to 21 and 7 to 10 knobs. **min failing = 7 for every one.** |E| grows
by exactly what the new knobs buy, every time.

`exact.py` upgrades greedy (an upper bound on minfail) to proof by enumerating every
sacrificed set of size <= 6. **Deliverable's own support (|E|=12, 7 knobs): PROVED
minfail >= 7** — all C(12,<=6) subsystems are integrally unsolvable.

Bug found and fixed on the way: `intsolve.col_hnf` had its column operation roles
swapped (`col_piv` updated instead of `col_dst`), so it never terminated. With it fixed,
my solver independently reproduces the optimum 7 AND the optimal sacrificed set
`[12231,12270,12350,14584,18673,22044,29125]` — the deliverable's own failing lines,
recovered from scratch by integer linear algebra.

## 7c. Probe of the off-branch structure — negative
`1,853` atoms occur in exactly one equation: **926 pairs plus one**. Each of 926
equations carries exactly two of them, always a product/difference atom `P` and a
bare-variable atom `X`. Hypothesis: `X` is free slack that cancels `P` at zero cost.
**Refuted by direct computation**: every one of those 926 variables occurs in exactly 3
atoms (the bare atom, a definition `X - Y*Z`, and the pairing atom `P + X`), so the
trade is local and self-cancelling. Equations with a genuinely free absorber (bare
variable, unit coefficient, variable in no other atom): **0**.

## 7d. Cascade (two-level) closure
`cascade.py`. A variable with one atom outside the support is usable if that atom is
held at zero by re-solving it for another variable whose atoms are all inside; such
atoms cost NO equations. Fixed-point closure from the deliverable's support absorbs
**1,817 atoms** as repairable, and the knob count stays at **exactly 7** — the closure
consumes one dependent variable for every knob it frees. minfail >= 7 re-proved
exhaustively with those knobs. The ripple / two-level move class buys nothing.

## 7e. The eq8680 hunt (`eq8680.py`, `hunt.py`, `cascade_rand.py`)
eq8680 is the one row that pins the knob `d28730` to zero and so costs the deliverable
its seventh equation.

**Complete candidate census.** To change eq8680's core, some atom IN eq8680 must move;
eq8680 has exactly 18 atoms; to move an atom you need a knob; a knob is a variable all
of whose atoms are in the support. So the candidate set is exactly `{v2a[x]}` over the
variables of eq8680's atoms: **43 groups, 30 with nonzero net effect on eq8680** (13
move it by exactly 0 and provably cannot compensate). This is an enumeration, not a
sample. Supporting fact: **only 28 atoms in the whole instance share an equation with
E(S)**; imported-equation counts 1,2,2,3,3,4,4,5,...

The cheapest compensator is `X19964`, group `{a1631, a23434}` — ONE new atom, net effect
exactly -1 on eq8680, the exact counterweight to `d28730`. My earlier adjacency search
missed it because `a1631`'s own equations do not touch E(S) at all.

**Every candidate tested gives minfail > 6.** The reason is uniform: every knob that
moves eq8680 also moves an atom living in 11-16 equations, all of which sit at base 0
and are pushed off zero by that same knob. Each candidate buys one row and pays 5-27.

**Cascade CHOICE is not a free parameter** (`cascade_rand.py`): randomising which
dependent variable each absorbed atom consumes gives **8/8 seeds identical** — 1,817
atoms, 1,817 dependents, the same 7 knobs. The closure is a confluent fixed point.

**Second bug found and fixed**: `minfail_bnb` ordered already-satisfied rows FIRST, so
the DFS spent its sacrifice budget on trivially keepable rows before reaching the hard
ones. Reversing it gave ~20x: 300 s timeouts became 2-36 s exact answers. That, not
instance hardness, is why the first run appeared to hang.

## 8. Final status
Best verified anywhere, re-verified by me with `solve_lab/checker.py`: **39,026 / 39,033**,
failing `[12231,12270,12350,14584,18673,22044,29125]`. I produced nothing above it and
wrote no new assignment.

Every realisable move class I could parameterise exactly — one-level knobs, cascade
knobs, supports extended by up to 3 adjacent atoms, and the complete knob-group census
for eq8680 — gives min failing = 7, with exhaustive proofs on six supports and on every
eq8680 candidate. Combined with `dim ker(M) = 0` (section 0), a full solution requires
all atoms exactly zero, which reduces to two explicit polynomials vanishing mod p under
a 256-bit boolean selector choice; nothing I measured gives a handle on that choice.

Not covered, stated plainly: atom vectors realisable only through simultaneous moves of
many circuit variables. My knob parameterisation is realisable by construction but
local; my global re-solves reach that space but land at cost >= 13 equations. A method
that is both global and exact over it is the only unexplored direction I can name.

## 9. Framing note
An earlier version of this log carried a geometric reading of items 4-6. It has been
removed; nothing above depends on it. All statements are congruences, integer linear
relations, and rank/Smith-normal-form facts about the coefficient matrices.
