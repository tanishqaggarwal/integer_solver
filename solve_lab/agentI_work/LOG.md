# Agent I — LOG

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

## 8. Open
`beam.py` (support search) is correct but `intsolve.solve_int` suffers HNF coefficient
blow-up on 300-bit entries. Needs: Bareiss rank filter over Q, then SNF integrality on
survivors, testing "can we fail <= 6" by enumerating the sacrificed set of size <= 6.
Then search for a compensator for **eq8680** whose own equations already lie in `E(S)`.

## 9. Framing note
An earlier version of this log carried a geometric reading of items 4-6. It has been
removed; nothing above depends on it. All statements are congruences, integer linear
relations, and rank/Smith-normal-form facts about the coefficient matrices.
