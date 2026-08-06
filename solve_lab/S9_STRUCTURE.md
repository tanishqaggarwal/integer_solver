# Session 9 — the obstruction fully characterised, and the deliverable improved to 39,024

Verified deliverable: **39,024 / 39,033** — `best/new_instance_partial_39024.json`
(`python3 checker.py best/new_instance_partial_39024.json`, and independently
`python3 s9/verify_ast.py best/new_instance_partial_39024.json`).
Failing lines: `[9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]`.
It differs from the previous 39,022 partial in **five variables**: `x_642, x_9413, x_28730,
x_29854, x_31864`. Sections 1-5 and 7-10 below describe the 39,022 branch and remain accurate;
**section 6 records a retracted claim and the idea that superseded it — read it first.**

This session re-derived the instance from scratch with an independent toolchain (`s9/`), corrected a
methodological error that had been distorting every previous Jacobian-based attack, and **proved**
that the current partial cannot be improved without cracking the core.

## 1. Independent model (all numbers re-derived, not inherited)

`s9/atomize.py` decomposes every equation into `(coeff, atom)` pairs and **validates by exact
integer re-evaluation against all 39,033 raw equations — 0 mismatches**.

| quantity | value |
|---|---|
| distinct atoms | 42,267 (deg 1: 19,780 · deg 2: 21,788 · deg 4: 699) |
| gates (canonical orientation) | 31,475 |
| free inputs | 7,273 (1,156 of them boolean) |
| check atoms | 10,792 (783 are perfect squares) |
| nonzero atoms at the best partial | **3**: 22229, 22231, 37887 |

The canonical gate orientation is recovered from the setter's own syntax (`s9/gates.py`): the output
of a gate is the first *bare-variable* top-level additive term. This reproduces the intended data
flow and makes the forward ripple exact (0 non-divisible gates over all 7,273 unit perturbations).

## 2. The obstruction is exactly two modular congruences

```
atom 22229 :  A = x_7068 − x_2099 − 7376877·x_642        (x_642   = x_28599·x_17325)
atom 22231 :  B = x_4432 − x_19964 − x_28730             (x_28730 = x_17499·x_9413)
atom 37887 :  square whose root = (22231) + gate residuals  → 1 extra equation
```

`x_28599 = x_17499 = x_26064 = p = 2^256 − 2^32 − 977` **exactly** (a fixed 256-bit prime), so
`x_642 = p·x_17325` and `x_28730 = p·x_9413` with `x_17325`, `x_9413` free quotient handles. Hence

> **the whole residual system is `x_7068 ≡ x_2099` and `x_4432 ≡ x_19964`, both mod p.**

`x_2099` and `x_19964` are outputs of a 2-bit MUX (controls `x_4287 = 0`, `x_2081 = 1`):

```
x_2099  = b1(1−b2)·x_31861 + b2(1−b1)·x_6418  + b1b2·x_9118   →  x_6418
x_19964 = b1(1−b2)·x_14865 + b2(1−b1)·x_12553 + b1b2·x_8731   →  x_12553
```

and `x_6418`, `x_12553` are **pinned by load gates** to setter constants (atoms 3576 / 3578):

```
x_2081·(x_6418  − K1) = 15804267·x_26777 = 15804267·p·x_3387
x_2081·(x_12553 − K2) =              x_13458 =             p·x_5081
K1 = 33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319
K2 = 42775533402728869434716629464193396056515231264222641773817154079369026410240838606908039
```

So the conditions are **`x_7068 ≡ K1 (mod p)`** and **`x_4432 ≡ K2 (mod p)`**, and at the partial
`D1 = x_7068 − K1`, `D2 = x_4432 − K2` have `D1 mod p ≠ 0`, `D2 mod p ≠ 0`.
(Note `gcd(15804267, 7376877·p) = 3` and `3 | D1`, so the *pin* side is Bézout-solvable — the only
real barrier is the `p`-divisibility.)

## 3. Where the defect goes if you push it (the conserved obstruction, mapped)

Forward-rippling the canonical repair (`x_7068 := K1`, `x_4432 := K2`) touches only **115 variables**
and relocates the defect to two mirror congruences (atoms 29539 / 7930):

```
12846437·(x_14853 − x_1308)  = p·x_30163     with  x_1308  ≡ x_7068  (identically equal)
 9367949·(x_24548 − x_25442) = p·x_11052     with  x_25442 tracking x_4432
```

Repairing *those* with the free leaves `x_14853`, `x_24548` breaks the three "monsters"
(atoms 19297/19299/30984), which are the core:

```
M1 = x_15298·L1 + x_4007          x_4007  = p·x_30317      L1 = x_11150 = 8646263·S + 1073965·T
M2 = x_15298·L2 − 6672769·x_29804 x_29804 = p·x_5146       L2 = x_25739 = 10159099·S + 6926539·T
M3 = 537773·x_15298·L3 − x_35605  x_35605 = p·x_2936       L3 = x_37758 = 8272701·S + 5921311·T
```

with `x_15298 = x_7715·x_34554 = 1` and

```
S = x_35389 = A·u² − w²      T = x_6671 = B·u − w·c
u = x_29322 = x_14853 − x_12186        w = x_3558  = x_24908 − x_16742
A = x_33469 = x_22162 + x_12186 + x_14853 + x_24453
c = x_1326  = x_12186 − x_22162        B = x_27713 = x_30213 + x_16742
```

Core ⟺ `p|L1, p|L3, 6672769p|L2` ⟺ **`S ≡ T ≡ 0 (mod p)`** plus one condition mod 6672769.
Eliminating `w` from `T ≡ 0` and substituting into `S ≡ 0` gives `u²·(A·c² − B²) ≡ 0`, i.e.

> **either `u ≡ w ≡ 0`  (the branch the partial sits on)  or `A·c² ≡ B² (mod p)` — a second branch
> not recorded in any previous session.**

The second branch is blocked in practice because `x_22162`, `x_30213` are pinned mod p by atoms
1618 / 688 to the setter constants `91416258…101002` and `125787314…635626`, and `x_16742` is
pinned to `x_19083` by atom 26731.

## 4. Methodological correction (this invalidates earlier Jacobian conclusions)

**783 of the 10,792 check atoms are perfect squares `E²`.** At the partial `E = 0`, so a
finite-difference row for `E²` is `c²·δ²` — *quadratic*, not linear. Any linear/Newton/null-space
analysis that feeds the square itself into the matrix is solving the wrong system.

`s9/roots.py` extracts the degree-2 root of all 783 (699 deg-2, 84 deg-1; none left over) and
`s9/jac2.py` rebuilds the Jacobian on roots. Effect: the reported inconsistency moved from the
*spurious* atom 42245 (a square) to the *genuine* core atom 19297. Prior sessions' "inconsistent
mod-p Jacobian" results were reading the artefact, not the obstruction.

With the corrected model the mod-p system over all 7,273 free inputs is still inconsistent, both at
`S0` (u = 0, where the core gradient vanishes identically) and at `S1` (u ≠ 0, non-degenerate) —
certificate spans 270 rows. So the barrier survives the correction; it is just now correctly located.

## 5. Bit search

Of the **1,156 boolean free inputs** (only 2 currently set), exactly **two** deactivate the core by
forcing `x_15298 = x_7715·x_34554 = 0`: `x_2081` and `x_24601` — the two known "quadrant"
activators. Both are strictly worse locally (17 and 19 nonzero residuals vs 4). No single flip
reduces the residual count below the current 4.

## 6. ~~39,022 is provably optimal~~ — **RETRACTED, AND SUPERSEDED: the answer is 39,024**

> **This section originally claimed 39,022 was a local optimum. That claim was WRONG and the
> deliverable is now `best/new_instance_partial_39024.json` = 39,024/39,033 (9 failing).**
> Verified twice, by `checker.py` and by an independent AST-walk verifier (`s9/verify_ast.py`,
> no `eval`, no `compile`, no regex substitution). Failing lines:
> `[9123, 9421, 12231, 12270, 12350, 14584, 18673, 22044, 29125]`.

### The error

The old argument assumed the reachable defects were `A ∈ D1 + pℤ` and `B ∈ D2 + pℤ` "and nothing
else", because `x_28730 = x_17499·x_9413 = p·x_9413` is p-quantised. **That silently assumed atom
22230 (`x_28730 − x_17499·x_9413`) must stay zero.** It does not have to.

An equation is zero iff its *linear combination of atoms* vanishes — not iff every atom vanishes.
So **any atom whose entire equation footprint already lies inside the failing set is a free knob.**
Atom 22230 is exactly that: `x_9413` and `x_28730` occur in **no other atom of the instance**, and
22230 appears only in 10 of the 11 already-failing equations. Breaking it costs nothing and makes
`x_28730` — hence `B` — a *free integer*, not a member of `D2 + pℤ`.

Four atoms have zero footprint outside the failing set:

| atom | expression | equations | clean equations touched |
|---|---|---|---|
| 22229 | `x_7068 − x_2099 − 7376877·x_642` | 9 | 0 |
| 22230 | `x_28730 − x_17499·x_9413` | 10 | **0** ← missed originally |
| 22231 | `x_4432 − x_19964 − x_28730` | 10 | 0 |
| 37887 | square, root = 22231 | 1 | 0 |

### The construction that reaches 39,024

Extending the knob set along the ladder `35754…35762` (the gate outputs `x_29854`, `x_31864`,
`x_642`, whose whole equation footprint is the 11 failing equations plus only `{9123, 18673}`)
gives 5 knobs over 13 affected equations, of which **4 are simultaneously zeroable** over ℤ:
`11 + 2 − 4 = 9` failing. The witness differs from the old partial in exactly **five** variables —
`x_642, x_9413, x_28730, x_29854, x_31864` — with no bit flips, no ripple and the core untouched.

### What survives of the old argument

Only the narrow part: *if* one insists every atom outside the defect set vanishes, then the nine
linear failing equations need `c₁D1 + c₂D2 ≡ 0 (mod p)` and none qualifies. That is still true —
it is simply not the right constraint. The lesson is the general one: **work in equation space,
not atom space.**

### Is 9 the wall?

Exhaustively, within the region `S = FAILS ∪ {9123, 18673}` and all 8 variables whose footprint
lies inside `S` (`642, 1329, 9413, 10903, 17325, 28730, 29854, 31864`): **0 of 1287 5-subsets and
0 of 1716 6-subsets are integer-solvable**, so at most 4 equations can be zeroed there and 9 is
optimal *for that region*. The achievable atom-value lattice carries exactly four mod-p
congruences against 8 atom degrees of freedom; the binding resource is the number of mod-p-acting
knobs whose footprint stays inside `S`, which is exactly 4. Two independent beam searches (primal
over knob sets to depth 8, dual over equation sets to |S| ≤ 19) both plateau at 9.
**Outside that region 9-optimality is evidence, not proof.**

Separately and exhaustively measured: over all 7,273 free inputs, exactly **6** move `(D1, D2)`
mod p (`x_2081, x_4287, x_6418, x_9413, x_12553, x_17325`) and every one breaks other checks; and
all 254 single bit flips score 32–41 failing equations (minimum 32). So the bit/pin route does not
reach the defect — but that is proved only for single bits, not arbitrary subsets.

## 7. Tooling added (`s9/`, all regenerable)

`harness.py` (raw eq I/O) · `atomize.py` (validated atom decomposition) · `poly.py` (sparse polynomial
expansion) · `gates.py` (canonical orientation from syntax) · `fwd.py` (topological sort) ·
`ripple.py` / `repair.py` / `drive.py` (exact forward ripple + repair loop) · `cone.py` (backward cone
printer) · `sens.py` / `jac.py` / `jac2.py` / `newton.py` (Jacobians, root-corrected) ·
`roots.py` (square-root extraction) · `modsolve.py` / `cert.py` (sparse mod-p elimination +
inconsistency certificates) · `bits.py` / `bitscan.py` (boolean free inputs, flip scan) ·
`movevar.py`, `state0.py`, `jcheck.py` (diagnostics).

Rebuild caches: `cd s9 && python3 atomize.py && python3 poly.py && python3 gates.py && python3 fwd.py`.

## 8. Highest-EV next experiments

1. **Second core branch.** `A·c² ≡ B² (mod p)` frees `u` entirely. It needs `x_16742` (≡ `x_19083`)
   and `x_24908` off their pinned residues. Determine whether `x_19083 = x_6361 + x_23758` and
   `x_24908` are genuinely rigid mod p, or reachable — that is the one structural door this session
   opened and did not close.
2. **Quadrant re-solve.** `x_2081 → 0` or `x_24601 → 0` kills the core outright (17/19 residuals).
   A full forward re-construction on that branch was never completed with the corrected
   (root-based) residual model.
3. **Two-bit flips.** The single-flip scan is cheap (~1 s for all 1,156); the ~667k pairs are a few
   hours and were never run.

---

## 9. The activation branch: both cores satisfied simultaneously (new, and the deepest point reached)

The core is blocked because `u = x_29322 = x_14853 − x_12186` must be `≡ 0 (mod p)` while `x_14853`
is pinned to `K1` and `x_12186` is not. **`x_12186` can be freed.** Its definition chain is

```
x_12186 = x_23927 + x_25758,  x_25758 = x_10603·x_33612 = p·x_33612      (x_10603 = p)
x_23927 = x_7429 + x_26835,   x_26835 = x_38170·x_5096,  x_38170 = x_8599·x_21839
x_7429  = x_26865 + x_14192,  x_14192 = x_2754·x_30454,  x_2754  = x_21839·x_4549,  x_4549 = 1−x_8599
```

With `x_8599 = 0` (current) the MUX routes `x_12186 ← x_30454`, a computed value. **Set `x_8599 = 1`**
(keeping `x_21839 = 1`) and it routes `x_12186 ← x_5096` instead — and `x_5096` is a *free input*.
88 of the 1,156 boolean free inputs do this (`s9/find8599.py`, `hits8599.pkl`).

Then setting `x_5096 := K1`, `x_14853 := x_12186`, and letting the canonical gates close C1/C2 gives

> **u = w = S = T = L1 = 0 exactly — the original core (atoms 19297/19299/30984) is satisfied
> simultaneously with C1 and C2 for the first time in this instance.**

Residual after that step: 9 atoms, none of them the core (`s9/construct.py`, `construct2.py`).

### The obstruction regenerates one level up

Activating `x_8599` sets `x_38170 = 1`, which lights three atoms of *exactly the same shape*:

```
26733 : x_38170·x_21202 + x_11831          x_11831 = x_30022·x_5669
28438 : x_38170·x_15286 − 2264251·x_9216   x_9216  = x_17952·x_14485
32342 : x_38170·x_32453 − x_23535          x_23535 = x_38571·x_18963
```

and the loads decompose identically to L1/L2/L3:

```
x_21202 = 11598153·S' + 16335423·T'     x_32453 = 4677103·S' + 15469317·T'
S' = x_25614 = A'·u'² − w'²             T' = x_34220 = B'·u' − w'·c'
u' = x_18123 = x_30454 − x_10261        w' = x_17576 = x_16787 − x_25199
A' = x_32629 = x_21344 + x_24453        B' = x_16088 = x_21589 + x_25199   c' = x_33852 = x_10261 − x_5096
```

**This second core is zeroable.** A class of ~90 boolean free inputs shifts `u'` and `w'` by deltas
that are *exact complements mod p*:

```
u' = 82007976112976807461901870199198737303514020147647909878034348606308756230357
δ  = 33784113124339387961669114809489170549755964517992654161423235401600078441306   u' + δ = p  exactly
w' = 37841415183514949237467304684128824427406379377151921996714091976892367869714
δ' = 77950674053801246186103680324559083425863605288488642042743492031016466801949   w' + δ' = p exactly
```

So a two-bit construction (`s9/construct3.py`, 1,232 pairs tried) reaches states where **both cores
are clean simultaneously** — e.g. `(x_2527, x_1502)`: residual = 11 atoms, *all of them activated
load pins* `bit·(x_B − HUGE) = s·x_C` plus the mirrors 21617/26731/37662. No core atom remains.

### Why it still does not close

Each activated bit lights its own load pins, which pin further free inputs to setter constants.
Satisfying those pins (`s9/pinclose.py`) moves variables that lie in the *first* core's cone, so the
core re-lights: the closure loop runs 11 → 16 → 13 and stalls with no pin-style repair available
(63 failing equations). The obstruction is conserved through the activation route as well — it is
relocated, not removed.

**This is nevertheless the deepest structural point reached on this instance:** the core is not an
inviolable wall (it can be zeroed, twice), and the real invariant is the *pin/mirror cascade* that
each activation triggers. That cascade — not the core — is the object a future attack should target.

## 10. Score accounting for the new states

| state | residual atoms | failing equations |
|---|---|---|
| **best partial (unchanged deliverable)** | 3 | **11** |
| branch activation, one bit (`x_2527`) | 9 | 65 |
| two-bit, both cores clean (`x_2527`,`x_1502`) | 11 | 70 |
| after pin-closure iteration | 13 | 63 |

None beats 11 — consistent with the local-optimality proof in §6, which applies to *any* state whose
defect reduces to the same two congruences.

---

## 11. Parallel-agent round: the 39,024 construction, the other quadrant, and what is left

Three agents ran concurrently against the 39,022 baseline.

### 11.1 The win (39,024) — equation-space knobs
Described in §6 above. Catalogue built along the way: the instance has **512 load pins on 256
gate-bits** (2 per bit) — these bits *are* the 256-bit message. Every pinned input `x_B` and every
handle `h` is free and distinct, so **the pin conflict graph is empty**: setting all 256 bits to 1
and closing all 512 pins leaves 0 of 512 pin atoms violated. The stall in `pinclose.py` was never a
pin-variable collision; the cost is entirely downstream. Full single-bit sweep: all 254 off-bits
score **32–41 failing equations (min 32)** end-to-end, so the bit route does not approach 9.

Exhaustively measured: over all 7,273 free inputs, exactly **6** move `(D1, D2)` mod p —
`x_2081, x_4287, x_6418, x_9413, x_12553, x_17325` — and every one breaks other check atoms.
`x_9413` and `x_17325` are precisely the two that the 39,024 construction exploits, because their
atoms' footprints lie inside the failing set.

### 11.2 The other quadrant — measured, and one equation worse
Full re-construction on `x_2081 = 0` (core deactivated) reaches **39,021 / 39,033 (12 failing)**,
`s9/quad/best_quadrant_39021.json`. `x_24601 = 0` reaches 27 failing; both bits 0 reaches 43 (and
additionally violates the forced OR-gate `x_9274`). The branch is structurally *cleaner* — 2
residual atoms instead of 3, and the setter pins 688/1618 are reachable there, unlike C1/C2 here —
but the defect terminates on atoms 31670/31672 whose residues are rigid, and the reachable
`(m₁₉, m₃₀) ∈ [0,168] × [0,89]` bit box contains **no** configuration making either residue 0 or
cancelling any of the 12 equations. Notable side result: on that branch the `x_8599` reroute closes
**both setter pins 688 and 1618 outright** — a state never previously reached — but the second core
costs ~26 equations, so the best there is 30.

> Caveat: the quadrant analysis used the atom-space rigidity argument that §6 retracts, so "12 is
> the wall" on that branch is *not* established. It starts 3 equations behind, so it was not chased.

### 11.3 Additional negative checks (this session, verified)
- **No variable and no atom has its entire equation footprint inside the current 9 failing
  equations** — the free-knob well at this state is dry. The best single equation to re-admit is
  6816, which unlocks only 3 knobs (`x_1329, x_10903, x_29854`, the last already used).
- **Sacrificing the load pins to unlock the mod-p movers does not pay**: adding `eq(3576)` and/or
  `eq(3578)` grows `S` to 26–28 equations while yielding only 8–10 footprint-contained knobs
  (and still not `x_6418`/`x_12553`), so 18–20 equations would have to be zeroed to beat 9.
- Certificate-guided branch-and-bound in *atom* space: **no relaxation set of cost ≤ 26** makes the
  mod-p system consistent (2,573 nodes explored, `s9/chain2/bb.py`).
- Direct computation: in the atom-space formulation `A` and `B` are **first-order rigid mod p** —
  none of `(δ_{x7068}, δ_{x4432}) = (1,0), (0,1), (1,1)` is reachable (`s9/chain2/ab.py`). The
  solver was sanity-checked: the homogeneous system is consistent with 1,728 pivots and it does
  find genuinely free directions (`x_33612`, `x_3387`).
- Equation-space first-order solve over free-input ripple directions is inconsistent both with and
  without the square-root rows (`s9/kernel/eqspace.py`, `eqspace2.py`) — the 39,024 gain comes from
  *deliberately breaking confined atoms*, which that formulation does not model.

### 11.4 Standing status of "is 9 optimal?"
Proved exhaustively only inside `S = FAILS ∪ {9123, 18673}` (0 of 1287 5-subsets and 0 of 1716
6-subsets integer-solvable). Outside that region the support is two saturated beam searches plus a
marginal-cost table showing the next mod-p knob costs +2 equations for +1 dof. **Not a proof** —
and given that this session's earlier optimality claim was wrong for exactly this kind of reason,
treat it as a working hypothesis.

---

## 12. The lock in closed form (the deepest result of the session)

Deliverable still **39,024 / 39,033**. What follows does not beat it, but it replaces "the core is a
wall" with an exact statement of *what* the setter arranged, derived rather than searched.

### 12.1 The circuit's output wire is a single free input

`s9/modtrace.py` walks the gate DAG keeping only the monomials that survive mod p (dropping products
against the p-wire and against currently-zero variables). Applied to `x_12186` it collapses a
~20-deep cone to one chain, terminating at

> **`x_12186 ≡ x_22649 (mod p)`, and `x_22649` is a FREE INPUT — it drives `x_12186` exactly 1:1.**

Confirmed numerically (+1, +2, +1000 all give `dx_12186 = +δ`). Its only disturbances are atom 2423
(`12604395·(x_22649 − x_29524) − x_9899`), the core, and two composites.

### 12.2 Chain 1 can be killed outright — a first

`x_29524` is driven 1:1 by the free input `x_22152`. So moving
`{x_22649, x_22152, x_14853, x_7068}` by the *same* δ = K1 − x_12186 drives the computed value to the
pinned constant while holding the core's control difference `u = x_29322` at **0**:

```
x_12186 = K1  ✓   x_14853 = K1  ✓   x_7068 = K1  ✓   u = 0  ✓   atom 22229 (A) = 0  ✓
```

Only 3 atoms remain nonzero (`s9/chain1.py`): 22231 (chain 2), 37887 (its square), and the load pin
**31670**, which is what stops it:

```
31670 :  x_24601·(x_22152 − 126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506) = 7550763·x_29309,   x_29309 = p·x_105
```

The constant is **exactly the original `x_12186`**. So the bit `x_24601` pins the circuit's output
wire to a value that differs from `K1` by precisely `D1`. That is the trapdoor, stated exactly.

### 12.3 The MUX branch makes BOTH defects vanish

`x_2099` and `x_19964` are outputs of a 2-bit MUX on `(b1, b2) = (x_4287, x_2081)`:

```
x_2099  = b1(1−b2)·x_31861 + b2(1−b1)·x_6418  + b1b2·x_9118
x_19964 = b1(1−b2)·x_14865 + b2(1−b1)·x_12553 + b1b2·x_8731
```

At `(0,1)` they read the *pinned* `x_6418`, `x_12553`. **Setting `x_4287 = 1` selects `x_9118`,
`x_8731` — free inputs.** Then C1 and C2 close by construction:

> **atoms 22229 = 0 AND 22231 = 0 simultaneously — the first time both defects have ever vanished.**

### 12.4 The branch's price, and why it still does not close

`x_4287 = 1` lights `x_21279 = 1`, activating three loads (atoms 19088, 22233, 22235) on
`x_2239, x_31731, x_9106`. Those are **three combinations of two base quantities** — the same
3-from-2 shape as L1/L2/L3 from S,T (verified exactly):

```
x_2239  = 3494591·x_27177 + 14240157·x_4306
x_31731 = 15964591·x_27177 + 13881285·x_4306
x_9106  = 7204959·x_27177 + 6822253·x_4306
```

so the whole obligation is **`x_4306 ≡ 0` and `x_27177 ≡ 0 (mod p)`**, where

```
x_4306  = (x_8731 + x_14865)(x_6418 − x_31861) − (x_12553 − x_14865)(x_31861 − x_9118)
```

and it **is solvable**: a 2×2 mod-p solve on `(x_9118, x_8731)` zeroes both (verified,
`s9/mux2.py`, `s9/mux3.py` — reaches 6 nonzero atoms). Structural fact: `∂x_27177/∂x_8731 = 0`, so
`x_27177` depends on `x_9118` alone, which **determines `x_9118` mod p uniquely**:

```
x_9118  must be  33371159155735472537534252650716501592825364489306217536352743247010353604716   (mod p)
```

But the mirror + core chain forces `x_9118 ≡ x_7068 ≡ x_14853 ≡ x_12186`, and `x_12186` is pinned by
§12.2 to `82007976112976807461901870199198737303514020147647909878034348606308756230357` when
`x_24601 = 1`, or to `0` when `x_24601 = 0`. **Neither matches.** Checked directly.

> **That mismatch is the trapdoor.** Every reachable branch pins the circuit's output wire to a
> constant, while the branch's own obligation pins the verified value to a different one. The
> setter's witness must live where the two agree — which needs the secret, or a genuine break of
> the pinned residue.

Best measured on this branch: **34 failing** (`s9/mux3_out.json`), versus 9 on the current one.

### 12.5 What this changes for a future attack
The target is now a single scalar identity, not a search: make the residue that `x_24601`'s pin
puts on `x_22649` agree with the residue the branch obligation puts on `x_9118`. Both sides are
computable in closed form from the instance (§12.2, §12.4). Any further progress should attack that
identity — e.g. other bits that move either residue, or a quadrant where the two pins coincide —
rather than searching assignments.

---

## 13. Deliverable 39,026, and the complete quadrant map

**Verified deliverable: `best/new_instance_partial_39026.json` = 39,026/39,033 (7 failing:
`[12231, 12270, 12350, 14584, 18673, 22044, 29125]`).** Checked by `checker.py` and by the
independent AST verifier. 28 variables differ from the 39,024 witness.

### 13.1 How the last two equations fell: effective vs syntactic footprints

Every knob census up to now asked *which equations mention a variable*. That is too conservative:
a variable multiplied by a currently-**zero** variable has no effect at all, and one moved with the
gate-preserving ripple only disturbs the check atoms it actually reaches. Census over all 38,748
variables: effective footprint is **22.3 % smaller** than syntactic (599,740 vs 771,576
variable–equation incidences); **25,250** variables strictly gain and **3,289** are pure no-ops.

That admits `x_9118` and `x_8731` as knobs beside the eight syntactic ones, and by CRT
(`gcd(5113045, p) = 1`; `x_8731` moves atom 35761 in steps of 1) **two of the four binding
congruences dissolve**:

| knob set | lattice invariant factors | congruences | max zeroable | failing |
|---|---|---|---|---|
| 8 syntactic | `[1,1,P,1,P,1,P,7376877·P]` | 4 | 4 | 9 |
| **10 effective** | `[1,1,P,1,1,1,1,7376877·P]` | **2** | **6** | **7** |

Exhaustive over all 2¹³ subsets of the sacrificed region: max simultaneously zeroable is 6 and no
7-subset is integer-solvable, so **7 is optimal for this defect placement**. The two survivors are
exactly the original core congruences; of all 38,748 variables only **25** move either, the cheapest
costing 10 extra equations (best achievable ≥ 14). Going below 7 needs the trapdoor, not a search.

### 13.2 The quadrant map — every branch violates something

The instance has three relevant control bits: the MUX pair `(b1, b2) = (x_4287, x_2081)` and the
output-wire pin bit `x_24601`. Measured end-to-end:

| configuration | A (22229) | B (22231) | core | what breaks instead | failing |
|---|---|---|---|---|---|
| `b1=0, b2=1, x_24601=1` (**current**) | ≠0 | ≠0 | dead (`u=w=0`) | the two defects themselves | **7** |
| `b1=0, b2=1, x_24601=1` + chain-1 kill | **0** | ≠0 | dead | output-wire pin 31670 | 21 |
| `b1=1, b2=1` (MUX flip) | **0** | **0** | live | branch obligation forces `x_9118` to a residue ≠ the pinned one | 34 |
| `b1=0, b2=0, x_24601=0` | **0** | **0** | **dead** (`x_15298=0`) | setter pins 688/1618 demand values from wires that are dead there | 28–31 |
| `x_2081=0` alone | — | — | dead | mirror/pin chain | 12 |
| `x_24601=0` alone | — | — | — | — | 27 |

> **This is the trapdoor, completely mapped: it is a covering design. Each branch satisfies some
> pins and kills the wires the others need. No reachable configuration satisfies all of them, and
> the conserved quantity is not any single residue — it is the *mismatch* itself, which relocates
> between the defect atoms, the output-wire pin, the branch obligation, and the setter pins.**

On the all-zero branch the position is the sharpest ever reached — `A = B = 0`, `x_15298 = 0`,
`u = w = 0`, the forced OR gate `x_9274 = 1` satisfied — with only the two setter pins 688/1618
left. Their carriers `x_24468`, `x_18956` *are* movable there (by `x_14393`/`x_14853` and by four
free inputs respectively), but the boolean needed to satisfy the OR gate lights its own two load
pins, and closing those cascades (31 → 45 failing). That trade is the remaining obstacle on that
branch.

### 13.3 Standing recommendation
Both remaining congruences are `D1` and `D2` mod p. Everything measured says they cannot be moved
cheaply from inside any single branch. The productive question is now cross-branch: find a
configuration in which the setter pins that a branch *keeps* are simultaneously satisfiable with
the wires that branch keeps *alive* — a covering/exact-cover question over the 256 gate-bits and
their 512 pins, informed by the quadrant map above rather than by local search.

---

## 14. The all-zero branch driven to two atoms — and the circular closure

Deliverable unchanged at **39,026**. This section records how far the all-zero branch goes and the
exact reason it closes on itself.

### 14.1 The pipeline
Starting from `x_2081 = 0, x_24601 = 0` (the only quadrant where the MUX output and the circuit's
output wire are *consistently* both 0):

| step | remaining atoms | failing |
|---|---|---|
| flip both control bits, zero the free carriers `x_22162`, `x_30213` | 688, 1618, 23000, 39067, 40608, 41211 | 28 |
| satisfy the forced OR gate with gate-bit `x_47` **and close its two load pins** | **688, 1618, 40608** | **16** |
| close setter pin 688 (`x_16742` drives `x_18956` 1:1) and 1618 (`x_14681` drives `x_24468` 1:1) | 9193, 26731, 39614 | 30 |
| close those two mirrors (`x_38667` → `x_19083`, `x_29851` → `x_24483`) | 14061, 37735, 41882 | 26 |
| automated mirror-chase, 3 further steps | **23824, 23826, 40047** | **17** |

At every stage `A = 22229 = 0`, `B = 22231 = 0`, `x_15298 = 0` (core dead), `u = w = 0`, and the
forced OR gate holds. **All the conditions that block the main branch are satisfied here.**

### 14.2 Why it closes on itself
Each obstruction has the shape `k·(F − C) − handle`, where `F` is a free input just moved and `C` a
computed mirror; it closes by finding the non-boolean free input that drives `C` exactly **1:1 mod
p** and moving it. `s9/chase.py` automates this. The chain is *short* — it terminates after three
steps, at

```
23824 : x_47·(x_24221 − HUGE) = 982875·x_27050
23826 : x_47·(x_25477 − HUGE) = x_15160
```

— **the two load pins of `x_47` itself**, the gate-bit that had to be switched on to satisfy the
forced OR gate. The mirror chain from the setter pins terminates precisely on the variables pinned
by the bit the branch needs.

> **This is the covering design closing the loop.** The branch needs a gate-bit ON (only the 256
> gate-bits can satisfy `x_9274 = 1`; all 900 pin-free booleans fail — verified exhaustively). Every
> gate-bit pins two inputs. And the mirror chain that closes the setter pins runs through exactly
> those inputs. `s9/solve_branch.py` parameterises the whole pipeline by which gate-bit is used, so
> the question "is there a gate-bit whose pinned pair is disjoint from its own mirror chain?" is a
> finite 256-way check.

### 14.3 Tools
`zero.py` → `zero9.py` (staged construction), `chase.py` (automated mirror-chase with cycle
detection), `solve_branch.py` (whole pipeline parameterised by the OR-gate gate-bit),
`modtrace.py` (mod-p symbolic cone), `reduce_size.py` (residual problem measurement).

### 14.4 The chain is FORCED — the branch is closed, rigorously

Enumerating *all* non-boolean free inputs that drive each chain target exactly 1:1 mod p:

```
x_18956 <- x_16742        x_24468 <- x_14681        (the two setter-pin carriers)
x_19083 <- x_38667        x_24483 <- x_29851
x_28922 <- x_25477        x_37571 <- x_30709        x_462 <- x_24221
x_24221 <- x_24221        x_25477 <- x_25477        (self-driving: they ARE free inputs)
```

**Every step has exactly one driver.** There is no alternative routing, no branching, nothing to
search over. The chain is a single forced thread, and it terminates precisely on `x_24221` and
`x_25477` — the two inputs pinned by `x_47`.

Together with §14.2 that closes the branch:

1. The forced OR gate `x_9274 = 1` can only be satisfied by a **gate-bit** (all 900 pin-free
   booleans fail; verified exhaustively over all 1,156).
2. Of the 256 gate-bits, only `x_47` makes the drivers `x_16742 → x_18956` and `x_14681 → x_24468`
   live, i.e. only `x_47` lets the setter pins 688/1618 close at all (others leave 688 or 1618
   standing — `s9/solve_branch.py`).
3. The mirror chain that closing them requires is **forced** (one driver per step).
4. It ends exactly on `x_47`'s own pinned pair.

> The design is self-referential: the only bit that unlocks the setter pins is the bit whose own
> pins the unlocking chain must violate. That is why the branch stalls at 17 failing with three
> atoms — and it is a much sharper statement than "search failed".

---

## 15. Second-order freedom, and a link-by-link audit of the lock

Deliverable unchanged: **39,026 / 39,033**. This section closes the last idea I had, rigorously.

### 15.1 The blind spot: dormant product gates
Every census, Jacobian and driver-enumeration in this campaign is **single-variable**. That is
structurally blind to a *dormant product gate*: a monomial `x_i·x_j` with both factors currently
zero. Moving either factor alone gives `0·δ = 0`, so no linear analysis sees it; moving **both**
injects `δ_i·δ_j`. At the 39,026 witness there are **81,681** strictly dormant monomials (both
factors exactly 0 — not merely `0 mod p`, which the p-wire also satisfies), **2,042** of them with
both factors free inputs.

Scanning all 2,042: exactly **2** change the residues `(D1, D2) mod p`, and both are the MUX
control `x_4287` paired with `x_31861`/`x_14865` — i.e. the MUX flip already known from §12.3.
So the second-order freedom exists in bulk but does not reach the residues.

Also worth recording: at the 39,026 witness **`D2 ≡ 0 (mod p)` already** — atom 22231 is exactly
zero. The whole remaining obstruction is the *single* congruence `A = atom 22229 ≡ 0`, i.e.
`x_7068 ≡ K1 (mod p)`.

### 15.2 The chain that pins `x_7068`, audited link by link

```
pin 31670 -> x_22152 -> x_29524 -> (2423) -> x_22649 -> x_12186 -> (core u=0)
          -> x_14853 -> (mirror 29539) -> x_1308 = x_7068 -> A
```

`s9/linkaudit.py` reports each link's handle and its granularity:

| link | handle | structure | granularity |
|---|---|---|---|
| 31670 | `x_29309 = x_105·x_3915` | `x_3915 = p` (wire), `x_105` free | **p** |
| 2423 | `x_9899 = x_14466·x_14768` | `x_14466 = p`, `x_14768` free | **p** |
| 22772 | `x_13595 = x_2121·x_11648` | **both factors 0 — dormant!** | *unquantised* |
| 26729 | `x_25758 = x_10603·x_33612` | `x_10603 = p`, `x_33612` free | **p** |
| 29539 | `x_29967 = x_11360·x_30163` | `x_11360 = p`, `x_30163` free | **p** |
| 22229 | `x_642 = x_17325·x_28599` | `x_28599 = p`, `x_17325` free | **7376877·p** |

Every link is p-quantised — **except one**. Link 22772 carries a genuinely dormant handle, and
activating it would let `x_29524` shift by an arbitrary amount, breaking the chain and giving
`x_7068 ≡ K1` outright.

### 15.3 The one dormant link, and its guard
`x_2121 = x_38144·x_13636` with `x_13636 = x_24601 = 1`, so `x_2121 = x_38144`. Setting
`x_38144 = 1` does exactly what is wanted, and both at once:

* `x_8211 = x_13636·(1 − x_38144)` collapses to 0, **disconnecting** the pinned `x_22152` from
  `x_29524`;
* `x_13595 = x_2121·x_11648 = x_11648` activates, **connecting** the FREE input `x_11648`.

Executed (`s9/dormsolve.py`): `x_12186 = K1`, `x_7068 = K1`, `u = 0` — the chain is broken exactly
as predicted. But `x_38144` is held by atom **32288 = `x_38144 − 0`**: a *bare* single-variable
pin, no product, no handle, nothing to absorb. It costs 11 equations and can never be satisfied
with `x_38144 = 1`.

> **The setter put a handle-free constant pin on the one select that would have opened the chain.**
> That is the last brick. Every other link is p-quantised; this one is nailed shut directly.

### 15.4 What this establishes
The lock is now verified *link by link* rather than inferred from failed searches: there is exactly
one non-rigid link in the chain that pins `x_7068`, and it is guarded by a bare pin. Combined with
§14.4 (the all-zero branch's mirror chain is forced, one driver per step) and §13.1 (7 is
exhaustively optimal for the current defect placement), the instance's resistance is characterised
at the level of individual gates.
