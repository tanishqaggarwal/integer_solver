# Session 9 — the obstruction, fully characterised (and a proof that 39,022 is a local optimum)

Verified deliverable unchanged: **39,022 / 39,033**
(`best/new_instance_partial_39022.json`, `python3 checker.py best/new_instance_partial_39022.json`).

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

`x_28599 = x_17499 = x_26064 = p = 2^256 − 2^32 − 977` **exactly** (the secp256k1 field prime), so
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

## 6. NEW: 39,022 is provably optimal for this defect structure

The handles absorb exactly the `p`-multiples, so the reachable defect values are
**`A ∈ D1 + pℤ`** and **`B ∈ D2 + pℤ`** — and nothing else:

* `x_17325` shifts `A` by `7376877p`; `x_9413` shifts `B` by `p`.
* `x_7068` may additionally move by any multiple of `p` (the mirror handle `x_30163` absorbs it);
  likewise `x_4432` via `x_11052`. Non-multiples of `p` break a 12–16-equation mirror atom.

Every failing equation evaluates to `m·(c₁A + c₂B)` (all other atoms vanish), so it can be rescued
only if `c₁·D1 + c₂·D2 ≡ 0 (mod p)`. Computed for all nine linear failing equations:

| eq | c₁ | c₂ | `c₁D1+c₂D2 mod p` = 0? |
|---|---|---|---|
| 2554 | 1 | −16 | no |
| 6816 | −15 | −12 | no |
| 8124 | 36 | −38 | no |
| 9421 | 13 | 1 | no |
| 12231 | 18 | 8 | no |
| 12270 | −31 | 6 | no |
| 12350 | −23 | 29 | no |
| 14584 | 17 | −33 | no |
| 22044 | −24 | −29 | no |

None can vanish: the required ratio is `−c₁/c₂` (a rational with numerator and denominator ≤ 40)
while the actual `D2/D1 mod p = 29086885819837044927165644576879200888791656444895144487473726012333934270214`
is a full 256-bit residue. The remaining two equations (8680, 29125) need `B = 0` **exactly**,
impossible since `B ≡ D2 ≢ 0 (mod p)`.

And every alternative placement of the two defects is strictly worse:

| placement | failing equations |
|---|---|
| **22229 + 22231 (+37887)  — current** | **11** |
| 22229 + 3578 | 23 |
| 29539 (+40826) + 22231 (+37887) | 24 |
| 3576 + 22231 (+37887) | 24 |
| 22229 + 7930 (+41512) | 25 |
| 29539 (+40826) + 7930 (+41512) | 29 |

> **Conclusion: 39,022/39,033 is a local optimum, not merely the best found. Any improvement
> requires satisfying `x_7068 ≡ K1` and `x_4432 ≡ K2 (mod p)` outright — i.e. cracking the core.**

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
