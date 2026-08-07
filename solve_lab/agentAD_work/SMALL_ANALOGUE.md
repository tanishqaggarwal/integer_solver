# SMALL_ANALOGUE.md — §8 settled in the small, and it is DEAD

Agent AD. Everything below was computed in this directory. No git commands were run.
`PYTHONDONTWRITEBYTECODE=1` throughout. Footprint: ~350 KB, no shared table touched.

**Rule-2 compliance.** The small system was built from the *decoded mathematics* — the
five-atom integer law recorded in `agentW_work/RESUME_W.md`, the mux quadrants and the
merge-tree shape from the campaign brief. `EQUATIONS.txt` was never opened by any script
here; no ordering, coefficient template, or emission artefact was used or inferred.

---

# 0. THE VERDICT

> ## §8 is DEAD.
> The integer-lift discipline (`c·P | R`) and the off-pin conditions **cannot** produce a
> useful upper bound on `w`, and the direction of the weak dependence they do have is
> **the opposite of what the §8 hypothesis needs at the top of the weight range.**

Four findings, in the order of how much each one closes:

1. **The off-pin family carries no signal at all, and is not even a constraint.** Its
   variables are private to a dead block and enter no other condition (verified
   adversarially: 11 987 condition values re-evaluated after driving every dead block's
   private outputs to random nonzero multiples of `P`; **0 moved**). It closes for every
   subset at every `|S|`. Its *active-condition count* is `2(n−|S|)` exactly — it **decreases**
   with `|S|`, so if it constrained anything it would be a *lower*-bound mechanism.

2. **The congruence family does depend on `|S|`, but only through a per-live-block
   constant.** Live merge blocks number exactly `|S|−1` (closed form, checked against the
   count for every subset at every size: 0 mismatches). The closure rate behaves as
   `ρ^{|S|−1}` with `ρ` a per-block constant, measured `0.55 … 1.00` (robust fit, 340 draws)
   across every curve, size, tree shape, coefficient pool and draw tried — and **flat in `n`**
   over five sizes, which is what §5 says it must be.

3. **The ABSORBING-BLOCK THEOREM explains why, and caps how bad it can get.** If `A = i1−i2`
   is invertible modulo the block's small modulus `c`, the block's *own* output lifts can
   satisfy all three conditions for **any** input lifts — the block constrains nothing.
   Since `A = a + P·δ` with `δ` a free lift difference and `P` invertible mod any small `c`,
   even a bad `a` is usually repairable. Failure needs a per-block coincidence; it never
   accumulates into a cutoff.

4. **At the real instance's own parameters the effect vanishes.** The real instance has a
   nontrivial modulus on `927 / 3707 = 25 %` of its lift conditions. At `frac = 0.25` the
   majority of coefficient draws close **all `2^n` subsets**, and in the draws that do not,
   the failures sit in the **middle** weights while the top weights are untouched — the
   closure-rate curve is **U-shaped, returning to 1.0000 at `|S| = n`.** Large `|S|` is the
   *easy* end.

**Extrapolation to `n = 256`** (model validated first, see §6 — over **340 draws** it predicts
`Bmax` exactly in 74 % and to within 1 in 99.4 %): taking the most constraining `ρ` seen
anywhere — including modulus pools chosen adversarially, four times more constraining than the
real instance — and holding it fixed, the implied bound is **`w ≤ 205`**. The free unconditional bound is `w ≤ 255` (AB §7). §8 needs `B ≲ 56` to beat
rho and `B ≲ 24` to be actionable. **It is short by ~170.**

**What would have to be true instead** (§7): a bound `w ≤ 56` needs `ρ ≤ 0.093`, i.e. every
live block binding with probability `> 90 %`. That is self-defeating: at `ρ = 0.093` the
*true* subset, of weight 56, would itself close with probability `2^−188`. An instance
tight enough to bound `w` this way is an instance that cannot be satisfiable by accident.

**Correction to the fleet's reading, both ways.** AB is right that the closures at
`|S| = 1..64` do not rule out an upper bound — but the reason they do not is that
**closure is not a function of `|S|` at all in the way the fleet has been reading it**: it is a
product over live blocks, and the number of *subsets* at weight `w` grows as `C(n,w)`, which
swamps `ρ^{w−1}` by hundreds of orders of magnitude at `n = 256`. And T's `|S| = 128` stall
is, in the analogue, **not** what an upper bound looks like: an upper bound would show up as
the top weights emptying, and the top weights are the fullest.

---

# 1. The analogue, and why it is the analogue

## 1.1 The law

`agentW_work/RESUME_W.md` records the integer form of every one of the 383 blocks:

```
congruence k=1,2,3 :   a_k  * L     * ( c_k1*N1 + c_k2*N2 )  =  c_k * P * u_k
off-pin    j=5,6   :   a'_j * (1-L) *   i_j                  =  c'_j* P * u'_j
A  = i1-i2   B = i4-i3   E = i1+i2+i5+Q
N1 = E*A^2 - B^2         N2 = A*(i3+i6) - B*(i2-i5)
```

Identifying `i1 = x_R, i2 = x_L, i3 = y_L, i4 = y_R, i5 = x_out − Q, i6 = y_out` turns
`N1 = N2 = 0 (mod P)` into exactly the chord addition law:

```
N1 = 0  <=>  x_L + x_R + x_out = lambda^2 ,  lambda = B/A = (y_R - y_L)/(x_R - x_L)
N2 = 0  <=>  y_L + y_out = lambda*(x_L - x_out)
```

`Q` is a global additive offset on the x-slot. Consistency of the ladder **forces `Q = 0`**:
with `Q ≠ 0`, `N1` still vanishes (the `Q` cancels inside `E`) but `N2` picks up `−b·Q`.
This is not a modelling choice, it is checked (`ad_validate.py` V1b: 400 random chord
additions, `Q ∈ {0,1,7}`, zero spurious consistencies).

Mux: `cA = a(1−b)`, `cB = b(1−a)`, `cC = ab`; the node output is the arithmetic combination
`X_v = cA·X_L + cB·X_R + cC·i5`. The gate `L` multiplying the gadget output is `cC = ab`
(W: 383/383 blocks, output pair consumed by exactly that multiplier), so the gadget law is
enforced iff **both** children are live, and the off-pins fire otherwise.

## 1.2 The integer-lift discipline — the part the experiment is about

`u_k` is a free integer and `gcd(a_k,P) = gcd(c_k,P) = 1`, so the atom's exact integer
content is

```
c_k*P | a_k*L*Z_k    <=>    P | L*Z_k       (the mod-P law: ordinary curve arithmetic)
                     AND    c_k | a_k*L*Z_k (the small-modulus condition, invisible mod P)
```

`Z_k = c_k1*N1 + c_k2*N2`. On a live block the mod-P law holds, so `Z_k = P·M_k`, and the
second condition is `m_k | M_k`, `m_k = c_k/gcd(a_k,c_k)`.

**`M_k` is not determined by the curve arithmetic.** Every slot is an integer variable fixed
only modulo `P`: `slot = alpha + P·t` with `alpha ∈ [0,P)` the reduced coordinate and `t ∈ ℤ`
a free **lift**. Writing `a = a1−a2`, `b = a4−a3`, `e = a1+a2+a5`, `g = a3+a6`, `h = a2−a5`,
`n1 = (e a² − b²)/P`, `n2 = (a g − b h)/P`, `δ = u1−u2`, `β = u4−u3`, `ε = u1+u2+u5`:

```
nu1 = N1/P = n1 + (2ae*d + eps*a^2 - 2b*bt) + P*(e*d^2 + 2a*eps*d - bt^2) + P^2*(eps*d^2)
nu2 = N2/P = n2 + (a*gam + d*g - b*eta - bt*h) + P*(d*gam - bt*eta)
```

(all terms kept — the `P²` term and the cubic `eps·d²` term included; see the trap list in
§2), and the condition is `m_k | c_k1·nu1 + c_k2·nu2`, `k = 1,2,3`.

Differentiating recovers **W's Jacobian independently**:

```
d(nu1,nu2)/d(u5,u6) = [[A^2, 0], [B, A]] ,  A = a + P*d ,  B = b + P*bt
i.e.   nu1 = nu1a + u5*A^2 ,   nu2 = nu2a + u5*B + u6*A
```

Lifts couple along the tree — a node's output lift **is** its parent's input lift — so closure
of the whole family is a tree DP over lift residues. By CRT the DP factors over the prime
powers dividing the moduli; each factor is an exact DP over `(ℤ/q)²` states. Exact, not
heuristic.

## 1.3 Curves and trees

`ad_curves.py`, `ad_curves2.py`. Short Weierstrass over prime fields, **prime order `N` of
exactly `n` bits**, with `n` leaf selectors carrying `2^i·G`, `i = 0..n−1` — matching the real
instance's `n = 256 = bitlen(N)`. Independently verified for every curve used: `G` on curve,
`N·G = O`, `N` prime, `bitlen(N) = n`.

| n | p | curves used |
|---|---|---|
| 8 | 241 | `a=0,b=112` and `a=0,b=110` (**j=0**), `a=126`, `a=105`, `a=141` |
| 10 | 1009 | `a=0,b=354` (**j=0**), `a=540,b=606` |
| 12 | 4093 | `a=0,b=335`, `a=0,b=306` (**j=0**), `a=3338`, `a=2343`, `a=1584` |
| 14 | 16363 | `a=0,b=14183` (**j=0**), `a=11097,b=3966` |
| 16 | 65521 | `a=0,b=345`, `a=0,b=51` (**j=0**), `a=14033`, `a=65072`, `a=38869` |

At least one `j = 0` (CM by `√−3`, the real instance's class) at every size, plus generic
curves at every size. **`j = 0` makes no measurable difference** — §5.

Trees: balanced (matching the real depth-9-over-256 shape) and **skew 178:78** (the real root
split), both tested everywhere. Leaf indices are assigned by a greedy sum-balancing split so
that every proper subtree's maximum partial sum stays below `N`; this is what makes the
encoding faithful, and the real instance gets it for free because `N/2^256 = 1 − 2^−128`.

---

# 2. VALIDATION BEFORE MEASUREMENT

`ad_validate.py`, `ad_v3b.py`, `ad_theorem.py`. All pass.

| check | what it rules out | result |
|---|---|---|
| **V1a** law `N1=N2=0 mod P` on 400 real chord additions | wrong slot identification | 0 violations |
| **V1b** the same with `Q ∈ {1,7}` | `Q = 0` being a free choice | 0 spurious consistencies |
| **V2a** `P | N1, P | N2` for random integer lifts | the lift model being ill-posed | 0 bad |
| **V2b** closed form for `nu1,nu2` vs **raw big-integer** `N1/P, N2/P` | **agent N's failure mode: discarding the quadratic part or truncating a fourth power.** The `P²` term and the cubic `eps·d²` term are compared term-for-term | 0 bad |
| **V2c** the Jacobian form `nu1 = nu1a + u5·A²`, `nu2 = nu2a + u5·B + u6·A` | an algebra slip in the DP's fast path | 0 bad |
| **V3** tree DP vs **exhaustive raw-integer enumeration of every lift assignment**, `n = 4,5`, both tree shapes, 3 seeds, `q ∈ {2,3,4}`, every subset | the DP being a leaky abstraction (AB's `d_reg` failure mode) | **3 456 comparisons, 0 mismatches** |
| **V3b** the same at `n = 7`, where the shape of the rate-vs-`|S|` curve is already visible; the **per-weight rate** is compared, not just the verdict | the headline curve being a DP artefact | **3 072 comparisons, 0 mismatches**; the brute-force per-weight counts equal the DP's |
| **V4** `n = 8` and `n = 12` **in full**: solution set from the equations vs `{S : k(S)·G = T}` from curve arithmetic, 4 curves × 2 tree shapes each | the model not being the problem it claims to model | **match in 16/16 configurations** |
| **V5a** awkward plant: a target with **two** valid subsets (`k₀ = 22` and `k₀+N = 233` at `p = 241`, weights **3 and 5**) | a solver that stops at the first hit | both recovered |
| **V5b** awkward plant: rigged coefficients (`2 | nu1` at the root block only) | agent X's vacuity failure: a plant every scan point hits | failure set is **3 of 64** subsets — proper and non-empty — and the DP reproduces it **exactly** (symmetric difference 0) against raw-integer brute force |
| **T1** absorbing-block theorem, brute force on real curve points, random `C`, random lifts | the theorem in §5 | 2 515 `gcd(A,c)=1` cases, **0** with no satisfying `(u5,u6)` |
| **T2** an independent, DP-free **greedy constructive certifier** (bottom-up, raw integers) | the DP over-reporting closure | thousands of certificates, **0** contradicted by the DP |

**V4's exceptional-block census, checked against a closed-form expectation.** The gadget law
has no doubling case: `A ≡ 0` forces `B ≡ 0` (W §2), so a merge of two *equal* points
(`k_L ≡ k_R mod N`) leaves the output completely free — DEGENERACY — and a merge of two
*opposite* points (`k_L ≡ −k_R`) is integrally infeasible — CANCEL. Prediction: with
sum-balanced leaf assignment every proper subtree sums below `N`, so DEGEN is impossible and
CANCEL can happen only at the root, at the unique subset with `k(S) = N`. **Measured: DEGEN
subsets 0, CANCEL subsets 1, in all 16 configurations.** That single subset has
`k(S)·G = O ≠ T`, so it is excluded from both sides and the solution sets coincide exactly.

> Without that census the model would have been wrong and every number after it noise: a
> degeneracy makes the equation solution set a **strict superset** of `{S : k(S)G = T}`.

---

# 3. THE MEASUREMENT — closure rate vs `|S|`, full curves

`ad_measure.py`. Exhaustive over all `2^n` subsets. The mod-P-infeasible subset (the CANCEL
one) is excluded so that what is measured is purely the **integer-lift** family.

**Live-block count is `|S|−1` exactly** — checked for every subset at every size, 0 mismatches.
So `|S|` enters only through the number of live blocks, which is what makes the whole
question a per-block question.

## 3.1 At the real instance's own parameter (`frac = 0.25` nontrivial moduli, `927/3707`)

`n = 12`, `p = 4093`, pool `{2,3,4,5,7,8,9,11,13}`, 8 draws; counts of **closing subsets** by
weight, against `C(12,w) = 1,12,66,220,495,792,924,792,495,220,66,12,1`:

```
draw 0 : 1,12,63,202,450,732,879,774,492,220,66,11,1
draw 1 : 1,12,61,190,420,692,849,762,490,220,66,11,1
draw 2 : 1,12,66,220,495,792,924,792,495,220,66,11,1     <- every subset closes
draw 3 : 1,12,64,208,465,752,894,780,493,220,66,11,1
```

(the `11` at `w = 11` is `12 − 1`, the excluded CANCEL subset).

> **The deficit is entirely in the middle weights. `w = 9,10,11,12` are FULL in every draw.**
> The closure-rate curve is U-shaped and returns to `1.0000` at `|S| = n`. **High `|S|` is the
> easy end.** This is the exact opposite of what §8 requires.

`n = 16`, same pool, one draw, rate by weight (`w = 0..16`):

```
1.0000 1.0000 .9750 .9357 .8923 .8526 .8216 .8021 .7955 .8021 .8216 .8525 .8923 .9357 .9750 1.0000 1.0000
```

Symmetric about `w = 8` to four decimals, minimum at the centre, `1.0000` at both ends.

At `n = 8` with the same setting, **11 of 12 draws close all 256 subsets outright.**

## 3.2 Under stress (`frac = 1.0`: a nontrivial modulus on *every* condition)

This is far more constrained than the real instance. `n = 12`, `p = 4093`, pool
`{2,3,4,5,7,8,9}`, aggregate over 6 draws:

```
 w    :  0     1     2     3     4     5     6     7     8     9    10    11    12
 rate : 1.000 1.000 .8409 .7826 .7037 .6420 .5979 .5644 .5391 .5159 .4798 .4394 .3333
```

Now monotone decreasing — but the decay is a *constant factor per live block*: the ratio of
successive rates is `.931 .899 .912 .931 .944 .955 .957 .930 .916`, i.e. `ρ ≈ 0.93` per extra
live block, stable across the whole range. Per-draw `Bmax` (largest weight with any closing
subset): `11,12,11,11,11,12` — **`n − Bmax ∈ {0,1}`.**

At `n = 8` under stress, `Bmax ∈ {6,7,8}`, `n − Bmax ∈ {0,1,2}`.

## 3.3 Adversarial modulus pools

`ad_scale.py` was swept over pools chosen to make the discipline as tight as the family
allows, all at `frac = 1.0`, `n = 8..14` exhaustive (`ad_scale_pools.log`):

| pool | robust `ρ` (min over draws) | worst `n − Bmax` |
|---|---|---|
| `{2}` | 0.66 | 0–3 |
| `{2,3}` | 0.67 | 0–3 |
| `{4}` | 0.64 | 0–2 |
| `{8}` | **0.55** | 4–5 |
| `{2,4,8}` | 0.57 | 0–4 |
| `{3,9}` | 0.55 | 1–3 |
| `{2,3,5,7,11,13}` | 0.64 | 0–2 |

**Concentrated small prime powers are the worst case; a wide pool is weaker,** because each
condition draws only one modulus and a block is then rarely bound by all three. Nothing in the
family reaches the `ρ ≤ 0.093` that §7 shows would be needed, and §5 says nothing can:
`M = 2^k` binds only when `A` is even, so `ρ ≥ 1/2` however large `k` gets — which is exactly
the 0.55 floor measured at `{8}`.

## 3.4 Tree shape and curve

Balanced vs skew-178:78 changes the individual numbers and nothing structural (both appear
throughout §3.1–3.2; e.g. at `n = 8` stress, balanced `Bmax = 6..8`, skew `Bmax = 7..8`).

**`j = 0` (CM by `√−3`, the real instance's class) makes no measurable difference**
(`ad_report.py` §3): over 340 draws, `ρ` median **0.7881** on `j = 0` curves against **0.7980**
on generic curves, with overlapping ranges. The real instance being `j = 0` is not what is
holding §8 up.
Five curves per size at `n = 8,12,16`, two at `n = 10,14`; the spread across curves is smaller than the spread
across coefficient draws.

---

# 4. THE OFF-PIN FAMILY, SEPARATELY

`ad_offpin.py`. Two claims, both tested rather than asserted.

**(1) Active off-pin conditions = `2(n − |S|)`** (and `2(n−1)` at `|S| = 0`, since there are
only `n−1` blocks). Checked for every subset at `n = 8` and `n = 12`, 2 curves × 2 tree
shapes: **0 mismatches**. Measured, `n = 12`:

```
|S|      : 0  1  2  3  4  5  6  7  8  9 10 11 12
#off-pins:22 22 20 18 16 14 12 10  8  6  4  2  0
```

**It decreases with `|S|`.** A constraint family that thins out as `|S|` grows cannot bound
`|S|` from above; if it constrained anything at all it would be a lower-bound mechanism.

**(2) It constrains nothing.** On a dead block the condition is `c'_j·P | a'_j·i_j`, i.e.
`P | i_j` (true value 0, the identity slot) and `m'_j | (i_j/P)`. The only consumer of
`(i5,i6)` is the mux multiplier `cC = ab = L`, which is **0 exactly when the off-pin fires**.
So `i_j/P` occurs in no other condition and `w = 0` satisfies it.

Checked adversarially rather than argued: the private outputs of every dead block were driven
to random nonzero multiples of `P` and **every other condition in the system re-evaluated**.
**11 987 condition values compared; 0 moved.**

> **The off-pin family closes for every subset at every `|S|`. It carries no signal.
> The whole of §8's hope lives in the congruence family.**

*(Scope: this is a statement about a full assignment. W's observation that the 39 026 witness
breaks two off-pins of a dead block is about a partial assignment whose gates are not at
their true values; nothing here contradicts it.)*

---

# 5. THE ABSORBING-BLOCK THEOREM — why the congruence family cannot make a cutoff

**Theorem.** At a live merge block let `A = i1 − i2` (the integer difference of the merged
x-slots) and let `c` be the block's small modulus. If `gcd(A,c) = 1` then for **any** values
of the four input lifts, the block's own output lifts `(u5,u6)` can be chosen so that all
three integer-lift conditions hold. The block then constrains nothing.

*Proof.* `nu1 = nu1a + u5·A²`, `nu2 = nu2a + u5·B + u6·A` (V2c), so `(u5,u6) ↦ (nu1,nu2)` is
affine with matrix `[[A²,0],[B,A]]`, determinant `A³`, invertible mod `c`. Take the preimage
of `(0,0)`; every condition `m_k | c_k1·nu1 + c_k2·nu2` is homogeneous, so all three hold. ∎

Machine-checked (T1): 2 515 cases with `gcd(A,c)=1`, real curve points, random `3×2` matrices,
random input lifts, brute-force search over `(u5,u6)` — **0 failures**.

**Corollary, and it is the whole of §8 in one line.** `A = a + P·δ` where `δ = t_R − t_L` is a
free lift difference and `P` is invertible modulo any small `c`. So `A mod c` is **not fixed by
the curve arithmetic** — the input lifts alone can usually make it invertible. A block can
only bite on a coincidence, the coincidence has probability `≈ 1/ℓ` per prime `ℓ` dividing its
modulus, and the resulting closure rate is `≈ ρ^{|S|−1}` with `ρ` a per-block constant.
**There is no mechanism here that can produce a threshold in `|S|`.**

How weak the binding is, measured: the "all live blocks have `gcd(a,c)=1` in reduced
coordinates" rate is a *large under-estimate* of the closure rate, precisely because the lifts
repair `A` (T3, `n = 12`, stress):

```
 w                    : 2      3      4      5      6      7
 all-absorbing (proxy): .3788  .0909  .0242  .0025  .0000  .0000
 actual closure rate  : .7727  .6045  .5051  .4192  .3528  .2980
```

---

# 6. SCALING, AND THE EXTRAPOLATION TO n = 256

**Model, stated before fitting:** `#closing subsets of weight w  ≈  C(n,w)·ρ^{w−1}`, hence
`Bmax(n) = max{ w : C(n,w)·ρ^{w−1} ≥ 1 }`. (`ρ` is fitted as the slope of `log(rate)` vs `w`,
independently per draw.)

**Model validated on the measured points before being used** (`ad_report.py` §1): over
**340 draws** spanning every configuration, curve, size, tree shape and modulus pool,

```
model-Bmax - measured-Bmax = -1 :  46 draws (13.5%)
model-Bmax - measured-Bmax =  0 : 252 draws (74.1%)
model-Bmax - measured-Bmax = +1 :  40 draws (11.8%)
model-Bmax - measured-Bmax = +2 :   2 draws ( 0.6%)
|error| <= 1                    : 99.4% of 340 draws
```

That is the licence to extrapolate, and it is the only licence claimed.

**`ρ` by size** — the question that decides whether the extrapolation is safe. Two fits: the
raw one (every weight with rate in `(0,1)`) and a **robust** one that requires a weight to
carry `≥ 20` subsets and `≥ 5` closing ones and at least 4 such weights to survive — the raw
fit at `n = 8` can rest on five points with single-digit counts, and that is exactly where the
outliers live.

| `n` | draws | raw `ρ` min / med / max | **robust `ρ`** min / med / max | `n − Bmax` |
|---|---|---|---|---|
| 8 | 93 | 0.2657 / 0.7756 / 1.0612 | **0.5450** / 0.8281 / 1.0612 | 0–4 |
| 10 | 53 | 0.7797 / 0.9710 / 1.0610 | **0.7797** / 0.9378 / 1.0610 | 0–2 |
| 12 | 100 | 0.6350 / 0.8024 / 1.0139 | **0.6493** / 0.8173 / 1.0139 | 0–3 |
| 14 | 55 | 0.5494 / 0.7479 / 1.0180 | **0.5707** / 0.7723 / 1.0180 | 0–4 |
| 16 | 39 | 0.6688 / 0.7980 / 1.0089 | **0.6688** / 0.8362 / 1.0089 | 0–4 |

**`ρ` is flat in `n`** over five sizes and 340 draws — it is a per-block quantity, as §5 says
it must be. `n − Bmax` stays in `0..4` while `n` doubles; it is not `c·n`, not `n/2`, not `n − c` with
`c` growing, and nothing that reaches `n − 200`.

**Extrapolation** (`ad_report.py` §4), taking the most constraining `ρ` observed *anywhere*
and holding it fixed — conservative, since the median is far higher:

| `ρ` | implied bound at `n = 256` |
|---|---|
| **0.5450** (most constraining ROBUST fit anywhere) | **`w ≤ 205`** |
| 0.5762 (5th percentile, raw) | `w ≤ 211` |
| 0.7980 (median, raw) | `w ≤ 241` |
| ≈1.00 (typical at the real instance's `frac = 0.25`) | `w ≤ 256`, i.e. nothing |
| *0.2657* (single raw outlier: one `n = 8` draw, 5-point fit, counts down to 1) | *`w ≤ 132`* |

**The outlier is reported, not hidden, and it does not change the verdict.** Even if that
5-point `n = 8` fit were taken at face value and held to `n = 256`, the implied bound is
`w ≤ 132` — half a σ below the null mean, excluding essentially none of the null mass, and
still `2.4×` above the `w ≲ 56` needed for §8 to beat rho. Its robust re-fit is 0.66.

The free unconditional bound is `w ≤ 255` (AB §7); the null puts `w ∈ [104,152]` with
probability 0.998. **`w ≤ 205` is 9.6σ above the null mean — it excludes a region of null mass ≈ `2^−100`,
i.e. it is vacuous, and it is worse than agent Y's complement mechanism at `W = 50`.**

> **How much three-to-five points can support — stated the way AB had to state `d_reg`.**
> What is *measured* is: `ρ` is a per-live-block constant, flat over `n = 8..16`, and the
> `C(n,w)·ρ^{w−1}` model predicts `Bmax` exactly at every measured `n`. What is *extrapolated*
> is only that `ρ` stays flat out to `n = 256`. That extrapolation is supported by a proof
> (§5: `ρ` is a per-block probability, with no `n` in it) and by five measured sizes, which is
> a better footing than a bare three-point fit — **but it is still an extrapolation, and the
> conclusion "§8 is dead" rests on it.** What would overturn it is a mechanism making `ρ`
> shrink with `n`; §5 says there is none, because a block sees only its own two children.

---

# 7. WHAT WOULD HAVE TO BE TRUE — and why it is self-defeating

`ad_report.py` §4 inverts the model: the per-live-block closure probability needed at
`n = 256` for a bound of `B`, and what that same `ρ` costs the instance's own satisfiability
(the true subset, of weight `B`, closes with probability `ρ^{B−1}`):

| target bound | `ρ` needed | probability the true weight-`B` subset itself closes |
|---|---|---|
| `w ≤ 198` (AB's break-even) | ≤ 0.511 | `2^−190.8` |
| `w ≤ 148` | ≤ 0.315 | `2^−245.0` |
| `w ≤ 106` (AB's crossover) | ≤ 0.199 | — |
| **`w ≤ 56`** (beats rho) | **≤ 0.093** | `2^−188.3` |
| **`w ≤ 24`** (actionable) | **≤ 0.037** | `2^−109.8` |

The smallest `ρ` observed anywhere in this study, across every curve, size, tree shape,
coefficient pool and draw — including settings **four times more constraining than the real
instance** (`frac = 1.0` vs `927/3707 = 0.25`) — is **0.658** (robust fit; **0.266** for a
single small-sample `n = 8` raw fit, which even taken at face value only reaches `w ≤ 132`);
the sweep over adversarial modulus pools (`{4}`, `{8}`, `{2,4,8}`, `{3,9}`,
`{2,3,5,7,11,13}`, all at `frac = 1.0`) pushes the robust minimum to **0.545**, which implies
`w ≤ 205` at `n = 256`.  Note which pools do it: **concentrated** small prime powers
(`{8}`, `{3,9}`) are the worst case; a *wide* pool `{2,3,5,7,11,13}` is much weaker
(`ρ ≈ 0.78–1.00`), because each condition then draws only one modulus and a block is rarely
bound by all of them.

**And there is a floor, from the theorem in §5.** A block is absorbing whenever `A` is
invertible modulo its own modulus `M_b`, so

```
    rho_b  >=  prod_{l | M_b} (1 - 1/l)
```

— and this is only a lower bound, because the input lifts can additionally repair a bad `A`.
By Mertens, `∏_{ℓ≤X}(1−1/ℓ) ≈ e^{−γ}/ln X = 0.5615/ln X`, so `ρ ≤ 0.093` requires the prime
support of **every block's** modulus to reach `X ≈ 420`. Three moduli per block cannot span
all primes below 420 unless each has magnitude `≈ 2^170`. **That is a checkable criterion on
the real instance, and it is nowhere near met:** W's census reports ordinary small `|c| > 1`
on 288 of the 1 149 congruences. For pure prime-power moduli it is starker still — with
`M = 2^k`, a block binds only when `A` is even, so `ρ ≥ 1/2` however large `k` is. That is
exactly the ≈0.55 floor measured at pool `{8}`.

And the right-hand column is the sting: **an instance tight enough to bound `w` this way is an
instance that is satisfiable only by design, not by accident.** The mechanism eats itself.

---

# 8. THE THREE QUESTIONS, ANSWERED

**Q. Do the lift conditions close for the true `S` only, or for many `S`?**
For many. **The lift conditions do not mention `T` at all** — they are per-block conditions on
partial-sum coordinates, and `T` enters the system only through the root pin. At the real
instance's `frac = 0.25`, most coefficient draws close **all `2^n` subsets**. Closure is
therefore not a filter on candidate solutions in any useful sense; it is not even correlated
with being a solution.

**Q. Is the set of `S` for which the lift closes correlated with `|S|` at all?**
Weakly, through one channel only: the number of live blocks, which is `|S|−1` exactly. The
full curves are in §3. Under stress the rate falls geometrically at `≈0.93` per live block; at
the real instance's parameter the curve is U-shaped with the **minimum in the middle** and
rate `1.0000` at `|S| = n`. In neither regime does the top of the range empty out faster than
`C(n,w)` fills it.

**Q. Same question for the off-pins?**
No correlation, and no constraint: §4. The active-condition count runs the other way
(`2(n−|S|)`), and the conditions are on variables that occur nowhere else.

---

# 9. SCOPE — what is NOT claimed

* This is a statement about **the small analogue**, built from the decoded law shapes. It is
  faithful in: the law (V1, V2), the mux/gate alignment, the integer-lift discipline, the
  live-block count `|S|−1`, prime group order with `n = bitlen(N)` and an `n`-bit ladder,
  `j = 0` presence, and both tree shapes including the real 178:78 root split. It does **not**
  reproduce the real instance's specific coefficients, which are not knowable to me under
  rule 2 — that is why every measurement is over a *distribution* of coefficient draws and why
  §7 states the threshold the real coefficients would have to meet.
* The real system has 383 blocks for 255 merges (101 pass-throughs, 27 dead) and 3 707 lift
  conditions against `5 × 383 = 1 915` block atoms, so there are lift conditions I have not
  modelled (leaf pins, aliases). The obvious one — leaf slots carrying their own lift modulus
  rather than being pinned exactly — is implemented (`leaf_free` in `ad_model.LiftDP`) and
  **measured, not assumed**: `ad_measure.py leaffree` gives closure rate `1.0000` at **every**
  weight, `Bmax = n` in 6/6 draws, in every configuration. Extra lift freedom only pushes `ρ`
  up, i.e. it can only strengthen the negative verdict.
* **Equation-level cancellation is not modelled.** W's boundary applies here too: everything
  above classifies solutions of `atoms = 0`, and the real checker requires *equations* — sums
  of ~12 atoms — to vanish, a strictly larger solution set. Larger again means more closure,
  again in the direction of the negative verdict.
* Nothing here is an infeasibility claim about the real instance, and nothing here bounds `w`.
  It says only that **this particular mechanism cannot bound `w`**.

---

# 10. FILES

| file | what |
|---|---|
| `ad_curves.py`, `ad_curves2.py` | small curves, prime order, `n` bits; `j=0` and generic |
| `ad_curves_partial.json`, `ad_curves2.json` | the curves used, each independently re-verified |
| `ad_model.py` | the analogue: law, tree, mux, lift algebra, exact CRT tree DP |
| `ad_validate.py` | V1–V5 (run this first; it exits nonzero on any failure) |
| `ad_v3b.py` | V3b: the rate-vs-`\|S\|` curve re-derived from raw integers at `n=7` |
| `ad_theorem.py` | T1 absorbing-block theorem, T2 greedy certifier, T3 the proxy bound |
| `ad_offpin.py` | the off-pin family, count law and adversarial independence check |
| `ad_measure.py` | closure rate vs `\|S\|`, exhaustive; modes `quick main main16 stress perprime leaffree n20` |
| `ad_scale.py` | the scaling study across `n`, with the model prediction alongside |
| `ad_report.py` | model validation, `ρ` by size, `j=0` vs generic, extrapolation to 256 |
| `ad_*.log`, `ad_*.json` | raw outputs |

Reproduce: `PYTHONDONTWRITEBYTECODE=1 AD_CURVES=ad_curves_partial.json python3 ad_validate.py`
then `ad_theorem.py`, `ad_offpin.py`, `ad_measure.py stress`, `ad_scale.py 2:1.0`, `ad_report.py`.
