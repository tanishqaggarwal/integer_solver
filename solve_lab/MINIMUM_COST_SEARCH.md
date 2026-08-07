# Minimum-cost search to solve the instance

Coordinator synthesis over the whole campaign. Every input is a measured result recorded in
`FLEET.md` or an agent's RESUME file; the arithmetic below was computed, not estimated.

## 1. What the problem reduces to

Established by **four independent parses sharing no code** (agents F, K, L, Q, corroborated by P,
R, T, U, W):

> The instance is satisfiable **iff** some subset `S ⊆ {0..255}` of leaf selectors satisfies
> `k·G = T`, where `k = Σ_{i∈S} 2^i`.

Supporting facts, each verified at instance level:

- **The lift is not an obstacle.** Agent T closed the 927 integer conditions at
  `|S| = 2, 3, 5, 6, 7, 8, 17`, each dumped and checked. So **finding `S` yields 39,033/39,033**;
  the remaining work after `k` is minutes.
- **The classification closes at exactly two families** (agent W): the Jacobian
  `d(N1,N2)/d(i5,i6)` has `det = A³`, so `A ≡ 0` forces `B ≡ 0` (degeneracy, output free) and
  `A ≢ 0` determines the output uniquely (chord). **No third case; the mixed case is impossible.**
  *Scope: this is at atom level. Equation-level cancellation is a larger set the theorem does not
  cover.*
- **The degeneracy route is closed** (agent U, exhaustive): all 510 proper slot supports have masked
  value strictly below `N`, largest ≈ 0.7987·N, so no two gadget inputs can coincide by
  configuration.
- **The solution is essentially unique.** `2²⁵⁶ − N ≈ 2¹²⁸`, so there are one or two valid `S`, and
  two only with probability ≈ 2⁻¹²⁸. No multi-target advantage exists.

## 2. The curve is secp256k1 — verified bit-for-bit

    p == 2^256 - 2^32 - 977                      : True
    N == 0xFFFF...FFFEBAAEDCE6AF48A03BBFD25E8CD0364141 : True

This is why every structural shortcut checked in this campaign failed, and it is not a coincidence
of this instance: prime order, **not** anomalous (`N ≠ p`), embedding degree past 24 (verified to
k ≤ 24), and the only endomorphism is the standard order-3 `λ`. Fifteen years of public
cryptanalysis have produced no better-than-generic attack on it.

## 3. The minimum-cost known attack

**Parallel Pollard rho with distinguished points, negation map, and the GLV endomorphism.**

| variant | cost |
|---|---|
| plain rho, `√(πN/4)` | 2^127.8 |
| + negation map (`√2`) | 2^127.3 |
| + GLV endomorphism (`√3`) | 2^127.0 |
| **+ both (`√6`)** | **2^126.5 group operations** |

Engineering specification, so this is a real answer rather than a shrug:

- **Iteration**: r-adding walk, r ≈ 20 partitions, walking on equivalence classes under the
  order-6 automorphism group `{±1, ±λ, ±λ²}`.
- **Fruitless cycles**: the negation map's known pathology. Detect by cycle-length check and escape
  by a deterministic perturbation; without this the `√2` is lost.
- **Distinguished points**: `x`-coordinate with `d` leading zero bits. Storage ≈ `2^126.5 / 2^d`
  points; per-worker memory `O(1)`.
- **Parallelism**: linear speedup in workers (van Oorschot–Wiener), no communication except DP
  reports.
- **Termination**: a DP collision yields `k` directly; then run agent T's `t_close2wj.py` lift and
  verify with `checker.py`.

### What it costs in the real world

| machine | wall time for 2^126.5 |
|---|---|
| one CPU core (10⁶ ops/s) | 3.9×10²⁴ years |
| one GPU (10⁹ ops/s) | 3.9×10²¹ years |
| one million GPUs (10¹⁵ ops/s) | 3.9×10¹⁵ years |
| **Bitcoin-network scale (6×10¹⁸ ops/s)** | **6.5×10¹¹ years ≈ 47× the age of the universe** |

**This is not a resource question.** No budget reaches it.

## 4. The only thing that beats rho: bounded-complexity `k`

If `k` is drawn from a small structured class, meet-in-the-middle over that class beats `2^126.5`.

**Unsigned Hamming weight** — `k` has `w` one-bits. MITM cost ≈ `C(256, w/2)`:

| bound | cost | vs rho |
|---|---|---|
| w ≤ 10 | 2^33.0 | beats |
| w ≤ 20 | 2^58.0 | beats |
| w ≤ 30 | 2^79.1 | beats |
| w ≤ 40 | 2^97.8 | beats |
| w ≤ 50 | 2^114.6 | beats |
| **w ≤ 56** | **2^123.8** | **crossover** |
| w ≤ 60 | 2^129.7 | loses |

**Signed-digit weight** — `k = Σ ε_j·2^{e_j}`, `ε ∈ {±1}`, `m` terms. **Strictly more general**: it
contains every low-Hamming-weight `k`, *and* every low-**run-length** `k` (a run of ones is
`2^{a+1} − 2^b`, two signed terms), *and* anything from a short addition-subtraction chain. Cost
≈ `C(256, m/2)·2^{m/2}`:

| bound | cost |
|---|---|
| m ≤ 8 | 2^31.4 |
| m ≤ 10 | 2^38.0 |
| **m ≤ 12** | **2^44.4** |
| m ≤ 14 | 2^50.6 |
| m ≤ 16 | 2^56.5 |

**Signed-digit `m ≤ 10` costs about what unsigned `w ≤ 12` costs and covers a strictly larger,
structurally different set.** Nobody in this campaign has searched it. Agent X has been tasked with
it after its plain-weight sweep.

### Budget → reachable bound

| budget | unsigned weight | signed-digit weight |
|---|---|---|
| 2^47 (this box, ~8 h) | w ≤ 14 | m ≤ 12 |
| 2^60 | w ≤ 20 | m ≤ 16 |
| 2^70 | w ≤ 24 | m ≤ 20 |
| 2^80 | w ≤ 30 | m ≤ 24 |
| 2^90 | w ≤ 34 | m ≤ 28 |
| 2^126.5 | — | **full solve by rho** |

## 5. Already excluded (all mod p; see the standing note on Q's withdrawal and restoration)

`w ≤ 6` exhausted; `w ≤ 7` partially; `k` not in the bottom or top `2^44`; ON-bits not confined to
any 34-bit window; `k` not a small multiple of a ladder point (`m ≤ 10⁷`); not in the endomorphism
orbit; no `k = a + bλ` with `|a|,|b| < 2^21`.

## 6. The bottom line

**If the instance was generated by choosing a random 256-bit `k` and publishing `T = k·G`, it is not
solvable by any known method, and no search this campaign or any successor can run will change
that.** The reduction is complete, the classification is closed, and the residual is a generic
ECDLP on the most-studied curve in existence.

**The minimum-cost search that can actually be run** is therefore bounded-complexity MITM —
signed-digit weight first, since it dominates unsigned at equal cost. It wins **only** if `k` lies
in a low-complexity class, which requires a prior the mathematics cannot supply.

**The one place such a prior could come from is how the instance was constructed** — the
leaf-to-exponent assignment, the ordering of the 512 pin constants, anything reflecting a generator
rather than the mathematics. **That direction was closed by user instruction at the start of the
campaign and the fleet has respected it throughout.** It is the only remaining lever, and it is the
user's to open.

---

## 7. Can `w` be bounded? — the objective is flat

**Measured, not argued.** Agent T closed the integer lift at `|S| = 2, 3, 5, 6, 7, 8, 17` and every
one scores **exactly 39,018 with the identical 15-equation failing set** — precisely the footprint of
the two target congruences, the only equations in the instance that depend on `k`.

> **The objective is FLAT: every correctly-lifted configuration scores 39,018 regardless of its
> weight, and 39,033 only at the exact answer. A cliff, not a gradient.**

Consequences: **no measurement of the instance leaks any information about `w`**; and no local
search, annealing, or gradient method can work, which is why every search in this campaign
plateaued. Seven independent configurations, same score, same failing lines.

### The null

If `k` is uniform, `w ~ Binomial(256, ½)`: **mean 128, σ = 8**, so `w ∈ [104,152]` at 99.7% and
`[120,136]` at 68%. `C(256,52) ≈ 2^117` — that range is unreachable.

### What exhaustion buys, against a structured hypothesis

`H_A`: the designer picked `w` uniform on `[1, Wmax]`. Missing at `W₀` gives likelihood ratio:

| Wmax | W₀=7 | W₀=10 | W₀=12 | W₀=16 | W₀=20 |
|---|---|---|---|---|---|
| 12 | 2.4× | 6.0× | **killed** | — | — |
| 20 | 1.5× | 2.0× | 2.5× | 5.0× | **killed** |
| 30 | 1.3× | 1.5× | 1.7× | 2.1× | 3.0× |
| 60 | 1.1× | 1.2× | 1.3× | 1.4× | 1.5× |

**Exhaustion is only informative when it approaches `Wmax`.** Finishing `w ≤ 12` moves a `Wmax = 30`
hypothesis by 1.7×. **The value is all in crossing the designer's bound, where the hypothesis dies
outright rather than being nudged.**

### The decisive threshold

Costs: `w ≤ 12` → 2^38 · **`w ≤ 20` → 2^58** · `w ≤ 24` → 2^67 · `w ≤ 30` → 2^79.

**2^58 is the number to aim at** — roughly 3 days on a thousand-GPU cluster; out of reach on this
box, which does ~2^47 in eight hours. It is decisive both ways: it **finds `k`**, or it **refutes
every "small designed weight" hypothesis up to 20 outright.** Signed-digit `m ≤ 20` costs 2^68 and
covers strictly more (low weight, low run-length, short addition-subtraction chains).

### Summary

- **Hard lower bound:** `w ≥ 8` once the current level completes. The only rigorous statement
  available; it grows about one level per 4× budget.
- **Upper bound:** none exists and none is obtainable from the instance.
- **The one real inference available:** clearing 2^58 licenses *"if this was designed with a
  low-weight key, the weight exceeds 20"* — a bound on the designer's choice, not on the mathematics.

**Everything below `w ≈ 20` is cheap enough to be worth doing and too weak to conclude much from.**

---

# CORRECTION (check-in 103) — the cost model in §4 was wrong

Agent AB, auditing its own model before leaning on it further, found that a radius-`W` Hamming ball
was priced at `C(256, W/2)`. **The correct per-ball cost is the cumulative half-volume
`√W · Vol₁₂₈(W/2)`** (Coppersmith / Stinson splitting systems). Wrong by up to **2^65** at large `W`.

The sanity check that catches it: **at `W = 256` one ball is the whole key space, so the cost must be
2^128. The old model returned 2^251.7.**

**~~Break-even at B = 198~~ → break-even at `B = 148`. Do not circulate 198.**

**§4's tables are pessimistic.** Time-only, unbounded memory:

| budget | ~~old~~ | corrected |
|---|---|---|
| 2^47 | ~~w ≤ 14~~ | **w ≤ 18** |
| 2^58 | ~~w ≤ 20~~ | **w ≤ 24** |
| 2^80 | ~~w ≤ 30~~ | **w ≤ 40** |
| rho crossover | ~~w ≈ 56~~ | **w ≈ 104** |

**MEMORY BINDS, and the corrected figures describe a machine nobody has.** At ~2^30 entries on this
box the real cap is **`w ≤ 10`**. A memory-aware costing has never been done here; AB is doing it,
and the deliverable is achievable `w` as a function of **(time, memory)**, not time alone.

**The qualitative conclusion survives and is cleaner.** `w ≤ 148` is +2.6σ on the null, excludes
~0.4% of null mass (~0.006 bits), and costs **2^126.4 against 2^126.5 to solve outright and learn
`w` exactly.** The covering optimum degenerates toward *"one ball = the whole space = solve it"*, and
the gap to AB's independently derived generic lower bound is now **≤2^3 everywhere**, against a
spurious 2^57 before. **A corrected model agreeing with an independent bound is stronger evidence
than the original was.**

## Two barriers now stand where §4 previously had only an absence

**Theorem C.** The only weight-preserving affine self-map of `Z_N` is the identity — forced by
`popcount(2^256 mod N) = 65 ≠ 1`, verified across all 255 `j`. In the generic model affine is all an
algorithm can realise, so **no weight-preserving randomised self-reduction exists**, and the
hardcore-bit machinery has nothing to run on.

**Theorem D.** Generic group model: every held element is `σ(α_i + β_i k)`, so a collision is one
affine equation over the **field** `Z_N` and has one root. For **any** predicate `P` and `m` queries,
`Adv ≤ m² / (2·min(|D₀|,|D₁|))`. For `P = [w ≤ B]`:

- **`B = 128` ⇒ `m ≥ 2^127.5`**, against 2^126.5 to solve. **Deciding the weight predicate is, to
  within the automorphism speedup, exactly as hard as solving the instance.**
- **`B = 20` ⇒ `m ≥ 2^49`**, and corrected MITM achieves **2^50.0 — within 2^1 of optimal.**

Knobs: generic model (coordinate encoding excluded), average-case over a distribution on `k`.
**HGJ / BCJ representation techniques are the one route that could beat the ball cost, and they need
`k₀ mod M` — the same obstruction as §3.**

## The instance side is closed by measurement (check-in 102, agent Z)

§7's asserted *"no upper bound is obtainable from the instance"* is now a **measured** result:
of **819,975 monomials, 0 contain two distinct selectors**; booleanity-reduced affine elimination
over the complete instance leaves **3,980 rows, all identically `0 = 0`, zero genuine linear
constraints on the selectors**, run both mod `2^61−1` and exactly over ℚ; and **0 adder-shaped atoms
among 9,527 all-boolean atoms — `Σ s_i` is never formed anywhere.**

**So the Hamming-weight angle is closed from both ends: bounds must be assumed, never derived.**

---

# AMENDMENT (check-in 105) — agent Z's audit of the corrected model

**Both AB results survive; neither survives unamended.** Independently re-derived in exact integers.

**CONFIRMED:** the corrected covering cost `poly(W)·Vol₁₂₈(⌈W/2⌉)` is right — plus a **third**
disqualifier for the round-1 model AB did not state: **`C(256,W/2)` is not monotone in `W`**, and a
ball cost must be. **Theorem B's qualitative claim confirmed independently**: scanning all
`W ∈ [0,256]`, the optimum is the degenerate one-ball cover for every `B` below ≈247. **The table is
confirmed exactly and filled in:**

| budget | 2^30 | 2^40 | 2^47 | 2^58 | 2^70 | 2^80 | 2^90 | 2^126 |
|---|---|---|---|---|---|---|---|---|
| reachable `w` (time-only) | 10 | 14 | **18** | **24** | 32 | **40** | 48 | 104 |

**AMENDED:**

1. **The `W = 256` certificate fails by 2^4** — required 2^128, AB's corrected model returns
   **2^132.0**. `√W` is spurious exactly there (`W/2 = 128` equals the half size, so every split is
   already balanced); without the poly factor it is **2^128.0 exactly**. **This is the very test AB
   used to certify its own correction.**
2. **Latent floor/ceil bug**: `V128[W//2]` should be `⌈W/2⌉`. **Odd radii underpriced by up to 4.6
   bits** (`W=9`: 2^25.0 vs 2^29.6). AB's published numbers are unaffected (even-only scans, floor =
   ceil) **but it fires on first reuse** — Z's own first pass produced spurious optima at
   `W = 9, 15, 35, 55, 75, 103, 127` and a spurious `+1` on every budget row. **Recorded so the next
   agent does not rediscover it as a finding.**
3. **Two numbers reconciled:** the largest affordable complement radius and the rho crossover are the
   **same function at the same budget** (AB's RESUME says 104, its script 107). **Reproducible:
   crossover `w = 106`, break-even `B = 149`** (AB says 148).
4. **The memory caveat is the binding one and is understated:** at 2^47 the half-list is **2^44.2
   entries**, so on this box (~2^30) the real reach is **`w ≤ 10`, not 18.** The time-only row
   describes a machine nobody has.

**PROVED, not spot-checked:** the MITM never dips below the generic floor —
`(Vol₁₂₈(w/2))² ≤ Vol₂₅₆(w)` by restricted Vandermonde, exact for every even `w`. **So MITM is
optimal to within `√W ≤ 2^4`**, and AB's "within 2^1 at `B = 20`" is exact at 1.50. **There is no
room left in the algorithm; class size is the whole story.**

## Theorem D — amended, and one claim STRUCK

- **The normaliser `min(|D₀|,|D₁|)` is right; the constant is not.**
  `Adv ≤ |Bad|·(1/|D₀| + 1/|D₁|) ≤ m²/min(|D₀|,|D₁|)`. **AB's `m²/(2·min)` drops the second side, is
  a factor 2 tighter than the argument supports, and the error direction OVERSTATES the barrier** by
  0.5 bits.
- **The single-root argument survives the automorphism group** (`λ` and negation are multiplication
  by fixed scalars, so held elements stay `σ(α′+β′k)`; one affine equation, one root over prime `N`)
  — **but not the encoding**: under `x`-coordinate + GLV the order-6 orbit collapses, giving `AUT`
  equations per pair and degrading the bound by `√6 = 1.29` bits.
- **`B = 128 ⇒ m ≥ 2^127.5` does follow**, recomputed by exact digit-DP over `k ∈ [0,N)`
  (self-checked `#{k<N} = N`): `min = 2^254.93`, `√min = 2^127.46`.

> ## ⚠ STRUCK: "deciding the weight predicate is HARDER than solving"
>
> **Taken literally, `m ≥ 2^127.5` against 2^126.5 to solve says deciding is harder than solving —
> which is impossible, since any solver decides.** The first figure **excludes** the automorphisms;
> the second **includes** them. **In one model, with the corrected constant and `AUT = 6`, the bound
> is 2^125.7 ≤ 2^126.5** and the inequality points the right way.
>
> **The qualitative conclusion is unharmed: there is no generic shortcut for the weight predicate,
> and deciding costs the same as solving to within the same `√6`.** Only the strict phrasing is
> withdrawn — from this file, from `FLEET.md`, and from the report to the user.

---

# SECOND AMENDMENT (check-in 106) — the crossover correction is WITHDRAWN

Agent AB, auditing its own round-2 correction, withdrew it. **Two self-retractions in three rounds,
both found by re-checking a model before leaning on it further.**

> ## ⚠ WITHDRAWN: ~~"the crossover moves `w ≈ 56 → 104`, doubling §8's payoff band"~~
>
> **An unbounded-memory artefact. Memory enters as a CUBE ROOT** — `2^30 → 2^60` moves the crossover
> by only 12 — **so the realistic crossover is `w ≈ 52–64`. The campaign's ORIGINAL `w ≈ 56` stands,
> and §8's payoff band is NOT doubled.**

**The reconciliation is the part worth keeping:** *the original table's time model was pessimistic
and it ignored memory, and the two errors cancelled.* **This box's reach is `w ≤ 14` — exactly what
this document said before any of the corrections.** A number derived twice from opposite errors is
worth more than one never challenged.

## The reference table is now memory-aware (vOW golden collision)

`T = rep·c·L^1.5/√M`, so reach `L ≤ (T²M/(rep·c)²)^{1/3}`. vOW dominates chunked rescanning for all
`M < L`.

| time \ memory | **2^30 (this box)** | 2^40 | 2^50 | unbounded |
|---|---|---|---|---|
| 2^47 | **14** | 16 | 18 | ~~18~~ |
| 2^58 | 18 | 20 | 22 | ~~24~~ |
| 2^80 | 26 | 30 | 32 | ~~40~~ |
| 2^126.5 | **52** | 56 | 60 | ~~106~~ |

**The unbounded-memory column is unreachable and is struck.** **Disk is not a way out:** random
access costs ~2^20 slowdown for 2^1 of memory.

## The `W = 256` certificate now passes exactly

`√W` was an asymptotic standing in for the exact hypergeometric reciprocal

    rep(W) = C(256,128) / ( C(W,⌈W/2⌉) · C(256−W, 128−⌈W/2⌉) )

which tends to `√(πW/2)` for small `W` (`rep(10) = 3.98` vs 3.96) but is **`rep(256) = 1.0000`
exactly**, since at full radius every split is already balanced. With floor→ceil fixed and
monotonicity repaired by a **suffix minimum** — principled, since a radius-`W` ball is searchable by
any `W' ≥ W` procedure — the certificate is **re-certified at 2^128.0000**. Crossover **106**
(matches Z); break-even **148** (Z: 149, a rounding difference).

## `d_reg` — the WEAK form, which is what the measurement supports

Measured by Macaulay/XL over `GF(q)` in python-flint: **`n=2 → 4`, `n=3 → 5`, `n=4 > 4`** (unfinished).
**Two points cannot support extrapolation to `n = 256`, and AB said so rather than quoting its
≈258.** What they establish: **`d_reg` increases with `n` at all — which was the only way §9.12
could have been wrong.** The verdict hardens on that basis alone.

**Two model fixes AB needed first, either of which alone would have produced a confident wrong
verdict:** the initial model was **leaky** (when `R_j = P_j` the slope is unconstrained, giving a
spurious positive-dimensional component — why it reported max-GB-degree 2 and failed to pin
selectors), and **reduced-GB degree is the wrong statistic**, being 1 for any unique-solution ideal;
the right one is **solving degree**.

## Tooling correction — the coordinator's error

**Singular is NOT installed.** No binary on the filesystem, nothing in dpkg. **Only `sympy 1.14.0`
and `python-flint 0.9.0` are present; no Sage, no msolve, no Macaulay2, no Magma, no PARI/gp.**
The earlier inventory naming "Singular 4.3.2" came from one agent and **was broadcast by the
coordinator without verification.** AB worked around it by building Macaulay/XL itself.
