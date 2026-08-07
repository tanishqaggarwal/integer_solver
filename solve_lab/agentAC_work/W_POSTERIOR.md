# W_POSTERIOR — the exact probabilistic bound on `w`

**Agent AC.** Everything below was computed in this directory. Scripts:
`ac_prior.py`, `ac_cond.py`, `ac_post.py`, `ac_extra.py`, `ac_point.py`, `ac_point2.py`.
Outputs kept as `ac_post.out`, `ac_extra.out`, `ac_point.out`, `ac_point2.out`,
`ac_thirdcheck.out`, `ac_table.out`. Total footprint 200 KB. No git commands were run;
no file outside this directory was written, renamed or deleted.

**What is mine and what is theirs.** The combinatorics, the digit-DPs, the quantiles, the tail
exponents and every test on `T` in §5 are mine, computed here from scratch and each verified by a
second independent route. Two *exhaustion* facts are taken on the record's authority and are
named inside every claim that uses them: **X's** "unsigned weight ≤ 9 exhausted, 0 hits" and
**Y's** "complement weight ≤ 9 exhausted, `ALLDONE`, per-`i0` lines summing to exactly
`C(256,5)`". I did not re-run either sweep — they are 2^33 candidates each and re-running them was
not my task. `p`, `N`, `G`, `T` are from `agentX_work/xdata.json` and `agentY_work/ydata.json`,
which are independent parses and agree bit-for-bit; every *structural* property of them
(`p` prime, `N` prime, `T` on curve, `N·T = O`, `a = 0`, `b/7` a sixth power) is re-verified here.

---

# AUTHORITATIVE SUMMARY

**Under the model `k₀ ~ uniform[0,N)`, conditioned on everything the fleet has established
(`10 ≤ w ≤ 246`):**

| quantity | exact value |
|---|---|
| `P(w ≤ 14)` — reach of this box (AB §S3) | **`2^-180.780`** = 3.7998 × 10⁻⁵⁵ |
| `P(w ≤ 24)` — "actionable" band (AB §S5) | **`2^-144.487`** = 3.1995 × 10⁻⁴⁴ |
| `P(w ≤ 56)` — beats rho at 2^40 memory (AB §S5) | **`2^-65.570`** = 1.8254 × 10⁻²⁰ |
| 90 % interval | `w ∈ [115, 141]` |
| 1 − 2^-80 interval | `w ∈ [49, 207]` |
| total movement of the distribution by the entire campaign | TV = **`2^-201.623`**; **2.9 × 10⁻⁶¹ bits** of information |

> **The one-sentence conclusion, unsoftened: under the uniformity model the probability that `w` is
> small enough for any of AB's mechanisms to help is `2^-65.6` at the most generous band (`w ≤ 56`,
> which only ties rho) and `2^-144.5` at the band that would actually be actionable (`w ≤ 24`) —
> so the honest answer is that the probability `w` is small enough to help is astronomically small,
> and no amount of further search of the kind this campaign has been running will change that,
> because the entire campaign so far has moved the distribution of `w` by `2^-201.6`.**

**The corollary that is actually useful, and it is not a number about `w`.** Under uniformity the
searches buy nothing. They buy something only under a *designer* hypothesis, and there the leverage
is computable and not small: **if one models "the instance's scalar was chosen with weight uniform
on `{1..W}`", then exhausting weight ≤ 9 has already eliminated `9/W` of that hypothesis** — 56 % at
`W = 16`, 28 % at `W = 32`, 14 % at `W = 64`. That, and only that, is what the fleet's 2^33-candidate
sweeps have been purchasing. Stated as a decision rule: **the sweeps are worth continuing exactly to
the extent that you hold a non-uniform-designer prior, and worth nothing otherwise** — which is the
same conclusion AB and Y reached about the complement, now with the exchange rate attached.

**Struck, mine, shown struck:**

1. ~~"`[0,N)` versus `[0,2^256)` is a `2^-128`-scale perturbation of every popcount class."~~
   **STRUCK.** It is a `2^-128`-scale perturbation *in total*, but its distribution over classes is
   nothing like uniform: **for every `b ≤ 127` the effect is EXACTLY ZERO** — `#{k < N : popcount = b}
   = C(256,b)` on the nose — and the first class where they differ is `b = 128`, by exactly 1. The
   entire deficit sits at `b ≥ 128`. Proof and consequence in §1.3. This matters: every headline
   number in §4 lives at `b ≤ 56`, so **`N` versus `2^256` is not an approximation in any of them,
   it is an identity.**
2. ~~"the fleet's per-bit figure is `< 2^-200.8` (X)"~~ — **corrected upward, not struck.** X's
   arithmetic is right on X's own inputs; X's table predates Y's complement sweep. Adding Y's family
   gives **`2^-200.12`** of a single-bit fibre. Still `< 2^-200`; the qualitative claim survives
   intact. Details and the two input-size discrepancies in §2.4.

---

# 1. THE EXACT PRIOR

## 1.1 The model, stated inside the claim

Everything in §§1–4 is conditional on: **`k₀ := log_G T` is uniformly distributed on `[0,N)`, with
`N` the secp256k1 group order**, and **`w := popcount(k₀)`**. Nothing in `EQUATIONS.txt` forces this
prior; it is the maximum-entropy stand-in for "I know nothing about how the scalar was chosen".
§5 is the attempt to find evidence against it. Nothing else in this document depends on the
uniformity assumption being *true*, only on it being the stated model.

## 1.2 The computation

`ac_prior.py` / `ac_post.py`: a digit-DP over the binary expansion of `N`. For each bit position `i`
where `bit_i(N) = 1`, every `k` agreeing with `N` above `i`, carrying 0 at `i`, and free below,
is `< N`; that contributes `C(i, j)` to class `ones_above(i) + j`. A few hundred `comb` calls.

**Verification, as required:**

* `Σ_b cnt[b] == N` — **exact integer equality, True.**
* A **second, structurally different DP** (MSB→LSB, carrying a `tight`/`free` pair) reproduces the
  vector **elementwise**. `ac_post.py` asserts `cnt == cnt2`.
* `Σ_b C(256,b) == 2^256` — True.
* `Σ_b (C(256,b) − cnt[b]) == 2^256 − N` — True, and every term is `≥ 0`.

## 1.3 `[0,N)` versus `[0,2^256)` — quantified, and it is not what I expected

`2^256 − N = 432420386565659656852420866394968145599 = 2^128.3457`, so the total deficit is
`2^-127.6543` of the space. But:

> **`gap[b] := C(256,b) − #{k < N : popcount(k) = b}` is exactly 0 for every `b ≤ 127`, and
> `gap[128] = 1`.**

Proof (and the DP agrees): `2^256 − N ∈ [2^128, 2^129)`, so `N ∈ (2^256 − 2^129, 2^256 − 2^128]`,
so bits 255…129 of `N` are all 1 and bit 128 is 0. Any `k ∈ [N, 2^256)` therefore has bits 255…129
all 1, giving `popcount(k) ≥ 127`; equality would force `k = 2^256 − 2^129 < N`, contradiction.
Hence **every `k ∈ [N, 2^256)` has `popcount ≥ 128`.**

Two consequences used later:

* **Every tail probability at `B ≤ 127` computed under `[0,2^256)` is exact under `[0,N)`.** All of
  §4 is in that range. The worst relative deficit anywhere in the live range `[10,246]` is at
  `b = 246` and is `2^-10.16`; at `b = 128` it is `2^-252`.
* **The "essentially unique solution" caveat is harmless in the direction it could have hurt.**
  A second legal ON-set exists iff `k₀ < 2^256 − N` (probability `2^-127.65`, coordinator's row 9,
  re-derived here). When it exists, the second representative is `k₀ + N ∈ [N, 2^256)` and therefore
  has **popcount ≥ 128 always**. So the alternative ON-set can never be the sparser one unless
  `w ≥ 128` already, and *defining* `w` as `popcount(k₀)` rather than as the minimum over legal
  ON-sets changes nothing at any `B ≤ 127`.

## 1.4 The prior itself

Mode `b = 128`, `P(w = 128) = 2^-4.3272`; mean 128, sd 8 exactly under `[0,2^256)` and to within
`2^-127` under `[0,N)`. Full vector in `ac_prior.json` as exact integers.

---

# 2. THE POSTERIOR

## 2.1 Which facts are genuine conditioning events on `w`

I went through `K_CONSTRAINTS.md §1` row by row. **Exactly two of the ten rows are events about `w`:**

| | event | source | effect on the support |
|---|---|---|---|
| **A** | `w ≥ 10` | X, unsigned weight ≤ 9 exhausted from cold, 0 hits, range totals exact | kills `b ∈ [0,9]` |
| **B** | `w ≤ 246` | Y, complement weight ≤ 9 exhausted, `ALLDONE`, per-`i0` sum `= C(256,5)` exactly | kills `b ∈ [247,256]` |

**A and B are neither nested nor overlapping** — they restrict disjoint ends of the support. The
posterior is the prior restricted to `[10,246]` and renormalised. **Nothing is multiplied**; treating
A and B as independent events to be multiplied would be a category error, since they are not events
that can co-occur in a way that compounds.

**Row 1b (signed weight ≥ 8) is NESTED INSIDE A and contributes no new bound on `w`**: a `k` of
unsigned weight `v` has signed weight `≤ v`, so signed ≥ 8 ⟹ unsigned ≥ 8, which A already
subsumes. It is *not* redundant as a set-exclusion — it removes low-run-length `k` of high unsigned
weight — and that residual effect is priced in §2.2. The fleet's own note that "signed ≤ m and
unsigned ≤ 9 are separate citable statements" is right about the *classes* and irrelevant to the
*bound on `w`*.

Rows 2–8 (BSGS interval, 34-bit window, small multiples, `a+bλ` box, endomorphism orbit) are
exclusions of explicit sets of `k` that cut across every weight class. Row 9 is arithmetic, not an
observation. Row 10 is a property of `N`. **AB's §8 instance-side constraint is UNSETTLED and is
therefore not a conditioning event** — I have not used it.

## 2.2 What rows 1b–8 do to the class weights, computed rather than waved away

`ac_cond.py` computes, per popcount class, an **upper bound on the number of `k` removed** by the
non-`w` families (a union bound, which is the right direction for bounding the movement):

* **signed/NAF weight ≤ 7:** exact joint census by DP, using the identity
  `NAFweight(v) = popcount(3v ⊕ v)` (verified against Reitwiesner's algorithm on 5,000 values and 5
  256-bit values; the DP verified against brute force over all `2^14`). **The set has 733,018,571,531,264
  = 2^49.38 elements**, against X's enumeration size `Σ_{m≤7} C(256,m)2^m = 2^50.60` — the enumeration
  over-counts because signed representations are not unique, and it is the *set* that conditions.
* **the wrapped half of that family** (`k = N − u` with `u` representable) **cannot touch any class
  below `b = 36`**: if `popcount(k) = b` and `NAFwt(u) ≤ 7` then `N = k + u` is a signed
  representation of `N` with `≤ b + 7` nonzero digits, and the minimum such is `NAFweight(N) = 43`
  (**recomputed here: `popcount(N) = 192`, `NAFweight(N) = 43`**, confirming X's row 10), so `b ≥ 36`.
* **`k < 2^52`:** `C(52,b)` exactly. **`N − k < 2^52`:** exact, by differencing the digit-DP at the
  two endpoints (sums to `2^52 − 1`, checked).
* **34-bit window:** `223·C(33,b−1)` exactly (each such `k` is `a·2^s` with `a` odd, uniquely).
* **`a+bλ` box, small multiples, orbit:** no useful per-class structure; bounded by total size.

**Result — the worst per-class relative distortion anywhere in `[10,246]` is `7.12 × 10⁻⁵` (at
`b = 10`), and it falls below `2^-100` for every `b ≥ 36`.** That bound is itself loose by orders of
magnitude, because it charges the whole `a+bλ`/small-multiple/orbit mass to every class. §4 carries
the resulting rigorous bracket on each headline tail; the bracket is **3.4 × 10⁻⁹ bits wide at
`B = 14`** and narrower above. **Rows 1b–8 do not move any number in this document.**

## 2.3 How much the entire campaign moved the distribution of `w`

Exact, from `ac_post.py`:

```
P_prior(w ≤ 9)     = 11711713815280289 / N            = 2^-202.6212
P_prior(w ≥ 247)   = 11689700375120489 / N            = 2^-202.6239
P_prior(excluded)  = 23401414190400778 / N            = 2^-201.6226
```

(The two ends are **not** equal — under `[0,N)` the high end is very slightly lighter. Under
`[0,2^256)` they are exactly equal, `2^-202.6212` each, as AB's symmetry check reports.)

> **TV(prior, posterior) = `P_prior(excluded)` = `2^-201.6226`.**
> **Information gained about `w` by the entire campaign = `−log₂(1 − P_excl)` = `2^-201.0938` bits
> = `2.916 × 10⁻⁶¹` bits.** `KL(posterior ‖ prior)` agrees to leading order.

For scale: this is the campaign's *whole* yield on the question it has been organised around. One
flip of a coin biased `0.5 + 10^-20` carries `2.9 × 10^-40` bits — **`10^21` times more**
(`ac_xcheck.out`).

## 2.4 Cross-check against X's `< 2^-200.8` per-bit figure — consistent, with a small correction

Recomputed here as a union over every excluded family (`ac_post.py`):

| family | size | log₂ |
|---|---|---|
| unsigned weight ≤ 9 (X) | 11,711,713,815,280,289 | 53.38 |
| complement weight ≤ 9 (Y) | 11,711,713,815,280,289 | 53.38 |
| BSGS both ends, `2·2^52` (X) | 9,007,199,254,740,992 | 53.00 |
| NAF weight ≤ 7 **as a set** (X row 1b) | 733,018,571,531,264 | 49.38 |
| `a+bλ` box (Q) | 17,592,186,044,416 | 44.00 |
| 34-bit window (Q) | 1,915,555,414,016 | 40.80 |
| small multiple `m ≤ 10^7` (Q) | 2,560,000,000 | 31.25 |
| endomorphism orbit (Q/X) | 1,536 | 10.58 |
| **union (upper bound)** | **33,183,155,758,292,802** | **54.88** |

* **X's list, with X's own family sizes: I get `2^54.2037 → 2^-200.7963` (`ac_xcheck.out`), against
  X's stated `2^54.20 → 2^-200.80`. X's arithmetic reproduces exactly.** My recomputation of the
  same list with my own family sizes gives `2^54.32 → 2^-200.68`; the 0.12-bit
  difference is two input sizes, not an error in X's sum: I use `(2·2^21)² = 2^44.00` for the `a+bλ`
  box where X has `2^43.58`, and `223·2^33 = 2^40.80` for the deduplicated window family where X has
  `2^41.80`. Neither matters.
* **Adding Y's complement family — which postdates X's table — gives `2^54.88 → `2^-200.12`.**

**Are the two figures consistent?** Yes, and they measure different things, so the comparison has to
be made carefully:

* X's is *(union of excluded sets) / 2^255*, a statement about one **bit** of `k`.
* Mine is *(excluded weight-classes) / N*, a statement about the **distribution of `w`**.
* The consistency requirement is `TV(w) ≤ union/2^256`, since the classes A and B remove are a
  subset of the union. **`2^-201.62 ≤ 2^-201.12` — holds.** The slack is exactly the families that
  remove `k` without removing whole weight classes.

**Neither figure is wrong. X's headline should now read `< 2^-200.1` rather than `< 2^-200.8`, and
the claim it supports — that no search this campaign has run moves any single bit of `k` — is
untouched.**

---

# 3. THE TABLE

Counts and CDF values are **exact rationals** (integer numerator over integer denominator; the three
headline ones are printed in full in `ac_table.out`). The `log₂` renderings are floating-point,
accurate to `~10^-13`, and were **cross-checked against a 60-digit `Decimal` recomputation
(agreeing to `2.8 × 10^-14`) and against `mpmath`'s regularised incomplete beta
`P(Bin(256,½) ≤ B) = I_{1/2}(256−B, B+1)` at 60 dps (agreeing to 6 decimal places in the exponent —
`ac_thirdcheck.out`).** Three independent routes; no exponent in this document rests on one sum.

**Prior and posterior give the same table to every digit shown** — conditioning on `10 ≤ w ≤ 246`
changes the quantiles by `2^-201`, which is the whole point of §2.3.

| `ε` | one-sided: smallest `B` with `P(w ≤ B) ≥ 1−ε` | achieved `1 − P(w≤B)` | two-sided equal-tailed `[L,U]`, each tail `≤ ε/2` | achieved coverage |
|---|---|---|---|---|
| `10^-1` | **138** | `2^-3.40` | **[115, 141]** | `1 − 2^-3.45` |
| `10^-2` | **147** | `2^-7.10` | **[107, 149]** | `1 − 2^-7.14` |
| `10^-3` | **153** | `2^-10.49` | **[102, 154]** | `1 − 2^-10.14` |
| `10^-6` | **166** | `2^-20.62` | **[89, 167]** | `1 − 2^-20.56` |
| `2^-20` | **166** | `2^-20.62` | **[89, 167]** | `1 − 2^-20.56` |
| `2^-40` | **183** | `2^-40.04` | **[72, 184]** | `1 − 2^-40.42` |
| `2^-80` | **207** | `2^-81.37` | **[49, 207]** | `1 − 2^-80.37` |

(`10^-6` and `2^-20` coincide because `2^-20 = 9.54 × 10^-7`; that is not a copy-paste error.)

The CDF, for reference (posterior; `ac_table.out`):

| `B` | `P(w ≤ B)` | `B` | `P(w ≤ B)` | `B` | `P(w ≤ B)` |
|---|---|---|---|---|---|
| 10 | `2^-198.05` | 45 | `2^-87.88` | 90 | `2^-19.70` |
| 15 | `2^-176.76` | 50 | `2^-77.21` | 100 | `2^-11.80` |
| 20 | `2^-158.05` | 55 | `2^-67.43` | 104 | `2^-9.27` |
| 25 | `2^-141.27` | 60 | `2^-58.47` | 110 | `2^-6.13` |
| 30 | `2^-126.06` | 65 | `2^-50.28` | 120 | `2^-2.52` |
| 35 | `2^-112.19` | 70 | `2^-42.83` | 128 | `2^-0.93` |
| 40 | `2^-99.51` | 80 | `2^-29.98` | 148 | `2^-0.0074` |

`P(w ≥ 256−B)` matches `P(w ≤ B)` to three decimal places in the exponent at every `B` in the table
— the residual difference is the `[0,N)` truncation, visible only in the third decimal.

---

# 4. THE COLLISION WITH AB'S COST BANDS

AB's bands, as stated in `UPPER_BOUND_MAP.md` §S3/§S5 and used here **as AB's numbers, not
re-derived**: `w ≲ 56` to beat rho at 2^40 memory (`w ≲ 52` at this box's 2^30); `w ≲ 24` to be
actionable; `w ≤ 14` reachable on this box at 2^47 time / 2^30 memory.

| band | `P` under the posterior | rigorous bracket incl. families 1b–8 | odds |
|---|---|---|---|
| `w ≤ 14` (this box) | **`2^-180.7801`** = 3.7998 × 10⁻⁵⁵ | width `3.4 × 10^-9` bits | 1 in `2^180.78` |
| `w ≤ 24` (actionable) | **`2^-144.4870`** = 3.1995 × 10⁻⁴⁴ | width `< 10^-13` bits | 1 in `2^144.49` |
| `w ≤ 52` (rho, 2^30 mem) | `2^-73.1941` | — | 1 in `2^73.19` |
| `w ≤ 56` (rho, 2^40 mem) | **`2^-65.5704`** = 1.8254 × 10⁻²⁰ | width `< 10^-13` bits | 1 in `2^65.57` |
| `w ≤ 64` (rho, 2^60 mem) | `2^-51.8589` | — | 1 in `2^51.86` |
| `w ≤ 104` (AB's struck unbounded-memory crossover) | `2^-9.2739` | — | 1 in `2^9.27` |

**Each exponent verified a second and third way** (`ac_post.py`, `ac_thirdcheck.out`):

| `B` | exact sum | single-term lower `C(256,B)/2^256` | geometric-ratio upper | entropy/Chernoff upper `2^{256(H(B/256)−1)}` | mpmath `I_{1/2}` | brackets hold |
|---|---|---|---|---|---|---|
| 14 | `2^-180.7801` | `2^-180.8653` | `2^-180.7797` | `2^-177.6679` | `2^-180.780135` | ✔ |
| 24 | `2^-144.4870` | `2^-144.6430` | `2^-144.4861` | `2^-141.0907` | `2^-144.487001` | ✔ |
| 56 | `2^-65.5704` | `2^-66.0367` | `2^-65.5656` | `2^-61.9831` | `2^-65.570354` | ✔ |

The Chernoff bound brackets the exact sum from above in every row, as it must; the single-term
bound brackets from below; the geometric-ratio bound is tight to 0.005 bits. A factor of 2 lost in a
tail would show up as a 1.0 shift in one of these columns and does not.

> ### The conclusion, in one sentence, plain language
>
> **Under the model `k₀ ~ uniform[0,N)`, the probability that `w` is small enough for any known
> mechanism to help is `1.8 × 10^-20` for the most generous band (`w ≤ 56`, which merely ties rho at
> a memory budget nobody here has) and `3.2 × 10^-44` for the band that would actually be actionable
> (`w ≤ 24`) — the probability that `w` is small enough to help is astronomically small, and the
> single fact that carries the whole result is that a uniform 256-bit scalar has weight `128 ± 8`
> and every useful band is more than nine standard deviations below that.**

Restated for planning: `w ≤ 24` is a 13σ event; `w ≤ 56` is a 9σ event; `w ≤ 14` is a 14.25σ event.
**With probability `1 − 2^-65.57`, `w` sits above every band in AB's table.** And by §2.3, running
more sweeps of the kind the fleet has been running changes that by `2^-201`.

---

# 5. THE ONE PLACE THE UNIFORMITY MODEL COULD BE WRONG

The model is the load-bearing assumption. Rule 2 forbids investigating how the instance was
generated, and nothing below does: **every test is a question about a given point on a given curve**,
answerable with no reference to a generator, a seed, an emission order, or a coefficient template.
`ac_point.py`, `ac_point2.py`; raw output in `ac_point.out`, `ac_point2.out`, `ac_point.json`.

## 5.1 Structure, re-verified here from scratch

| test | result |
|---|---|
| `p` prime; `p = 2^256 − 2^32 − 977` | **True** |
| `N` prime | **True** |
| `a = 0` (so `j = 0`, `Aut(E) = μ₆`) | **True** |
| `T` on curve; `N·T = O`; `T ≠ O`; `T ≠ G` | **True** |
| `b/7` is a 6th power mod `p` (⟹ `F_p`-isomorphic to secp256k1) | **True** |

## 5.2 Every test, with its computed answer

| # | question | answer | value tested |
|---|---|---|---|
| 1 | does `T` lie in a small-index subgroup? | **NO — impossible** | `N` prime ⟹ the only subgroups are `{O}` and `E(F_p)`; `T ≠ O` ⟹ `ord(T) = N`, index 1. Nothing to check beyond `N` prime, which is **True**. |
| 2 | are there order-3 points at all (`x = 0`)? | **NO** | `b` is a non-residue mod `p`; consistent with `3 ∤ N` |
| 3 | is `T` fixed by any non-trivial automorphism of `E`? | **NO** | all 5 non-trivial elements of `μ₆` checked as `(x,y) ↦ (ζ₃^i x, ±y)`; fixing `T` would need `x(T) = 0` (**False**) or `y(T) = 0` (**False**) |
| 4 | is `T = ±ζ₃^j·(c·G)` for a tiny scalar `c`? | **NO** | swept **all `c ≤ 2^20`** against all 3 values of `ζ₃^j·x(T)`; `x`-equality covers both signs since `x(−P) = x(P)`. This subsumes `T = ±G, ±λG, ±λ²G` and the "few hundred points" the brief names. Engine cross-checked against a naive affine engine for `c = 1..299` |
| 5 | is `T` in the 1,536-point endomorphism orbit `{±λ^j 2^i G}`? | **NO** | all 256 doublings × 3 `ζ₃` images. Independently confirms Q/X row 6 |
| 6 | is `G = ±ζ₃^j·(m·T)`, i.e. `T = G/m` for tiny `m`? | **NO** | swept **all `m ≤ 2^20`** |
| 7 | is `x(T)` small, or near `0`, `p`, `2^255`, `2^256`? | **NO** | `x(T)` is 255 bits; distances `2^254.06`, `2^255.57`, `2^253.94`, `2^255.57` — all `≫ 2^40` |
| 8 | is `y(T)` small, or near `0`, `p`, `2^255`, `2^256`? | **NO** | `y(T)` is 245 bits; distances `2^244.67`, `2^256.00`, `2^255.00`, `2^256.00` — all `≫ 2^40` |
| 9 | is `x(T)` or `y(T)` of low Hamming weight? | **NO** | `popcount(x) = 130`, `popcount(y) = 123` (uniform: 128 ± 8, i.e. +0.25σ and −0.63σ); `NAFweight(x) = 87`, `NAFweight(y) = 92` (uniform ≈ 85.3) |
| 10 | is `x(T)` or `y(T)` a perfect square over ℤ? a perfect cube? a power of two? | **NO** to all six | integer `nthroot` |
| 11 | `T = ±G`? `x(T) = y(T)`? `gcd(x(T), p) ≠ 1`? | **NO** to all three | — |
| 12 | is `x(T)` a QR mod `p`? is `y(T)`? is `y` the even root? | `x`: **no**; `y`: **no**; `y < p/2`: **yes** | recorded for completeness; each is a 50/50 under **any** model and is not evidence of anything |

**Clean sweep of "no".** That is the expected outcome and it is worth recording: **it removes the
last cheap thing anyone can point at.**

## 5.3 The one observation that is not perfectly bland, reported rather than buried

**`y(T)` is 245 bits, not 256 — its top 11 bits are zero.** For a uniform point,
`P(min(y, p−y) < 2^245) = 2·2^245/p ≈ 2^-10.0`, so this is roughly a 1-in-1000 coincidence. I ran
about 30 tests; the chance that *at least one* of them landed in a 1-in-1000 cell is a few percent,
so **as evidence this is worth essentially nothing**, and I would not have mentioned it if the brief
had not asked for every value.

It is also worth being precise about why it could not help even if it were real:

> **A property of `T` bounds `w` only if the hit exhibits `k` explicitly.** Tests 4, 5 and 6 have
> that character — a hit there hands you `k = ±λ^j c` or `k = ±λ^j/m` and hence `w` outright. Tests
> 7–12 do not: a property of the *point* with prior probability `q` cuts the scalar space down to a
> set of density `q` that you cannot enumerate, which is exactly the mechanism AB marked **DEAD**
> as "character sums / analytic" (#6). So the small `y` is not a lead, it is a coincidence with no
> handle on it.

## 5.4 What actually remains, stated plainly

The uniformity model survives every cheap test that could have refuted it from the instance file
alone. What it cannot survive, and what nobody here is permitted to investigate, is the possibility
that the scalar was **chosen** rather than drawn. That possibility is not testable from `T`; it is a
prior, and it belongs to the user, not to the fleet. The honest decomposition, with the model inside
the sentence:

> Writing `π` for your prior that the instance's scalar was deliberately made sparse, and `1−π` for
> uniform,
> `P(w ≤ 24) = π · P(w ≤ 24 | sparse-by-design) + (1−π) · 2^-144.49`.
> **The second term is nothing. Your answer to "is `w` small?" is `π`, and `π` is not something any
> computation in this repository can move.**

And the campaign's searches touch exactly that term and nothing else. **Under the model "a designer
picks `w` uniform on `{1..W}`", exhausting weight ≤ 9 has removed `9/W` of `π`** — 56 % at `W = 16`,
28 % at `W = 32`, 14 % at `W = 64`, 7 % at `W = 128 `— and pushing this box to its `w ≤ 14` limit
(AB §S3) would take that to `14/W`. **That exchange rate is the entire case for continuing to
search, and it should be stated in those terms rather than as "progress on bounding `w`."**

---

# 6. REPRODUCTION

```
cd solve_lab/agentAC_work
PYTHONDONTWRITEBYTECODE=1 python3 ac_prior.py     # digit DP, sum == N check
PYTHONDONTWRITEBYTECODE=1 python3 ac_cond.py      # NAF DP self-test + per-class removals
PYTHONDONTWRITEBYTECODE=1 python3 ac_post.py      # prior, posterior, quantiles, tails, brackets
PYTHONDONTWRITEBYTECODE=1 python3 ac_extra.py     # operational CDF + Decimal cross-check
PYTHONDONTWRITEBYTECODE=1 python3 ac_point.py     # point tests, sections 0-2,4,5
PYTHONDONTWRITEBYTECODE=1 python3 ac_point2.py    # the two 2^20 reachability sweeps (~3 min)
```

Every script asserts its own invariants and exits non-zero on failure; none writes outside this
directory; none writes a `.pyc` anywhere.
