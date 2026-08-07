# THEOREM_B_AUDIT — red-team report on agent AB's Theorem B

Agent AG. Adversarial audit, not a confirmation exercise. Everything below was recomputed in this
directory from scratch (`ag_recompute.py`, `ag_attack.py`, `ag_verify.py`, `ag_disk.py`); AB's
scripts were read only *after* my own model was written, and where I agree I say whether the
derivation was independent or a re-check of AB's formula.

Reference point recomputed here, not quoted: **rho = √(πN/4)/√6 = 2^126.5333**.

---

## S0. VERDICT IN ONE PARAGRAPH

**Theorem B is false in the exact form it is stated, and its two headline numbers are both wrong —
but its qualitative conclusion survives, in a corrected form that is stronger than AB's.** The
covering-code attack the coordinator asked for **fails**: AB's per-ball pricing sits at most **2^2.66**
above an amortisation-proof information-theoretic floor, so no code structure, shared half-list or
overlap can move the break-even by more than **12 in B**. The quantifier attack **succeeds**: not
every search-based upper bound is a covering of `{wt > B}`, and the counterexample is the low-weight
sweep the campaign is already running — AB's model overprices it by up to **2^98**. The memory
attack **succeeds**: `B = 148` and `w = 106` are unbounded-memory numbers that AB corrected
elsewhere in round 3 and never propagated; **memory-aware at 2^30 the dead band is `[54, 200]`**, not
`[107,147]`. Three further arithmetic defects are recorded in §4, all of them in the conservative
direction.

**Round 2 (§6a): four of my own claims are struck or corrected, two of them errors of exactly the
class I convicted AB of.** The one that did not survive is §4.4's *"`cover(B) = 2^128.000` exactly
for every `B ≤ 148`"* — true only down to `B = 142`, and contradicted by my own §1.3 table. Attacks
1, 2 and 3 stand as reported; AB conceded 2 and 3 in full and independently reproduced both the
`rep(W)` fix and the `√Z` floor.

---

## S1. WHAT I ATTACKED, AND WHAT HAPPENED

| # | attack | outcome |
|---|---|---|
| 1 | covering **code** instead of covering; amortisation across overlapping half-lists | **FAILS** — headroom ≤ 2^2.66 (§1) |
| 2 | is the quantifier "every search-based upper bound is a ball covering" honest? | **SUCCEEDS** — theorem false as stated (§2) |
| 3 | memory-aware break-even and crossover | **SUCCEEDS** — 148 → 201 at 2^30 (§3) |
| 4 | arithmetic audit of `rep`, volumes, self-certificate, monotonicity, tables, disk | **3 defects found**, all conservative (§4) |
| 5 | algebraic / non-generic certificates against Theorems B ∧ D jointly | **gap confirmed and named** (§2.3) |

---

## 1. ATTACK 1 — covering codes. FAILS, with a bound on how much it could ever have bought.

### 1.1 A different model from AB's

AB priced *balls*. I priced *rectangles*, and the rectangle model contains every amortisation the
coordinator asked about.

Because the exponents are **distinct powers of two**, splitting the 256 positions in half gives
`k = k_H·2^128 + k_L` **over ℤ with no carries**, hence `wt(k) = wt(k_H) + wt(k_L)` exactly. Every
meet-in-the-middle — ball, code, design or ad hoc — stores a set `A ⊆ F_2^128` of high halves and
streams a set `B ⊆ F_2^128` of low halves, and certifies precisely the **combinatorial rectangle
`A × B`**. Two balls "sharing most of their half-lists" are two rectangles with overlapping sides;
if they share the same `A` they **merge into a single rectangle `A × (B₁ ∪ B₂)` at no extra cost**.
The merge is free, so it is already in the model, and after merging WLOG every rectangle has a
distinct `A`.

### 1.2 The floor, twice, from two unrelated directions

**(a) Combinatorial.** A rectangle of area `z` costs `≥ √z` to search by *any* algorithm (classical
MITM is `|A|+|B| ≥ 2√z`; a memory-free vOW/rho pass over the same rectangle is `≥ √z`). To prove
`w ≤ B` the rectangles must cover `{(a,b) : wt(a)+wt(b) > B}`, of size `Z = |{wt > B}|`. Then

> `total ≥ Σᵢ √zᵢ ≥ √(Σᵢ zᵢ) ≥ √Z`   (subadditivity of `√`)

and this is **indifferent to how the covering is structured** — code, design, random pile, arbitrary
overlap. That is the whole answer to "does a code structure amortise work across balls".

**(b) Generic-group counting, with no rectangles in it at all.** After `m` generic queries the
algorithm holds `m` elements `σ(αᵢ+βᵢk)`; it can rule out a scalar `x` only if `x` is the (unique)
root of one of the `≤ C(m,2)` affine collision equations or of one of the `m` direct tests. Ruling
out `Z` scalars needs `m ≥ √(2Z)`. **Same floor.**

### 1.3 Measured headroom

| B | AB's covering cost | `√Z` floor | headroom |
|---|---|---|---|
| 245 | 2^30.03 | 2^29.01 | 2^1.03 |
| 220 | 2^74.75 | 2^71.90 | 2^2.84 |
| 200 | 2^96.97 | 2^94.29 | 2^2.68 |
| 180 | 2^112.56 | 2^109.96 | 2^2.60 |
| 152 | 2^125.70 | 2^123.06 | 2^2.63 |
| 148 | 2^126.85 | 2^124.20 | 2^2.66 |

**Max headroom over `B ∈ [120,251]`: 2^2.66, a factor of 6.3.** Translating that into the numbers the
campaign plans against: even a perfect, information-theoretically optimal certifier moves

* **break-even `B`: 148 → 136** (best possible, 12 lower)
* **crossover `w`: 106 → 118** (best possible, 12 higher)
* **the dead band where nothing beats solving: `[107,147]` → `[119,135]` — it narrows but does not close.**

### 1.4 Two structural reasons the amortisation was never there

* **The Hamming symmetry is not an algorithmic symmetry.** `{wt > B}` is invariant under all 256!
  coordinate permutations, which is exactly what makes good covering codes cheap in coding theory.
  But the `2^i·G` are unrelated group elements, so a permutation of bit positions is not a map the
  algorithm can apply to its lists. **Orbit-based amortisation, the main engine of covering-code
  constructions, is unavailable here.** AB never says this; it is the real reason the coding-theory
  intuition does not transfer.
* **HGJ/BCJ representations die for a stronger reason than AB gives.** AB says the technique needs
  `s mod M` and `s = k₀` is unknown. Sharper: **even if `k₀ mod M` were given free, the technique
  would still fail**, because the representation filter must be *enumerable* — you must be able to
  list the filtered sublist directly — and a constraint on a **group element** is only *testable*,
  never enumerable. You would pay the full unfiltered enumeration to build the filtered list.

**VERDICT: Attack 1 fails. AB's per-ball, no-amortisation pricing is conservative, and the
conservatism is bounded by 2^2.66.** This is the strongest part of Theorem B and it should be
restated in the rectangle form, which proves more with less machinery.

---

## 2. ATTACK 2 — the quantifier. SUCCEEDS: the theorem is false as stated.

> Theorem B: *"Every search-based upper bound on `w` is a Hamming-ball covering of `{wt > B}`."*

### 2.1 The counterexample is the algorithm the campaign is already running

To certify `w ≤ B` you may **search `{wt ≤ B}`** instead of covering `{wt > B}`. A **hit** hands you
`k`, hence a certified (indeed exact) `w ≤ B`. It is search-based; it is unconditional on success;
and it is a covering of the **complement** of the set Theorem B says must be covered.

| B | Theorem-B covering cost | search-`{wt ≤ B}` cost | AB overprices by | `P(w ≤ B \| null)` |
|---|---|---|---|---|
| 10 | 2^128.00 | 2^30.03 | **2^97.97** | 2.5e−60 |
| 20 | 2^128.00 | 2^50.26 | **2^77.74** | 2.6e−48 |
| 40 | 2^128.00 | 2^79.83 | 2^48.18 | 1.1e−30 |
| 56 | 2^128.00 | 2^96.97 | 2^31.03 | 1.8e−20 |
| 106 | 2^128.00 | 2^126.32 | 2^1.68 | 3.5e−03 |

**Consequence for a published claim.** AB's round-2 line *"the gap to the generic lower bound is now
≤ 2^3 everywhere, nothing large is left on the table"* is **false outside the range AB tabulated
(`B ≥ 128`)**. At `B = 20` AB's model says 2^128.0 while Theorem D's generic floor is 2^49.0 — and
**the floor is the truth**, because the low-weight sweep achieves 2^50.3.

**What actually rescues the qualitative conclusion is not cost — it is success probability.** This
certifier returns an upper bound only in the branch where the bound is true, and in that branch it
has produced `k₀` itself. **Theorem B needs the hypothesis "zero-error / correct for every `w`",
which appears nowhere in its statement.**

### 2.2 A second non-covering certifier: intervals

Kangaroo decides `k₀ < 2^t` in `≈2^{t/2+1}`; a **hit** proves `w ≤ t`. It is search-based and is not
a Hamming-ball covering — an interval is not a ball.

| t | cost | proves on a hit | `P(hit \| null)` |
|---|---|---|---|
| 128 | 2^65 | `w ≤ 128` | 2^−128 |
| 200 | 2^101 | `w ≤ 200` | 2^−56 |

At `t = 128` that is **2^63 below AB's curve** for the same conclusion. AB kills it in §9.13 by the
right argument (a miss says nothing) — but that argument is *success probability*, not Theorem B.

### 2.3 What Theorems B and D together do and do not exclude

**Excluded:** (i) any covering/rectangle MITM, by §1's floor; (ii) any **generic-group** algorithm
deciding `w ≤ B` on average, by Theorem D.

**NOT excluded — and this is the honest gap:**

* **Non-generic algebraic certificates.** A Nullstellensatz/Positivstellensatz refutation of
  `{ladder equations} ∧ {Σ sᵢ > B}` is a proof of `w ≤ B` that is not a covering (so Theorem B
  misses it) and works in the coordinate ring of `E/F_p` (so Theorem D, which explicitly excludes
  the encoding, misses it too). **Note the certificate *size* is never the barrier: `k₀` is itself a
  256-bit certificate verifiable by one scalar multiplication, so a short proof of `w ≤ B` always
  exists.** The only possible barrier is the cost of *finding* one — which is exactly AB's `d_reg`
  question, and `d_reg` has been measured at `n = 2, 3` only.
* **Any non-generic decision procedure returning `w` without locating `S`.** Generic case covered by
  D; non-generic case covered by nothing on the page.
* **Non-search arithmetic bounds** (`k₀ < N` + digit-DP ⇒ `w ≤ 255`). AB's own §7 exhibits one, which
  is direct evidence that the word "search-based" in the quantifier is load-bearing.

### 2.4 The repair: a trichotomy and a unified law

> **Restated (AG).** For every `B`, the cost of *deciding* `[w ≤ B]` is
> **`Θ( √( min(|{w ≤ B}|, |{w > B}|) ) )`**, achieved to within 2^2.7 by searching whichever side is
> smaller, and matched from below by Theorem D's generic bound **at every `B`, not only `B ≥ 128`**.

| B | `√min(·,·)` | best known | which side |
|---|---|---|---|
| 20 | 2^48.98 | 2^50.26 | search `{w ≤ B}` |
| 60 | 2^98.77 | 2^100.56 | search `{w ≤ B}` |
| 106 | 2^123.93 | 2^126.32 | search `{w ≤ B}` |
| 120–148 | ≈2^126–127 | **2^126.53** | **rho (just solve)** |
| 180 | 2^109.96 | 2^112.56 | cover `{w > B}` |
| 245 | 2^29.01 | 2^30.03 | cover `{w > B}` |

and the dichotomy becomes a **trichotomy**:

> Every search-based upper bound on `w` is **(i)** vacuous (`B ≳ 148`, ≥ 2.5σ above the null mean),
> **(ii)** priced at what solving costs (`B` in the dead band), or **(iii)** cheap but one-sided —
> it produces an upper bound only with the null probability that the bound is true, and in that
> branch it has produced `k₀`, so it *is* solving restricted to a subclass.

Branch (iii) is missing from AB's statement, and it is the branch the campaign's whole low-weight
programme lives in. **Nothing strategic changes; the classification does.**

---

## 3. ATTACK 3 — memory. SUCCEEDS: the two headline numbers are unbounded-memory artefacts.

AB's round 3 corrected the `(time, memory) → w` reach table with vOW golden-collision search
(`T = rep·c·L^1.5/√M`, `c = 2.5`) and correctly withdrew "crossover moves to 104". **It did not
propagate the same correction to the break-even or to the AUTHORITATIVE SUMMARY's `w = 106 / B = 148`.**
Those are the unbounded-memory column of the very table AB struck through.

Re-derived in **one model** (AB's own vOW model, memory `M`, rho at 2^126.533, per-ball cost capped
at rho since any ball can be searched by solving outright):

| memory | crossover `w` | **break-even `B`** | dead band |
|---|---|---|---|
| **2^30 (this box)** | ~~52~~ → **53** | **201** | **`[54, 200]`** |
| 2^35 | 54 | 200 | `[55,199]` |
| 2^40 | 56 | 198 | `[57,197]` |
| 2^50 | 60 | 194 | `[61,193]` |
| 2^60 | 64 | 190 | `[65,189]` |
| 2^80 | 72 | 181 | `[73,180]` |
| unbounded | 106 | 148 | `[107,147]` |

> **CORRECTED (round 2): the crossover at 2^30 is 53, not 52. AB is right, and the reason is mine.**
> I derived the odd-`W` `rep` fix in §4.1 and then **ran my own memory-aware table with the
> un-fixed `rep`**. Recomputed: `rep_ab` → crossover 52, `rep_tight` → crossover **53**; break-even
> is 201 under both. At `M = 2^30` the ball time is 2^124.486 at `w = 52` and **2^126.424 at
> `w = 53`**, still under rho = 2^126.533; `w = 54` costs 2^127.417.
> **This is a propagation failure of exactly the kind I opened §3 by accusing AB of** — find a fix
> in one section, fail to carry it into the table in the next. It is recorded rather than quietly
> patched. **The band is `[54, 200]`.**

* **`B = 201` is +9.1σ on the null.** The corrected barrier is therefore *much* stronger than AB
  published: the cheapest non-vacuous ceiling a real machine can prove is far further out.
* The **dead band widens from 41 wide to 148 wide** at 2^30. For any realistic machine, essentially
  **no** `B` in the interesting range admits a certificate cheaper than solving.
* **Caution, stated because it is exactly the trap this fleet keeps falling into:** AB's *struck*
  round-1 break-even was **198**, and my memory-aware answer is **201**. These are numerically close
  for entirely unrelated reasons (round 1 was a wrong time-only cost model; 201 is a correct
  memory-aware one). **This is a coincidence and must not be read as vindicating the struck claim.**

**Re-check of AB's S3 reach table** (my code, same model): **all 35 cells MATCH**, including
`w ≤ 14` at 2^47/2^30 and `w ≤ 52` at rho/2^30. That part of round 3 is sound.

### 3.1 AB's §5/§8 payoff table is still running the *retracted* round-1 model

`UPPER_BOUND_MAP.md` §8's table (`w ≤ 20 → 2^58.0`, `w ≤ 56 → 2^123.8`, `w ≥ 60 → 2^129.7 worthless`)
uses `C(256,B/2)`, which AB retracted in round 2 and never propagated here. Corrected:

| B | AB §8 (retracted model) | time-only corrected | **M = 2^30** | M = 2^40 | verdict at 2^30 |
|---|---|---|---|---|---|
| 14 | 2^43.6 | 2^38.8 | **2^43.4** | 2^38.8 | beats rho |
| 20 | 2^58.0 | 2^50.3 | **2^60.5** | 2^55.5 | beats rho |
| 30 | 2^79.1 | 2^66.4 | **2^84.6** | 2^79.6 | beats rho |
| 40 | 2^97.8 | 2^79.8 | **2^104.6** | 2^99.6 | beats rho |
| 52 | 2^117.7 | 2^93.1 | **2^124.5** | 2^119.5 | beats rho |
| 56 | 2^123.8 | 2^97.0 | **2^130.2** | 2^125.2 | **worthless** |
| 60 | 2^129.7 | 2^100.6 | **2^135.6** | 2^130.6 | **worthless** |

AB's §8 *conclusion* ("only `B ≲ 56` changes anything") is right **memory-aware** while every number
in its table is wrong in the time-only model it claims to be using — the third instance in this
document of two errors cancelling. The corrected column is the one to plan against.

---

## 4. ARITHMETIC AUDIT — three further defects, all conservative

### 4.1 `rep(W)` is exactly 2× too large for every odd `W`

AB's `rep(W) = C(256,128)/(C(W,⌈W/2⌉)·C(256−W,128−⌈W/2⌉))` counts only the split
`(⌈W/2⌉, ⌊W/2⌋)`. For odd `W` the split `(⌊W/2⌋, ⌈W/2⌉)` is equally admissible — both sides are
`≤ ⌈W/2⌉`, which is the only constraint the half-lists impose.

*Proof.* Write `W = 2c−1`. The two admissible terms are `C(W,c)C(256−W,128−c)` and
`C(W,c−1)C(256−W,129−c)`. `C(W,c) = C(W,c−1)` since `W = 2c−1`; and `256−W = 257−2c` is odd with
`(128−c)+(129−c) = 257−2c`, so the two binomials are the twin central ones of an odd-`n` row and are
**equal**. Hence the honest denominator is exactly twice AB's. ∎

Verified for all `W = 1..255` at `n = 256`: **0 exceptions** (`ag_verify.py` T2). Independently,
`1/rep_honest` was checked against a **measured** balance frequency at `n = 32` over 200 000 trials
per `W`, max relative error 0.0087 (T1) — a test that could have failed and did not.

**Effect:** crossover `106 → 109`, break-even `148 → 145`. AB's error **overcharges the attacker**,
i.e. it errs toward overstating the barrier — the same direction as the Theorem-D constant Z caught.

### 4.2 The `W = 256` self-certificate is asymmetric: it can refute, it cannot confirm

> ~~"It is structurally incapable of failing."~~ **AMENDED (round 2) — AB's precision is fair and I
> adopt it.** The certificate **did** fail, informatively, against round 2's model, which had
> `rep(256) = 16` and returned 2^132.0. It is not vacuous.

The correct statement is asymmetric:

> The `W = 256` check **can refute any model with `rep(256) ≠ 1`** (round 2's, and it did), but it
> **cannot confirm a model with `rep(256) = 1`**, because every model of the form
> `rep(W)·Vol₁₂₈(⌈W/2⌉)` with `rep(256) = 1` returns 2^128 identically — `Vol₁₂₈(128) = 2^128`.

AB used it in the second, unavailable direction: as the certificate that round 3's model is *right*.
In that direction it is blind to §4.1's odd-`W` error (256 is even) and to the floor/ceil family
generally. A check that *can* fail in the confirming direction is §1.3's comparison against `√Z`,
which is a different function of `W` entirely — and §4.4 above is what happens when a boundary is
sampled rather than scanned.

### 4.3 A zero-error covering proof costs more than AB charges

`rep(W)` is the **expected** number of *random* splits until one balances the (unknown) true flip
set — correct for a **Las Vegas** search. A *proof* must be zero-error, and needs a genuine
**splitting system**: a fixed family balancing **every** `W`-subset.

Such a system exists and is small: arrange the 256 positions on a cycle and take the contiguous
windows of length 128. For a fixed `W`-set `D`, `f(i) = |D ∩ window_i|` changes by at most 1 per step
and satisfies `f(i) + f(i+128) = W`, so by a discrete intermediate-value argument some window holds
exactly `⌊W/2⌋` or `⌈W/2⌉`. **Exhaustively verified over every `W`-subset for `n = 12, 14, 16`**
(`ag_verify.py` T3) — a test that could have failed.

> ~~"128 windows"~~ → **129. AB is right.** The IVT path must be evaluated at
> `i = 0 … 128`, which is **129 window positions**; those induce **128 distinct partitions**, since
> window `i` and window `i+128` are complementary. Quote 129 evaluations / 128 partitions.

So the deterministic factor is `≤ 128` against AB's `rep ≈ 10`: **AB underprices a zero-error proof
by at most 2^3.7.** Again conservative for the negative conclusion, but it means the true cost sits
in `[√Z, 2^3.7·AB]`, a total modelling band of about 2^6.4 — which is still far too small to matter.

### 4.4 The covering curve is FLAT below the break-even — there is no cliff at 148

> ~~"The minimising `W` for every `B ≤ 148` is the saturated one (one ball = the whole space), so
> **`cover(B) = 2^128.000` exactly for every `B ≤ 148`**."~~
> **STRUCK (round 2). AB is right and my own §1.3 table already said so.** The largest `B` with
> `cover(B) = 2^128.000` exactly is **142**, not 148. `cover(148) = 2^126.854` — the very value
> printed in §1.3 — and the minimiser there is `W = 106`, not the saturated radius.
>
> **How I got it wrong, since it is the same failure mode I charged AB with in §4.2.** My scan
> printed `B = 148` and then jumped to `B = 140`; I saw 2^128.000 at 140 and below and generalised
> over the six values 143–147 that the sample skipped. **A sampled scan hiding a boundary is
> exactly AB's even-`W`-only scan hiding the floor/ceil bug.** I criticised it and then did it.

**Corrected, recomputed over every `B` (not a sample):**

| B | 143 | 144 | 145 | 146 | 147 | **148** | 149 |
|---|---|---|---|---|---|---|---|
| `cover(B)` | 2^127.882 | 2^127.851 | 2^127.401 | 2^127.373 | 2^126.881 | **2^126.854** | 2^126.320 |
| minimiser `W` | 112 | 110 | 110 | 108 | 108 | 106 | 106 |

`cover(B)` is non-increasing in `B` on `[0,148]` (verified for all 149 values), so

> **`cover(B) ∈ [2^126.854, 2^128.000]` for every `B ≤ 148` — a band of `2^1.146`; quoted against
> rho as the floor, `[2^126.533, 2^128.000]`, a band of `2^1.467` (AB's phrasing, and the one to
> use, since rho is the thing being compared against).**

**The substance is unaffected.** Proving `w ≤ 0` (producing `k₀`) is only 2^1.15 dearer than proving
`w ≤ 148`; the word "break-even" still invites the reading that something changes at 148, and
nothing does. The conclusion "no cliff" stands; the exact-value phrasing does not.

### 4.5 Disk — conclusion stands, arithmetic off by ~2^5, premise off by 3×

AB: *"~30 GB ≈ 2^31 entries; vOW is random-access so disk runs at seek rate ~10²/s vs memory
~10⁸/s: a 2^20 slowdown for a 2^1 memory gain."*

Measured here (`ag_disk.py`: 200 MB scratch file, `O_DIRECT` so the page cache is bypassed, 4000
random 4 KiB reads, box under load 16.6, file deleted afterwards):

* **4.92 × 10³ random 4 KiB reads/s**, not 10²/s. AB is **2^5.6 too pessimistic** on the seek rate.
* Honest slowdown against a C-speed in-RAM hash table (10⁸–10⁹ probes/s): **2^14.3 – 2^17.6**, not 2^20.
* **The premise is also wrong: there are 11.5 GB free now, not ~30 GB** (`df`), and the shared tables
  may not be deleted. At 32 B per vOW distinguished-point entry that is 2^28.4 entries — **fewer
  than the 2^30 RAM figure.** Disk currently buys **negative** memory.

**Conclusion "disk is not a way out" therefore stands and in fact hardens** (vOW gains only `√M`, so
a ≤2^1 memory gain against a ≥2^14 slowdown is a catastrophic trade), but neither number in AB's
sentence is right and the sentence should be replaced.

---

## 5. WHICH HYPOTHESES ARE LOAD-BEARING

The coordinator asked precisely this. Dropping one at a time:

| hypothesis | drop it and… |
|---|---|
| **zero-error / correct for every `w`** | **THEOREM FALSE.** Cost collapses to 2^30 at `B = 20` via the low-weight sweep (§2.1). **Most load-bearing hypothesis in the theorem, and it is not stated.** |
| **unbounded memory** | Numbers move hard, conclusion **strengthens**: break-even 148 → 201, dead band 41 → 148 wide (§3). |
| **"search-based" = ball covering** | Generalise to rectangles: **no change** (≤2^2.66, §1). Generalise to intervals: classification breaks, conclusion survives via success probability (§2.2). This hypothesis is **cosmetic**, and can be dropped in favour of the rectangle model at a profit. |
| **Hamming structure only (no group-structure exploitation)** | **REAL GAP.** Neither B nor D covers a non-generic algebraic certificate (§2.3). The only evidence against that route is `d_reg` measured at `n = 2, 3`. |
| **rho at 2^126.5 as the comparison** | Everything shifts by the same amount; this is the ECDLP assumption itself and is not Theorem B's business. |
| **`rep(W)` exact** | 2× error for odd `W`; moves 148 → 145 (§4.1). Not load-bearing. |
| **Las Vegas vs deterministic splitting** | ≤2^3.7 the other way (§4.3). Not load-bearing. |

---

## 6. RECOMMENDED RESTATEMENT

> **THEOREM B (AG restatement).** *Model: any meet-in-the-middle certifier, i.e. any covering of
> `{(k_H,k_L) : wt(k_H)+wt(k_L) > B}` by combinatorial rectangles over any splits of the 256
> exponent positions, with a rectangle of area `z` costing `≥ √z`; memory `M`; comparison point rho
> at 2^126.53.*
>
> 1. **(Floor.)** Any such certifier costs `≥ √|{wt>B}|`, independently of overlap, code structure
>    or shared lists; the best known construction is within **2^2.66** of it.
> 2. **(Saturation.)** With unbounded memory `cover(B) ∈ [2^126.854, 2^128.000]` for every
>    `B ≤ 148`, saturating at exactly 2^128.000 from `B = 142` down — within 2^1.47 of solving
>    outright, which returns `w` exactly. **There is no cliff at 148.**
> 3. **(Memory.)** With `M = 2^30` the break-even is `B = 201` (+9.1σ) and the crossover `w = 53`;
>    no `B ∈ [54,200]` admits a certificate cheaper than solving.
> 4. **(Trichotomy.)** Every search-based upper bound is vacuous, or priced at solving, **or** a
>    zero-error decider of `[w ≤ B]` that is cheap for `B` far from 128 — and which returns an
>    *upper* bound only with the null probability that the bound holds, in which branch it has
>    produced `k₀`. *(Zero-error, per AB: the cheap branch is a decider, not merely one-sided.)*
> 5. **(Unified law.)** For every `B`, deciding `[w ≤ B]` costs
>    `Θ(√min(|{w≤B}|,|{w>B}|))`, matching Theorem D's generic bound at **every** `B`.
>
> **Not covered by this theorem or by Theorem D:** non-generic certificates that work in the
> coordinate ring of `E/F_p` rather than through the group law or the Hamming metric.

---

## 6a. ROUND 2 — adjudication against AB

AB replied item by item. **Three concessions to me in full, one with a precision, one rebuttal
against me that lands.** Everything below was re-run in this directory before being recorded.

### What I withdraw

| my claim | status |
|---|---|
| §4.4 "`cover(B) = 2^128.000` exactly for every `B ≤ 148`" | **STRUCK.** Largest such `B` is **142**; `cover(148) = 2^126.854`, which my own §1.3 table printed. Corrected band `[2^126.854, 2^128.000]`. Conclusion "no cliff at 148" survives. |
| §4.2 "the `W=256` certificate cannot fail" | **AMENDED.** It refuted round 2's model (`rep(256)=16` → 2^132.0). Correct form: **can refute `rep(256) ≠ 1`, cannot confirm `rep(256) = 1`** — and AB used it in the second direction. |
| §3 crossover **52** at `M = 2^30` | **CORRECTED to 53**, band `[54,200]`. My own odd-`W` fix, which I failed to propagate into my own table. |
| §4.3 "128 windows" | **CORRECTED to 129 window positions** (128 distinct partitions). |

Two of the four are the *same* error class I convicted AB of — a sampled scan hiding a boundary
(§4.4) and a fix not propagated one section later (§3). Recorded, not patched over.

### What AB conceded, and what it independently confirmed

* **Memory (Attack 3): conceded in full** — AB: *"I struck the unbounded-memory column in S3 and
  then quoted `106/148` — values from that column — as the headline one section later."*
* **Quantifier (Attack 2): conceded, and AB says my concession was too small.** Exhausting
  `{w ≤ B}` is not merely one-sided — it is a **zero-error decider** (hit ⇒ `w ≤ B`, miss ⇒
  `w > B`). **So the 2^77.7 overprice at `B = 20` is not attributable to one-sidedness at all.**
  AB is right; my §2.1 framing understated my own result. The zero-error hypothesis is now inside
  AB's statement, together with my trichotomy and my `√min(·,·)` law.
* **Disk: conceded cleanly**, conclusion hardens.
* **`rep(W)` odd-`W` factor 2: reproduced exactly** from AB's independent code (106 → 109, 148 → 145).
* **My no-carry step was checked first and then generalised by AB**: disjoint bit support makes
  `wt(k) = wt(k_H)+wt(k_L)` hold for **any** split, contiguous or not — 20,000 random
  (split, `S`) pairs including non-contiguous ones, 0 failures. My rectangle model is therefore not
  tied to the `2^128` split.
* **My `√Z` floor re-derived by AB via my generic-query route, matching to 0.01 bits at four values
  of `Z`. Attack 1 fails, and I reported it as failing.**

### Where this landed

Agent AC's exact posterior on `w` gives a 90 % interval `[115, 141]`; the coordinator's join with
the corrected band is **`P(w ∉ [54,200]) = 2^-67.327`**. My 52-vs-53 dispute moves that by 0.07 bits.

> **With probability `1 − 2^-67.3`, `w` lies inside the band where — by the corrected theorem — no
> upper bound on `w` is cheaper than solving the instance outright.**

That statement is the product of Attack 3, and it is what the memory correction was worth. It is
also, note, a statement whose force comes from a *prior* (AC's posterior), not from a proof about
this instance: it is not an infeasibility claim, and §8 remains open and untouched by any of this.

### The remaining gap, endorsed by AB as the sharpest stated so far

> A **non-generic algebraic certificate** is missed by Theorem B (it is not a covering) and by
> Theorem D (which excludes the encoding). **Certificate size is never the barrier — `k₀` is itself
> a 256-bit certificate verifiable by one scalar multiplication — so only *finding* cost can be.**

That is now the `d_reg` question; `n = 4` is computing under agent AI's custody with AB's read-off
written in advance. **Nothing in this document bears on it**, and it is the one hypothesis of
Theorem B (§5, row 4) whose removal is not covered by anything the fleet has proved.

---

## 7. FILES

`ag_recompute.py` (arithmetic audit, `rep`, volumes, curve, crossover, break-even) ·
`ag_attack.py` (the three attacks, memory-aware tables) ·
`ag_verify.py` (T1–T5 falsifiable checks) ·
`ag_disk.py` (measured random-read rate; scratch file deleted).
No process launched; no shared table touched; ~0 MB persistent footprint.
