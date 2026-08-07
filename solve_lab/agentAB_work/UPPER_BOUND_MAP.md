# UPPER_BOUND_MAP — every mechanism that could bound `w` FROM ABOVE

Agent AB. Theory task, not search. Every number below was recomputed in this directory
(`ab_facts.py`, `ab_cost.py`, `ab_rank.py`, and round 2's `ab_barrier.py`, `ab_costfix.py`,
`ab_soft.py`); nothing is quoted from another agent without
re-derivation, except where explicitly attributed as *their measurement*.

**Verified before starting:** `python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
→ `satisfied 39026/39033 (7 failing) [12231,12270,12350,14584,18673,22044,29125]`.


> ## ⚠ CORRECTION NOTICE — READ BEFORE §10
>
> **Theorem B's numeric table in §10 was WRONG and is retracted below in §13.0. The qualitative
> conclusion survives; the numbers do not.** I priced a radius-`W` ball at `C(256, W/2)`. The correct
> per-ball cost is the *cumulative* half-volume `Vol₁₂₈(W/2) = Σ_{j≤W/2} C(128,j)`, which is smaller
> — by up to **2^65** at large `W`. The break-even moves from **`B = 198` to `B = 148`**. This is the
> lab's documented failure mode (a number computed under one configuration reported as a property),
> committed by me, and the coordinator has already circulated the wrong figure. **§13.0 is the
> retraction; §13.1–13.2 are the new barrier that replaces it.**

---

## 0. The object, restated exactly (this reframing does work later)

From `agentX_work/xdata.json`, re-verified with my own field/curve arithmetic (`ab_facts.py` §0):

* `p = 2^256 − 2^32 − 977` ✔, `N = 0xFFFF…364141` (secp256k1 order) ✔, both prime (sympy) ✔
* `a = 0`, so `j = 0`, `Aut = μ₆`, CM by `√−3` (`(4p − t²)/3` is a perfect square ✔)
* `b/7` **is a 6th power mod p**, so this curve is `F_p`-**isomorphic** to secp256k1 `y² = x³ + 7` —
  not merely "a sextic twist of the same order". Recomputed here.
* `ladder[i] == 2^i·G` for **all** `i = 0..255` by an independent doubling chain ✔; the exponent
  set is exactly `{0,…,255}` with no gaps.
* `N·G = O`, `N·T = O`, `G`,`T` on curve ✔.

Now the restatement that matters:

> `k₀ := log_G T ∈ [0, N)`. Since `k₀ < N < 2^256`, the set `S = {i : bit_i(k₀) = 1}` is **always**
> a legal subset. **The instance is therefore always satisfiable, and `w = popcount(k₀)`.**
> A second solution `k₀ + N` exists iff `k₀ < 2^256 − N = 2^128.35`, probability `2^−127.7`.

Two consequences used repeatedly below:

1. **`bit_i(k₀)` *is* the predicate `i ∈ S`.** "Extract one bit" and "decide one membership" are the
   same operation. 256 of them reconstruct `k₀` outright.
2. **There is no existence question**, only a search question. Every statement about `w` is a
   statement about `popcount` of one specific unknown integer.

**Null:** if `k₀` is uniform, `w ~ Bin(256, ½)`, mean 128, σ = 8; `P(w ∈ [104,152]) = 0.99787`
(recomputed exactly).

---

## 1. The complement identity — Agent Y's mechanism

**Mechanism.** `fold(S) + fold(S̄) = 2^256 − 1` as integers, hence
`T' := (2^256 − 1)·G − T` has `log_G T' ≡ 2^256 − 1 − k`. Exhaust weight `≤ W` on `T'`; a miss
proves `256 − w > W`, i.e. **`w ≤ 255 − W`**.

**Algebra checked, not assumed** (`ab_facts.py` §1):

* 200 random `S`: `fold(S) + fold(S̄) = 2^256 − 1` and `fold(S)G + fold(S̄)G = (2^256−1)G` — 0 failures.
* End-to-end planted test: `k` of weight **250** planted, `(2^256−1)G − kG` recovered as exactly
  `(2^256−1−k)·G` with complement weight **6**. The mechanism recovers it.
* `T'` computed and on-curve, `N·T' = O`:
  `x = 34393883340176920870250405895813312293662019612542515124677151347541157631736`
  `y = 113211777249390963896039371927650543610689519859184504098294945786752299991906`

**Soundness under the mod-N wrap** (the failure mode worth naming): the true complement `S̄` is a
genuine subset of `{0..255}`, so an exhaustive weight-`≤W` sweep on `T'` *would* find it. The wrap
can only manufacture **false positives** (some `S'` with `fold(S') ≡ log T'` but `≠` as integers),
never false negatives. **The one-sided inference "miss ⟹ `w ≤ 255 − W`" is valid in every case,
including the `k = k₀ + N` branch.** ✔

**Strength, stated honestly** (`ab_rank.py`):

| W | cost `C(256,W/2)` | proves | null mass of the excluded region | σ above mean |
|---|---|---|---|---|
| 10 | 2^33.0 | `w ≤ 245` | 2^−198.0 | +14.8σ |
| 20 | 2^58.0 | `w ≤ 235` | 2^−158.0 | +13.5σ |
| 30 | 2^79.1 | `w ≤ 225` | 2^−126.1 | +12.2σ |

Information gained by excluding a region of null mass `q` is `≈ q/ln2` bits: at `W = 10` that is
**3.6×10⁻⁶⁰ bits**. Exactly symmetric to the low-weight side — `P(w ≤ W) = P(w ≥ 256−W)` exactly,
verified for `W = 7,10,12,14,16,20`.

> **VERDICT: LIVE — it is the *only* mechanism in this document that produces an unconditional upper
> bound at affordable cost — but it is vacuous against the null and is a bet on a designer
> hypothesis, precisely the mirror of the low-weight bet the fleet is already making.**

The honest framing: the complement search is worth running **not** because `w ≤ 245` is
informative, but because a *hit* solves the instance, and "designer chose `k = 2^256 − small`" or
"`k = all-ones-ish`" is about as plausible a priori as "designer chose `k` sparse". Its value is
entirely in the hit branch, identically to the low-weight sweep. **Do not report `w ≤ 245` as
progress on bounding `w`.**

---

## 2. Partial information / bit security

**Mechanism.** Extract one bit of `k` for less than full DLP cost, then iterate or use it to prune.

**Settled here, without needing to trust a citation.** By §0(1), `bit_i(k₀)` = `[i ∈ S]`, so a bit
oracle answers the instance in 256 queries: **any single-bit extractor costs ≥ 2^126.5/256 = 2^118.5**,
a self-contained reduction requiring no literature. (This bounds *exact* extractors; it says nothing
about biased predictors — see the adversarial note.)

**Why bits are easy in `Z_p^*` and not here — the actual reason.** In `Z_p^*` the order `p−1` is
even, so the index-2 subgroup gives a *computable* quadratic character (the Legendre symbol), which
is exactly the LSB of the discrete log. **`N` is prime.** The character group of a cyclic group of
prime order `N` is cyclic of order `N`: every nontrivial character has order `N`, and evaluating one
requires the DL itself. **There is no proper subgroup to project onto and no low-order character to
evaluate.** Verified: `N` prime; `E(F_p) ≅ Z/N` with `N·G = O` and no smaller order.

**Boneh–Venkatesan / HNP points the wrong way.** HNP says: an oracle giving `≈√(log N)` MSBs of
`α·tᵢ mod N` for *known random* `tᵢ` lets you recover `α` by lattice/CVP. That is a *hardness
amplification for bits given an oracle*, not a way to manufacture the oracle. To build HNP samples
here we would need MSBs of `k·tᵢ mod N`; from `T` we can compute `tᵢ·T = (k tᵢ)·G` but never any
integer digits of `k tᵢ`. The reduction is unusable in this direction. Boneh–Shparlinski (2001) is
about the bit security of the *ECDH shared value's x-coordinate*, not of the discrete log, and also
assumes an oracle.

**Adversarial note on my own verdict.** The reduction above kills *exact* bit oracles. It does not
kill a *biased predictor* — an algorithm computing some predicate `P(k)` with advantage ε at low
cost. Random self-reducibility (`k ↦ k + r`) would normally convert a biased predictor into an exact
one, but the predicate we actually want, `w(k) ≤ B`, is **not compatible with that self-reduction**:
`k ↦ k + r` scrambles Hamming weight completely. So the standard hardcore machinery neither proves
nor refutes hardness of the *weight* predicate. I can say **no mechanism is known**; I cannot say
"provably hard".

> **VERDICT: DEAD as a route to an upper bound.** Exact bit extraction is DLP-hard by a 256-query
> reduction proved here; the HNP machinery runs in the opposite direction; prime order removes the
> character shortcut that makes bits easy in `Z_p^*`. *Least-secure part:* the weight predicate
> itself is outside the reach of both the reduction and the hardcore-bit theory — see §11.

---

## 3. Lattice / LLL on the subset-sum

**Density, recomputed.** `n = 256` items `{2^i}` modulo `N`: density `n / log₂N = 256/256 = 1.000`.
Low-density attacks (Lagarias–Odlyzko; Coster–Joux–LaMacchia–Odlyzko–Schnorr–Stern) require
density `< 0.9408`. **Not met, and not marginally** — this is the worst possible density.

**But the density is a red herring, and this is the settling argument.** Powers-of-two subset-sum is
*trivial* at any density: **given `k₀`, you read off `S` from the binary expansion in O(256) bit
operations**. There is no subset-sum problem here at all. **100% of the hardness is "we do not know
`k₀`."** A lattice attack needs an integer (or rational) target vector to reduce. The only integers
obtainable from `T` are `x_T, y_T ∈ F_p`; their relation to `k₀` runs through division polynomials
`φ_{k₀}/ψ_{k₀}²` of degree `≈ k₀²/2 ≈ 2^511`. There is **no linear or near-linear structure over ℤ
connecting anything computable from `T` to `k₀`**, so there is nothing to put in a basis.

*Knob set for this claim:* lattices whose basis is built from instance-derivable integers
(`x_T`, `y_T`, `p`, `N`, `2^i`, and integer combinations thereof) under LLL/BKZ/CVP. Within that
knob set, no formulation reaches `k₀`, because `k₀` never appears as an integer in any of them.

*Adversarial check:* lattice methods **do** finish ECDLP when partial information exists (biased
nonces, HNP-style). That is §2's problem, not §3's — the lattice is downstream of information we do
not have. Also checked and rejected: index calculus / Semaev summation polynomials / Gaudry–Diem
decomposition require a factor base from a proper subfield; `p` is prime, there is no subfield.
Weil descent/GHS requires `F_{q^n}`. Both inapplicable.

> **VERDICT: DEAD.** Not because of density, but because a lattice needs an integer target and the
> instance supplies none. The subset-sum framing contributes *nothing* to the difficulty.

---

## 4. 2-adic / valuation arguments

**Mechanism.** `v₂(k) = min(S)` = trailing-zero count of `k₀`. Does `T` constrain it?

**Killed outright, exactly as suspected.** `E(F_p)` has odd prime order `N`, so `[2]` is a
**bijection** on the group: `2^{-1} mod N` exists (verified `2·2^{-1} ≡ 1`), every point has a unique
halving, and iterating halving from `T` 40 times stays consistent (verified). Multiplying `T` by
`2^{-1}` produces `k₀/2` if `k₀` even and `(k₀+N)/2` if odd — **both are legitimate points and are
indistinguishable without knowing `k₀`.** No 2-adic invariant of `T` exists to read.

There is also no formal-group / reduction-of-a-global-point argument available: the points live in
`E(F_p)`, not `E(ℚ)`, so there is no `ℤ_p`-valued logarithm to take valuations in.

*Adversarial check:* is anything 2-adic hiding in `N − 1` (relevant to `F_N` arithmetic) or `p − 1`?
Neither enters DLP in a group of prime order `N` — Pohlig–Hellman uses the factorisation of `N`
itself, and `N` is prime.

> **VERDICT: DEAD, unconditionally. Odd group order kills it; there is no residual angle.**

---

## 5. The endomorphism `λ`

**Mechanism.** `λ³ ≡ 1 (mod N)`, `λ ≠ 1`; `[λ](x,y) = (βx, y)`. Orbit `{k, λk, λ²k}`, and with `−1`
the full `μ₆` orbit `{±k, ±λk, ±λ²k}`.

**Computed** (`ab_facts.py` §5): both cube roots found, the matching `(λ, β)` pair identified against
the actual point operation, so the endomorphism is the standard one and is real on this curve.

**Weight information: none.** For 300 random weight-4 keys `k`:

* `popcount(λk mod N)`: mean **127.41**, sd **8.31**, range [105,151] — indistinguishable from the
  `Bin(256,½)` null. **`λ` destroys sparsity completely.**
* `popcount(−k mod N)` = `popcount(N−k)`: mean **189.98**, sd **2.30**. Sharply concentrated near
  `popcount(N) = 192`, but that is a fact about `N`, not about `k` — and it goes the wrong way: a
  miss on `−T` bounds `popcount(N − k₀)` from **below**, which constrains `popcount(k₀)` not at all.

`1 + λ + λ² ≡ 0`, so `T + λT + λ²T = O` — a relation among the *targets*, giving no handle on `k`.

What `λ` actually buys: `√3` in rho (already in the 2^126.5 figure), 3× (6× with negation) target
coverage in MITM. Both are **lower-bound** machinery. GLV decomposition `k = a + bλ`, `|a|,|b| ≈ 2^128`,
defines a *different* low-complexity class (sparse `a`,`b`); sparse `(a,b)` gives a **dense** `k`, so
it says nothing about `w`. (Campaign already excluded `|a|,|b| < 2^21`.)

> **VERDICT: DEAD for upper bounds. `λ` is a √3 search accelerator and a 3× coverage multiplier and
> nothing else. Measured, not argued.**

---

## 6. Character sums / analytic equidistribution

**Mechanism.** Exponential-sum bounds (Shparlinski, Konyagin, Lange–Shparlinski on distribution of
sequences of points on elliptic curves) show that `{kG : popcount(k) = w}` equidistributes on `E(F_p)`
once `w` is not tiny.

**Why it cannot give a single-instance upper bound.** Equidistribution is a statement about the
*image measure* of a weight class. It licenses "for a **random** `T`, `popcount(log T) ~ Bin(256,½)`"
— which is **exactly the null already in use**, not a bound. Reversed, it says: any fixed `T` is hit
by `#{k : w(k) ≤ W} = 2^{198.0}` keys at `W = 10` out of `2^256`, so a random `T` fails to have small
weight — again average-case. `w` for *this* `T` is a fixed integer; no analytic statement about a
family bounds a fixed member. To conclude anything about this `T` you would need the sum restricted
to a set defined by `T`, which is the DLP.

*Adversarial check:* is there any *effective* character-sum inequality that is nonvacuous for a
single point? All such bounds are of the form `|Σ_{k∈W-class} ψ(x(kG))| ≤ B` with `B ≫ 1`; to isolate
one `T` you need the full sum's cancellation to fail, and the bounds are far too weak (`B` exceeds 1
by many orders for any class of size < 2^128). Verified by size: the smallest weight class that could
be resolved by a square-root-cancellation bound would need `|class| > √#E ≈ 2^128`, i.e. `w ≥ 44` —
and at that size the bound tells you the class is equidistributed, i.e. *nothing* about membership.

> **VERDICT: DEAD as an upper bound. LIVE only as the source of the null prior — which the campaign
> already has.** State it that way and stop.

---

## 7. Counting / uniqueness

**Mechanism.** `2^256 − N = 2^128.35` (recomputed), so for a given `T` there are 1 or 2 valid `S`.

Uniqueness holds for **every** target regardless of weight — it is a statement about the map
`S ↦ fold(S) mod N` being 2^256 → N, not about `w`. It carries zero weight information.

The only weight-relevant crumb: a **second** solution exists iff `k₀ < 2^256 − N`, which would force
`k₀ < 2^128.35`, hence `popcount(k₀) ≤ 129`. But (a) `P = 2^−127.7`, and (b) detecting two-ness
requires knowing `k₀`, i.e. solving. Circular.

**One genuinely free, genuinely unconditional bound falls out of counting**, and it is the only
non-search upper bound in this document. `k₀ < N`, and `N`'s binary form is: bits 255..129 set
(127 ones), **bit 128 = 0**, low 128 bits `= 0xBAAE…4141`; `popcount(N) = 192`, 64 zeros. An exact
digit-DP for `max{popcount(x) : 0 ≤ x < N}` (validated against brute force on 7 small moduli) gives:

> **`w = popcount(k₀) ≤ 255`, free, unconditional.** (Attained by `x` = bits 255..130 set, bit 129
> clear, bits 128..0 set.) With the `k₀ + N` branch, `w ≤ 256`.

> **VERDICT: DEAD for anything useful.** Uniqueness gives nothing; counting gives `w ≤ 255` at zero
> cost, which is 15.9σ above the null mean and worse than the complement mechanism at `W = 6`.

---

## 8. Instance-side constraints (Agent Z's angle) — **THE ONE THAT IS NOT SETTLED**

I have not duplicated Z. What follows is (a) what a positive result buys and (b) a correction to the
evidence the fleet currently believes settles it.

**What a positive result buys — the largest payoff on this page.** An instance-level constraint
`w ≤ B` is the *only* upper bound that does not cost a search (see §10's theorem). Its value is
sharply threshold-shaped:

| instance forces | MITM cost `C(256,B/2)` | verdict |
|---|---|---|
| `w ≤ 20` | 2^58.0 | **solvable on a cluster in days** |
| `w ≤ 30` | 2^79.1 | national-scale |
| `w ≤ 40` | 2^97.8 | beats rho, still infeasible |
| `w ≤ 56` | 2^123.8 | crossover with rho |
| `w ≥ 60` | 2^129.7 | worthless — rho is cheaper |

**So only a bound `B ≲ 56` changes anything at all, and only `B ≲ 24` changes anything the fleet can
act on.**

**Correction to the fleet's evidence — this matters and I believe it is currently being read
backwards.** `MINIMUM_COST_SEARCH.md` §7 argues from agent T's closures at `|S| = 2,3,5,6,7,8,17`
(all scoring 39,018 with the identical 15-equation footprint) that "no measurement of the instance
leaks any information about `w`". Re-reading T's own `RESUME_T.md`:

* Those closures, plus `|S| = 32` (`close_T32f.json`) and `|S| = 64` (`close_T64.json`), all close
  with exactly **2 nonzero atoms** — the two target congruences — and the identical 15-equation
  failing set.
* **`|S| = 128` was started and did not finish.** `agentT_work/t_close2wj_T128.log` ends at
  `outer 8: global nonzero 3, violated c-conditions 3, nonzero handle-less 0`, then
  `no addable collateral -- giving up on ((x2820-x17195)-(8271997*x17079))`. **It stalled at three
  nonzero atoms — one more than every closure — and gave up.**

Two things follow that the current synthesis does not say:

1. **All confirmed closures are at `|S| ≤ 64`, i.e. below the null mean.** They are evidence against
   a *lower*-bound-style instance constraint. **They cannot rule out an upper-bound constraint of
   the form `w ≤ B` for any `B ≥ 64`** — such a constraint is perfectly consistent with every data
   point the fleet holds.
2. The single high-`|S|` probe is **open, and it stalled in the direction an upper bound would
   predict.** T's own standard applies (it applied it to the `|S| = 8` stall, which later closed with
   the joint two-wire pass): this is *"this solver did not close it"*, not *"it cannot be closed"*.
   T's measured trend argues for a solver gap — nonzero-atom count grows as `≈|S|/2 + 2`
   (18 / 34 / 70 at `|S| = 32/64/128`) while the **handle-less count is flat at 2, 3, 3** — but the
   probe is unresolved.

> **VERDICT: LIVE, and the only live non-search mechanism. Currently UNSETTLED, with the decisive
> experiment cheap and already tooled.** Recipe in §12.

Caveat to carry into that experiment: T's `|S| = 32/64/128` ON-sets are **nested prefixes of one
`random.Random(7)` chain**, so they are one correlated sample, not three. Independent draws are
required before any conclusion either way.

---

## 9. Everything else I could find

Each of these was considered as an upper-bound mechanism and each is dead; the ones with a computed
input have it recorded.

| # | Mechanism | Result | Verdict |
|---|---|---|---|
| 9.1 | **Pohlig–Hellman** | `N` prime (verified) | DEAD |
| 9.2 | **Smart / anomalous** | `N ≠ p` (verified) | DEAD |
| 9.3 | **MOV / Frey–Rück** | mult. order of `p` mod `N` **> 20 000** (computed here, extending the lab's verified `k ≤ 24`); the MOV target field `F_{N^k}` therefore has **> 5 million bits** | DEAD |
| 9.4 | **Weil descent / GHS** | `p` prime, no subfield | DEAD |
| 9.5 | **Index calculus / Semaev / Gaudry–Diem** | needs a subfield factor base; base field is prime | DEAD |
| 9.6 | **CM by `√−3`, `j = 0`, `Aut = μ₆`** | verified `(4p−t²)/3` is a perfect square; CM curves are **not** weaker — the only cash value is GLV (§5) | DEAD |
| 9.7 | **`F_p`-isomorphism class** | `b/7` **is** a 6th power mod `p` — the instance curve *is* secp256k1 up to `F_p`-isomorphism, not just an order-matching twist. Removes any hope that the twist has extra structure | DEAD (and closes a loose end) |
| 9.8 | **Cheon's algorithm** | needs `k^d·G` for auxiliary `d | N ± 1`; instance supplies only `2^i·G` and `T`. No `k`-dependent auxiliary point exists | DEAD |
| 9.9 | **Structure of `N`'s expansion** (192 ones, 64 zeros — recomputed) | affects only `popcount(N − k)` (§5, mean 190, sd 2.3) and the free `w ≤ 255` (§7). No handle on `popcount(k)` | DEAD |
| 9.10 | **Doubling as a weight-preserving map**: `k ↦ 2k mod N` is a left shift, hence **exactly weight-preserving whenever `2k < N`** | Real but useless: searching targets `2^j·T` covers only shifted low-weight keys, which are themselves low-weight keys. Zero new coverage | DEAD |
| 9.11 | **Division polynomials / algebraic elimination** | `deg ψ_{k₀} ≈ k₀²/2 ≈ 2^511` | DEAD by size |
| 9.12 | **Gröbner / polynomial-system encoding of the ladder** | This is *literally the instance file*: 39 033 equations encoding 256 selector-gated point additions. The **entire campaign of 26 agents has been attacking exactly this system** and has produced a 39 026 plateau, `ker(M) = 0`, and a proven-tight equation-level bound of 7. Empirically tested at a scale no fresh Gröbner attempt would match | DEAD, with the strongest empirical support of anything here |
| 9.13 | **Interval / kangaroo bounds on `k₀`** (`k₀ < 2^t ⟹ w ≤ t`) | BSGS decides `k₀ ∈ [a, a+L)` in `√L`, but a **miss** gives a *lower* bound on `k₀`, which bounds `w` not at all; a **hit** *is* the solution. Upper bounds via intervals require solving | DEAD |
| 9.14 | **Multi-target / batch DLP** | one target; `μ₆` orbit gives 6, already counted as √6 in rho | DEAD |
| 9.15 | **Shor / quantum** | would give `k₀` and hence `w` exactly, in poly time | LIVE in principle, **no hardware exists**; excluded on availability, not on mathematics |
| 9.16 | **Masked complement `c_A`, `A ⊊ {0..255}`** | see §10 | DEAD (conditional bounds only) |

---

## 10. Two structural theorems I can prove, which organise the whole page

### Theorem A (the complement is the *unique* unconditional shifted centre)

For a search at centre `c` (i.e. on target `cG − T`) to yield an unconditional upper bound on `w`,
we need `wt(c − k)` to determine `wt(k)` downward for **all** `k`. Now `c − k = c ⊕ k` **iff there
are no borrows iff `supp(k) ⊆ supp(c)`**; and then `wt(c−k) = wt(c) − wt(k)`. Requiring this for all
`k` forces `supp(c) = {0,…,255}`, i.e. **`c = 2^256 − 1`, uniquely.**

Measured confirmation over 4000 random `k` (`ab_cost.py` §B), `corr(wt(k), wt((c−k) mod N))`:

| `c` | correlation | exact `wt(k)+wt(c−k)=256`? |
|---|---|---|
| `2^256−1` | **−1.0000** | **yes** |
| `N−1` | −0.4954 | no |
| `2^255` | −0.7073 | no |
| `2^128−1` | −0.5048 | no |
| random | +0.0386 | no |
| `0` (plain search) | −0.4964 | no |

For a proper mask `A`, the relation holds **conditionally on `supp(k) ⊆ A`** (verified 2000/2000 at
`|A| = 255, 224, 192, 128`) — but `P(supp(k) ⊆ A) = 2^{|A|−256}`, i.e. `2^−32` at `|A| = 224` and
`2^−128` at `|A| = 128`. **Conditional bounds on an event of null mass 2^−32 are worthless.**

**Consequence: agent Y's mechanism has no cheaper or stronger sibling. The complement is it.**

### Theorem B (the cost floor for *any* search-based upper bound) — the most useful number here

*Knob set:* MITM over a Hamming ball `B(c, W)` around any **known** centre `c`. This is fully
general for search: `k = c + Σ_{i∈D} ε_i 2^i` with signs `ε_i` fixed by `c`, `|D| ≤ W`, cost
`C(256, W/2)` time and memory. `c = 0` is the plain low-weight sweep; `c = 2^256−1` is the complement
sweep; general `c` is a signed-digit sweep about `c`.

To **prove** `w ≤ B` you must exclude every `k` with `wt(k) > B`, i.e. **cover** that set with balls.
`#balls ≥ |{wt > B}| / Vol(256,W)`. Total cost ≥ `(|{wt>B}| / Vol(W)) · C(256,W/2)`, minimised over
`W` (`ab_cost.py` §A):

| prove | `|{wt>B}|` | best `W` | #balls | cost/ball | **TOTAL** | vs solving (2^126.5) |
|---|---|---|---|---|---|---|
| `w ≤ 245` | 2^58.0 | 10 | 2^0.0 | 2^33.0 | **2^33.0** | affordable |
| `w ≤ 240` | 2^79.2 | 14 | 2^4.0 | 2^43.6 | **2^47.6** | affordable |
| `w ≤ 220` | 2^143.8 | 34 | 2^2.7 | 2^86.9 | **2^89.6** | affordable |
| `w ≤ 200` | 2^188.6 | 54 | 2^1.9 | 2^120.8 | **2^122.7** | affordable, barely |
| `w ≤ 180` | 2^219.9 | 74 | 2^1.3 | 2^148.8 | **2^150.0** | 2^23.5× worse than solving |
| `w ≤ 152` (+3σ) | 2^246.1 | 76 | 2^24.9 | 2^151.3 | **2^176.2** | 2^49.7× worse |
| `w ≤ 128` (null median) | 2^254.9 | 76 | 2^33.7 | 2^151.3 | **2^185.0** | 2^58.5× worse |

Break-even, computed exactly: **`B = 198`**. `w ≤ 198` costs 2^126.3; `w ≤ 197` costs 2^126.8.

> **THEOREM B. No search-based upper bound below `w ≤ 198` is cheaper than solving the instance
> outright by rho — and solving returns `w` exactly. Every search-based upper bound is therefore
> either vacuous (`B ≥ 198`, `≥ 8.75σ` above the null mean) or strictly dominated.**

This is what makes §8 the only thing that matters: an instance-side constraint is the **only** route
to an upper bound that escapes Theorem B, because it is not a search.

---

## 11. Where my own DEAD verdicts are least secure

Ranked, most-suspect first. This lab has retracted five barriers; these are the three of mine most
likely to become the sixth.

1. **§2, the *weight* predicate specifically.** My 256-query reduction kills exact **bit** oracles.
   It does **not** cover a low-cost algorithm that decides `w ≤ B` directly, and the usual
   hardcore-bit machinery cannot be applied to it either, because the weight predicate is not
   preserved by the DLP's random self-reduction `k ↦ k + r`. I claim **"no mechanism known"**, and I
   am explicitly *not* claiming "provably hard". If any mechanism on this page is wrong, it is this
   one. **Nobody in the fleet has looked for a direct weight predicate**; that gap is real.
2. **§6, character sums.** My size argument ("classes below 2^44 are too small for square-root
   cancellation to say anything") is a *heuristic* application of a square-root-cancellation
   threshold, not a theorem I derived. A genuine effective bound in the Lange–Shparlinski line,
   applied to a class defined using `x_T` itself, is the shape of thing I would want an expert to
   rule out rather than me. I still expect it to be dead — every such bound is average-case by
   construction — but my reason is weaker than in §3–§5.
3. **§9.12, "the Gröbner route is dead because the fleet already attacked it."** This is empirical,
   not mathematical, and the fleet attacked the system as *equation repair and coset decoding*, not
   as *elimination with a term order*. I would not resist someone spending a day on genuine
   elimination in Singular over the reduced core. I would resist spending a week.

Verdicts I consider **secure**: §4 (odd order makes `[2]` a bijection — this is a proof), §5 (measured,
mean 127.4/sd 8.3), §7 (arithmetic), §3 (the "no integer target exists" argument is structural),
Theorem A (proof + 4000-sample confirmation), Theorem B (exact combinatorics).

---

## 12. Ranked live mechanisms, and the one experiment to run

### Ranked by cost-to-informativeness

| rank | mechanism | cost | what it yields | why ranked here |
|---|---|---|---|---|
| **1** | **§8 instance-side constraint on `\|S\|`** — finish the high-`\|S\|` lift probe | **~10–60 min CPU per probe** (T's measured curve: ~2× CPU per doubling of `\|S\|` + fixed tail) | Either kills the last non-search upper-bound mechanism, or produces the campaign's first real one | Only mechanism escaping Theorem B; decisive in **both** directions; cheapest thing on the page by ~2^25 |
| **2** | **§1 complement search** (agent Y) | 2^33 at `W=10`, 2^43.6 at `W=14` | `w ≤ 245` / `w ≤ 241`, **or a hit that solves the instance** | Sound, cheap, unconditional — but the bound is 14.8σ out and worth ~10⁻⁶⁰ bits. **Its value is the hit branch only**, exactly like the low-weight sweep |
| **3** | **Extend agent X's signed search to exponent 256** | ~free (one more table entry) | Adds complement-sparse keys to the signed class | See below — X's table currently **cannot** reach them |
| — | §9.15 quantum | n/a | exact `w` | no hardware |

**Concrete note for agent X.** `agentX_work/xsigned.c` line 107 reads exactly **256** ladder points
into 512 signed entries `±2^i, i ∈ 0..255`. Since `2^256 ≡ 2^256 − N (mod N)` and that residue has
popcount 65, a complement-sparse key `k = 2^256 − 1 − k'` has **no short signed representation over
exponents 0..255** — X's search **does not** cover Y's class, and Y's does not cover X's. They are
complementary, not redundant; and appending a 257th point `(2^256 − N)·G` to `data_real.txt` extends
the signed class to include complement-sparse keys at `m = wt(k') + 2` terms for one extra index.
(Dedicated complement MITM at `C(256,W/2)` remains far cheaper for that class — `2^33` vs `2^44.4` —
so this is coverage, not a replacement for Y.)

### The single highest-value next experiment

> **Run agent T's `t_close2wj.py` integer-lift closer at HIGH `|S|` — `|S| = 250` first, then `192`,
> then finish the stalled `|S| = 128` — on `random.Random`-seeds *independent of the seed-7 chain*
> that produced the correlated 32/64/128 probes.**

* **Command shape:** `python3 t_close2wj.py <tag> <n>` (arg 2 = `n` or an explicit comma-list of
  selectors; `outer_max` is arg 3, default 16). Owner is agent T's tooling; agent Z's angle.
  Output goes to `agentT_work/close_<tag>.json` + `t_close2wj_<tag>.log`; every result must be
  confirmed with `checker.py` (expect **39 018 / 39 033** with the 15-equation footprint).
* **Why `|S| = 250` and not `128`:** `250` is the *complement regime* — it is the direct test of the
  upper-bound question, and it is the regime in which **not one data point exists**. Every confirmed
  closure in the campaign is at `|S| ≤ 64`, i.e. below the null mean.
* **Both outcomes are decisive:**
  * **It closes** → the instance imposes no upper bound on `|S|`; §8 dies; combined with Theorem B,
    **no affordable upper bound on `w` exists by any mechanism on this page**, and the campaign's
    entire remaining lever is the hit branch of §1 and the low-weight sweeps. That is a clean, final
    answer to the question the fleet is spending compute on.
  * **It does not close, on independent seeds, after the joint two-wire pass** → the first
    instance-side evidence for an upper bound on `w`, and the highest-value finding available to this
    campaign. It must then be pushed down (`|S| = 200, 160, 128, 96`) to find where closure breaks,
    because only `B ≲ 56` has cash value (§8's table).
* **Cost:** T measured ~4 min CPU at `|S| = 32` and ~6 min at `|S| = 64`; budget ~30–60 min CPU per
  high-`|S|` probe. **This is ~2^25 times cheaper than the next-cheapest mechanism on the page.**
* **Discipline for it:** one process at a time, PID recorded, and the correlated-seed caveat above
  respected — a single ON-set that fails to close proves nothing about the instance.

### Standing caveat, stated as required

Nothing in this document is an infeasibility claim about the instance, and none should be drawn from
it. Every "nothing can do X" above carries its knob set explicitly (§3 lattice bases from
instance-derivable integers; §10 Theorem B MITM over Hamming balls with known centres; §2 exact bit
oracles, *not* biased weight predicates). `w` remains unknown, with no non-vacuous upper bound
established by anyone, and §8 unsettled.


---

# ROUND 2 — the §2 barrier, and the retraction of my own Theorem B

Round-2 scripts: `ab_barrier.py` (Theorems C and D), `ab_costfix.py` (the retraction),
`ab_soft.py` (§6 and §9.12).

## 13.0 RETRACTION — Theorem B's cost model was wrong

**What was wrong.** §10 priced a MITM over a Hamming ball of radius `W` at `C(256, W/2)`. Two errors:

1. The standard low-Hamming-weight DLP algorithm (Coppersmith; **Stinson's splitting systems**)
   splits the **256 positions** into halves and enumerates half-patterns, giving a half-list of size
   `Σ_{j≤W/2} C(128,j) = Vol₁₂₈(W/2)`, with a `√W`-ish factor of repeated partitions. That is
   cheaper than `C(256,W/2)` — the naive split is redundant by `C(W,W/2)`.
2. More seriously, for large `W` the list is the **cumulative volume**, not the single binomial
   `C(128,W/2)`. That is where the 2^65 lived.

**Sanity check the round-1 model failed.** At `W = 256` a single ball is the whole key space, so the
cost must be the birthday bound `2^128`. Round-1's model returned `C(256,128) = 2^251.7`. The
corrected model returns **2^132.0** (`√256 · Vol₁₂₈(128) = 16 · 2^128`) ✔. *I should have run that
check before publishing the table.*

**Corrected Theorem B** (`ab_costfix.py`; `ball_cost(W) = √W · Vol₁₂₈(W/2)`):

| prove | `|{wt>B}|` | best `W` | **corrected total** | round-1 said | Theorem-D generic LB |
|---|---|---|---|---|---|
| `w ≤ 245` | 2^58.0 | 10 | **2^29.7** | 2^33.0 | 2^29.0 |
| `w ≤ 230` | 2^114.7 | 24 | **2^60.1** | 2^70.0 | 2^57.4 |
| `w ≤ 220` | 2^143.8 | 34 | **2^74.5** | 2^89.6 | 2^71.9 |
| `w ≤ 200` | 2^188.6 | 54 | **2^96.8** | 2^122.7 | 2^94.3 |
| `w ≤ 180` | 2^219.9 | 74 | **2^112.5** | 2^150.0 | 2^110.0 |
| `w ≤ 152` | 2^246.1 | 102 | **2^125.7** | 2^176.2 | 2^123.1 |
| `w ≤ 128` | 2^254.9 | 126 | **2^130.5** | 2^185.0 | 2^127.5 |

**Corrected break-even: `B = 148`, not 198.** The largest affordable complement radius is `W = 107`
(cost 2^126.4 < 2^126.5), proving `w ≤ 148`.

**What survives, and it is the whole qualitative claim.** `w ≤ 148` is `+2.6σ` on the null — still
above the mean, still excluding only ≈0.4% of the null mass (≈0.006 bits), and costing 2^126.4 when
**solving outright costs 2^126.5 and returns `w` exactly**. The corrected table also shows *why*:
the covering optimum degenerates toward "one ball = the whole space = solve it."

> **THEOREM B, CORRECTED. No search-based upper bound below `w ≤ 148` is cheaper than solving the
> instance, and `w ≤ 148` itself costs 99.5% of what solving costs while returning ~0.006 bits
> instead of the answer. Every search-based upper bound is vacuous or dominated.**

**A consequence the fleet must act on — the campaign's MITM cost table is pessimistic.**
`MINIMUM_COST_SEARCH.md` §4 uses `C(256,w/2)` throughout. Corrected, **time-only, unbounded memory**:

| budget | campaign's table | **corrected** | half-list memory |
|---|---|---|---|
| 2^47 | `w ≤ 14` | **`w ≤ 18`** | 2^44.2 entries |
| 2^58 | `w ≤ 20` | **`w ≤ 24`** | 2^54.6 |
| 2^70 | `w ≤ 24` | **`w ≤ 32`** | 2^66.6 |
| 2^80 | `w ≤ 30` | **`w ≤ 40`** | 2^77.0 |
| 2^126 | `w ≤ 56` (rho crossover) | **`w ≤ 104`** | 2^122.4 |

**The rho crossover moves from `w ≈ 56` to `w ≈ 104`** — which nearly doubles the payoff band for §8:
an instance-side constraint would only have to reach **104**, three σ *below* the null mean, rather
than 56, to make MITM beat rho.

**Two caveats I am stating rather than burying, having just been burned by not stating one.**
(i) These are **time-only** figures with unbounded memory. **Memory binds hard**: on this box
~2^30 entries ⇒ `w ≤ 10`. Beyond that a low-memory variant (van Oorschot–Wiener parallel collision
search) is needed and costs *more* than the table. **The `w ≤ 104` crossover is a time-only upper
limit, not a realizable plan.** (ii) A proper memory-aware costing of low-weight MITM has not been
done by anyone in this campaign and is a concrete, cheap follow-up.

## 13.1 THEOREM C — there is no weight-preserving self-reduction. *Proved, not observed.*

The coordinator asked for the crux: additive shifts scramble weight, so is there **any**
weight-preserving group action?

**In the generic model every map on scalars an algorithm can realise is affine**: from `G` and `T` it
can produce `σ(α + βk)` for known `α, β`, and nothing else. So the question is exactly: which affine
`ψ(k) = ak + b (mod N)` preserve `popcount`?

> **THEOREM C. The only weight-preserving affine self-map of `Z_N` is the identity.**
>
> *Proof.* (1) `k = 0` forces `popcount(b) = 0`; since `0 ≤ b < N`, `b = 0`.
> (2) `k = 1` forces `popcount(a) = 1`, so `a = 2^j`, `j ∈ [0,255]`.
> (3) For any `j ≥ 1` take `k = 2^{256−j}`. Then `k ≤ 2^255 < N` and `popcount(k) = 1`, while
> `a·k = 2^256 mod N`, and **`popcount(2^256 mod N) = 65 ≠ 1`**. So `j = 0` and `a = 1`. ∎

Step (3) verified for **all 255** values of `j` — no exceptions (`ab_barrier.py`).

**No approximate version either** (20 000 random `k`, measured):

| `a` | `Pr[wt(ak+b) = wt(k)]` | `corr(wt(k), wt(ak+b))` | `corr` of the predicate `[w ≤ 128]` |
|---|---|---|---|
| `1` | 1.0000 | +1.0000 | +1.0000 |
| `2` | 0.5212 | +0.7451 | +0.6643 |
| `4` | 0.2893 | +0.6252 | +0.5002 |
| `2^8` | 0.0504 | +0.4533 | +0.2937 |
| `2^{-1}` | 0.5213 | +0.7420 | +0.6596 |
| `−1` | 0.0296 | −0.4934 | −0.3174 |
| `λ` | 0.0375 | −0.0025 | −0.0040 |
| `2`, `b=1` | 0.0240 | +0.7451 | +0.6231 |
| random `a` | 0.0350 | +0.0033 | +0.0018 |
| random `a,b` | 0.0351 | +0.0104 | +0.0021 |

*Read the baseline correctly:* two independent 256-bit integers share a popcount with probability
≈ 0.035 (visible in the random rows), so anything near 0.035 is **pure coincidence, not
preservation**. `a = 2^j` preserves exactly on `{k : 2^j k < N}`, a fraction `2^{−j}` (predicted
0.500/0.250/0.125/0.0039 for `j = 1,2,3,8`; measured 0.521/0.289/0.169/0.050 = prediction +
coincidence baseline ✔).

**Consequence.** The only maps with *any* weight preservation are `k ↦ 2^j k`, and their preservation
set shrinks like `2^{−j}` — so they cannot randomise. **There is no weight-preserving randomised
self-reduction.** Round 1 observed that additive shifts scramble weight; this proves nothing else
works. The standard hardcore-bit amplification machinery therefore has **no self-reduction to run
on** — not "we could not find one".

## 13.2 THEOREM D — the weight predicate has no generic shortcut. *This is the barrier.*

> **THEOREM D.** In the generic group model, every element an algorithm holds after `m` queries is
> `σ(α_i + β_i k)` with `(α_i, β_i)` known, and its entire view is the equality pattern among them.
> Two collide iff `(β_i − β_j)k = α_j − α_i` — **one affine equation over the field `Z_N`**
> (`N` prime), hence **at most one root**. For `k` drawn from any distribution `D`, each nontrivial
> pair collides with probability `≤ 1/|supp D|`, so
> `Pr[the view depends on k at all] ≤ m²/(2|supp D|)`. Conditioned on no collision the view — hence
> the output — is **independent of `k`**. Therefore for **any** predicate `P` and any generic
> algorithm making `m` queries,
> `Adv ≤ m² / (2·min(|D₀|,|D₁|))`, i.e. `m ≥ √(2ε·min(|D₀|,|D₁|))`.

Applied to `P = [w(k) ≤ B]` with `D₀` uniform on `{w ≤ B}` and `D₁` uniform on `{w > B}`:

| `B` | `|{w≤B}|` | `|{w>B}|` | generic lower bound (ε = ½) |
|---|---|---|---|
| 20 | 2^98.0 | 2^256.0 | `m ≥ 2^49.0` |
| 30 | 2^129.9 | 2^256.0 | `m ≥ 2^65.0` |
| 56 | 2^190.4 | 2^256.0 | `m ≥ 2^95.2` |
| **128** | 2^255.1 | 2^254.9 | **`m ≥ 2^127.5`** |
| 152 | 2^256.0 | 2^246.1 | `m ≥ 2^123.1` |
| 198 | 2^256.0 | 2^192.3 | `m ≥ 2^96.1` |
| 245 | 2^256.0 | 2^58.0 | `m ≥ 2^29.0` |

**Read both ways — this is the whole point.**

* **`B = 128` (the null median): `m ≥ 2^127.5`.** Solving costs 2^126.5. **Deciding whether `w ≤ 128`
  is, to within the automorphism speedup, exactly as hard as solving the instance.** The weight
  predicate has **no generic shortcut** — no algorithm that decides it "without producing `k`" can
  exist in this model, whatever its internal structure.
* **`B = 20`: `m ≥ 2^49` only.** Small classes really are cheap to decide, and MITM (2^50.0
  corrected) is within **2^1** of the generic optimum. **This is why every *lower* bound in this
  campaign was affordable and every *upper* bound is not: the asymmetry is class size, and it is now
  matched by a lower bound, not merely by a cost model.**

**Theorem D also repairs the gap Theorem B left.** Corrected best-known vs generic lower bound:
`w ≤ 245`: 2^29.7 vs 2^29.0 · `w ≤ 200`: 2^96.8 vs 2^94.3 · `w ≤ 128`: 2^130.5 vs 2^127.5. **The gap
is now ≤ 2^3 everywhere** (round-1's model showed a spurious 2^57 gap). Nothing large is left on the
table.

**Knob set and honest limits — stated, not buried.**
* Theorem D is a **generic group model** bound: it constrains algorithms that touch the group only
  through the group law, inversion and equality, and **not** through the `F_p` coordinate encoding.
  Non-generic attacks exist for other groups (index calculus in `F_p^*`; ECDLP over extension
  fields). For elliptic curves over **prime** fields none is known — §9.5, §3.
* It is an **average-case** bound over a distribution on `k`. For a single fixed `k` no lower bound
  of this shape can hold (an algorithm can hardwire the answer). The statement is about algorithms
  that work for a non-negligible fraction of instances.
* The HGJ/BCJ **representation technique** is the one route that could beat the corrected ball cost.
  It filters intermediate lists by `s mod M` for an auxiliary modulus `M` — and `s = k₀` is exactly
  what we do not have. Guessing `s mod M` multiplies cost by `M` and cancels the gain. **Same
  obstruction as §3, and it is why the group version is harder than the integer subset-sum version.**

> **§2 VERDICT, UPGRADED: DEAD, and now a barrier rather than a survey.** Theorem C proves no
> weight-preserving self-reduction exists; Theorem D proves no generic algorithm decides `w ≤ B` with
> advantage `ε` in fewer than `√(2ε·min(|D₀|,|D₁|))` operations, which at the only thresholds that
> would matter is the cost of solving. **The "nobody has looked" gap is closed.**

## 13.3 The other two soft verdicts, re-ranked

**§6 analytic — upgraded from heuristic to structural.** Round 1 leaned on a square-root-cancellation
threshold I had not derived. The real argument needs no threshold. Any equidistribution statement
reads `#{k ∈ W : kG ∈ A} = |W||A|/N + E`. To decide membership for one target (`|A| = 1`) you must
resolve 0 vs 1, so you need `E < 1`; every bound in this literature has `E ≥ √p ≈ 2^128` or
`E ≥ |W|^{1−δ}`. And decisively: **a *perfect* equidistribution theorem (`E = 0`) says the weight
class spreads uniformly over targets — i.e. that this `T` is exactly as likely as any other. The
stronger the equidistribution result, the less it says about one instance.** Only a *bias* could
help, and it would have to be resolvable at a single point, which is the opposite of what these
theorems assert. **DEAD, structurally.**

**§9.12 Gröbner — DEAD on complexity, independently of the campaign's empirical failure.** The
"empirical" objection was fair, so here is the arithmetic (`ab_soft.py`). Macaulay matrix width at
degree `d` over 256 boolean variables: `d=8 → 2^48.6`, `d=16 → 2^83.2`, `d=32 → 2^135.6` columns.
At `ω ≈ 2.37`, `d = 8` costs 2^115.2 (affordable), `d = 16` costs **2^197.1** (dead), `d = 24` costs
2^264.3. **So elimination is affordable only if `d_reg ≲ 11`**, while the standard semi-regular
estimate for a boolean system with a unique solution puts `d_reg` at a constant fraction of `n`
(≈ 0.09·256 ≈ 23 for the quadratic-over-`F_2` analogue, larger for the mixed `F_p`/boolean system
actually present). **And a term order cannot be the missing ingredient: `d_reg` is a property of the
ideal, not of the order — a different order changes the constant, not the exponent.**
*Worth exactly one bounded experiment:* degree-truncated linear algebra in Singular on the reduced
core to **measure** `d_reg` rather than estimate it. Cheap, replaces an estimate with a number. **Not
worth a week**, and it cannot change the verdict unless `d_reg ≤ 11`, which would itself be the
finding.

## 13.4 Where round 2 leaves the map

| verdict | status after round 2 |
|---|---|
| §2 bit security / **weight predicate** | **DEAD, upgraded to a barrier** (Theorems C + D) — was my least-secure verdict, is now my most secure after §4 |
| §6 analytic | DEAD, upgraded from heuristic to structural |
| §9.12 Gröbner | DEAD on complexity; one cheap `d_reg` measurement would close it fully |
| §10 Theorem B numbers | **RETRACTED and corrected** (§13.0); qualitative claim survives, break-even 198 → 148 |
| §10 Theorem A | unchanged, unaffected by the cost-model error (it is about centres, not costs) |
| §8 instance-side | unchanged: **LIVE, UNSETTLED, still rank 1** — and now *more* valuable, since the rho crossover moved from `w ≈ 56` to `w ≈ 104` |
| everything else | unchanged |

**My least-secure verdict is now §9.12** (one estimate, `d_reg`, not measured), then §6 (rests on the
shape of the literature rather than a computation I ran), then the **memory-aware costing of MITM**,
which is not a verdict at all but an unmeasured input that my own corrected table depends on.

**Standing caveat, restated because a retraction is exactly when it matters:** none of this is an
infeasibility claim about the instance. Theorem D is a generic-model, average-case bound with the
encoding knob explicitly excluded; Theorem B (corrected) prices ball-covering MITM with unbounded
memory; Theorem C covers affine maps, which is all the generic model can realise. `w` remains
unknown, no non-vacuous upper bound is established, and §8 is open.
