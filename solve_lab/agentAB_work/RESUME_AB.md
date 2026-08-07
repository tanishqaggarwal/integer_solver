# RESUME_AB — agent AB, upper-bound theory. Self-contained.

**Mandate:** enumerate every mechanism that could bound the solution's Hamming weight `w` FROM
ABOVE, settle each DEAD or LIVE. Theory, not search. Deliverable: `UPPER_BOUND_MAP.md` (this dir).

**Status: COMPLETE (round 2). ROUND-1 THEOREM B NUMBERS ARE RETRACTED — see below.** No compute process launched, no long job running, nothing to resume mid-flight.
Everything below is reproducible from cold in < 3 minutes.

## Files in this directory
| file | what it is | run time |
|---|---|---|
| `UPPER_BOUND_MAP.md` | **the deliverable** — 9 numbered mechanisms + 7 extras, 2 theorems, ranked live list | — |
| `ab_facts.py` | independent curve/G/T/ladder re-verification; complement identity incl. a planted weight-250 end-to-end test; N's binary structure + exact digit-DP for `max popcount(k<N)`; null tail probabilities; λ and negation weight measurements; halving bijection; counting | ~50 s |
| `ab_cost.py` | **Theorem B** ball-covering cost floor; Theorem A masked-complement correlations (4000 samples); residual curve checks (primality, anomalous, CM, embedding degree, `b/7` a 6th power) | ~60 s |
| `ab_rank.py` | round-1 break-even + cost/informativeness ledger (**numbers superseded by `ab_costfix.py`**) | ~20 s |
| `ab_barrier.py` | **round 2** — Theorem C (no weight-preserving self-map, proved) and Theorem D (generic-group bound for an arbitrary predicate of `k`) | ~40 s |
| `ab_costfix.py` | **round 2** — the retraction: corrected ball cost `√W·Vol₁₂₈(W/2)`, corrected Theorem B, corrected budget→weight table | ~15 s |
| `ab_soft.py` | **round 2** — §6 upgraded to a structural argument; §9.12 costed by degree of regularity | ~5 s |

Reproduce: `cd solve_lab/agentAB_work && PYTHONDONTWRITEBYTECODE=1 python3 ab_facts.py && python3 ab_cost.py && python3 ab_rank.py`.
Only external file read: `../agentX_work/xdata.json` (read-only, for `G`, `T`, ladder). Nothing
outside this directory was modified. No git commands were run.

Deliverable re-verified at start: `checker.py best/new_instance_partial_39026.json` → **39026/39033**,
failing `[12231,12270,12350,14584,18673,22044,29125]`.

## ⚠ RETRACTION (round 2) — my own Theorem B numbers

I priced a radius-`W` Hamming ball at `C(256,W/2)`. The correct per-ball cost is the **cumulative
half-volume** `√W · Vol₁₂₈(W/2)` (Coppersmith / Stinson splitting systems). Wrong by up to **2^65**
at large `W`. The sanity check I failed to run: at `W = 256` one ball is the whole key space, so the
cost must be 2^128; round-1's model returned 2^251.7, corrected returns 2^132.0.

* **Break-even moves `B = 198` → `B = 148`.** Largest affordable complement radius `W = 107`.
* **Qualitative claim survives intact:** `w ≤ 148` is +2.6σ on the null, excludes ~0.4% of null mass
  (~0.006 bits), and costs 2^126.4 while **solving costs 2^126.5 and returns `w` exactly**.
* **Consequence the fleet must act on:** `MINIMUM_COST_SEARCH.md` §4 uses `C(256,w/2)` and is
  pessimistic. Corrected (**time-only, unbounded memory**): 2^47 → `w ≤ 18` (not 14); 2^58 → `w ≤ 24`
  (not 20); 2^80 → `w ≤ 40` (not 30). **The rho crossover moves from `w ≈ 56` to `w ≈ 104`**, which
  nearly doubles the payoff band for §8.
* **Caveat stated, not buried:** memory binds. ~2^30 entries on this box ⇒ `w ≤ 10`. A memory-aware
  costing of low-weight MITM has not been done by anyone and is a cheap concrete follow-up.

## The three results worth carrying forward

1. **THEOREM B (new, exact combinatorics).** Every search-based upper bound is a Hamming-ball
   covering of `{wt > B}`. Cost floor `= min_W (|{wt>B}|/Vol(256,W))·C(256,W/2)`. **Break-even at
   `B = 198`:** proving `w ≤ 198` costs 2^126.3, and solving outright by rho costs 2^126.5 and
   returns `w` exactly. Proving `w ≤ 128` (the null median) costs **2^185.0**.
   ⇒ *Every search-based upper bound is vacuous (`B ≥ 198`, ≥ 8.75σ out) or dominated by solving.*
2. **THEOREM A (proof + measurement).** `c − k = c ⊕ k` iff `supp(k) ⊆ supp(c)`, so
   `c = 2^256 − 1` is the **unique** centre giving an unconditional upper bound. Agent Y's mechanism
   has no cheaper or stronger sibling. Masked variants give only bounds conditional on events of
   null mass `2^{|A|−256}`.
3. **§8 is the only live non-search mechanism and it is UNSETTLED — and the fleet is currently
   reading its evidence backwards.** Every confirmed integer-lift closure is at `|S| ≤ 64`
   (T: 1,2,3,5,6,7,8,17,32,64 → 39,018, identical 15-eq footprint). Those rule out *lower*-bound
   constraints; they are fully consistent with an upper-bound constraint `w ≤ B` for any `B ≥ 64`.
   The one high-`|S|` probe, `|S| = 128`, **stalled and gave up** —
   `agentT_work/t_close2wj_T128.log`: `outer 8: global nonzero 3 ... no addable collateral --
   giving up on ((x2820-x17195)-(8271997*x17079))`. Also: T's 32/64/128 ON-sets are nested prefixes
   of one `random.Random(7)` chain, i.e. **one correlated sample, not three.**

## Round 2: §2 upgraded from "no mechanism known" to a real barrier

* **THEOREM C — proved.** The only weight-preserving affine self-map of `Z_N` is the identity.
  `k=0 ⇒ b=0`; `k=1 ⇒ a=2^j`; then `k = 2^{256−j}` gives `ak = 2^256 mod N` with
  **popcount 65 ≠ 1**, forcing `j=0`. Verified for all 255 `j`. In the generic model affine is all
  an algorithm can realise, so **no weight-preserving randomised self-reduction exists** — the
  hardcore-bit machinery has nothing to run on. Measured: the only maps with *any* preservation are
  `k ↦ 2^j k`, preserving on a fraction `2^{−j}` (coincidence baseline 0.035).
* **THEOREM D — the barrier.** Generic group model: every held element is `σ(α_i+β_i k)`, collisions
  are single affine equations over the field `Z_N`, so for any predicate `P` and `m` queries,
  `Adv ≤ m²/(2·min(|D₀|,|D₁|))`. For `P = [w ≤ B]`: **`B = 128` needs `m ≥ 2^127.5`** (solving costs
  2^126.5) — **deciding the weight predicate is as hard as solving**; `B = 20` needs only `m ≥ 2^49`,
  and corrected MITM achieves 2^50.0, within 2^1 of optimal. **The lower/upper asymmetry is class
  size, now matched by a lower bound rather than a cost model.**
  *Knobs:* generic model (coordinate encoding excluded), average-case over a distribution on `k`.
  HGJ/BCJ representations are blocked because they need `k₀ mod M` — the same obstruction as §3.

## Verdict table (full detail in UPPER_BOUND_MAP.md)
| # | mechanism | verdict |
|---|---|---|
| 1 | complement identity (agent Y) | **LIVE** — sound, `w ≤ 255−W` at cost `C(256,W/2)`; vacuous vs the null; value is the hit branch only |
| 2 | bit security / HNP | DEAD — exact bit oracle ⇒ 256-query full solve (≥2^118.5); prime `N` kills the Legendre-symbol shortcut; HNP points the wrong way |
| 3 | lattice / LLL on the subset-sum | DEAD — density 1.000, but the real reason is that no integer target exists; given `k₀` the subset-sum is *trivial* |
| 4 | 2-adic / `v₂(k)` | DEAD — odd prime order ⇒ `[2]` is a bijection; proof, not heuristic |
| 5 | endomorphism `λ` | DEAD — `popcount(λk)` for weight-4 `k`: mean 127.41, sd 8.31. √3 accelerator only |
| 6 | character sums / analytic | DEAD as a bound; LIVE only as the null prior the campaign already has |
| 7 | counting / uniqueness | DEAD — but yields the free unconditional `w ≤ 255` (exact digit-DP over `k < N`) |
| 8 | instance-side constraint (agent Z) | **LIVE, UNSETTLED — rank 1** |
| 9 | 16 further mechanisms (PH, Smart, MOV, GHS, index calculus, CM/`j=0`, Cheon, `N`'s expansion, weight-preserving doubling, division polys, Gröbner, kangaroo intervals, multi-target, quantum, masked complement) | all DEAD except quantum (no hardware) |

Also settled in passing: **`b/7` is a 6th power mod `p`, so the instance curve is `F_p`-isomorphic to
secp256k1 `y² = x³ + 7`** — not merely an order-matching sextic twist. And **embedding degree
> 20 000**, extending the lab's verified `k ≤ 24`.

## Least-secure verdicts AFTER round 2 (re-ranked)
1. **§9.12 Gröbner** — now dead on complexity (`d = 16` Macaulay costs 2^197.1 at `ω=2.37`; affordable
   only if `d_reg ≲ 11`, estimate ≈ 23), and a term order cannot help since `d_reg` is a property of
   the ideal. But `d_reg` is **estimated, not measured**. One bounded Singular experiment
   (degree-truncated linear algebra on the reduced core) would replace the estimate with a number.
2. **§6 analytic** — upgraded from heuristic to structural (*a perfect equidistribution theorem says
   the class spreads uniformly over targets, i.e. says nothing about this `T`; the stronger the
   result, the less it says*), but it rests on the shape of the literature rather than on a
   computation I ran.
3. **The memory-aware costing of MITM** — not a verdict but an unmeasured input my own corrected
   table depends on. Nobody in the campaign has done it.

*(Round 1's #1, §2's weight predicate, is now closed by Theorems C and D above.)*

## THE next experiment (handed to agent T's tooling / agent Z's angle — I did NOT run it)
> Run `agentT_work/t_close2wj.py <tag> <n>` at **`|S| = 250`** first, then 192, then finish the
> stalled 128 — on ON-set seeds **independent** of the seed-7 chain. Confirm with `checker.py`
> (expect 39,018/39,033, 15-equation footprint). Budget ~30–60 min CPU per probe (T measured ~4 min
> at 32, ~6 min at 64). **`|S| = 250` is the complement regime, where the campaign holds zero data
> points and where an upper-bound constraint would first show.**
> * **Closes** ⇒ §8 dies; with Theorem B, no affordable upper bound on `w` exists at all — a clean
>   final answer to the question the fleet is spending compute on.
> * **Fails to close on independent seeds** ⇒ the campaign's first real upper bound; then bisect
>   downward (200/160/128/96), since only `B ≲ 56` beats rho and only `B ≲ 24` is actionable.

## One free improvement for agent X
`agentX_work/xsigned.c:107` loads exactly **256** ladder points → signed exponents 0..255 only.
`2^256 ≡ 2^256 − N (mod N)` has popcount 65, so **complement-sparse keys have no short signed
representation in X's table** — X's class and Y's class are disjoint, not redundant. Appending a
257th point `(2^256 − N)·G` to `data_real.txt` adds them at `m = wt(k')+2`. (Y's dedicated MITM stays
far cheaper for that class: 2^33 vs 2^44.4, so this is added coverage, not a replacement.)
