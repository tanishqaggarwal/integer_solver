# RESUME_AB — agent AB, upper-bound theory. Self-contained.

**Mandate:** enumerate every mechanism that could bound the solution's Hamming weight `w` FROM
ABOVE, settle each DEAD or LIVE. Theory, not search. Deliverable: `UPPER_BOUND_MAP.md` (this dir).

**Status: COMPLETE.** No compute process launched, no long job running, nothing to resume mid-flight.
Everything below is reproducible from cold in < 3 minutes.

## Files in this directory
| file | what it is | run time |
|---|---|---|
| `UPPER_BOUND_MAP.md` | **the deliverable** — 9 numbered mechanisms + 7 extras, 2 theorems, ranked live list | — |
| `ab_facts.py` | independent curve/G/T/ladder re-verification; complement identity incl. a planted weight-250 end-to-end test; N's binary structure + exact digit-DP for `max popcount(k<N)`; null tail probabilities; λ and negation weight measurements; halving bijection; counting | ~50 s |
| `ab_cost.py` | **Theorem B** ball-covering cost floor; Theorem A masked-complement correlations (4000 samples); residual curve checks (primality, anomalous, CM, embedding degree, `b/7` a 6th power) | ~60 s |
| `ab_rank.py` | break-even `B` for Theorem B; cost/informativeness ledger for the complement mechanism | ~20 s |

Reproduce: `cd solve_lab/agentAB_work && PYTHONDONTWRITEBYTECODE=1 python3 ab_facts.py && python3 ab_cost.py && python3 ab_rank.py`.
Only external file read: `../agentX_work/xdata.json` (read-only, for `G`, `T`, ladder). Nothing
outside this directory was modified. No git commands were run.

Deliverable re-verified at start: `checker.py best/new_instance_partial_39026.json` → **39026/39033**,
failing `[12231,12270,12350,14584,18673,22044,29125]`.

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

## Least-secure DEAD verdicts (my own, ranked)
1. **§2 for the *weight* predicate specifically** — the 256-query reduction covers exact **bit**
   oracles, not a direct low-cost decider for `w ≤ B`; and standard hardcore-bit machinery does not
   apply because `k ↦ k+r` scrambles Hamming weight. "No mechanism known", **not** "provably hard".
   Nobody has looked for a direct weight predicate.
2. **§6** — my "class too small for square-root cancellation" threshold is heuristic, not derived.
3. **§9.12 Gröbner** — the argument is empirical (the fleet attacked this system as equation repair
   and coset decoding, not as elimination with a term order).

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
