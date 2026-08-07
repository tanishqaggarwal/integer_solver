# STRUCTURED_KEYS — the cheap structural families of `k₀`, priced and swept

**Agent AE.** Everything below was measured in this directory. Scripts: `ae_setup.py`,
`ae_lib.py`, `aekang.c`, `ae_plant.py`, `ae_families.py`, `aequot.c`, `ae_quot.py`,
`ae_report.py`, `ae_verify.py`, `ae_support.py`. Nothing is quoted from another agent without
re-derivation except where the row says otherwise.

**If any run in this directory ever produces a `CAND` line, do not believe it until
`python3 ae_verify.py <k>` agrees.** That verifier shares no code with the search: Jacobian
coordinates instead of affine, a fixed-window scalar multiplication instead of binary
double-and-add, extended-Euclid inversion instead of Fermat, and curve parameters re-read from
`agentX_work/xdata.json` rather than from this directory's `ae_data.json`. `ae_verify.py
--selftest` proves it against the instance ladder and a planted scalar and **passes**; it was
built before any hit could occur, not after one.

---

# 0. THE ONE-PARAGRAPH VERSION

Every family here is a **magnitude or quotient** family, not a weight family. A **hit** in any of
them solves the instance outright and hands you `w` exactly. A **miss** excludes a set of scalars
and says **nothing whatever about `popcount(k₀)`** — the exclusions below are statements about
where `k₀` is **not**, not about how heavy it is. That is the whole honest content, and it is
what agent AB's §9.13 already said before I started. The reason to run them anyway is that they
are cheap and that nobody had. **Result: a clean sweep of misses. No hit anywhere, in any family,
at any radius.**

Three numbers the campaign can take away:

1. **`k₀ ∉ [0, 2^58)` at 98.2 % confidence**, for `2^32` group operations and **12.6 MB**, in 56
   minutes on one core. That extends agent X's deterministic `k₀ > 2^52` by six bits and is the
   best keys-per-operation ratio any instrument on this fleet has bought. It is still `2^198`
   short of mattering.
2. **`k₀ ≢ ±λ^{-e}·a·b^{-1} (mod N)` for every `1 ≤ a,b ≤ 2^26` and `e ∈ {0,1,2}` — exhaustively,
   not probabilistically.** That family has `2^53.87` members, which is **larger than agent X's
   entire weight-`≤9` class (`2^53.38`)**, and it cost `2^27` point operations against X's
   `2^33.1`. It is the cheapest large exhaustion this campaign has run, and it is disjoint from
   every weight-based family.
3. **On this box more threads make the search slower in total throughput, not merely
   sub-linearly** (§7). Anyone budgeting compute here from a per-core figure × 4 will be wrong by
   ~1.7× in the wrong direction.

None of it moves `w`. Under agent AC's posterior the total prior mass removed by everything in
this file is `2^-197.9`.

---

# 1. VERIFIED BEFORE STARTING

Deliverable baseline, re-checked at the end of this session, unchanged by anything here:
`python3 solve_lab/checker.py solve_lab/best/new_instance_partial_39026.json`
→ `satisfied 39026/39033 (7 failing) [12231, 12270, 12350, 14584, 18673, 22044, 29125]`.
**I did not beat it and nothing in `agentAE_work/` is a partial assignment.**

Curve data, re-derived and re-verified here (`ae_setup.py`):

| check | result |
|---|---|
| `p == 2^256 − 2^32 − 977` | **True** |
| `N ==` secp256k1 group order | **True** |
| `a == 0`; `G`, `T` on the curve | **True** |
| `N·G == O`, `N·T == O` | **True** |
| `ladder[i] == 2^i·G` for **all** `i = 0..255`, independent doubling chain | **256/256, 0 bad** |
| all 256 ladder points on the curve | **True** |
| `G`, `T` identical to agent **Y**'s independently extracted `ydata.json` | **True** |
| GLV: `β³ ≡ 1 mod p`, `λ³ ≡ 1 mod N`, `φ(x,y) = (βx, y) == [λ]` on 12 points | **True** |

`ae_data.json` is written only if every check passes; the script exits non-zero otherwise.

---

# 2. THE INSTRUMENT — a distinguished-point kangaroo, and what it does and does not prove

`aekang.c`: secp256k1 field arithmetic on 4×64-bit limbs, affine walks, one Montgomery batched
inversion per lock-step round over `K` kangaroos, jump index = low 5 bits of `x`, distinguished
point = low `d` bits of `x` zero, shared open-addressed DP table. **Memory is the DP table alone
and does not grow with the range**: `#DP = jumps / 2^d` at 48 B each, which was **12.6 MB** in
every production run here including `R = 58`, because `d` is chosen per run to hold `#DP` near
`2^18`. Half the herd is tame (start `u·G`, `u` uniform in `[0,L)`), half wild (start `Q + v·G`,
`v` uniform in `[0,L/8)`); mean jump `2^{R/2}`; `K = 1024` kangaroos.

## 2.1 It is **not** an exhaustion. This matters and is stated everywhere below.

BSGS decides `k₀ ∈ [a, a+L)` **deterministically**. A kangaroo does not: a miss means "no
tame–wild collision within the jump budget". The honest form of every miss below is therefore
**"excluded at confidence 1 − ε"**, with ε computed, not asserted.

## 2.2 Calibration — 34 planted keys, every one recovered **exactly**

A plant is a uniform `k` in `[0,2^R)`; PASS requires the engine's candidate to be confirmed by an
**independent Python bignum recomputation** of `k·G` *and* to equal the planted `k`.

| plant set | n | recovered exactly | mean `jumps/√L` | max |
|---|---|---|---|---|
| `R = 32` (`ae_plant.py 32 12 2 32`) | 12 | **12/12** | 2.406 | 7.29 |
| `R = 40` | 12 | **12/12** | 1.946 | 4.72 |
| `R = 48` | 6 | **6/6** | 1.546 | 3.13 |
| **shifted interval** `k = c + j`, random 256-bit `c`, run on `T − c·G` | 4 | **4/4** | 1.42 | 2.85 |
| **total (engine)** | **34** | **34/34** | **1.976** | 7.29 |
| **family-builder plants** — planted inside the real `c_N_over_2`, `c_pow2_128`, `c_two256_modN`, `c_lam` windows, run through `ae_families.build()` and its `lo`/reconstruction path | 4 | **4/4** | — | — |
| quotient-engine plants (`ae_quot.py plant`), incl. 4 deliberate lane-0 plants | 8 | **8/8** | — | — |
| **grand total** | **46** | **46/46** | | |

The family-builder row matters separately from the engine row: it exercises the code that
computes each centre's `lo`, shifts the target, and reconstructs `k = lo + j`. An engine that
works perfectly against a wrong `lo` would produce a confident, silent, universal miss.

**Closed-form check #1 (rule 5).** The theoretical mean is `2√L`. Measured over 34 plants:
**1.976** — this is the number that says the engine is not merely running but running *correctly
at the right cost*.

**Closed-form check #2 (rule 5).** Every run reports `dps` and the closed form `jumps/2^d`.
Across every run in this document the ratio is **0.99–1.01**. A run whose DP count does not match
its own jump count over `2^d` has not searched what it claims to have searched. A crashed or
table-less run cannot produce this ratio; it produces no `DONE` line at all.

**Closed-form check #3.** The engine self-tests before every run: it re-derives `2·L_i == L_{i+1}`
for 8 ladder indices and `jd_i·G == J_i` for all 32 jump points from the ladder, and exits **3**
with `SELFTEST_FAIL` otherwise. Every run below printed `SELFTEST_OK`.

## 2.2a A bug the plants missed and a counter caught — worth the space

The first version of the quotient engine (`aequot.c`) reported `deg1 = 1023` and `deg2 = 1023`
degenerate steps on a table of `W = 1024` lanes with `chunk = 1024`. **Four random planted
quotients passed anyway** — `(a₀,b₀)` were large, and the damage was confined to lane 0.

`1023 = chunk − 1` is not a plausible count for a `2^-128` event. Reading it against its closed
form: lane 0 starts at `1·G`, so its first step computes `dx = x(G) − x(G) = 0` — an addition
that is really a **doubling** — the batch-inversion path skipped it, the lane never advanced, and
**every entry of lane 0's chunk was recorded against `x(G)`** — `a ∈ [1, 2^10]` in the `m = 20`
plant runs where it was caught, and it would have been `a ∈ [1, 2^15]` (and likewise `b`) in the
`m = 25` production run. About 0.2 % of the pairs — but including `a = 1`, i.e. the **entire
`k₀ = 1/b` sub-family**, which is the most interesting part of the family and the one that
generalises agent X's row 5.

Fixed by supplying `2G` and `2T` and resolving the `P == base` case explicitly (the `P == −base`
case aborts with a non-zero exit; it is unreachable here). After the fix the counter reads
**`deg1 = deg2 = 1` exactly** — the one genuine doubling — which is itself the closed form. Four
**deliberate lane-0 plants** (`(a,b) = (1,3), (2,1), (1,2^{19}), (5,7)`) were added and are part of
the standing validation.

> **The random plants could not have found this and did not. The counter found it.** That is the
> case for rule 5 stated as concretely as I can state it.

## 2.3 The failure model, measured

| `c` | plants needing `> c√L` | empirical `P` | `exp(−c/2)` |
|---|---|---|---|
| 2 | 13/34 | 0.382 | 0.368 |
| 4 | 4/34 | 0.118 | 0.135 |
| 6 | 1/34 | 0.029 | 0.050 |
| 8 | 0/34 | 0.000 | **0.018** |

The exponential model with mean `2√L` fits (sd 1.53 vs mean 1.98). **A run capped at `8√L`
therefore misses a genuine hit with probability ≈ 1.8 %; capped at `4√L`, ≈ 13.5 %.** Each family
row states its cap, so each exclusion carries its own ε.

## 2.4 Build hygiene

Both binaries are compiled with the exit status checked directly (`gcc …; rc=$?`), never through
a pipe — a pipeline's status is the last stage's, which is how a failed link on this fleet once
left a stale binary in place. `aekang` and `aequot` were both rebuilt this way from unchanged
source after the first plant runs and re-verified to execute.

---

# 3. WHAT A MISS ON A MAGNITUDE FAMILY GIVES — the direction of every implication

This section exists because getting it wrong is the main way this thread could produce a
falsehood.

| family | a **hit** proves | a **miss** proves | does the miss bound `w`? |
|---|---|---|---|
| `k₀ < 2^R` | `k₀` known ⇒ `w` known exactly; and `w ≤ R` | `k₀ ≥ 2^R`, i.e. `max(S) ≥ R` | **NO** |
| `N − k₀ < 2^R` | `k₀` known ⇒ `w` known | `N − k₀ ≥ 2^R` | **NO** |
| `\|k₀ − c\| < 2^{R−1}` | `k₀` known ⇒ `w` known | `k₀` is not within `2^{R−1}` of `c` | **NO** |
| `k₀ = ±λ^{-e}·a·b^{-1}` | `k₀` known ⇒ `w` known | no such representation with `a,b ≤ 2^m` | **NO** |
| **NAF weight ≤ m** (agent X) | `k₀` known | `NAF-wt(k₀) > m`, and since `w ≥ NAF-wt`, only `w > m` | **only from BELOW** |

**The NAF row is the trap the brief warned about, and the check comes out the way the brief said.**
The non-adjacent form is the *minimal* signed-digit representation, so `NAF-wt(k₀) ≤ w` always.
A signed-digit exhaustion therefore pushes `w` **up**, never down: agent X's `m ≤ 7` exhausted
yields `w ≥ 8`, which is strictly **weaker** than X's own unsigned `w ≥ 10`. That is exactly why
the two numbers differ, and it is why **no signed/NAF experiment was run here** — the implication
points away from the upper bound this thread was asked for. Priced at zero, correctly.

**So: every row in §4 is a lottery ticket whose only payoff is the hit branch.** Under agent AC's
posterior (`w ∈ [115,141]` at 90 %), the hit branch is where all of the value is, and there is
almost none of it.

---

# 4. THE FAMILY SWEEP — priors, costs, outcomes, exclusions

Priors are `|F| / N` with `log₂ N = 256.0000`. Every kangaroo row ran to its **full jump cap**
(a miss always costs the whole budget; only a hit stops early), so `jumps == cap` in every row is
itself a check that the run was not truncated.

## 4.1 The rows

| family | keys covered | prior P(k0 in F) | ops (jumps) | DP/closed-form | exclusion confidence | outcome |
|---|---|---|---|---|---|---|
| magnitude k0 < 2^58 (headline) | 2^58.00 | 2^-198.0 | 2^32.0 | 0.9986 | 98.2% | miss |
| magnitude k0 < 2^60 (abandoned) | 2^60.00 | 2^-196.0 | 2^28.1 | None | 12.9% | NOT RUN (partial, 13% conf) -- reported as not-run |
| magnitude k0 < 2^64 (abandoned) | 2^64.00 | 2^-192.0 | 2^27.0 | None | 1.5% | NOT RUN (partial, 2% conf) -- reported as not-run |
| const tier (37 families, R=44) | 2^49.21 | 2^-206.8 | 2^30.2 | 0.9954-1.0031 | 98.2% | miss |
| orbit tier (5 families, R=48) | 2^50.32 | 2^-205.7 | 2^29.3 | 0.9973-1.0028 | 98.2% | miss |
| quotient a*b^-1, m=25 (EXHAUSTIVE) | 2^51.87 | 2^-204.1 | 2^26.0 | n/a | 100% (exhaustive) | exhausted-no-hit |
| quotient a*b^-1, m=26 (EXHAUSTIVE) | 2^53.87 | 2^-202.1 | 2^27.0 | n/a | 100% (exhaustive) | exhausted-no-hit |

TOTAL keys excluded by rows that completed: 2^58.11  (prior mass 2^-197.9)
For comparison: agent X weight<=9 covered 2^53.38; agent Y complement 2^53.38.

Regenerate with `python3 ae_report.py`; it reads only the engines' own `DONE` lines and the tier
JSONs, and prints `NOT RUN` for anything without one. The `m = 25` quotient row is nested inside
`m = 26` and is kept only because it is an independent second run of the same instrument.

**Every kangaroo row shows `jumps` exactly equal to its cap** (`2^32`, `37 × 2^25 = 2^30.21`,
`5 × 2^27 = 2^29.32`) — a miss spends the whole budget, so this equality is the check that no run
was cut short. **Every `DP/closed-form` is within 0.5 % of 1.** **Every `dxzero` is 0.**

## 4.2 What each miss excludes — as a statement about `k₀`, not about the search

1. **`k₀ ≥ 2^58`.** Equivalently `max(S) ≥ 58`: the highest ON selector index is at least 58.
   Confidence 98.2 %. This supersedes agent X's row 3 (`k₀ > 2^52`, BSGS, deterministic) by six
   bits — but X's is an exhaustion and mine is not, so **both should be cited, not just mine.**
2. **`λ·k₀`, `−λ·k₀`, `λ²·k₀`, `−λ²·k₀` are each `≥ 2^48`** (98.2 % each) — i.e. `k₀` is not
   `±λ^{∓1}` times a scalar below `2^48`. This is the **magnitude** version of the orbit question;
   agent Y's orbit work was the **weight** version, and the two exclude different sets.

   **The fifth family in that tier, `N − k₀ ≥ 2^48`, adds nothing and I am not claiming it.**
   Agent X's row 3 already has `N − k₀ > 2^52` **deterministically**, which is four bits stronger
   *and* exhaustive. I ran it because it was one line in the same tier and because it is the
   control that proves the tier's `lo`-shifting works against a target whose answer is known to
   be "no"; **cite X, not me, for the top end.** Symmetrically, my `R = 58` low-end result is the
   only magnitude row here that improves on X, and only at the bottom.
3. **`k₀` is not within `2^43` of any of 37 named constants** — `N/d` for
   `d ∈ {2,3,4,5,6,7,8,10,16,100}`, `2N/3`, `3N/4`, `d^{-1} mod N` for `d ∈ {5,7,11,13}`, `2^e`
   for `e ∈ {16,64,96,128,160,192,224,250,254,255}`, `2^256 mod N`, `λ`, `λ²`, `β mod N`, and
   seven rational multiples of `2^256 − 1` (`(2^256−1)/d` for `d ∈ {3,5,15,17,257}` and
   `2(2^256−1)/3`, `4(2^256−1)/5` — the `0x5555…`, `0x3333…`, `0xAAAA…`, `0xCCCC…` patterns).
   Confidence 98.2 % each. Six further centres were **skipped as duplicates** because they fall
   inside another centre's window, and two of those are worth recording as arithmetic facts:
   `3^{-1} mod N` lies within `2^44` of `⌊2N/3⌋`, and `−λ` lies within `2^44` of `λ²` — the latter
   because `1 + λ + λ² ≡ 0`, so `−λ = λ² + 1` exactly.
4. **`k₀ ≢ ±λ^{-e}·a·b^{-1} (mod N)` for any `1 ≤ a,b ≤ 2^26`, `e ∈ {0,1,2}`.** This one is
   **exhaustive, not probabilistic**: the engine stores the low 64 bits of `x`, and truncation can
   only manufacture false positives, never suppress a true match — so **0 matches means 0 true
   matches, exactly.** It contains `b = 1` (`k₀` small), `a = 1` (`k₀ = 1/b`, generalising agent
   X's row 5 from `2^i/m` to `a/b`), and every "nice rational" scalar.

## 4.3 What none of it excludes

**Nothing about `w`.** Under agent AC's posterior (`w ∈ [115,141]` at 90 %) and agent AG's
corrected dead band `[53,200]`, the total prior mass removed by everything in §4.1 is
**`2^-197.9`** — which is, to the precision anyone cares about, zero. For scale, agent AC put the
whole campaign's total-variation movement at `2^-201.6`; this file roughly triples that and it
still rounds to nothing.

The campaign gains the sentence *"the obvious cheap structural shortcuts were checked and are not
there"*, and that is the entire return. It is worth having and it is not progress.

---

# 5. THE INSTRUMENT COMPARISON — why kangaroo, and where it stops

Coverage per operation, on this box (higher is better):

Kangaroo memory is `#DP × 48 B`, and `#DP = jumps / 2^d` is a **tunable** — raise `d` and the
table shrinks while only the end-of-run tail lengthens. The figures below are the DP storage the
runs actually needed; I allocated `2^21` slots (100 MB) out of caution.

| instrument | memory | coverage | ops | keys/op | exhaustive? |
|---|---|---|---|---|---|
| agent X's weight `≤ 9` sweep | 1.9 GB table | `2^53.4` | `2^33.1` | `2^20.3` | **yes** |
| kangaroo `R = 44` (run, 37×) | 12.6 MB | `2^44` | `2^25` | `2^19.0` | no (ε = 1.8 %) |
| kangaroo `R = 48` (run, 5×) | 12.6 MB | `2^48` | `2^27` | `2^21.0` | no (ε = 1.8 %) |
| **kangaroo `R = 58`** (run) | **12.6 MB** | **`2^58`** | **`2^32`** | **`2^26.0`** | no (ε = 1.8 %) |
| kangaroo `R = 60` | 12.6 MB | `2^60` | `2^33` | `2^27.0` | no |
| kangaroo `R = 64` | 12.6 MB | `2^64` | `2^35` | `2^29.0` | no |
| **quotient sweep `m = 25`** (run) | **512 MB** | **`2^51.9`** | **`2^26`** | **`2^25.9`** | **yes** |
| **quotient sweep `m = 26`** (run) | **1 GB** | **`2^53.9`** | **`2^27`** | **`2^26.9`** | **yes** |
| quotient sweep `m = 27` | 2 GB | `2^55.9` | `2^28` | `2^27.9` | yes — **memory wall here** |
| BSGS, whole 2 GB budget | 2 GB | `2^56` | `2^29` | `2^27.0` | **yes** |

**The memory column is the point.** The kangaroo's is flat in `R`; every other instrument's grows
with coverage until it hits the 2 GB wall.

Four things this table settles:

1. **The brief's headline observation is right in direction and optimistic in size.**
   Plain BSGS with `2^30` entries would give `R = 60` at `~2^30` work — but `2^30` entries is
   **8 GB** at the 8 bytes/entry the brief assumed (and 8 bytes cannot hold a usable fingerprint
   *plus* a 30-bit index, so it is really 16). This box has ~8 GB *available* shared across the
   whole fleet against a 2 GB budget for me. **At 2 GB the deterministic reach is `R = 56`, not
   60**, and the `2^30` work figure ignores that the run must also *build* and *sort* the table.
2. **The kangaroo reaches further than any table method here precisely because its memory does
   not grow with the range** — 12.6 MB against 2 GB, and the same 12.6 MB would serve `R = 64`.
   That is the trade the brief anticipated, and it is the right one on a box whose binding
   constraint is memory and whose second constraint is that *more threads make it slower* (§7).
3. **The quotient sweep is the best buy on the board at low cost** and it is genuinely
   exhaustive. Both it and the kangaroo scale as coverage ∝ cost²; the quotient sweep's constant
   is ~3.6× better (it gets `±λ^e` — a factor 6 — for two field multiplications per probe), but it
   is capped by memory at `m ≈ 26–27` while the kangaroo is capped only by time. So:
   **quotient sweep first, kangaroo for reach.** Measured: `m = 26` covered `2^53.87` scalars
   exhaustively in 1,184 s; matching that coverage with the kangaroo needs `R = 54`, `2^31` jumps,
   ~1,700 s, and only at 98.2 % rather than exactly.
4. **Disk does not rescue BSGS.** Its access pattern is random; agent AB measured the penalty at
   ~`2^20` slowdown for a `2^1` memory gain. Not attempted.

---

# 6. FAMILIES PRICED AND **NOT** RUN — with the reason

Reporting these is the point of the exercise as much as the runs are.

| family | why not run |
|---|---|
| `k₀ = a·b`, `a ≤ 2^20`, `b ≤ 2^B` | **Dominated, and provably so.** `a·b < 2^{20+B}`, so the whole family is a **subset** of the magnitude family `k₀ < 2^{20+B}`, which the kangaroo decides at `2^{(20+B)/2+4}` ops. Reaching `B = 36` by multi-target BSGS costs `2^32` and covers `2^56` keys — the kangaroo covers those same `2^56` keys for `2^32` and covers `2^60` besides. The brief's suggestion is right that it is cheap; it is cheap **and redundant**. The non-redundant version — `a` small, `b` *not* magnitude-bounded — is the quotient family, which **was** run (§4). |
| `k₀` with low NAF / signed weight | Implication points the wrong way (§3), and agent X exhausted `m ≤ 7` already. |
| window family `k₀ ≡ a·2^s (mod N)`, `\|a\| < 2^{R−1}`, 256 shifts | **Built, priced, NOT run.** It would extend agent X's row 4 (`a < 2^34`, integer enumeration, no wraparound) to `a < 2^39` **including** the mod-`N` wraparound cases X could not see — genuinely new, but dominated: 256 shifts at `R = 40` cover `2^48` keys for `2^30` ops, whereas the same `2^30` spent on one kangaroo buys `R = 52`, i.e. `2^52`, sixteen times more. Its only argument is disjointness, and at `2^-208` prior that argument does not pay for an hour of a contended box. The driver is in `ae_families.py` (`tier = 'window'`) and runs with one command if anyone wants it. |
| magnitude analogue of agent Y's complement sweep | **The obvious version is empty, and this is worth recording.** If the ON-set folds to `k₀` itself, its integer complement `2^256 − 1 − k₀` is `≥ 2^256 − N = 2^128.346` for every `k₀ ∈ [0,N)` — so "the complement of `k₀` is small in magnitude" is **impossible in that branch, at any `R < 128.346`**, recomputed exactly. It is reachable only when the ON-set folds to `k₀ + N`, which needs `k₀ < 2^256 − N`. Working that through: a complement below `2^R` forces `k₀ ∈ (2^256 − 1 − N − 2^R, 2^256 − N − 1]` — an interval of width `2^R` sitting immediately below `2^256 − N`. **That is (up to its single lowest point) the `c_two256_modN` centre in §4**, so the family is covered — but it must not be called "the complement family": it is a magnitude window living entirely inside the `2^-127.7` two-solution branch, and it has nothing to do with Hamming distance. |
| deterministic BSGS at `R = 56` | Subsumed in coverage by the completed `R = 58` kangaroo; its only advantage over that is exhaustiveness (100 % vs 98.2 %), bought at 2 GB and heavy random-access pressure on a box that spent the session between load 14 and load 45. Priced (§5), declined. **If the campaign wants the `[0, 2^56)` interval closed deterministically rather than at 98.2 %, this is the way and it is about 30 minutes of CPU plus 2 GB.** |
| any test of `k₀` against hashes, timestamps, PRNG outputs or file emission order | **Forbidden by rule 2 and not attempted.** See §9. |

---

# 7. AN OPERATIONAL MEASUREMENT WORTH CARRYING

Thread scaling of the same 2^26-jump workload, measured back-to-back on this box (load average
14–20 on 4 cores):

| threads | wall | rate |
|---|---|---|
| **1** | 35.1 s | **1.915 M jumps/s** |
| 2 | 57.3 s | 1.172 M jumps/s |
| 4 | 61.1 s | 1.100 M jumps/s |

**More threads made it strictly slower in total throughput, not merely sub-linear.** Every
production run below is therefore single-threaded. Any agent planning a compute budget on this
box from a per-core figure times four will be out by a factor of ~1.7 in the wrong direction.

---

# 7a. THE RUN LEDGER — everything started, including what was thrown away

A status marker is a claim. This is the full list, so nothing in this directory can be mistaken
for evidence it is not.

| artefact | what it is | standing |
|---|---|---|
| `head58.out` / `head58.err` | the headline `R = 58` run | **evidence** |
| `const44.log` → `res_const.json` | 37-family named-constant tier at `R = 44` | **evidence** |
| `orbit48.log` → `res_orbit.json` | 5-family orbit/top tier at `R = 48` | **evidence** |
| `res_quot_25.json`, `res_quot_26.json` | exhaustive quotient sweeps | **evidence** |
| `calib_R32.json`, `calib_R40.json`, `quot_plant.log` | the plant records | **evidence** |
| `head64_PARTIAL_ABANDONED.txt` | `R = 64` started 21:37, stopped for budget at `2^27.0` jumps (**1.5 %** of the exclusion it would need) | **NOT evidence — reported as not-run** |
| `head60_PARTIAL_ABANDONED.txt` | `R = 60` started 21:52, stopped for budget at `2^28.1` jumps (**12.9 %**) | **NOT evidence — reported as not-run** |
| `head60b.out` / `head60b.err` (seed 424242) | a **second** `R = 60` attempt, launched 23:36 once the box emptied, budget `2^33` jumps (`8√L`, ε = 1.8 %), ~85 min at the 1.7 M/s the free box delivers | **in flight at hand-off.** Read it the same way as every other row: a `DONE … cands=0` line with `jumps = 4294967296×2` and `dps/dpexp ≈ 1` means `k₀ ∉ [0,2^60)` at 98.2 %; a `CAND` line means run `python3 ae_verify.py <k>`. **A `STATUS` line is not a result** — if it never reached `DONE`, quote §4's `R = 58` row and nothing more. `finish_R60.sh` (watcher pid recorded in `AE_PIDS.txt`, waiting on the engine **by PID**) writes **`R60_VERDICT.txt`** when the engine exits: it re-checks `jumps == 2^33`, `dps/dpexp ≈ 1` and `cands == 0` before writing a verdict, and refuses to write one otherwise. |
| `DISCARDED_const_R44_double_launch.log` | two `ae_families.py` processes wrote to it concurrently — my launch error, the file is interleaved | **NOT evidence**, kept only so the mistake is on the record |
| `DISCARDED_const_R48_partial.log`, `const_R48_PARTIAL_5of37.txt` | 5 of 37 constants completed at `R = 48` before the box load tripled; superseded by the complete `R = 44` tier | **not cited**; the 5 rows are real but the tier is not, and a mixed-`R` tier is not a tier |
| the first `aequot` binary and its 4-plant run | had the lane-0 bug of §2.2a | **destroyed and rebuilt**; no result from it is quoted |
| `res_quot_25.json` as it existed at 22:21 (`rc = -15`) | output of the buggy engine, killed | **deleted**, so it could not be picked up by the report generator |

Two abandoned runs and one interleaved log are the price of re-planning a compute budget three
times on a box whose load went from 14 to 45 and back to 20 during the session. **The re-planning
was right; publishing a 1.5 %-confidence run as an exclusion would not have been.**

---

# 8. THE COORDINATOR'S ITEM 3 — U's partition theorem and the 927 lift conditions

Full working in **`AE_ITEM3.md`**; scripts `ae_support.py`, outputs `res_support.json` and
`res_927_supports.json`. Headline, so it is not lost behind a pointer:

* The conjectured chain *"no wraparound ⇒ exact integer addition ⇒ the lift conditions are
  `|S|`-blind"* **runs backwards at the second step.** Exactness makes the fold value at a node an
  **injective** encoding of `S` restricted to that node (distinct powers of two ⇒ unique binary
  representation), so the conditions there see `S` *perfectly*. U's theorem enables the mechanism;
  it does not refute it.
* What does limit §8 is **support locality, measured**: of the 927 `c > 1` conditions, **925 have
  a selector support that is a proper subset of `{0..255}`**, **793 see ≤ 4 selectors**, **48 see
  none at all**, and **only 2 see all 256**. The 879 with non-empty support cover **exactly** U's
  511-node tree, set-equal.
* Locality alone does **not** close §8 — a conjunction of subtree-local conditions over a laminar
  family can still bound `|S|` — so the clean "kills it for all `|S|` at once" result **does not
  follow**, and I do not claim it.
* It does collapse the search to **two named conditions** (`u = x5146, h = x29804, c = 6672769`
  and `u = x14393, h = x34243, c = 12354891`). Whether those two vary with `|S|` is a **probe**,
  which is agent T's and agent AH's ground; route it there.
* **Scope flag:** this question referenced a crux and an object that were never handed to me
  (my brief was the curve-side sweeps). I answered it because it was short, but the measurement
  is a **join of F's parse, L's handle list and U's support closure** — each re-checked here, none
  re-derived from `EQUATIONS.txt`.

---

# 8a. IF SOMEONE PICKS THIS UP

* **Cheapest thing left that is not dominated:** the quotient sweep at `m = 27` — `2^28` point
  operations, ~40 min, **2 GB of RAM** (transient, not disk), covering `2^55.9` scalars
  *exhaustively*. I did not take that 2 GB unilaterally. One command:
  `python3 ae_quot.py real 27`.
* **If the campaign wants `[0, 2^56)` closed deterministically** rather than at 98.2 %: BSGS, 2 GB,
  ~30 min (§5). Not built here.
* **The window tier is built and priced and I recommend against it** (§6):
  `python3 ae_families.py window R_window=40 log2max=23`.
* **Everything scales as coverage ∝ cost².** At the 1.7 M jumps/s this box delivers when it is
  quiet, `R = 64` at `8√L` is `2^35` jumps ≈ 5.6 h on one core. There is no cleverness left in the
  instrument; only hours.

---

# 9. RULE-2 BOUNDARY — where I stopped

Everything in this file is arithmetic on a point of an elliptic curve: is `T` a small multiple of
`G`; is it near a named constant; does it factor through a small rational. Every centre used in
§4 is definable from the curve and field data alone (`N`, `p`, `λ`, `β`, powers of two, and
elementary rational multiples of `2^256 − 1`, the constant agent Y's complement identity already
put in play). **No test of provenance was run and none was designed:** no PRNG traces, no
hashes of strings or timestamps, no examination of emission order or coefficient templates, no
hypothesis about who wrote the file or how. Where a family looked like it was drifting from
"low-complexity integer" toward "plausible authoring choice", it was dropped rather than run.
